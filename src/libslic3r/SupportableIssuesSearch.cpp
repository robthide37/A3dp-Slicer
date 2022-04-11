#include "SupportableIssuesSearch.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <stack>

#include "libslic3r/Layer.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include "libslic3r/ClipperUtils.hpp"

#define DEBUG_FILES

#ifdef DEBUG_FILES
#include <boost/nowide/cstdio.hpp>
#endif

namespace Slic3r {
namespace SupportableIssues {

void Issues::add(const Issues &layer_issues) {
    supports_nedded.insert(supports_nedded.end(),
            layer_issues.supports_nedded.begin(), layer_issues.supports_nedded.end());
    curling_up.insert(curling_up.end(), layer_issues.curling_up.begin(),
            layer_issues.curling_up.end());
}

bool Issues::empty() const {
    return supports_nedded.empty() && curling_up.empty();
}

namespace Impl {

struct LayerDescriptor {
    Vec2f centroid { 0.0f, 0.0f };
    size_t segments_count { 0 };
    float perimeter_length { 0.0f };
};

struct EdgeGridWrapper {
    EdgeGridWrapper(coord_t edge_width, ExPolygons ex_polys) :
            ex_polys(ex_polys), edge_width(edge_width) {

        grid.create(this->ex_polys, edge_width);
        grid.calculate_sdf();
    }

    bool signed_distance(const Point &point, coordf_t point_width, coordf_t &dist_out) const {
        coordf_t tmp_dist_out;
        bool found = grid.signed_distance(point, edge_width, tmp_dist_out);
        dist_out = tmp_dist_out - edge_width / 2 - point_width / 2;
        return found;

    }

    EdgeGrid::Grid grid;
    ExPolygons ex_polys;
    coord_t edge_width;
};

#ifdef DEBUG_FILES
void debug_export(Issues issues, std::string file_name) {
    Slic3r::CNumericLocalesSetter locales_setter;

    {
        FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_supports.obj").c_str()).c_str(), "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Debug files: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < issues.supports_nedded.size(); ++i) {
            fprintf(fp, "v %f %f %f  %f %f %f\n",
                    issues.supports_nedded[i](0), issues.supports_nedded[i](1), issues.supports_nedded[i](2),
                    1.0, 0.0, 0.0
                    );
        }

        fclose(fp);
    }
    {
        FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_curling.obj").c_str()).c_str(), "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Debug files: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < issues.curling_up.size(); ++i) {
            fprintf(fp, "v %f %f %f  %f %f %f\n",
                    issues.curling_up[i](0), issues.curling_up[i](1), issues.curling_up[i](2),
                    0.0, 1.0, 0.0
                    );
        }
        fclose(fp);
    }

}
#endif

EdgeGridWrapper compute_layer_edge_grid(const Layer *layer) {
    float min_region_flow_width { 1.0f };
    for (const auto *region : layer->regions()) {
        min_region_flow_width = std::min(min_region_flow_width, region->flow(FlowRole::frExternalPerimeter).width());
    }
    ExPolygons ex_polygons;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                if (perimeter->role() == ExtrusionRole::erExternalPerimeter
                        || perimeter->role() == ExtrusionRole::erOverhangPerimeter) {
                    Points perimeter_points { };
                    perimeter->collect_points(perimeter_points);
                    assert(perimeter->is_loop());
                    perimeter_points.pop_back();
                    ex_polygons.push_back(ExPolygon { perimeter_points });
                } // ex_perimeter
            } // perimeter
        } // ex_entity
    }

    return EdgeGridWrapper(scale_(min_region_flow_width), ex_polygons);
}

coordf_t get_max_allowed_distance(ExtrusionRole role, coord_t flow_width, const Params &params) { // <= distance / flow_width (can be larger for perimeter, if not external perimeter first)
    if (!params.external_perimeter_first
            && (role == ExtrusionRole::erExternalPerimeter || role == ExtrusionRole::erOverhangPerimeter)) {
        return params.max_ex_perim_unsupported_distance_factor * flow_width;
    } else {
        return params.max_unsupported_distance_factor * flow_width;
    }
}

Issues check_extrusion_entity_stability(const ExtrusionEntity *entity,
        float slice_z,
        coordf_t flow_width,
        const EdgeGridWrapper &supported_grid,
        const Params &params) {

    Issues issues { };
    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            issues.add(check_extrusion_entity_stability(e, slice_z, flow_width, supported_grid, params));
        }
    } else { //single extrusion path, with possible varying parameters
        std::stack<Point> points { };
        for (const auto &p : entity->as_polyline().points) {
            points.push(p);
        }

        float unsupported_distance = params.bridge_distance + 1.0f;
        float curvature = 0;
        float max_curvature = 0;
        Vec2f tmp = unscale(points.top()).cast<float>();
        Vec3f prev_fpoint = Vec3f(tmp.x(), tmp.y(), slice_z);

        const coordf_t max_allowed_dist_from_prev_layer = get_max_allowed_distance(entity->role(), flow_width,
                params);

        while (!points.empty()) {
            Point point = points.top();
            points.pop();
            Vec2f tmp = unscale(point).cast<float>();
            Vec3f fpoint = Vec3f(tmp.x(), tmp.y(), slice_z);

            coordf_t dist_from_prev_layer { 0 };
            if (!supported_grid.signed_distance(point, flow_width, dist_from_prev_layer)) {
                issues.supports_nedded.push_back(fpoint);
                continue;
            }

            if (dist_from_prev_layer > max_allowed_dist_from_prev_layer) { //unsupported
                unsupported_distance += (fpoint - prev_fpoint).norm();

                if (!points.empty()) {
                    const Vec2f v1 = (fpoint - prev_fpoint).head<2>();
                    const Vec2f v2 = unscale(points.top()).cast<float>() - fpoint.head<2>();
                    float dot = v1(0) * v2(0) + v1(1) * v2(1);
                    float cross = v1(0) * v2(1) - v1(1) * v2(0);
                    float angle = float(atan2(float(cross), float(dot)));

                    curvature += angle;
                    max_curvature = std::max(abs(curvature), max_curvature);
                }

                if (unsupported_distance
                        > params.bridge_distance
                                / (1.0f + (max_curvature * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    issues.supports_nedded.push_back(fpoint);

                    std::cout << "SUPP: " << "udis: " << unsupported_distance << " curv: " << curvature << " max curv: "
                            << max_curvature << std::endl;
                    std::cout << "max dist from layer: " << max_allowed_dist_from_prev_layer << "  measured dist: "
                            << dist_from_prev_layer << "  FW: " << flow_width << std::endl;

                    unsupported_distance = 0;
                    curvature = 0;
                    max_curvature = 0;
                }
            } else {
                unsupported_distance = 0;
                curvature = 0;
                max_curvature = 0;
            }

            if (max_curvature / (PI * unsupported_distance) > params.limit_curvature) {
                issues.curling_up.push_back(fpoint);
            }

            prev_fpoint = fpoint;

            if (!points.empty()) {
                Vec2f next = unscale(points.top()).cast<float>();
                Vec2f reverse_v = fpoint.head<2>() - next;
                float dist_to_next = reverse_v.norm();
                reverse_v.normalize();
                int new_points_count = dist_to_next / params.bridge_distance;
                float step_size = dist_to_next / (new_points_count + 1);
                for (int i = 1; i <= new_points_count; ++i) {
                    points.push(Point::new_scale(Vec2f(next + reverse_v * (i * step_size))));
                }
            }

        }
    }
    return issues;
}

//TODO needs revision
coordf_t get_flow_width(const LayerRegion *region, ExtrusionRole role) {
    switch (role) {
        case ExtrusionRole::erBridgeInfill:
            return region->flow(FlowRole::frExternalPerimeter).scaled_width();
        case ExtrusionRole::erExternalPerimeter:
            return region->flow(FlowRole::frExternalPerimeter).scaled_width();
        case ExtrusionRole::erGapFill:
            return region->flow(FlowRole::frInfill).scaled_width();
        case ExtrusionRole::erPerimeter:
            return region->flow(FlowRole::frPerimeter).scaled_width();
        case ExtrusionRole::erSolidInfill:
            return region->flow(FlowRole::frSolidInfill).scaled_width();
        default:
            return region->flow(FlowRole::frPerimeter).scaled_width();
    }
}

Issues check_layer_stability(const PrintObject *po, size_t layer_idx, bool full_check, const Params &params) {
    std::cout << "Checking: " << layer_idx << std::endl;
    if (layer_idx == 0) {
        // first layer is usually ok
        return {};
    }
    const Layer *layer = po->get_layer(layer_idx);
    EdgeGridWrapper supported_grid = compute_layer_edge_grid(layer->lower_layer);

    Issues issues { };
    if (full_check) {
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    issues.add(check_extrusion_entity_stability(perimeter,
                            layer->slice_z, get_flow_width(layer_region, perimeter->role()),
                            supported_grid, params));
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill || fill->role() == ExtrusionRole::erBridgeInfill) {
                        issues.add(check_extrusion_entity_stability(fill,
                                layer->slice_z, get_flow_width(layer_region, fill->role()),
                                supported_grid, params));
                    }
                } // fill
            } // ex_entity
        } // region
    } else { //check only external perimeters
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter
                            || perimeter->role() == ExtrusionRole::erOverhangPerimeter) {
                        issues.add(check_extrusion_entity_stability(perimeter,
                                layer->slice_z, get_flow_width(layer_region, perimeter->role()),
                                supported_grid, params));
                    }; // ex_perimeter
                } // perimeter
            } // ex_entity
        } //region
    }

    return issues;
}

} //Impl End

std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
    using namespace Impl;
    std::vector<LayerDescriptor> descriptors(po->layer_count());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, po->layer_count()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    const Layer *layer = po->get_layer(layer_idx);

                    LayerDescriptor &descriptor = descriptors[layer_idx];
                    size_t point_count { 0 };

                    for (const LayerRegion *layer_region : layer->regions()) {
                        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
        if (perimeter->role() == ExtrusionRole::erExternalPerimeter || perimeter->role() == ExtrusionRole::erOverhangPerimeter) {
            assert(perimeter->is_loop());
            descriptor.segments_count++;
            const ExtrusionLoop *loop = static_cast<const ExtrusionLoop*>(perimeter);
            for (const ExtrusionPath& path : loop->paths) {
                Vec2f prev_pos = unscale(path.polyline.last_point()).cast<float>();
                for (size_t p_idx = 0; p_idx < path.polyline.points.size(); ++p_idx) {
                    point_count++;
                    Vec2f point_pos = unscale(path.polyline.points[p_idx]).cast<float>();
                    descriptor.centroid += point_pos;
                    descriptor.perimeter_length += (point_pos - prev_pos).norm();
                    prev_pos = point_pos;
                } //point
            } //path
        } // ex_perimeter
    } // perimeter
} // ex_entity
} // region

                    descriptor.centroid /= float(point_count);

                } // layer
            } // thread
            );

    std::vector<size_t> suspicious_layers_indices { };

    for (size_t desc_idx = 1; desc_idx < descriptors.size(); ++desc_idx) {
        const LayerDescriptor &prev = descriptors[desc_idx - 1];
        const LayerDescriptor &descriptor = descriptors[desc_idx];
        if (descriptor.segments_count - prev.segments_count != 0
                ||
                std::abs(descriptor.perimeter_length - prev.perimeter_length)
                        > params.perimeter_length_diff_tolerance ||
                (descriptor.centroid - prev.centroid).norm() > params.centroid_offset_tolerance
                        ) {
            suspicious_layers_indices.push_back(desc_idx);
        }
#ifdef DEBUG_FILES
        std::cout << "SIS layer idx: " << desc_idx << " reg count: " << descriptor.segments_count << " len: "
                << descriptor.perimeter_length <<
                " centroid: " << descriptor.centroid.x() << " | " << descriptor.centroid.y() << std::endl;
        std::cout << "SIS diff: " << desc_idx << " reg count: "
                << (int(descriptor.segments_count) - int(prev.segments_count)) <<
                " len: " << (descriptor.perimeter_length - prev.perimeter_length) <<
                " centroid: " << (descriptor.centroid - prev.centroid).norm() << std::endl;
#endif
    }

    std::vector<bool> layer_needs_supports(suspicious_layers_indices.size(), false);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, suspicious_layers_indices.size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t suspicious_index = r.begin(); suspicious_index < r.end(); ++suspicious_index) {
                    auto layer_issues = check_layer_stability(po, suspicious_layers_indices[suspicious_index],
                            false, params);
                    if (!layer_issues.supports_nedded.empty()) {
                        layer_needs_supports[suspicious_index] = true;
                    }
                }
            });

    std::vector<size_t> problematic_layers;

    for (size_t index = 0; index < suspicious_layers_indices.size(); ++index) {
        if (layer_needs_supports[index]) {
            problematic_layers.push_back(suspicious_layers_indices[index]);
        }
    }
    return problematic_layers;
}

Issues full_search(const PrintObject *po, const Params &params) {
    using namespace Impl;
    size_t layer_count = po->layer_count();
    Issues found_issues = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, layer_count), Issues { },
            [&](tbb::blocked_range<size_t> r, const Issues &init) {
                Issues issues = init;
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    auto layer_issues = check_layer_stability(po, layer_idx, true, params);
                    if (!layer_issues.empty()) {
                        issues.add(layer_issues);
                    }
                }
                return issues;
            },
            [](Issues left, const Issues &right) {
                left.add(right);
                return left;
            }
    );

#ifdef DEBUG_FILES
    Impl::debug_export(found_issues, "issues");
#endif

    return found_issues;
}

}
}
