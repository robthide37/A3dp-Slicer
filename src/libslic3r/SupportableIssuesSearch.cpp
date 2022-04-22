#include "SupportableIssuesSearch.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <stack>

#include "libslic3r/Layer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "PolygonPointTest.hpp"

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

struct Cell {
    float weight;
    char last_extrusion_id;
};

struct WeightDistributionMatrix {
    // Lets make Z coord half the size of X (and Y).
    // This corresponds to angle of ~26 degrees between center of one cell and other one up and sideways
    // which is approximately a limiting printable angle.

    WeightDistributionMatrix(const PrintObject *po, size_t layer_idx_begin, size_t layer_idx_end) {
        Vec3crd object_origin = scaled(po->trafo_centered() * Vec3d::Zero());
        Vec3crd min = object_origin - po->size() / 2 - Vec3crd::Ones();
        Vec3crd max = object_origin + po->size() / 2 + Vec3crd::Ones();

        cell_size = Vec3crd { int(cell_height * 2), int(cell_height * 2), int(cell_height) };

        global_origin = min;
        global_size = max - min;
        global_cell_count = global_size.cwiseQuotient(cell_size);

        coord_t local_min_z = scale_(po->layers()[layer_idx_begin]->slice_z);
        coord_t local_max_z = scale_(po->layers()[layer_idx_end]->slice_z);
        coord_t local_min_z_index = local_min_z / cell_size.z();
        coord_t local_max_z_index = local_max_z / cell_size.z();

        local_z_index_offset = local_min_z_index;
        local_z_cell_count = local_max_z_index - local_min_z_index + 1;

        cells.resize(local_z_cell_count * global_cell_count.y() * global_cell_count.x());
    }

    Vec3i to_global_cell_coords(const Point &p, float slice_z) const {
        Vec3crd position = Vec3crd { p.x(), p.y(), coord_t(scale_(slice_z)) };
        Vec3i cell_coords = position.cwiseQuotient(cell_size);
        return cell_coords;
    }

    Vec3i to_local_cell_coords(const Point &p, float slice_z) const {
        Vec3i cell_coords = to_global_cell_coords(p, slice_z);
        Vec3i local_cell_coords = cell_coords - local_z_index_offset * Vec3i::UnitZ();
        return local_cell_coords;
    }

    size_t to_cell_index(const Vec3i &local_cell_coords) {
        assert(local_cell_coords.x() >= 0);
        assert(local_cell_coords.x() < global_cell_count.x());
        assert(local_cell_coords.y() >= 0);
        assert(local_cell_coords.y() < global_cell_count.y());
        assert(local_cell_coords.z() >= 0);
        assert(local_cell_coords.z() < local_z_cell_count);
        return local_cell_coords.z() * global_cell_count.x() * global_cell_count.y()
                + local_cell_coords.y() * global_cell_count.x() +
                local_cell_coords.x();
    }

    Vec3crd cell_center(const Vec3i &global_cell_coords) {
        return global_origin + global_cell_coords.cwiseProduct(cell_size);
    }

    Cell& access_cell(const Point &p, float slice_z) {
        return cells[to_cell_index(to_local_cell_coords(p, slice_z))];
    }

    Cell& access_cell(const Vec3i& local_cell_coords) {
        return cells[to_cell_index(local_cell_coords)];
    }


#ifdef DEBUG_FILES
    void debug_export(std::string file_name) {
        Slic3r::CNumericLocalesSetter locales_setter;
        {
            FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_matrix.obj").c_str()).c_str(), "w");
            if (fp == nullptr) {
                BOOST_LOG_TRIVIAL(error)
                << "Debug files: Couldn't open " << file_name << " for writing";
                return;
            }

            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    for (int z = 0; z < local_z_cell_count; ++z) {
                        Vec3f center = unscale(cell_center(Vec3i(x, y, z + local_z_index_offset))).cast<float>();
                        Cell &cell = access_cell(Vec3i(x, y, z));
                        fprintf(fp, "v %f %f %f  %f %f %f\n",
                                center(0), center(1),
                                center(2),
                                cell.weight, 0.0, 0.0
                                );
                    }
                }
            }

            fclose(fp);
        }
    }
#endif

    static constexpr float cell_height = scale_(0.15f);

    Vec3crd cell_size;

    Vec3crd global_origin;
    Vec3crd global_size;
    Vec3i global_cell_count;

    int local_z_index_offset;
    int local_z_cell_count;
    std::vector<Cell> cells;

};

namespace Impl {

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
    std::vector<Points> lines;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            lines.push_back(Points { });
            ex_entity->collect_points(lines.back());
        } // ex_entity

        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            lines.push_back(Points { });
            ex_entity->collect_points(lines.back());
        } // ex_entity
    }

    return EdgeGridWrapper(scale_(min_region_flow_width), lines);
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

coordf_t get_max_allowed_distance(ExtrusionRole role, coord_t flow_width, bool external_perimeters_first,
        const Params &params) { // <= distance / flow_width (can be larger for perimeter, if not external perimeter first)
    if ((role == ExtrusionRole::erExternalPerimeter || role == ExtrusionRole::erOverhangPerimeter)
            && (external_perimeters_first)
            ) {
        return params.max_first_ex_perim_unsupported_distance_factor * flow_width;
    } else {
        return params.max_unsupported_distance_factor * flow_width;
    }
}

struct SegmentAccumulator {
    float distance = 0; //accumulated distance
    float curvature = 0; //accumulated signed ccw angles
    float max_curvature = 0; //max absolute accumulated value

    void add_distance(float dist) {
        distance += dist;
    }

    void add_angle(float ccw_angle) {
        curvature += ccw_angle;
        max_curvature = std::max(max_curvature, std::abs(curvature));
    }

    void reset() {
        distance = 0;
        curvature = 0;
        max_curvature = 0;
    }

};

Issues check_extrusion_entity_stability(const ExtrusionEntity *entity,
        float slice_z,
        const LayerRegion *layer_region,
        const EdgeGridWrapper &supported_grid,
        const Params &params) {

    Issues issues { };
    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            issues.add(check_extrusion_entity_stability(e, slice_z, layer_region, supported_grid, params));
        }
    } else { //single extrusion path, with possible varying parameters
        //prepare stack of points on the extrusion path. If there are long segments, additional points might be pushed onto the stack during the algorithm.
        std::stack<Point> points { };
        for (const auto &p : entity->as_polyline().points) {
            points.push(p);
        }

        SegmentAccumulator supports_acc { };
        supports_acc.add_distance(params.bridge_distance + 1.0f); // initialize unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter start and short loops into air.
        SegmentAccumulator curling_acc { };

        const auto to_vec3f = [slice_z](const Point &point) {
            Vec2f tmp = unscale(point).cast<float>();
            return Vec3f(tmp.x(), tmp.y(), slice_z);
        };

        Vec3f prev_fpoint = to_vec3f(points.top()); // prev point of the path. Initialize with first point.
        coordf_t flow_width = get_flow_width(layer_region, entity->role());
        bool external_perimters_first = layer_region->region().config().external_perimeters_first;
        const coordf_t max_allowed_dist_from_prev_layer = get_max_allowed_distance(entity->role(), flow_width,
                external_perimters_first, params);

        while (!points.empty()) {
            Point point = points.top();
            points.pop();
            Vec2f tmp = unscale(point).cast<float>();
            Vec3f fpoint = Vec3f(tmp.x(), tmp.y(), slice_z);
            float edge_len = (fpoint - prev_fpoint).norm();

            coordf_t dist_from_prev_layer { 0 };
            if (!supported_grid.signed_distance(point, flow_width, dist_from_prev_layer)) { // dist from prev layer not found, assume empty layer
                issues.supports_nedded.push_back(fpoint);
                supports_acc.reset();
            }

            float angle = 0;
            if (!points.empty()) {
                const Vec2f v1 = (fpoint - prev_fpoint).head<2>();
                const Vec2f v2 = unscale(points.top()).cast<float>() - fpoint.head<2>();
                float dot = v1(0) * v2(0) + v1(1) * v2(1);
                float cross = v1(0) * v2(1) - v1(1) * v2(0);
                angle = float(atan2(float(cross), float(dot))); // ccw angle, TODO replace with angle func, once it gets into master
            }

            supports_acc.add_angle(angle);
            curling_acc.add_angle(angle);

            if (dist_from_prev_layer > max_allowed_dist_from_prev_layer) { //extrusion point is unsupported
                supports_acc.add_distance(edge_len); // for algorithm simplicity, expect that the whole line between prev and current point is unsupported

                if (supports_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                > params.bridge_distance
                        / (1.0f
                                + (supports_acc.max_curvature
                                        * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    issues.supports_nedded.push_back(fpoint);
                    supports_acc.reset();
                }
            } else {
                supports_acc.reset();
            }

            // Estimation of short curvy segments which are not supported -> problems with curling
            if (dist_from_prev_layer > 0.0f) { //extrusion point is unsupported or poorly supported
                curling_acc.add_distance(edge_len);
                if (curling_acc.max_curvature / (PI * curling_acc.distance) > params.limit_curvature) {
                    issues.curling_up.push_back(fpoint);
                    curling_acc.reset();
                }
            } else {
                curling_acc.reset();
            }

            prev_fpoint = fpoint;

            if (!points.empty()) { //oversampling if necessary
                Vec2f next = unscale(points.top()).cast<float>();
                Vec2f reverse_v = fpoint.head<2>() - next; // vector from next to current
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

Issues check_layer_stability(const PrintObject *po, size_t layer_idx, bool full_check, const Params &params) {
    std::cout << "Checking: " << layer_idx << std::endl;
    if (layer_idx == 0) {
        // first layer is usually ok
        return {};
    }
    const Layer *layer = po->get_layer(layer_idx);
    //Prepare edge grid of previous layer, will be used to check if the extrusion path is supported
    EdgeGridWrapper supported_grid = compute_layer_edge_grid(layer->lower_layer);

    Issues issues { };
    if (full_check) { // If full checkm check stability of perimeters, gap fills, and bridges.
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    issues.add(check_extrusion_entity_stability(perimeter,
                            layer->slice_z, layer_region,
                            supported_grid, params));
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill || fill->role() == ExtrusionRole::erBridgeInfill) {
                        issues.add(check_extrusion_entity_stability(fill,
                                layer->slice_z, layer_region,
                                supported_grid, params));
                    }
                } // fill
            } // ex_entity
        } // region

    } else { // If NOT full check, check only external perimeters
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter
                            || perimeter->role() == ExtrusionRole::erOverhangPerimeter) {
                        issues.add(check_extrusion_entity_stability(perimeter,
                                layer->slice_z, layer_region,
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

    size_t layer_count = po->layer_count();
    std::vector<bool> layer_needs_supports(layer_count, false);
    tbb::parallel_for(tbb::blocked_range<size_t>(1, layer_count),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    auto layer_issues = check_layer_stability(po, layer_idx,
                            false, params);
                    if (!layer_issues.supports_nedded.empty()) {
                        layer_needs_supports[layer_idx] = true;
                    }
                }
            });

    std::vector<size_t> problematic_layers;
    for (size_t index = 0; index < layer_needs_supports.size(); ++index) {
        if (layer_needs_supports[index]) {
            problematic_layers.push_back(index);
        }
    }
    return problematic_layers;
}

Issues full_search(const PrintObject *po, const Params &params) {
    using namespace Impl;

    WeightDistributionMatrix matrix { po, 0, po->layers().size() };
    matrix.debug_export("matrix");

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
