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
    int island_id;
};

struct WeightDistributionMatrix {
    // Lets make Z coord half the size of X (and Y).
    // This corresponds to angle of ~26 degrees between center of one cell and other one up and sideways
    // which is approximately a limiting printable angle.

    WeightDistributionMatrix() = default;

    void init(const PrintObject *po, size_t layer_idx_begin, size_t layer_idx_end) {
        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec3crd min = Vec3crd(-size_half.x(), -size_half.y(), 0);
        Vec3crd max = Vec3crd(size_half.x(), size_half.y(), po->height());

        cell_size = Vec3crd { int(cell_height * 2), int(cell_height * 2), int(cell_height) };
        assert(cell_size.x() == cell_size.y());

        global_origin = min;
        global_size = max - min;
        global_cell_count = global_size.cwiseQuotient(cell_size) + Vec3i::Ones();

        coord_t local_min_z = scale_(po->layers()[layer_idx_begin]->print_z);
        coord_t local_max_z = scale_(po->layers()[layer_idx_end > 0 ? layer_idx_end - 1 : 0]->print_z);
        int local_min_z_index = local_min_z / cell_size.z();
        int local_max_z_index = local_max_z / cell_size.z() + 1;

        local_z_index_offset = local_min_z_index;
        local_z_cell_count = local_max_z_index + 1 - local_min_z_index;

        cells.resize(local_z_cell_count * global_cell_count.y() * global_cell_count.x());
    }

    Vec3i to_global_cell_coords(const Vec3i &local_cell_coords) const {
        return local_cell_coords + local_z_index_offset * Vec3i::UnitZ();
    }

    Vec3i to_local_cell_coords(const Vec3i &global_cell_coords) const {
        return global_cell_coords - local_z_index_offset * Vec3i::UnitZ();
    }

    Vec3i to_global_cell_coords(const Point &p, float print_z) const {
        Vec3i position = Vec3crd { p.x(), p.y(), int(scale_(print_z)) };
        Vec3i cell_coords = (position - this->global_origin).cwiseQuotient(this->cell_size);
        return cell_coords;
    }

    Vec3i to_local_cell_coords(const Point &p, float print_z) const {
        Vec3i cell_coords = this->to_global_cell_coords(p, print_z);
        return this->to_local_cell_coords(cell_coords);
    }

    size_t to_cell_index(const Vec3i &local_cell_coords) const {
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

    Vec3crd get_cell_center(const Vec3i &global_cell_coords) const {
        return global_origin + global_cell_coords.cwiseProduct(this->cell_size)
                + this->cell_size.cwiseQuotient(Vec3crd(2, 2, 2));
    }

    Cell& access_cell(const Point &p, float print_z) {
        return cells[this->to_cell_index(to_local_cell_coords(p, print_z))];
    }

    Cell& access_cell(const Vec3i &local_cell_coords) {
        return cells[this->to_cell_index(local_cell_coords)];
    }

    const Cell& access_cell(const Vec3i &local_cell_coords) const {
        return cells[this->to_cell_index(local_cell_coords)];
    }

    void ditribute_edge_weight(const Point &p1, const Point &p2, float print_z, coordf_t width) {
        Vec2d dir = (p2 - p1).cast<double>();
        double length = dir.norm();
        if (length < 0.01) {
            return;
        }
        dir /= length;
        double step_size = this->cell_size.x() / 2.0;

        double distributed_length = 0;
        while (distributed_length < length) {
            double next_len = std::min(length, distributed_length + step_size);
            double current_dist_payload = next_len - distributed_length;

            Point location = p1 + ((next_len / length) * dir).cast<coord_t>();
            double payload = current_dist_payload * width;

            Vec3i local_index = this->to_local_cell_coords(location, print_z);

            if (this->to_cell_index(local_index) >= this->cells.size() || this->to_cell_index(local_index) < 0) {
                std::cout << "loc: " << local_index.x() << "  " << local_index.y() << "  " << local_index.z()
                        << "   globals: " << this->global_cell_count.x() << "  "
                        << this->global_cell_count.y() << "   " << this->local_z_cell_count <<
                        "+" << this->local_z_cell_count << std::endl;
                return;
            }
            this->access_cell(location, print_z).weight += payload;

            distributed_length = next_len;
        }
    }

    void merge(const WeightDistributionMatrix &other) {
        int z_start = std::max(local_z_index_offset, other.local_z_index_offset);
        int z_end = std::min(local_z_index_offset + local_z_cell_count,
                other.local_z_index_offset + other.local_z_cell_count);

        for (int x = 0; x < global_cell_count.x(); ++x) {
            for (int y = 0; y < global_cell_count.y(); ++y) {
                for (int z = z_start; z < z_end; ++z) {
                    Vec3i global_coords { x, y, z };
                    Vec3i local_coords = this->to_local_cell_coords(global_coords);
                    Vec3i other_local_coords = other.to_local_cell_coords(global_coords);
                    this->access_cell(local_coords).weight += other.access_cell(other_local_coords).weight;
                }
            }
        }
    }

    void distribute_top_down() {
        const auto validate_xy_coords = [&](const Vec2i &local_coords) {
            return local_coords.x() >= 0 && local_coords.y() >= 0 &&
                    local_coords.x() < this->global_cell_count.x() && local_coords.y() < this->global_cell_count.y();
        };

        Vec2i valid_coords[9];

        for (int x = 0; x < global_cell_count.x(); ++x) {
            for (int y = 0; y < global_cell_count.y(); ++y) {
                for (int z = local_z_cell_count - 1; z > local_z_index_offset; --z) {
                    Cell &current = this->access_cell(Vec3i(x, y, z));
                    size_t valid_coords_count = 0;
                    if (current.weight > 0) {
                        for (int y_offset = -1; y_offset <= 1; ++y_offset) {
                            for (int x_offset = -1; x_offset <= 1; ++x_offset) {
                                Vec2i xy_coords { x + x_offset, y + y_offset };
                                if (validate_xy_coords(xy_coords)
                                        &&
                                        this->access_cell(Vec3i(xy_coords.x(), xy_coords.y(), z - 1)).weight != 0) {
                                    valid_coords[valid_coords_count] = xy_coords;
                                    valid_coords_count++;
                                }
                            }
                        }

                        float distribution = current.weight / valid_coords_count;
                        for (size_t index = 0; index < valid_coords_count; ++index) {
                            this->access_cell(Vec3i(valid_coords[index].x(), valid_coords[index].y(), z - 1)).weight +=
                                    distribution;
                        }

                        if (valid_coords_count > 0) {
                            current.weight = 0;
                        }
                    }
                }

            }
        }
    }

#ifdef DEBUG_FILES
    void debug_export(std::string file_name) const {
        Slic3r::CNumericLocalesSetter locales_setter;
        {
            FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_matrix.obj").c_str()).c_str(), "w");
            if (fp == nullptr) {
                BOOST_LOG_TRIVIAL(error)
                << "Debug files: Couldn't open " << file_name << " for writing";
                return;
            }

            float max_weight = 0;
            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    for (int z = 0; z < local_z_cell_count; ++z) {
                        const Cell &cell = access_cell(Vec3i(x, y, z));
                        max_weight = std::max(max_weight, cell.weight);
                    }
                }
            }

            max_weight *= 0.8;

            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    for (int z = 0; z < local_z_cell_count; ++z) {
                        Vec3f center = unscale(get_cell_center(to_global_cell_coords(Vec3i { x, y, z }))).cast<float>();
                        const Cell &cell = access_cell(Vec3i(x, y, z));
                        if (cell.weight != 0) {
                            fprintf(fp, "v %f %f %f  %f %f %f\n",
                                    center(0), center(1),
                                    center(2),
                                    cell.weight / max_weight, 0.0, 0.0
                                    );
                        }
                    }
                }
            }

            fclose(fp);
        }
    }
#endif

    static constexpr float cell_height = scale_(0.3f);

    Vec3crd cell_size { };

    Vec3crd global_origin { };
    Vec3crd global_size { };
    Vec3i global_cell_count { };

    int local_z_index_offset { };
    int local_z_cell_count { };
    std::vector<Cell> cells { };

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

coordf_t get_max_allowed_distance(ExtrusionRole role, coordf_t flow_width, bool external_perimeters_first,
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
        float print_z,
        const LayerRegion *layer_region,
        const EdgeGridWrapper &supported_grid,
        WeightDistributionMatrix &weight_matrix,
        const Params &params) {

    Issues issues { };
    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            issues.add(
                    check_extrusion_entity_stability(e, print_z, layer_region, supported_grid, weight_matrix, params));
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

        const auto to_vec3f = [print_z](const Point &point) {
            Vec2f tmp = unscale(point).cast<float>();
            return Vec3f(tmp.x(), tmp.y(), print_z);
        };

        Point prev_point = points.top(); // prev point of the path. Initialize with first point.
        Vec3f prev_fpoint = to_vec3f(prev_point);
        coordf_t flow_width = get_flow_width(layer_region, entity->role());
        bool external_perimters_first = layer_region->region().config().external_perimeters_first;
        const coordf_t max_allowed_dist_from_prev_layer = get_max_allowed_distance(entity->role(), flow_width,
                external_perimters_first, params);

        while (!points.empty()) {
            Point point = points.top();
            points.pop();
            Vec2f tmp = unscale(point).cast<float>();
            Vec3f fpoint = Vec3f(tmp.x(), tmp.y(), print_z);
            float edge_len = (fpoint - prev_fpoint).norm();

            weight_matrix.ditribute_edge_weight(prev_point, point, print_z, flow_width);

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

            prev_point = point;
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

Issues check_layer_stability(const PrintObject *po, size_t layer_idx, bool full_check,
        WeightDistributionMatrix &weight_matrix, const Params &params) {
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
                            layer->print_z, layer_region,
                            supported_grid, weight_matrix, params));
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill || fill->role() == ExtrusionRole::erBridgeInfill) {
                        issues.add(check_extrusion_entity_stability(fill,
                                layer->print_z, layer_region,
                                supported_grid, weight_matrix, params));
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
                                layer->print_z, layer_region,
                                supported_grid, weight_matrix, params));
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

    WeightDistributionMatrix matrix { };
    matrix.init(po, 0, po->layers().size());
    std::mutex matrix_mutex;

    size_t layer_count = po->layer_count();
    std::vector<bool> layer_needs_supports(layer_count, false);
    tbb::parallel_for(tbb::blocked_range<size_t>(1, layer_count),
            [&](tbb::blocked_range<size_t> r) {
                WeightDistributionMatrix weight_matrix { };
                weight_matrix.init(po, r.begin(), r.end());

                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    auto layer_issues = check_layer_stability(po, layer_idx,
                            false, weight_matrix, params);
                    if (!layer_issues.supports_nedded.empty()) {
                        layer_needs_supports[layer_idx] = true;
                    }
                }

                matrix_mutex.lock();
                matrix.merge(weight_matrix);
                matrix_mutex.unlock();
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

    WeightDistributionMatrix matrix { };
    matrix.init(po, 0, po->layers().size());
    std::mutex matrix_mutex;

    size_t layer_count = po->layer_count();
    Issues found_issues = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, layer_count), Issues { },
            [&](tbb::blocked_range<size_t> r, const Issues &init) {
                WeightDistributionMatrix weight_matrix { };
                weight_matrix.init(po, r.begin(), r.end());
                Issues issues = init;
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    auto layer_issues = check_layer_stability(po, layer_idx, true, weight_matrix, params);
                    if (!layer_issues.empty()) {
                        issues.add(layer_issues);
                    }
                }

                matrix_mutex.lock();
                matrix.merge(weight_matrix);
                matrix_mutex.unlock();

                return issues;
            },
            [](Issues left, const Issues &right) {
                left.add(right);
                return left;
            }
    );

    matrix.distribute_top_down();

    matrix.debug_export("weight");

#ifdef DEBUG_FILES
    Impl::debug_export(found_issues, "issues");
#endif

    return found_issues;
}

}
}
