#include "SupportableIssuesSearch.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <unordered_set>
#include <stack>

#include "AABBTreeLines.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "Geometry/ConvexHull.hpp"

#define DEBUG_FILES

#ifdef DEBUG_FILES
#include <boost/nowide/cstdio.hpp>
#include "libslic3r/Color.hpp"
#endif

namespace Slic3r {

static const size_t NULL_ACC_ID = std::numeric_limits<size_t>::max();

class ExtrusionLine
{
public:
    ExtrusionLine() :
            a(Vec2f::Zero()), b(Vec2f::Zero()), len(0.0f) {
    }
    ExtrusionLine(const Vec2f &_a, const Vec2f &_b) :
            a(_a), b(_b), len((_a - _b).norm()) {
    }

    float length() {
        return (a - b).norm();
    }

    Vec2f a;
    Vec2f b;
    float len;

    size_t stability_accumulator_id = NULL_ACC_ID;

    static const constexpr int Dim = 2;
    using Scalar = Vec2f::Scalar;
};

auto get_a(ExtrusionLine &&l) {
    return l.a;
}
auto get_b(ExtrusionLine &&l) {
    return l.b;
}

namespace SupportableIssues {

void Issues::add(const Issues &layer_issues) {
    supports_nedded.insert(supports_nedded.end(), layer_issues.supports_nedded.begin(),
            layer_issues.supports_nedded.end());
    curling_up.insert(curling_up.end(), layer_issues.curling_up.begin(), layer_issues.curling_up.end());
}

bool Issues::empty() const {
    return supports_nedded.empty() && curling_up.empty();
}

SupportPoint::SupportPoint(const Vec3f &position, float weight) :
        position(position), weight(weight) {
}

CurledFilament::CurledFilament(const Vec3f &position, float estimated_height) :
        position(position), estimated_height(estimated_height) {
}

CurledFilament::CurledFilament(const Vec3f &position) :
        position(position), estimated_height(0.0f) {
}

class LayerLinesDistancer {
private:
    std::vector<ExtrusionLine> lines;
    AABBTreeIndirect::Tree<2, float> tree;

public:
    explicit LayerLinesDistancer(std::vector<ExtrusionLine> &&lines) :
            lines(lines) {
        tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(lines);
    }

    // negative sign means inside
    float signed_distance_from_lines(const Vec2f &point, size_t &nearest_line_index_out,
            Vec2f &nearest_point_out) const {
        auto distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, point, nearest_line_index_out,
                nearest_point_out);
        if (distance < 0)
            return std::numeric_limits<float>::infinity();

        distance = sqrt(distance);
        const ExtrusionLine &line = lines[nearest_line_index_out];
        Vec2f v1 = line.b - line.a;
        Vec2f v2 = point - line.a;
        if ((v1.x() * v2.y()) - (v1.y() * v2.x()) > 0.0) {
            distance *= -1;
        }
        return distance;
    }

    const ExtrusionLine& get_line(size_t line_idx) const {
        return lines[line_idx];
    }

    const std::vector<ExtrusionLine>& get_lines() const {
        return lines;
    }
};

class StabilityAccumulator {
private:
    Polygon base_convex_hull { };
    Points support_points { };
    Vec3f centroid_accumulator = Vec3f::Zero();
    float accumulated_volume { };
    float base_area { };
    float base_height { };

public:
    explicit StabilityAccumulator(float base_height) :
            base_height(base_height) {
    }

    void add_base_extrusion(const ExtrusionLine &line, float width, float print_z, float mm3_per_mm) {
        base_area += line.len * width;
        support_points.push_back(Point::new_scale(line.a));
        support_points.push_back(Point::new_scale(line.b));
        base_convex_hull.clear();
        add_extrusion(line, print_z, mm3_per_mm);
    }

    void add_support_point(const Point &position, float area) {
        support_points.push_back(position);
        base_convex_hull.clear();
        base_area += area;
    }

    void add_extrusion(const ExtrusionLine &line, float print_z, float mm3_per_mm) {
        float volume = line.len * mm3_per_mm;
        accumulated_volume += volume;
        Vec2f center = (line.a + line.b) / 2.0f;
        centroid_accumulator += volume * Vec3f(center.x(), center.y(), print_z);
    }

    Vec3f get_centroid() const {
        return centroid_accumulator / accumulated_volume;
    }

    float get_base_area() const {
        return base_area;
    }
    float get_base_height() const {
        return base_height;
    }
    float get_accumulated_volume() const {
        return accumulated_volume;
    }

    const Polygon& segment_base_hull() {
        if (this->base_convex_hull.empty()) {
            this->base_convex_hull = Geometry::convex_hull(this->support_points);
        }
        return this->base_convex_hull;
    }

    const Points& get_support_points() {
        return support_points;
    }

    void add_from(const StabilityAccumulator &acc) {
        this->support_points.insert(this->support_points.end(), acc.support_points.begin(),
                acc.support_points.end());
        base_convex_hull.clear();
        this->centroid_accumulator += acc.centroid_accumulator;
        this->accumulated_volume += acc.accumulated_volume;
        this->base_area += acc.base_area;
    }
};

struct StabilityAccumulators {
private:
    size_t next_id = 0;
    std::unordered_map<size_t, size_t> mapping;
    std::vector<StabilityAccumulator> accumulators;

    void merge_to(size_t from_id, size_t to_id) {
        StabilityAccumulator &from_acc = this->access(from_id);
        StabilityAccumulator &to_acc = this->access(to_id);
        if (&from_acc == &to_acc) {
            return;
        }
        to_acc.add_from(from_acc);
        mapping[from_id] = mapping[to_id];
        from_acc = StabilityAccumulator { 0.0f };

    }

public:
    StabilityAccumulators() = default;

    int create_accumulator(float base_height) {
        size_t id = next_id;
        next_id++;
        mapping[id] = accumulators.size();
        accumulators.push_back(StabilityAccumulator { base_height });
        return id;
    }

    StabilityAccumulator& access(size_t id) {
        return accumulators[mapping[id]];
    }

    void merge_accumulators(size_t from_id, size_t to_id) {
        if (from_id == NULL_ACC_ID || to_id == NULL_ACC_ID) {
            return;
        }
        StabilityAccumulator &from_acc = this->access(from_id);
        StabilityAccumulator &to_acc = this->access(to_id);
        if (&from_acc == &to_acc) {
            return;
        }
        to_acc.add_from(from_acc);
        mapping[from_id] = mapping[to_id];
        from_acc = StabilityAccumulator { 0.0f };
    }

#ifdef DEBUG_FILES
    Vec3f get_accumulator_color(size_t id) {
        if (mapping.find(id) == mapping.end()) {
            std::cerr << " ERROR: uknown accumulator ID: " << id << std::endl;
            return Vec3f(1.0f, 1.0f, 1.0f);
        }

        size_t pseudornd = ((mapping[id] + 127) * 33331 + 6907) % 987;
        return value_to_rgbf(0.0f, float(987), float(pseudornd));
    }
#endif
};

float get_flow_width(const LayerRegion *region, ExtrusionRole role) {
    switch (role) {
        case ExtrusionRole::erBridgeInfill:
            return region->flow(FlowRole::frExternalPerimeter).width();
        case ExtrusionRole::erExternalPerimeter:
            return region->flow(FlowRole::frExternalPerimeter).width();
        case ExtrusionRole::erGapFill:
            return region->flow(FlowRole::frInfill).width();
        case ExtrusionRole::erPerimeter:
            return region->flow(FlowRole::frPerimeter).width();
        case ExtrusionRole::erSolidInfill:
            return region->flow(FlowRole::frSolidInfill).width();
        case ExtrusionRole::erInternalInfill:
            return region->flow(FlowRole::frInfill).width();
        case ExtrusionRole::erTopSolidInfill:
            return region->flow(FlowRole::frTopSolidInfill).width();
        default:
            return region->flow(FlowRole::frPerimeter).width();
    }
}

struct ExtrusionPropertiesAccumulator {
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

void check_extrusion_entity_stability(const ExtrusionEntity *entity,
        StabilityAccumulators &stability_accs,
        Issues &issues,
        std::vector<ExtrusionLine> &checked_lines,
        float print_z,
        const LayerRegion *layer_region,
        const LayerLinesDistancer &prev_layer_lines,
        const Params &params) {

    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            check_extrusion_entity_stability(e, stability_accs, issues, checked_lines, print_z, layer_region,
                    prev_layer_lines,
                    params);
        }
    } else { //single extrusion path, with possible varying parameters
        const auto to_vec3f = [print_z](const Point &point) {
            Vec2f tmp = unscale(point).cast<float>();
            return Vec3f(tmp.x(), tmp.y(), print_z);
        };
        Points points { };
        entity->collect_points(points);
        std::vector<ExtrusionLine> lines;
        lines.reserve(points.size() * 1.5);
        lines.emplace_back(unscaled(points[0]).cast<float>(), unscaled(points[0]).cast<float>());
        for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
            Vec2f start = unscaled(points[point_idx]).cast<float>();
            Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
            Vec2f v = next - start; // vector from next to current
            float dist_to_next = v.norm();
            v.normalize();
            int lines_count = int(std::ceil(dist_to_next / params.bridge_distance));
            float step_size = dist_to_next / lines_count;
            for (int i = 0; i < lines_count; ++i) {
                Vec2f a(start + v * (i * step_size));
                Vec2f b(start + v * ((i + 1) * step_size));
                lines.emplace_back(a, b);
            }
        }

        size_t current_stability_acc = NULL_ACC_ID;
        ExtrusionPropertiesAccumulator bridging_acc { };
        bridging_acc.add_distance(params.bridge_distance + 1.0f); // Initialise unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter start and short loops into air.
        const float flow_width = get_flow_width(layer_region, entity->role());
        const float max_allowed_dist_from_prev_layer = flow_width;
        float distance_from_last_support_point = params.min_distance_between_support_points * 2.0f;

        for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx) {
            ExtrusionLine &current_line = lines[line_idx];
            Point current = Point::new_scale(current_line.b);
            distance_from_last_support_point += current_line.len;
            float mm3_per_mm = float(entity->min_mm3_per_mm());

            float curr_angle = 0;
            if (line_idx + 1 < lines.size()) {
                const Vec2f v1 = current_line.b - current_line.a;
                const Vec2f v2 = lines[line_idx + 1].b - lines[line_idx + 1].a;
                curr_angle = angle(v1, v2);
            }
            bridging_acc.add_angle(curr_angle);

            size_t nearest_line_idx;
            Vec2f nearest_point;
            float dist_from_prev_layer = prev_layer_lines.signed_distance_from_lines(current_line.b, nearest_line_idx,
                    nearest_point);
            if (dist_from_prev_layer < max_allowed_dist_from_prev_layer) {
                const ExtrusionLine &nearest_line = prev_layer_lines.get_line(nearest_line_idx);
                size_t acc_id = nearest_line.stability_accumulator_id;
                stability_accs.merge_accumulators(std::max(acc_id, current_stability_acc),
                        std::min(acc_id, current_stability_acc));
                current_stability_acc = std::min(acc_id, current_stability_acc);
                current_line.stability_accumulator_id = current_stability_acc;
                stability_accs.access(current_stability_acc).add_extrusion(current_line, print_z, mm3_per_mm);
                bridging_acc.reset();
                // TODO curving here
            } else {
                bridging_acc.add_distance(current_line.len);
                if (current_stability_acc == NULL_ACC_ID) {
                    current_stability_acc = stability_accs.create_accumulator(print_z);
                }
                StabilityAccumulator &current_segment = stability_accs.access(current_stability_acc);
                current_line.stability_accumulator_id = current_stability_acc;
                current_segment.add_extrusion(current_line, print_z, mm3_per_mm);
                if (distance_from_last_support_point > params.min_distance_between_support_points &&
                        bridging_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                        > params.bridge_distance
                                / (1.0f + (bridging_acc.max_curvature
                                        * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    current_segment.add_support_point(current, params.extrusion_support_points_area);
                    issues.supports_nedded.emplace_back(to_vec3f(current), 1.0);
                    bridging_acc.reset();
                    distance_from_last_support_point = 0.0f;
                }
            }
        }
        checked_lines.insert(checked_lines.end(), lines.begin(), lines.end());
    }
}

void check_layer_global_stability(StabilityAccumulators &stability_accs,
        Issues &issues,
        const std::vector<ExtrusionLine> &checked_lines,
        float print_z,
        const Params &params) {
    std::unordered_map<StabilityAccumulator*, std::vector<size_t>> layer_accs_lines;
    for (size_t i = 0; i < checked_lines.size(); ++i) {
        layer_accs_lines[&stability_accs.access(checked_lines[i].stability_accumulator_id)].push_back(i);
    }

    for (auto &acc_lines : layer_accs_lines) {
        StabilityAccumulator *acc = acc_lines.first;
        Vec3f centroid = acc->get_centroid();
        Vec2f hull_centroid = unscaled(acc->segment_base_hull().centroid()).cast<float>();
        std::vector<ExtrusionLine> hull_lines;
        for (const Line &line : acc->segment_base_hull().lines()) {
            Vec2f start = unscaled(line.a).cast<float>();
            Vec2f next = unscaled(line.b).cast<float>();
            hull_lines.push_back( { start, next });
        }
        if (hull_lines.empty()) {
            if (acc->get_support_points().empty()) {
                acc->add_support_point(Point::new_scale(checked_lines[acc_lines.second[0]].a),
                        params.support_points_interface_area);
                issues.supports_nedded.emplace_back(to_3d(checked_lines[acc_lines.second[0]].a, print_z), 1.0);
            }
            hull_lines.push_back( { unscaled(acc->get_support_points()[0]).cast<float>(),
                    unscaled(acc->get_support_points()[0]).cast<float>() });
            hull_centroid = unscaled(acc->get_support_points()[0]).cast<float>();
        }

        LayerLinesDistancer hull_distancer(std::move(hull_lines));

        float sticking_force = acc->get_base_area()
                * (acc->get_base_height() == 0 ? params.base_adhesion : params.support_adhesion);
        float mass = acc->get_accumulated_volume() * params.filament_density;
        float weight = mass * params.gravity_constant;

        float distance_from_last_support_point = params.min_distance_between_support_points * 2.0f;
        for (size_t line_idx : acc_lines.second) {
            const ExtrusionLine &line = checked_lines[line_idx];
            distance_from_last_support_point += line.len;

            Vec3f extruder_pressure_direction = to_3d(Vec2f(line.b - line.a), 0.0f).normalized();
            Vec2f pivot_site_search = line.b + extruder_pressure_direction.head<2>() * 1000.0f;
            extruder_pressure_direction.z() = -0.1f;
            extruder_pressure_direction.normalize();

            size_t nearest_line_idx;
            Vec2f pivot;
            hull_distancer.signed_distance_from_lines(pivot_site_search, nearest_line_idx, pivot);

            float sticking_arm = (pivot - hull_centroid).norm();
            float sticking_torque = sticking_arm * sticking_force;

            std::cout << "sticking_arm: " << sticking_arm << std::endl;
            std::cout << "sticking_torque: " << sticking_torque << std::endl;

            float weight_arm = (pivot - centroid.head<2>()).norm();
            float weight_torque = weight_arm * weight;

            std::cout << "weight_arm: " << sticking_arm << std::endl;
            std::cout << "weight_torque: " << weight_torque << std::endl;

            float bed_movement_arm = centroid.z() - acc->get_base_height();
            float bed_movement_force = params.max_acceleration * mass;
            float bed_movement_torque = bed_movement_force * bed_movement_arm;

            std::cout << "bed_movement_arm: " << bed_movement_arm << std::endl;
            std::cout << "bed_movement_torque: " << bed_movement_torque << std::endl;

            float conflict_torque_arm = (to_3d(Vec2f(pivot - line.b), print_z).cross(
                    extruder_pressure_direction)).norm();
            float extruder_conflict_torque = params.tolerable_extruder_conflict_force * conflict_torque_arm;

            std::cout << "conflict_torque_arm: " << conflict_torque_arm << std::endl;
            std::cout << "extruder_conflict_torque: " << extruder_conflict_torque << std::endl;

            float total_torque = bed_movement_torque + extruder_conflict_torque - weight_torque - sticking_torque;

            std::cout << "total_torque: " << total_torque << "   printz: " << print_z << std::endl;

            if (total_torque > 0) {
                acc->add_support_point(Point::new_scale(line.b), std::min(params.support_points_interface_area, (pivot - line.b).squaredNorm()));
                issues.supports_nedded.emplace_back(to_3d(line.b, print_z), extruder_conflict_torque - sticking_torque);
                distance_from_last_support_point = 0.0f;
            }

        }
    }
}

Issues check_object_stability(const PrintObject *po, const Params &params) {
#ifdef DEBUG_FILES
    FILE *debug_acc = boost::nowide::fopen(debug_out_path("accumulators.obj").c_str(), "w");
#endif
    StabilityAccumulators stability_accs;
    LayerLinesDistancer prev_layer_lines { { } };
    Issues issues { };
    std::vector<ExtrusionLine> checked_lines;

    // PREPARE BASE LAYER
    float max_flow_width = 0.0f;
    const Layer *layer = po->layers()[0];
    float base_print_z = layer->print_z;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                const float flow_width = get_flow_width(layer_region, perimeter->role());
                max_flow_width = std::max(flow_width, max_flow_width);
                const float mm3_per_mm = float(perimeter->min_mm3_per_mm());
                int id = stability_accs.create_accumulator(base_print_z);
                StabilityAccumulator &acc = stability_accs.access(id);
                Points points { };
                perimeter->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    acc.add_base_extrusion(line, flow_width, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
                if (perimeter->is_loop()) {
                    Vec2f start = unscaled(points[points.size() - 1]).cast<float>();
                    Vec2f next = unscaled(points[0]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    acc.add_base_extrusion(line, flow_width, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
            } // perimeter
        } // ex_entity
        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                const float flow_width = get_flow_width(layer_region, fill->role());
                max_flow_width = std::max(flow_width, max_flow_width);
                const float mm3_per_mm = float(fill->min_mm3_per_mm());
                int id = stability_accs.create_accumulator(base_print_z);
                StabilityAccumulator &acc = stability_accs.access(id);
                Points points { };
                fill->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    acc.add_base_extrusion(line, flow_width, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
            } // fill
        } // ex_entity
    } // region

#ifdef DEBUG_FILES
    for (const auto &line : checked_lines) {
        Vec3f color = stability_accs.get_accumulator_color(line.stability_accumulator_id);
        fprintf(debug_acc, "v %f %f %f  %f %f %f\n", line.b[0],
                line.b[1], base_print_z, color[0], color[1], color[2]);
    }
#endif

    //MERGE BASE LAYER STABILITY ACCS
    prev_layer_lines = LayerLinesDistancer { std::move(checked_lines) };
    for (const ExtrusionLine &l : prev_layer_lines.get_lines()) {
        size_t nearest_line_idx;
        Vec2f nearest_pt;
        float dist = prev_layer_lines.signed_distance_from_lines(l.a, nearest_line_idx, nearest_pt);
        if (std::abs(dist) < max_flow_width*1.1f) {
            size_t other_line_acc_id = prev_layer_lines.get_line(nearest_line_idx).stability_accumulator_id;
            size_t from_id = std::max(other_line_acc_id, l.stability_accumulator_id);
            size_t to_id = std::min(other_line_acc_id, l.stability_accumulator_id);
            stability_accs.merge_accumulators(from_id, to_id);
        }
    }

    //CHECK STABILITY OF ALL LAYERS
    for (size_t layer_idx = 1; layer_idx < po->layer_count(); ++layer_idx) {
        const Layer *layer = po->layers()[layer_idx];
        checked_lines = std::vector<ExtrusionLine> { };
        std::vector<std::pair<Vec2f, size_t>> fill_points;
        float max_fill_flow_width = 0.0f;

        float print_z = layer->print_z;
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    check_extrusion_entity_stability(perimeter, stability_accs, issues, checked_lines, print_z,
                            layer_region,
                            prev_layer_lines, params);
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill
                            || fill->role() == ExtrusionRole::erBridgeInfill) {
                        check_extrusion_entity_stability(fill, stability_accs, issues, checked_lines, print_z,
                                layer_region,
                                prev_layer_lines, params);
                    } else {
                        const float flow_width = get_flow_width(layer_region, fill->role());
                        max_fill_flow_width = std::max(max_fill_flow_width, flow_width);
                        Vec2f start = unscaled(fill->first_point()).cast<float>();
                        size_t nearest_line_idx;
                        Vec2f nearest_pt;
                        float dist = prev_layer_lines.signed_distance_from_lines(start, nearest_line_idx, nearest_pt);
                        if (dist < flow_width) {
                            size_t acc_id = prev_layer_lines.get_line(nearest_line_idx).stability_accumulator_id;
                            StabilityAccumulator &acc = stability_accs.access(acc_id);
                            Points points { };
                            const float mm3_per_mm = float(fill->min_mm3_per_mm());
                            fill->collect_points(points);
                            for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                                Vec2f start = unscaled(points[point_idx]).cast<float>();
                                Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                                ExtrusionLine line { start, next };
                                line.stability_accumulator_id = acc_id;
                                acc.add_extrusion(line, print_z, mm3_per_mm);
                            }
                            fill_points.emplace_back(start, acc_id);
                        } else {
                            std::cout << " SUPPORTS POINT GEN, start infill in the air? on printz: " << print_z
                                    << std::endl;
                        }
                    }
                } // fill
            } // ex_entity
        } // region

        prev_layer_lines = LayerLinesDistancer { std::move(checked_lines) };

        for (const std::pair<Vec2f, size_t> &fill_point : fill_points) {
            size_t nearest_line_idx;
            Vec2f nearest_pt;
            float dist = prev_layer_lines.signed_distance_from_lines(fill_point.first, nearest_line_idx, nearest_pt);
            if (dist < max_fill_flow_width) {
                size_t other_line_acc_id = prev_layer_lines.get_line(nearest_line_idx).stability_accumulator_id;
                size_t from_id = std::max(other_line_acc_id, fill_point.second);
                size_t to_id = std::min(other_line_acc_id, fill_point.second);
                stability_accs.merge_accumulators(from_id, to_id);
            } else {
                std::cout << " SUPPORTS POINT GEN, no connection on current layer for infill? on printz: " << print_z
                        << std::endl;
            }
        }

        check_layer_global_stability(stability_accs, issues, prev_layer_lines.get_lines(), print_z, params);

#ifdef DEBUG_FILES
        for (const auto &line : prev_layer_lines.get_lines()) {
            Vec3f color = stability_accs.get_accumulator_color(line.stability_accumulator_id);
            fprintf(debug_acc, "v %f %f %f  %f %f %f\n", line.b[0],
                    line.b[1], print_z, color[0], color[1], color[2]);
        }
#endif
    }

#ifdef DEBUG_FILES
    fclose(debug_acc);
#endif

    std::cout << " SUPP: " << issues.supports_nedded.size() << std::endl;
    return issues;
}

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
            fprintf(fp, "v %f %f %f  %f %f %f\n", issues.supports_nedded[i].position(0),
                    issues.supports_nedded[i].position(1),
                    issues.supports_nedded[i].position(2), 1.0, 0.0, 1.0);
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
            fprintf(fp, "v %f %f %f  %f %f %f\n", issues.curling_up[i].position(0),
                    issues.curling_up[i].position(1),
                    issues.curling_up[i].position(2), 0.0, 1.0, 0.0);
        }
        fclose(fp);
    }
}
#endif

std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
    check_object_stability(po, params);
    return {};
}

Issues full_search(const PrintObject *po, const Params &params) {
    auto issues = check_object_stability(po, params);
    debug_export(issues, "issues");
    return issues;

}
} //SupportableIssues End
}

