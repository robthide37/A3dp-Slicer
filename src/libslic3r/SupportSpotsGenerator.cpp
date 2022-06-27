#include "SupportSpotsGenerator.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <unordered_set>
#include <stack>

#include "AABBTreeLines.hpp"
#include "KDTreeIndirect.hpp"
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

    float malformation = 0.0f;
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

namespace SupportSpotsGenerator {

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

struct VoxelGrid {
private:
    Vec3f cell_size;
    Vec3f origin;
    Vec3f size;
    Vec3i cell_count;

    std::unordered_set<size_t> taken_cells { };

public:
    VoxelGrid(const PrintObject *po, float voxel_size) {
        cell_size = Vec3f(voxel_size, voxel_size, voxel_size);

        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec3f min = unscale(Vec3crd(-size_half.x(), -size_half.y(), 0)).cast<float>() - cell_size;
        Vec3f max = unscale(Vec3crd(size_half.x(), size_half.y(), po->height())).cast<float>() + cell_size;

        origin = min;
        size = max - min;
        cell_count = size.cwiseQuotient(cell_size).cast<int>() + Vec3i::Ones();
    }

    Vec3i to_cell_coords(const Vec3f &position) const {
        Vec3i cell_coords = (position - this->origin).cwiseQuotient(this->cell_size).cast<int>();
        return cell_coords;
    }

    size_t to_cell_index(const Vec3i &cell_coords) const {
        assert(cell_coords.x() >= 0);
        assert(cell_coords.x() < cell_count.x());
        assert(cell_coords.y() >= 0);
        assert(cell_coords.y() < cell_count.y());
        assert(cell_coords.z() >= 0);
        assert(cell_coords.z() < cell_count.z());

        return cell_coords.z() * cell_count.x() * cell_count.y()
                + cell_coords.y() * cell_count.x()
                + cell_coords.x();
    }

    Vec3f get_cell_center(const Vec3i &cell_coords) const {
        return origin + cell_coords.cast<float>().cwiseProduct(this->cell_size)
                + this->cell_size.cwiseQuotient(Vec3f(2.0f, 2.0f, 2.0));
    }

    void take_position(const Vec3f &position) {
        taken_cells.insert(to_cell_index(to_cell_coords(position)));
    }

    bool position_taken(const Vec3f &position) const {
        return taken_cells.find(to_cell_index(to_cell_coords(position))) != taken_cells.end();
    }

};

class LayerLinesDistancer {
private:
    std::vector<ExtrusionLine> lines;
    AABBTreeIndirect::Tree<2, float> tree;

public:
    explicit LayerLinesDistancer(std::vector<ExtrusionLine> &&lines) :
            lines(lines) {
        tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(this->lines);
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

// StabilityAccumulator accumulates extrusions for each connected model part from bed to current printed layer.
//  If the originaly disconected parts meet in the layer, their stability accumulators get merged and continue as one.
// (think legs of table, which get connected when the top desk is being printed).
// The class gathers mass, centroid mass, sticking force (bed extrusions, support points) and sticking centroid for the
// connected part. These values are then used to check global part stability.
class StabilityAccumulator {
private:
    std::vector<Vec2f> support_points { };
    Vec3f centroid_accumulator = Vec3f::Zero();
    float accumulated_volume { };
    Vec2f sticking_centroid_accumulator = Vec2f::Zero();
    float accumulated_sticking_force { };

public:
    StabilityAccumulator() = default;

    void add_base_extrusion(const ExtrusionLine &line, float sticking_force, float print_z, float mm3_per_mm) {
        accumulated_sticking_force += sticking_force;
        sticking_centroid_accumulator += sticking_force * ((line.a + line.b) / 2.0f);
        support_points.push_back(line.a);
        support_points.push_back(line.b);
        add_extrusion(line, print_z, mm3_per_mm);
    }

    void add_support_point(const Vec2f &position, float sticking_force) {
        support_points.push_back(position);
        accumulated_sticking_force += sticking_force;
        sticking_centroid_accumulator += sticking_force * position;
    }

    void add_extrusion(const ExtrusionLine &line, float print_z, float mm3_per_mm) {
        float volume = line.len * mm3_per_mm;
        accumulated_volume += volume;
        Vec2f center = (line.a + line.b) / 2.0f;
        centroid_accumulator += volume * Vec3f(center.x(), center.y(), print_z);
    }

    Vec3f get_centroid() const {
        if (accumulated_volume <= 0.0f) {
            return Vec3f::Zero();
        }
        return centroid_accumulator / accumulated_volume;
    }

    float get_sticking_force() const {
        return accumulated_sticking_force;
    }

    float get_accumulated_volume() const {
        return accumulated_volume;
    }

    const std::vector<Vec2f>& get_support_points() const {
        return support_points;
    }

    Vec2f get_sticking_centroid() const {
        if (accumulated_sticking_force <= 0.0f) {
            return Vec2f::Zero();
        }
        return sticking_centroid_accumulator / accumulated_sticking_force;
    }

    void add_from(const StabilityAccumulator &acc) {
        this->support_points.insert(this->support_points.end(), acc.support_points.begin(),
                acc.support_points.end());
        this->centroid_accumulator += acc.centroid_accumulator;
        this->accumulated_volume += acc.accumulated_volume;
        this->accumulated_sticking_force += acc.accumulated_sticking_force;
        this->sticking_centroid_accumulator += acc.sticking_centroid_accumulator;
    }
};

// StabilityAccumulators class is wrapper over the vector of stability accumualtors. It provides a level of indirection
// between accumulator ID and the accumulator instance itself. While each extrusion line has one id, which is asigned
// when algorithm reaches the line's layer, the accumulator this id points to can change due to merging.
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
        from_acc = StabilityAccumulator { };

    }

public:
    StabilityAccumulators() = default;

    int create_accumulator() {
        size_t id = next_id;
        next_id++;
        mapping[id] = accumulators.size();
        accumulators.push_back(StabilityAccumulator { });
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
        from_acc = StabilityAccumulator { };
    }

#ifdef DEBUG_FILES
    Vec3f get_accumulator_color(size_t id) {
        if (mapping.find(id) == mapping.end()) {
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: ERROR: uknown accumulator ID: " << id;
            return Vec3f(1.0f, 1.0f, 1.0f);
        }

        size_t pseudornd = ((mapping[id] + 127) * 33331 + 6907) % 987;
        return value_to_rgbf(0.0f, float(987), float(pseudornd));
    }

    void log_accumulators() {
        for (size_t i = 0; i < accumulators.size(); ++i) {
            const auto &acc = accumulators[i];
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: accumulator POS: " << i << "\n"
                    << "SSG: get_accumulated_volume: " << acc.get_accumulated_volume() << "\n"
                    << "SSG: get_sticking_force: " << acc.get_sticking_force() << "\n"
                    << "SSG: support points count: " << acc.get_support_points().size() << "\n";

        }
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

// Accumulator of current extruion path properties
// It remembers unsuported distance and maximum accumulated curvature over that distance.
// Used to determine local stability issues (too long bridges, extrusion curves into air)
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

// check_extrusion_entity_stability checks each extrusion for local issues, appends the extrusion
// into checked lines, and gives it a stability accumulator id. If support is needed it pushes it
// into issues as well.
// Rules for stability accumulator id assigment:
// If there is close extrusion under, use min extrusion id between the id of the previous line,
//      and id of line under. Also merge the accumulators of those two ids!
// If there is no close extrusion under, use id of the previous extrusion line.
// If there is no previous line, create new stability accumulator.
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
        const auto to_vec3f = [print_z](const Vec2f &point) {
            return Vec3f(point.x(), point.y(), print_z);
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
        ExtrusionPropertiesAccumulator malformation_acc { };
        bridging_acc.add_distance(params.bridge_distance + 1.0f); // Initialise unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter starts and short loops into air.
        const float flow_width = get_flow_width(layer_region, entity->role());

        for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx) {
            ExtrusionLine &current_line = lines[line_idx];
            float mm3_per_mm = float(entity->min_mm3_per_mm());

            float curr_angle = 0;
            if (line_idx + 1 < lines.size()) {
                const Vec2f v1 = current_line.b - current_line.a;
                const Vec2f v2 = lines[line_idx + 1].b - lines[line_idx + 1].a;
                curr_angle = angle(v1, v2);
            }
            bridging_acc.add_angle(curr_angle);
            malformation_acc.add_angle(curr_angle);

            size_t nearest_line_idx;
            Vec2f nearest_point;
            float dist_from_prev_layer = prev_layer_lines.signed_distance_from_lines(current_line.b, nearest_line_idx,
                    nearest_point);

            if (dist_from_prev_layer < flow_width) {
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
                    current_stability_acc = stability_accs.create_accumulator();
                }
                StabilityAccumulator &current_segment = stability_accs.access(current_stability_acc);
                current_line.stability_accumulator_id = current_stability_acc;
                current_segment.add_extrusion(current_line, print_z, mm3_per_mm);
                if (bridging_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                > params.bridge_distance
                        / (1.0f + bridging_acc.max_curvature
                                * params.bridge_distance_decrease_by_curvature_factor / PI)) {
                    current_segment.add_support_point(current_line.b, 0.0f); // Do not count extrusion supports into the sticking force. They can be very densely placed, causing algorithm to overestimate stickiness.
                    issues.supports_nedded.emplace_back(to_vec3f(current_line.b), 1.0);
                    bridging_acc.reset();
                }
            }

            //malformation
            if (fabs(dist_from_prev_layer) < flow_width * 2.0f) {
                const ExtrusionLine &nearest_line = prev_layer_lines.get_line(nearest_line_idx);
                current_line.malformation += 0.7 * nearest_line.malformation;
            }
            if (dist_from_prev_layer > flow_width * 0.3) {
                current_line.malformation += 0.6 + 0.4 * malformation_acc.max_curvature / PI;
            } else {
                malformation_acc.reset();
            }
        }
        checked_lines.insert(checked_lines.end(), lines.begin(), lines.end());
    }
}

//check_layer_global_stability checks stability of the accumulators that are still present on the current line
// ( this is determined from the gathered checked_lines vector)
// For each accumulator and each its extrusion, forces and torques (weight, bed movement, extruder pressure, stickness to bed)
// are computed and if stability is not sufficient, support points are added
// accumualtors are filtered by their pointer address, since one accumulator can have multiple IDs due to merging
void check_layer_global_stability(StabilityAccumulators &stability_accs,
        VoxelGrid &supports_presence_grid,
        Issues &issues,
        float flow_width,
        const std::vector<ExtrusionLine> &checked_lines,
        float print_z,
        const Params &params,
        std::mt19937_64& generator) {
    std::unordered_map<StabilityAccumulator*, std::vector<ExtrusionLine>> layer_accs_w_lines;
    for (size_t i = 0; i < checked_lines.size(); ++i) {
        layer_accs_w_lines[&stability_accs.access(checked_lines[i].stability_accumulator_id)].push_back(
                checked_lines[i]);
    }

    for (auto &accumulator : layer_accs_w_lines) {
        StabilityAccumulator *acc = accumulator.first;
        std::shuffle(accumulator.second.begin(), accumulator.second.end(), generator);
        LayerLinesDistancer acc_lines(std::move(accumulator.second));

        if (acc->get_support_points().empty()) {
            // acc_lines cannot be empty - if the accumulator has no extrusion in the current layer, it is not considered in stability computation
            acc->add_support_point(acc_lines.get_line(0).a, 0.0f);
            issues.supports_nedded.emplace_back(to_3d(acc_lines.get_line(0).a, print_z), 0.0);
        }
        const std::vector<Vec2f> &support_points = acc->get_support_points();

        auto coord_fn = [&support_points](size_t idx, size_t dim) {
            return support_points[idx][dim];
        };
        KDTreeIndirect<2, float, decltype(coord_fn)> supports_tree(coord_fn, support_points.size());

        for (const ExtrusionLine &line : acc_lines.get_lines()) {
            Vec2f line_dir = (line.b - line.a).normalized();
            Vec2f pivot_site_search_point = line.b + line_dir * 300.0f;
            size_t pivot_idx = find_closest_point(supports_tree, pivot_site_search_point);
            const Vec2f &pivot = support_points[pivot_idx];

            const Vec2f &sticking_centroid = acc->get_sticking_centroid();
            float sticking_arm = (pivot - sticking_centroid).norm();
            float sticking_torque = sticking_arm * acc->get_sticking_force();

            float mass = acc->get_accumulated_volume() * params.filament_density;
            const Vec3f &mass_centorid = acc->get_centroid();
            float weight = mass * params.gravity_constant;
            float weight_arm = (pivot - mass_centorid.head<2>()).norm();
            float weight_torque = weight_arm * weight;

            float bed_movement_arm = mass_centorid.z();
            float bed_movement_force = params.max_acceleration * mass;
            float bed_movement_torque = bed_movement_force * bed_movement_arm;

            Vec3f extruder_pressure_direction = to_3d(line_dir, 0.0f);
            extruder_pressure_direction.z() = -0.2 - line.malformation * 0.5;
            extruder_pressure_direction.normalize();
            float conflict_torque_arm = (to_3d(Vec2f(pivot - line.b), print_z).cross(
                    extruder_pressure_direction)).norm();
            float extruder_conflict_force = params.tolerable_extruder_conflict_force +
                    line.malformation * params.malformations_additive_conflict_extruder_force;
            float extruder_conflict_torque = extruder_conflict_force * conflict_torque_arm;

            float total_torque = bed_movement_torque + extruder_conflict_torque - weight_torque - sticking_torque;

            if (total_torque > 0) {
                Vec2f target_point;
                size_t _idx;
                acc_lines.signed_distance_from_lines(pivot_site_search_point, _idx, target_point);
                if (!supports_presence_grid.position_taken(to_3d(target_point, print_z))) {
                    float area = params.support_points_interface_radius * params.support_points_interface_radius
                            * float(PI);
                    float sticking_force = area * params.support_adhesion;
                    acc->add_support_point(target_point, sticking_force);
                    issues.supports_nedded.emplace_back(to_3d(target_point, print_z),
                            extruder_conflict_torque - sticking_torque);
                    supports_presence_grid.take_position(to_3d(target_point, print_z));
                }
            }
#if 1
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: sticking_arm: " << sticking_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: sticking_torque: " << sticking_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: weight_arm: " << sticking_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: weight_torque: " << weight_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_movement_arm: " << bed_movement_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_movement_torque: " << bed_movement_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conflict_torque_arm: " << conflict_torque_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: extruder_conflict_torque: " << extruder_conflict_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: total_torque: " << total_torque << "   printz: " << print_z;
#endif
        }
    }
}

Issues check_object_stability(const PrintObject *po, const Params &params) {
#ifdef DEBUG_FILES
    FILE *debug_acc = boost::nowide::fopen(debug_out_path("accumulators.obj").c_str(), "w");
    FILE *malform_f = boost::nowide::fopen(debug_out_path("malformations.obj").c_str(), "w");
#endif
    StabilityAccumulators stability_accs;
    LayerLinesDistancer prev_layer_lines { { } };
    Issues issues { };
    std::vector<ExtrusionLine> checked_lines;
    VoxelGrid supports_presence_grid { po, params.min_distance_between_support_points };
    std::mt19937_64 generator { 27644437 };

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
                int id = stability_accs.create_accumulator();
                StabilityAccumulator &acc = stability_accs.access(id);
                Points points { };
                perimeter->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    float line_sticking_force = line.len * flow_width * params.base_adhesion;
                    acc.add_base_extrusion(line, line_sticking_force, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
                if (perimeter->is_loop()) {
                    Vec2f start = unscaled(points[points.size() - 1]).cast<float>();
                    Vec2f next = unscaled(points[0]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    float line_sticking_force = line.len * flow_width * params.base_adhesion;
                    acc.add_base_extrusion(line, line_sticking_force, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
            } // perimeter
        } // ex_entity
        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                const float flow_width = get_flow_width(layer_region, fill->role());
                max_flow_width = std::max(flow_width, max_flow_width);
                const float mm3_per_mm = float(fill->min_mm3_per_mm());
                int id = stability_accs.create_accumulator();
                StabilityAccumulator &acc = stability_accs.access(id);
                Points points { };
                fill->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next };
                    line.stability_accumulator_id = id;
                    float line_sticking_force = line.len * flow_width * params.base_adhesion;
                    acc.add_base_extrusion(line, line_sticking_force, base_print_z, mm3_per_mm);
                    checked_lines.push_back(line);
                }
            } // fill
        } // ex_entity
    } // region

    //MERGE BASE LAYER STABILITY ACCS
    prev_layer_lines = LayerLinesDistancer { std::move(checked_lines) };
    for (const ExtrusionLine &l : prev_layer_lines.get_lines()) {
        size_t nearest_line_idx;
        Vec2f nearest_pt;
        Vec2f line_dir = (l.b - l.a).normalized();
        Vec2f site_search_location = l.a + Vec2f(line_dir.y(), -line_dir.x()) * max_flow_width;
        float dist = prev_layer_lines.signed_distance_from_lines(site_search_location, nearest_line_idx, nearest_pt);
        if (std::abs(dist) < max_flow_width) {
            size_t other_line_acc_id = prev_layer_lines.get_line(nearest_line_idx).stability_accumulator_id;
            size_t from_id = std::max(other_line_acc_id, l.stability_accumulator_id);
            size_t to_id = std::min(other_line_acc_id, l.stability_accumulator_id);
            stability_accs.merge_accumulators(from_id, to_id);
        }
    }

#ifdef DEBUG_FILES
    for (const auto &line : prev_layer_lines.get_lines()) {
        Vec3f color = stability_accs.get_accumulator_color(line.stability_accumulator_id);
        fprintf(debug_acc, "v %f %f %f  %f %f %f\n", line.b[0],
                line.b[1], base_print_z, color[0], color[1], color[2]);
    }

    stability_accs.log_accumulators();
#endif

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
                            BOOST_LOG_TRIVIAL(debug)
                            << "SSG: ERROR: seem that infill starts in the air? on printz: " << print_z;
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
                BOOST_LOG_TRIVIAL(debug)
                << "SSG: ERROR: seem that infill starts in the air? on printz: " << print_z;
            }
        }

        check_layer_global_stability(stability_accs,
                supports_presence_grid,
                issues,
                max_flow_width,
                prev_layer_lines.get_lines(),
                print_z,
                params,
                generator);

#ifdef DEBUG_FILES
        for (const auto &line : prev_layer_lines.get_lines()) {
            Vec3f color = value_to_rgbf(0, 5.0f, line.malformation);
            fprintf(malform_f, "v %f %f %f  %f %f %f\n", line.b[0],
                    line.b[1], print_z, color[0], color[1], color[2]);
        }
        for (const auto &line : prev_layer_lines.get_lines()) {
            Vec3f color = stability_accs.get_accumulator_color(line.stability_accumulator_id);
            fprintf(debug_acc, "v %f %f %f  %f %f %f\n", line.b[0],
                    line.b[1], print_z, color[0], color[1], color[2]);
        }
        stability_accs.log_accumulators();
#endif
    }

#ifdef DEBUG_FILES
    fclose(debug_acc);
    fclose(malform_f);
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
#ifdef DEBUG_FILES
    debug_export(issues, "issues");
#endif
    return issues;

}
} //SupportableIssues End
}

