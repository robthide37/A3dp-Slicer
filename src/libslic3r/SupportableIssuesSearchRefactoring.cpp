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

    size_t supported_segment_accumulator_id = NULL_ACC_ID;

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
    float signed_distance_from_lines(const Point &point, size_t &nearest_line_index_out,
            Vec2f &nearest_point_out) const {
        Vec2f p = unscaled(point).cast<float>();
        auto distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, p, nearest_line_index_out,
                nearest_point_out);
        if (distance < 0)
            return std::numeric_limits<float>::infinity();

        distance = sqrt(distance);
        const ExtrusionLine &line = lines[nearest_line_index_out];
        Vec2f v1 = line.b - line.a;
        Vec2f v2 = p - line.a;
        if ((v1.x() * v2.y()) - (v1.y() * v2.x()) > 0.0) {
            distance *= -1;
        }
        return distance;
    }

    const ExtrusionLine& get_line(size_t line_idx) const {
        return lines[line_idx];
    }
};

class SupportedSegmentAccumulator {
private:
    Polygon base_convex_hull { };
    Points supported_points { };
    Vec3f centroid_accumulator = Vec3f::Zero();
    float accumulated_volume { };
    float base_area { };
    float base_height { };

public:
    explicit SupportedSegmentAccumulator(float base_height) :
            base_height(base_height) {
    }

    void add_base_extrusion(const ExtrusionLine &line, float width, float print_z, float cross_section) {
        base_area += line.len * width;
        supported_points.push_back(Point::new_scale(line.a));
        supported_points.push_back(Point::new_scale(line.b));
        base_convex_hull.clear();
        add_extrusion(line, print_z, cross_section);
    }

    void add_support_point(const Point &position, float area) {
        supported_points.push_back(position);
        base_convex_hull.clear();
        base_area += area;
    }

    void add_extrusion(const ExtrusionLine &line, float print_z, float cross_section) {
        float volume = line.len * cross_section;
        accumulated_volume += volume;
        Vec2f center = (line.a + line.b) / 2.0f;
        centroid_accumulator += volume * Vec3f(center.x(), center.y(), print_z);
    }

    const Polygon& segment_base_hull() {
        if (this->base_convex_hull.empty()) {
            this->base_convex_hull = Geometry::convex_hull(this->supported_points);
        }
        return this->base_convex_hull;
    }

    void add_from(const SupportedSegmentAccumulator &acc) {
        this->supported_points.insert(this->supported_points.end(), acc.supported_points.begin(),
                acc.supported_points.end());
        base_convex_hull.clear();
        this->centroid_accumulator += acc.centroid_accumulator;
        this->accumulated_volume += acc.accumulated_volume;
        this->base_area += acc.base_area;
    }

    bool check_stability() {
        return true;
    }
};

struct SupportedSegmentAccumulators {
private:
    size_t next_id = 0;
    std::unordered_map<size_t, size_t> mapping;
    std::vector<SupportedSegmentAccumulator> acccumulators;

    void merge_to(size_t from_id, size_t to_id) {
        SupportedSegmentAccumulator &from_acc = this->access(from_id);
        SupportedSegmentAccumulator &to_acc = this->access(to_id);
        if (&from_acc == &to_acc) {
            return;
        }
        to_acc.add_from(from_acc);
        mapping[from_id] = mapping[to_id];
        from_acc = SupportedSegmentAccumulator { 0.0f };

    }

public:
    SupportedSegmentAccumulators() = default;

    int create_accumulator(float base_height) {
        size_t id = next_id;
        next_id++;
        mapping[id] = acccumulators.size();
        acccumulators.push_back(SupportedSegmentAccumulator { base_height });
        return id;
    }

    SupportedSegmentAccumulator& access(size_t id) {
        return acccumulators[mapping[id]];
    }

    void merge_accumulators(size_t from_id, size_t to_id) {
        if (from_id == NULL_ACC_ID || to_id == NULL_ACC_ID) {
            return;
        }
        SupportedSegmentAccumulator &from_acc = this->access(from_id);
        SupportedSegmentAccumulator &to_acc = this->access(to_id);
        if (&from_acc == &to_acc) {
            return;
        }
        to_acc.add_from(from_acc);
        mapping[from_id] = mapping[to_id];
        from_acc = SupportedSegmentAccumulator { 0.0f };
    }

    std::unordered_set<size_t> get_active_acc_indices() const {
        std::unordered_set<size_t> result;
        for (const auto &pair : mapping) {
            result.insert(pair.second);
        }
        return result;
    }
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

float get_max_allowed_distance(ExtrusionRole role, float flow_width, bool external_perimeters_first,
        const Params &params) { // <= distance / flow_width (can be larger for perimeter, if not external perimeter first)
    if ((role == ExtrusionRole::erExternalPerimeter || role == ExtrusionRole::erOverhangPerimeter)
            && (external_perimeters_first)) {
        return params.max_first_ex_perim_unsupported_distance_factor * flow_width;
    } else {
        return params.max_unsupported_distance_factor * flow_width;
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
        SupportedSegmentAccumulators &stability_accs,
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
            Vec2f next = unscaled(points[point_idx]).cast<float>();
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

        checked_lines.insert(checked_lines.end(), lines.begin(), lines.end());

        size_t current_stability_acc = NULL_ACC_ID;
        ExtrusionPropertiesAccumulator bridging_acc { };
        bridging_acc.add_distance(params.bridge_distance + 1.0f); // Initialise unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter start and short loops into air.
        const float flow_width = get_flow_width(layer_region, entity->role());
        const float region_height = layer_region->layer()->height;
        const float max_allowed_dist_from_prev_layer = get_max_allowed_distance(entity->role(), flow_width,
                layer_region->region().config().external_perimeters_first, params);

        for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx) {
            ExtrusionLine current_line = lines[line_idx];
            Point current = Point::new_scale(current_line.b);
            float cross_section = region_height * flow_width * 0.7071f;

            float angle = 0;
            if (line_idx + 1 < lines.size()) {
                const Vec2f v1 = current_line.b - current_line.a;
                const Vec2f v2 = lines[line_idx + 1].b - lines[line_idx + 1].a;
                float dot = v1(0) * v2(0) + v1(1) * v2(1);
                float cross = v1(0) * v2(1) - v1(1) * v2(0);
                angle = float(atan2(float(cross), float(dot))); // ccw angle, TODO replace with angle func, once it gets into master
            }
            bridging_acc.add_angle(angle);

            size_t nearest_line_idx;
            Vec2f nearest_point;
            float dist_from_prev_layer = prev_layer_lines.signed_distance_from_lines(current, nearest_line_idx,
                    nearest_point);
            if (dist_from_prev_layer - flow_width < max_allowed_dist_from_prev_layer) {
                const ExtrusionLine &nearest_line = prev_layer_lines.get_line(nearest_line_idx);
                size_t acc_id = nearest_line.supported_segment_accumulator_id;
                stability_accs.merge_accumulators(std::max(acc_id, current_stability_acc),
                        std::min(acc_id, current_stability_acc));
                current_stability_acc = std::min(acc_id, current_stability_acc);
                current_line.supported_segment_accumulator_id = current_stability_acc;
                stability_accs.access(current_stability_acc).add_extrusion(current_line, print_z, cross_section);
                bridging_acc.reset();
                // TODO curving here
            } else {
                bridging_acc.add_distance(current_line.len);
                if (current_stability_acc < 0) {
                    size_t acc_id = stability_accs.create_accumulator(print_z);
                    stability_accs.merge_accumulators(std::max(acc_id, current_stability_acc),
                            std::min(acc_id, current_stability_acc));
                    current_stability_acc = std::min(acc_id, current_stability_acc);
                }
                SupportedSegmentAccumulator &current_segment = stability_accs.access(current_stability_acc);
                current_segment.add_extrusion(current_line, print_z, cross_section);
                if (bridging_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                > params.bridge_distance
                        / (1.0f + (bridging_acc.max_curvature
                                * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    current_segment.add_support_point(current, 5.0f);
                    issues.supports_nedded.emplace_back(to_vec3f(current), 1.0);
                    bridging_acc.reset();
                }
            }
        }
    }
}

Issues check_object_stability(const PrintObject *po, const Params &params) {
    SupportedSegmentAccumulators stability_accs;
    LayerLinesDistancer prev_layer_lines { { } };
    Issues issues { };
    std::vector<ExtrusionLine> checked_lines;

    const Layer *layer = po->layers()[0];
    float base_print_z = layer->print_z;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                const float flow_width = get_flow_width(layer_region, perimeter->role());
                const float region_height = layer_region->layer()->height;
                const float cross_section = region_height * flow_width * 0.7071f;
                int id = stability_accs.create_accumulator(base_print_z);
                SupportedSegmentAccumulator &acc = stability_accs.access(id);
                Points points { };
                perimeter->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx]).cast<float>();
                    ExtrusionLine line{start, next};
                    line.supported_segment_accumulator_id = id;
                    acc.add_base_extrusion( line, flow_width, base_print_z, cross_section);
                    checked_lines.push_back(line);
                }
            } // perimeter
        } // ex_entity
        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                const float flow_width = get_flow_width(layer_region, fill->role());
                const float region_height = layer_region->layer()->height;
                const float cross_section = region_height * flow_width * 0.7071f;
                int id = stability_accs.create_accumulator(base_print_z);
                SupportedSegmentAccumulator &acc = stability_accs.access(id);
                Points points { };
                fill->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx]).cast<float>();
                    acc.add_base_extrusion( { start, next }, flow_width, base_print_z, cross_section);
                }
            } // fill
        } // ex_entity
    } // region

    for (size_t layer_idx = 1; layer_idx < po->layer_count(); ++layer_idx) {
        const Layer *layer = po->layers()[layer_idx];
        prev_layer_lines = LayerLinesDistancer{std::move(checked_lines)};
        checked_lines = std::vector<ExtrusionLine>{};

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
                    }
                } // fill
            } // ex_entity
        } // region
    }

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

