#ifndef slic3r_ExtrusionProcessor_hpp_
#define slic3r_ExtrusionProcessor_hpp_

#include "../AABBTreeLines.hpp"
#include "../SupportSpotsGenerator.hpp"
#include "../libslic3r.h"
#include "../ExtrusionEntity.hpp"
#include "../Layer.hpp"
#include "../Point.hpp"
#include "../SVG.hpp"
#include "../BoundingBox.hpp"
#include "../Polygon.hpp"
#include "../ClipperUtils.hpp"
#include "../Flow.hpp"
#include "../Config.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Slic3r {

class SlidingWindowCurvatureAccumulator
{
    float        window_size;
    float        total_distance  = 0; // accumulated distance
    float        total_curvature = 0; // accumulated signed ccw angles
    deque<float> distances;
    deque<float> angles;

public:
    SlidingWindowCurvatureAccumulator(float window_size) : window_size(window_size) {}

    void add_point(float distance, float angle)
    {
        total_distance += distance;
        total_curvature += std::abs(angle);
        distances.push_back(distance);
        angles.push_back(std::abs(angle));

        while (distances.size() > 1 && total_distance > window_size) {
            total_distance -= distances.front();
            total_curvature -= angles.front();
            distances.pop_front();
            angles.pop_front();
        }
    }

    float get_curvature() const
    {
        if (total_distance < EPSILON) { return 0.0; }

        return total_curvature / window_size;
    }

    void reset()
    {
        total_curvature = 0;
        total_distance  = 0;
        distances.clear();
        angles.clear();
    }
};

class CurvatureEstimator
{
    static const size_t               sliders_count          = 2;
    SlidingWindowCurvatureAccumulator sliders[sliders_count] = {{1.5}, {3.0}};

public:
    void add_point(float distance, float angle)
    {
        if (distance < EPSILON) return;
        for (SlidingWindowCurvatureAccumulator &slider : sliders) { slider.add_point(distance, angle); }
    }
    float get_curvature()
    {
        float max_curvature = std::numeric_limits<float>::min();
        for (const SlidingWindowCurvatureAccumulator &slider : sliders) { max_curvature = std::max(max_curvature, slider.get_curvature()); }
        return max_curvature;
    }
    void reset()
    {
        for (SlidingWindowCurvatureAccumulator &slider : sliders) { slider.reset(); }
    }
};

struct ExtendedPoint
{
    ExtendedPoint(Vec2d position, float distance = 0.0, size_t nearest_prev_layer_line = size_t(-1), float curvature = 0.0)
        : position(position), distance(distance), nearest_prev_layer_line(nearest_prev_layer_line), curvature(curvature)
    {}

    Vec2d  position;
    float  distance;
    size_t nearest_prev_layer_line;
    float  curvature;
};

template<bool SCALED_INPUT, bool ADD_INTERSECTIONS, bool PREV_LAYER_BOUNDARY_ONLY, bool CONCAVITY_RESETS_CURVATURE, typename P, typename L>
std::vector<ExtendedPoint> estimate_points_properties(const std::vector<P>                   &extrusion_points,
                                                      const AABBTreeLines::LinesDistancer<L> &unscaled_prev_layer,
                                                      float                                   flow_width)
{
    if (extrusion_points.empty()) return {};
    float              boundary_offset = PREV_LAYER_BOUNDARY_ONLY ? 0.5 * flow_width : 0.0f;
    CurvatureEstimator cestim;
    float              min_malformation_dist = 0.55 * flow_width;

    std::vector<ExtendedPoint> points;
    points.reserve(extrusion_points.size() * (ADD_INTERSECTIONS ? 1.5 : 1));
    auto maybe_unscale = [](const P &p) { return SCALED_INPUT ? unscaled(p) : p.template cast<double>(); };

    {
        ExtendedPoint start_point{maybe_unscale(extrusion_points.front())};
        auto [distance, nearest_line, x]    = unscaled_prev_layer.signed_distance_from_lines_extra(start_point.position);
        start_point.distance                = distance + boundary_offset;
        start_point.nearest_prev_layer_line = nearest_line;
        points.push_back(start_point);
    }
    for (size_t i = 1; i < extrusion_points.size(); i++) {
        ExtendedPoint next_point{maybe_unscale(extrusion_points[i])};
        auto [distance, nearest_line, x]   = unscaled_prev_layer.signed_distance_from_lines_extra(next_point.position);
        next_point.distance                = distance + boundary_offset;
        next_point.nearest_prev_layer_line = nearest_line;

        if (ADD_INTERSECTIONS) {
            const ExtendedPoint &prev_point = points.back();
            if ((prev_point.distance < min_malformation_dist) != (next_point.distance < min_malformation_dist)) { // one in air, one not
                auto intersections = unscaled_prev_layer.template intersections_with_line<true>(L{prev_point.position, next_point.position});
                for (const auto &intersection : intersections) { points.emplace_back(intersection, boundary_offset, 1.0); }
            }

            if (PREV_LAYER_BOUNDARY_ONLY && prev_point.distance > min_malformation_dist &&
                next_point.distance > min_malformation_dist) { // both in air
                double line_len = (prev_point.position - next_point.position).norm();
                if (line_len > 3.0) {
                    double a0 = std::clamp((boundary_offset + prev_point.distance) / line_len, 0.0, 1.0);
                    double a1 = std::clamp((boundary_offset + next_point.distance) / line_len, 0.0, 1.0);
                    double t0 = std::min(a0, a1);
                    double t1 = std::max(a0, a1);

                    auto p0                         = prev_point.position + t0 * (next_point.position - prev_point.position);
                    auto [p0_dist, p0_near_l, p0_x] = unscaled_prev_layer.signed_distance_from_lines_extra(p0);
                    points.push_back(ExtendedPoint{p0, float(p0_dist + boundary_offset), p0_near_l});

                    auto p1                         = prev_point.position + t1 * (next_point.position - prev_point.position);
                    auto [p1_dist, p1_near_l, p1_x] = unscaled_prev_layer.signed_distance_from_lines_extra(p1);
                    points.push_back(ExtendedPoint{p1, float(p1_dist + boundary_offset), p1_near_l});
                }
            }
        }
        points.push_back(next_point);
    }

    for (int point_idx = 0; point_idx < int(points.size()); ++point_idx) {
        ExtendedPoint &a    = points[point_idx];
        ExtendedPoint &prev = points[point_idx > 0 ? point_idx - 1 : point_idx];

        int prev_point_idx = point_idx;
        while (prev_point_idx > 0) {
            prev_point_idx--;
            if ((a.position - points[prev_point_idx].position).squaredNorm() > EPSILON) { break; }
        }

        int next_point_index = point_idx;
        while (next_point_index < int(points.size()) - 1) {
            next_point_index++;
            if ((a.position - points[next_point_index].position).squaredNorm() > EPSILON) { break; }
        }

        if (prev_point_idx != point_idx && next_point_index != point_idx) {
            float distance = (prev.position - a.position).norm();
            float alfa     = angle(a.position - points[prev_point_idx].position, points[next_point_index].position - a.position);
            cestim.add_point(distance, alfa);
            if (CONCAVITY_RESETS_CURVATURE && alfa < 0.0) { cestim.reset(); }
        }

        if (a.distance < min_malformation_dist) { cestim.reset(); }
        a.curvature = cestim.get_curvature();
    }

    return points;
}

struct ProcessedPoint
{
    Point p;
    float speed = 1.0f;
};

class ExtrusionQualityEstimator
{
    std::unordered_map<const PrintObject *, AABBTreeLines::LinesDistancer<Linef>> prev_layer_boundaries;
    std::unordered_map<const PrintObject *, AABBTreeLines::LinesDistancer<Linef>> next_layer_boundaries;
    const PrintObject                                                            *current_object;

public:
    void set_current_object(const PrintObject *object) { current_object = object; }

    void prepare_for_new_layer(const Layer *layer)
    {
        if (layer == nullptr) return;
        const PrintObject *object     = layer->object();
        prev_layer_boundaries[object] = next_layer_boundaries[object];
        next_layer_boundaries[object] = AABBTreeLines::LinesDistancer<Linef>{to_unscaled_linesf(layer->lslices)};
    }

    std::vector<ProcessedPoint> estimate_extrusion_quality(const ExtrusionPath                &path,
                                                           const ConfigOptionPercents         &overlaps,
                                                           const ConfigOptionFloatsOrPercents &speeds,
                                                           float                               ext_perimeter_speed,
                                                           float                               original_speed)
    {
        size_t                               speed_sections_count = std::min(overlaps.values.size(), speeds.values.size());
        std::vector<std::pair<float, float>> speed_sections;
        for (size_t i = 0; i < speed_sections_count; i++) {
            float distance = path.width * (1.0 - (overlaps.get_at(i) / 100.0));
            float speed    = speeds.get_at(i).percent ? (ext_perimeter_speed * speeds.get_at(i).value / 100.0) : speeds.get_at(i).value;
            speed_sections.push_back({distance, speed});
        }
        std::sort(speed_sections.begin(), speed_sections.end(),
                  [](const std::pair<float, float> &a, const std::pair<float, float> &b) { return a.first < b.first; });

        std::vector<ExtendedPoint> extended_points =
            estimate_points_properties<true, true, true, false>(path.polyline.points, prev_layer_boundaries[current_object], path.width);

        std::vector<ProcessedPoint> processed_points;
        processed_points.reserve(extended_points.size());
        for (size_t i = 0; i < extended_points.size(); i++) {
            const ExtendedPoint &curr = extended_points[i];
            const ExtendedPoint &next = extended_points[i + 1 < extended_points.size() ? i + 1 : i];

            auto calculate_speed = [&speed_sections, &original_speed](float distance) {
                float final_speed;
                if (distance <= speed_sections.front().first) {
                    final_speed = original_speed;
                } else if (distance >= speed_sections.back().first) {
                    final_speed = speed_sections.back().second;
                } else {
                    size_t section_idx = 0;
                    while (distance > speed_sections[section_idx + 1].first) {
                        section_idx++;
                    }
                    float t = (distance - speed_sections[section_idx].first) /
                              (speed_sections[section_idx + 1].first - speed_sections[section_idx].first);
                    t           = std::clamp(t, 0.0f, 1.0f);
                    final_speed = (1.0f - t) * speed_sections[section_idx].second + t * speed_sections[section_idx + 1].second;
                }
                return final_speed;
            };

            float extrusion_speed = (calculate_speed(curr.distance) + calculate_speed(next.distance)) / 2.0f;

            processed_points.push_back({scaled(curr.position), extrusion_speed});
        }
        return processed_points;
    }
};

} // namespace Slic3r

#endif // slic3r_ExtrusionProcessor_hpp_
