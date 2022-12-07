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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <unordered_map>
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
    static const size_t               sliders_count          = 3;
    SlidingWindowCurvatureAccumulator sliders[sliders_count] = {{2.0}, {4.0}, {8.0}};

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
    ExtendedPoint(Vec2d position, float distance = 0.0, size_t nearest_prev_layer_line = size_t(-1), float curvature = 0.0, float quality = 1.0)
        : position(position), distance(distance), nearest_prev_layer_line(nearest_prev_layer_line), curvature(curvature), quality(quality)
    {}

    Vec2d  position;
    float  distance;
    size_t nearest_prev_layer_line;
    float  curvature;
    float  quality;
};

template<bool SCALED_INPUT, bool ADD_INTERSECTIONS, bool PREV_LAYER_BOUNDARY_ONLY, bool CONCAVITY_RESETS_CURVATURE, typename P, typename L>
std::vector<ExtendedPoint> estimate_points_properties(const std::vector<P>                   &extrusion_points,
                                                      const AABBTreeLines::LinesDistancer<L> &unscaled_prev_layer,
                                                      float                                   flow_width)
{
    if (extrusion_points.empty()) return {};
    float              boundary_offset = PREV_LAYER_BOUNDARY_ONLY ? 0.5 * flow_width : 0.0f;
    CurvatureEstimator cestim;
    float              min_malformation_dist  = 0.2 * flow_width;
    float              peak_malformation_dist = 0.6 * flow_width;

    std::vector<ExtendedPoint> points;
    points.reserve(extrusion_points.size() * (ADD_INTERSECTIONS ? 1.5 : 1));
    auto maybe_unscale = [](const P &p) { return SCALED_INPUT ? unscaled(p) : p.template cast<double>(); };

    {
        ExtendedPoint start_point{maybe_unscale(extrusion_points.begin())};
        auto [distance, nearest_line, x]    = unscaled_prev_layer.signed_distance_from_lines_extra(start_point.position) + boundary_offset;
        start_point.distance                = distance;
        start_point.nearest_prev_layer_line = nearest_line;
        points.push_back(start_point);
    }
    for (size_t i = 1; i < extrusion_points.size(); i++) {
        ExtendedPoint next_point{maybe_unscale(extrusion_points[i])};
        auto [distance, nearest_line, x]   = unscaled_prev_layer.signed_distance_from_lines_extra(next_point.position) + boundary_offset;
        next_point.distance                = distance;
        next_point.nearest_prev_layer_line = nearest_line;

        if (ADD_INTERSECTIONS) {
            const ExtendedPoint &prev_point = points.back();
            if ((prev_point.distance < min_malformation_dist) != (next_point.distance < min_malformation_dist)) { // one in air, one not
                auto intersections = unscaled_prev_layer.intersections_with_line<true>(L{prev_point.position, next_point.position});
                for (const auto &intersection : intersections) { points.push_back({intersection, boundary_offset, 1.0}); }
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
                    auto [p0_dist, p0_near_l, p0_x] = unscaled_prev_layer.signed_distance_from_lines(p0) + boundary_offset;
                    points.push_back(ExtendedPoint{p0, p0_dist, p0_near_l});

                    auto p1                         = prev_point.position + t1 * (next_point.position - prev_point.position);
                    auto [p1_dist, p1_near_l, p1_x] = unscaled_prev_layer.signed_distance_from_lines(p1) + boundary_offset;
                    points.push_back(ExtendedPoint{p1, p1_dist, p1_near_l});
                }
            }
        }
        points.push_back(next_point);
    }

    for (int point_idx = 0; point_idx < points.size(); ++point_idx) {
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
            if (alfa > 0.95 * 0.5 * PI) {
                alfa = 0; // Ignore very sharp corners.. The curling problem happens mostly on rounded surfaces, not sudden sharp turns
            }
            cestim.add_point(distance, alfa);
            if (CONCAVITY_RESETS_CURVATURE && alfa < 0.0) { cestim.reset(); }
        }

        if (a.distance < min_malformation_dist) {
            a.quality = 1.0;
            cestim.reset();
        } else {
            float distance_quality = std::min(1.0f, std::abs(a.distance - peak_malformation_dist) /
                                                        (peak_malformation_dist - min_malformation_dist));
            distance_quality       = distance_quality * distance_quality;

            float curvature_penalty = 0.0f;
            a.curvature             = cestim.get_curvature();
            float curvature         = std::abs(a.curvature);
            if (curvature > 1.0f) {
                curvature_penalty = 1.0f;
            } else if (curvature > 0.1f) {
                curvature_penalty = sqrt(1.0 - distance_quality) * curvature;
            }

            a.quality = std::clamp(distance_quality - curvature_penalty, 0.0f, 1.0f);
        }
    }

    return points;
}

struct ProcessedPoint
{
    Point p;
    float speed_factor = 1.0f;
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

    std::vector<ProcessedPoint> estimate_extrusion_quality(const ExtrusionPath &path)
    {
        std::vector<ExtendedPoint> extended_points =
            estimate_points_properties<true, true, true, false>(path.polyline.points, prev_layer_boundaries[current_object], path.width);
        
        std::vector<ProcessedPoint> processed_points;
        processed_points.reserve(extended_points.size());
        for (size_t i = 0; i < extended_points.size(); i++) {
            Point position = scaled(extended_points[i].position);
            float speed_factor = std::min(extended_points[i].quality, extended_points[i+1].quality);
            processed_points.push_back({position, speed_factor});
        }
        return processed_points;
    }
};

} // namespace Slic3r

#endif // slic3r_ExtrusionProcessor_hpp_
