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

struct ProcessedPoint
{
    Point p;
    float speed_factor = 1.0f;
};

class ExtrusionQualityEstimator
{
    std::unordered_map<const PrintObject *, AABBTreeLines::LinesDistancer<Linef>> prev_layer_boundaries;
    std::unordered_map<const PrintObject *, AABBTreeLines::LinesDistancer<Linef>> next_layer_boundaries;
    CurvatureEstimator                                                            cestim;
    const PrintObject                                                            *current_object;

public:
    void reset_for_next_extrusion() { cestim.reset(); }

    void set_current_object(const PrintObject *object) { current_object = object; }

    void prepare_for_new_layer(const Layer *layer)
    {
        if (layer == nullptr) return;
        const PrintObject *object   = layer->object();
        prev_layer_boundaries[object] = next_layer_boundaries[object];
        next_layer_boundaries[object] = AABBTreeLines::LinesDistancer<Linef>{to_unscaled_linesf(layer->lslices)};
    }

    std::vector<ProcessedPoint> estimate_extrusion_quality(const ExtrusionPath &path)
    {
        struct ExtendedPoint
        {
            ExtendedPoint(const Vec2d &pos, float dist, float quality) : position(pos), distance(dist), quality(quality) {}

            Vec2d position;
            float distance; // in multiples of flow_width
            float quality;
        };

        float flow_width             = path.width;
        float min_malformation_dist  = 0.1 * flow_width;
        float peak_malformation_dist = 0.65 * flow_width;

        const Points              &original_points = path.polyline.points;
        std::vector<ExtendedPoint> points;

        const auto& prev_layer_boundary = prev_layer_boundaries[current_object];

        float distance = prev_layer_boundary.signed_distance_from_lines(unscaled(original_points[0])) + 0.5 * flow_width;
        points.push_back({unscaled(original_points[0]), distance, 1.0f});
        for (size_t i = 1; i < original_points.size(); i++) {
            Vec2d next_point_pos   = unscaled(original_points[i]);
            float distance_of_next = prev_layer_boundary.signed_distance_from_lines(next_point_pos) + 0.5 * flow_width;
            if ((points.back().distance < min_malformation_dist) != (distance_of_next < min_malformation_dist)) { // one in air, one not
                auto intersections = prev_layer_boundary.intersections_with_line<true>({points.back().position, next_point_pos});
                for (const auto &intersection : intersections) { points.push_back({intersection, 0.5f * flow_width, 1.0}); }
                points.push_back({next_point_pos, distance_of_next, 1.0});
            }

            if (points.back().distance > min_malformation_dist && distance_of_next > min_malformation_dist) { // both in air
                double line_len = (points.back().position - next_point_pos).norm();
                if (line_len > 3.0) {
                    double a0 = std::clamp((0.5 * flow_width + points.back().distance) / line_len, 0.0, 1.0);
                    double a1 = std::clamp((0.5 * flow_width + distance_of_next) / line_len, 0.0, 1.0);
                    double t0 = std::min(a0, a1);
                    double t1 = std::max(a0, a1);

                    auto p0      = points.back().position + t0 * (next_point_pos - points.back().position);
                    auto p0_dist = prev_layer_boundary.signed_distance_from_lines(p0) + 0.5 * flow_width;
                    points.push_back({p0, float(p0_dist), 1.0});
                    auto p1      = points.back().position + t1 * (next_point_pos - points.back().position);
                    auto p1_dist = prev_layer_boundary.signed_distance_from_lines(p1) + 0.5 * flow_width;
                    points.push_back({p1, float(p1_dist), 1.0});
                }
            }

            points.push_back({next_point_pos, distance_of_next, 1.0});
        }

        for (int point_idx = 0; point_idx < int(points.size()) - 1; ++point_idx) {
            ExtendedPoint &a = points[point_idx];
            ExtendedPoint &b = points[point_idx + 1];

            float distance = std::min(a.distance, b.distance);

            int prev_point_idx = point_idx;
            while (prev_point_idx > 0) {
                prev_point_idx--;
                if ((b.position - points[prev_point_idx].position).squaredNorm() > EPSILON) { break; }
            }

            int next_point_index = point_idx;
            while (next_point_index < int(points.size()) - 1) {
                next_point_index++;
                if ((b.position - points[next_point_index].position).squaredNorm() > EPSILON) { break; }
            }

            if (prev_point_idx != point_idx && next_point_index != point_idx) {
                float distance = (b.position - a.position).norm();
                float alfa     = angle(b.position - points[prev_point_idx].position, points[next_point_index].position - b.position);
                cestim.add_point(distance, alfa);
            }

            if (distance < min_malformation_dist) {
                a.quality = 1.0;
                cestim.reset();
            } else if (distance < peak_malformation_dist) {
                a.quality               = 1.0 - (distance - min_malformation_dist) / (peak_malformation_dist - min_malformation_dist);
                float curvature_penalty = 0.0f;
                float curvature         = std::abs(cestim.get_curvature());
                if (curvature > 1.0f) {
                    curvature_penalty = 1.0f;
                } else if (curvature > 0.1f) {
                    curvature_penalty = fmin(1.0, (distance - min_malformation_dist) / flow_width) * curvature;
                }
                a.quality -= curvature_penalty;

            } else {
                a.quality = 0.0f;
            }
        }

        if (points.size() > 3) {
            points[points.size() - 2].quality = points[points.size()-3].quality;
        }

        std::vector<ProcessedPoint> result;
        result.reserve(points.size());
        for (const ExtendedPoint &p : points) { result.push_back({Point::new_scale(p.position), std::clamp(p.quality, 0.0f, 1.0f)}); }

        return result;
    }
};

} // namespace Slic3r

#endif // slic3r_ExtrusionProcessor_hpp_
