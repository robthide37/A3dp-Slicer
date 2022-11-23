#ifndef slic3r_ExtrusionProcessor_hpp_
#define slic3r_ExtrusionProcessor_hpp_

#include "../AABBTreeLines.hpp"
#include "../SupportSpotsGenerator.hpp"
#include "../libslic3r.h"
#include "../ExtrusionEntity.hpp"
#include "../Layer.hpp"

#include <cstddef>
#include <limits>
#include <numeric>
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
        total_curvature += angle;
        distances.push_back(distance);
        angles.push_back(angle);

        while (distances.size() > 1 && total_distance > window_size) {
            total_distance -= distances.front();
            total_curvature -= angles.front();
            distances.pop_front();
            angles.pop_front();
        }
    }

    float get_curvature() const
    {
        if (total_distance <= 0.0) { return 0.0; }

        return total_curvature / std::min(total_distance, window_size);
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
    static const size_t               sliders_count          = 4;
    SlidingWindowCurvatureAccumulator sliders[sliders_count] = {{2.0}, {4.0}, {8.0}, {16.0}};

public:
    void add_point(float distance, float angle)
    {
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

class ExtrusionQualityEstimator
{
    AABBTreeLines::LinesDistancer<Linef> prev_layer_boundary;
    AABBTreeLines::LinesDistancer<Linef> next_layer_boundary;
    CurvatureEstimator                   cestim{};

public:
    void prepare_for_new_layer(const std::vector<const Layer *> &layers)
    {
        std::vector<Linef> layer_lines;
        for (const Layer *layer : layers) {
            if (layer == nullptr) continue;
            std::vector<Linef> object_lines = to_unscaled_linesf(layer->lslices);
            layer_lines.insert(layer_lines.end(), object_lines.begin(), object_lines.end());
        }
        prev_layer_boundary = next_layer_boundary;
        next_layer_boundary = AABBTreeLines::LinesDistancer<Linef>{std::move(layer_lines)};
    }

    void reset_for_next_extrusion() { cestim.reset(); }

    std::vector<float> estimate_extrusion_quality(const ExtrusionPath &path)
    {
        float flow_width              = path.width;
        float min_malformation_dist   = 0.2 * flow_width;
        float max_malformation_dist   = 1.1 * flow_width;
        float worst_malformation_dist = 0.5 * (min_malformation_dist + max_malformation_dist);

        std::vector<Vec2f> points;
        Polyline           pl = path.as_polyline();
        points.reserve(pl.size());
        for (const Point &p : pl) { points.push_back(unscaled(p).cast<float>()); }

        std::vector<float> point_qualities(points.size(), 1.0);
        for (size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
            Vec2f b = points[point_idx];

            double dist_from_prev_layer = prev_layer_boundary.signed_distance_from_lines(b.cast<double>()) + flow_width * 0.5f;
            if (dist_from_prev_layer < min_malformation_dist) continue;

            Vec2f a = points[point_idx > 0 ? point_idx - 1 : point_idx];
            Vec2f c = points[point_idx < points.size() - 1 ? point_idx + 1 : point_idx];

            const Vec2f v1         = b - a;
            const Vec2f v2         = c - b;
            float       curr_angle = angle(v1, v2);

            cestim.add_point(v1.norm(), curr_angle);

            float distance_quality = std::min(1.0, std::abs(dist_from_prev_layer - worst_malformation_dist) /
                                                       (worst_malformation_dist - min_malformation_dist));

            // Curvature is 1 / R, where is radius of the touching sphere
            // if the radius of the touching sphere is greater than 10 mm, dont lower quality, for sharper corners do lower the quality of the point
            float curvature_value = std::abs(cestim.get_curvature()) * 10.0f;
            curvature_value       = std::max(curvature_value, 1.0f);
            distance_quality /= curvature_value;

            point_qualities[point_idx] = distance_quality;
        }

        if (points.size() > 1) { point_qualities[0] = point_qualities[1]; }

        return point_qualities;
    }
};

} // namespace Slic3r

#endif // slic3r_ExtrusionProcessor_hpp_
