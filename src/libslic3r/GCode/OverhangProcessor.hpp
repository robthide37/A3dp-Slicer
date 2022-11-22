#ifndef slic3r_OverhangProcessor_hpp_
#define slic3r_OverhangProcessor_hpp_

#include "../AABBTreeLines.hpp"
#include "../SupportSpotsGenerator.hpp"
#include "../libslic3r.h"

#include <limits>
#include <numeric>

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

inline float estimate_overhang_quality(const ExtrusionPath                  &entity,
                                       float                                 flow_width,
                                       AABBTreeLines::LinesDistancer<Linef> &prev_layer_boundary)
{
    // value of 1 is for nice straigth lines that are either in air or mostly lying on the prev layer.
    float quality = 1.0;

    float min_malformation_dist = 0.0f;
    float max_malformation_dist = 0.7 * flow_width;

    CurvatureEstimator cestim{};
    std::vector<Vec2f> points;
    Polyline           pl = entity.as_polyline();
    points.reserve(pl.size());
    for (const Point &p : pl) { points.push_back(unscaled(p).cast<float>()); }

    for (size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
        Vec2f a = points[point_idx > 0 ? point_idx - 1 : point_idx];
        Vec2f b = points[point_idx];
        Vec2f c = points[point_idx < points.size() - 1 ? point_idx + 1 : point_idx];

        const Vec2f v1         = b - a;
        const Vec2f v2         = c - b;
        float       curr_angle = angle(v1, v2);

        cestim.add_point(v1.norm(), curr_angle);
        // malformation in concave angles does not happen
        if (curr_angle < -20.0 * PI / 180.0) { cestim.reset(); }

        double dist_from_prev_layer = prev_layer_boundary.signed_distance_from_lines(b.cast<double>());

        float distance_quality  = std::abs(dist_from_prev_layer - (max_malformation_dist + min_malformation_dist) * 0.5);
        float curvature_quality = std::abs(cestim.get_curvature()) * 10.0f;
        curvature_quality       = std::max(curvature_quality, 1.0f);
        distance_quality /= curvature_quality;

        if (distance_quality < quality) { quality = 0.8 * quality + 0.2 * distance_quality; }
    }

    return quality;
}

}; // namespace Slic3r

#endif // slic3r_OverhangProcessor_hpp_
