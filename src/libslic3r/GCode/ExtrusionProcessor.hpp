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

#include <algorithm>
#include <cmath>
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
        if (total_distance < EPSILON) { return 0.0; }

        return total_curvature / std::max(total_distance, window_size);
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

class ExtrusionQualityEstimator
{
    AABBTreeLines::LinesDistancer<Linef> prev_layer_boundary;
    AABBTreeLines::LinesDistancer<Linef> next_layer_boundary;
    CurvatureEstimator                   cestim;

public:
    void reset_for_next_extrusion() { cestim.reset(); }

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

#if 0 // EXPORT DEBUG FILES
        Lines scaled_lines;
        for (const Linef &lf : layer_lines) { scaled_lines.push_back({Point::new_scale(lf.a), Point::new_scale(lf.b)}); }
        BoundingBox bb = get_extents(scaled_lines);

        Points inside;
        for (const Layer *layer : layers) {
            if (layer == nullptr) continue;
            auto in = to_points(to_polygons(offset_ex(layer->lslices, -scale_(0.4))));
            inside.insert(inside.end(), in.begin(), in.end());
        }

        ::Slic3r::SVG svg(debug_out_path(("path_jps" + std::to_string(rand() % 1000)).c_str()).c_str(), bb);
        svg.draw(scaled_lines, "black", scale_(0.10));
        for (Point p : inside) {
            auto [distance, line_idx, nearest_point] = next_layer_boundary.signed_distance_from_lines_extra(unscaled(p));
            if (distance > 0) {
                svg.draw(p, "red", scale_(0.2));
                svg.draw(Point::new_scale(nearest_point.x(), nearest_point.y()), "blue", scale_(0.2));
                auto li = next_layer_boundary.get_line(line_idx);
                Line ls{Point::new_scale(li.a), Point::new_scale(li.b)};
                svg.draw(ls, "yellow", scale_(0.2));


            }
        }
#endif
    }

    std::vector<float> estimate_extrusion_quality(const ExtrusionPath &path)
    {
        float flow_width              = path.width;
        float min_malformation_dist   = 0.0 * flow_width;
        float max_malformation_dist   = 1.1 * flow_width;
        float worst_malformation_dist = 0.7 * flow_width;

        std::vector<Vec2f> points;
        Polyline           pl = path.as_polyline();
        points.reserve(pl.size());
        for (const Point &p : pl) { points.push_back(unscaled(p).cast<float>()); }

        std::vector<float> point_qualities(points.size(), 1.0);
        for (int point_idx = 0; point_idx < int(points.size()); ++point_idx) {
            const Vec2f &p = points[point_idx];

            double dist_from_prev_layer = prev_layer_boundary.signed_distance_from_lines(p.cast<double>()) + flow_width * 0.5f;

            float default_dist_quality = 0.3f;
            float distance_quality     = 1.0f;
            if (dist_from_prev_layer < min_malformation_dist) {
                distance_quality = 1.0f;
                cestim.reset();
                continue;
            } else if (dist_from_prev_layer < worst_malformation_dist) {
                distance_quality = (worst_malformation_dist - dist_from_prev_layer) / (worst_malformation_dist - min_malformation_dist);
            } else if (dist_from_prev_layer < max_malformation_dist) {
                distance_quality = default_dist_quality * (1.0f - (max_malformation_dist - dist_from_prev_layer) /
                                                                      (max_malformation_dist - worst_malformation_dist));
            } else { // completely in the air. use the default value in that case
                distance_quality = default_dist_quality;
            }

            int prev_point_idx = point_idx;
            while (prev_point_idx > 0) {
                prev_point_idx--;
                if ((p - points[prev_point_idx]).squaredNorm() > EPSILON) { break; }
            }

            int next_point_index = point_idx;
            while (next_point_index < int(points.size()) - 1) {
                next_point_index++;
                if ((p - points[next_point_index]).squaredNorm() > EPSILON) { break; }
            }

            float curvature_penalty = 0.0f;
            if (prev_point_idx != point_idx && next_point_index != point_idx) {
                float distance = (p - points[prev_point_idx]).norm();
                float alfa     = angle(p - points[prev_point_idx], points[next_point_index] - p);
                cestim.add_point(distance, alfa);

                float curvature = std::abs(cestim.get_curvature());
                if (curvature > 1.0f) {
                    curvature_penalty = 1.0f;
                } else if (curvature > 0.1f) {
                    curvature_penalty = std::min(1.0, dist_from_prev_layer - min_malformation_dist) * curvature;
                }
            }

            point_qualities[point_idx] = std::clamp(distance_quality - curvature_penalty, 0.0f, 1.0f);
        }

        for (size_t point_idx = 1; point_idx < points.size(); ++point_idx) {
            point_qualities[point_idx - 1] = std::max(point_qualities[point_idx - 1], point_qualities[point_idx]);
        }

        return point_qualities;
    }
};

} // namespace Slic3r

#endif // slic3r_ExtrusionProcessor_hpp_
