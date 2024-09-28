///|/ Copyright (c) Prusa Research 2022 - 2023 Pavel Mikuš @Godrak, Vojtěch Bubník @bubnikv
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
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
#include "../Line.hpp"
#include "../Exception.hpp"
#include "../PrintConfig.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <ostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Slic3r { namespace ExtrusionProcessor {

struct ExtendedPoint
{
    Vec2d position;
    float distance;
    float curvature;
};

template<bool SCALED_INPUT, bool ADD_INTERSECTIONS, bool PREV_LAYER_BOUNDARY_OFFSET, bool SIGNED_DISTANCE, typename POINTS, typename L>
std::vector<ExtendedPoint> estimate_points_properties(const POINTS                           &input_points,
                                                      const AABBTreeLines::LinesDistancer<L> &unscaled_prev_layer,
                                                      float                                   flow_width,
                                                      float                                   max_line_length = -1.0f)
{
    bool   looped     = input_points.front() == input_points.back();
    std::function<size_t(size_t,size_t)> get_prev_index = [](size_t idx, size_t count) {
        if (idx > 0) {
            return idx - 1;
        } else
            return idx;
    };
    if (looped) {
        get_prev_index = [](size_t idx, size_t count) {
            if (idx == 0)
                idx = count;
            return --idx;
        };
    };
    std::function<size_t(size_t,size_t)> get_next_index = [](size_t idx, size_t size) {
        if (idx + 1 < size) {
            return idx + 1;
        } else
            return idx;
    };
    if (looped) {
        get_next_index = [](size_t idx, size_t count) {
            if (++idx == count)
                idx = 0;
            return idx;
        };
    };

    using P = typename POINTS::value_type;

    using AABBScalar = typename AABBTreeLines::LinesDistancer<L>::Scalar;
    if (input_points.empty())
        return {};
    float boundary_offset = PREV_LAYER_BOUNDARY_OFFSET ? 0.5 * flow_width : 0.0f;
    auto  maybe_unscale   = [](const P &p) { return SCALED_INPUT ? unscaled(p) : p.template cast<double>(); };

    std::vector<ExtendedPoint> points;
    points.reserve(input_points.size() * (ADD_INTERSECTIONS ? 1.5 : 1));

    {
        ExtendedPoint start_point{maybe_unscale(input_points.front())};
        auto [distance, nearest_line,
              x] = unscaled_prev_layer.template distance_from_lines_extra<SIGNED_DISTANCE>(start_point.position.cast<AABBScalar>());
        start_point.distance = distance + boundary_offset;
        points.push_back(start_point);
    }
    for (size_t i = 1; i < input_points.size(); i++) {
        ExtendedPoint next_point{maybe_unscale(input_points[i])};
        auto [distance, nearest_line,
              x] = unscaled_prev_layer.template distance_from_lines_extra<SIGNED_DISTANCE>(next_point.position.cast<AABBScalar>());
        next_point.distance = distance + boundary_offset;

        if (ADD_INTERSECTIONS &&
            ((points.back().distance > boundary_offset + EPSILON) != (next_point.distance > boundary_offset + EPSILON))) {
            const ExtendedPoint &prev_point    = points.back();
            auto                 intersections = unscaled_prev_layer.template intersections_with_line<true>(
                L{prev_point.position.cast<AABBScalar>(), next_point.position.cast<AABBScalar>()});
            for (const auto &intersection : intersections) {
                ExtendedPoint p{};
                p.position = intersection.first.template cast<double>();
                p.distance = boundary_offset;
                points.push_back(p);
            }
        }
        points.push_back(next_point);
    }

    if (PREV_LAYER_BOUNDARY_OFFSET && ADD_INTERSECTIONS) {
        std::vector<ExtendedPoint> new_points;
        new_points.reserve(points.size() * 2);
        new_points.push_back(points.front());
        for (int point_idx = 0; point_idx < int(points.size()) - 1; ++point_idx) {
            const ExtendedPoint &curr = points[point_idx];
            const ExtendedPoint &next = points[point_idx + 1];

            if ((curr.distance > -boundary_offset && curr.distance < boundary_offset + 2.0f) ||
                (next.distance > -boundary_offset && next.distance < boundary_offset + 2.0f)) {
                double line_len = (next.position - curr.position).norm();
                if (line_len > 4.0f) {
                    double a0 = std::clamp((curr.distance + 3 * boundary_offset) / line_len, 0.0, 1.0);
                    double a1 = std::clamp(1.0f - (next.distance + 3 * boundary_offset) / line_len, 0.0, 1.0);
                    double t0 = std::min(a0, a1);
                    double t1 = std::max(a0, a1);

                    if (t0 < 1.0) {
                        auto p0     = curr.position + t0 * (next.position - curr.position);
                        auto [p0_dist, p0_near_l,
                              p0_x] = unscaled_prev_layer.template distance_from_lines_extra<SIGNED_DISTANCE>(p0.cast<AABBScalar>());
                        ExtendedPoint new_p{};
                        new_p.position = p0;
                        new_p.distance = float(p0_dist + boundary_offset);
                        new_points.push_back(new_p);
                    }
                    if (t1 > 0.0) {
                        auto p1     = curr.position + t1 * (next.position - curr.position);
                        auto [p1_dist, p1_near_l,
                              p1_x] = unscaled_prev_layer.template distance_from_lines_extra<SIGNED_DISTANCE>(p1.cast<AABBScalar>());
                        ExtendedPoint new_p{};
                        new_p.position = p1;
                        new_p.distance = float(p1_dist + boundary_offset);
                        new_points.push_back(new_p);
                    }
                }
            }
            new_points.push_back(next);
        }
        points = std::move(new_points);
    }

    if (max_line_length > 0) {
        std::vector<ExtendedPoint> new_points;
        new_points.reserve(points.size() * 2);
        {
            for (size_t i = 0; i + 1 < points.size(); i++) {
                const ExtendedPoint &curr = points[i];
                const ExtendedPoint &next = points[i + 1];
                new_points.push_back(curr);
                double len             = (next.position - curr.position).squaredNorm();
                double t               = sqrt((max_line_length * max_line_length) / len);
                size_t new_point_count = 1.0 / t;
                for (size_t j = 1; j < new_point_count + 1; j++) {
                    Vec2d pos  = curr.position * (1.0 - j * t) + next.position * (j * t);
                    auto [p_dist, p_near_l,
                          p_x] = unscaled_prev_layer.template distance_from_lines_extra<SIGNED_DISTANCE>(pos.cast<AABBScalar>());
                    ExtendedPoint new_p{};
                    new_p.position = pos;
                    new_p.distance = float(p_dist + boundary_offset);
                    new_points.push_back(new_p);
                }
            }
            new_points.push_back(points.back());
        }
        points = std::move(new_points);
    }

    float accumulated_distance = 0;
    std::vector<float> distances_for_curvature(points.size());
    for (size_t point_idx = 0; point_idx < points.size(); ++point_idx) {
        const ExtendedPoint &a = points[point_idx];
        const ExtendedPoint &b = points[get_prev_index(point_idx, points.size())];

        distances_for_curvature[point_idx] = (b.position - a.position).norm();
        accumulated_distance += distances_for_curvature[point_idx];
    }

    if (accumulated_distance > EPSILON)
        for (float window_size : {3.0f, 9.0f, 16.0f}) {
            for (int point_idx = 0; point_idx < int(points.size()); ++point_idx) {
                ExtendedPoint &current = points[point_idx];

                Vec2d back_position = current.position;
                {
                    size_t back_point_index = point_idx;
                    float  dist_backwards   = 0;
                    while (dist_backwards < window_size * 0.5 && back_point_index != get_prev_index(back_point_index, points.size())) {
                        float line_dist = distances_for_curvature[get_prev_index(back_point_index, points.size())];
                        if (dist_backwards + line_dist > window_size * 0.5) {
                            back_position = points[back_point_index].position +
                                            (window_size * 0.5 - dist_backwards) *
                                                (points[get_prev_index(back_point_index, points.size())].position -
                                                 points[back_point_index].position)
                                                    .normalized();
                            dist_backwards += window_size * 0.5 - dist_backwards + EPSILON;
                        } else {
                            dist_backwards += line_dist;
                            back_point_index = get_prev_index(back_point_index, points.size());
                        }
                    }
                }

                Vec2d front_position = current.position;
                {
                    size_t front_point_index = point_idx;
                    float  dist_forwards     = 0;
                    while (dist_forwards < window_size * 0.5 && front_point_index != get_next_index(front_point_index, points.size())) {
                        float line_dist = distances_for_curvature[front_point_index];
                        if (dist_forwards + line_dist > window_size * 0.5) {
                            front_position = points[front_point_index].position +
                                             (window_size * 0.5 - dist_forwards) *
                                                 (points[get_next_index(front_point_index, points.size())].position -
                                                  points[front_point_index].position)
                                                     .normalized();
                            dist_forwards += window_size * 0.5 - dist_forwards + EPSILON;
                        } else {
                            dist_forwards += line_dist;
                            front_point_index = get_next_index(front_point_index, points.size());
                        }
                    }
                }
                float new_curvature = angle_ccw(current.position - back_position, front_position - current.position) / window_size;
                if (abs(current.curvature) < abs(new_curvature)) {
                    current.curvature = new_curvature;
                }
            }
        }

    return points;
}

//ExtrusionPaths calculate_and_split_overhanging_extrusions(const ExtrusionPath                             &path,
//                                                          const AABBTreeLines::LinesDistancer<Linef>      &unscaled_prev_layer,
//                                                          const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines);

ExtrusionEntityCollection calculate_and_split_overhanging_extrusions(
    const ExtrusionEntityCollection                 *ecc,
    const AABBTreeLines::LinesDistancer<Linef>      &unscaled_prev_layer,
    const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines,
    const double max_width);

std::pair<float, float> calculate_overhang_speed(const ExtrusionAttributes &attributes,
                                                 const FullPrintConfig     &config,
                                                 size_t                     extruder_id);

}} // namespace Slic3r::ExtrusionProcessor

#endif // slic3r_ExtrusionProcessor_hpp_
