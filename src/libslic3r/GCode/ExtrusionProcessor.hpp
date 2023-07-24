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

    std::vector<float> angles_for_curvature(points.size());
    std::vector<float> distances_for_curvature(points.size());

    for (int point_idx = 0; point_idx < int(points.size()); ++point_idx) {
        ExtendedPoint &a    = points[point_idx];
        ExtendedPoint &prev = points[point_idx > 0 ? point_idx - 1 : point_idx];

        int prev_point_idx = point_idx;
        while (prev_point_idx > 0) {
            prev_point_idx--;
            if ((a.position - points[prev_point_idx].position).squaredNorm() > EPSILON) {
                break;
            }
        }

        int next_point_index = point_idx;
        while (next_point_index < int(points.size()) - 1) {
            next_point_index++;
            if ((a.position - points[next_point_index].position).squaredNorm() > EPSILON) {
                break;
            }
        }

        distances_for_curvature[point_idx] = (prev.position - a.position).norm();
        if (prev_point_idx != point_idx && next_point_index != point_idx) {
            float alfa = angle(a.position - points[prev_point_idx].position, points[next_point_index].position - a.position);
            angles_for_curvature[point_idx] = alfa;
        } // else keep zero
    }

    for (float window_size : {3.0f, 9.0f, 16.0f}) {
        size_t tail_point      = 0;
        float  tail_window_acc = 0;
        float  tail_angle_acc  = 0;

        size_t head_point      = 0;
        float  head_window_acc = 0;
        float  head_angle_acc  = 0;

        for (int point_idx = 0; point_idx < int(points.size()); ++point_idx) {
            if (point_idx > 0) {
                tail_window_acc += distances_for_curvature[point_idx - 1];
                tail_angle_acc += angles_for_curvature[point_idx - 1];
                head_window_acc -= distances_for_curvature[point_idx - 1];
                head_angle_acc -= angles_for_curvature[point_idx - 1];
            }
            while (tail_window_acc > window_size * 0.5 && int(tail_point) < point_idx) {
                tail_window_acc -= distances_for_curvature[tail_point];
                tail_angle_acc -= angles_for_curvature[tail_point];
                tail_point++;
            }

            while (head_window_acc < window_size * 0.5 && int(head_point) < int(points.size()) - 1) {
                head_window_acc += distances_for_curvature[head_point];
                head_angle_acc += angles_for_curvature[head_point];
                head_point++;
            }

            float curvature = (tail_angle_acc + head_angle_acc) / (tail_window_acc + head_window_acc);
            if (std::abs(curvature) > std::abs(points[point_idx].curvature)) {
                points[point_idx].curvature = curvature;
            }
        }
    }

    return points;
}

ExtrusionPaths calculate_and_split_overhanging_extrusions(const ExtrusionPath                             &path,
                                                          const AABBTreeLines::LinesDistancer<Linef>      &unscaled_prev_layer,
                                                          const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines);

ExtrusionEntityCollection calculate_and_split_overhanging_extrusions(
    const ExtrusionEntityCollection                 *ecc,
    const AABBTreeLines::LinesDistancer<Linef>      &unscaled_prev_layer,
    const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines);

std::pair<float, float> calculate_overhang_speed(const ExtrusionAttributes &attributes,
                                                 const FullPrintConfig     &config,
                                                 size_t                     extruder_id,
                                                 float                      external_perim_reference_speed,
                                                 float                      default_speed);

}} // namespace Slic3r::ExtrusionProcessor

#endif // slic3r_ExtrusionProcessor_hpp_
