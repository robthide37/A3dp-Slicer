// The following code for merging circles into arches originates from https://github.com/FormerLurker/ArcWelderLib

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Arc Welder: Anti-Stutter Library
//
// Compresses many G0/G1 commands into G2/G3(arc) commands where possible, ensuring the tool paths stay within the specified resolution.
// This reduces file size and the number of gcodes per second.
//
// Uses the 'Gcode Processor Library' for gcode parsing, position processing, logging, and other various functionality.
//
// Copyright(C) 2021 - Brad Hochgesang
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This program is free software : you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU Affero General Public License for more details.
//
//
// You can contact the author at the following email address: 
// FormerLurker@pm.me
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ArcWelder.hpp"

#include "../MultiPoint.hpp"
#include "../Polygon.hpp"

#include <numeric>
#include <random>
#include <boost/log/trivial.hpp>

namespace Slic3r { namespace Geometry { namespace ArcWelder {

// Area of a parallelogram of two vectors to be considered collinear.
static constexpr const double  Parallel_area_threshold = 0.0001;
static constexpr const auto    Parallel_area_threshold_scaled = int64_t(Parallel_area_threshold / sqr(SCALING_FACTOR));
// FIXME do we want to use EPSILON here?
static constexpr const double  epsilon = 0.000005;

struct Circle
{
    Point  center;
    double radius;
};

// Interpolate three points with a circle.
// Returns false if the three points are collinear or if the radius is bigger than maximum allowed radius.
//FIXME unit test!
static std::optional<Circle> try_create_circle(const Point &p1, const Point &p2, const Point &p3, const double max_radius)
{
    // Use area of triangle to judge whether three points are considered collinear.
    Vec2i64 v2 = (p2 - p1).cast<int64_t>();
    Vec2i64 v3 = (p3 - p2).cast<int64_t>();
    if (std::abs(cross2(v2, v3)) <= Parallel_area_threshold_scaled)
        return {};

    int64_t det = cross2(p2.cast<int64_t>(), p3.cast<int64_t>()) - cross2(p1.cast<int64_t>(), v3);
    if (std::abs(det) < int64_t(SCALED_EPSILON))
        return {};

    Point center = ((1. / 2.0 * double(det)) *
        (double(p1.cast<int64_t>().squaredNorm()) * perp(v3).cast<double>() + 
         double(p2.cast<int64_t>().squaredNorm()) * perp(p1 - p3).cast<double>() + 
         double(p3.cast<int64_t>().squaredNorm()) * perp(v2).cast<double>())).cast<coord_t>();
    double r = sqrt(double((center - p1).squaredNorm()));
    return r > max_radius ? std::make_optional<Circle>() : std::make_optional<Circle>({ center, r });
}

// Returns a closest point on the segment.
// Returns false if the closest point is not inside the segment, but at its boundary.
static bool foot_pt_on_segment(const Point &p1, const Point &p2, const Point &c, Point &out)
{
    Vec2i64 v21   = (p2 - p1).cast<int64_t>();
    int64_t denom = v21.squaredNorm();
    if (denom > epsilon) {
        if (double t = double((c - p1).cast<int64_t>().dot(v21)) / double(denom);
            t >= epsilon && t < 1. - epsilon) {
            out = p1 + (t * v21.cast<double>()).cast<coord_t>();
            return true;
        }
    }
    // The segment is short or the closest point is an end point.
    return false;
}

static inline bool circle_approximation_sufficient(const Circle &circle, const Points::const_iterator begin, const Points::const_iterator end, const double tolerance)
{
    // The circle was calculated from the 1st and last point of the point sequence, thus the fitting of those points does not need to be evaluated.
    assert(std::abs((*begin - circle.center).cast<double>().norm() - circle.radius) < epsilon);
    assert(std::abs((*std::prev(end) - circle.center).cast<double>().norm() - circle.radius) < epsilon);
    assert(end - begin >= 3);

    for (auto it = begin; std::next(it) != end; ++ it) {
        if (it != begin) {
            if (double distance_from_center = (*it - circle.center).cast<double>().norm();
                std::abs(distance_from_center - circle.radius) > tolerance)
                return false;
        }
        Point closest_point;
        if (foot_pt_on_segment(*it, *std::next(it), circle.center, closest_point)) {
            if (double distance_from_center = (closest_point - circle.center).cast<double>().norm();
                std::abs(distance_from_center - circle.radius) > tolerance)
                return false;
        }
    }
    return true;
}

static inline bool get_deviation_sum_squared(const Circle &circle, const Points::const_iterator begin, const Points::const_iterator end, const double tolerance, double &total_deviation)
{
    // The circle was calculated from the 1st and last point of the point sequence, thus the fitting of those points does not need to be evaluated.
    assert(std::abs((*begin - circle.center).cast<double>().norm() - circle.radius) < epsilon);
    assert(std::abs((*std::prev(end) - circle.center).cast<double>().norm() - circle.radius) < epsilon);
    assert(end - begin >= 3);

    total_deviation = 0;

    const double tolerance2 = sqr(tolerance);
    for (auto it = std::next(begin); std::next(it) != end; ++ it)
        if (double deviation2 = sqr((*it - circle.center).cast<double>().norm() - circle.radius); deviation2 > tolerance2)
            return false;
        else
            total_deviation += deviation2;

    for (auto it = begin; std::next(it) != end; ++ it) {
        Point closest_point;
        if (foot_pt_on_segment(*it, *std::next(it), circle.center, closest_point)) {
            if (double deviation2 = sqr((closest_point - circle.center).cast<double>().norm() - circle.radius); deviation2 > tolerance2)
                return false;
            else
                total_deviation += deviation2;
        }
    }

    return true;
}

static std::optional<Circle> try_create_circle(const Points::const_iterator begin, const Points::const_iterator end, const double max_radius, const double tolerance)
{
    std::optional<Circle> out;
    size_t size = end - begin;
    if (size == 3) {
        out = try_create_circle(*begin, *std::next(begin), *std::prev(end), max_radius);
        if (! circle_approximation_sufficient(*out, begin, end, tolerance))
            out.reset();
    } else {
        size_t ipivot = size / 2;
        // Take a center difference of points at the center of the path.
        //FIXME does it really help? For short arches, the linear interpolation may be 
        Point pivot = (size % 2 == 0) ? (*(begin + ipivot) + *(begin + ipivot - 1)) / 2 :
            (*(begin + ipivot - 1) + *(begin + ipivot + 1)) / 2;
        if (std::optional<Circle> circle = try_create_circle(*begin, pivot, *std::prev(end), max_radius);
            circle_approximation_sufficient(*circle, begin, end, tolerance))
            return circle;

        // Find the circle with the least deviation, if one exists.
        double least_deviation;
        double current_deviation;
        for (auto it = std::next(begin); std::next(it) != end; ++ it)
            if (std::optional<Circle> circle = try_create_circle(*begin, *it, *std::prev(end), max_radius); 
                circle && get_deviation_sum_squared(*circle, begin, end, tolerance, current_deviation)) {
                if (! out || current_deviation < least_deviation) {
                    out = circle;
                    least_deviation = current_deviation;
                }
            }
    }
    return out;
}

// ported from ArcWelderLib/ArcWelder/segmented/shape.h class "arc"
class Arc {
public:

    Arc() {}
#if 0
    Arc(Point center, double radius, Point start, Point end, Orientation dir) :
        center(center),
        radius(radius),
        start_point(start),
        end_point(end),
        direction(dir) {
        if (radius == 0.0 ||
            start_point == center ||
            end_point == center ||
            start_point == end_point) {
            is_arc = false;
            return;
        }
        is_arc = true;
    }
#endif

    Point  center;
    double radius { 0 };
    bool   is_arc { false };
    Point  start_point{ 0, 0 };
    Point  end_point{ 0, 0 };
    Orientation direction { Orientation::Unknown };

    static std::optional<Arc> try_create_arc(
        const Points::const_iterator begin,
        const Points::const_iterator end,
        double max_radius = default_scaled_max_radius,
        double tolerance = default_scaled_resolution,
        double path_tolerance_percent = default_arc_length_percent_tolerance);

    bool is_valid() const { return is_arc; }
};

static inline int sign(const int64_t i)
{
    return i > 0 ? 1 : i < 0 ? -1 : 0;
}

static inline std::optional<Arc> try_create_arc_impl(
    const Circle                &circle,
    const Points::const_iterator begin,
    const Points::const_iterator end,
    double                       path_tolerance_percent)
{
    assert(end - begin >= 3);
    // Assumption: Two successive points of a single segment span an angle smaller than PI.
    Vec2i64 vstart  = (*begin - circle.center).cast<int64_t>();
    Vec2i64 vprev   = vstart;
    int     arc_dir = 0;
    for (auto it = std::next(begin); it != end; ++ it) {
        Vec2i64 v = (*it - circle.center).cast<int64_t>();
        int     dir = sign(cross2(vprev, v));
        if (dir == 0) {
            // Ignore radial segments.
        } else if (arc_dir * dir < 0) {
            // The path turns back and overextrudes. Such path is likely invalid, but the arc interpolation should not cover it.
            return {};
        } else {
            // Success, moving in the same direction.
            arc_dir = dir;
            vprev = v;
        }
    }
    
    if (arc_dir == 0)
        // All points were radial, this should not happen.
        return {};

    Vec2i64 vend  = (*std::prev(end) - circle.center).cast<int64_t>();
    double  angle = atan2(double(cross2(vstart, vend)), double(vstart.dot(vend)));
    if (arc_dir > 0) {
        if (angle < 0)
            angle += 2. * M_PI;
    } else {
        if (angle > 0)
            angle -= 2. * M_PI;
    }

    // Check the length against the original length.
    // This can trigger simply due to the differing path lengths
    // but also could indicate that the vector calculation above
    // got wrong direction
    const double arc_length                     = std::abs(circle.radius * angle);
    const double approximate_length             = length(begin, end);
    assert(approximate_length > 0);
    const double arc_length_difference_relative = (arc_length - approximate_length) / approximate_length;
    if (std::fabs(arc_length_difference_relative) >= path_tolerance_percent)
        return {};

    Arc out;
    out.is_arc            = true;
    out.direction         = arc_dir > 0 ? Orientation::CCW : Orientation::CW;
    out.center            = circle.center;
    out.radius            = circle.radius;
    out.start_point       = *begin;
    out.end_point         = *std::prev(end);
    return std::make_optional<Arc>(out);
}

std::optional<Arc> Arc::try_create_arc(
    const Points::const_iterator begin,
    const Points::const_iterator end,
    double                       max_radius,
    double                       tolerance,
    double                       path_tolerance_percent)
{
    std::optional<Circle> circle = try_create_circle(begin, end, max_radius, tolerance);
    if (! circle)
        return {};
    return try_create_arc_impl(*circle, begin, end, path_tolerance_percent);
}

float arc_angle(const Vec2f &start_pos, const Vec2f &end_pos, Vec2f &center_pos, bool is_ccw)
{
    if ((end_pos - start_pos).squaredNorm() < sqr(1e-6)) {
        // If start equals end, full circle is considered.
        return float(2. * M_PI);
    } else {
        Vec2f v1 = start_pos - center_pos;
        Vec2f v2 = end_pos   - center_pos;
        if (! is_ccw)
            std::swap(v1, v2);
        float radian = atan2(cross2(v1, v2), v1.dot(v2));
        return radian < 0 ? float(2. * M_PI) + radian : radian;
    }
}

float arc_length(const Vec2f &start_pos, const Vec2f &end_pos, Vec2f &center_pos, bool is_ccw)
{
    return (center_pos - start_pos).norm() * arc_angle(start_pos, end_pos, center_pos, is_ccw);
}

// Reduces polyline in the <begin, end) range in place,
// returns the new end iterator.
static inline Segments::iterator douglas_peucker_in_place(Segments::iterator begin, Segments::iterator end, const double tolerance)
{
    return douglas_peucker<int64_t>(begin, end, begin, tolerance, [](const Segment &s) { return s.point; });
}

Path fit_path(const Points &src, double tolerance, double fit_circle_percent_tolerance)
{
    assert(tolerance >= 0);
    assert(fit_circle_percent_tolerance >= 0);

    Path out;
    out.reserve(src.size());
    if (tolerance <= 0 || src.size() <= 2) {
        // No simplification, just convert.
        std::transform(src.begin(), src.end(), std::back_inserter(out), [](const Point &p) -> Segment { return { p }; });
    } else if (fit_circle_percent_tolerance <= 0) {
        // Convert and simplify to a polyline.
        std::transform(src.begin(), src.end(), std::back_inserter(out), [](const Point &p) -> Segment { return { p }; });
        out.erase(douglas_peucker_in_place(out.begin(), out.end(), tolerance), out.end());
    } else {
        // Perform simplification & fitting.
        int begin_pl_idx = 0;
        for (auto begin = src.begin(); begin < src.end();) {
            // Minimum 3 points required for circle fitting.
            auto end = begin + 3;
            std::optional<Arc> arc;
            while (end <= src.end()) {
                if (std::optional<Arc> this_arc = ArcWelder::Arc::try_create_arc(
                                                        begin, end,
                                                        ArcWelder::default_scaled_max_radius,
                                                        tolerance, fit_circle_percent_tolerance);
                    this_arc) {
                    arc = this_arc;
                    ++ end;
                } else
                    break;
            }
            if (arc) {
                // If there is a trailing polyline, decimate it first before saving a new arc.
                if (out.size() - begin_pl_idx > 2)
                    out.erase(douglas_peucker_in_place(out.begin() + begin_pl_idx, out.end(), tolerance), out.end());
                // Save the end of the last circle segment, which may become the first point of a possible future polyline.
                begin_pl_idx = int(out.size());
                -- end;
                out.push_back({ arc->end_point, float(arc->direction == Orientation::CCW ? arc->radius : - arc->radius) });
            } else
                out.push_back({ arc->end_point, 0.f });
        }
    }

    return out;
}

void reverse(Path &path)
{
    if (path.size() > 1) {
        std::reverse(path.begin(), path.end());
        auto prev = path.begin();
        for (auto it = std::next(prev); it != path.end(); ++ it) {
            it->radius      = prev->radius;
            it->orientation = prev->orientation == Orientation::CCW ? Orientation::CW : Orientation::CCW;
            prev = it;
        }
        path.front().radius = 0;
    }
}

double clip_start(Path &path, const double len)
{
    reverse(path);
    double remaining = clip_end(path, len);
    reverse(path);
    // Return remaining distance to go.
    return remaining;
}

double clip_end(Path &path, double distance)
{
    while (distance > 0) {
        Segment &last = path.back();
        path.pop_back();
        if (path.empty())
            break;
        if (last.linear()) {
            // Linear segment
            Vec2d  v    = (path.back().point - last.point).cast<double>();
            double lsqr = v.squaredNorm();
            if (lsqr > sqr(distance)) {
                path.push_back({ last.point + (v * (distance / sqrt(lsqr))).cast<coord_t>(), 0.f, Orientation::CCW });
                return 0;
            }
            distance -= sqrt(lsqr);
        } else {
            // Circular segment
            float angle = arc_angle(path.back().point.cast<float>(), last.point.cast<float>(), last.radius);
            double len = std::abs(last.radius) * angle;
            if (len > distance) {
                // Rotate the segment end point in reverse towards the start point.
                path.push_back({ 
                    last.point.rotated(- angle * (distance / len), 
                        arc_center(path.back().point.cast<float>(), last.point.cast<float>(), last.radius, last.ccw()).cast<coord_t>()),
                    last.radius, last.orientation });
                return 0;
            }
            distance -= len;
        }
    }

    // Return remaining distance to go.
    assert(distance >= 0);
    return distance;
}

PathSegmentProjection point_to_path_projection(const Path &path, const Point &point, double search_radius2)
{
    assert(path.size() != 1);
    PathSegmentProjection out;
    out.distance2 = search_radius2;
    if (path.size() < 2 || path.front().point == point) {
        // First point is the closest point.
        if (path.empty()) {
        } else if (const Point p0 = path.front().point; p0 == point) {
            out.segment_id = 0;
            out.point      = p0;
            out.distance2  = 0;
        } else if (double d2 = (p0 - point).cast<double>().squaredNorm(); d2 < out.distance2) {
            out.segment_id = 0;
            out.point      = p0;
            out.distance2  = d2;
        }
    } else {
        auto  min_point_it = path.cbegin();
        Point prev         = path.front().point;
        for (auto it = path.cbegin() + 1; it != path.cend(); ++ it) {
            if (it->linear()) {
                // Linear segment
                Point proj;
                if (double d2 = line_alg::distance_to_squared(Line(prev, it->point), point, &proj); d2 < out.distance2) {
                    out.point     = proj;
                    out.distance2 = d2;
                    min_point_it  = it;
                }
            } else {
                // Circular arc
                Vec2i64 center = arc_center(prev.cast<float>(), it->point.cast<float>(), it->radius, it->ccw()).cast<int64_t>();
                // Test whether point is inside the wedge.
                Vec2i64 v1 = prev.cast<int64_t>() - center;
                Vec2i64 v2 = it->point.cast<int64_t>() - center;
                Vec2i64 vp = point.cast<int64_t>() - center;
                bool inside = it->radius > 0 ?
                    // Smaller (convex) wedge.
                    (it->ccw() ?
                        cross2(v1, vp) > 0 && cross2(vp, v2) > 0 :
                        cross2(v1, vp) < 0 && cross2(vp, v2) < 0) :
                    // Larger (concave) wedge.
                    (it->ccw() ?
                        cross2(v2, vp) < 0 || cross2(vp, v1) < 0 :
                        cross2(v2, vp) > 0 || cross2(vp, v1) > 0);
                if (inside) {
                    // Distance of the radii.
                    if (double d2 = sqr(std::abs(it->radius) - sqrt(double(v1.squaredNorm()))); d2 < out.distance2) {
                        out.distance2 = d2;
                        min_point_it  = it;
                    }
                } else {
                    // Distance to the start point.
                    if (double d2 = double((v1 - vp).squaredNorm()); d2 < out.distance2) {
                        out.point     = prev;
                        out.distance2 = d2;
                        min_point_it  = it;
                    }
                }
            }
            prev = it->point;
        }
        if (! path.back().linear()) {
            // Calculate distance to the end point.
            if (double d2 = (path.back().point - point).cast<double>().norm(); d2 < out.distance2) {
                out.point     = path.back().point;
                out.distance2 = d2;
                min_point_it  = std::prev(path.end());
            }
        }
        out.segment_id = min_point_it - path.begin();
    }

    return out;
}

std::pair<Path, Path> split_at(const Path &path, PathSegmentProjection proj, const double min_segment_length)
{
    std::pair<Path, Path> out;
    if (proj.segment_id == 0 && proj.point == path.front().point)
        out.second = path;
    else if (proj.segment_id + 1 == path.size() || (proj.segment_id + 2 == path.size() && proj.point == path.back().point))
        out.first = path;
    else {
        const Segment &start = path[proj.segment_id];
        const Segment &end   = path[proj.segment_id + 1];
        bool           split_segment = true;
        if (int64_t d = (proj.point - start.point).cast<int64_t>().squaredNorm(); d < sqr(min_segment_length)) {
            split_segment = false;
        } else if (int64_t d = (proj.point - end.point).cast<int64_t>().squaredNorm(); d < sqr(min_segment_length)) {
            ++ proj.segment_id;
            split_segment = false;
        }
        if (split_segment) {
            out.first.assign(path.begin(), path.begin() + proj.segment_id + 2);
            out.second.assign(path.begin() + proj.segment_id, path.end());
            out.first.back().point = proj.point;
            out.second.front().point = proj.point;
        } else {
            out.first.assign(path.begin(), path.begin() + proj.segment_id + 1);
            out.second.assign(path.begin() + proj.segment_id, path.end());
        }
        out.second.front().radius = 0;
    }

    return out;
}

std::pair<Path, Path> split_at(const Path &path, const Point &point, const double min_segment_length)
{
    return split_at(path, point_to_path_projection(path, point), min_segment_length);
}

} } } // namespace Slic3r::Geometry::ArcWelder
