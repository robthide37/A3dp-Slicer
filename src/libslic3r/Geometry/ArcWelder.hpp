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

#ifndef slic3r_Geometry_ArcWelder_hpp_
#define slic3r_Geometry_ArcWelder_hpp_

#include <optional>

#include "../Point.hpp"

namespace Slic3r { namespace Geometry { namespace ArcWelder {

// Calculate center point of an arc given two points and a radius.
// positive radius: take shorter arc
// negative radius: take longer arc
// radius must NOT be zero!
template<typename Derived, typename Derived2>
inline Eigen::Matrix<typename Derived::Scalar, 2, 1, Eigen::DontAlign> arc_center(
    const Eigen::MatrixBase<Derived>   &start_pos,
    const Eigen::MatrixBase<Derived2>  &end_pos, 
    const typename Derived::Scalar      radius, 
    const bool                          is_ccw)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "arc_center(): first parameter is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_center(): second parameter is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value, "arc_center(): Both vectors must be of the same type.");
    assert(radius != 0);
    using Float = typename Derived::Scalar;
    using Vector = Eigen::Matrix<Float, 2, 1, Eigen::DontAlign>;
    auto  v  = end_pos - start_pos;
    Float q2 = v.squaredNorm();
    assert(q2 > 0);
    Float t = sqrt(sqr(radius) / q2 - Float(.25f));
    auto mid = Float(0.5) * (start_pos + end_pos);
    Vector vp{ -v.y() * t, v.x() * t };
    return (radius > Float(0)) == is_ccw ? (mid + vp).eval() : (mid - vp).eval();
}

// Calculate angle of an arc given two points and a radius.
// Returned angle is in the range <0, 2 PI)
// positive radius: take shorter arc
// negative radius: take longer arc
// radius must NOT be zero!
template<typename Derived, typename Derived2>
inline typename Derived::Scalar arc_angle(
    const Eigen::MatrixBase<Derived>   &start_pos,
    const Eigen::MatrixBase<Derived2>  &end_pos, 
    const typename Derived::Scalar      radius)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "arc_angle(): first parameter is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_angle(): second parameter is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value, "arc_angle(): Both vectors must be of the same type.");
    assert(radius != 0);
    using Float = typename Derived::Scalar;
    Float a     = Float(2.) * asin(Float(0.5) * (end_pos - start_pos).norm() / radius);
    return radius > Float(0) ? a : Float(2. * M_PI) + a;
}

// Calculate positive length of an arc given two points and a radius.
// positive radius: take shorter arc
// negative radius: take longer arc
// radius must NOT be zero!
template<typename Derived, typename Derived2>
inline typename Derived::Scalar arc_length(
    const Eigen::MatrixBase<Derived>   &start_pos,
    const Eigen::MatrixBase<Derived2>  &end_pos,
    const typename Derived::Scalar      radius)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "arc_length(): first parameter is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_length(): second parameter is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value, "arc_length(): Both vectors must be of the same type.");
    assert(radius != 0);
    return arc_angle(start_pos, end_pos, radius) * std::abs(radius);
}

// 1.2m diameter, maximum given by coord_t
static constexpr const double default_scaled_max_radius = scaled<double>(600.);
// 0.05mm
static constexpr const double default_scaled_resolution = scaled<double>(0.05);
// 5 percent
static constexpr const double default_arc_length_percent_tolerance = 0.05;

enum class Orientation : unsigned char {
    Unknown,
    CCW,
    CW,
};

// Single segment of a smooth path.
struct Segment
{
    // End point of a linear or circular segment.
    // Start point is provided by the preceding segment.
    Point       point;
    // Radius of a circular segment. Positive - take the shorter arc. Negative - take the longer arc. Zero - linear segment.
    float       radius{ 0.f };
    // CCW or CW. Ignored for zero radius (linear segment).
    Orientation orientation{ Orientation::CCW };

    bool    linear() const { return radius == 0; }
    bool    ccw() const { return orientation == Orientation::CCW; }
    bool    cw() const { return orientation == Orientation::CW; }
};

using Segments = std::vector<Segment>;
using Path = Segments;

// Interpolate polyline path with a sequence of linear / circular segments given the interpolation tolerance.
// Only convert to polyline if zero tolerance.
// Convert to polyline and decimate polyline if zero fit_circle_percent_tolerance.
Path fit_path(const Points &points, double tolerance, double fit_circle_percent_tolerance);
inline Path fit_polyline(const Points &points, double tolerance) { return fit_path(points, tolerance, 0.); }

inline double segment_length(const Segment &start, const Segment &end)
{
    return end.linear() ?
        (end.point - start.point).cast<double>().norm() :
        arc_length(start.point.cast<float>(), end.point.cast<float>(), end.radius);
}

// Estimate minimum path length of a segment cheaply without having to calculate center of an arc and it arc length.
// Used for caching a smooth path chunk that is certainly longer than a threshold.
inline int64_t estimate_min_segment_length(const Segment &start, const Segment &end)
{
    if (end.linear() || end.radius > 0) {
        // Linear segment or convex wedge, take the larger X or Y component.
        Point v = (end.point - start.point).cwiseAbs();
        return std::max(v.x(), v.y());
    } else {
        // Arc with angle > PI.
        // Returns estimate of PI * r
        return - 3 * int64_t(end.radius);
    }
}

// Estimate minimum path length cheaply without having to calculate center of an arc and it arc length.
// Used for caching a smooth path chunk that is certainly longer than a threshold.
inline int64_t estimate_path_length(const Path &path)
{
    int64_t len = 0;
    for (size_t i = 1; i < path.size(); ++ i)
        len += Geometry::ArcWelder::estimate_min_segment_length(path[i - 1], path[i]);
    return len;
}

void reverse(Path &path);

// Clip start / end of a smooth path by len.
// If path is shorter than len, remaining path length to trim will be returned.
double clip_start(Path &path, const double len);
double clip_end(Path &path, const double len);

struct PathSegmentProjection
{
    // Start segment of a projection on the path.
    size_t segment_id { std::numeric_limits<size_t>::max() };
    Point  point      { 0, 0 };
    // Square of a distance of the projection.
    double distance2  { std::numeric_limits<double>::max() };

    bool   valid() const { return this->segment_id != std::numeric_limits<size_t>::max(); }
};
// Returns closest segment and a parameter along the closest segment of a path to a point.
PathSegmentProjection point_to_path_projection(const Path &path, const Point &point, double search_radius2 = std::numeric_limits<double>::max());
// Split a path into two paths at a segment point. Snap to an existing point if the projection of "point is closer than min_segment_length.
std::pair<Path, Path> split_at(const Path &path, PathSegmentProjection proj, const double min_segment_length);
// Split a path into two paths at a point closest to "point". Snap to an existing point if the projection of "point is closer than min_segment_length.
std::pair<Path, Path> split_at(const Path &path, const Point &point, const double min_segment_length);

} } } // namespace Slic3r::Geometry::ArcWelder

#endif // slic3r_Geometry_ArcWelder_hpp_
