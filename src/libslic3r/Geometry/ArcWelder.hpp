#ifndef slic3r_Geometry_ArcWelder_hpp_
#define slic3r_Geometry_ArcWelder_hpp_

#include <optional>

#include "../Point.hpp"

namespace Slic3r { namespace Geometry { namespace ArcWelder {

// Calculate center point of an arc given two points and a radius.
// positive radius: take shorter arc
// negative radius: take longer arc
// radius must NOT be zero!
template<typename Derived, typename Derived2, typename Float>
inline Eigen::Matrix<Float, 2, 1, Eigen::DontAlign> arc_center(
    const Eigen::MatrixBase<Derived>   &start_pos,
    const Eigen::MatrixBase<Derived2>  &end_pos, 
    const Float                         radius,
    const bool                          is_ccw)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "arc_center(): first parameter is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_center(): second parameter is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value, "arc_center(): Both vectors must be of the same type.");
    static_assert(std::is_same<typename Derived::Scalar, Float>::value, "arc_center(): Radius must be of the same type as the vectors.");
    assert(radius != 0);
    using Vector = Eigen::Matrix<Float, 2, 1, Eigen::DontAlign>;
    auto  v  = end_pos - start_pos;
    Float q2 = v.squaredNorm();
    assert(q2 > 0);
    Float t2 = sqr(radius) / q2 - Float(.25f);
    // If the start_pos and end_pos are nearly antipodal, t2 may become slightly negative.
    // In that case return a centroid of start_point & end_point.
    Float t = t2 > 0 ? sqrt(t2) : Float(0);
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
    Float a = Float(0.5) * (end_pos - start_pos).norm() / radius;
    return radius > Float(0) ?
        // acute angle:
        (a > Float( 1.) ? Float(M_PI) : Float(2.) * std::asin(a)) :
        // obtuse angle:
        (a < Float(-1.) ? Float(M_PI) : Float(2. * M_PI) + Float(2.) * std::asin(a));
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

// Calculate positive length of an arc given two points, center and orientation.
template<typename Derived, typename Derived2, typename Derived3>
inline typename Derived::Scalar arc_length(
    const Eigen::MatrixBase<Derived>   &start_pos,
    const Eigen::MatrixBase<Derived2>  &end_pos,
    const Eigen::MatrixBase<Derived3>  &center_pos,
    const bool                          ccw)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "arc_length(): first parameter is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_length(): second parameter is not a 2D vector");
    static_assert(Derived3::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "arc_length(): third parameter is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value &&
                  std::is_same<typename Derived::Scalar, typename Derived3::Scalar>::value, "arc_length(): All third points must be of the same type.");
    using Float = typename Derived::Scalar;
    auto  vstart = start_pos - center_pos;
    auto  vend   = end_pos - center_pos;
    Float radius = vstart.norm();
    Float angle  = atan2(double(cross2(vstart, vend)), double(vstart.dot(vend)));
    if (! ccw)
        angle *= Float(-1.);
    if (angle < 0)
        angle += Float(2. * M_PI);
    assert(angle >= Float(0.) && angle < Float(2. * M_PI + EPSILON));
    return angle * radius;
}

// Test whether a point is inside a wedge of an arc.
template<typename Derived, typename Derived2, typename Derived3>
inline bool inside_arc_wedge_vectors(
    const Eigen::MatrixBase<Derived>    &start_vec,
    const Eigen::MatrixBase<Derived2>   &end_vec,
    const bool                           shorter_arc,
    const bool                           ccw,
    const Eigen::MatrixBase<Derived3>   &query_vec)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "inside_arc_wedge_vectors(): start_vec is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "inside_arc_wedge_vectors(): end_vec is not a 2D vector");
    static_assert(Derived3::IsVectorAtCompileTime && int(Derived3::SizeAtCompileTime) == 2, "inside_arc_wedge_vectors(): query_vec is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value &&
                  std::is_same<typename Derived::Scalar, typename Derived3::Scalar>::value, "inside_arc_wedge_vectors(): All vectors must be of the same type.");
    return shorter_arc ?
        // Smaller (convex) wedge.
        (ccw ?
            cross2(start_vec, query_vec) > 0 && cross2(query_vec, end_vec) > 0 :
            cross2(start_vec, query_vec) < 0 && cross2(query_vec, end_vec) < 0) :
        // Larger (concave) wedge.
        (ccw ?
            cross2(end_vec, query_vec) < 0 || cross2(query_vec, start_vec) < 0 :
            cross2(end_vec, query_vec) > 0 || cross2(query_vec, start_vec) > 0);
}

template<typename Derived, typename Derived2, typename Derived3, typename Derived4>
inline bool inside_arc_wedge(
    const Eigen::MatrixBase<Derived>    &start_pt,
    const Eigen::MatrixBase<Derived2>   &end_pt,
    const Eigen::MatrixBase<Derived3>   &center_pt,
    const bool                           shorter_arc,
    const bool                           ccw,
    const Eigen::MatrixBase<Derived4>   &query_pt)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "inside_arc_wedge(): start_pt is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "inside_arc_wedge(): end_pt is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived3::SizeAtCompileTime) == 2, "inside_arc_wedge(): center_pt is not a 2D vector");
    static_assert(Derived3::IsVectorAtCompileTime && int(Derived4::SizeAtCompileTime) == 2, "inside_arc_wedge(): query_pt is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value &&
        std::is_same<typename Derived::Scalar, typename Derived3::Scalar>::value &&
        std::is_same<typename Derived::Scalar, typename Derived4::Scalar>::value, "inside_arc_wedge(): All vectors must be of the same type.");
    return inside_arc_wedge_vectors(start_pt - center_pt, end_pt - center_pt, shorter_arc, ccw, query_pt - center_pt);
}

template<typename Derived, typename Derived2, typename Derived3, typename Float>
inline bool inside_arc_wedge(
    const Eigen::MatrixBase<Derived>    &start_pt,
    const Eigen::MatrixBase<Derived2>   &end_pt,
    const Float                          radius,
    const bool                           ccw,
    const Eigen::MatrixBase<Derived3>   &query_pt)
{
    static_assert(Derived::IsVectorAtCompileTime && int(Derived::SizeAtCompileTime) == 2, "inside_arc_wedge(): start_pt is not a 2D vector");
    static_assert(Derived2::IsVectorAtCompileTime && int(Derived2::SizeAtCompileTime) == 2, "inside_arc_wedge(): end_pt is not a 2D vector");
    static_assert(Derived3::IsVectorAtCompileTime && int(Derived3::SizeAtCompileTime) == 2, "inside_arc_wedge(): query_pt is not a 2D vector");
    static_assert(std::is_same<typename Derived::Scalar, typename Derived2::Scalar>::value &&
        std::is_same<typename Derived::Scalar, typename Derived3::Scalar>::value &&
        std::is_same<typename Derived::Scalar, Float>::value, "inside_arc_wedge(): All vectors + radius must be of the same type.");
    return inside_arc_wedge(start_pt, end_pt,
        arc_center(start_pt, end_pt, radius, ccw), 
        radius > 0, ccw, query_pt);
}

// Return number of linear segments necessary to interpolate arc of a given positive radius and positive angle to satisfy
// maximum deviation of an interpolating polyline from an analytic arc.
template<typename FloatType>
size_t arc_discretization_steps(const FloatType radius, const FloatType angle, const FloatType deviation)
{
    assert(radius > 0);
    assert(angle > 0);
    assert(angle <= FloatType(2. * M_PI));
    assert(deviation > 0);

    FloatType d = radius - deviation;
    return d < EPSILON ?
        // Radius smaller than deviation.
        (   // Acute angle: a single segment interpolates the arc with sufficient accuracy.
            angle < M_PI || 
            // Obtuse angle: Test whether the furthest point (center) of an arc is closer than deviation to the center of a line segment.
            radius * (FloatType(1.) + cos(M_PI - FloatType(.5) * angle)) < deviation ?
            // Single segment is sufficient
            1 :
            // Two segments are necessary, the middle point is at the center of the arc.
            2) :
        size_t(ceil(angle / (2. * acos(d / radius))));
}

// Discretize arc given the radius, orientation and maximum deviation from the arc.
// Returned polygon starts with p1, ends with p2 and it is discretized to guarantee the maximum deviation.
Points arc_discretize(const Point &p1, const Point &p2, const double radius, const bool ccw, const double deviation);

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

// Returns orientation of a polyline with regard to the center.
// Successive points are expected to take less than a PI angle step.
// Returns Orientation::Unknown if the orientation with regard to the center 
// is not monotonous.
Orientation arc_orientation(
    const Point                 &center,
    const Points::const_iterator begin,
    const Points::const_iterator end);

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

inline bool operator==(const Segment &lhs, const Segment &rhs) {
    return lhs.point == rhs.point && lhs.radius == rhs.radius && lhs.orientation == rhs.orientation;
}

using Segments = std::vector<Segment>;
using Path = Segments;

// Interpolate polyline path with a sequence of linear / circular segments given the interpolation tolerance.
// Only convert to polyline if zero tolerance.
// Convert to polyline and decimate polyline if zero fit_circle_percent_tolerance.
// Path fitting is inspired with the arc fitting algorithm in
//      Arc Welder: Anti-Stutter Library by Brad Hochgesang FormerLurker@pm.me
//      https://github.com/FormerLurker/ArcWelderLib 
Path fit_path(const Points &points, double tolerance, double fit_circle_percent_tolerance);

// Decimate polyline into a smooth path structure using Douglas-Peucker polyline decimation algorithm.
inline Path fit_polyline(const Points &points, double tolerance) { return fit_path(points, tolerance, 0.); }

template<typename FloatType>
inline FloatType segment_length(const Segment &start, const Segment &end)
{
    return end.linear() ?
        (end.point - start.point).cast<FloatType>().norm() :
        arc_length(start.point.cast<FloatType>(), end.point.cast<FloatType>(), FloatType(end.radius));
}

template<typename FloatType>
inline FloatType path_length(const Path &path)
{
    FloatType len = 0;
    for (size_t i = 1; i < path.size(); ++ i)
        len += segment_length<FloatType>(path[i - 1], path[i]);
    return len;
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
        return - int64_t(3) * int64_t(end.radius);
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
    // Projection of the point on the segment.
    Point  point      { 0, 0 };
    // If the point lies on an arc, the arc center is cached here.
    Point  center     { 0, 0 };
    // Square of a distance of the projection.
    double distance2  { std::numeric_limits<double>::max() };

    bool   valid() const { return this->segment_id != std::numeric_limits<size_t>::max(); }
};
// Returns closest segment and a parameter along the closest segment of a path to a point.
PathSegmentProjection point_to_path_projection(const Path &path, const Point &point, double search_radius2 = std::numeric_limits<double>::max());
// Split a path into two paths at a segment point. Snap to an existing point if the projection of "point is closer than min_segment_length.
std::pair<Path, Path> split_at(const Path &path, const PathSegmentProjection &proj, const double min_segment_length);
// Split a path into two paths at a point closest to "point". Snap to an existing point if the projection of "point is closer than min_segment_length.
std::pair<Path, Path> split_at(const Path &path, const Point &point, const double min_segment_length);

} } } // namespace Slic3r::Geometry::ArcWelder

#endif // slic3r_Geometry_ArcWelder_hpp_
