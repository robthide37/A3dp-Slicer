///|/ Copyright (c) Prusa Research 2016 - 2023 Tomáš Mészáros @tamasmeszaros, Vojtěch Bubník @bubnikv, Lukáš Hejl @hejllukas, Enrico Turri @enricoturri1966
///|/ Copyright (c) Slic3r 2013 - 2016 Alessandro Ranellucci @alranel
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_MultiPoint_hpp_
#define slic3r_MultiPoint_hpp_

#include "libslic3r.h"
#include <algorithm>
#include <vector>
#include "Line.hpp"
#include "Point.hpp"

namespace Slic3r {

class BoundingBox;
class BoundingBox3;
#ifdef _DEBUG
// Reduces polyline in the <begin, end) range, outputs into the output iterator.
// Output iterator may be equal to input iterator as long as the iterator value type move operator supports move at the same input / output address.
template<typename SquareLengthType, typename InputIterator, typename OutputIterator, typename PointGetter>
inline OutputIterator douglas_peucker_old(InputIterator begin, InputIterator end, OutputIterator out, const coord_t tolerance, PointGetter point_getter)
{
    using InputIteratorCategory = typename std::iterator_traits<InputIterator>::iterator_category;
    static_assert(std::is_base_of_v<std::input_iterator_tag, InputIteratorCategory>);
    using Vector = Eigen::Matrix<SquareLengthType, 2, 1, Eigen::DontAlign>;
    if (begin != end) {
        // Supporting in-place reduction and the data type may be generic, thus we are always making a copy of the point value before there is a chance
        // to override input by moving the data to the output.
        auto a = point_getter(*begin);
        *out ++ = std::move(*begin);
        if (auto next = std::next(begin); next == end) {
            // Single point input only.
        } else if (std::next(next) == end) {
            // Two points input.
            *out ++ = std::move(*next);
        } else {
            const SquareLengthType tolerance_sq = (SquareLengthType(tolerance*tolerance));
            InputIterator anchor  = begin;
            InputIterator floater = std::prev(end);
            std::vector<InputIterator> dpStack;
            if constexpr (std::is_base_of_v<std::random_access_iterator_tag, InputIteratorCategory>)
                dpStack.reserve(end - begin);
            dpStack.emplace_back(floater);
            auto f = point_getter(*floater);
            for (;;) {
                assert(anchor != floater);
                bool            take_floater = false;
                InputIterator   furthest     = anchor;
                if (std::next(anchor) == floater) {
                    // Two point segment. Accept the floater.
                    take_floater = true;
                } else {
                    SquareLengthType max_dist_sq = 0;
                    // Find point furthest from line seg created by (anchor, floater) and note it.
                    const Vector v = (f - a).template cast<SquareLengthType>();
                    if (const SquareLengthType l2 = v.squaredNorm(); l2 == 0) {
                        // Zero length segment, find the furthest point between anchor and floater.
                        for (auto it = std::next(anchor); it != floater; ++ it)
                            if (SquareLengthType dist_sq = (point_getter(*it) - a).template cast<SquareLengthType>().squaredNorm(); 
                                dist_sq > max_dist_sq) {
                                max_dist_sq  = dist_sq;
                                furthest = it;
                            }
                    } else {
                        // Find Find the furthest point from the line <anchor, floater>.
                        const double dl2 = double(l2);
                        const Vec2d  dv  = v.template cast<double>();
                        for (auto it = std::next(anchor); it != floater; ++ it) {
                            const auto   p  = point_getter(*it);
                            const Vector va = (p - a).template cast<SquareLengthType>();
                            const SquareLengthType t = va.dot(v);
                            SquareLengthType dist_sq;
                            if (t <= 0) {
                                dist_sq = va.squaredNorm();
                            } else if (t >= l2) {
                                dist_sq = (p - f).template cast<SquareLengthType>().squaredNorm();
                            } else if (double dt = double(t) / dl2; dt <= 0) {
                                // is this case useful? seems not
                                assert(false);
                                dist_sq = va.squaredNorm();
                            } else if (dt >= 1.) {
                                // is this case useful? seems not
                                assert(false);
                                dist_sq = (p - f).template cast<SquareLengthType>().squaredNorm();
                            } else {
                                const Vector w = (dt * dv).cast<SquareLengthType>();
                                const Vec2crd vtemp = (w - va);
                                const double normal_dist = vtemp.cast<double>().norm();
                                dist_sq = (w - va).squaredNorm();
                            }
                            if (dist_sq > max_dist_sq) {
                                max_dist_sq  = dist_sq;
                                furthest     = it;
                            }
                        }                        
                    }
                    // remove point if less than tolerance
                    take_floater = max_dist_sq <= tolerance_sq;
                }
                if (take_floater) {
                    // The points between anchor and floater are close to the <anchor, floater> line.
                    // Drop the points between them.
                    a = f;
                    *out ++ = std::move(*floater);
                    anchor = floater;
                    assert(dpStack.back() == floater);
                    dpStack.pop_back();
                    if (dpStack.empty())
                        break;
                    floater = dpStack.back();
                    f = point_getter(*floater);
                } else {
                    // The furthest point is too far from the segment <anchor, floater>. 
                    // Divide recursively.
                    floater = furthest;
                    f = point_getter(*floater);
                    dpStack.emplace_back(floater);
                }
            }
        }
    }
    return out;
}
#endif


// Reduces polyline in the <begin, end) range, outputs into the output iterator.
// Output iterator may be equal to input iterator as long as the iterator value type move operator supports move at the same input / output address.
// note: SquareLengthType is int64 becasue it's not like it won't be enough, we're looking at deviation, not path length.
template<typename InputIterator, typename OutputIterator, typename PointGetter>
inline OutputIterator douglas_peucker_double(InputIterator begin, InputIterator end, OutputIterator out, const coord_t tolerance, PointGetter point_getter)
{
    using InputIteratorCategory = typename std::iterator_traits<InputIterator>::iterator_category;
    static_assert(std::is_base_of_v<std::input_iterator_tag, InputIteratorCategory>);
    if (begin != end) {
        // Supporting in-place reduction and the data type may be generic, thus we are always making a copy of the point value before there is a chance
        // to override input by moving the data to the output.
        Point pt_start = point_getter(*begin);
        *out ++ = std::move(*begin);
        if (InputIterator next = std::next(begin); next == end) {
            // Single point input only.
        } else if (std::next(next) == end) {
            // Two points input.
            *out ++ = std::move(*next);
        } else {
            const double tolerance_sq_d = std::max(100.*100., sqr(double(tolerance)));
            InputIterator anchor  = begin;
            InputIterator floater = std::prev(end);
            std::vector<InputIterator> dpStack;
            if constexpr (std::is_base_of_v<std::random_access_iterator_tag, InputIteratorCategory>)
                dpStack.reserve(end - begin);
            dpStack.emplace_back(floater);
            Point pt_floater = point_getter(*floater);
            for (;;) {
                assert(anchor != floater);
                bool            take_floater = false;
                InputIterator   furthest     = anchor;
                if (std::next(anchor) == floater) {
                    // Two point segment. Accept the floater.
                    take_floater = true;
                } else {
                    double max_dist_sq_d = 0;
                    // Find point furthest from line seg created by (anchor, floater) and note it.
                    const Vec2d vec_af_d = (pt_floater - pt_start).cast<double>();
                    const double length_squared_d = vec_af_d.squaredNorm();
                    if ( length_squared_d == 0) {
                        // Zero length segment, find the furthest point between anchor and floater.
                        for (InputIterator it = std::next(anchor); it != floater; ++it) {
                            const Point ptit = point_getter(*it);
                            const Vec2d pttemp = (ptit - pt_start).cast<double>();
                            double dist_sq_d = pttemp.squaredNorm();
                            if (dist_sq_d > max_dist_sq_d) {
                                max_dist_sq_d = dist_sq_d;
                                furthest = it;
                            }
                        }
                    } else {
                        // Find Find the furthest point from the line <anchor, floater>.
                        for (InputIterator it = std::next(anchor); it != floater; ++ it) {
                            const Point &pt_check  = point_getter(*it);
                            const Vec2d vec_ap_d = (pt_check - pt_start).cast<double>();
                            const double dot_d = vec_ap_d.dot(vec_af_d);
                            double dist_sq_d;
                            if (dot_d <= 0) {
                                // point below pt_start, or at 90°
                                dist_sq_d = vec_ap_d.squaredNorm();
                            } else if (dot_d >= length_squared_d) {
                                // point after pt_floater, or at 90°
                                dist_sq_d = (pt_check - pt_floater).cast<double>().squaredNorm();
                            } else {
                                // point between pt_start and pt_floater
                                dist_sq_d = (((dot_d / length_squared_d) * vec_af_d) - vec_ap_d).squaredNorm();
                            }
                            if (dist_sq_d > max_dist_sq_d) {
                                max_dist_sq_d  = dist_sq_d;
                                furthest     = it;
                            }
                        }
                    }
                    // remove point if less than tolerance
                    take_floater = max_dist_sq_d <= tolerance_sq_d;
                }
                if (take_floater) {
                    // The points between anchor and floater are close to the <anchor, floater> line.
                    // Drop the points between them.
                    pt_start = pt_floater;
                    *out ++ = std::move(*floater);
                    anchor = floater;
                    assert(dpStack.back() == floater);
                    dpStack.pop_back();
                    if (dpStack.empty()) {
                        break;
                    }
                    floater = dpStack.back();
                    pt_floater = point_getter(*floater);
                } else {
                    // The furthest point is too far from the segment <anchor, floater>. 
                    // Divide recursively.
                    floater = furthest;
                    pt_floater = point_getter(*floater);
                    dpStack.emplace_back(floater);
                }
            }
        }
    }
    return out;
}

// Reduces polyline in the <begin, end) range, outputs into the output iterator.
// Output iterator may be equal to input iterator as long as the iterator value type move operator supports move at the same input / output address.
// note: SquareLengthType is int64 becasue it's not like it won't be enough, we're looking at deviation, not path length.
// note: loss of precision in dot() is too much and the result is too different.
template<typename InputIterator, typename OutputIterator, typename PointGetter>
inline OutputIterator douglas_peucker_int(InputIterator begin, InputIterator end, OutputIterator out, const coord_t tolerance, PointGetter point_getter)
{
    using InputIteratorCategory = typename std::iterator_traits<InputIterator>::iterator_category;
    static_assert(std::is_base_of_v<std::input_iterator_tag, InputIteratorCategory>);
    if (begin != end) {
        // Supporting in-place reduction and the data type may be generic, thus we are always making a copy of the point value before there is a chance
        // to override input by moving the data to the output.
        Point pt_start = point_getter(*begin);
        *out ++ = std::move(*begin);
        if (InputIterator next = std::next(begin); next == end) {
            // Single point input only.
        } else if (std::next(next) == end) {
            // Two points input.
            *out ++ = std::move(*next);
        } else {
            const lengthsqr_t tolerance_sq = std::max(lengthsqr_t(1), Slic3r::coord_int_sqr(tolerance));
            const double tolerance_sq_d = std::max(128*128., sqr(double(tolerance)));
            InputIterator anchor  = begin;
            InputIterator floater = std::prev(end);
            std::vector<InputIterator> dpStack;
            if constexpr (std::is_base_of_v<std::random_access_iterator_tag, InputIteratorCategory>)
                dpStack.reserve(end - begin);
            dpStack.emplace_back(floater);
            Point pt_floater = point_getter(*floater);
            for (;;) {
                assert(anchor != floater);
                bool            take_floater = false;
                InputIterator   furthest     = anchor;
                if (std::next(anchor) == floater) {
                    // Two point segment. Accept the floater.
                    take_floater = true;
                } else {
                    lengthsqr_t max_dist_sq = 0;
                    double max_dist_sq_d = 0;
                    // Find point furthest from line seg created by (anchor, floater) and note it.
                    const Vec2crd vec_af = (pt_floater - pt_start);
                    const Vec2d vec_af_d = vec_af.cast<double>();
                    const lengthsqr_t length_squared = squared_int_norm(vec_af);
                    const double length_squared_d = vec_af_d.squaredNorm();
                    assert(length_squared_d !=0 || length_squared == 0);
                    if ( length_squared == 0) {
                        // Zero length segment, find the furthest point between anchor and floater.
                        for (InputIterator it = std::next(anchor); it != floater; ++it) {
                            lengthsqr_t dist_sq = squared_int_norm(point_getter(*it) - pt_start);
                            const Vec2crd tempcrd = (point_getter(*it) - pt_start);
                            const Vec2d temppt = tempcrd.cast<double>();
                            const double dist_sq_d = temppt.squaredNorm();
                            if ( dist_sq > max_dist_sq) {
                                assert(dist_sq_d > max_dist_sq_d);
                                max_dist_sq_d = dist_sq_d;
                                max_dist_sq = dist_sq;
                                furthest = it;
                            } else {
                                assert(dist_sq_d <= max_dist_sq_d);
                            }
                        }
                    } else {
                        // Find Find the furthest point from the line <anchor, floater>.
                        const double d_length_squared = double(length_squared);
                        const Vec2d  d_vec_af  = vec_af.cast<double>();
                        for (InputIterator it = std::next(anchor); it != floater; ++ it) {
                            const Point  &pt_check  = point_getter(*it);
                            const Vec2crd vec_ap = (pt_check - pt_start);
                            const Vec2d vec_ap_d = (pt_check - pt_start).cast<double>();
                            const int64_t dot = dot_int(vec_ap, vec_af);
                            const double dot_d = vec_ap_d.dot(vec_af_d);
                            lengthsqr_t dist_sq;
                            double dist_sq_d;
                            if (dot <= 0) {
                                assert(dot_d <= 0);
                                // vec_ap and vec_af are in opposite direction
                                dist_sq = squared_int_norm(vec_ap);
                                dist_sq_d = vec_ap_d.squaredNorm();
                                if (dist_sq > max_dist_sq) {
                                    assert(dist_sq <= tolerance_sq == (dist_sq_d <= tolerance_sq_d));
                                }
                            } else
                                // dot is >0, so can be cast to uint64_t if needed
                                if (lengthsqr_t(dot) >= length_squared) {
                                assert(dot_d >= length_squared_d);
                                // vec_ap and vec_af are in same direction, and the angle is smaller than almost 90°
                                dist_sq = squared_int_norm(pt_check - pt_floater);
                                dist_sq_d = (pt_check - pt_floater).cast<double>().squaredNorm();
                                if (dist_sq > max_dist_sq) {
                                    assert(dist_sq <= tolerance_sq == (dist_sq_d <= tolerance_sq_d));
                                }
                            } else {
                                assert(dot_d < length_squared_d && dot > 0);
                                // vec_ap and vec_af are in same direction, angle is (or almost) 90°
                                const Vec2d w_d = ((dot_d / length_squared_d) * d_vec_af);
                                const Vec2d w = ((double(dot)/ length_squared) * d_vec_af);
                                const Vec2crd vtemp = w.cast<coord_t>() - vec_ap;
                                const Vec2d vtemp_d = w_d - vec_ap_d;
                                const double normal_dist = vtemp.cast<double>().norm();
                                const double normal_dist_d = vtemp_d.norm();
                                dist_sq_d = vtemp_d.squaredNorm();
                                if (normal_dist > std::numeric_limits<int32_t>::max()) {
                                    dist_sq = std::numeric_limits<int64_t>::max();
                                    assert(false);
                                } else {
                                    dist_sq = 1 + squared_int_norm(vtemp);
                                    if (dist_sq > max_dist_sq) {
                                        assert(dist_sq <= tolerance_sq == (dist_sq_d <= tolerance_sq_d));
                                    }
                                }
                            }
                            if (dist_sq > max_dist_sq) {
                                assert(dist_sq_d >= max_dist_sq_d);
                                max_dist_sq_d  = dist_sq_d;
                                max_dist_sq  = dist_sq;
                                furthest     = it;
                                assert(max_dist_sq <= tolerance_sq == (max_dist_sq_d <= tolerance_sq_d));
                            } else {
                                assert(dist_sq_d <= max_dist_sq_d);
                            }
                        }
                    }
                    // remove point if less than tolerance
                    take_floater = max_dist_sq <= tolerance_sq;
                    assert(take_floater == (max_dist_sq_d <= tolerance_sq_d));
                }
                if (take_floater) {
                    // The points between anchor and floater are close to the <anchor, floater> line.
                    // Drop the points between them.
                    pt_start = pt_floater;
                    *out ++ = std::move(*floater);
                    anchor = floater;
                    assert(dpStack.back() == floater);
                    dpStack.pop_back();
                    if (dpStack.empty()) {
                        break;
                    }
                    floater = dpStack.back();
                    pt_floater = point_getter(*floater);
                } else {
                    // The furthest point is too far from the segment <anchor, floater>. 
                    // Divide recursively.
                    floater = furthest;
                    pt_floater = point_getter(*floater);
                    dpStack.emplace_back(floater);
                }
            }
        }
    }
    return out;
}
template<typename InputIterator, typename OutputIterator, typename PointGetter>
inline OutputIterator douglas_peucker_impl(
    InputIterator begin, InputIterator end, OutputIterator out, const coord_t tolerance, PointGetter point_getter) {
    return douglas_peucker_double(begin, end, out, tolerance, point_getter);
}

// Reduces polyline in the <begin, end) range, outputs into the output iterator.
// Output iterator may be equal to input iterator as long as the iterator value type move operator supports move at the same input / output address.
template<typename OutputIterator>
inline OutputIterator douglas_peucker(Points::const_iterator begin, Points::const_iterator end, OutputIterator out, const coord_t tolerance)
{
    return douglas_peucker_impl(begin, end, out, tolerance, [](const Point &p) { return p; });
}

inline Points douglas_peucker(const Points &src, const coord_t tolerance) 
{
    Points out;
    out.reserve(src.size());
    douglas_peucker(src.begin(), src.end(), std::back_inserter(out), tolerance);
    return out;
}

class MultiPoint
{
public:
    // TODO: makes that private?
    Points points;

    MultiPoint() = default;
    MultiPoint(const MultiPoint &other) : points(other.points) {}
    MultiPoint(MultiPoint &&other) : points(std::move(other.points)) {}
    MultiPoint(std::initializer_list<Point> list) : points(list) {}
    explicit MultiPoint(const Points &_points) : points(_points) {}
    explicit MultiPoint(Points &&_points) : points(std::move(_points)) {}
    MultiPoint &operator=(const MultiPoint &other) {
        points = other.points;
        return *this;
    }
    MultiPoint &operator=(MultiPoint &&other) {
        points = std::move(other.points);
        return *this;
    }
    void scale(double factor);
    void scale(double factor_x, double factor_y);
    void translate(double x, double y) { this->translate(Point(coord_t(x), coord_t(y))); }
    void translate(const Vector &vector);
    void rotate(double angle) { this->rotate(cos(angle), sin(angle)); }
    void rotate(double cos_angle, double sin_angle);
    void rotate(double angle, const Point &center);
    virtual void reverse() { std::reverse(this->points.begin(), this->points.end()); }

    const Point &front() const { return this->points.front(); }
    const Point &back() const { return this->points.back(); }
    const Point &first_point() const { return this->front(); }
    virtual bool is_loop() const { return size() <= 1 || front() == back(); }
    size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }
    bool is_valid() const { return this->points.size() >= 2; }

    // Return index of a polygon point exactly equal to point.
    // Return -1 if no such point exists.
    int find_point(const Point &point) const;
    // Return index of the closest point to point closer than scaled_epsilon.
    // Return -1 if no such point exists.
    int find_point(const Point &point, const coordf_t scaled_epsilon) const;
    int closest_point_index(const Point &point) const {
        int idx = -1;
        if (!this->points.empty()) {
            idx = 0;
            double dist_min = (point - this->points.front()).cast<double>().norm();
            for (int i = 1; i < int(this->points.size()); ++i) {
                double d = (this->points[i] - point).cast<double>().norm();
                if (d < dist_min) {
                    dist_min = d;
                    idx = i;
                }
            }
        }
        return idx;
    }
    const Point *closest_point(const Point &point) const {
        return this->points.empty() ? nullptr : &this->points[this->closest_point_index(point)];
    }
    BoundingBox bounding_box() const;
    // Return true if there are exact duplicates.
    bool has_duplicate_points() const;
    // Remove exact duplicates, return true if any duplicate has been removed.
    bool remove_duplicate_points();
    virtual void douglas_peucker(coord_t tolerance = SCALED_EPSILON) {
        auto it_end = Slic3r::douglas_peucker(points.begin(), this->points.end(), this->points.begin(), double(tolerance));
        assert(it_end <= points.end());
        points.resize(std::distance(points.begin(), it_end));
    }
    virtual void clear() { this->points.clear(); }
    void append(const Point &point) { this->points.push_back(point); }
    void append(const Points &src) { this->append(src.begin(), src.end()); }
    void append(const Points::const_iterator &begin, const Points::const_iterator &end) {
        this->points.insert(this->points.end(), begin, end);
    }
    void append(Points &&src) {
        if (this->points.empty()) {
            this->points = std::move(src);
        } else {
            this->points.insert(this->points.end(), src.begin(), src.end());
            src.clear();
        }
    }

    static Points douglas_peucker(const Points &src, const coord_t tolerance) {
        return Slic3r::douglas_peucker(src, tolerance);
    }
    static Points _douglas_peucker_plus(const Points& pts, const double tolerance, const double min_length);
    static Points visivalingam(const Points &src, const double tolerance);

    // Projection of a point onto the lines defined by the points.
    virtual std::pair<Point, size_t> point_projection(const Point &point) const;

    inline auto begin() { return points.begin(); }
    inline auto begin() const { return points.begin(); }
    inline auto end() { return points.end(); }
    inline auto end() const { return points.end(); }
    inline auto cbegin() const { return points.begin(); }
    inline auto cend() const { return points.end(); }
    inline auto rbegin() { return points.rbegin(); }
    inline auto rbegin() const { return points.rbegin(); }
    inline auto rend() { return points.rend(); }
    inline auto rend() const { return points.rend(); }
    inline auto crbegin() const { return points.crbegin(); }
    inline auto crend() const { return points.crend(); }

#ifdef _DEBUGINFO
    virtual void assert_valid() const;
    // to create a cpp multipoint to create test units.
    std::string to_debug_string();
#else
    void assert_valid() const;
#endif
};

class MultiPoint3
{
public:
    Points3 points;

    void append(const Vec3crd& point) { this->points.push_back(point); }

    void translate(double x, double y);
    void translate(const Point& vector);
    bool is_valid() const { return this->points.size() >= 2; }

    BoundingBox3 bounding_box() const;

    // Remove exact duplicates, return true if any duplicate has been removed.
    bool remove_duplicate_points();
};

extern BoundingBox get_extents(const MultiPoint &mp);
extern BoundingBox get_extents_rotated(const Points &points, double angle);
extern BoundingBox get_extents_rotated(const MultiPoint &mp, double angle);

inline double length(const Points::const_iterator begin, const Points::const_iterator end) {
    double total = 0;
    if (begin != end) {
        auto it = begin;
        for (auto it_prev = it ++; it != end; ++ it, ++ it_prev)
            total += (*it - *it_prev).cast<double>().norm();
    }
    return total;
}

inline double length(const Points &pts) {
    return length(pts.begin(), pts.end());
}

inline double area(const Points &polygon) {
    double area = 0.;
    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i ++)
		area += double(polygon[i](0) + polygon[j](0)) * double(polygon[i](1) - polygon[j](1));
    return area;
}

} // namespace Slic3r

#endif
