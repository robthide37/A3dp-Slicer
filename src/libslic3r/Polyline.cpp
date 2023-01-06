#include "BoundingBox.hpp"
#include "Polyline.hpp"
#include "Exception.hpp"
#include "ExPolygon.hpp"
#include "ExPolygonCollection.hpp"
#include "Line.hpp"
#include "Polygon.hpp"
#include <iostream>
#include <utility>
#include <algorithm>

namespace Slic3r {

const Point& Polyline::leftmost_point() const
{
    const Point *p = &this->points.front();
    for (Points::const_iterator it = this->points.begin() + 1; it != this->points.end(); ++ it) {
        if (it->x() < p->x()) 
        	p = &(*it);
    }
    return *p;
}

Lines Polyline::lines() const
{
    Lines lines;
    if (this->points.size() >= 2) {
        lines.reserve(this->points.size() - 1);
        for (Points::const_iterator it = this->points.begin(); it != this->points.end()-1; ++it) {
            lines.push_back(Line(*it, *(it + 1)));
        }
    }
    return lines;
}

// removes the given distance from the end of the polyline
void Polyline::clip_end(coordf_t distance)
{
    while (distance > 0) {
        Vec2d  last_point = this->last_point().cast<coordf_t>();
        this->points.pop_back();
        if (this->points.empty())
            break;
        Vec2d  v    = this->last_point().cast<coordf_t>() - last_point;
        coordf_t lsqr = v.squaredNorm();
        if (lsqr > distance * distance) {
            this->points.emplace_back((last_point + v * (distance / sqrt(lsqr))).cast<coord_t>());
            return;
        }
        distance -= sqrt(lsqr);
    }
}

// removes the given distance from the start of the polyline
void Polyline::clip_start(coordf_t distance)
{
    this->reverse();
    this->clip_end(distance);
    if (this->points.size() >= 2)
        this->reverse();
}

void Polyline::extend_end(coordf_t distance)
{
    // relocate last point by extending the last segment by the specified length
    Vec2d v = (this->points.back() - *(this->points.end() - 2)).cast<coordf_t>().normalized();
    this->points.back() += (v * distance).cast<coord_t>();
}

void Polyline::extend_start(coordf_t distance)
{
    // relocate first point by extending the first segment by the specified length
    Vec2d v = (this->points.front() - this->points[1]).cast<coordf_t>().normalized();
    this->points.front() += (v * distance).cast<coord_t>();
}

/* this method returns a collection of points picked on the polygon contour
   so that they are evenly spaced according to the input distance */
Points Polyline::equally_spaced_points(coordf_t distance) const
{
    Points points;
    points.emplace_back(this->first_point());
    double len = 0;
    
    for (Points::const_iterator it = this->points.begin() + 1; it != this->points.end(); ++it) {
        Vec2d  p1 = (it-1)->cast<double>();
        Vec2d  v  = it->cast<double>() - p1;
        coordf_t segment_length = v.norm();
        len += segment_length;
        if (len < distance)
            continue;
        if (len == distance) {
            points.emplace_back(*it);
            len = 0;
            continue;
        }
        coordf_t take = segment_length - (len - distance);  // how much we take of this segment
        points.emplace_back((p1 + v * (take / v.norm())).cast<coord_t>());
        -- it;
        len = - take;
    }
    return points;
}

void Polyline::simplify(coordf_t tolerance)
{
    this->points = MultiPoint::_douglas_peucker(this->points, tolerance);
}

#if 0
// This method simplifies all *lines* contained in the supplied area
template <class T>
void Polyline::simplify_by_visibility(const T &area)
{
    Points &pp = this->points;
    
    size_t s = 0;
    bool did_erase = false;
    for (size_t i = s+2; i < pp.size(); i = s + 2) {
        if (area.contains(Line(pp[s], pp[i]))) {
            pp.erase(pp.begin() + s + 1, pp.begin() + i);
            did_erase = true;
        } else {
            ++s;
        }
    }
    if (did_erase)
        this->simplify_by_visibility(area);
}
template void Polyline::simplify_by_visibility<ExPolygon>(const ExPolygon &area);
template void Polyline::simplify_by_visibility<ExPolygonCollection>(const ExPolygonCollection &area);
#endif

void Polyline::split_at(const Point &point, Polyline* p1, Polyline* p2) const
{
    if (this->size() < 2) {
        *p1 = *this;
        p2->clear();
        return;
    }

    if (this->points.front() == point) {
        *p1 = { point };
        *p2 = *this;
    }

    double min_dist2    = std::numeric_limits<double>::max();
    auto   min_point_it = this->points.cbegin();
    Point  prev         = this->points.front();
    for (auto it = this->points.cbegin() + 1; it != this->points.cend(); ++ it) {
        Point  proj = point.projection_onto(Line(prev, *it));
        double d2   = (proj - point).cast<double>().squaredNorm();
        if (d2 < min_dist2) {
            min_dist2    = d2;
            min_point_it = it;
        }
        prev = *it;
    }

    p1->points.assign(this->points.cbegin(), min_point_it);
    if (p1->points.back() != point)
        p1->points.emplace_back(point);
    
    p2->points = { point };
    if (*min_point_it == point)
        ++ min_point_it;
    p2->points.insert(p2->points.end(), min_point_it, this->points.cend());
}

bool Polyline::is_straight() const
{
    // Check that each segment's direction is equal to the line connecting
    // first point and last point. (Checking each line against the previous
    // one would cause the error to accumulate.)
    double dir = Line(this->first_point(), this->last_point()).direction();
    for (const auto &line: this->lines())
        if (! line.parallel_to(dir))
            return false;
    return true;
}

BoundingBox get_extents(const Polyline &polyline)
{
    return polyline.bounding_box();
}

BoundingBox get_extents(const Polylines &polylines)
{
    BoundingBox bb;
    if (! polylines.empty()) {
        bb = polylines.front().bounding_box();
        for (size_t i = 1; i < polylines.size(); ++ i)
            bb.merge(polylines[i].points);
    }
    return bb;
}

const Point& leftmost_point(const Polylines &polylines)
{
    if (polylines.empty())
        throw Slic3r::InvalidArgument("leftmost_point() called on empty PolylineCollection");
    Polylines::const_iterator it = polylines.begin();
    const Point *p = &it->leftmost_point();
    for (++ it; it != polylines.end(); ++it) {
        const Point *p2 = &it->leftmost_point();
        if (p2->x() < p->x())
            p = p2;
    }
    return *p;
}

bool remove_degenerate(Polylines &polylines)
{
    bool modified = false;
    size_t j = 0;
    for (size_t i = 0; i < polylines.size(); ++ i) {
        if (polylines[i].points.size() >= 2) {
            if (j < i) 
                std::swap(polylines[i].points, polylines[j].points);
            ++ j;
        } else
            modified = true;
    }
    if (j < polylines.size())
        polylines.erase(polylines.begin() + j, polylines.end());
    return modified;
}

std::pair<int, Point> foot_pt(const Points &polyline, const Point &pt)
{
    if (polyline.size() < 2)
        return std::make_pair(-1, Point(0, 0));

    auto  d2_min  = std::numeric_limits<double>::max();
    Point foot_pt_min;
    Point prev = polyline.front();
    auto  it = polyline.begin();
    auto  it_proj = polyline.begin();
    for (++ it; it != polyline.end(); ++ it) {
        Point  foot_pt = pt.projection_onto(Line(prev, *it));
        double d2 = (foot_pt - pt).cast<double>().squaredNorm();
        if (d2 < d2_min) {
            d2_min      = d2;
            foot_pt_min = foot_pt;
            it_proj     = it;
        }
        prev = *it;
    }
    return std::make_pair(int(it_proj - polyline.begin()) - 1, foot_pt_min);
}

ThickLines ThickPolyline::thicklines() const
{
    ThickLines lines;
    if (this->points.size() >= 2) {
        lines.reserve(this->points.size() - 1);
        for (size_t i = 0; i < this->points.size() - 1; ++i) {
            ThickLine line(this->points[i], this->points[i + 1]);
            line.a_width = this->points_width[i];
            line.b_width = this->points_width[i + 1];
            lines.push_back(line);
        }
    }
    return lines;
}

// Removes the given distance from the end of the ThickPolyline
void ThickPolyline::clip_end(coordf_t distance)
{
    assert(this->points_width.size() == this->points.size());
    while (distance > 0) {
        Vec2d    last_point = this->last_point().cast<double>();
        coord_t last_width = this->points_width.back();
        assert(this->points_width.size() == this->points.size());
        this->points.pop_back();
        this->points_width.pop_back();
        if (this->points.empty())
            break;

        assert(this->points_width.size() == this->points.size());
        Vec2d    vec            = this->last_point().cast<double>() - last_point;
        coordf_t width_diff     = this->points_width.back() - last_width;
        coordf_t   vec_length_sqr = vec.squaredNorm();
        if (vec_length_sqr > distance * distance) {
            coordf_t t = (distance / std::sqrt(vec_length_sqr));
            this->points.emplace_back((last_point + vec * t).cast<coord_t>());
            this->points_width.emplace_back(last_width + width_diff * t);
            assert(this->points_width.size() == this->points.size());
            return;
        }

        distance -= std::sqrt(vec_length_sqr);
    }
    assert(this->points_width.size() == this->points.size());
}

Lines3 Polyline3::lines() const
{
    Lines3 lines;
    if (points.size() >= 2)
    {
        lines.reserve(points.size() - 1);
        for (Points3::const_iterator it = points.begin(); it != points.end() - 1; ++it)
        {
            lines.emplace_back(*it, *(it + 1));
        }
    }
    return lines;
}

void concatThickPolylines(ThickPolylines& pp) {
    bool changes = true;
    while (changes){
        changes = false;
        //concat polyline if only 2 polyline at a point
        for (size_t i = 0; i < pp.size(); ++i) {
            ThickPolyline *polyline = &pp[i];
            if (polyline->first_point().coincides_with_epsilon(polyline->last_point())) {
                polyline->endpoints.first = false;
                polyline->endpoints.second = false;
                continue;
            }

            size_t id_candidate_first_point = -1;
            size_t id_candidate_last_point = -1;
            size_t nbCandidate_first_point = 0;
            size_t nbCandidate_last_point = 0;
            // find another polyline starting here
            for (size_t j = 0; j < pp.size(); ++j) {
                if (j == i) continue;
                ThickPolyline *other = &pp[j];
                if (polyline->last_point().coincides_with_epsilon(other->last_point())) {
                    id_candidate_last_point = j;
                    nbCandidate_last_point++;
                }
                if (polyline->last_point().coincides_with_epsilon(other->first_point())) {
                    id_candidate_last_point = j;
                    nbCandidate_last_point++;
                }
                if (polyline->first_point().coincides_with_epsilon(other->last_point())) {
                    id_candidate_first_point = j;
                    nbCandidate_first_point++;
                }
                if (polyline->first_point().coincides_with_epsilon(other->first_point())) {
                    id_candidate_first_point = j;
                    nbCandidate_first_point++;
                }
            }
            if (id_candidate_last_point == id_candidate_first_point && nbCandidate_first_point == 1 && nbCandidate_last_point == 1) {
                if (polyline->first_point().coincides_with_epsilon(pp[id_candidate_first_point].first_point())) pp[id_candidate_first_point].reverse();
                // it's a trap! it's a  loop!
                polyline->points.insert(polyline->points.end(), pp[id_candidate_first_point].points.begin() + 1, pp[id_candidate_first_point].points.end());
                polyline->points_width.insert(polyline->points_width.end(), pp[id_candidate_first_point].points_width.begin() + 1, pp[id_candidate_first_point].points_width.end());
                pp.erase(pp.begin() + id_candidate_first_point);
                changes = true;
                polyline->endpoints.first = false;
                polyline->endpoints.second = false;
            } else {

                if (nbCandidate_first_point == 1) {
                    if (polyline->first_point().coincides_with_epsilon(pp[id_candidate_first_point].first_point())) pp[id_candidate_first_point].reverse();
                    //concat at front
                    polyline->points_width[0] = std::max(polyline->points_width.front(), pp[id_candidate_first_point].points_width.back());
                    polyline->points.insert(polyline->points.begin(), pp[id_candidate_first_point].points.begin(), pp[id_candidate_first_point].points.end() - 1);
                    polyline->points_width.insert(polyline->points_width.begin(), pp[id_candidate_first_point].points_width.begin(), pp[id_candidate_first_point].points_width.end() - 1);
                    polyline->endpoints.first = pp[id_candidate_first_point].endpoints.first;
                    pp.erase(pp.begin() + id_candidate_first_point);
                    changes = true;
                    if (id_candidate_first_point < i) {
                        i--;
                        polyline = &pp[i];
                    }
                    if (id_candidate_last_point > id_candidate_first_point) {
                        id_candidate_last_point--;
                    }
                } else if (nbCandidate_first_point == 0) {
                    //update endpoint
                    polyline->endpoints.first = true;
                }
                if (nbCandidate_last_point == 1) {
                    if (polyline->last_point().coincides_with_epsilon(pp[id_candidate_last_point].last_point())) pp[id_candidate_last_point].reverse();
                    //concat at back
                    polyline->points_width[polyline->points_width.size() - 1] = std::max(polyline->points_width.back(), pp[id_candidate_last_point].points_width.front());
                    polyline->points.insert(polyline->points.end(), pp[id_candidate_last_point].points.begin() + 1, pp[id_candidate_last_point].points.end());
                    polyline->points_width.insert(polyline->points_width.end(), pp[id_candidate_last_point].points_width.begin() + 1, pp[id_candidate_last_point].points_width.end());
                    polyline->endpoints.second = pp[id_candidate_last_point].endpoints.second;
                    pp.erase(pp.begin() + id_candidate_last_point);
                    changes = true;
                    if (id_candidate_last_point < i) {
                        i--;
                        polyline = &pp[i];
                    }
                } else if (nbCandidate_last_point == 0) {
                    //update endpoint
                    polyline->endpoints.second = true;
                }

                if (polyline->last_point().coincides_with_epsilon(polyline->first_point())) {
                    //the concat has created a loop : update endpoints
                    polyline->endpoints.first = false;
                    polyline->endpoints.second = false;
                }
            }
        }
    }
}

}
