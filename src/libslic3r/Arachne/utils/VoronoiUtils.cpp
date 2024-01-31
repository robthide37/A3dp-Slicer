//Copyright (c) 2021 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include <stack>
#include <optional>
#include <boost/log/trivial.hpp>

#include "linearAlg2D.hpp"
#include "VoronoiUtils.hpp"
#include "libslic3r/Geometry/VoronoiUtils.hpp"

namespace Slic3r::Arachne
{

Vec2i64 VoronoiUtils::p(const vd_t::vertex_type *node)
{
    const double x = node->x();
    const double y = node->y();
    assert(std::isfinite(x) && std::isfinite(y));
    assert(x <= double(std::numeric_limits<int64_t>::max()) && x >= std::numeric_limits<int64_t>::lowest());
    assert(y <= double(std::numeric_limits<int64_t>::max()) && y >= std::numeric_limits<int64_t>::lowest());
    return {int64_t(x + 0.5 - (x < 0)), int64_t(y + 0.5 - (y < 0))}; // Round to the nearest integer coordinates.
}

Point VoronoiUtils::getSourcePoint(const vd_t::cell_type& cell, const std::vector<Segment>& segments)
{
    assert(cell.contains_point());
    if(!cell.contains_point())
        BOOST_LOG_TRIVIAL(debug) << "Voronoi cell doesn't contain a source point!";

    switch (cell.source_category()) {
        case boost::polygon::SOURCE_CATEGORY_SINGLE_POINT:
            assert(false && "Voronoi diagram is always constructed using segments, so cell.source_category() shouldn't be SOURCE_CATEGORY_SINGLE_POINT!\n");
            BOOST_LOG_TRIVIAL(error) << "Voronoi diagram is always constructed using segments, so cell.source_category() shouldn't be SOURCE_CATEGORY_SINGLE_POINT!";
            break;
        case boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT:
            assert(cell.source_index() < segments.size());
            return boost::polygon::segment_traits<Segment>::get(segments[cell.source_index()], boost::polygon::LOW);
            break;
        case boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT:
            assert(cell.source_index() < segments.size());
            return boost::polygon::segment_traits<Segment>::get(segments[cell.source_index()], boost::polygon::HIGH);
            break;
        default:
            assert(false && "getSourcePoint should only be called on point cells!\n");
            break;
    }

    assert(false && "cell.source_category() is equal to an invalid value!\n");
    BOOST_LOG_TRIVIAL(error) << "cell.source_category() is equal to an invalid value!";
    return {};
}

PolygonsPointIndex VoronoiUtils::getSourcePointIndex(const vd_t::cell_type& cell, const std::vector<Segment>& segments)
{
    assert(cell.contains_point());
    if(!cell.contains_point())
        BOOST_LOG_TRIVIAL(debug) << "Voronoi cell doesn't contain a source point!";

    assert(cell.source_category() != boost::polygon::SOURCE_CATEGORY_SINGLE_POINT);
    switch (cell.source_category()) {
        case boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT: {
            assert(cell.source_index() < segments.size());
            return segments[cell.source_index()];
            break;
        }
        case boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT: {
            assert(cell.source_index() < segments.size());
            return segments[cell.source_index()].next();
            break;
        }
        default:
            assert(false && "getSourcePoint should only be called on point cells!\n");
            break;
    }
    PolygonsPointIndex ret = segments[cell.source_index()];
    return ++ret;
}

const VoronoiUtils::Segment &VoronoiUtils::getSourceSegment(const vd_t::cell_type &cell, const std::vector<Segment> &segments)
{
    assert(cell.contains_segment());
    if (!cell.contains_segment())
        BOOST_LOG_TRIVIAL(debug) << "Voronoi cell doesn't contain a source segment!";

    return segments[cell.source_index()];
}

Points VoronoiUtils::discretizeParabola(const Point &source_point, const Segment &source_segment, const Point &start, const Point &end, const coord_t approximate_step_size, float transitioning_angle)
{
    Points discretized;
    // x is distance of point projected on the segment ab
    // xx is point projected on the segment ab
    const Point   a       = source_segment.from();
    const Point   b       = source_segment.to();
    const Point   ab      = b - a;
    const Point   as      = start - a;
    const Point   ae      = end - a;
    const coord_t ab_size = ab.cast<int64_t>().norm();
    const coord_t sx      = as.cast<int64_t>().dot(ab.cast<int64_t>()) / ab_size;
    const coord_t ex      = ae.cast<int64_t>().dot(ab.cast<int64_t>()) / ab_size;
    const coord_t sxex    = ex - sx;

    const Point   ap = source_point - a;
    const coord_t px = ap.cast<int64_t>().dot(ab.cast<int64_t>()) / ab_size;

    Point pxx;
    Line(a, b).distance_to_infinite_squared(source_point, &pxx);
    const Point   ppxx = pxx - source_point;
    const coord_t d    = ppxx.cast<int64_t>().norm();

    const Vec2d  rot           = perp(ppxx).cast<double>().normalized();
    const double rot_cos_theta = rot.x();
    const double rot_sin_theta = rot.y();

    if (d == 0) {
        discretized.emplace_back(start);
        discretized.emplace_back(end);
        return discretized;
    }

    const double marking_bound = atan(transitioning_angle * 0.5);
    int64_t      msx           = -marking_bound * int64_t(d); // projected marking_start
    int64_t      mex           = marking_bound * int64_t(d);  // projected marking_end

    const coord_t marking_start_end_h = msx * msx / (2 * d) + d / 2;
    Point         marking_start       = Point(coord_t(msx), marking_start_end_h).rotated(rot_cos_theta, rot_sin_theta) + pxx;
    Point         marking_end         = Point(coord_t(mex), marking_start_end_h).rotated(rot_cos_theta, rot_sin_theta) + pxx;
    const int     dir                 = (sx > ex) ? -1 : 1;
    if (dir < 0) {
        std::swap(marking_start, marking_end);
        std::swap(msx, mex);
    }

    bool add_marking_start = msx * int64_t(dir) > int64_t(sx - px) * int64_t(dir) && msx * int64_t(dir) < int64_t(ex - px) * int64_t(dir);
    bool add_marking_end   = mex * int64_t(dir) > int64_t(sx - px) * int64_t(dir) && mex * int64_t(dir) < int64_t(ex - px) * int64_t(dir);

    const Point apex     = Point(0, d / 2).rotated(rot_cos_theta, rot_sin_theta) + pxx;
    bool        add_apex = int64_t(sx - px) * int64_t(dir) < 0 && int64_t(ex - px) * int64_t(dir) > 0;

    assert(!add_marking_start || !add_marking_end || add_apex);
    if (add_marking_start && add_marking_end && !add_apex)
        BOOST_LOG_TRIVIAL(warning) << "Failing to discretize parabola! Must add an apex or one of the endpoints.";

    const coord_t step_count = lround(static_cast<double>(std::abs(ex - sx)) / approximate_step_size);
    discretized.emplace_back(start);
    for (coord_t step = 1; step < step_count; ++step) {
        const int64_t x = int64_t(sx) + int64_t(sxex) * int64_t(step) / int64_t(step_count) - int64_t(px);
        const int64_t y = int64_t(x) * int64_t(x) / int64_t(2 * d) + int64_t(d / 2);

        if (add_marking_start && msx * int64_t(dir) < int64_t(x) * int64_t(dir)) {
            discretized.emplace_back(marking_start);
            add_marking_start = false;
        }

        if (add_apex && int64_t(x) * int64_t(dir) > 0) {
            discretized.emplace_back(apex);
            add_apex = false; // only add the apex just before the
        }

        if (add_marking_end && mex * int64_t(dir) < int64_t(x) * int64_t(dir)) {
            discretized.emplace_back(marking_end);
            add_marking_end = false;
        }

        assert(Geometry::VoronoiUtils::is_in_range<coord_t>(x) && Geometry::VoronoiUtils::is_in_range<coord_t>(y));
        const Point result = Point(x, y).rotated(rot_cos_theta, rot_sin_theta) + pxx;
        discretized.emplace_back(result);
    }

    if (add_apex)
        discretized.emplace_back(apex);

    if (add_marking_end)
        discretized.emplace_back(marking_end);

    discretized.emplace_back(end);
    return discretized;
}

}//namespace Slic3r::Arachne
