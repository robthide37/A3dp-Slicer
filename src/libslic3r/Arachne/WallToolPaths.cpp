// Copyright (c) 2020 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include <algorithm> //For std::partition_copy and std::min_element.
#include <unordered_set>

#include "WallToolPaths.hpp"

#include "SkeletalTrapezoidation.hpp"
#include "../ClipperUtils.hpp"
#include "Arachne/utils/linearAlg2D.hpp"
#include "EdgeGrid.hpp"
#include "utils/SparseLineGrid.hpp"
#include "Geometry.hpp"

namespace Slic3r::Arachne
{

WallToolPaths::WallToolPaths(const Polygons& outline, const coord_t nominal_bead_width, const size_t inset_count, const coord_t wall_0_inset,
                             const PrintObjectConfig &print_object_config)
    : outline(outline)
    , bead_width_0(nominal_bead_width)
    , bead_width_x(nominal_bead_width)
    , inset_count(inset_count)
    , wall_0_inset(wall_0_inset)
    , strategy_type(print_object_config.beading_strategy_type.value)
    , print_thin_walls(Slic3r::Arachne::fill_outline_gaps)
    , min_feature_size(scaled<coord_t>(print_object_config.min_feature_size.value))
    , min_bead_width(scaled<coord_t>(print_object_config.min_bead_width.value))
    , small_area_length(static_cast<double>(nominal_bead_width) / 2.)
    , toolpaths_generated(false)
    , print_object_config(print_object_config)
{
}

WallToolPaths::WallToolPaths(const Polygons& outline, const coord_t bead_width_0, const coord_t bead_width_x,
                             const size_t inset_count, const coord_t wall_0_inset, const PrintObjectConfig &print_object_config)
    : outline(outline)
    , bead_width_0(bead_width_0)
    , bead_width_x(bead_width_x)
    , inset_count(inset_count)
    , wall_0_inset(wall_0_inset)
    , strategy_type(print_object_config.beading_strategy_type.value)
    , print_thin_walls(Slic3r::Arachne::fill_outline_gaps)
    , min_feature_size(scaled<coord_t>(print_object_config.min_feature_size.value))
    , min_bead_width(scaled<coord_t>(print_object_config.min_bead_width.value))
    , small_area_length(static_cast<double>(bead_width_0) / 2.)
    , toolpaths_generated(false)
    , print_object_config(print_object_config)
{
}

void simplify(Polygon &thiss, const int64_t smallest_line_segment_squared, const int64_t allowed_error_distance_squared)
{
    if (thiss.size() < 3)
    {
        thiss.points.clear();
        return;
    }
    if (thiss.size() == 3)
    {
        return;
    }

    Polygon new_path;
    Point previous = thiss.points.back();
    Point previous_previous = thiss.points.at(thiss.points.size() - 2);
    Point current = thiss.points.at(0);

    /* When removing a vertex, we check the height of the triangle of the area
     being removed from the original polygon by the simplification. However,
     when consecutively removing multiple vertices the height of the previously
     removed vertices w.r.t. the shortcut path changes.
     In order to not recompute the new height value of previously removed
     vertices we compute the height of a representative triangle, which covers
     the same amount of area as the area being cut off. We use the Shoelace
     formula to accumulate the area under the removed segments. This works by
     computing the area in a 'fan' where each of the blades of the fan go from
     the origin to one of the segments. While removing vertices the area in
     this fan accumulates. By subtracting the area of the blade connected to
     the short-cutting segment we obtain the total area of the cutoff region.
     From this area we compute the height of the representative triangle using
     the standard formula for a triangle area: A = .5*b*h
     */
    int64_t accumulated_area_removed = int64_t(previous.x()) * int64_t(current.y()) - int64_t(previous.y()) * int64_t(current.x()); // Twice the Shoelace formula for area of polygon per line segment.

    for (size_t point_idx = 0; point_idx < thiss.points.size(); point_idx++)
    {
        current = thiss.points.at(point_idx % thiss.points.size());

        //Check if the accumulated area doesn't exceed the maximum.
        Point next;
        if (point_idx + 1 < thiss.points.size())
        {
            next = thiss.points.at(point_idx + 1);
        }
        else if (point_idx + 1 == thiss.points.size() && new_path.size() > 1)
        { // don't spill over if the [next] vertex will then be equal to [previous]
            next = new_path[0]; //Spill over to new polygon for checking removed area.
        }
        else
        {
            next = thiss.points.at((point_idx + 1) % thiss.points.size());
        }
        const int64_t removed_area_next = int64_t(current.x()) * int64_t(next.y()) - int64_t(current.y()) * int64_t(next.x()); // Twice the Shoelace formula for area of polygon per line segment.
        const int64_t negative_area_closing = int64_t(next.x()) * int64_t(previous.y()) - int64_t(next.y()) * int64_t(previous.x()); // area between the origin and the short-cutting segment
        accumulated_area_removed += removed_area_next;

        const int64_t length2 = (current - previous).cast<int64_t>().squaredNorm();
        if (length2 < scaled<int64_t>(25.))
        {
            // We're allowed to always delete segments of less than 5 micron.
            continue;
        }

        const int64_t area_removed_so_far = accumulated_area_removed + negative_area_closing; // close the shortcut area polygon
        const int64_t base_length_2 = (next - previous).cast<int64_t>().squaredNorm();

        if (base_length_2 == 0) //Two line segments form a line back and forth with no area.
        {
            continue; //Remove the vertex.
        }
        //We want to check if the height of the triangle formed by previous, current and next vertices is less than allowed_error_distance_squared.
        //1/2 L = A           [actual area is half of the computed shoelace value] // Shoelace formula is .5*(...) , but we simplify the computation and take out the .5
        //A = 1/2 * b * h     [triangle area formula]
        //L = b * h           [apply above two and take out the 1/2]
        //h = L / b           [divide by b]
        //h^2 = (L / b)^2     [square it]
        //h^2 = L^2 / b^2     [factor the divisor]
        const int64_t height_2 = double(area_removed_so_far) * double(area_removed_so_far) / double(base_length_2);
        if ((height_2 <= Slic3r::sqr(scaled<coord_t>(0.005)) //Almost exactly colinear (barring rounding errors).
             && Line::distance_to_infinite(current, previous, next) <= scaled<double>(0.005))) // make sure that height_2 is not small because of cancellation of positive and negative areas
        {
            continue;
        }

        if (length2 < smallest_line_segment_squared
            && height_2 <= allowed_error_distance_squared) // removing the vertex doesn't introduce too much error.)
        {
            const int64_t next_length2 = (current - next).cast<int64_t>().squaredNorm();
            if (next_length2 > smallest_line_segment_squared)
            {
                // Special case; The next line is long. If we were to remove this, it could happen that we get quite noticeable artifacts.
                // We should instead move this point to a location where both edges are kept and then remove the previous point that we wanted to keep.
                // By taking the intersection of these two lines, we get a point that preserves the direction (so it makes the corner a bit more pointy).
                // We just need to be sure that the intersection point does not introduce an artifact itself.
                Point intersection_point;
                bool has_intersection = Line(previous_previous, previous).intersection_infinite(Line(current, next), &intersection_point);
                if (!has_intersection
                    || Line::distance_to_infinite_squared(intersection_point, previous, current) > double(allowed_error_distance_squared)
                    || (intersection_point - previous).cast<int64_t>().squaredNorm() > smallest_line_segment_squared  // The intersection point is way too far from the 'previous'
                    || (intersection_point - next).cast<int64_t>().squaredNorm() > smallest_line_segment_squared)     // and 'next' points, so it shouldn't replace 'current'
                {
                    // We can't find a better spot for it, but the size of the line is more than 5 micron.
                    // So the only thing we can do here is leave it in...
                }
                else
                {
                    // New point seems like a valid one.
                    current = intersection_point;
                    // If there was a previous point added, remove it.
                    if(!new_path.empty())
                    {
                        new_path.points.pop_back();
                        previous = previous_previous;
                    }
                }
            }
            else
            {
                continue; //Remove the vertex.
            }
        }
        //Don't remove the vertex.
        accumulated_area_removed = removed_area_next; // so that in the next iteration it's the area between the origin, [previous] and [current]
        previous_previous = previous;
        previous = current; //Note that "previous" is only updated if we don't remove the vertex.
        new_path.points.push_back(current);
    }

    thiss = new_path;
}

/*!
     * Removes vertices of the polygons to make sure that they are not too high
     * resolution.
     *
     * This removes points which are connected to line segments that are shorter
     * than the `smallest_line_segment`, unless that would introduce a deviation
     * in the contour of more than `allowed_error_distance`.
     *
     * Criteria:
     * 1. Never remove a vertex if either of the connceted segments is larger than \p smallest_line_segment
     * 2. Never remove a vertex if the distance between that vertex and the final resulting polygon would be higher than \p allowed_error_distance
     * 3. The direction of segments longer than \p smallest_line_segment always
     * remains unaltered (but their end points may change if it is connected to
     * a small segment)
     *
     * Simplify uses a heuristic and doesn't neccesarily remove all removable
     * vertices under the above criteria, but simplify may never violate these
     * criteria. Unless the segments or the distance is smaller than the
     * rounding error of 5 micron.
     *
     * Vertices which introduce an error of less than 5 microns are removed
     * anyway, even if the segments are longer than the smallest line segment.
     * This makes sure that (practically) colinear line segments are joined into
     * a single line segment.
     * \param smallest_line_segment Maximal length of removed line segments.
     * \param allowed_error_distance If removing a vertex introduces a deviation
     * from the original path that is more than this distance, the vertex may
     * not be removed.
 */
void simplify(Polygons &thiss, const int64_t smallest_line_segment = scaled<coord_t>(0.01), const int64_t allowed_error_distance = scaled<coord_t>(0.005))
{
    const int64_t allowed_error_distance_squared = int64_t(allowed_error_distance) * int64_t(allowed_error_distance);
    const int64_t smallest_line_segment_squared = int64_t(smallest_line_segment) * int64_t(smallest_line_segment);
    for (size_t p = 0; p < thiss.size(); p++)
    {
        simplify(thiss[p], smallest_line_segment_squared, allowed_error_distance_squared);
        if (thiss[p].size() < 3)
        {
            thiss.erase(thiss.begin() + p);
            p--;
        }
    }
}

/*!
 * Locator to extract a line segment out of a \ref PolygonsPointIndex
 */
struct PolygonsPointIndexSegmentLocator
{
    std::pair<Point, Point> operator()(const PolygonsPointIndex &val) const
    {
        const Polygon &poly           = (*val.polygons)[val.poly_idx];
        const Point    start          = poly[val.point_idx];
        unsigned int   next_point_idx = (val.point_idx + 1) % poly.size();
        const Point    end            = poly[next_point_idx];
        return std::pair<Point, Point>(start, end);
    }
};

typedef SparseLineGrid<PolygonsPointIndex, PolygonsPointIndexSegmentLocator> LocToLineGrid;
std::unique_ptr<LocToLineGrid>                                               createLocToLineGrid(const Polygons &polygons, int square_size)
{
    unsigned int n_points = 0;
    for (const auto &poly : polygons)
        n_points += poly.size();

    auto ret = std::make_unique<LocToLineGrid>(square_size, n_points);

    for (unsigned int poly_idx = 0; poly_idx < polygons.size(); poly_idx++)
        for (unsigned int point_idx = 0; point_idx < polygons[poly_idx].size(); point_idx++)
            ret->insert(PolygonsPointIndex(&polygons, poly_idx, point_idx));
    return ret;
}

/* Note: Also tries to solve for near-self intersections, when epsilon >= 1
 */
void fixSelfIntersections(const coord_t epsilon, Polygons &thiss)
{
    if (epsilon < 1) {
        ClipperLib::SimplifyPolygons(ClipperUtils::PolygonsProvider(thiss));
        return;
    }

    const int64_t half_epsilon = (epsilon + 1) / 2;

    // Points too close to line segments should be moved a little away from those line segments, but less than epsilon,
    //   so at least half-epsilon distance between points can still be guaranteed.
    constexpr coord_t grid_size  = scaled<coord_t>(2.);
    auto              query_grid = createLocToLineGrid(thiss, grid_size);

    const auto    move_dist         = std::max<int64_t>(2L, half_epsilon - 2);
    const int64_t half_epsilon_sqrd = half_epsilon * half_epsilon;

    const size_t n = thiss.size();
    for (size_t poly_idx = 0; poly_idx < n; poly_idx++) {
        const size_t pathlen = thiss[poly_idx].size();
        for (size_t point_idx = 0; point_idx < pathlen; ++point_idx) {
            Point &pt = thiss[poly_idx][point_idx];
            for (const auto &line : query_grid->getNearby(pt, epsilon)) {
                const size_t line_next_idx = (line.point_idx + 1) % thiss[line.poly_idx].size();
                if (poly_idx == line.poly_idx && (point_idx == line.point_idx || point_idx == line_next_idx))
                    continue;

                const Line segment(thiss[line.poly_idx][line.point_idx], thiss[line.poly_idx][line_next_idx]);
                Point      segment_closest_point;
                segment.distance_to_squared(pt, &segment_closest_point);

                if (half_epsilon_sqrd >= (pt - segment_closest_point).cast<int64_t>().squaredNorm()) {
                    const Point  &other = thiss[poly_idx][(point_idx + 1) % pathlen];
                    const Vec2i64 vec   = (LinearAlg2D::pointIsLeftOfLine(other, segment.a, segment.b) > 0 ? segment.b - segment.a : segment.a - segment.b).cast<int64_t>();
                    assert(Slic3r::sqr(double(vec.x())) < double(std::numeric_limits<int64_t>::max()));
                    assert(Slic3r::sqr(double(vec.y())) < double(std::numeric_limits<int64_t>::max()));
                    const int64_t len   = vec.norm();
                    pt.x() += (-vec.y() * move_dist) / len;
                    pt.y() += (vec.x() * move_dist) / len;
                }
            }
        }
    }

    ClipperLib::SimplifyPolygons(ClipperUtils::PolygonsProvider(thiss));
}

/*!
     * Removes overlapping consecutive line segments which don't delimit a positive area.
 */
void removeDegenerateVerts(Polygons &thiss)
{
    for (unsigned int poly_idx = 0; poly_idx < thiss.size(); poly_idx++) {
        Polygon &poly = thiss[poly_idx];
        Polygon  result;

        auto isDegenerate = [](const Point &last, const Point &now, const Point &next) {
            Vec2i64 last_line = (now - last).cast<int64_t>();
            Vec2i64 next_line = (next - now).cast<int64_t>();
            return last_line.dot(next_line) == -1 * last_line.norm() * next_line.norm();
        };
        bool isChanged = false;
        for (unsigned int idx = 0; idx < poly.size(); idx++) {
            const Point &last = (result.size() == 0) ? poly.back() : result.back();
            if (idx + 1 == poly.size() && result.size() == 0) { break; }
            Point &next = (idx + 1 == poly.size()) ? result[0] : poly[idx + 1];
            if (isDegenerate(last, poly[idx], next)) { // lines are in the opposite direction
                // don't add vert to the result
                isChanged = true;
                while (result.size() > 1 && isDegenerate(result[result.size() - 2], result.back(), next)) { result.points.pop_back(); }
            } else {
                result.points.emplace_back(poly[idx]);
            }
        }

        if (isChanged) {
            if (result.size() > 2) {
                poly = result;
            } else {
                thiss.erase(thiss.begin() + poly_idx);
                poly_idx--; // effectively the next iteration has the same poly_idx (referring to a new poly which is not yet processed)
            }
        }
    }
}

void removeSmallAreas(Polygons &thiss, const double min_area_size, const bool remove_holes)
{
    auto to_path = [](const Polygon &poly) -> ClipperLib::Path {
        ClipperLib::Path out;
        for (const Point &pt : poly.points)
            out.emplace_back(ClipperLib::cInt(pt.x()), ClipperLib::cInt(pt.y()));
        return out;
    };

    auto new_end = thiss.end();
    if(remove_holes)
    {
        for(auto it = thiss.begin(); it < new_end; it++)
        {
            // All polygons smaller than target are removed by replacing them with a polygon from the back of the vector
            if(fabs(ClipperLib::Area(to_path(*it))) < min_area_size)
            {
                new_end--;
                *it = std::move(*new_end);
                it--; // wind back the iterator such that the polygon just swaped in is checked next
            }
        }
    }
    else
    {
        // For each polygon, computes the signed area, move small outlines at the end of the vector and keep pointer on small holes
        std::vector<Polygon> small_holes;
        for(auto it = thiss.begin(); it < new_end; it++) {
            double area = ClipperLib::Area(to_path(*it));
            if (fabs(area) < min_area_size)
            {
                if(area >= 0)
                {
                    new_end--;
                    if(it < new_end) {
                        std::swap(*new_end, *it);
                        it--;
                    }
                    else
                    { // Don't self-swap the last Path
                        break;
                    }
                }
                else
                {
                    small_holes.push_back(*it);
                }
            }
        }

        // Removes small holes that have their first point inside one of the removed outlines
        // Iterating in reverse ensures that unprocessed small holes won't be moved
        const auto removed_outlines_start = new_end;
        for(auto hole_it = small_holes.rbegin(); hole_it < small_holes.rend(); hole_it++)
        {
            for(auto outline_it = removed_outlines_start; outline_it < thiss.end() ; outline_it++)
            {
                if(Polygon(*outline_it).contains(*hole_it->begin())) {
                    new_end--;
                    *hole_it = std::move(*new_end);
                    break;
                }
            }
        }
    }
    thiss.resize(new_end-thiss.begin());
}

void removeColinearEdges(Polygon &poly, const double max_deviation_angle)
{
    // TODO: Can be made more efficient (for example, use pointer-types for process-/skip-indices, so we can swap them without copy).
    size_t num_removed_in_iteration = 0;
    do {
        num_removed_in_iteration = 0;
        std::vector<bool> process_indices(poly.points.size(), true);

        bool go = true;
        while (go) {
            go = false;

            const auto  &rpath   = poly;
            const size_t pathlen = rpath.size();
            if (pathlen <= 3)
                return;

            std::vector<bool> skip_indices(poly.points.size(), false);

            Polygon new_path;
            for (size_t point_idx = 0; point_idx < pathlen; ++point_idx) {
                // Don't iterate directly over process-indices, but do it this way, because there are points _in_ process-indices that should nonetheless
                // be skipped:
                if (!process_indices[point_idx]) {
                    new_path.points.push_back(rpath[point_idx]);
                    continue;
                }

                // Should skip the last point for this iteration if the old first was removed (which can be seen from the fact that the new first was skipped):
                if (point_idx == (pathlen - 1) && skip_indices[0]) {
                    skip_indices[new_path.size()] = true;
                    go                            = true;
                    new_path.points.push_back(rpath[point_idx]);
                    break;
                }

                const Point &prev = rpath[(point_idx - 1 + pathlen) % pathlen];
                const Point &pt   = rpath[point_idx];
                const Point &next = rpath[(point_idx + 1) % pathlen];

                float angle = LinearAlg2D::getAngleLeft(prev, pt, next); // [0 : 2 * pi]
                if (angle >= float(M_PI)) { angle -= float(M_PI); }                    // map [pi : 2 * pi] to [0 : pi]

                // Check if the angle is within limits for the point to 'make sense', given the maximum deviation.
                // If the angle indicates near-parallel segments ignore the point 'pt'
                if (angle > max_deviation_angle && angle < M_PI - max_deviation_angle) {
                    new_path.points.push_back(pt);
                } else if (point_idx != (pathlen - 1)) {
                    // Skip the next point, since the current one was removed:
                    skip_indices[new_path.size()] = true;
                    go                            = true;
                    new_path.points.push_back(next);
                    ++point_idx;
                }
            }
            poly = new_path;
            num_removed_in_iteration += pathlen - poly.points.size();

            process_indices.clear();
            process_indices.insert(process_indices.end(), skip_indices.begin(), skip_indices.end());
        }
    } while (num_removed_in_iteration > 0);
}

void removeColinearEdges(Polygons &thiss, const double max_deviation_angle = 0.0005)
{
    for (int p = 0; p < int(thiss.size()); p++) {
        removeColinearEdges(thiss[p], max_deviation_angle);
        if (thiss[p].size() < 3) {
            thiss.erase(thiss.begin() + p);
            p--;
        }
    }
}

const VariableWidthPaths& WallToolPaths::generate()
{
    const coord_t smallest_segment = Slic3r::Arachne::meshfix_maximum_resolution;
    const coord_t allowed_distance = Slic3r::Arachne::meshfix_maximum_deviation;
    const coord_t epsilon_offset = (allowed_distance / 2) - 1;
    const double  transitioning_angle = Geometry::deg2rad(this->print_object_config.wall_transition_angle.value);
    constexpr coord_t discretization_step_size = scaled<coord_t>(0.8);

    // Simplify outline for boost::voronoi consumption. Absolutely no self intersections or near-self intersections allowed:
    // TODO: Open question: Does this indeed fix all (or all-but-one-in-a-million) cases for manifold but otherwise possibly complex polygons?
    Polygons prepared_outline = offset(offset(offset(outline, -epsilon_offset), epsilon_offset * 2), -epsilon_offset);
    simplify(prepared_outline, smallest_segment, allowed_distance);
    fixSelfIntersections(epsilon_offset, prepared_outline);
    removeDegenerateVerts(prepared_outline);
    removeColinearEdges(prepared_outline, 0.005);
    // Removing collinear edges may introduce self intersections, so we need to fix them again
    fixSelfIntersections(epsilon_offset, prepared_outline);
    removeDegenerateVerts(prepared_outline);
    removeSmallAreas(prepared_outline, small_area_length * small_area_length, false);

    if (area(prepared_outline) > 0)
    {
        const coord_t wall_transition_length = scaled<coord_t>(this->print_object_config.wall_transition_length.value);
        const double wall_split_middle_threshold = this->print_object_config.wall_split_middle_threshold.value / 100.;  // For an uneven nr. of lines: When to split the middle wall into two.
        const double wall_add_middle_threshold = this->print_object_config.wall_add_middle_threshold.value / 100.;      // For an even nr. of lines: When to add a new middle in between the innermost two walls.
        const int wall_distribution_count = this->print_object_config.wall_distribution_count.value;
        const size_t max_bead_count = (inset_count < std::numeric_limits<coord_t>::max() / 2) ? 2 * inset_count : std::numeric_limits<coord_t>::max();
        const auto beading_strat = BeadingStrategyFactory::makeStrategy
            (
                strategy_type,
                bead_width_0,
                bead_width_x,
                wall_transition_length,
                transitioning_angle,
                print_thin_walls,
                min_bead_width,
                min_feature_size,
                wall_split_middle_threshold,
                wall_add_middle_threshold,
                max_bead_count,
                wall_0_inset,
                wall_distribution_count
            );
        const coord_t transition_filter_dist = scaled<coord_t>(this->print_object_config.wall_transition_filter_distance.value);
        SkeletalTrapezoidation wall_maker
        (
            prepared_outline,
            *beading_strat,
            beading_strat->getTransitioningAngle(),
            discretization_step_size,
            transition_filter_dist,
            wall_transition_length
        );
        wall_maker.generateToolpaths(toolpaths);
        computeInnerContour();
    }
    simplifyToolPaths(toolpaths);

    removeEmptyToolPaths(toolpaths);
    assert(std::is_sorted(toolpaths.cbegin(), toolpaths.cend(),
                          [](const VariableWidthLines& l, const VariableWidthLines& r)
                          {
                              return l.front().inset_idx < r.front().inset_idx;
                          }) && "WallToolPaths should be sorted from the outer 0th to inner_walls");
    toolpaths_generated = true;
    return toolpaths;
}

void WallToolPaths::simplifyToolPaths(VariableWidthPaths& toolpaths/*, const Settings& settings*/)
{
    for (size_t toolpaths_idx = 0; toolpaths_idx < toolpaths.size(); ++toolpaths_idx)
    {
        const int64_t maximum_resolution = Slic3r::Arachne::meshfix_maximum_resolution;
        const int64_t maximum_deviation = Slic3r::Arachne::meshfix_maximum_deviation;
        const int64_t maximum_extrusion_area_deviation = Slic3r::Arachne::meshfix_maximum_extrusion_area_deviation; // unit: μm²
        for (auto& line : toolpaths[toolpaths_idx])
        {
            line.simplify(maximum_resolution * maximum_resolution, maximum_deviation * maximum_deviation, maximum_extrusion_area_deviation);
        }
    }
}

const VariableWidthPaths& WallToolPaths::getToolPaths()
{
    if (!toolpaths_generated)
    {
        return generate();
    }
    return toolpaths;
}

void WallToolPaths::computeInnerContour()
{
    //We'll remove all 0-width paths from the original toolpaths and store them separately as polygons.
    VariableWidthPaths actual_toolpaths;
    actual_toolpaths.reserve(toolpaths.size()); //A bit too much, but the correct order of magnitude.
    VariableWidthPaths contour_paths;
    contour_paths.reserve(toolpaths.size() / inset_count);
    std::partition_copy(toolpaths.begin(), toolpaths.end(), std::back_inserter(actual_toolpaths), std::back_inserter(contour_paths),
        [](const VariableWidthLines& path)
        {
            for(const ExtrusionLine& line : path)
            {
                for(const ExtrusionJunction& junction : line.junctions)
                {
                    return junction.w != 0; //On the first actual junction, decide: If it's got 0 width, this is a contour. Otherwise it is an actual toolpath.
                }
            }
            return true; //No junctions with any vertices? Classify it as a toolpath then.
        });
    if (! actual_toolpaths.empty())
    {
        toolpaths = std::move(actual_toolpaths); //Filtered out the 0-width paths.
    }
    else
    {
        toolpaths.clear();
    }

    //Now convert the contour_paths to Polygons to denote the inner contour of the walled areas.
    inner_contour.clear();

    //We're going to have to stitch these paths since not all walls may be closed contours.
    //Since these walls have 0 width they should theoretically be closed. But there may be rounding errors.
    const coord_t minimum_line_width = bead_width_0 / 2;
    stitchContours(contour_paths, minimum_line_width, inner_contour);

    //The output walls from the skeletal trapezoidation have no known winding order, especially if they are joined together from polylines.
    //They can be in any direction, clockwise or counter-clockwise, regardless of whether the shapes are positive or negative.
    //To get a correct shape, we need to make the outside contour positive and any holes inside negative.
    //This can be done by applying the even-odd rule to the shape. This rule is not sensitive to the winding order of the polygon.
    //The even-odd rule would be incorrect if the polygon self-intersects, but that should never be generated by the skeletal trapezoidation.
    inner_contour = union_(inner_contour);
}

const Polygons& WallToolPaths::getInnerContour()
{
    if (!toolpaths_generated && inset_count > 0)
    {
        generate();
    }
    else if(inset_count == 0)
    {
        return outline;
    }
    return inner_contour;
}

bool WallToolPaths::removeEmptyToolPaths(VariableWidthPaths& toolpaths)
{
    toolpaths.erase(std::remove_if(toolpaths.begin(), toolpaths.end(), [](const VariableWidthLines& lines)
                                   {
                                       return lines.empty();
                                   }), toolpaths.end());
    return toolpaths.empty();
}

void WallToolPaths::stitchContours(const VariableWidthPaths& input, const coord_t stitch_distance, Polygons& output)
{
    // Create a bucket grid to find endpoints that are close together.
    struct ExtrusionLineStartLocator
    {
        const Point *operator()(const ExtrusionLine *line) { return &line->junctions.front().p; }
    };
    struct ExtrusionLineEndLocator
    {
        const Point *operator()(const ExtrusionLine *line) { return &line->junctions.back().p; }
    };

    // Only find endpoints closer than minimum_line_width, so we can't ever accidentally make crossing contours.
    ClosestPointInRadiusLookup<const ExtrusionLine*, ExtrusionLineStartLocator> line_starts(coord_t(stitch_distance * std::sqrt(2.)));
    ClosestPointInRadiusLookup<const ExtrusionLine*, ExtrusionLineEndLocator> line_ends(coord_t(stitch_distance * std::sqrt(2.)));

    auto get_search_bbox = [](const Point &pt, const coord_t radius) -> BoundingBox {
        const Point min_grid((pt - Point(radius, radius)) / radius);
        const Point max_grid((pt + Point(radius, radius)) / radius);
        return {min_grid * radius, (max_grid + Point(1, 1)) * radius - Point(1, 1)};
    };

    for (const VariableWidthLines &path : input) {
        for (const ExtrusionLine &line : path) {
            line_starts.insert(&line);
            line_ends.insert(&line);
        }
    }
    //Then go through all lines and construct chains of polylines if the endpoints are nearby.
    std::unordered_set<const ExtrusionLine*> processed_lines; //Track which lines were already processed to not process them twice.
    for(const VariableWidthLines& path : input)
    {
        for(const ExtrusionLine& line : path)
        {
            if(processed_lines.find(&line) != processed_lines.end()) //We already added this line before. It got added as a nearby line.
            {
                continue;
            }
            //We'll create a chain of polylines that get joined together. We can add polylines on both ends!
            std::deque<const ExtrusionLine*> chain;
            std::deque<bool> is_reversed; //Lines could need to be inserted in reverse. Must coincide with the `chain` deque.
            const ExtrusionLine* nearest = &line; //At every iteration, add the polyline that joins together most closely.
            bool nearest_reverse = false; //Whether the next line to insert must be inserted in reverse.
            bool nearest_before = false; //Whether the next line to insert must be inserted in the front of the chain.
            while(nearest)
            {
                if(processed_lines.find(nearest) != processed_lines.end())
                {
                    break; //Looping. This contour is already processed.
                }
                processed_lines.insert(nearest);
                if(nearest_before)
                {
                    chain.push_front(nearest);
                    is_reversed.push_front(nearest_reverse);
                }
                else
                {
                    chain.push_back(nearest);
                    is_reversed.push_back(nearest_reverse);
                }

                //Find any nearby lines to attach. Look on both ends of our current chain and find both ends of polylines.
                const Point chain_start = is_reversed.front() ? chain.front()->junctions.back().p : chain.front()->junctions.front().p;
                const Point chain_end = is_reversed.back() ? chain.back()->junctions.front().p : chain.back()->junctions.back().p;

                std::vector<std::pair<const ExtrusionLine *const *, double>> starts_near_start = line_starts.find_all(chain_start);
                std::vector<std::pair<const ExtrusionLine *const *, double>> ends_near_start   = line_ends.find_all(chain_start);
                std::vector<std::pair<const ExtrusionLine *const *, double>> starts_near_end   = line_starts.find_all(chain_end);
                std::vector<std::pair<const ExtrusionLine *const *, double>> ends_near_end     = line_ends.find_all(chain_end);

                nearest = nullptr;
                int64_t nearest_dist2 = std::numeric_limits<int64_t>::max();
                for (const auto &candidate_ptr : starts_near_start) {
                    const ExtrusionLine* candidate = *candidate_ptr.first;
                    if(processed_lines.find(candidate) != processed_lines.end())
                        continue; //Already processed this line before. It's linked to something else.

                    if (const int64_t dist2 = (candidate->junctions.front().p - chain_start).cast<int64_t>().squaredNorm(); dist2 < nearest_dist2) {
                        nearest         = candidate;
                        nearest_dist2   = dist2;
                        nearest_reverse = true;
                        nearest_before  = true;
                    }
                }
                for (const auto &candidate_ptr : ends_near_start) {
                    const ExtrusionLine* candidate = *candidate_ptr.first;
                    if(processed_lines.find(candidate) != processed_lines.end())
                        continue;

                    if (const int64_t dist2 = (candidate->junctions.back().p - chain_start).cast<int64_t>().squaredNorm(); dist2 < nearest_dist2) {
                        nearest         = candidate;
                        nearest_dist2   = dist2;
                        nearest_reverse = false;
                        nearest_before  = true;
                    }
                }
                for (const auto &candidate_ptr : starts_near_end) {
                    const ExtrusionLine* candidate = *candidate_ptr.first;
                    if(processed_lines.find(candidate) != processed_lines.end())
                        continue; //Already processed this line before. It's linked to something else.

                    if (const int64_t dist2 = (candidate->junctions.front().p - chain_start).cast<int64_t>().squaredNorm(); dist2 < nearest_dist2) {
                        nearest         = candidate;
                        nearest_dist2   = dist2;
                        nearest_reverse = false;
                        nearest_before  = false;
                    }
                }
                for (const auto &candidate_ptr : ends_near_end) {
                    const ExtrusionLine* candidate = *candidate_ptr.first;
                    if (processed_lines.find(candidate) != processed_lines.end())
                        continue;

                    if (const int64_t dist2 = (candidate->junctions.back().p - chain_start).cast<int64_t>().squaredNorm(); dist2 < nearest_dist2) {
                        nearest         = candidate;
                        nearest_dist2   = dist2;
                        nearest_reverse = true;
                        nearest_before  = false;
                    }
                }
            }

            //Now serialize the entire chain into one polygon.
            output.emplace_back();
            for (size_t i = 0; i < chain.size(); ++i) {
                if(!is_reversed[i])
                    for (const ExtrusionJunction& junction : chain[i]->junctions)
                        output.back().points.emplace_back(junction.p);
                else
                    for (auto junction = chain[i]->junctions.rbegin(); junction != chain[i]->junctions.rend(); ++junction)
                        output.back().points.emplace_back(junction->p);
            }
        }
    }
}

size_t getOuterRegionId(const Arachne::VariableWidthPaths& toolpaths, size_t& out_max_region_id)
{
    // Polygons show up here one by one, so there are always only a) the outer lines and b) the lines that are part of the holes.
    // Therefore, the outer-regions' lines will always have the region-id that is larger then all of the other ones.

    // First, build the bounding boxes:
    std::map<size_t, BoundingBox> region_ids_to_bboxes; // Use a sorted map, ordered by region_id, so that we can find the largest region_id quickly.
    for (const Arachne::VariableWidthLines &path : toolpaths) {
        for (const Arachne::ExtrusionLine &line : path) {
            BoundingBox &aabb =
                region_ids_to_bboxes[line.region_id]; // Empty AABBs are default initialized when region_ids are encountered for the first time.
            for (const auto &junction : line.junctions) aabb.merge(junction.p);
        }
    }

    // Then, the largest of these will be the one that's needed for the outer region, the others' all belong to hole regions:
    BoundingBox outer_bbox;
    size_t      outer_region_id = 0; // Region-ID 0 is reserved for 'None'.
    for (const auto &region_id_bbox_pair : region_ids_to_bboxes) {
        if (region_id_bbox_pair.second.contains(outer_bbox)) {
            outer_bbox      = region_id_bbox_pair.second;
            outer_region_id = region_id_bbox_pair.first;
        }
    }

    // Maximum Region-ID (using the ordering of the map)
    out_max_region_id = region_ids_to_bboxes.empty() ? 0 : region_ids_to_bboxes.rbegin()->first;
    return outer_region_id;
}

Arachne::BinJunctions variableWidthPathToBinJunctions(const Arachne::VariableWidthPaths& toolpaths, const bool pack_regions_by_inset, const bool center_last, std::set<size_t>* p_bins_with_index_zero_insets)
{
    // Find the largest inset-index:
    size_t max_inset_index = 0;
    for (const Arachne::VariableWidthLines &path : toolpaths)
        max_inset_index = std::max(path.front().inset_idx, max_inset_index);

    // Find which regions are associated with the outer-outer walls (which region is the one the rest is holes inside of):
    size_t       max_region_id   = 0;
    const size_t outer_region_id = getOuterRegionId(toolpaths, max_region_id);

    //Since we're (optionally!) splitting off in the outer and inner regions, it may need twice as many bins as inset-indices.
    //Add two extra bins for the center-paths, if they need to be stored separately. One bin for inner and one for outer walls.
    const size_t max_bin = (pack_regions_by_inset ? (max_region_id * 2) + 2 : (max_inset_index + 1) * 2) + center_last * 2;
    Arachne::BinJunctions insets(max_bin + 1);
    for (const Arachne::VariableWidthLines &path : toolpaths) {
        if (path.empty()) // Don't bother printing these.
            continue;

        const size_t inset_index = path.front().inset_idx;

        // Convert list of extrusion lines to vectors of extrusion junctions, and add those to the binned insets.
        for (const Arachne::ExtrusionLine &line : path) {
            // Sort into the right bin, ...
            size_t     bin_index;
            const bool in_hole_region = line.region_id != outer_region_id && line.region_id != 0;
            if (center_last && line.is_odd) {
                bin_index = inset_index > 0;
            } else if (pack_regions_by_inset) {
                bin_index = std::min(inset_index, static_cast<size_t>(1)) + 2 * (in_hole_region ? line.region_id : 0) + center_last * 2;
            } else {
                bin_index = inset_index + (in_hole_region ? (max_inset_index + 1) : 0) + center_last * 2;
            }
            insets[bin_index].emplace_back(line.junctions.begin(), line.junctions.end());

            // Collect all bins that have zero-inset indices in them, if needed:
            if (inset_index == 0 && p_bins_with_index_zero_insets != nullptr)
                p_bins_with_index_zero_insets->insert(bin_index);
        }
    }
    return insets;
}

BinJunctions WallToolPaths::getBinJunctions(std::set<size_t> &bins_with_index_zero_insets)
{
    if (!toolpaths_generated)
        generate();

    return variableWidthPathToBinJunctions(toolpaths, true, true, &bins_with_index_zero_insets);
}

} // namespace Slic3r::Arachne
