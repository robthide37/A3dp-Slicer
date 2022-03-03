//Copyright (c) 2020 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.


#ifndef UTILS_EXTRUSION_LINE_H
#define UTILS_EXTRUSION_LINE_H

#include "ExtrusionJunction.hpp"
#include "../../Polyline.hpp"
#include "../../Polygon.hpp"

namespace Slic3r {
class ThickPolyline;
}

namespace Slic3r::Arachne
{

/*!
 * Represents a polyline (not just a line) that is to be extruded with variable
 * line width.
 *
 * This polyline is a sequence of \ref ExtrusionJunction, with a bit of metadata
 * about which inset it represents.
 */
struct ExtrusionLine
{
    /*!
     * Which inset this path represents, counted from the outside inwards.
     *
     * The outer wall has index 0.
     */
    size_t inset_idx;

    /*!
     * If a thin piece needs to be printed with an odd number of walls (e.g. 5
     * walls) then there will be one wall in the middle that is not a loop. This
     * field indicates whether this path is such a line through the middle, that
     * has no companion line going back on the other side and is not a closed
     * loop.
     */
    bool is_odd;

    /*!
     * Which region this line is part of. A solid polygon without holes has only one region.
     * A polygon with holes has 2. Disconnected parts of the polygon are also separate regions.
     * Will be 0 if no region was given.
     */
    size_t region_id;

    /*!
     * The list of vertices along which this path runs.
     *
     * Each junction has a width, making this path a variable-width path.
     */
    std::vector<ExtrusionJunction> junctions;

    ExtrusionLine(const size_t inset_idx, const bool is_odd, const size_t region_id = 0);

    /*!
     * Sum the total length of this path.
     */
    coord_t getLength() const;

    /*!
     * Get the minimal width of this path
     */
    coord_t getMinimalWidth() const;

    /*!
     * Export the included junctions as vector.
     */
    void appendJunctionsTo(LineJunctions& result) const;

    /*!
     * Removes vertices of the ExtrusionLines to make sure that they are not too high
     * resolution.
     *
     * This removes junctions which are connected to line segments that are shorter
     * than the `smallest_line_segment`, unless that would introduce a deviation
     * in the contour of more than `allowed_error_distance`.
     *
     * Criteria:
     * 1. Never remove a junction if either of the connected segments is larger than \p smallest_line_segment
     * 2. Never remove a junction if the distance between that junction and the final resulting polygon would be higher
     *    than \p allowed_error_distance
     * 3. The direction of segments longer than \p smallest_line_segment always
     *    remains unaltered (but their end points may change if it is connected to
     *    a small segment)
     * 4. Never remove a junction if it has a distinctively different width than the next junction, as this can
     *    introduce unwanted irregularities on the wall widths.
     *
     * Simplify uses a heuristic and doesn't necessarily remove all removable
     * vertices under the above criteria, but simplify may never violate these
     * criteria. Unless the segments or the distance is smaller than the
     * rounding error of 5 micron.
     *
     * Vertices which introduce an error of less than 5 microns are removed
     * anyway, even if the segments are longer than the smallest line segment.
     * This makes sure that (practically) co-linear line segments are joined into
     * a single line segment.
     * \param smallest_line_segment Maximal length of removed line segments.
     * \param allowed_error_distance If removing a vertex introduces a deviation
     *         from the original path that is more than this distance, the vertex may
     *         not be removed.
     * \param maximum_extrusion_area_deviation The maximum extrusion area deviation allowed when removing intermediate
     *        junctions from a straight ExtrusionLine
     */
    void simplify(int64_t smallest_line_segment_squared, int64_t allowed_error_distance_squared, int64_t maximum_extrusion_area_deviation);

    /*!
     * Computes and returns the total area error (in μm²) of the AB and BC segments of an ABC straight ExtrusionLine
     * when the junction B with a width B.w is removed from the ExtrusionLine. The area changes due to the fact that the
     * new simplified line AC has a uniform width which equals to the weighted average of the width of the subsegments
     * (based on their length).
     *
     * \param A Start point of the 3-point-straight line
     * \param B Intermediate point of the 3-point-straight line
     * \param C End point of the 3-point-straight line
     * \param weighted_average_width The weighted average of the widths of the two colinear extrusion segments
     * */
    static int64_t calculateExtrusionAreaDeviationError(ExtrusionJunction A, ExtrusionJunction B, ExtrusionJunction C, coord_t& weighted_average_width);
};

using VariableWidthLines = std::vector<ExtrusionLine>; //<! The ExtrusionLines generated by libArachne for each Path
using VariableWidthPaths = std::vector<VariableWidthLines>; //<! The toolpaths generated by libArachne

static inline Slic3r::ThickPolyline to_thick_polyline(const Arachne::LineJunctions &line_junctions)
{
    assert(line_junctions.size() >= 2);
    Slic3r::ThickPolyline out;
    out.points.emplace_back(line_junctions.front().p);
    out.width.emplace_back(line_junctions.front().w);
    out.points.emplace_back(line_junctions[1].p);
    out.width.emplace_back(line_junctions[1].w);

    auto it_prev = line_junctions.begin() + 1;
    for (auto it = line_junctions.begin() + 2; it != line_junctions.end(); ++it) {
        out.points.emplace_back(it->p);
        out.width.emplace_back(it_prev->w);
        out.width.emplace_back(it->w);
        it_prev = it;
    }

    return out;
}

static inline Polygon to_polygon(const ExtrusionLine& line) {
    Polygon out;
    assert(line.junctions.size() >= 3);
    assert(line.junctions.front().p == line.junctions.back().p);
    out.points.reserve(line.junctions.size() - 1);
    for (auto it = line.junctions.begin(); it != line.junctions.end() - 1; ++it)
        out.points.emplace_back(it->p);
    return out;
}

} // namespace Slic3r::Arachne
#endif // UTILS_EXTRUSION_LINE_H
