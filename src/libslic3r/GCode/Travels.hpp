/**
 * @file
 * @brief Utility functions for travel gcode generation.
 */

#ifndef slic3r_GCode_Travels_hpp_
#define slic3r_GCode_Travels_hpp_

#include <vector>
#include <tcbspan/span.hpp>
#include <functional>
#include <optional>

#include "libslic3r/Line.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/PrintConfig.hpp"

namespace Slic3r::GCode::Impl::Travels {
struct DistancedPoint
{
    Point point;
    double distance_from_start;
};

/**
 * @brief Takes a path described as a list of points and adds points to it.
 *
 * @param xy_path A list of points describing a path in xy.
 * @param sorted_distances A sorted list of distances along the path.
 * @return Sliced path.
 *
 * The algorithm travels along the path segments and adds points to
 * the segments in such a way that the points have specified distances
 * from the xy_path start. **Any distances over the xy_path end will
 * be simply ignored.**
 *
 * Example usage - simplified for clarity:
 * @code
 * std::vector<double> distances{0.5, 1.5};
 * std::vector<Points> xy_path{{0, 0}, {1, 0}};
 * // produces
 * {{0, 0}, {0, 0.5}, {1, 0}}
 * // notice that 1.5 is omitted
 * @endcode
 */
std::vector<DistancedPoint> slice_xy_path(
    tcb::span<const Point> xy_path, tcb::span<const double> sorted_distances
);

/**
 * @brief Simply return the xy_path with z coord set to elevation.
 */
Points3 generate_flat_travel(tcb::span<const Point> xy_path, const float elevation);

/**
 * @brief Take xy_path and genrate a travel acording to elevation.
 *
 * @param xy_path A list of points describing a path in xy.
 * @param ensure_points_at_distances See slice_xy_path sorted_distances.
 * @param elevation  A function taking current distance in mm as input and returning elevation in mm
 * as output.
 *
 * **Be aweare** that the elevation function operates in mm, while xy_path and returned travel are
 * in scaled coordinates.
 */
Points3 generate_elevated_travel(
    const tcb::span<const Point> xy_path,
    const std::vector<double> &ensure_points_at_distances,
    const double initial_elevation,
    const std::function<double(double)> &elevation
);

/**
 * @brief Given a AABB tree over lines find intersection with xy_path closest to the xy_path start.
 *
 * @param xy_path A path in 2D.
 * @param distancer AABB Tree over lines.
 * @return Distance to the first intersection if there is one.
 *
 * **Ignores intersection with xy_path starting point.**
 */
std::optional<double> get_first_crossed_line_distance(
    tcb::span<const Line> xy_path, const AABBTreeLines::LinesDistancer<Linef> &distancer
);

/**
 * @brief Extract parameters and decide wheather the travel can be elevated.
 * Then generate the whole travel 3D path - elevated if possible.
 */
Points3 generate_travel_to_extrusion(
    const Polyline &xy_path,
    const FullPrintConfig &config,
    const unsigned extruder_id,
    const double initial_elevation,
    const std::optional<AABBTreeLines::LinesDistancer<Linef>> &previous_layer_distancer,
    const Point &xy_path_coord_origin
);
} // namespace Slic3r::GCode::Impl::Travels

#endif // slic3r_GCode_Travels_hpp_
