/**
 * @file
 * @brief Utility functions for layer change gcode generation.
 */

#ifndef slic3r_GCode_LayerChanges_hpp_
#define slic3r_GCode_LayerChanges_hpp_

#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"

namespace Slic3r::GCode::Impl::LayerChanges {
/**
 * Generates a regular polygon - all angles are the same (e.g. typical hexagon).
 *
 * @param centroid Central point.
 * @param start_point The polygon point are ordered. This is the first point.
 * @param points_count Amount of nodes of the polygon (e.g. 6 for haxagon).
 *
 * Distance between centroid and start point sets the scale of the polygon.
 */
Polygon generate_regular_polygon(
    const Point &centroid, const Point &start_point, const unsigned points_count
);

/**
 * @brief A representation of the bed shape with inner padding.
 *
 * Its purpose is to facilitate the bed boundary checking.
 */
class Bed
{
private:
    Polygon inner_offset;
    static Polygon get_inner_offset(const std::vector<Vec2d> &shape, const double padding);

public:
    /**
     * Bed shape with inner padding.
     */
    Bed(const std::vector<Vec2d> &shape, const double padding);

    Vec2d centroid;

    /**
     * Returns true if the point is within the bed shape including inner padding.
     */
    bool contains_within_padding(const Vec2d &point) const;
};
} // namespace Slic3r::GCode::Impl::LayerChanges

#endif // slic3r_GCode_LayerChanges_hpp_
