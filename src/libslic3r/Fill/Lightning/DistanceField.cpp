//Copyright (c) 2021 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include "DistanceField.hpp" //Class we're implementing.
#include "../FillRectilinear.hpp"
#include "../../ClipperUtils.hpp"

namespace Slic3r::FillLightning
{

constexpr coord_t radius_per_cell_size = 6;  // The cell-size should be small compared to the radius, but not so small as to be inefficient.

DistanceField::DistanceField(const coord_t& radius, const Polygons& current_outline, const BoundingBox& current_outlines_bbox, const Polygons& current_overhang) :
    m_cell_size(radius / radius_per_cell_size),
    m_supporting_radius(radius),
    m_unsupported_points_bbox(current_outlines_bbox)
{
    m_supporting_radius2 = Slic3r::sqr(int64_t(radius));
    // Sample source polygons with a regular grid sampling pattern.
    for (const ExPolygon &expoly : union_ex(current_overhang)) {
        for (const Point &p : sample_grid_pattern(expoly, m_cell_size)) {
            // Find a squared distance to the source expolygon boundary.
            double d2 = std::numeric_limits<double>::max();
            for (size_t icontour = 0; icontour <= expoly.holes.size(); ++icontour) {
                const Polygon &contour = icontour == 0 ? expoly.contour : expoly.holes[icontour - 1];
                if (contour.size() > 2) {
                    Point prev = contour.points.back();
                    for (const Point &p2 : contour.points) {
                        d2   = std::min(d2, Line::distance_to_squared(p, prev, p2));
                        prev = p2;
                    }
                }
            }
            m_unsupported_points.emplace_back(p, sqrt(d2));
            assert(m_unsupported_points_bbox.contains(p));
        }
    }
    m_unsupported_points.sort([&radius](const UnsupportedCell &a, const UnsupportedCell &b) {
        constexpr coord_t prime_for_hash = 191;
        return std::abs(b.dist_to_boundary - a.dist_to_boundary) > radius ?
               a.dist_to_boundary < b.dist_to_boundary :
               (PointHash{}(a.loc) % prime_for_hash) < (PointHash{}(b.loc) % prime_for_hash);
        });
    for (auto it = m_unsupported_points.begin(); it != m_unsupported_points.end(); ++it) {
        UnsupportedCell& cell = *it;
        m_unsupported_points_grid.emplace(this->to_grid_point(cell.loc), it);
    }
    // Because the distance between two points is at least one axis equal to m_cell_size, every cell
    // in m_unsupported_points_grid contains exactly one point.
    assert(m_unsupported_points.size() == m_unsupported_points_grid.size());
}

void DistanceField::update(const Point& to_node, const Point& added_leaf)
{
    Vec2d       v  = (added_leaf - to_node).cast<double>();
    auto        l2 = v.squaredNorm();
    Vec2d       extent = Vec2d(-v.y(), v.x()) * m_supporting_radius / sqrt(l2);

    BoundingBox grid;
    {
        Point diagonal(m_supporting_radius, m_supporting_radius);
        Point iextent(extent.cast<coord_t>());
        grid = BoundingBox(added_leaf - diagonal, added_leaf + diagonal);
        grid.merge(to_node - iextent);
        grid.merge(to_node + iextent);
        grid.merge(added_leaf - iextent);
        grid.merge(added_leaf + iextent);

        // Clip grid by m_unsupported_points_bbox. Mainly to ensure that grid.min is a non-negative value.
        grid.min.x() = std::max(grid.min.x(), m_unsupported_points_bbox.min.x());
        grid.min.y() = std::max(grid.min.y(), m_unsupported_points_bbox.min.y());
        grid.max.x() = std::min(grid.max.x(), m_unsupported_points_bbox.max.x());
        grid.max.y() = std::min(grid.max.y(), m_unsupported_points_bbox.max.y());

        grid.min = this->to_grid_point(grid.min);
        grid.max = this->to_grid_point(grid.max);
    }

    Point grid_addr;
    Point grid_loc;
    for (grid_addr.y() = grid.min.y(); grid_addr.y() <= grid.max.y(); ++grid_addr.y()) {
        for (grid_addr.x() = grid.min.x(); grid_addr.x() <= grid.max.x(); ++grid_addr.x()) {
            grid_loc = this->from_grid_point(grid_addr);
            // Test inside a circle at the new leaf.
            if ((grid_loc - added_leaf).cast<int64_t>().squaredNorm() > m_supporting_radius2) {
                // Not inside a circle at the end of the new leaf.
                // Test inside a rotated rectangle.
                Vec2d  vx = (grid_loc - to_node).cast<double>();
                double d  = v.dot(vx);
                if (d >= 0 && d <= l2) {
                    d = extent.dot(vx);
                    if (d < -1. || d > 1.)
                        // Not inside a rotated rectangle.
                        continue;
                }
            }
            // Inside a circle at the end of the new leaf, or inside a rotated rectangle.
            // Remove unsupported leafs at this grid location.
            if (auto it = m_unsupported_points_grid.find(grid_addr); it != m_unsupported_points_grid.end()) {
                std::list<UnsupportedCell>::iterator& list_it = it->second;
                UnsupportedCell& cell = *list_it;
                if ((cell.loc - added_leaf).cast<int64_t>().squaredNorm() <= m_supporting_radius2) {
                    m_unsupported_points.erase(list_it);
                    m_unsupported_points_grid.erase(it);
                }
            }
        }
    }
}

#if 0
void DistanceField::update(const Point &to_node, const Point &added_leaf)
{
    const Point supporting_radius_point(m_supporting_radius, m_supporting_radius);
    const BoundingBox grid(this->to_grid_point(added_leaf - supporting_radius_point), this->to_grid_point(added_leaf + supporting_radius_point));

    for (coord_t grid_y = grid.min.y(); grid_y <= grid.max.y(); ++grid_y) {
        for (coord_t grid_x = grid.min.x(); grid_x <= grid.max.x(); ++grid_x) {
            if (auto it = m_unsupported_points_grid.find({grid_x, grid_y}); it != m_unsupported_points_grid.end()) {
                std::list<UnsupportedCell>::iterator &list_it = it->second;
                UnsupportedCell                      &cell    = *list_it;
                if ((cell.loc - added_leaf).cast<int64_t>().squaredNorm() <= m_supporting_radius2) {
                    m_unsupported_points.erase(list_it);
                    m_unsupported_points_grid.erase(it);
                }
            }
        }
    }
}
#endif

} // namespace Slic3r::FillLightning
