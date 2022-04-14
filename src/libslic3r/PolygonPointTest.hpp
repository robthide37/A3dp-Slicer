#ifndef SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_
#define SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_

#include "libslic3r/Point.hpp"
#include "libslic3r/EdgeGrid.hpp"

namespace Slic3r {

struct EdgeGridWrapper {
    EdgeGridWrapper(coord_t resolution, ExPolygons ex_polys) :
            ex_polys(ex_polys) {

        grid.create(this->ex_polys, resolution);
        grid.calculate_sdf();
    }

    bool signed_distance(const Point &point, coordf_t point_width, coordf_t &dist_out) const {
        coordf_t tmp_dist_out;
        bool found = grid.signed_distance(point, point_width, tmp_dist_out);
        // decrease the distance by half of edge width of previous layer and half of flow width of current layer
        dist_out = tmp_dist_out - point_width / 2;
        return found;

    }

    EdgeGrid::Grid grid;
    ExPolygons ex_polys;
};

class PolygonPointTest {
public:
    PolygonPointTest(const ExPolygons &ex_polygons) {
        std::vector<Line> lines;
        for (const auto &exp : ex_polygons) {
            Lines contour = exp.contour.lines();
            lines.insert(lines.end(), contour.begin(), contour.end());
            for (const auto &hole : exp.holes) {
                Lines hole_lines = hole.lines();
                for (Line &line : hole_lines) {
                    line.reverse(); // reverse hole lines, so that we can use normal to deduce where the object is
                }
                lines.insert(lines.end(), hole_lines.begin(), hole_lines.end());
            }
        }

        std::vector<std::pair<size_t, bool>> sweeping_data(lines.size());
        sweeping_data.reserve(lines.size() * 2);
        for (int line_index = 0; line_index < lines.size(); ++line_index) {
            sweeping_data[line_index].first = line_index;
            sweeping_data[line_index].second = true;
        }

        const auto data_comparator = [&lines](const std::pair<size_t, bool> &left,
                const std::pair<size_t, bool> &right) {
            return std::min(lines[left.first].a.x(), lines[left.first].b.x())
                    < std::min(lines[right.first].a.x(), lines[right.first].b.x());
        };

        std::make_heap(sweeping_data.begin(), sweeping_data.end(), data_comparator);
        std::set<size_t> active_lines;

        //TODO continue




    }

};

}

#endif /* SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_ */
