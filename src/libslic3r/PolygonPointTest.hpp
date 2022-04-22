#ifndef SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_
#define SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_

#include "libslic3r/Point.hpp"
#include "libslic3r/EdgeGrid.hpp"

namespace Slic3r {

struct EdgeGridWrapper {
    EdgeGridWrapper(coord_t edge_width, std::vector<Points> lines) :
        lines(lines), edge_width(edge_width) {

        grid.create(this->lines, edge_width, true);
        grid.calculate_sdf();
    }

    bool signed_distance(const Point &point, coordf_t point_width, coordf_t &dist_out) const {
        coordf_t tmp_dist_out;
        bool found = grid.signed_distance(point, edge_width, tmp_dist_out);
        dist_out = tmp_dist_out - edge_width / 2 - point_width / 2;
        return found;

    }

    EdgeGrid::Grid grid;
    std::vector<Points> lines;
    coord_t edge_width;
};

namespace TODO {

class PolygonPointTest {

    struct Segment {
        coord_t start;
        std::vector<size_t> lines;
    };

    std::vector<Segment> x_coord_segments;

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

        std::vector<std::pair<size_t, bool>> sweeping_data(lines.size() * 2);
        for (size_t line_index = 0; line_index < lines.size(); ++line_index) {
            sweeping_data[line_index].first = line_index;
            sweeping_data[line_index].second = true;
            sweeping_data[line_index * 2 + 1].first = line_index;
            sweeping_data[line_index * 2 + 1].second = false;
        }

        const auto data_comparator = [&lines](const std::pair<size_t, bool> &left,
                const std::pair<size_t, bool> &right) {
            const auto left_x =
                    left.second ?
                                  std::min(lines[left.first].a.x(), lines[left.first].b.x()) :
                                  std::max(lines[left.first].a.x(), lines[left.first].b.x());
            const auto right_x =
                    right.second ?
                                   std::min(lines[right.first].a.x(), lines[right.first].b.x()) :
                                   std::max(lines[right.first].a.x(), lines[right.first].b.x());

            return left_x < right_x;
        };

        std::sort(sweeping_data.begin(), sweeping_data.end(), data_comparator);
        std::set<size_t> active_lines;

        //TODO continue

    }

};
}

}

#endif /* SRC_LIBSLIC3R_POLYGONPOINTTEST_HPP_ */
