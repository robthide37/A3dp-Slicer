#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Arachne/WallToolPaths.hpp"

#include "AABBTreeLines.hpp"
#include "ExPolygon.hpp"
#include "FillEnsuring.hpp"
#include "Line.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "SVG.hpp"
#include "libslic3r.h"

#include <algorithm>
#include <boost/log/trivial.hpp>
#include <functional>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace Slic3r {

ThickPolylines FillEnsuring::fill_surface_arachne(const Surface *surface, const FillParams &params)
{
    assert(params.use_arachne);
    assert(this->print_config != nullptr && this->print_object_config != nullptr && this->print_region_config != nullptr);

    auto rotate_thick_polylines = [](ThickPolylines &tpolylines, double cos_angle, double sin_angle) {
        for (ThickPolyline &tp : tpolylines) {
            for (auto &p : tp.points) {
                double px = double(p.x());
                double py = double(p.y());
                p.x()     = coord_t(round(cos_angle * px - sin_angle * py));
                p.y()     = coord_t(round(cos_angle * py + sin_angle * px));
            }
        }
    };

    auto segments_overlap = [](coord_t alow, coord_t ahigh, coord_t blow, coord_t bhigh) {
        return (alow >= blow && alow <= bhigh) || (ahigh >= blow && ahigh <= bhigh) || (blow >= alow && blow <= ahigh) ||
               (bhigh >= alow && bhigh <= ahigh);
    };

    const coord_t scaled_spacing                      = scaled<coord_t>(this->spacing);
    double        squared_distance_limit_reconnection = 4 * double(scaled_spacing) * double(scaled_spacing);
    Polygons      filled_area                         = to_polygons(surface->expolygon);
    double        aligning_angle                      = -this->angle + PI * 0.5;
    polygons_rotate(filled_area, aligning_angle);
    Polygons    internal_area = shrink(filled_area, 0.5 * scaled_spacing - scale_(this->overlap));
    BoundingBox bb            = get_extents(filled_area);

    const size_t      n_vlines = (bb.max.x() - bb.min.x() + scaled_spacing - 1) / scaled_spacing;
    std::vector<Line> vertical_lines(2 * n_vlines + 1);
    coord_t           y_min = bb.min.y();
    coord_t           y_max = bb.max.y();
    for (size_t i = 0; i < n_vlines; i++) {
        coord_t x0                  = bb.min.x() + i * double(scaled_spacing) - scaled_spacing * 0.5;
        coord_t x1                  = bb.min.x() + i * double(scaled_spacing);
        vertical_lines[i * 2].a     = Point{x0, y_min};
        vertical_lines[i * 2].b     = Point{x0, y_max};
        vertical_lines[i * 2 + 1].a = Point{x1, y_min};
        vertical_lines[i * 2 + 1].b = Point{x1, y_max};
    }
    vertical_lines.back().a = Point{coord_t(bb.min.x() + n_vlines * double(scaled_spacing) + scaled_spacing * 0.5), y_min};
    vertical_lines.back().b = Point{vertical_lines.back().a.x(), y_max};

    auto                                                         area_walls = AABBTreeLines::LinesDistancer<Line>{to_lines(internal_area)};
    std::vector<std::vector<std::pair<Vec<2, coord_t>, size_t>>> vertical_lines_intersections(vertical_lines.size());
    for (int i = 0; i < vertical_lines.size(); i++) {
        vertical_lines_intersections[i] = area_walls.intersections_with_line<true>(vertical_lines[i]);
    }
    std::vector<std::vector<Line>> polygon_sections(n_vlines);

    for (size_t i = 0; i < n_vlines; i++) {
        const auto &central_intersections = vertical_lines_intersections[i * 2 + 1];
        const auto &left_intersections    = vertical_lines_intersections[i * 2];
        const auto &right_intersections   = vertical_lines_intersections[i * 2 + 2];

        for (int intersection_idx = 0; intersection_idx < int(central_intersections.size()) - 1; intersection_idx++) {
            const auto &a = central_intersections[intersection_idx];
            const auto &b = central_intersections[intersection_idx + 1];
            if (area_walls.outside((a.first + b.first) / 2) < 0) {
                // central part is inside. Now check for reasonable side distances

                auto get_closest_intersection_squared_dist =
                    [](const std::pair<Vec<2, coord_t>, size_t>              &point,
                       const std::vector<std::pair<Vec<2, coord_t>, size_t>> &sorted_intersections) {
                        if (sorted_intersections.empty()) {
                            return 0.0;
                        }
                        auto closest_higher = std::upper_bound(sorted_intersections.begin(), sorted_intersections.end(), point,
                                                               [](const std::pair<Vec<2, coord_t>, size_t> &left,
                                                                  const std::pair<Vec<2, coord_t>, size_t> &right) {
                                                                   return left.first.y() < right.first.y();
                                                               });
                        if (closest_higher == sorted_intersections.end()) {
                            return (point.first - sorted_intersections.back().first).cast<double>().squaredNorm();
                        }
                        double candidate_dist = (point.first - closest_higher->first).cast<double>().squaredNorm();
                        if (closest_higher != sorted_intersections.begin()) {
                            double closest_lower_dist = (point.first - (closest_higher--)->first).cast<double>().squaredNorm();
                            candidate_dist            = std::min(candidate_dist, closest_lower_dist);
                        }
                        return candidate_dist;
                    };
                Point section_a = a.first;
                Point section_b = b.first;

                double max_a_squared_dist = std::max(get_closest_intersection_squared_dist(a, left_intersections),
                                             get_closest_intersection_squared_dist(a, right_intersections));

                double max_b_squared_dist = std::max(get_closest_intersection_squared_dist(b, left_intersections),
                                             get_closest_intersection_squared_dist(b, right_intersections));

                if (max_a_squared_dist > 0.4 * squared_distance_limit_reconnection) {
                    section_a.y() += std::min(2.0 * scaled_spacing, sqrt(max_a_squared_dist));
                }

                if (max_b_squared_dist > 0.4 * squared_distance_limit_reconnection) {
                    section_b.y() -= std::min(2.0 * scaled_spacing, sqrt(max_b_squared_dist));
                }

                section_a.y() = std::min(section_a.y(), section_b.y());
                section_b.y() = std::max(section_a.y(), section_b.y());
                polygon_sections[i].emplace_back(section_a, section_b);
            }
        }
    }

    struct Node
    {
        int section_idx;
        int line_idx;
        int skips_taken = 0;
        bool neighbours_explored = false;
        std::vector<std::pair<int,int>> neighbours{};
    };

    coord_t length_filter = scale_(4);
    size_t skips_allowed = 2;
    size_t min_removal_conut = 4;
    for (int section_idx = 0; section_idx < polygon_sections.size(); section_idx++) {
        for (int line_idx = 0; line_idx < polygon_sections[section_idx].size(); line_idx++) {
            if (const Line &line = polygon_sections[section_idx][line_idx]; line.a != line.b && line.length() < length_filter) {
                std::set<std::pair<int, int>> to_remove{{section_idx, line_idx}};
                std::vector<Node>             to_visit{{section_idx, line_idx}};

                bool initial_touches_long_lines = false;
                if (section_idx > 0) {
                    for (int prev_line_idx = 0; prev_line_idx < polygon_sections[section_idx - 1].size(); prev_line_idx++) {
                        if (const Line &nl = polygon_sections[section_idx - 1][prev_line_idx];
                            nl.a != nl.b && segments_overlap(line.a.y(), line.b.y(), nl.a.y(), nl.b.y())) {
                            initial_touches_long_lines = true;
                        }
                    }
                }

                while (!to_visit.empty()) {
                    Node        curr   = to_visit.back();
                    const Line &curr_l = polygon_sections[curr.section_idx][curr.line_idx];
                    if (curr.neighbours_explored) {
                        bool is_valid_for_removal = (curr_l.length() < length_filter) &&
                                                    ((int(to_remove.size()) - curr.skips_taken > min_removal_conut) ||
                                                     (curr.neighbours.empty() && !initial_touches_long_lines));
                        if (!is_valid_for_removal) {
                            for (const auto &n : curr.neighbours) {
                                if (to_remove.find(n) != to_remove.end()) {
                                    is_valid_for_removal = true;
                                    break;
                                }
                            }
                        }
                        if (!is_valid_for_removal) {
                            to_remove.erase({curr.section_idx, curr.line_idx});
                        }
                        to_visit.pop_back();
                    } else {
                        to_visit.back().neighbours_explored = true;
                        int  curr_index                     = to_visit.size() - 1;
                        bool can_use_skip                   = curr_l.length() <= length_filter && curr.skips_taken < skips_allowed;
                        if (curr.section_idx + 1 < polygon_sections.size()) {
                            for (int lidx = 0; lidx < polygon_sections[curr.section_idx + 1].size(); lidx++) {
                                if (const Line &nl = polygon_sections[curr.section_idx + 1][lidx];
                                    nl.a != nl.b && segments_overlap(curr_l.a.y(), curr_l.b.y(), nl.a.y(), nl.b.y()) &&
                                    (nl.length() < length_filter || can_use_skip)) {
                                    to_visit[curr_index].neighbours.push_back({curr.section_idx + 1, lidx});
                                    to_remove.insert({curr.section_idx + 1, lidx});
                                    Node next_node{curr.section_idx + 1, lidx, curr.skips_taken + (nl.length() >= length_filter)};
                                    to_visit.push_back(next_node);
                                }
                            }
                        }
                    }
                }

                for (const auto &pair : to_remove) {
                    Line &l = polygon_sections[pair.first][pair.second];
                    l.a     = l.b;
                }
            }
        }
    }

    for (size_t section_idx = 0; section_idx < polygon_sections.size(); section_idx++) {
        polygon_sections[section_idx].erase(std::remove_if(polygon_sections[section_idx].begin(), polygon_sections[section_idx].end(),
                                                           [](const Line &s) { return s.a == s.b; }),
                                            polygon_sections[section_idx].end());
        std::sort(polygon_sections[section_idx].begin(), polygon_sections[section_idx].end(),
                  [](const Line &a, const Line &b) { return a.a.y() < b.b.y(); });
    }

    ThickPolylines thick_polylines_out;
    {
        ThickPolylines current_traced_paths;
        for (const auto &polygon_slice : polygon_sections) {
            std::unordered_set<const Line *> used_segments;
            for (ThickPolyline &traced_path : current_traced_paths) {
                Point max_y = traced_path.last_point();
                Point min_y = traced_path.points[traced_path.size() - 2];

                if (max_y.y() < min_y.y())
                    std::swap(max_y, min_y);

                auto candidates_begin = std::upper_bound(polygon_slice.begin(), polygon_slice.end(), min_y,
                                                         [](const Point &low, const Line &seg) { return seg.b.y() > low.y(); });
                auto candidates_end   = std::upper_bound(polygon_slice.begin(), polygon_slice.end(), max_y,
                                                         [](const Point &high, const Line &seg) { return seg.a.y() > high.y(); });
                bool segment_added    = false;
                for (auto candidate = candidates_begin; candidate != candidates_end && !segment_added; candidate++) {
                    if (used_segments.find(&(*candidate)) != used_segments.end()) {
                        continue;
                    }
                    if ((traced_path.last_point() - candidate->a).cast<double>().squaredNorm() < squared_distance_limit_reconnection) {
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.points.push_back(candidate->a);
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.points.push_back(candidate->b);
                        traced_path.width.push_back(scaled_spacing);
                        used_segments.insert(&(*candidate));
                        segment_added = true;
                    } else if ((traced_path.last_point() - candidate->b).cast<double>().squaredNorm() <
                               squared_distance_limit_reconnection) {
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.points.push_back(candidate->b);
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.width.push_back(scaled_spacing);
                        traced_path.points.push_back(candidate->a);
                        traced_path.width.push_back(scaled_spacing);
                        used_segments.insert(&(*candidate));
                        segment_added = true;
                    }
                }

                if (!segment_added) {
                    // Zero overlapping segments. Finish the polyline.
                    thick_polylines_out.push_back(std::move(traced_path));
                    traced_path.clear();
                }
            }

            current_traced_paths.erase(std::remove_if(current_traced_paths.begin(), current_traced_paths.end(),
                                                      [](const ThickPolyline &tp) { return tp.empty(); }),
                                       current_traced_paths.end());

            for (const Line &segment : polygon_slice) {
                if (used_segments.find(&segment) == used_segments.end()) {
                    ThickPolyline &new_path = current_traced_paths.emplace_back();
                    new_path.points.push_back(segment.a);
                    new_path.width.push_back(scaled_spacing);
                    new_path.points.push_back(segment.b);
                    new_path.width.push_back(scaled_spacing);
                    new_path.endpoints = {true, true};
                }
            }
        }

        thick_polylines_out.insert(thick_polylines_out.end(), current_traced_paths.begin(), current_traced_paths.end());
    }

    Polygons reconstructed_area{};
    // reconstruct polygon from polygon sections
    {
        struct TracedPoly
        {
            std::vector<Point> lows;
            std::vector<Point> highs;
        };

        std::vector<TracedPoly> current_traced_polys;
        for (const auto &polygon_slice : polygon_sections) {
            std::unordered_set<const Line *> used_segments;
            for (TracedPoly &traced_poly : current_traced_polys) {
                auto candidates_begin = std::upper_bound(polygon_slice.begin(), polygon_slice.end(), traced_poly.lows.back(),
                                                    [](const Point &low, const Line &seg) { return seg.b.y() > low.y(); });
                auto candidates_end = std::upper_bound(polygon_slice.begin(), polygon_slice.end(), traced_poly.highs.back(),
                                                  [](const Point &high, const Line &seg) { return seg.a.y() > high.y(); });

                bool segment_added = false;
                for (auto candidate = candidates_begin; candidate != candidates_end && !segment_added; candidate++) {
                    if (used_segments.find(&(*candidate)) != used_segments.end()) {
                        continue;
                    }
                    if ((traced_poly.lows.back() - candidates_begin->a).cast<double>().squaredNorm() < squared_distance_limit_reconnection) {
                        traced_poly.lows.push_back(candidates_begin->a);
                    } else {
                        traced_poly.lows.push_back(traced_poly.lows.back() + Point{scaled_spacing / 2, 0});
                        traced_poly.lows.push_back(candidates_begin->a - Point{scaled_spacing / 2, 0});
                        traced_poly.lows.push_back(candidates_begin->a);
                    }

                    if ((traced_poly.highs.back() - candidates_begin->b).cast<double>().squaredNorm() <
                        squared_distance_limit_reconnection) {
                        traced_poly.highs.push_back(candidates_begin->b);
                    } else {
                        traced_poly.highs.push_back(traced_poly.highs.back() + Point{scaled_spacing / 2, 0});
                        traced_poly.highs.push_back(candidates_begin->b - Point{scaled_spacing / 2, 0});
                        traced_poly.highs.push_back(candidates_begin->b);
                    }
                    segment_added = true;
                    used_segments.insert(&(*candidates_begin));
                }

                if (!segment_added) {
                    // Zero or multiple overlapping segments. Resolving this is nontrivial,
                    // so we just close this polygon and maybe open several new. This will hopefully happen much less often
                    traced_poly.lows.push_back(traced_poly.lows.back() + Point{scaled_spacing / 2, 0});
                    traced_poly.highs.push_back(traced_poly.highs.back() + Point{scaled_spacing / 2, 0});
                    Polygon &new_poly = reconstructed_area.emplace_back(std::move(traced_poly.lows));
                    new_poly.points.insert(new_poly.points.end(), traced_poly.highs.rbegin(), traced_poly.highs.rend());
                    traced_poly.lows.clear();
                    traced_poly.highs.clear();
                }
            }

            current_traced_polys.erase(std::remove_if(current_traced_polys.begin(), current_traced_polys.end(),
                                                      [](const TracedPoly &tp) { return tp.lows.empty(); }),
                                       current_traced_polys.end());

            for (const auto &segment : polygon_slice) {
                if (used_segments.find(&segment) == used_segments.end()) {
                    TracedPoly &new_tp = current_traced_polys.emplace_back();
                    new_tp.lows.push_back(segment.a - Point{scaled_spacing / 2, 0});
                    new_tp.lows.push_back(segment.a);
                    new_tp.highs.push_back(segment.b - Point{scaled_spacing / 2, 0});
                    new_tp.highs.push_back(segment.b);
                }
            }
        }

        // add not closed polys
        for (TracedPoly &traced_poly : current_traced_polys) {
            Polygon &new_poly = reconstructed_area.emplace_back(std::move(traced_poly.lows));
            new_poly.points.insert(new_poly.points.end(), traced_poly.highs.rbegin(), traced_poly.highs.rend());
        }
    }

    reconstructed_area                     = closing(reconstructed_area, float(SCALED_EPSILON), float(SCALED_EPSILON));
    ExPolygons gaps_for_additional_filling = diff_ex(filled_area, reconstructed_area);
     if (this->overlap != 0) {
        gaps_for_additional_filling = offset_ex(gaps_for_additional_filling, scaled<float>(this->overlap));
    }
    gaps_for_additional_filling            = opening_ex(gaps_for_additional_filling, 0.3 * scaled_spacing);

    BoundingBox bbox = get_extents(filled_area);
    bbox.offset(scale_(1.));
    ::Slic3r::SVG svg(debug_out_path(("surface" + std::to_string(surface->area())).c_str()).c_str(), bbox);
    svg.draw(to_lines(filled_area), "red", scale_(0.4));
    svg.draw(to_lines(reconstructed_area), "blue", scale_(0.3));
    svg.draw(to_lines(gaps_for_additional_filling), "green", scale_(0.2));
    svg.draw(vertical_lines, "black", scale_(0.1));
    svg.Close();

    for (ExPolygon &ex_poly : gaps_for_additional_filling) {
        Point    bbox_size   = ex_poly.contour.bounding_box().size();
        coord_t  loops_count = std::max(bbox_size.x(), bbox_size.y()) / scaled_spacing + 1;
        Polygons polygons    = to_polygons(ex_poly);
        Arachne::WallToolPaths wall_tool_paths(polygons, scaled_spacing, scaled_spacing, loops_count, 0, params.layer_height, *this->print_object_config, *this->print_config);
        if (std::vector<Arachne::VariableWidthLines> loops = wall_tool_paths.getToolPaths(); !loops.empty()) {
            std::vector<const Arachne::ExtrusionLine *> all_extrusions;
            for (Arachne::VariableWidthLines &loop : loops) {
                if (loop.empty())
                    continue;
                for (const Arachne::ExtrusionLine &wall : loop)
                    all_extrusions.emplace_back(&wall);
            }

            // Split paths using a nearest neighbor search.
            size_t firts_poly_idx = thick_polylines_out.size();
            Point  last_pos(0, 0);
            for (const Arachne::ExtrusionLine *extrusion : all_extrusions) {
                if (extrusion->empty())
                    continue;

                ThickPolyline thick_polyline = Arachne::to_thick_polyline(*extrusion);
                if (thick_polyline.length() == 0.)
                    //FIXME this should not happen.
                    continue;
                assert(thick_polyline.size() > 1);
                assert(thick_polyline.length() > 0.);
                //assert(thick_polyline.points.size() == thick_polyline.width.size());
                if (extrusion->is_closed)
                    thick_polyline.start_at_index(nearest_point_index(thick_polyline.points, last_pos));

                assert(thick_polyline.size() > 1);
                //assert(thick_polyline.points.size() == thick_polyline.width.size());
                thick_polylines_out.emplace_back(std::move(thick_polyline));
                last_pos = thick_polylines_out.back().last_point();
            }

            // clip the paths to prevent the extruder from getting exactly on the first point of the loop
            // Keep valid paths only.
            size_t j = firts_poly_idx;
            for (size_t i = firts_poly_idx; i < thick_polylines_out.size(); ++i) {
                assert(thick_polylines_out[i].size() > 1);
                assert(thick_polylines_out[i].length() > 0.);
                //assert(thick_polylines_out[i].points.size() == thick_polylines_out[i].width.size());
                thick_polylines_out[i].clip_end(this->loop_clipping);
                assert(thick_polylines_out[i].size() > 1);
                if (thick_polylines_out[i].is_valid()) {
                    if (j < i)
                        thick_polylines_out[j] = std::move(thick_polylines_out[i]);
                    ++j;
                }
            }
            if (j < thick_polylines_out.size())
                thick_polylines_out.erase(thick_polylines_out.begin() + int(j), thick_polylines_out.end());
        }
    }

    rotate_thick_polylines(thick_polylines_out, cos(-aligning_angle), sin(-aligning_angle));

    return thick_polylines_out;
}

} // namespace Slic3r
