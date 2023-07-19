#include "ExtrusionProcessor.hpp"

namespace Slic3r {

ExtrusionPaths calculate_and_split_overhanging_extrusions(const ExtrusionPath                             &path,
                                                          const AABBTreeLines::LinesDistancer<Linef>      &unscaled_prev_layer,
                                                          const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines)
{
    std::vector<ExtendedPoint>           extended_points = estimate_points_properties<true, true, true, true>(path.polyline.points,
                                                                                                    unscaled_prev_layer, path.width());
    std::vector<std::pair<float, float>> calculated_distances(extended_points.size());

    for (size_t i = 0; i < extended_points.size(); i++) {
        const ExtendedPoint &curr = extended_points[i];
        const ExtendedPoint &next = extended_points[i + 1 < extended_points.size() ? i + 1 : i];

        // The following code artifically increases the distance to provide slowdown for extrusions that are over curled lines
        float        proximity_to_curled_lines = 0.0;
        const double dist_limit                = 10.0 * path.width();
        {
            Vec2d middle       = 0.5 * (curr.position + next.position);
            auto  line_indices = prev_layer_curled_lines.all_lines_in_radius(Point::new_scale(middle), scale_(dist_limit));
            if (!line_indices.empty()) {
                double len = (next.position - curr.position).norm();
                // For long lines, there is a problem with the additional slowdown. If by accident, there is small curled line near the middle
                // of this long line
                //  The whole segment gets slower unnecesarily. For these long lines, we do additional check whether it is worth slowing down.
                // NOTE that this is still quite rough approximation, e.g. we are still checking lines only near the middle point
                // TODO maybe split the lines into smaller segments before running this alg? but can be demanding, and GCode will be huge
                if (len > 8) {
                    Vec2d dir   = Vec2d(next.position - curr.position) / len;
                    Vec2d right = Vec2d(-dir.y(), dir.x());

                    Polygon box_of_influence = {
                        scaled(Vec2d(curr.position + right * dist_limit)),
                        scaled(Vec2d(next.position + right * dist_limit)),
                        scaled(Vec2d(next.position - right * dist_limit)),
                        scaled(Vec2d(curr.position - right * dist_limit)),
                    };

                    double projected_lengths_sum = 0;
                    for (size_t idx : line_indices) {
                        const CurledLine &line   = prev_layer_curled_lines.get_line(idx);
                        Lines             inside = intersection_ln({{line.a, line.b}}, {box_of_influence});
                        if (inside.empty())
                            continue;
                        double projected_length = abs(dir.dot(unscaled(Vec2d((inside.back().b - inside.back().a).cast<double>()))));
                        projected_lengths_sum += projected_length;
                    }
                    if (projected_lengths_sum < 0.4 * len) {
                        line_indices.clear();
                    }
                }

                for (size_t idx : line_indices) {
                    const CurledLine &line                 = prev_layer_curled_lines.get_line(idx);
                    float             distance_from_curled = unscaled(line_alg::distance_to(line, Point::new_scale(middle)));
                    float proximity = (1.0 - (distance_from_curled / dist_limit)) * (1.0 - (distance_from_curled / dist_limit)) *
                                      (line.curled_height / (path.height() * 10.0f)); // max_curled_height_factor from SupportSpotGenerator
                    proximity_to_curled_lines = std::max(proximity_to_curled_lines, proximity);
                }
            }
        }
        calculated_distances[i].first  = std::max(curr.distance, next.distance);
        calculated_distances[i].second = proximity_to_curled_lines;
    }

    ExtrusionPaths      result;
    ExtrusionAttributes new_attrs = path.attributes();
    new_attrs.overhang_attributes = std::optional<OverhangAttributes>({calculated_distances[0].first, calculated_distances[0].second});
    result.emplace_back(new_attrs);
    result.back().polyline.append(Point::new_scale(extended_points[0].position));
    for (size_t i = 1; i < extended_points.size(); i++) {
        if (std::abs(calculated_distances[i].first - calculated_distances[i - 1].first) < path.width() * 0.1 &&
            std::abs(calculated_distances[i].second - calculated_distances[i - 1].second) < 0.1) {
            // do not start new path, the attributes are similar enough
        } else {
            result.back().polyline.append(Point::new_scale(extended_points[i].position));
            new_attrs.overhang_attributes->max_distance_from_prev_layer = calculated_distances[i].first;
            new_attrs.overhang_attributes->proximity_to_curled_lines    = calculated_distances[i].second;
            result.emplace_back(new_attrs);
        }
        result.back().polyline.append(Point::new_scale(extended_points[i].position));
    }

    return result;
};

ExtrusionEntityCollection calculate_and_split_overhanging_extrusions(const ExtrusionEntityCollection            *ecc,
                                                                     const AABBTreeLines::LinesDistancer<Linef> &unscaled_prev_layer,
                                                                     const AABBTreeLines::LinesDistancer<CurledLine> &prev_layer_curled_lines)
{
    ExtrusionEntityCollection result{};
    result.no_sort = ecc->no_sort;
    for (const auto *e : ecc->entities) {
        if (e == nullptr) {
             BOOST_LOG_TRIVIAL(debug) << "perimeters collections contain nullptr entities for some reason";
        } else if (auto *col = static_cast<const ExtrusionEntityCollection *>(e)) {
            result.append(calculate_and_split_overhanging_extrusions(col, unscaled_prev_layer, prev_layer_curled_lines));
        } else if (auto *loop = static_cast<const ExtrusionLoop *>(e)) {
            ExtrusionLoop new_loop = *loop;
            new_loop.paths.clear();
            for (const ExtrusionPath &p : loop->paths) {
                auto paths = calculate_and_split_overhanging_extrusions(p, unscaled_prev_layer, prev_layer_curled_lines);
                new_loop.paths.insert(new_loop.paths.end(), paths.begin(), paths.end());
            }
            result.append(new_loop);
        } else if (auto *mp = static_cast<const ExtrusionMultiPath *>(e)) {
            ExtrusionMultiPath new_mp = *mp;
            new_mp.paths.clear();
            for (const ExtrusionPath &p : mp->paths) {
                auto paths = calculate_and_split_overhanging_extrusions(p, unscaled_prev_layer, prev_layer_curled_lines);
                new_mp.paths.insert(new_mp.paths.end(), paths.begin(), paths.end());
            }
            result.append(new_mp);
        } else if (auto *op = static_cast<const ExtrusionPathOriented *>(e)) {
            auto paths = calculate_and_split_overhanging_extrusions(*op, unscaled_prev_layer, prev_layer_curled_lines);
            for (const ExtrusionPath &p : paths) {
                result.append(ExtrusionPathOriented(p.polyline, p.attributes()));
            }
        } else if (auto *p = static_cast<const ExtrusionPath *>(e)) {
            auto paths = calculate_and_split_overhanging_extrusions(*p, unscaled_prev_layer, prev_layer_curled_lines);
            result.append(paths);
        } else {
            throw Slic3r::InvalidArgument("Unknown extrusion entity type");
        }
    }
    return result;
};

} // namespace Slic3r