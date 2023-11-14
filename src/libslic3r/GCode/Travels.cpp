#include "Travels.hpp"

namespace Slic3r::GCode::Impl::Travels {

Points3 generate_flat_travel(tcb::span<const Point> xy_path, const float elevation) {
    Points3 result;
    result.reserve(xy_path.size() - 1);
    for (const Point &point : xy_path.subspan(1)) {
        result.emplace_back(point.x(), point.y(), scaled(elevation));
    }
    return result;
}

Vec2d place_at_segment(
    const Vec2d &current_point, const Vec2d &previous_point, const double distance
) {
    Vec2d direction = (current_point - previous_point).normalized();
    return previous_point + direction * distance;
}

std::vector<DistancedPoint> slice_xy_path(
    tcb::span<const Point> xy_path, tcb::span<const double> sorted_distances
) {
    assert(xy_path.size() >= 2);
    std::vector<DistancedPoint> result;
    result.reserve(xy_path.size() + sorted_distances.size());
    double total_distance{0};
    result.emplace_back(DistancedPoint{xy_path.front(), 0});
    Point previous_point = result.front().point;
    std::size_t offset{0};
    for (const Point &point : xy_path.subspan(1)) {
        Vec2d unscaled_point{unscaled(point)};
        Vec2d unscaled_previous_point{unscaled(previous_point)};
        const double current_segment_length = (unscaled_point - unscaled_previous_point).norm();
        for (const double distance_to_add : sorted_distances.subspan(offset)) {
            if (distance_to_add <= total_distance + current_segment_length) {
                Point to_place = scaled(place_at_segment(
                    unscaled_point, unscaled_previous_point, distance_to_add - total_distance
                ));
                if (to_place != previous_point && to_place != point) {
                    result.emplace_back(DistancedPoint{to_place, distance_to_add});
                }
                ++offset;
            } else {
                break;
            }
        }
        total_distance += current_segment_length;
        result.emplace_back(DistancedPoint{point, total_distance});
        previous_point = point;
    }
    return result;
}

struct ElevatedTravelParams
{
    double lift_height{};
    double slope_end{};
};

struct ElevatedTravelFormula
{
    double operator()(double distance_from_start) const {
        if (distance_from_start < this->params.slope_end) {
            const double lift_percent = distance_from_start / this->params.slope_end;
            return lift_percent * this->params.lift_height;
        } else {
            return this->params.lift_height;
        }
    }

    ElevatedTravelParams params{};
};

Points3 generate_elevated_travel(
    const tcb::span<const Point> xy_path,
    const std::vector<double> &ensure_points_at_distances,
    const double initial_elevation,
    const std::function<double(double)> &elevation
) {
    Points3 result{};

    std::vector<DistancedPoint> extended_xy_path = slice_xy_path(xy_path, ensure_points_at_distances);
    result.reserve(extended_xy_path.size());

    for (const DistancedPoint &point : extended_xy_path) {
        result.emplace_back(
            point.point.x(), point.point.y(),
            scaled(initial_elevation + elevation(point.distance_from_start))
        );
    }

    return result;
}

std::optional<double> get_first_crossed_line_distance(
    tcb::span<const Line> xy_path, const AABBTreeLines::LinesDistancer<Linef> &distancer
) {
    assert(!xy_path.empty());
    if (xy_path.empty()) {
        return {};
    }

    double traversed_distance = 0;
    for (const Line &line : xy_path) {
        const Linef unscaled_line = {unscaled(line.a), unscaled(line.b)};
        auto intersections = distancer.intersections_with_line<true>(unscaled_line);
        if (!intersections.empty()) {
            const Vec2d intersection = intersections.front().first;
            const double distance = traversed_distance + (unscaled_line.a - intersection).norm();
            if (distance > EPSILON) {
                return distance;
            } else if (intersections.size() >= 2) { // Edge case
                const Vec2d second_intersection = intersections[1].first;
                return traversed_distance + (unscaled_line.a - second_intersection).norm();
            }
        }
        traversed_distance += (unscaled_line.a - unscaled_line.b).norm();
    }

    return {};
}

std::optional<double> get_obstacle_adjusted_slope_end(
    const Lines &xy_path,
    const std::optional<AABBTreeLines::LinesDistancer<Linef>> &previous_layer_distancer
) {
    if (!previous_layer_distancer) {
        return std::nullopt;
    }
    std::optional<double> first_obstacle_distance =
        get_first_crossed_line_distance(xy_path, *previous_layer_distancer);
    if (!first_obstacle_distance) {
        return std::nullopt;
    }
    return *first_obstacle_distance;
}

ElevatedTravelParams get_elevated_traval_params(
    const Lines &xy_path,
    const FullPrintConfig &config,
    const unsigned extruder_id,
    const std::optional<AABBTreeLines::LinesDistancer<Linef>> &previous_layer_distancer
) {
    ElevatedTravelParams elevation_params{};
    if (!config.travel_ramping_lift.get_at(extruder_id)) {
        elevation_params.slope_end = 0;
        elevation_params.lift_height = config.retract_lift.get_at(extruder_id);
        return elevation_params;
    }
    elevation_params.lift_height = config.travel_max_lift.get_at(extruder_id);

    const double slope_deg = config.travel_slope.get_at(extruder_id);

    if (slope_deg >= 90 || slope_deg <= 0) {
        elevation_params.slope_end = 0;
    } else {
        const double slope_rad = slope_deg * (M_PI / 180); // rad
        elevation_params.slope_end = elevation_params.lift_height / std::tan(slope_rad);
    }

    std::optional<double> obstacle_adjusted_slope_end{
        get_obstacle_adjusted_slope_end(xy_path, previous_layer_distancer)};

    if (obstacle_adjusted_slope_end && obstacle_adjusted_slope_end < elevation_params.slope_end) {
        elevation_params.slope_end = *obstacle_adjusted_slope_end;
    }

    return elevation_params;
}

Points3 generate_travel_to_extrusion(
    const Polyline &xy_path,
    const FullPrintConfig &config,
    const unsigned extruder_id,
    const double initial_elevation,
    const std::optional<AABBTreeLines::LinesDistancer<Linef>> &previous_layer_distancer,
    const Point &xy_path_coord_origin
) {
    const double upper_limit = config.retract_lift_below.get_at(extruder_id);
    const double lower_limit = config.retract_lift_above.get_at(extruder_id);
    if ((lower_limit > 0 && initial_elevation < lower_limit) ||
        (upper_limit > 0 && initial_elevation > upper_limit)) {
        return generate_flat_travel(xy_path.points, initial_elevation);
    }

    Lines global_xy_path;
    for (const Line &line : xy_path.lines()) {
        global_xy_path.emplace_back(line.a + xy_path_coord_origin, line.b + xy_path_coord_origin);
    }

    ElevatedTravelParams elevation_params{
        get_elevated_traval_params(global_xy_path, config, extruder_id, previous_layer_distancer)};

    const std::vector<double> ensure_points_at_distances{elevation_params.slope_end};

    Points3 result{generate_elevated_travel(
        xy_path.points, ensure_points_at_distances, initial_elevation,
        ElevatedTravelFormula{elevation_params}
    )};

    result.emplace_back(xy_path.back().x(), xy_path.back().y(), scaled(initial_elevation));
    return result;
}
} // namespace Slic3r::GCode::Impl::Travels
