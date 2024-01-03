#include "Travels.hpp"

namespace Slic3r::GCode::Impl::Travels {

ElevatedTravelFormula::ElevatedTravelFormula(const ElevatedTravelParams &params)
    : smoothing_from(params.slope_end - params.blend_width / 2.0)
    , smoothing_to(params.slope_end + params.blend_width / 2.0)
    , blend_width(params.blend_width)
    , lift_height(params.lift_height)
    , slope_end(params.slope_end) {
    if (smoothing_from < 0) {
        smoothing_from = params.slope_end;
        smoothing_to = params.slope_end;
    }
}

double parabola(const double x, const double a, const double b, const double c) {
    return a * x * x + b * x + c;
}

double ElevatedTravelFormula::slope_function(double distance_from_start) const {
    if (distance_from_start < this->slope_end) {
        const double lift_percent = distance_from_start / this->slope_end;
        return lift_percent * this->lift_height;
    } else {
        return this->lift_height;
    }
}

double ElevatedTravelFormula::operator()(const double distance_from_start) const {
    if (distance_from_start > this->smoothing_from && distance_from_start < this->smoothing_to) {
        const double slope = this->lift_height / this->slope_end;

        // This is a part of a parabola going over a specific
        // range and with specific end slopes.
        const double a = -slope / 2.0 / this->blend_width;
        const double b = slope * this->smoothing_to / this->blend_width;
        const double c = this->lift_height + a * boost::math::pow<2>(this->smoothing_to);
        return parabola(distance_from_start, a, b, c);
    }
    return slope_function(distance_from_start);
}

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

struct SmoothingParams
{
    double blend_width{};
    unsigned points_count{1};
};

SmoothingParams get_smoothing_params(
    const double lift_height,
    const double slope_end,
    unsigned extruder_id,
    const double travel_length,
    const FullPrintConfig &config
) {
    if (config.gcode_flavor != gcfMarlinFirmware)
        // Smoothing is supported only on Marlin.
        return {0, 1};

    const double slope = lift_height / slope_end;
    const double max_machine_z_velocity = config.machine_max_feedrate_z.get_at(extruder_id);
    const double max_xy_velocity =
        Vec2d{
            config.machine_max_feedrate_x.get_at(extruder_id),
            config.machine_max_feedrate_y.get_at(extruder_id)}
            .norm();

    const double xy_acceleration = config.machine_max_acceleration_travel.get_at(extruder_id);

    const double xy_acceleration_time = max_xy_velocity / xy_acceleration;
    const double xy_acceleration_distance = 1.0 / 2.0 * xy_acceleration *
        boost::math::pow<2>(xy_acceleration_time);

    if (travel_length < xy_acceleration_distance) {
        return {0, 1};
    }

    const double max_z_velocity = std::min(max_xy_velocity * slope, max_machine_z_velocity);
    const double deceleration_time = max_z_velocity /
        config.machine_max_acceleration_z.get_at(extruder_id);
    const double deceleration_xy_distance = deceleration_time * max_xy_velocity;

    const double blend_width = slope_end > deceleration_xy_distance / 2.0 ? deceleration_xy_distance :
                                                                          slope_end * 2.0;

    const unsigned points_count = blend_width > 0 ?
        std::ceil(max_z_velocity / config.machine_max_jerk_z.get_at(extruder_id)) :
        1;

    if (blend_width <= 0     // When there is no blend with, there is no need for smoothing.
        || points_count > 6  // That would be way to many points. Do not do it at all.
        || points_count <= 0 // Always return at least one point.
    )
        return {0, 1};

    return {blend_width, points_count};
}

ElevatedTravelParams get_elevated_traval_params(
    const Polyline &xy_path,
    const FullPrintConfig &config,
    const unsigned extruder_id,
    const std::optional<AABBTreeLines::LinesDistancer<Linef>> &previous_layer_distancer
) {
    ElevatedTravelParams elevation_params{};
    if (!config.travel_ramping_lift.get_at(extruder_id)) {
        elevation_params.slope_end = 0;
        elevation_params.lift_height = config.retract_lift.get_at(extruder_id);
        elevation_params.blend_width = 0;
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
        get_obstacle_adjusted_slope_end(xy_path.lines(), previous_layer_distancer)};

    if (obstacle_adjusted_slope_end && obstacle_adjusted_slope_end < elevation_params.slope_end) {
        elevation_params.slope_end = *obstacle_adjusted_slope_end;
    }

    SmoothingParams smoothing_params{get_smoothing_params(
        elevation_params.lift_height, elevation_params.slope_end, extruder_id,
        unscaled(xy_path.length()), config
    )};

    elevation_params.blend_width = smoothing_params.blend_width;
    elevation_params.parabola_points_count = smoothing_params.points_count;
    return elevation_params;
}

/**
 * @brief Generate regulary spaced points on 1 axis. Includes both from and to.
 *
 * If count is 1, the point is in the middle of the range.
 */
std::vector<double> linspace(const double from, const double to, const unsigned count) {
    if (count == 0) {
        return {};
    }
    std::vector<double> result;
    result.reserve(count);
    if (count == 1) {
        result.emplace_back((from + to) / 2.0);
        return result;
    }
    const double step = (to - from) / count;
    for (unsigned i = 0; i < count - 1; ++i) {
        result.emplace_back(from + i * step);
    }
    result.emplace_back(to); // Make sure the last value is exactly equal to the value of "to".
    return result;
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

    Points global_xy_path;
    for (const Point &point : xy_path.points) {
        global_xy_path.emplace_back(point + xy_path_coord_origin);
    }

    ElevatedTravelParams elevation_params{get_elevated_traval_params(
        Polyline{std::move(global_xy_path)}, config, extruder_id, previous_layer_distancer
    )};

    const std::vector<double> ensure_points_at_distances = linspace(
        elevation_params.slope_end - elevation_params.blend_width / 2.0,
        elevation_params.slope_end + elevation_params.blend_width / 2.0,
        elevation_params.parabola_points_count
    );

    Points3 result{generate_elevated_travel(
        xy_path.points, ensure_points_at_distances, initial_elevation,
        ElevatedTravelFormula{elevation_params}
    )};

    result.emplace_back(xy_path.back().x(), xy_path.back().y(), scaled(initial_elevation));
    return result;
}
} // namespace Slic3r::GCode::Impl::Travels
