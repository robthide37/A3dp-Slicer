#include "LayerChanges.hpp"
#include "libslic3r/ClipperUtils.hpp"

namespace Slic3r::GCode::Impl::LayerChanges {

Polygon generate_regular_polygon(
    const Point &centroid, const Point &start_point, const unsigned points_count
) {
    Points points;
    points.reserve(points_count);
    const double part_angle{2 * M_PI / points_count};
    for (unsigned i = 0; i < points_count; ++i) {
        const double current_angle{i * part_angle};
        points.emplace_back(scaled(std::cos(current_angle)), scaled(std::sin(current_angle)));
    }

    Polygon regular_polygon{points};
    const Vec2d current_vector{unscaled(regular_polygon.points.front())};
    const Vec2d expected_vector{unscaled(start_point) - unscaled(centroid)};

    const double current_scale = current_vector.norm();
    const double expected_scale = expected_vector.norm();
    regular_polygon.scale(expected_scale / current_scale);

    regular_polygon.rotate(angle(current_vector, expected_vector));

    regular_polygon.translate(centroid);

    return regular_polygon;
}

Bed::Bed(const std::vector<Vec2d> &shape, const double padding)
    : inner_offset(get_inner_offset(shape, padding)), centroid(unscaled(inner_offset.centroid())) {}

bool Bed::contains_within_padding(const Vec2d &point) const {
    return inner_offset.contains(scaled(point));
}

Polygon Bed::get_inner_offset(const std::vector<Vec2d> &shape, const double padding) {
    Points shape_scaled;
    shape_scaled.reserve(shape.size());
    using std::begin, std::end, std::back_inserter, std::transform;
    transform(begin(shape), end(shape), back_inserter(shape_scaled), [](const Vec2d &point) {
        return scaled(point);
    });
    const Polygons inner_offset{shrink({Polygon{shape_scaled}}, scaled(padding))};
    if (inner_offset.empty()) {
        return Polygon{};
    }
    return inner_offset.front();
}

} // namespace Slic3r::GCode::Impl::LayerChanges
