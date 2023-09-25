#include "libslic3r/Point.hpp"
#include <catch2/catch.hpp>
#include <libslic3r/SupportSpotsGenerator.hpp>

using namespace Slic3r;
using namespace SupportSpotsGenerator;


TEST_CASE("Numerical integral calculation compared with exact solution.", "[SupportSpotsGenerator]") {
    const float width = 10;
    const float height = 20;
    const Polygon polygon = {
        scaled(Vec2f{-width / 2, -height / 2}),
        scaled(Vec2f{width / 2, -height / 2}),
        scaled(Vec2f{width / 2, height / 2}),
        scaled(Vec2f{-width / 2, height / 2})
    };

    const Integrals integrals{{polygon}};
    CHECK(integrals.area == Approx(width * height));
    CHECK(integrals.x_i.x() == Approx(0));
    CHECK(integrals.x_i.y() == Approx(0));
    CHECK(integrals.x_i_squared.x() == Approx(std::pow(width, 3) * height / 12));
    CHECK(integrals.x_i_squared.y() == Approx(width * std::pow(height, 3) / 12));
}

TEST_CASE("Moment values and ratio check.", "[SupportSpotsGenerator]") {
    const float width = 40;
    const float height = 2;

    // Moments are calculated at centroid.
    // Polygon centroid must not be (0, 0).
    const Polygon polygon = {
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{width, 0}),
        scaled(Vec2f{width, height}),
        scaled(Vec2f{0, height})
    };

    const Integrals integrals{{polygon}};

    const Vec2f x_axis{1, 0};
    const float x_axis_moment = compute_second_moment(integrals, x_axis);

    const Vec2f y_axis{0, 1};
    const float y_axis_moment = compute_second_moment(integrals, y_axis);

    const float moment_ratio = std::pow(width / height, 2);

    // Ensure the object transaltion has no effect.
    CHECK(x_axis_moment == Approx(width * std::pow(height, 3) / 12));
    CHECK(y_axis_moment == Approx(std::pow(width, 3) * height / 12));
    // If the object is "wide" the y axis moments should be large compared to x axis moment.
    CHECK(y_axis_moment / x_axis_moment == Approx(moment_ratio));
}

TEST_CASE("Moments calculation for rotated axis.", "[SupportSpotsGenerator]") {

    Polygon polygon = {
        scaled(Vec2f{6.362284076172198, 138.9674202217155}),
        scaled(Vec2f{97.48779843751677, 106.08136606617076}),
        scaled(Vec2f{135.75221821532384, 66.84428834668765}),
        scaled(Vec2f{191.5308049852741, 45.77905628725614}),
        scaled(Vec2f{182.7525148049201, 74.01799041087513}),
        scaled(Vec2f{296.83210979283473, 196.80022572637228}),
        scaled(Vec2f{215.16434429179148, 187.45715418834143}),
        scaled(Vec2f{64.64574271229334, 284.293883209721}),
        scaled(Vec2f{110.76507036894843, 174.35633141113783}),
        scaled(Vec2f{77.56229640885199, 189.33057746591336})
    };

    Integrals integrals{{polygon}};

    std::mt19937 generator{std::random_device{}()};
    std::uniform_real_distribution<float> angle_distribution{0.f, float(2*M_PI)};

    // Meassured counterclockwise from (1, 0)
    const float angle = angle_distribution(generator);
    Vec2f axis{std::cos(angle), std::sin(angle)};

    float moment_calculated_then_rotated = compute_second_moment(
        integrals,
        axis
    );

    // We want to rotate the object clockwise by angle to align the axis with (1, 0)
    // Method .rotate is counterclockwise for positive angle
    polygon.rotate(-angle);

    Integrals integrals_rotated{{polygon}};
    float moment_rotated_polygon = compute_second_moment(
        integrals_rotated,
        Vec2f{1, 0}
    );

    CHECK(moment_calculated_then_rotated == Approx(moment_rotated_polygon));
}

TEST_CASE("TODO", "[SupportSpotsGenerator]") {
    const Polyline polyline{
        Point{scaled(Vec2f{0, 0})},
        Point{scaled(Vec2f{1, 0})},
    };
    ExtrusionAttributes attributes;
    attributes.width = 0.1;
    const ExtrusionPath path{polyline, attributes};
    ExtrusionEntityCollection collection;
    collection.append(path);
    std::vector<const ExtrusionEntityCollection*> collections{&collection};

    Polygons polygons = path.polygons_covered_by_width();
    for (const Polygon& polygon : polygons) {
        std::cout << "Polygon: " << std::endl;
        for (const Line& line : polygon.lines()) {
            std::cout << "(" << line.a.x() << ", " << line.a.y() << ")" << std::endl;
        }
    }
}
