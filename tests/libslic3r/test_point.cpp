
/**
 * Ported from xs/t/03_point.t
 *  - it used to check ccw() but it does not exist anymore
 *  and cross product uses doubles
 */

#include <catch2/catch.hpp>
#include <libslic3r/Point.hpp>
#include <libslic3r/MultiPoint.hpp>
#include "test_utils.hpp"

using namespace Slic3r;

TEST_CASE("Nearest point", "[Point]") {
    const Point point{10, 15};
    const Point point2{30, 15};

    const Point nearest{nearest_point({point2, Point{100, 200}}, point).first};
    CHECK(nearest == point2);
}

TEST_CASE("Distance to line", "[Point]") {
    const Line line{{0, 0}, {100, 0}};
    CHECK(line.distance_to(Point{0, 0}) == Approx(0));
    CHECK(line.distance_to(Point{100, 0}) == Approx(0));
    CHECK(line.distance_to(Point{50, 0}) == Approx(0));
    CHECK(line.distance_to(Point{150, 0}) == Approx(50));
    CHECK(line.distance_to(Point{0, 50}) == Approx(50));
    CHECK(line.distance_to(Point{50, 50}) == Approx(50));
    CHECK(line.perp_distance_to(Point{50, 50}) == Approx(50));
    CHECK(line.perp_distance_to(Point{150, 50}) == Approx(50));
}

TEST_CASE("Distance to diagonal line", "[Point]") {
    const Line line{{50, 50}, {125, -25}};
    CHECK(std::abs(line.distance_to(Point{100, 0})) == Approx(0));
}

TEST_CASE("Perp distance to line does not overflow", "[Point]") {
    const Line line{
        {18335846, 18335845},
        {18335846, 1664160},
    };

    CHECK(line.distance_to(Point{1664161, 18335848}) == Approx(16671685));
}

TEST_CASE("douglas_peucker", "[MultiPoint]") {
    Points pt_in{Point{-15395527,-2111143},Point{-13895526,-2111143},Point{-13895526,2888857},Point{-13895526,7888857},Point{-15395527,7888857},Point{-16895527,7888857},Point{-16895527,2888857},Point{-16895527,-2111143},Point{-15395527,-2111143}};
    Points pt_check{Point{-15395527,-2111143},Point{-13895526,-2111143},Point{-13895526,7888857},Point{-16895527,7888857},Point{-16895527,-2111143},Point{-15395527,-2111143}};
    Points pt_out;
    douglas_peucker_impl(pt_in.begin(), pt_in.end(), std::back_inserter(pt_out), 3125, [](const Point &p) { return p; });
    CHECK(pt_check == pt_out);
}

