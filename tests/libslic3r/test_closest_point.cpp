#include <catch2/catch.hpp>

#include <libslic3r/ClosestPoint.hpp>

using namespace Slic3r;
TEST_CASE("Find the closest point from 2", "[ClosestPoint]")
{
    Points pts = {{0, 1}, {0, 2}};
    CHECK(std::is_sorted(pts.begin(), pts.end(), closestPoint::sort_fnc));
    CHECK(0 == find_closest_in_sorted(Point{0, 0}, pts));
    CHECK(0 == find_closest_in_sorted(Point{1, 1}, pts));
    CHECK(1 == find_closest_in_sorted(Point{1, 2}, pts));
}

TEST_CASE("Find the closest point from 9", "[ClosestPoint]")
{
    //  0 - 3 - 6
    //  |   |   |
    //  1 - 4 - 7
    //  |   |   |
    //  2 - 5 - 8
    Points pts = {{-3, 3}, {-3, 0}, {-3, -3}, {0, 3}, {0, 0},
                  {0, -3}, {3, 3}, {3, 0}, {3, -3}};
    CHECK(std::is_sorted(pts.begin(), pts.end(), closestPoint::sort_fnc));

    CHECK(0 == find_closest_in_sorted(Point{-4, 4}, pts));
    CHECK(0 == find_closest_in_sorted(Point{-2, 2}, pts));
    // check center
    CHECK(4 == find_closest_in_sorted(Point{-1, 1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 0, 1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 1, 1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{-1, 0}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 0, 0}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 1, 0}, pts));
    CHECK(4 == find_closest_in_sorted(Point{-1,-1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 0,-1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 1,-1}, pts));
    CHECK(4 == find_closest_in_sorted(Point{ 0, 0}, pts));

    CHECK(7 == find_closest_in_sorted(Point{2, 1}, pts));
    CHECK(8 == find_closest_in_sorted(Point{2,-2}, pts));
}
