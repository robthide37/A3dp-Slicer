#include <catch2/catch.hpp>

#include <libslic3r/ClosestPoint.hpp>
#include <libslic3r/Point.hpp>

using namespace Slic3r;
TEST_CASE("Find the closest point from 2", "[ClosestPoint]")
{
    Points pts = {{0, 1}, {0, 2}};
    CHECK(std::is_sorted(pts.begin(), pts.end(), 
        closestPoint::sort_fnc<Point>));
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
    CHECK(std::is_sorted(pts.begin(), pts.end(), 
        closestPoint::sort_fnc<Point>));

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

TEST_CASE("Find the closest point from 9 unsorted", "[ClosestPoint]")
{
    //  4 - 3 - 0
    //  |   |   |
    //  1 - 6 - 5
    //  |   |   |
    //  2 - 7 - 8
    Points pts = {
        /*0*/ {3, 3}, 
        /*1*/ {-3, 0}, 
        /*2*/ {-3, -3}, 
        /*3*/ {0, 3},
        /*4*/ {-3, 3},
        /*5*/ {3, 0},
        /*6*/ {0, 0},
        /*7*/ {0, -3},
        /*8*/ {3, -3}
    };
    CHECK(4 == find_closest(Point{-4, 4}, pts));
    CHECK(4 == find_closest(Point{-2, 2}, pts));
    // check center
    CHECK(6 == find_closest(Point{-1, 1}, pts));
    CHECK(6 == find_closest(Point{ 0, 1}, pts));
    CHECK(6 == find_closest(Point{ 1, 1}, pts));
    CHECK(6 == find_closest(Point{-1, 0}, pts));
    CHECK(6 == find_closest(Point{ 0, 0}, pts));
    CHECK(6 == find_closest(Point{ 1, 0}, pts));
    CHECK(6 == find_closest(Point{-1,-1}, pts));
    CHECK(6 == find_closest(Point{ 0,-1}, pts));
    CHECK(6 == find_closest(Point{ 1,-1}, pts));
    CHECK(6 == find_closest(Point{ 0, 0}, pts));

    CHECK(5 == find_closest(Point{2, 1}, pts));
    CHECK(8 == find_closest(Point{2,-2}, pts));
}