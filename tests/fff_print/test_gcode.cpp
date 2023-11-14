#include <catch2/catch.hpp>

#include <memory>

#include "libslic3r/GCode.hpp"

using namespace Slic3r;
using namespace Slic3r::GCode::Impl;

SCENARIO("Origin manipulation", "[GCode]") {
	Slic3r::GCodeGenerator gcodegen;
	WHEN("set_origin to (10,0)") {
    	gcodegen.set_origin(Vec2d(10,0));
    	REQUIRE(gcodegen.origin() == Vec2d(10, 0));
    }
	WHEN("set_origin to (10,0) and translate by (5, 5)") {
		gcodegen.set_origin(Vec2d(10,0));
		gcodegen.set_origin(gcodegen.origin() + Vec2d(5, 5));
		THEN("origin returns reference to point") {
    		REQUIRE(gcodegen.origin() == Vec2d(15,5));
    	}
    }
}

TEST_CASE("Generate regular polygon", "[GCode]") {
    const unsigned points_count{32};
    const Point centroid{scaled(Vec2d{5, -2})};
    const Polygon result{generate_regular_polygon(centroid, scaled(Vec2d{0, 0}), points_count)};
    const Point oposite_point{centroid * 2};

    REQUIRE(result.size() == 32);
    CHECK(result[16].x() == Approx(oposite_point.x()));
    CHECK(result[16].y() == Approx(oposite_point.y()));

    std::vector<double> angles;
    angles.reserve(points_count);
    for (unsigned index = 0; index < points_count; index++) {
        const unsigned previous_index{index == 0 ? points_count - 1 : index - 1};
        const unsigned next_index{index == points_count - 1 ? 0 : index + 1};

        const Point previous_point = result.points[previous_index];
        const Point current_point = result.points[index];
        const Point next_point = result.points[next_index];

        angles.emplace_back(angle(Vec2crd{previous_point - current_point}, Vec2crd{next_point - current_point}));
    }

    std::vector<double> expected;
    angles.reserve(points_count);
    std::generate_n(std::back_inserter(expected), points_count, [&](){
        return angles.front();
    });

    CHECK_THAT(angles, Catch::Matchers::Approx(expected));
}

TEST_CASE("Square bed with padding", "[GCode]") {
    const Bed bed{
        {
            Vec2d{0, 0},
            Vec2d{100, 0},
            Vec2d{100, 100},
            Vec2d{0, 100}
        },
        10.0
    };

    CHECK(bed.centroid.x() == 50);
    CHECK(bed.centroid.y() == 50);
    CHECK(bed.contains_within_padding(Vec2d{10, 10}));
    CHECK_FALSE(bed.contains_within_padding(Vec2d{9, 10}));

}
