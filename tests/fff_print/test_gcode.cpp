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

struct ApproxEqualsPoints : public Catch::MatcherBase<Points> {
    ApproxEqualsPoints(const Points& expected, unsigned tolerance): expected(expected), tolerance(tolerance) {}
    bool match(const Points& points) const override {
        if (points.size() != expected.size()) {
            return false;
        }
        for (auto i = 0u; i < points.size(); ++i) {
            const Point& point = points[i];
            const Point& expected_point = this->expected[i];
            if (
                std::abs(point.x() - expected_point.x()) > this->tolerance
                || std::abs(point.y() - expected_point.y()) > this->tolerance
            ) {
                return false;
            }
        }
        return true;
    }
    std::string describe() const override {
        std::stringstream ss;
        ss << std::endl;
        for (const Point& point : expected) {
            ss << "(" << point.x() << ", " << point.y() << ")" << std::endl;
        }
        ss << "With tolerance: " << this->tolerance;

        return "Equals " + ss.str();
    }

private:
    Points expected;
    unsigned tolerance;
};

Points get_points(const std::vector<DistancedPoint>& result) {
    Points result_points;
    std::transform(
        result.begin(),
        result.end(),
        std::back_inserter(result_points),
        [](const DistancedPoint& point){
            return point.point;
        }
    );
    return result_points;
}

std::vector<double> get_distances(const std::vector<DistancedPoint>& result) {
    std::vector<double> result_distances;
    std::transform(
        result.begin(),
        result.end(),
        std::back_inserter(result_distances),
        [](const DistancedPoint& point){
            return point.distance_from_start;
        }
    );
    return result_distances;
}

TEST_CASE("Place points at distances - expected use", "[GCode]") {
    std::vector<Point> line{
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{1, 0}),
        scaled(Vec2f{2, 1}),
        scaled(Vec2f{2, 2})
    };
    std::vector<double> distances{0, 0.2, 0.5, 1 + std::sqrt(2)/2, 1 + std::sqrt(2) + 0.5, 100.0};
    std::vector<DistancedPoint> result = slice_xy_path(line, distances);

    REQUIRE_THAT(get_points(result), ApproxEqualsPoints(Points{
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{0.2, 0}),
        scaled(Vec2f{0.5, 0}),
        scaled(Vec2f{1, 0}),
        scaled(Vec2f{1.5, 0.5}),
        scaled(Vec2f{2, 1}),
        scaled(Vec2f{2, 1.5}),
        scaled(Vec2f{2, 2})
    }, 5));

    REQUIRE_THAT(get_distances(result), Catch::Matchers::Approx(std::vector<double>{
        distances[0], distances[1], distances[2], 1, distances[3], 1 + std::sqrt(2), distances[4], 2 + std::sqrt(2)
    }));
}

TEST_CASE("Place points at distances - edge case", "[GCode]") {
    std::vector<Point> line{
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{1, 0}),
        scaled(Vec2f{2, 0})
    };
    std::vector<double> distances{0, 1, 1.5, 2};
    Points result{get_points(slice_xy_path(line, distances))};
    CHECK(result == Points{
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{1, 0}),
        scaled(Vec2f{1.5, 0}),
        scaled(Vec2f{2, 0})
    });
}

TEST_CASE("Generate elevated travel", "[GCode]") {
    std::vector<Point> xy_path{
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{1, 0}),
    };
    std::vector<double> ensure_points_at_distances{0.2, 0.5};
    Points3 result{generate_elevated_travel(xy_path, ensure_points_at_distances, 2.0, [](double x){return 1 + x;})};

    CHECK(result == Points3{
        scaled(Vec3f{0, 0, 3.0}),
        scaled(Vec3f{0.2, 0, 3.2}),
        scaled(Vec3f{0.5, 0, 3.5}),
        scaled(Vec3f{1, 0, 4.0})
    });
}
