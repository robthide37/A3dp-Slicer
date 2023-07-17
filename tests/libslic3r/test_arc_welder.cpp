#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/Geometry/ArcWelder.hpp>
#include <libslic3r/libslic3r.h>

using namespace Slic3r;

TEST_CASE("arc basics", "[ArcWelder]") {
    using namespace Slic3r::Geometry;

    WHEN("arc from { 2000.f, 1000.f } to { 1000.f, 2000.f }") {
        Vec2f p1{ 2000.f, 1000.f };
        Vec2f p2{ 1000.f, 2000.f };
        float r{ 1000.f };
        THEN("90 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, true);
            REQUIRE(is_approx(c, Vec2f{ 1000.f, 1000.f }));
            REQUIRE(ArcWelder::arc_angle(p1, p2, r) == Approx(0.5 * M_PI));
            REQUIRE(ArcWelder::arc_length(p1, p2, r) == Approx(r * 0.5 * M_PI).epsilon(0.001));
        }
        THEN("90 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, false);
            REQUIRE(is_approx(c, Vec2f{ 2000.f, 2000.f }));
        }
        THEN("270 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, true);
            REQUIRE(is_approx(c, Vec2f{ 2000.f, 2000.f }));
            REQUIRE(ArcWelder::arc_angle(p1, p2, - r) == Approx(1.5 * M_PI));
            REQUIRE(ArcWelder::arc_length(p1, p2, - r) == Approx(r * 1.5 * M_PI).epsilon(0.001));
        }
        THEN("270 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, false);
            REQUIRE(is_approx(c, Vec2f{ 1000.f, 1000.f }));
        }
    }
    WHEN("arc from { 1707.11f, 1707.11f } to { 1000.f, 2000.f }") {
        Vec2f p1{ 1707.11f, 1707.11f };
        Vec2f p2{ 1000.f, 2000.f };
        float r{ 1000.f };
        Vec2f center1 = Vec2f{ 1000.f, 1000.f };
        // Center on the other side of the CCW arch.
        Vec2f center2 = center1 + 2. * (0.5 * (p1 + p2) - center1);
        THEN("45 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, true);
            REQUIRE(is_approx(c, center1, 1.f));
            REQUIRE(ArcWelder::arc_angle(p1, p2, r) == Approx(0.25 * M_PI));
            REQUIRE(ArcWelder::arc_length(p1, p2, r) == Approx(r * 0.25 * M_PI).epsilon(0.001));
        }
        THEN("45 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, false);
            REQUIRE(is_approx(c, center2, 1.f));
        }
        THEN("315 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, true);
            REQUIRE(is_approx(c, center2, 1.f));
            REQUIRE(ArcWelder::arc_angle(p1, p2, - r) == Approx((2. - 0.25) * M_PI));
            REQUIRE(ArcWelder::arc_length(p1, p2, - r) == Approx(r * (2. - 0.25) * M_PI).epsilon(0.001));
        }
        THEN("315 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, false);
            REQUIRE(is_approx(c, center1, 1.f));
        }
    }
    WHEN("arc from { 1866.f, 1500.f } to { 1000.f, 2000.f }") {
        Vec2f p1{ 1866.f, 1500.f };
        Vec2f p2{ 1000.f, 2000.f };
        float r{ 1000.f };
        Vec2f center1 = Vec2f{ 1000.f, 1000.f };
        // Center on the other side of the CCW arch.
        Vec2f center2 = center1 + 2. * (0.5 * (p1 + p2) - center1);
        THEN("60 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, true);
            REQUIRE(is_approx(c, center1, 1.f));
            REQUIRE(is_approx(ArcWelder::arc_angle(p1, p2, r), float(M_PI / 3.), 0.001f));
            REQUIRE(ArcWelder::arc_length(p1, p2, r) == Approx(r * M_PI / 3.).epsilon(0.001));
        }
        THEN("60 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, r, false);
            REQUIRE(is_approx(c, center2, 1.f));
        }
        THEN("300 degrees arc, CCW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, true);
            REQUIRE(is_approx(c, center2, 1.f));
            REQUIRE(is_approx(ArcWelder::arc_angle(p1, p2, - r), float((2. - 1./3.) * M_PI), 0.001f));
            REQUIRE(ArcWelder::arc_length(p1, p2, - r) == Approx(r * (2. - 1. / 3.) * M_PI).epsilon(0.001));
        }
        THEN("300 degrees arc, CW") {
            Vec2f c = ArcWelder::arc_center(p1, p2, - r, false);
            REQUIRE(is_approx(c, center1, 1.f));
        }
    }
}

TEST_CASE("arc discretization", "[ArcWelder]") {
    using namespace Slic3r::Geometry;
    WHEN("arc from { 2, 1 } to { 1, 2 }") {
        const Point p1         = Point::new_scale(2., 1.);
        const Point p2         = Point::new_scale(1., 2.);
        const Point center     = Point::new_scale(1., 1.);
        const float radius     = scaled<float>(1.);
        const float resolution = scaled<float>(0.002);
        auto test = [center, resolution, radius](const Point &p1, const Point &p2, const float r, const bool ccw) {
            Vec2f  c = ArcWelder::arc_center(p1.cast<float>(), p2.cast<float>(), r, ccw);
            REQUIRE((p1.cast<float>() - c).norm() == Approx(radius));
            REQUIRE((c - center.cast<float>()).norm() == Approx(0.));
            Points pts = ArcWelder::arc_discretize(p1, p2, r, ccw, resolution);
            REQUIRE(pts.size() >= 2);
            REQUIRE(pts.front() == p1);
            REQUIRE(pts.back() == p2);
            for (const Point &p : pts)
                REQUIRE(std::abs((p.cast<double>() - c.cast<double>()).norm() - double(radius)) < double(resolution + SCALED_EPSILON));
        };
        THEN("90 degrees arc, CCW") {
            test(p1, p2, radius, true);
        }
        THEN("270 degrees arc, CCW") {
            test(p2, p1, - radius, true);
        }
        THEN("90 degrees arc, CW") {
            test(p2, p1, radius, false);
        }
        THEN("270 degrees arc, CW") {
            test(p1, p2, - radius, false);
        }
    }
}

TEST_CASE("arc fitting", "[ArcWelder]") {
    using namespace Slic3r::Geometry;

    WHEN("arc from { 2, 1 } to { 1, 2 }") {
        const Point p1         = Point::new_scale(2., 1.);
        const Point p2         = Point::new_scale(1., 2.);
        const Point center     = Point::new_scale(1., 1.);
        const float radius     = scaled<float>(1.);
        const float resolution = scaled<float>(0.002);
        auto test = [center, resolution](const Point &p1, const Point &p2, const float r, const bool ccw) {
            Points pts = ArcWelder::arc_discretize(p1, p2, r, ccw, resolution);
            ArcWelder::Path path = ArcWelder::fit_path(pts, resolution + SCALED_EPSILON, ArcWelder::default_scaled_resolution);
            REQUIRE(path.size() == 2);
            REQUIRE(path.front().point == p1);
            REQUIRE(path.front().radius == 0.f);
            REQUIRE(path.back().point == p2);
            REQUIRE(path.back().radius == Approx(r));
            REQUIRE(path.back().ccw() == ccw);
        };
        THEN("90 degrees arc, CCW is fitted") {
            test(p1, p2, radius, true);
        }
        THEN("270 degrees arc, CCW is fitted") {
            test(p2, p1, - radius, true);
        }
        THEN("90 degrees arc, CW is fitted") {
            test(p2, p1, radius, false);
        }
        THEN("270 degrees arc, CW is fitted") {
            test(p1, p2, - radius, false);
        }
    }

    WHEN("arc from { 2, 1 } to { 1, 2 }, another arc from { 2, 1 } to { 0, 2 }, tangentially connected") {
        const Point p1 = Point::new_scale(2., 1.);
        const Point p2 = Point::new_scale(1., 2.);
        const Point p3 = Point::new_scale(0., 3.);
        const Point center1 = Point::new_scale(1., 1.);
        const Point center2 = Point::new_scale(1., 3.);
        const float radius = scaled<float>(1.);
        const float resolution = scaled<float>(0.002);
        auto test = [center1, center2, resolution](const Point &p1, const Point &p2, const Point &p3, const float r, const bool ccw) {
            Points pts = ArcWelder::arc_discretize(p1, p2, r, ccw, resolution);
            {
                Points pts2 = ArcWelder::arc_discretize(p2, p3, - r, ! ccw, resolution);
                REQUIRE(pts.back() == pts2.front());
                pts.insert(pts.end(), std::next(pts2.begin()), pts2.end());
            }
            ArcWelder::Path path = ArcWelder::fit_path(pts, resolution + SCALED_EPSILON, ArcWelder::default_scaled_resolution);
            REQUIRE(path.size() == 3);
            REQUIRE(path.front().point == p1);
            REQUIRE(path.front().radius == 0.f);
            REQUIRE(path[1].point == p2);
            REQUIRE(path[1].radius == Approx(r));
            REQUIRE(path[1].ccw() == ccw);
            REQUIRE(path.back().point == p3);
            REQUIRE(path.back().radius == Approx(- r));
            REQUIRE(path.back().ccw() == ! ccw);
        };
        THEN("90 degrees arches, CCW are fitted") {
            test(p1, p2, p3, radius, true);
        }
        THEN("270 degrees arc, CCW is fitted") {
            test(p3, p2, p1, -radius, true);
        }
        THEN("90 degrees arc, CW is fitted") {
            test(p3, p2, p1, radius, false);
        }
        THEN("270 degrees arc, CW is fitted") {
            test(p1, p2, p3, -radius, false);
        }
    }
}

TEST_CASE("arc wedge test", "[ArcWelder]") {
    using namespace Slic3r::Geometry;

    WHEN("test point inside wedge, arc from { 2, 1 } to { 1, 2 }") {
        const int64_t s  = 1000000;
        const Vec2i64 p1{ 2 * s, s };
        const Vec2i64 p2{ s, 2 * s };
        const Vec2i64 center{ s, s };
        const int64_t radius{ s };
        auto test = [center](
            // Arc data
            const Vec2i64 &p1, const Vec2i64 &p2, const int64_t r, const bool ccw,
            // Test data
            const Vec2i64 &ptest, const bool ptest_inside) {
            const Vec2d c = ArcWelder::arc_center(p1.cast<double>(), p2.cast<double>(), double(r), ccw);
            REQUIRE(is_approx(c, center.cast<double>()));
            REQUIRE(ArcWelder::inside_arc_wedge(p1, p2, center, r > 0, ccw, ptest) == ptest_inside);
            REQUIRE(ArcWelder::inside_arc_wedge(p1.cast<double>(), p2.cast<double>(), double(r), ccw, ptest.cast<double>()) == ptest_inside);
        };
        auto test_quadrants = [center, test](
            // Arc data
            const Vec2i64 &p1, const Vec2i64 &p2, const int64_t r, const bool ccw,
            // Test data
            const Vec2i64 &ptest1, const bool ptest_inside1,
            const Vec2i64 &ptest2, const bool ptest_inside2, 
            const Vec2i64 &ptest3, const bool ptest_inside3,
            const Vec2i64 &ptest4, const bool ptest_inside4) {
            test(p1, p2, r, ccw, ptest1 + center, ptest_inside1);
            test(p1, p2, r, ccw, ptest2 + center, ptest_inside2);
            test(p1, p2, r, ccw, ptest3 + center, ptest_inside3);
            test(p1, p2, r, ccw, ptest4 + center, ptest_inside4);
        };
        THEN("90 degrees arc, CCW") {
            test_quadrants(p1, p2, radius, true, 
                Vec2i64{   s,   s }, true,
                Vec2i64{   s, - s }, false,
                Vec2i64{ - s,   s }, false,
                Vec2i64{ - s, - s }, false);
        }
        THEN("270 degrees arc, CCW") {
            test_quadrants(p2, p1, -radius, true,
                Vec2i64{   s,   s }, false,
                Vec2i64{   s, - s }, true,
                Vec2i64{ - s,   s }, true,
                Vec2i64{ - s, - s }, true);
        }
        THEN("90 degrees arc, CW") {
            test_quadrants(p2, p1, radius, false,
                Vec2i64{   s,   s }, true,
                Vec2i64{   s, - s }, false,
                Vec2i64{ - s,   s }, false,
                Vec2i64{ - s, - s }, false);
        }
        THEN("270 degrees arc, CW") {
            test_quadrants(p1, p2, -radius, false,
                Vec2i64{   s,   s }, false,
                Vec2i64{   s, - s }, true,
                Vec2i64{ - s,   s }, true,
                Vec2i64{ - s, - s }, true);
        }
    }
}
