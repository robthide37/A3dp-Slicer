#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/Geometry/ArcWelder.hpp>

TEST_CASE("arc_center", "[ArcWelder]") {
    using namespace Slic3r;
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
