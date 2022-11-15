#include <catch2/catch.hpp>

#include "test_data.hpp"
#include "clipper/clipper_z.hpp"
#include "libslic3r/clipper.hpp"

using namespace Slic3r;

// tests for ExPolygon::overlaps(const ExPolygon &other)
SCENARIO("Clipper intersection with polyline", "[Clipper]")
{
    struct TestData { 
        ClipperLib::Path subject;
        ClipperLib::Path clip;
        ClipperLib::Paths result;
    };

    auto run_test = [](const TestData &data) {
        ClipperLib::Clipper clipper;
        clipper.AddPath(data.subject, ClipperLib::ptSubject, false);
        clipper.AddPath(data.clip, ClipperLib::ptClip, true);

        ClipperLib::PolyTree polytree;
        ClipperLib::Paths    paths;
        clipper.Execute(ClipperLib::ctIntersection, polytree, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
        ClipperLib::PolyTreeToPaths(polytree, paths);

        REQUIRE(paths == data.result);
    };

    WHEN("Open polyline completely inside stays inside") {
        run_test({
            { { 10, 0 }, { 20, 0 } },
            { { -1000, -1000 }, { -1000,  1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { { { 20, 0 }, { 10, 0 } } }
        });
    };
    WHEN("Closed polyline completely inside stays inside") {
        run_test({
            { { 10, 0 }, { 20, 0 }, { 20, 20 }, { 10, 20 }, { 10, 0 } },
            { { -1000, -1000 }, { -1000,  1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { { { 10, 0 }, { 20, 0 }, { 20, 20 }, { 10, 20 }, { 10, 0 } } }
        });
    };
    WHEN("Polyline which crosses right rectangle boundary is trimmed") {
        run_test({
            { { 10, 0 }, { 2000, 0 } },
            { { -1000, -1000 }, { -1000,  1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { { { 1000, 0 }, { 10, 0 } } }
        });
    };
    WHEN("Polyline which is outside clipping region is removed") {
        run_test({
            { { 1500, 0 }, { 2000, 0 } },
            { { -1000, -1000 }, { -1000,  1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { }
        });
    };

    WHEN("Polyline on left vertical boundary is kept") {
        run_test({
            { { -1000, -1000 }, { -1000, 1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { { { -1000, -1000 }, { -1000, 1000 } } }
        });
        run_test({
            { { -1000, 1000 }, { -1000, -1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { { { -1000, 1000 }, { -1000, -1000 } } }
        });
    };
    WHEN("Polyline on right vertical boundary is kept") {
        run_test({
            { { 1000, -1000 }, { 1000, 1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { { { 1000, -1000 }, { 1000, 1000 } } }
        });
        run_test({
            { { 1000, 1000 }, { 1000, -1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000,  1000 }, { 1000, -1000 } },
            { { { 1000, 1000 }, { 1000, -1000 } } }
        });
    };
    WHEN("Polyline on bottom horizontal boundary is removed") {
        run_test({
            { { -1000, -1000 }, { 1000, -1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { }
        });
        run_test({
            { { 1000, -1000 }, { -1000, -1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { }
            });
    };
    WHEN("Polyline on top horizontal boundary is removed") {
        run_test({
            { { -1000, 1000 }, { 1000, 1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { }
        });
        run_test({
            { { 1000, 1000 }, { -1000, 1000 } },
            { { -1000, -1000 }, { -1000, 1000 }, { 1000, 1000 }, { 1000, -1000 } },
            { }
            });
    };
}

SCENARIO("Clipper Z", "[ClipperZ]")
{
    ClipperLib_Z::Path subject { { -2000, -1000, 10 }, { -2000,  1000, 10 }, { 2000,  1000, 10 }, { 2000, -1000, 10 } };
    ClipperLib_Z::Path clip{ { -1000, -2000, -5 }, { -1000,  2000, -5 }, { 1000,  2000, -5 }, { 1000, -2000, -5 } };

    ClipperLib_Z::Clipper clipper;
    clipper.ZFillFunction([](const ClipperLib_Z::IntPoint &e1bot, const ClipperLib_Z::IntPoint &e1top, const ClipperLib_Z::IntPoint &e2bot,
                             const ClipperLib_Z::IntPoint &e2top, ClipperLib_Z::IntPoint &pt) {
        pt.z() = 1;
    });

    clipper.AddPath(subject, ClipperLib_Z::ptSubject, false);
    clipper.AddPath(clip, ClipperLib_Z::ptClip, true);

    ClipperLib_Z::PolyTree polytree;
    ClipperLib_Z::Paths    paths;
    clipper.Execute(ClipperLib_Z::ctIntersection, polytree, ClipperLib_Z::pftNonZero, ClipperLib_Z::pftNonZero);
    ClipperLib_Z::PolyTreeToPaths(polytree, paths);

    REQUIRE(paths.size() == 1);
    REQUIRE(paths.front().size() == 2);
    for (const ClipperLib_Z::IntPoint &pt : paths.front())
        REQUIRE(pt.z() == 1);
}

