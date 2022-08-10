#include <catch2/catch.hpp>

#include "libslic3r/Arachne/WallToolPaths.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/SVG.hpp"
#include "libslic3r/Utils.hpp"

using namespace Slic3r;
using namespace Slic3r::Arachne;

//#define ARACHNE_DEBUG_OUT

#ifdef ARACHNE_DEBUG_OUT
static void export_perimeters_to_svg(const std::string &path, const Polygons &contours, const std::vector<Arachne::VariableWidthLines> &perimeters, const ExPolygons &infill_area)
{
    coordf_t    stroke_width = scale_(0.03);
    BoundingBox bbox         = get_extents(contours);
    bbox.offset(scale_(1.));
    ::Slic3r::SVG svg(path.c_str(), bbox);

    svg.draw(infill_area, "cyan");

    for (const Arachne::VariableWidthLines &perimeter : perimeters)
        for (const Arachne::ExtrusionLine &extrusion_line : perimeter) {
            ThickPolyline thick_polyline = to_thick_polyline(extrusion_line);
            svg.draw({thick_polyline}, "green", "blue", stroke_width);
        }

    for (const Line &line : to_lines(contours))
        svg.draw(line, "red", stroke_width);
}
#endif

TEST_CASE("Arachne - Closed ExtrusionLine", "[ArachneClosedExtrusionLine]") {
    Polygon poly = {
        Point(62478540, -7411873),  Point(62478540, 9978540),   Point(-62478540, 9978540),  Point(-62478540, -7411873), Point(-58818049, -7411874),
        Point(-58639424, -7393054), Point(-58430204, -7325743), Point(-58317958, -7261069), Point(-58187096, -7150294), Point(-58032997, -6934055),
        Point(-57956770, -6723830), Point(-57927922, -6536131), Point(-57948215, -6249353), Point(-58038066, -5971432), Point(-58400146, -5419566),
        Point(-58720844, -4711417), Point(-58812247, -4429032), Point(-58945877, -3868129), Point(-58999054, -3458536), Point(-59021495, -3000104),
        Point(-58978088, -2345755), Point(-58945641, -2130557), Point(-58719408, -1284462), Point(-58350699, -492550),  Point(-57854519, 218384),
        Point(-57690839, 403070),   Point(-57242241, 834472),   Point(-56937894, 1068372),  Point(-56522699, 1341801),  Point(-56245930, 1488656),
        Point(-55633586, 1748152),  Point(-54872819, 1945077),  Point(-54279560, 2011550),  Point(-53999789, 2021482),  Point(-53550613, 1995780),
        Point(-53127364, 1945077),  Point(-52852060, 1886312),  Point(-52262142, 1711199),  Point(-51479386, 1343079),  Point(-50763215, 839141),
        Point(-50302640, 384003),   Point(-50150730, 224721),   Point(-49616021, -551391),  Point(-49279548, -1287513), Point(-49210805, -1492871),
        Point(-49054178, -2131559), Point(-49004097, -2510916), Point(-48978506, -2999863), Point(-49012690, -3563440), Point(-49054263, -3868963),
        Point(-49280289, -4714703), Point(-49649307, -5507428), Point(-49921685, -5909300), Point(-49982145, -6037568), Point(-50058406, -6311890),
        Point(-50071436, -6450122), Point(-50045218, -6724784), Point(-50006267, -6852153), Point(-49876677, -7082551), Point(-49687434, -7260679),
        Point(-49451449, -7373245), Point(-49184639, -7411873), Point(-13258854, -7411871), Point(-13258853, -7021460), Point(-44501672, -7021460),
        Point(-44816622, -6971390), Point(-45100558, -6826592), Point(-45326629, -6600516), Point(-45471725, -6315721), Point(-45521460, -6001689),
        Point(-45521460, 3000836),  Point(-45509082, 3159381),  Point(-45412385, 3459576),  Point(-45326445, 3600421),  Point(-45221777, 3722975),
        Point(-44964550, 3909861),  Point(-44815338, 3971666),  Point(-44598502, 4016264),  Point(58501687, 4021460),   Point(58814884, 3971856),
        Point(58964551, 3909863),   Point(59221777, 3722976),   Point(59326668, 3600164),   Point(59436707, 3406438),   Point(59508907, 3160338),
        Point(59521460, 3000842),   Point(59521460, -6001688),  Point(59471724, -6315713),  Point(59326597, -6600557),  Point(59100555, -6826595),
        Point(58814822, -6971976),  Point(58501662, -7021460),  Point(27258850, -7021460),  Point(27258851, -7411871),  Point(59755385, -7411874),
    };

    Polygons polygons    = {poly};
    coord_t  spacing     = 407079;
    coord_t  inset_count = 8;

    Arachne::WallToolPaths wallToolPaths(polygons, spacing, spacing, inset_count, 0, 0.2, PrintObjectConfig::defaults(), PrintConfig::defaults());
    wallToolPaths.generate();
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
    export_perimeters_to_svg(debug_out_path("arachne-closed-extrusion-line.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif

    for (VariableWidthLines &perimeter : perimeters)
        for (ExtrusionLine &el : perimeter)
            if (el.is_closed) {
//                REQUIRE(el.junctions.front().p == el.junctions.back().p);
            }
}

// This test case was distilled from GitHub issue #8472.
// Where for wall_distribution_count == 3 sometime middle perimeter was missing.
TEST_CASE("Arachne - Missing perimeter - #8472", "[ArachneMissingPerimeter8472]") {
    Polygon poly = {
        Point(-9000000,  8054793),
        Point( 7000000,  8054793),
        Point( 7000000, 10211874),
        Point(-8700000, 10211874),
        Point(-9000000,  9824444)
    };

    Polygons polygons    = {poly};
    coord_t  spacing     = 437079;
    coord_t  inset_count = 3;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
    print_object_config.wall_distribution_count.setInt(3);

    Arachne::WallToolPaths wallToolPaths(polygons, spacing, spacing, inset_count, 0, 0.2, print_object_config, PrintConfig::defaults());
    wallToolPaths.generate();
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
    export_perimeters_to_svg(debug_out_path("arachne-missing-perimeter-8472.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif

    REQUIRE(perimeters.size() == 3);
}

// This test case was distilled from GitHub issue #8593.
// Where on the symmetrical model, there were missing parts of extrusions in gear teeth based on model rotation.
TEST_CASE("Arachne - #8593 - Missing a part of the extrusion", "[ArachneMissingPartOfExtrusion8593]") {
    const Polygon poly_orig = {
        Point( 1800000, 28500000),
        Point( 1100000, 30000000),
        Point( 1000000, 30900000),
        Point(  600000, 32300000),
        Point( -600000, 32300000),
        Point(-1000000, 30900000),
        Point(-1100000, 30000000),
        Point(-1800000, 29000000),
    };

    coord_t  spacing     = 377079;
    coord_t  inset_count = 3;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
    print_object_config.min_bead_width.set(new ConfigOptionFloatOrPercent(0.315, false));
    print_object_config.wall_transition_angle.set(new ConfigOptionFloat(40.));
    print_object_config.wall_transition_length.set(new ConfigOptionFloatOrPercent(1., false));

    // This behavior seems to be related to the rotation of the input polygon.
    // There are specific angles in which this behavior is always triggered.
    for (const double angle : {0., -PI / 2., -PI / 15.}) {
        Polygon poly = poly_orig;
        if (angle != 0.)
            poly.rotate(angle);

        Polygons polygons    = {poly};
        Arachne::WallToolPaths wall_tool_paths(polygons, spacing, spacing, inset_count, 0, 0.2, print_object_config, PrintConfig::defaults());
        wall_tool_paths.generate();
        std::vector<Arachne::VariableWidthLines> perimeters = wall_tool_paths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
        {
            static int iRun = 0;
            export_perimeters_to_svg(debug_out_path("arachne-missing-part-of-extrusion-8593-%d.svg", iRun++), polygons, perimeters, union_ex(wall_tool_paths.getInnerContour()));
        }
#endif
    }
}

// This test case was distilled from GitHub issue #8573.
TEST_CASE("Arachne - #8573 - A gap in the perimeter - 1", "[ArachneGapInPerimeter8573_1]") {
    const Polygon poly = {
        Point(13960000,  500000),
        Point(13920000, 1210000),
        Point(13490000, 2270000),
        Point(12960000, 3400000),
        Point(12470000, 4320000),
        Point(12160000, 4630000),
        Point(12460000, 3780000),
        Point(12700000, 2850000),
        Point(12880000, 1910000),
        Point(12950000, 1270000),
        Point(13000000,  500000),
    };

    Polygons polygons    = {poly};
    coord_t  spacing     = 407079;
    coord_t  inset_count = 2;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
//    print_object_config.wall_transition_angle.set(new ConfigOptionFloat(20.));

    Arachne::WallToolPaths wallToolPaths(polygons, spacing, spacing, inset_count, 0, 0.2, print_object_config, PrintConfig::defaults());
    wallToolPaths.generate();
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
    export_perimeters_to_svg(debug_out_path("arachne-gap-in-perimeter-1-8573.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif
}

// This test case was distilled from GitHub issue #8444.
TEST_CASE("Arachne - #8444 - A gap in the perimeter - 2", "[ArachneGapInPerimeter8444_2]") {
    const Polygon poly = {
        Point(14413938, 3825902),
        Point(16817613,  711749),
        Point(19653030,   67154),
        Point(20075592,  925370),
        Point(20245428, 1339788),
        Point(20493219, 2121894),
        Point(20570295, 2486625),
        Point(20616559, 2835232),
        Point(20631964, 3166882),
        Point(20591800, 3858877),
        Point(19928267, 2153012),
        Point(19723020, 1829802),
        Point(19482017, 1612364),
        Point(19344810, 1542433),
        Point(19200249, 1500902),
        Point(19047680, 1487200),
        Point(18631073, 1520777),
        Point(18377524, 1567627),
        Point(18132517, 1641174),
        Point(17896307, 1741360),
        Point(17669042, 1868075),
        Point(17449999, 2021790),
    };

    Polygons polygons    = {poly};
    coord_t  spacing     = 594159;
    coord_t  inset_count = 2;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
    //    print_object_config.wall_transition_angle.set(new ConfigOptionFloat(20.));

    Arachne::WallToolPaths wallToolPaths(polygons, spacing, spacing, inset_count, 0, 0.4, print_object_config, PrintConfig::defaults());
    wallToolPaths.generate();
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
    export_perimeters_to_svg(debug_out_path("arachne-gap-in-perimeter-2-8444.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif
}

// This test case was distilled from GitHub issue #8528.
// There is a hole in the place where the number of perimeters is changing from 6 perimeters to 7 perimeters.
TEST_CASE("Arachne - #8528 - A hole when number of perimeters is changing", "[ArachneHoleOnPerimetersChange8528]") {
    const Polygon poly = {
        Point(-30000000, 27650000),
        Point(-30000000, 33500000),
        Point(-40000000, 33500000),
        Point(-40500000, 33500000),
        Point(-41100000, 33400000),
        Point(-41600000, 33200000),
        Point(-42100000, 32900000),
        Point(-42600000, 32600000),
        Point(-43000000, 32200000),
        Point(-43300000, 31700000),
        Point(-43600000, 31200000),
        Point(-43800000, 30700000),
        Point(-43900000, 30100000),
        Point(-43900000, 29600000),
        Point(-43957080, 25000000),
        Point(-39042920, 25000000),
        Point(-39042920, 27650000),
    };

    Polygons polygons    = {poly};
    coord_t  spacing     = 814159;
    coord_t  inset_count = 5;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
    print_object_config.min_bead_width.set(new ConfigOptionFloatOrPercent(0.68, false));

    // Changing min_bead_width to 0.66 seems that resolve this issue, at least in this case.
    print_object_config.min_bead_width.set(new ConfigOptionFloatOrPercent(0.66, false));

    Arachne::WallToolPaths wallToolPaths(polygons, spacing, spacing, inset_count, 0, 0.4, print_object_config, PrintConfig::defaults());
    wallToolPaths.generate();
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
    export_perimeters_to_svg(debug_out_path("arachne-hole-on-perimeters-change-8528.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif
}

// This test case was distilled from GitHub issue #8528.
// There is an inconsistency between layers in length of the single perimeters.
TEST_CASE("Arachne - #8555 - Inconsistent single perimeter", "[ArachneInconsistentSinglePerimeter8555]") {
    const Polygon poly_0 = {
        Point(5527411, -38490007),
        Point(11118814, -36631169),
        Point(13529600, -36167120),
        Point(11300145, -36114514),
        Point(10484024, -36113916),
        Point(5037323, -37985945),
        Point(4097054, -39978866)
    };
    const Polygon poly_1 = {
        Point(5566841, -38517205),
        Point(11185208, -36649404),
        Point(13462719, -36211009),
        Point(11357290, -36161329),
        Point(10583855, -36160763),
        Point(5105952, -38043516),
        Point(4222019, -39917031)
    };
    const Polygon poly_2 = {
        Point(5606269, -38544404),
        Point(11251599, -36667638),
        Point(13391666, -36255700),
        Point(10683552, -36207653),
        Point(5174580, -38101085),
        Point(4346981, -39855197)
    };
    const Polygon poly_3 = {
        Point(5645699, -38571603),
        Point(11317993, -36685873),
        Point(13324786, -36299588),
        Point(10783383, -36254499),
        Point(5243209, -38158655),
        Point(4471947, -39793362)
    };
    const Polygon poly_4 = {
        Point(5685128, -38598801),
        Point(11384385, -36704108),
        Point(13257907, -36343476),
        Point(10883211, -36301345),
        Point(5311836, -38216224),
        Point(4596909, -39731528)
    };
    const Polygon poly_5 = {
        Point(5724558, -38626000),
        Point(11450778, -36722343),
        Point(13191026, -36387365),
        Point(10983042, -36348191),
        Point(5380466, -38273795),
        Point(4721874, -39669693)
    };

    Polygons polygons    = {poly_0, poly_1, poly_2, poly_3, poly_4, poly_5};
    coord_t  spacing     = 417809;
    coord_t  inset_count = 2;

    for (size_t poly_idx = 0; poly_idx < polygons.size(); ++poly_idx) {
        Polygons input_polygons{polygons[poly_idx]};
        Arachne::WallToolPaths wallToolPaths(input_polygons, spacing, spacing, inset_count, 0, 0.15, PrintObjectConfig::defaults(), PrintConfig::defaults());
        wallToolPaths.generate();
        std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
        export_perimeters_to_svg(debug_out_path("arachne-inconsistent-single-perimeter-8555-%d.svg", poly_idx), input_polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif
    }
}

// This test case was distilled from GitHub issue #8633.
// Open perimeter extrusion is shorter on endpoints in comparison to closed perimeter.
TEST_CASE("Arachne - #8633 - Shorter open perimeter", "[ArachneShorterOpenPerimeter8633]") {
    const Polygon poly_0 = {
        Point(6507498, 4189461),
        Point(6460382, 3601960),
        Point(6390896, 3181097),
        Point(6294072, 2765838),
        Point(6170293, 2357794),

        Point(7090581, 2045388),
        Point(7232821, 2514293),
        Point(7344089, 2991501),
        Point(7423910, 3474969),
        Point(7471937, 3962592),
        Point(7487443, 4436235),
        Point(6515575, 4436235),
    };

    const Polygon poly_1 = {
        Point(6507498, 4189461),
        Point(6460382, 3601960),
        Point(6390896, 3181097),
        Point(6294072, 2765838),
        Point(6170293, 2357794),

        Point(6917958, 1586830),
        Point(7090552, 2045398),

        Point(7232821, 2514293),
        Point(7344089, 2991501),
        Point(7423910, 3474969),
        Point(7471937, 3962592),
        Point(7487443, 4436235),
        Point(6515575, 4436235),
    };

    Polygons polygons    = {poly_0, poly_1};
    coord_t  spacing     = 617809;
    coord_t  inset_count = 1;

    PrintObjectConfig print_object_config = PrintObjectConfig::defaults();
    print_object_config.min_bead_width.set(new ConfigOptionFloatOrPercent(0.51, false));
    print_object_config.min_feature_size.set(new ConfigOptionFloatOrPercent(0.15, false));
    print_object_config.wall_transition_length.set(new ConfigOptionFloatOrPercent(0.6, false));

    for (size_t poly_idx = 0; poly_idx < polygons.size(); ++poly_idx) {
        Polygons input_polygons{polygons[poly_idx]};
        Arachne::WallToolPaths wallToolPaths(input_polygons, spacing, spacing, inset_count, 0, 0.15, print_object_config, PrintConfig::defaults());
        wallToolPaths.generate();
        std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#ifdef ARACHNE_DEBUG_OUT
        export_perimeters_to_svg(debug_out_path("arachne-shorter-open-perimeter-8633-%d.svg", poly_idx), input_polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif
    }
}