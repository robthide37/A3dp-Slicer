#include <catch2/catch.hpp>

#include "libslic3r/Arachne/WallToolPaths.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/SVG.hpp"
#include "libslic3r/Utils.hpp"

using namespace Slic3r;
using namespace Slic3r::Arachne;

//#define ARACHNE_DEBUG_OUT

#ifdef ARACHNE_DEBUG_OUT
static void draw_extrusion(const std::string &path, const Polygons &polygons, const std::vector<VariableWidthLines> &vlines, const ExPolygons &contours)
{
    coordf_t    stroke_width = scale_(0.03);
    BoundingBox bbox         = get_extents(polygons);
    bbox.offset(scale_(1.));
    ::Slic3r::SVG svg(path.c_str(), bbox);

    svg.draw(contours, "cyan");

    for (const VariableWidthLines &vl : vlines)
        for (const ExtrusionLine &el : vl) {
            ThickPolyline thick_polyline = to_thick_polyline(el);
            svg.draw({thick_polyline}, "green", "blue", stroke_width);
        }

    for (const Line &line : to_lines(polygons))
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
    draw_extrusion(debug_out_path("arachne-closed-extrusion-line.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
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
    draw_extrusion(debug_out_path("arachne-missing-perimeter-8472.svg"), polygons, perimeters, union_ex(wallToolPaths.getInnerContour()));
#endif

    REQUIRE(perimeters.size() == 3);
}

