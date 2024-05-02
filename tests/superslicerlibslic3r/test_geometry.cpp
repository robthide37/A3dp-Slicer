
#include <catch_main.hpp>

#include <libslic3r/Point.hpp>
#include <libslic3r/BoundingBox.hpp>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/Polyline.hpp>
#include <libslic3r/Line.hpp>
#include <libslic3r/Geometry.hpp>
#include <libslic3r/ClipperUtils.hpp>

using namespace Slic3r;

TEST_CASE("Polygon::contains works properly", ""){
   // this test was failing on Windows (GH #1950)
    Slic3r::Polygon polygon( Points{
        Point{207802834,-57084522},
        Point{196528149,-37556190},
        Point{173626821,-25420928},
        Point{171285751,-21366123},
        Point{118673592,-21366123},
        Point{116332562,-25420928},
        Point{93431208,-37556191},
        Point{82156517,-57084523},
        Point{129714478,-84542120},
        Point{160244873,-84542120}
    } );
    Point point{ 95706562, -57294774 };
    REQUIRE(polygon.contains(point));
}

SCENARIO("Intersections of line segments"){
    GIVEN("Integer coordinates"){
        Line line1{ Point{5,15},Point{30,15} };
        Line line2{ Point{10,20}, Point{10,10} };
        THEN("The intersection is valid"){
            Point point;
            line1.intersection(line2,&point);
            REQUIRE(Point{ 10,15 } == point);
        }
    }

    GIVEN("Scaled coordinates"){
        Line line1{ Point::new_scale(73.6310778185108, 371.74239268924), Point::new_scale(73.6310778185108, 501.74239268924) };
        Line line2{ Point::new_scale(75, 437.9853), Point::new_scale(62.7484, 440.4223) };
        THEN("There is still an intersection"){
            Point point;
            REQUIRE(line1.intersection(line2,&point));
        }
    }
}

/*
Tests for unused methods still written in perl
{
    my $polygon = Slic3r::Polygon->new(
        [45919000, 515273900], [14726100, 461246400], [14726100, 348753500], [33988700, 315389800], 
        [43749700, 343843000], [45422300, 352251500], [52362100, 362637800], [62748400, 369577600], 
        [75000000, 372014700], [87251500, 369577600], [97637800, 362637800], [104577600, 352251500], 
        [107014700, 340000000], [104577600, 327748400], [97637800, 317362100], [87251500, 310422300], 
        [82789200, 309534700], [69846100, 294726100], [254081000, 294726100], [285273900, 348753500], 
        [285273900, 461246400], [254081000, 515273900],
    );
    
    # this points belongs to $polyline
    # note: it's actually a vertex, while we should better check an intermediate point
    my $point = Slic3r::Point->new(104577600, 327748400);
    
    local $Slic3r::Geometry::epsilon = 1E-5;
    is_deeply Slic3r::Geometry::polygon_segment_having_point($polygon, $point)->pp, 
        [ [107014700, 340000000], [104577600, 327748400] ],
        'polygon_segment_having_point';
}
{
        auto point = Point{736310778.185108, 5017423926.8924};
        auto line = Line(Point{(long int} 627484000, (long int) 3695776000), Point{(long int} 750000000, (long int)3720147000));
        //is Slic3r::Geometry::point_in_segment($point, $line), 0, 'point_in_segment';
}

// Possible to delete
{
        //my $p1 = [10, 10];
        //my $p2 = [10, 20];
        //my $p3 = [10, 30];
        //my $p4 = [20, 20];
        //my $p5 = [0,  20];
        
        THEN("Points in a line give the correct angles"){
            //is Slic3r::Geometry::angle3points($p2, $p3, $p1),  PI(),   'angle3points';
            //is Slic3r::Geometry::angle3points($p2, $p1, $p3),  PI(),   'angle3points';
        }
        THEN("Left turns give the correct angle"){
            //is Slic3r::Geometry::angle3points($p2, $p4, $p3),  PI()/2, 'angle3points';
            //is Slic3r::Geometry::angle3points($p2, $p1, $p4),  PI()/2, 'angle3points';
        }
        THEN("Right turns give the correct angle"){
            //is Slic3r::Geometry::angle3points($p2, $p3, $p4),  PI()/2*3, 'angle3points';
            //is Slic3r::Geometry::angle3points($p2, $p1, $p5),  PI()/2*3, 'angle3points';
        }
        //my $p1 = [30, 30];
        //my $p2 = [20, 20];
        //my $p3 = [10, 10];
        //my $p4 = [30, 10];
        
        //is Slic3r::Geometry::angle3points($p2, $p1, $p3), PI(),       'angle3points';
        //is Slic3r::Geometry::angle3points($p2, $p1, $p4), PI()/2*3,   'angle3points';
        //is Slic3r::Geometry::angle3points($p2, $p1, $p1), 2*PI(),     'angle3points';
}

SCENARIO("polygon_is_convex works"){
    GIVEN("A square of dimension 10"){
        //my $cw_square = [ [0,0], [0,10], [10,10], [10,0] ];
        THEN("It is not convex clockwise"){
            //is polygon_is_convex($cw_square), 0, 'cw square is not convex';
        }
        THEN("It is convex counter-clockwise"){
            //is polygon_is_convex([ reverse @$cw_square ]), 1, 'ccw square is convex';
        } 

    }
    GIVEN("A concave polygon"){
        //my $convex1 = [ [0,0], [10,0], [10,10], [0,10], [0,6], [4,6], [4,4], [0,4] ];
        THEN("It is concave"){
            //is polygon_is_convex($convex1), 0, 'concave polygon';
        }
    }
}*/


TEST_CASE("Creating a polyline generates the obvious lines"){
    auto polyline = Slic3r::Polyline();
    polyline.points = std::vector<Point>({Point{0, 0}, Point{10, 0}, Point{20, 0}});
    REQUIRE(polyline.lines().at(0).a == Point{0,0});
    REQUIRE(polyline.lines().at(0).b == Point{10,0});
    REQUIRE(polyline.lines().at(1).a == Point{10,0});
    REQUIRE(polyline.lines().at(1).b == Point{20,0});
}

TEST_CASE("Splitting a Polygon generates a polyline correctly"){
    auto polygon = Slic3r::Polygon(std::vector<Point>({Point{0, 0}, Point{10, 0}, Point{5, 5}}));
    auto split = polygon.split_at_index(1);
    REQUIRE(split.points[0]==Point{10,0});
    REQUIRE(split.points[1]==Point{5,5});
    REQUIRE(split.points[2]==Point{0,0});
    REQUIRE(split.points[3]==Point{10,0});
}


TEST_CASE("Bounding boxes are scaled appropriately"){
    auto bb = BoundingBox(std::vector<Point>({Point{0, 1}, Point{10, 2}, Point{20, 2}}));
    bb.scale(2);
    REQUIRE(bb.min == Point{0,2});
    REQUIRE(bb.max == Point{40,4});
}


TEST_CASE("Offseting a line generates a polygon correctly"){
    Slic3r::Polyline tmp(Points{{10,10},{20,10} });
    Slic3r::Polygon area = offset(tmp,5).at(0);
    REQUIRE(area.area() == Slic3r::Polygon(std::vector<Point>({Point{10,5},Point{20,5},Point{20,15},Point{10,15}})).area());
}

SCENARIO("Circle Fit, TaubinFit with Newton's method") {
    GIVEN("A vector of Pointfs arranged in a half-circle with approximately the same distance R from some point") {
        Vec2d expected_center(-6, 0);
        Pointfs sample {Vec2d{6.0, 0}, Vec2d{5.1961524, 3}, Vec2d{3 ,5.1961524}, Vec2d{0, 6.0}, Vec2d{-3, 5.1961524}, Vec2d{-5.1961524, 3}, Vec2d{-6.0, 0}};
        std::transform(sample.begin(), sample.end(), sample.begin(), [expected_center] (const Vec2d& a) { return a + expected_center;});

        WHEN("Circle fit is called on the entire array") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample);
            THEN("A center point of -6,0 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
        WHEN("Circle fit is called on the first four points") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin(), sample.cbegin()+4);
            THEN("A center point of -6,0 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
        WHEN("Circle fit is called on the middle four points") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin()+2, sample.cbegin()+6);
            THEN("A center point of -6,0 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
    }
    GIVEN("A vector of Pointfs arranged in a half-circle with approximately the same distance R from some point") {
        Vec2d expected_center(-3, 9);
        Vec2ds sample {Vec2d{6.0, 0}, Vec2d{5.1961524, 3}, Vec2d{3 ,5.1961524},
                        Vec2d{0, 6.0}, 
                        Vec2d{3, 5.1961524}, Vec2d{-5.1961524, 3}, Vec2d{-6.0, 0}};

        std::transform(sample.begin(), sample.end(), sample.begin(), [expected_center] (const Vec2d& a) { return a + expected_center;});


        WHEN("Circle fit is called on the entire array") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample);
            THEN("A center point of 3,9 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
        WHEN("Circle fit is called on the first four points") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin(), sample.cbegin()+4);
            THEN("A center point of 3,9 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
        WHEN("Circle fit is called on the middle four points") {
            Vec2d result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin()+2, sample.cbegin()+6);
            THEN("A center point of 3,9 is returned.") {
                REQUIRE((result_center - expected_center).norm() < EPSILON);
            }
        }
    }
    GIVEN("A vector of Points arranged in a half-circle with approximately the same distance R from some point") {
        Point expected_center { Point::new_scale(-3, 9)};
        Points sample {Point::new_scale(6.0, 0), Point::new_scale(5.1961524, 3), Point::new_scale(3 ,5.1961524), 
                        Point::new_scale(0, 6.0), 
                        Point::new_scale(3, 5.1961524), Point::new_scale(-5.1961524, 3), Point::new_scale(-6.0, 0)};

        std::transform(sample.begin(), sample.end(), sample.begin(), [expected_center] (const Point& a) { return a + expected_center;});


        WHEN("Circle fit is called on the entire array") {
            Point result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample);
            THEN("A center point of scaled 3,9 is returned.") {
                REQUIRE(result_center.coincides_with_epsilon(expected_center));
            }
        }
        WHEN("Circle fit is called on the first four points") {
            Point result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin(), sample.cbegin()+4);
            THEN("A center point of scaled 3,9 is returned.") {
                REQUIRE(result_center.coincides_with_epsilon(expected_center));
            }
        }
        WHEN("Circle fit is called on the middle four points") {
            Point result_center(0,0);
            result_center = Geometry::circle_center_taubin_newton(sample.cbegin()+2, sample.cbegin()+6);
            THEN("A center point of scaled 3,9 is returned.") {
                REQUIRE(result_center.coincides_with_epsilon(expected_center));
            }
        }
    }
}

// A PU
//TEST_CASE("Chained path working correctly"){
//    // if chained_path() works correctly, these points should be joined with no diagonal paths
//    // (thus 26 units long)
//    std::vector<Point> points = {Point{26,26},Point{52,26},Point{0,26},Point{26,52},Point{26,0},Point{0,52},Point{52,52},Point{52,0}};
//    std::vector<Points::size_type> indices;
//    Geometry::chained_path(points,indices);
//    for(Points::size_type i = 0; i < indices.size()-1;i++){
//        double dist = points.at(indices.at(i)).distance_to(points.at(indices.at(i+1)));
//        REQUIRE(abs(dist-26) <= EPSILON);
//    }
//}

SCENARIO("Line distances"){
    GIVEN("A line"){
        Line line{ Point{0, 0}, Point{20, 0} };
        THEN("Points on the line segment have 0 distance"){
            REQUIRE(Point{0, 0}.distance_to(line)  == 0);
            REQUIRE(Point{20, 0}.distance_to(line) == 0);
            REQUIRE(Point{10, 0}.distance_to(line) == 0);
        
        }
        THEN("Points off the line have the appropriate distance"){
            REQUIRE(Point{10, 10}.distance_to(line) == 10);
            REQUIRE(Point{50, 0}.distance_to(line) == 30);
        }
    }
}

SCENARIO("Polygon convex/concave detection"){
    GIVEN(("A Square with dimension 100")){
        Slic3r::Polygon square/*new_scale*/{ std::vector<Point>{
            Point{100,100},
            Point{200,100},
            Point{200,200},
            Point{100,200}}};
        THEN("It has 4 convex points counterclockwise"){
            REQUIRE(square.concave_points(PI*4/3).size() == 0);
            REQUIRE(square.convex_points(PI*2/3).size() == 4);
        }
        THEN("It has 4 concave points clockwise"){
            square.make_clockwise();
            REQUIRE(square.concave_points(PI*4/3).size() == 4);
            REQUIRE(square.convex_points(PI*2/3).size() == 0);
        }
    }
    GIVEN("A Square with an extra colinearvertex"){
        Slic3r::Polygon square /*new_scale*/{ std::vector<Point>{
            Point{150,100},
            Point{200,100},
            Point{200,200},
            Point{100,200},
            Point{100,100}} };
        THEN("It has 4 convex points counterclockwise"){
            REQUIRE(square.concave_points(PI*4/3).size() == 0);
            REQUIRE(square.convex_points(PI*2/3).size() == 4);
        }
    }
    GIVEN("A Square with an extra collinear vertex in different order"){
        Slic3r::Polygon square /*new_scale*/{ std::vector<Point>{
            Point{200,200},
            Point{100,200},
            Point{100,100},
            Point{150,100},
            Point{200,100}} };
        THEN("It has 4 convex points counterclockwise"){
            REQUIRE(square.concave_points(PI*4/3).size() == 0);
            REQUIRE(square.convex_points(PI*2/3).size() == 4);
        }
    }

    GIVEN("A triangle"){
        Slic3r::Polygon triangle{ std::vector<Point>{
            Point{16000170,26257364},
            Point{714223,461012},
            Point{31286371,461008}
        } };
        THEN("it has three convex vertices"){
            REQUIRE(triangle.concave_points(PI*4/3).size() == 0);
            REQUIRE(triangle.convex_points(PI*2/3).size() == 3);
        }
    }

    GIVEN("A triangle with an extra collinear point"){
        Slic3r::Polygon triangle{ std::vector<Point>{
            Point{16000170,26257364},
            Point{714223,461012},
            Point{20000000,461012},
            Point{31286371,461012}
        } };
        THEN("it has three convex vertices"){
            REQUIRE(triangle.concave_points(PI*4/3).size() == 0);
            REQUIRE(triangle.convex_points(PI*2/3).size() == 3);
        }
    }
    GIVEN("A polygon with concave vertices with angles of specifically 4/3pi"){
        // Two concave vertices of this polygon have angle = PI*4/3, so this test fails
        // if epsilon is not used.
        Slic3r::Polygon polygon{ std::vector<Point>{
            Point{60246458,14802768},Point{64477191,12360001},
            Point{63727343,11060995},Point{64086449,10853608},
            Point{66393722,14850069},Point{66034704,15057334},
            Point{65284646,13758387},Point{61053864,16200839},
            Point{69200258,30310849},Point{62172547,42483120},
            Point{61137680,41850279},Point{67799985,30310848},
            Point{51399866,1905506},Point{38092663,1905506},
            Point{38092663,692699},Point{52100125,692699}
        } };
        THEN("the correct number of points are detected"){
            REQUIRE(polygon.concave_points(PI*4/3).size() == 6);
            REQUIRE(polygon.convex_points(PI*2/3).size() == 10);
        }
    }
}

TEST_CASE("Triangle Simplification does not result in less than 3 points"){
    Slic3r::Polygon triangle{ std::vector<Point>{
        Point{16000170,26257364}, Point{714223,461012}, Point{31286371,461008}
    } };
    REQUIRE(triangle.simplify(250000).at(0).points.size() == 3);
}

    
