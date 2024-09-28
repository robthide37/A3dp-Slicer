
//#define CATCH_CONFIG_DISABLE

#include <catch_main.hpp>
#include "test_data.hpp"
#include <libslic3r/Fill/Fill.hpp>
#include <libslic3r/Print.hpp>
#include <libslic3r/ExtrusionEntity.hpp>
#include <libslic3r/Layer.hpp>
#include <libslic3r/Geometry.hpp>
#include <libslic3r/Flow.hpp>
#include <libslic3r/ClipperUtils.hpp>
#include <libslic3r/SVG.hpp>
#include <libslic3r/Format/3mf.hpp>

using namespace Slic3r;
using namespace Slic3r::Geometry;
using namespace Slic3r::Test;

bool test_if_solid_surface_filled(const ExPolygon& expolygon, double flow_spacing, double angle = 0, double density = 1.0);

//TEST_CASE("Fill: adjusted solid distance") {
//    Print print;
//    int surface_width {250};
//
//    int distance {Slic3r::Flow::solid_spacing(surface_width, 47)};
//
//    REQUIRE(distance == Approx(50));
//    REQUIRE(surface_width % distance == 0);
//}
Polylines test(const ExPolygon& poly, Fill &filler, const FillParams &params){
	Surface surface{ Slic3r::Surface((stPosTop | stDensSolid), poly) };
	return filler.fill_surface(&surface, params);
}

TEST_CASE("Fill: Pattern Path Length") {
    Fill* filler {Slic3r::Fill::new_from_type("rectilinear")};
	filler->angle = -(PI) / 2.0;
    FillParams params{};
    params.dont_adjust = true;
    params.density = 0.1; // 5/50
    filler->set_bounding_box(BoundingBox(Point(0, 0), Point::new_scale(Point(100, 100))));
    filler->init_spacing(5, params);
	//params.endpoints_overlap = false;



    SECTION("Square") {
        Points test_set{};
        test_set.reserve(4);
        Points points {Point{0,0}, Point{100,0}, Point{100,100}, Point{0,100}};
        for (size_t i = 0; i < 4; ++i) {
            std::transform(points.cbegin()+i, points.cend(),   std::back_inserter(test_set), [] (const Point& a) -> Point { return Point::new_scale(a); } ); 
            std::transform(points.cbegin(), points.cbegin()+i, std::back_inserter(test_set), [] (const Point& a) -> Point { return Point::new_scale(a); } ); 
			Slic3r::ExPolygon expoly{};
            expoly.contour = Slic3r::Polygon{ test_set };
            Polylines paths {test(expoly, *filler, params)};
            REQUIRE(paths.size() == 1); // one continuous path

            // TODO: determine what the "Expected length" should be for rectilinear fill of a 100x100 polygon. 
            // This check only checks that it's above scale(3*100 + 2*50) + scaled_epsilon.
            // ok abs($paths->[0]->length - scale(3*100 + 2*50)) - scaled_epsilon, 'path has expected length';
            REQUIRE(std::abs(paths[0].length() - static_cast<double>(scale_(3*100 + 2*50))) - SCALED_EPSILON > 0); // path has expected length

            test_set.clear();
        }
    }
    SECTION("Diamond with endpoints on grid") {
        Points points {Point{0,0}, Point{100,0}, Point{150,50}, Point{100,100}, Point{0,100}, Point{-50,50}};
        Points test_set{};
        test_set.reserve(6);
        std::transform(points.cbegin(), points.cend(),   std::back_inserter(test_set), [] (const Point& a) -> Point { return Point::new_scale(a); } );
		Slic3r::ExPolygon expoly;
		expoly.contour = Slic3r::Polygon(test_set);
        Polylines paths {test(expoly, *filler, params)};
        REQUIRE(paths.size() == 1); // one continuous path
    }

    SECTION("Square with hole") {
        Points square { Point{0,0}, Point{100,0}, Point{100,100}, Point{0,100}};
        Points hole {Point{25,25}, Point{75,25}, Point{75,75}, Point{25,75} };
        std::reverse(hole.begin(), hole.end());

        Points test_hole{};
        Points test_square{};

        std::transform(square.cbegin(), square.cend(), std::back_inserter(test_square), [] (const Point& a) -> Point { return Point::new_scale(a); } ); 
        std::transform(hole.cbegin(), hole.cend(), std::back_inserter(test_hole), [] (const Point& a) -> Point { return Point::new_scale(a); } ); 

        for (double angle : {-(PI/2.0), -(PI/4.0), -(PI), PI/2.0, PI}) {
            for (double spacing : {25.0, 5.0, 7.5, 8.5}) {
				FillParams params_local = params;
				params_local.density = filler->get_spacing() / spacing;
                filler->angle = angle;
				Slic3r::ExPolygon e{};
				e.contour = Slic3r::Polygon(test_square);
                e.holes = Slic3r::Polygons{Slic3r::Polygon(test_hole)};
                Polylines paths {test(e, *filler, params_local)};
				//std::cout << "paths.size="<<paths.size() << "\n";
				//{
				//    std::stringstream stri;
				//    stri << "squarewithhole.svg";
				//    SVG svg(stri.str());
				//    svg.draw(paths);
				//    svg.draw(e);
				//    svg.Close();
				//}
                //path CAN loop around the hole
                REQUIRE(paths.size() >= 1);
                REQUIRE(paths.size() <= 3);
                // paths don't cross hole
                REQUIRE(diff_pl(paths, offset(e, (float)(+SCALED_EPSILON * 10))).size() == 0);
            }
        }
    }
    SECTION("Regression: Missing infill segments in some rare circumstances") {
		FillParams params_local = params;
        params_local.density = 1;
        params_local.dont_adjust = false;
		Fill* filler_local = { Slic3r::Fill::new_from_type("rectilinear") };
		filler_local->angle = (PI/4.0);
        filler_local->set_bounding_box(BoundingBox(Point(0, 0), Point(2512749, 2512749)));
		filler_local->init_spacing(0.654498, params_local);
        //filler_local->endpoints_overlap = unscale(359974);
        filler_local->layer_id = 66;
        filler_local->z = 20.15;

        Points points {Point{25771516,14142125},Point{14142138,25771515},Point{2512749,14142131},Point{14142125,2512749}};
		Slic3r::ExPolygon expoly{};
		expoly.contour = Slic3r::Polygon(points);
        Polylines paths {test(expoly, *filler_local, params_local)};
        REQUIRE(paths.size() == 1); // one continuous path

        // TODO: determine what the "Expected length" should be for rectilinear fill of a 100x100 polygon. 
        // This check only checks that it's above scale(3*100 + 2*50) + scaled_epsilon.
        // ok abs($paths->[0]->length - scale(3*100 + 2*50)) - scaled_epsilon, 'path has expected length';
        REQUIRE(std::abs(paths[0].length() - static_cast<double>(scale_(3*100 + 2*50))) - SCALED_EPSILON > 0); // path has expected length
    }

    SECTION("Rotated Square") {
        Points square { Point::new_scale(0,0), Point::new_scale(50,0), Point::new_scale(50,50), Point::new_scale(0,50)};
		ExPolygon expolygon{};
		expolygon.contour = Slic3r::Polygon(square);
        auto filler {Slic3r::Fill::new_from_type("rectilinear")};
        filler->bounding_box = expolygon.contour.bounding_box();
        filler->angle = 0.F;
        
        Surface surface {(stPosTop|stDensSolid), expolygon};
        Flow flow = Flow::new_from_width(0.69f, 0.4f, 0.50f, 1.f, false); //width, height, nozzle_diameter, spacing_ratio, is_bridge

        params.density = 1.0;
        filler->init_spacing(flow.spacing(), params);

        for (auto angle : { 0.0, 45.0}) {
            surface.expolygon.rotate(angle, Point{0,0});
            Polylines paths = filler->fill_surface(&surface, params);
            REQUIRE(paths.size() == 1);
        }
    }
    SECTION("Solid surface fill") {
        Points points {
            Point(6883102, 9598327),
            Point(6883102, 20327272),
            Point(3116896, 20327272),
            Point(3116896, 9598327) 
        };
        Slic3r::ExPolygon expolygon{};
        expolygon.contour = Slic3r::Polygon{ points };
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        for (size_t i = 0; i <= 20; ++i)
        {
            expolygon.scale(1.05);
            //FIXME number overflow.
            REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        }
    }
    SECTION("Solid surface fill") {
        Points points {
                Point{59515297,5422499},Point{59531249,5578697},Point{59695801,6123186},
                Point{59965713,6630228},Point{60328214,7070685},Point{60773285,7434379},
                Point{61274561,7702115},Point{61819378,7866770},Point{62390306,7924789},
                Point{62958700,7866744},Point{63503012,7702244},Point{64007365,7434357},
                Point{64449960,7070398},Point{64809327,6634999},Point{65082143,6123325},
                Point{65245005,5584454},Point{65266967,5422499},Point{66267307,5422499},
                Point{66269190,8310081},Point{66275379,17810072},Point{66277259,20697500},
                Point{65267237,20697500},Point{65245004,20533538},Point{65082082,19994444},
                Point{64811462,19488579},Point{64450624,19048208},Point{64012101,18686514},
                Point{63503122,18415781},Point{62959151,18251378},Point{62453416,18198442},
                Point{62390147,18197355},Point{62200087,18200576},Point{61813519,18252990},
                Point{61274433,18415918},Point{60768598,18686517},Point{60327567,19047892},
                Point{59963609,19493297},Point{59695865,19994587},Point{59531222,20539379},
                Point{59515153,20697500},Point{58502480,20697500},Point{58502480,5422499}
        };
		Slic3r::ExPolygon expolygon;
        expolygon.contour = Slic3r::Polygon{ points };
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55, PI/2.0) == true);
    }
    SECTION("Solid surface fill") {
        Points points {
            Point::new_scale(0,0),Point::new_scale(98,0),Point::new_scale(98,10), Point::new_scale(0,10)
        };
		Slic3r::ExPolygon expolygon{};
        expolygon.contour = Slic3r::Polygon{ points };
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.5, 45.0, 0.99) == true);
    }

}

TEST_CASE("Fill area: check if periemter give the good values")
{
        Model model{};
        TriangleMesh sample_mesh = make_cube(5, 5, 0.2);
        double volume = (5 * 5 * 0.2);
        DynamicPrintConfig &config = Slic3r::DynamicPrintConfig::full_print_config();
        config.set_key_value("perimeters", new ConfigOptionInt(1));
        config.set_key_value("top_solid_layers", new ConfigOptionInt(1));
        config.set_key_value("bottom_solid_layers", new ConfigOptionInt(1));

        config.set_key_value("enforce_full_fill_volume", new ConfigOptionBool(false));
        config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0.1, false));

        config.set_key_value("skirts", new ConfigOptionInt(0));

        config.set_key_value("layer_height", new ConfigOptionFloat(0.2)); // get a known number of layers
        config.set_key_value("first_layer_height", new ConfigOptionFloatOrPercent(0.2, false));

        config.set_key_value("extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("first_layer_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("external_perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("solid_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("top_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        SECTION("no overlap")
        {
            config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0, false));
            config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(0));
            config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(0));
            Print print{};
            Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
            print.process();
            const LayerRegion *lr = print.get_object(0)->get_layer(0)->regions()[0];
            double area_infill = unscaled(unscaled(lr->fill_surfaces.surfaces.front().area()));
            REQUIRE(lr->fill_no_overlap_expolygons.empty());
            double area_computed = (5-0.5*2) * (5-0.5*2);
            REQUIRE(std::abs(area_computed - area_infill) < 0.001);
        }
        SECTION("only encroachment (0.2mm)")
        {
            config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0.2, false));
            config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(0));
            config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(0));
            Print print{};
            Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
            print.process();
            const LayerRegion *lr = print.get_object(0)->get_layer(0)->regions()[0];
            REQUIRE(1 == lr->fill_surfaces.surfaces.size());
            REQUIRE(1 == lr->fill_no_overlap_expolygons.size());
            double area_infill = unscaled(unscaled(lr->fill_surfaces.surfaces[0].area()));
            double area_infill_no_overlap = unscaled(unscaled(lr->fill_no_overlap_expolygons[0].area())); // note: don't need to intersect as there is only one fill_surfaces
            double area_no_encroach_computed = (5-0.5*2) * (5-0.5*2);
            double area_computed = (5-0.3*2) * (5-0.3*2);
            REQUIRE(area_infill_no_overlap < area_infill);
            REQUIRE(std::abs(area_computed - area_infill) < 0.001);
            REQUIRE(std::abs(area_no_encroach_computed - area_infill_no_overlap) < 0.001);
        }
        SECTION("only encroachment (40%)")
        {
            // % over (perimeter_spacing + solid_fill_spacing)/2, but no periemter overlap 
            // (note: here it's the external perimeter, as we have only one perimeter)
            config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(40, true)); // 40% -> 0.2 value
            config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(0));
            config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(0));
            Print print{};
            Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
            print.process();
            const LayerRegion *lr = print.get_object(0)->get_layer(0)->regions()[0];
            REQUIRE(1 == lr->fill_surfaces.surfaces.size());
            REQUIRE(1 == lr->fill_no_overlap_expolygons.size());
            double area_infill = unscaled(unscaled(lr->fill_surfaces.surfaces[0].area()));
            double area_infill_no_overlap = unscaled(unscaled(lr->fill_no_overlap_expolygons[0].area())); // note: don't need to intersect as there is only one fill_surfaces
            double area_no_encroach_computed = (5-0.5*2) * (5-0.5*2);
            double area_computed = (5-0.3*2) * (5-0.3*2);
            REQUIRE(area_infill_no_overlap < area_infill);
            REQUIRE(std::abs(area_computed - area_infill) < 0.001);
            REQUIRE(std::abs(area_no_encroach_computed - area_infill_no_overlap) < 0.001);
        }
        SECTION("only overlap")
        {
            config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0, false));
            config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(100));
            Print print{};
            Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
            print.process();
            const LayerRegion *lr = print.get_object(0)->get_layer(0)->regions()[0];
            REQUIRE(1 == lr->fill_surfaces.surfaces.size());
            REQUIRE(lr->fill_no_overlap_expolygons.empty());
            double area_infill = unscaled(unscaled(lr->fill_surfaces.surfaces[0].area()));
            double spacing_diff = (0.5f - Flow::rounded_rectangle_extrusion_spacing(0.5f, 0.2f, 1.f))/2;
            double area_computed = (5-(0.5-spacing_diff)*2) * (5-(0.5-spacing_diff)*2);
            REQUIRE(std::abs(area_computed - area_infill) < 0.001);
        }
        SECTION("both")
        {
            config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0.1, false));
            config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
            config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(100));
            Print print{};
            Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
            print.process();
            const LayerRegion *lr = print.get_object(0)->get_layer(0)->regions()[0];
            REQUIRE(1 == lr->fill_surfaces.surfaces.size());
            REQUIRE(1 == lr->fill_no_overlap_expolygons.size());
            double area_infill = unscaled(unscaled(lr->fill_surfaces.surfaces[0].area()));
            double area_infill_no_overlap = unscaled(unscaled(lr->fill_no_overlap_expolygons[0].area())); // note: don't need to intersect as there is only one fill_surfaces
            double spacing_diff = (0.5f - Flow::rounded_rectangle_extrusion_spacing(0.5f, 0.2f, 1.f))/2;
            double area_no_encroach_computed = (5-(0.5-spacing_diff)*2) * (5-(0.5-spacing_diff)*2);
            double area_computed = (5-(0.4-spacing_diff)*2) * (5-(0.4-spacing_diff)*2);
            REQUIRE(area_infill_no_overlap < area_infill);
            REQUIRE(std::abs(area_computed - area_infill) < 0.001);
            REQUIRE(std::abs(area_no_encroach_computed - area_infill_no_overlap) < 0.001);
        }
}

void test_all(DynamicPrintConfig &config, double& extrusion_width){

    
            SECTION("45°"){
                config.set_deserialize("fill_angle", "45");
                config.set_key_value("first_layer_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
                extrusion_width = 0.5;
                //test all solid fills
                SECTION("rectilinear") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinear)); }
                SECTION("rectilinear with gap fill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinearWGapFill)); }
                SECTION("ipMonotonic") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipMonotonic)); }
                SECTION("ipMonotonicWGapFill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipMonotonicWGapFill)); }
                SECTION("ipConcentric") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentric)); }
                SECTION("ipConcentricGapFill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentricGapFill)); }
                SECTION("ipHilbertCurve") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipHilbertCurve)); }
                SECTION("ipArchimedeanChords") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipArchimedeanChords)); }
                SECTION("ipOctagramSpiral") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipOctagramSpiral)); }
                SECTION("ipSmooth") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipSmooth)); }
            }
            SECTION("0° with bad spacing") {
                config.set_deserialize("fill_angle", "0");
                config.set_key_value("first_layer_extrusion_width", new ConfigOptionFloatOrPercent(0.415, false));
                extrusion_width = 0.415;
                SECTION("rectilinear") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinear)); }
                SECTION("rectilinear with gap fill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinearWGapFill)); }
                SECTION("ipMonotonic") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipMonotonic)); }
                SECTION("ipMonotonicWGapFill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipMonotonicWGapFill)); }
                SECTION("ipConcentric") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentric)); }
                SECTION("ipConcentricGapFill") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentricGapFill)); }
                SECTION("ipHilbertCurve") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipHilbertCurve)); }
                SECTION("ipArchimedeanChords") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipArchimedeanChords)); }
                SECTION("ipOctagramSpiral") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipOctagramSpiral)); }
                SECTION("ipSmooth") { config.set_key_value("bottom_fill_pattern", new ConfigOptionEnum<InfillPattern>(ipSmooth)); }
            }
}

#include "libslic3r/GCodeReader.hpp"
TEST_CASE("Fill: extrude gcode and check it")
{

    SECTION("simple square")
    {
        Model model{};
        TriangleMesh sample_mesh = make_cube(5, 5, 0.2);
        const double volume = (5 * 5 * 0.2);
        //sample_mesh.repair();

        DynamicPrintConfig &config = Slic3r::DynamicPrintConfig::full_print_config();
        config.set_key_value("perimeters", new ConfigOptionInt(1));
        config.set_key_value("top_solid_layers", new ConfigOptionInt(1));
        config.set_key_value("bottom_solid_layers", new ConfigOptionInt(1));

        config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0.1, false));
        config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(100));
        config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
        config.set_key_value("solid_infill_overlap", new ConfigOptionPercent(100));

        config.set_key_value("extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("first_layer_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("external_perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("solid_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("top_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        double extrusion_width = 0.5;
        config.set_deserialize("only_one_perimeter_top", "0");
        
        SECTION("classic"){
                config.set_key_value("perimeter_generator", new ConfigOptionEnum<PerimeterGeneratorType>(PerimeterGeneratorType::Classic));
                test_all(config, extrusion_width);
        }
        SECTION("arachne"){
                config.set_key_value("perimeter_generator", new ConfigOptionEnum<PerimeterGeneratorType>(PerimeterGeneratorType::Arachne));
                test_all(config, extrusion_width);
        }


        config.set_key_value("enforce_full_fill_volume", new ConfigOptionBool(true));
        config.set_key_value("skirts", new ConfigOptionInt(0));
        //simplier than auto opt = new ConfigOptionFloatsOrPercents{FloatOrPercent{0, false}}; opt.set_is_extruder_size(true);
        config.set_deserialize("seam_gap", "0");
        config.set_key_value("layer_height", new ConfigOptionFloat(0.2)); // get a known number of layers
        config.set_key_value("first_layer_height", new ConfigOptionFloatOrPercent(0.2, false));
        auto event_counter{ 0U };
        std::string stage;
        Print print{};
        Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
        print.process();

        std::string gcode_filepath{ "" };
        Slic3r::Test::gcode(gcode_filepath, print);
        //std::cout << "gcode generation done\n";
        std::string gcode_from_file = read_to_string(gcode_filepath);
        // model = print.model();
        // Slic3r::store_3mf("test.3mf", &model, &print.full_print_config(), OptionStore3mf{});

        //string[] lineArray = gcode_from_file
        GCodeReader parser;
        double volume_extruded = 0;
        int idx = 0;
        double volume_perimeter_extruded = 0;
        double volume_infill_extruded = 0;
        // add remaining time lines where needed
        parser.parse_buffer(gcode_from_file,
            [&](GCodeReader& reader, const GCodeReader::GCodeLine& line)
        {
            if (line.cmd_is("G1"))
            {
                if (line.dist_E(reader) > 0 && line.dist_XY(reader) > 0) {
                    //std::cout << "add " << line.dist_E(reader)<<" now "<< volume_extruded<<"=>";
                    volume_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                    //std::cout << volume_extruded << "\n";
                    if (idx<4)volume_perimeter_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                    else volume_infill_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                    idx++;
                }
            }
        });

        ExPolygons perimeter_center_line = offset_ex(print.get_object(0)->get_layer(0)->lslices[0], -scale_t(0.25f));

        //double perimeterRoundGapRemove = unscaled(print.get_object(0)->get_layer(0)->lslices[0].contour.length()) * 0.1*0.1 * (2 - (PI / 2));
        double perimeterRoundGapRemove = unscaled(perimeter_center_line[0].contour.length()) * 0.1*0.1 * (2 - (PI / 2));
        //double perimeterRoundGapAdd = unscaled(print.get_object(0)->get_layer(0)->lslices[0].contour.length()) * 0.1*0.1 * ((PI / 2));
        //for (Line &l : print.get_object(0)->get_layer(0)->slices.expolygons[0].contour.lines()) {

        //}
        //std::cout << "flow mm3permm: " << Flow{ 0.5f,0.2f,0.4f,false }.mm3_per_mm() << "\n";
        //std::cout << "perimeter : " << unscaled(print.get_object(0)->get_layer(0)->slices.expolygons[0].contour.length()) << " != " << (PI * 10) << "\n";

        //std::cout << "created a mesh of volume " << volume << " and i have extruded " << volume_extruded << " mm3.\n";
        //std::cout << "Note that if we remove the bits of the external extrusion, it's only a volume of " << (volume - perimeterRoundGapRemove) << " that needs to be filled\n";
        //std::cout << "Note that if we add the bits of the external extrusion, it's a volume of " << (volume + perimeterRoundGapAdd) << " that needs to be filled\n";

        double volumeExtrPerimeter = ExtrusionVolume{}.get(print.get_object(0)->get_layer(0)->regions()[0]->perimeters);
        double volumeExtrInfill = ExtrusionVolume{}.get(print.get_object(0)->get_layer(0)->regions()[0]->fills);

        double volumeInfill = 0;
        for (const ExPolygon & p : print.get_object(0)->get_layer(0)->regions()[0]->fill_no_overlap_expolygons) {
            volumeInfill += unscaled(unscaled(p.area()));
        }
        double spacing_diff = (extrusion_width - Flow::rounded_rectangle_extrusion_spacing(extrusion_width, 0.2f, 1.f))/2;
        double fill_raw_area_no_encroach = (5-(extrusion_width-spacing_diff)*2) * (5-(extrusion_width-spacing_diff)*2);
        double fill_raw_area = (5-(extrusion_width-0.1-spacing_diff)*2) * (5-(extrusion_width-0.1-spacing_diff)*2);

        double compute_perimeter_area = (5-extrusion_width)*4*Flow::rounded_rectangle_extrusion_spacing(extrusion_width, 0.2f, 1.f);

        std::cout << "area fill_no_overlap_expolygons= " << (unscaled(unscaled(print.get_object(0)->get_layer(0)->regions()[0]->fill_no_overlap_expolygons.front().contour.area()))) << "\n";
        volumeInfill *= 0.2;
        std::cout << "\nvolumeRealr=" << (volume_perimeter_extruded + volume_infill_extruded) << " volumeRealPerimeter= " << volume_perimeter_extruded << " and volumeRealInfill=" << volume_infill_extruded << " mm3." << "\n";
        std::cout << "volumeExtr=" << (volumeExtrPerimeter + volumeExtrInfill) << " volumeExtrPerimeter= " << volumeExtrPerimeter << " and volumeExtrInfill=" << volumeExtrInfill << " mm3." << "\n";
        std::cout << "volumePerimeter= " << (volume - volumeInfill) << " volumePerimeter(wo/bits)= " << (volume - volumeInfill- perimeterRoundGapRemove) << " and volumeInfill=" << volumeInfill << " mm3." << "\n";
        std::cout << "volume= " << (volume) << " raw_fill_volume="<<fill_raw_area*0.2<<" raw_fill_volume_no_encroach=" << fill_raw_area_no_encroach*0.2 << "\n";
        std::cout << "raw_fill_area="<<fill_raw_area<<" raw_fill_area_no_encroach=" << fill_raw_area_no_encroach << "\n";
        std::cout << "computed peri= " << (unscaled(perimeter_center_line[0].contour.length())*0.2*0.5 - 2*perimeterRoundGapRemove)<< " perimeterRoundGapRemove= " << (perimeterRoundGapRemove) << "\n";

        //Flow fl{0.5f, 0.2f, 0.4f, false};

        //{
        //    std::stringstream stri;
        //    stri << "extrusion_width_learning" << ".svg";
        //    SVG svg(stri.str());
        //    //svg.draw(bounds);
        //    svg.draw(print.get_object(0)->get_layer(0)->slices.expolygons[0].contour, "green");
        //    svg.draw(print.get_object(0)->get_layer(0)->regions()[0]->fill_no_overlap_expolygons, "black", scale_(0.01));
        //    svg.draw(print.get_object(0)->get_layer(0)->regions()[0]->perimeters.as_polylines(), "orange", fl.scaled_width());
        //    svg.draw(print.get_object(0)->get_layer(0)->regions()[0]->perimeters.as_polylines(), "red", fl.scaled_spacing());
        //    svg.draw(print.get_object(0)->get_layer(0)->regions()[0]->fills.as_polylines(), "cyan", fl.scaled_width());
        //    svg.draw(print.get_object(0)->get_layer(0)->regions()[0]->fills.as_polylines(), "blue", fl.scaled_spacing());
        //    svg.Close();
        //}
        REQUIRE(abs(fill_raw_area_no_encroach*0.2 - volumeInfill) < 0.01);
        REQUIRE(abs(compute_perimeter_area * 0.2 - volumeExtrPerimeter) < 0.01);

        //std::cout << gcode_from_file;
        if(abs(volumeInfill - volumeExtrInfill) > EPSILON*5) // *5 for archimean chords
            std::cout<<"stop";
        REQUIRE(abs(volumeInfill - volumeExtrInfill) < EPSILON * 5);// *5 for archimean chords
        REQUIRE(abs(volumeInfill - volume_infill_extruded) < 0.01);
        REQUIRE(abs((volume - volumeInfill - perimeterRoundGapRemove) - volumeExtrPerimeter) < 0.01);
        REQUIRE(abs((volume - volumeInfill - perimeterRoundGapRemove) - volume_perimeter_extruded) < 0.1); //there are a bit less for seam mitigation 
        // lower than the full volume because of the rounded extenral perimeter
        REQUIRE(volume_extruded < volume);
        //lower than the full volume - rounded external perimeter (because I used 4*5 as perimeter length instead of (5-offset)*4 )
        REQUIRE(volume_extruded > volume - ((5*4) * 0.2 * (0.5f-Flow::rounded_rectangle_extrusion_spacing(extrusion_width, 0.2f, 1.f))/2));
        clean_file(gcode_filepath, "gcode");

    }
    SECTION("simple disk") {
        Model model{};
        TriangleMesh sample_mesh = make_cylinder(5, 0.2);
        const double volume = (PI * 25 * 0.2);

        DynamicPrintConfig &config = Slic3r::DynamicPrintConfig::full_print_config();
        config.set_key_value("perimeters", new ConfigOptionInt(1));
        config.set_key_value("top_solid_layers", new ConfigOptionInt(1));
        config.set_key_value("bottom_solid_layers", new ConfigOptionInt(1));

        config.set_key_value("enforce_full_fill_volume", new ConfigOptionBool(true));
        config.set_key_value("infill_overlap", new ConfigOptionFloatOrPercent(0.1, true));
        config.set_key_value("perimeter_overlap", new ConfigOptionPercent(100));
        config.set_key_value("external_perimeter_overlap", new ConfigOptionPercent(100));
        config.set_deserialize("external_perimeter_cut_corners", "0");

        config.set_key_value("skirts", new ConfigOptionInt(0));

        config.set_key_value("layer_height", new ConfigOptionFloat(0.2)); // get a known number of layers
        config.set_key_value("first_layer_height", new ConfigOptionFloatOrPercent(0.2, false));

        config.set_key_value("extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("first_layer_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("external_perimeter_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("solid_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        config.set_key_value("top_infill_extrusion_width", new ConfigOptionFloatOrPercent(0.5, false));
        auto event_counter{ 0U };
        std::string stage;
        Print print{};
        Slic3r::Test::init_print(print, { sample_mesh }, model, &config);
        print.process();

        std::string gcode_filepath{ "" };
        Slic3r::Test::gcode(gcode_filepath, print);
        //std::cout << "gcode generation done\n";
        std::string gcode_from_file = read_to_string(gcode_filepath);

        //string[] lineArray = gcode_from_file
        GCodeReader parser;
        double volume_extruded = 0;
        //int idx = 0;
        int step = 0;
        double volume_perimeter_extruded = 0;
        double volume_infill_extruded = 0;
        // add remaining time lines where needed
        parser.parse_buffer(gcode_from_file,
            [&](GCodeReader& reader, const GCodeReader::GCodeLine& line)
        {
            if(line.comment() == "TYPE:External perimeter")
                step = 1;
            if(line.comment() == "TYPE:Solid infill")
                step = 2;
            if (line.cmd_is("G1") && step > 0)
            {
                if (line.dist_E(reader) > 0 && line.dist_XY(reader) > 0) {
                    //std::cout << "add " << line.dist_E(reader)<<" now "<< volume_extruded<<"=>";
                    volume_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                    //std::cout << volume_extruded << "\n";
                    if (step == 1) volume_perimeter_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                    else if (step == 2) volume_infill_extruded += line.dist_E(reader)*(PI*1.75*1.75 / 4.);
                }
            }
        });

        ExPolygons perimeter_center_line = offset_ex(print.get_object(0)->get_layer(0)->lslices[0], -scale_t(0.25f));

        //double perimeterRoundGapRemove = unscaled(print.get_object(0)->get_layer(0)->lslices[0].contour.length()) * 0.1*0.1 * (2 - (PI / 2));
        double perimeterRoundGapRemove = unscaled(perimeter_center_line[0].contour.length()) * 0.1*0.1 * (2 - (PI / 2));
        //double perimeterRoundGapAdd = unscaled(print.get_object(0)->get_layer(0)->lslices[0].contour.length()) * 0.1*0.1 * ((PI / 2));

        double volumeExtrPerimeter = ExtrusionVolume{}.get(print.get_object(0)->get_layer(0)->regions()[0]->perimeters);
        double volumeExtrInfill = ExtrusionVolume{}.get(print.get_object(0)->get_layer(0)->regions()[0]->fills);

        double volumeInfill = 0;
        ExPolygons infill_area = intersection_ex(print.get_object(0)->get_layer(0)->regions()[0]->fill_no_overlap_expolygons, print.get_object(0)->get_layer(0)->regions()[0]->fill_expolygons);
        for (const ExPolygon & p : infill_area) {
            volumeInfill += unscaled(unscaled(p.area()));
        }
        volumeInfill *= 0.2;
        std::cout << "volumeRealr=" << (volume_perimeter_extruded + volume_infill_extruded) << " volumeRealPerimeter= " << volume_perimeter_extruded << " and volumeRealInfill=" << volume_infill_extruded << " mm3." << "\n";
        std::cout << "volumeExtr=" << (volumeExtrPerimeter + volumeExtrInfill) << " volumeExtrPerimeter= " << volumeExtrPerimeter << " and volumeExtrInfill=" << volumeExtrInfill << " mm3." << "\n";
        std::cout << "volumePerimeter= " << (volume - volumeInfill) << " volumePerimeter(wo/bits)= " << (volume - volumeInfill - perimeterRoundGapRemove) << " and volumeInfill=" << volumeInfill << " mm3." << "\n";
        std::cout << "volume= " << (volume) << "\n";

        REQUIRE(abs(volumeInfill - volumeExtrInfill) < 0.001);
        REQUIRE(abs(volumeInfill - volume_infill_extruded) < 0.001);
        REQUIRE(abs((volume - volumeInfill - perimeterRoundGapRemove) - volumeExtrPerimeter) < 0.1); //there are a bit less for seam mitigation 
        REQUIRE(abs(volumeExtrPerimeter - volume_perimeter_extruded) < 0.01);
        clean_file(gcode_filepath, "gcode");

    }
}

/* 

{
    my $collection = Slic3r::Polyline::Collection->new(
            Slic3r::Polyline->new([0,15], [0,18], [0,20]),
            Slic3r::Polyline->new([0,10], [0,8], [0,5]),
            );
    is_deeply
        [ map $_->[Y], map @$_, @{$collection->chained_path_from(Slic3r::Point->new(0,30), 0)} ],
        [20, 18, 15, 10, 8, 5],
        'chained path';
}

{
    my $collection = Slic3r::Polyline::Collection->new(
            Slic3r::Polyline->new([4,0], [10,0], [15,0]),
            Slic3r::Polyline->new([10,5], [15,5], [20,5]),
            );
    is_deeply
        [ map $_->[X], map @$_, @{$collection->chained_path_from(Slic3r::Point->new(30,0), 0)} ],
        [reverse 4, 10, 15, 10, 15, 20],
        'chained path';
}

{
    my $collection = Slic3r::ExtrusionPath::Collection->new(
            map Slic3r::ExtrusionPath->new(polyline => $_, role => 0, mm3_per_mm => 1),
            Slic3r::Polyline->new([0,15], [0,18], [0,20]),
            Slic3r::Polyline->new([0,10], [0,8], [0,5]),
            );
    is_deeply
        [ map $_->[Y], map @{$_->polyline}, @{$collection->chained_path_from(Slic3r::Point->new(0,30), 0)} ],
        [20, 18, 15, 10, 8, 5],
        'chained path';
}

{
    my $collection = Slic3r::ExtrusionPath::Collection->new(
            map Slic3r::ExtrusionPath->new(polyline => $_, role => 0, mm3_per_mm => 1),
            Slic3r::Polyline->new([15,0], [10,0], [4,0]),
            Slic3r::Polyline->new([10,5], [15,5], [20,5]),
            );
    is_deeply
        [ map $_->[X], map @{$_->polyline}, @{$collection->chained_path_from(Slic3r::Point->new(30,0), 0)} ],
        [reverse 4, 10, 15, 10, 15, 20],
        'chained path';
}

for my $pattern (qw(rectilinear honeycomb hilbertcurve concentric)) {
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('fill_pattern', $pattern);
    $config->set('external_fill_pattern', $pattern);
    $config->set('perimeters', 1);
    $config->set('skirts', 0);
    $config->set('fill_density', 20);
    $config->set('layer_height', 0.05);
    $config->set('perimeter_extruder', 1);
    $config->set('infill_extruder', 2);
    my $print = Slic3r::Test::init_print('20mm_cube', config => $config, scale => 2);
    ok my $gcode = Slic3r::Test::gcode($print), "successful $pattern infill generation";
    my $tool = undef;
    my @perimeter_points = my @infill_points = ();
    Slic3r::GCode::Reader->new->parse($gcode, sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            if ($tool == $config->perimeter_extruder-1) {
            push @perimeter_points, Slic3r::Point->new_scale($args->{X}, $args->{Y});
            } elsif ($tool == $config->infill_extruder-1) {
            push @infill_points, Slic3r::Point->new_scale($args->{X}, $args->{Y});
            }
            }
            });
    my $convex_hull = convex_hull(\@perimeter_points);
    ok !(defined first { !$convex_hull->contains_point($_) } @infill_points), "infill does not exceed perimeters ($pattern)";
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('infill_only_where_needed', 1);
    $config->set('bottom_solid_layers', 0);
    $config->set('infill_extruder', 2);
    $config->set('infill_extrusion_width', 0.5);
    $config->set('fill_density', 40);
    $config->set('cooling', 0);                 # for preventing speeds from being altered
        $config->set('first_layer_speed', '100%');  # for preventing speeds from being altered

        my $test = sub {
            my $print = Slic3r::Test::init_print('pyramid', config => $config);

            my $tool = undef;
            my @infill_extrusions = ();  # array of polylines
                Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
                        my ($self, $cmd, $args, $info) = @_;

                        if ($cmd =~ /^T(\d+)/) {
                        $tool = $1;
                        } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
                        if ($tool == $config->infill_extruder-1) {
                        push @infill_extrusions, Slic3r::Line->new_scale(
                                [ $self->X, $self->Y ],
                                [ $info->{new_X}, $info->{new_Y} ],
                                );
                        }
                        }
                        });
            return 0 if !@infill_extrusions;  # prevent calling convex_hull() with no points

                my $convex_hull = convex_hull([ map $_->pp, map @$_, @infill_extrusions ]);
            return unscale unscale sum(map $_->area, @{offset([$convex_hull], scale(+$config->infill_extrusion_width/2))});
        };

    my $tolerance = 5;  # mm^2

        $config->set('solid_infill_below_area', 0);
    ok $test->() < $tolerance,
       'no infill is generated when using infill_only_where_needed on a pyramid';

    $config->set('solid_infill_below_area', 70);
    ok abs($test->() - $config->solid_infill_below_area) < $tolerance,
       'infill is only generated under the forced solid shells';
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('skirts', 0);
    $config->set('perimeters', 1);
    $config->set('fill_density', 0);
    $config->set('top_solid_layers', 0);
    $config->set('bottom_solid_layers', 0);
    $config->set('solid_infill_below_area', 20000000);
    $config->set('solid_infill_every_layers', 2);
    $config->set('perimeter_speed', 99);
    $config->set('external_perimeter_speed', 99);
    $config->set('cooling', 0);
    $config->set('first_layer_speed', '100%');

    my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
    my %layers_with_extrusion = ();
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd eq 'G1' && $info->{dist_XY} > 0 && $info->{extruding}) {
            if (($args->{F} // $self->F) != $config->perimeter_speed*60) {
            $layers_with_extrusion{$self->Z} = ($args->{F} // $self->F);
            }
            }
            });

    ok !%layers_with_extrusion,
       "solid_infill_below_area and solid_infill_every_layers are ignored when fill_density is 0";
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('skirts', 0);
    $config->set('perimeters', 3);
    $config->set('fill_density', 0);
    $config->set('layer_height', 0.2);
    $config->set('first_layer_height', 0.2);
    $config->set('nozzle_diameter', [0.35]);
    $config->set('infill_extruder', 2);
    $config->set('solid_infill_extruder', 2);
    $config->set('infill_extrusion_width', 0.52);
    $config->set('solid_infill_extrusion_width', 0.52);
    $config->set('first_layer_extrusion_width', 0);

    my $print = Slic3r::Test::init_print('A', config => $config);
    my %infill = ();  # Z => [ Line, Line ... ]
        my $tool = undef;
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            if ($tool == $config->infill_extruder-1) {
            my $z = 1 * $self->Z;
            $infill{$z} ||= [];
            push @{$infill{$z}}, Slic3r::Line->new_scale(
                    [ $self->X, $self->Y ],
                    [ $info->{new_X}, $info->{new_Y} ],
                    );
            }
            }
            });
    my $grow_d = scale($config->infill_extrusion_width)/2;
    my $layer0_infill = union([ map @{$_->grow($grow_d)}, @{ $infill{0.2} } ]);
    my $layer1_infill = union([ map @{$_->grow($grow_d)}, @{ $infill{0.4} } ]);
    my $diff = diff($layer0_infill, $layer1_infill);
    $diff = offset2_ex($diff, -$grow_d, +$grow_d);
    $diff = [ grep { $_->area > 2*(($grow_d*2)**2) } @$diff ];
    is scalar(@$diff), 0, 'no missing parts in solid shell when fill_density is 0';
}

{
    # GH: #2697
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('perimeter_extrusion_width', 0.72);
    $config->set('top_infill_extrusion_width', 0.1);
    $config->set('infill_extruder', 2);         # in order to distinguish infill
        $config->set('solid_infill_extruder', 2);   # in order to distinguish infill

        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
    my %infill = ();  # Z => [ Line, Line ... ]
        my %other  = ();  # Z => [ Line, Line ... ]
        my $tool = undef;
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            my $z = 1 * $self->Z;
            my $line = Slic3r::Line->new_scale(
                    [ $self->X, $self->Y ],
                    [ $info->{new_X}, $info->{new_Y} ],
                    );
            if ($tool == $config->infill_extruder-1) {
            $infill{$z} //= [];
            push @{$infill{$z}}, $line;
            } else {
            $other{$z} //= [];
            push @{$other{$z}}, $line;
            }
            }
            });
    my $top_z = max(keys %infill);
    my $top_infill_grow_d = scale($config->top_infill_extrusion_width)/2;
    my $top_infill = union([ map @{$_->grow($top_infill_grow_d)}, @{ $infill{$top_z} } ]);
    my $perimeters_grow_d = scale($config->perimeter_extrusion_width)/2;
    my $perimeters = union([ map @{$_->grow($perimeters_grow_d)}, @{ $other{$top_z} } ]);
    my $covered = union_ex([ @$top_infill, @$perimeters ]);
    my @holes = map @{$_->holes}, @$covered;
    ok sum(map unscale unscale $_->area*-1, @holes) < 1, 'no gaps between top solid infill and perimeters';
}
*/

class ExtrusionGetExtrusionArea : public ExtrusionVisitor {
    Polygons grown_paths;
public:
    ExtrusionGetExtrusionArea() {}
    void use(ExtrusionPath &path) override {
        polygons_append(grown_paths, offset(path.as_polyline().as_polyline(), scale_t(path.width)));
    }
    void use(ExtrusionPath3D &path3D) override { assert(false); }
    void use(ExtrusionMultiPath &multipath) override { for (ExtrusionPath path : multipath.paths) path.visit(*this);    }
    void use(ExtrusionMultiPath3D &multipath) override { for (ExtrusionPath path : multipath.paths) path.visit(*this);    }
    void use(ExtrusionLoop &loop) override { for (ExtrusionPath path : loop.paths) path.visit(*this); }
    void use(ExtrusionEntityCollection &collection) override { for (ExtrusionEntity *entity : collection.entities()) entity->visit(*this); }
    Polygons get(ExtrusionEntityCollection &coll) {
        for (ExtrusionEntity *entity : coll.entities()) entity->visit(*this);
        return grown_paths;
    }
};
//TODO: also check by volume extruded
//TODO: replace the simple area coverage check by one that takes into account the width of the path, not only the default flow spacing
//TODO: test more fills
bool test_if_solid_surface_filled(const ExPolygon& expolygon, double flow_width, double angle, double density) {
    auto* filler {Slic3r::Fill::new_from_type("concentricgapfill")};
    filler->bounding_box = expolygon.contour.bounding_box();
    filler->angle = angle;
	FillParams params;
	params.dont_adjust = false;
    FullPrintConfig config= FullPrintConfig::defaults();
    params.config = &config;

    Surface surface((stPosBottom | stDensSolid), expolygon);
	//note: here we do flow.width = flow_width , flow.gheight = 0.4, flow.nozzle_size = flow_width;
    params.flow = Flow::new_from_width( float(flow_width), 0.4, float(flow_width), 1.f, false);

    params.density = density;
    filler->init_spacing(params.flow.spacing(), params);

    // concentricgapfill can't output only Polylines, as it's a composed thing with gapfill
    //Polylines paths {filler->fill_surface(&surface, params)};
    ExtrusionEntityCollection coll;
    filler->fill_surface_extrusion(&surface, params, coll.set_entities());

    // check whether any part was left uncovered

// figure out what is actually going on here re: data types
    Polygons grown_paths = ExtrusionGetExtrusionArea{}.get(coll);
    //grown_paths.reserve(paths.size());
    //std::for_each(paths.begin(), paths.end(), [filler, &grown_paths] (const Slic3r::Polyline& p) {
    //    polygons_append(grown_paths, offset(p, scale_(filler->get_spacing() / 2.0)));
    //});
    
    ExPolygons uncovered = diff_ex(expolygon, grown_paths, Slic3r::ApplySafetyOffset::Yes);

    // ignore very small dots
    const auto scaled_flow_width { std::pow(scale_(flow_width), 2) };
    auto iter {std::remove_if(uncovered.begin(), uncovered.end(), [scaled_flow_width] (const ExPolygon& poly) {
        return poly.area() > scaled_flow_width;
    }) };
    uncovered.erase(iter, uncovered.end());

	double uncovered_area = 0;
	for (ExPolygon &p : uncovered) uncovered_area += unscaled(unscaled(p.area()));
	std::cout << "uncovered size =" << uncovered_area << " / "<< unscaled(unscaled(expolygon.area()))<<"\n";
    return uncovered_area < 0.05; // solid surface is (almost) fully filled

}
