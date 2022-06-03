#include <unordered_set>
#include <unordered_map>
#include <random>
#include <numeric>
#include <cstdint>

#include "sla_test_utils.hpp"

#include <libslic3r/TriangleMeshSlicer.hpp>
#include <libslic3r/SLA/SupportTreeMesher.hpp>
#include <libslic3r/BranchingTree/PointCloud.hpp>

namespace {

const char *const BELOW_PAD_TEST_OBJECTS[] = {
    "20mm_cube.obj",
    "V.obj",
};

const char *const AROUND_PAD_TEST_OBJECTS[] = {
    "20mm_cube.obj",
    "V.obj",
    "frog_legs.obj",
    "cube_with_concave_hole_enlarged.obj",
};

const char *const SUPPORT_TEST_MODELS[] = {
    "cube_with_concave_hole_enlarged_standing.obj",
    "A_upsidedown.obj",
    "extruder_idler.obj"
};

} // namespace

TEST_CASE("Pillar pairhash should be unique", "[SLASupportGeneration]") {
    test_pairhash<int, int>();
    test_pairhash<int, long>();
    test_pairhash<unsigned, unsigned>();
    test_pairhash<unsigned, unsigned long>();
}

TEST_CASE("Support point generator should be deterministic if seeded", 
          "[SLASupportGeneration], [SLAPointGen]") {
    TriangleMesh mesh = load_model("A_upsidedown.obj");
    
    AABBMesh emesh{mesh};
    
    sla::SupportTreeConfig supportcfg;
    sla::SupportPointGenerator::Config autogencfg;
    autogencfg.head_diameter = float(2 * supportcfg.head_front_radius_mm);
    sla::SupportPointGenerator point_gen{emesh, autogencfg, [] {}, [](int) {}};
        
    auto   bb      = mesh.bounding_box();
    double zmin    = bb.min.z();
    double zmax    = bb.max.z();
    double gnd     = zmin - supportcfg.object_elevation_mm;
    auto   layer_h = 0.05f;
    
    auto slicegrid = grid(float(gnd), float(zmax), layer_h);
    std::vector<ExPolygons> slices = slice_mesh_ex(mesh.its, slicegrid, CLOSING_RADIUS);
    
    point_gen.seed(0);
    point_gen.execute(slices, slicegrid);
    
    auto get_chksum = [](const std::vector<sla::SupportPoint> &pts){
        int64_t chksum = 0;
        for (auto &pt : pts) {
            auto p = scaled(pt.pos);
            chksum += p.x() + p.y() + p.z();
        }
        
        return chksum;
    };
    
    int64_t checksum = get_chksum(point_gen.output());
    size_t ptnum = point_gen.output().size();
    REQUIRE(point_gen.output().size() > 0);
    
    for (int i = 0; i < 20; ++i) {
        point_gen.output().clear();
        point_gen.seed(0);
        point_gen.execute(slices, slicegrid);
        REQUIRE(point_gen.output().size() == ptnum);
        REQUIRE(checksum == get_chksum(point_gen.output()));
    }
}

TEST_CASE("Flat pad geometry is valid", "[SLASupportGeneration]") {
    sla::PadConfig padcfg;
    
    // Disable wings
    padcfg.wall_height_mm = .0;
    
    for (auto &fname : BELOW_PAD_TEST_OBJECTS) test_pad(fname, padcfg);
}

TEST_CASE("WingedPadGeometryIsValid", "[SLASupportGeneration]") {
    sla::PadConfig padcfg;
    
    // Add some wings to the pad to test the cavity
    padcfg.wall_height_mm = 1.;
    
    for (auto &fname : BELOW_PAD_TEST_OBJECTS) test_pad(fname, padcfg);
}

TEST_CASE("FlatPadAroundObjectIsValid", "[SLASupportGeneration]") {
    sla::PadConfig padcfg;
    
    // Add some wings to the pad to test the cavity
    padcfg.wall_height_mm = 0.;
    // padcfg.embed_object.stick_stride_mm = 0.;
    padcfg.embed_object.enabled = true;
    padcfg.embed_object.everywhere = true;
    
    for (auto &fname : AROUND_PAD_TEST_OBJECTS) test_pad(fname, padcfg);
}

TEST_CASE("WingedPadAroundObjectIsValid", "[SLASupportGeneration]") {
    sla::PadConfig padcfg;
    
    // Add some wings to the pad to test the cavity
    padcfg.wall_height_mm = 1.;
    padcfg.embed_object.enabled = true;
    padcfg.embed_object.everywhere = true;
    
    for (auto &fname : AROUND_PAD_TEST_OBJECTS) test_pad(fname, padcfg);
}

TEST_CASE("DefaultSupports::ElevatedSupportGeometryIsValid", "[SLASupportGeneration]") {
    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 10.;
    
    for (auto fname : SUPPORT_TEST_MODELS) test_supports(fname, supportcfg);
}

TEST_CASE("DefaultSupports::FloorSupportGeometryIsValid", "[SLASupportGeneration]") {
    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 0;
    
    for (auto &fname: SUPPORT_TEST_MODELS) test_supports(fname, supportcfg);
}

TEST_CASE("DefaultSupports::ElevatedSupportsDoNotPierceModel", "[SLASupportGeneration]") {
    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 10.;

    for (auto fname : SUPPORT_TEST_MODELS)
        test_support_model_collision(fname, supportcfg);
}

TEST_CASE("DefaultSupports::FloorSupportsDoNotPierceModel", "[SLASupportGeneration]") {
    
    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 0;
    
    for (auto fname : SUPPORT_TEST_MODELS)
        test_support_model_collision(fname, supportcfg);
}

//TEST_CASE("BranchingSupports::ElevatedSupportGeometryIsValid", "[SLASupportGeneration][Branching]") {
//    sla::SupportTreeConfig supportcfg;
//    supportcfg.object_elevation_mm = 10.;
//    supportcfg.tree_type = sla::SupportTreeType::Branching;

//    for (auto fname : SUPPORT_TEST_MODELS) test_supports(fname, supportcfg);
//}

//TEST_CASE("BranchingSupports::FloorSupportGeometryIsValid", "[SLASupportGeneration][Branching]") {
//    sla::SupportTreeConfig supportcfg;
//    supportcfg.object_elevation_mm = 0;
//    supportcfg.tree_type = sla::SupportTreeType::Branching;

//    for (auto &fname: SUPPORT_TEST_MODELS) test_supports(fname, supportcfg);
//}

bool is_outside_support_cone(const Vec3f &supp, const Vec3f &pt, float angle)
{
    Vec3d D = (pt - supp).cast<double>();
    double dot_sq = -D.z() * std::abs(-D.z());

    return dot_sq <
           D.squaredNorm() * std::cos(angle) * std::abs(std::cos(angle));
}

TEST_CASE("BranchingSupports::MergePointFinder", "[SLASupportGeneration][Branching]") {
    SECTION("Identical points have the same merge point") {
        Vec3f a{0.f, 0.f, 0.f}, b = a;
        float slope = PI / 4.f;

        auto mergept = branchingtree::find_merge_pt(a, b, slope);

        REQUIRE(bool(mergept));
        REQUIRE((*mergept - b).norm() < EPSILON);
        REQUIRE((*mergept - a).norm() < EPSILON);
    }

    // ^ Z
    // | a *
    // |
    // | b * <= mergept
    SECTION("Points at different heights have the lower point as mergepoint") {
        Vec3f a{0.f, 0.f, 0.f}, b = {0.f, 0.f, -1.f};
        float slope = PI / 4.f;

        auto mergept = branchingtree::find_merge_pt(a, b, slope);

        REQUIRE(bool(mergept));
        REQUIRE((*mergept - b).squaredNorm() < EPSILON);
    }

    // -|---------> X
    //  a       b
    //  *       *
    //      * <= mergept
    SECTION("Points at different X have mergept in the middle at lower Z") {
        Vec3f a{0.f, 0.f, 0.f}, b = {1.f, 0.f, 0.f};
        float slope = PI / 4.f;

        auto mergept = branchingtree::find_merge_pt(a, b, slope);

        REQUIRE(bool(mergept));

           // Distance of mergept should be equal from both input points
        float D = std::abs((*mergept - b).squaredNorm() - (*mergept - a).squaredNorm());

        REQUIRE(D < EPSILON);
        REQUIRE(!is_outside_support_cone(a, *mergept, slope));
        REQUIRE(!is_outside_support_cone(b, *mergept, slope));
    }

    // -|---------> Y
    //  a       b
    //  *       *
    //      * <= mergept
    SECTION("Points at different Y have mergept in the middle at lower Z") {
        Vec3f a{0.f, 0.f, 0.f}, b = {0.f, 1.f, 0.f};
        float slope = PI / 4.f;

        auto mergept = branchingtree::find_merge_pt(a, b, slope);

        REQUIRE(bool(mergept));

        // Distance of mergept should be equal from both input points
        float D = std::abs((*mergept - b).squaredNorm() - (*mergept - a).squaredNorm());

        REQUIRE(D < EPSILON);
        REQUIRE(!is_outside_support_cone(a, *mergept, slope));
        REQUIRE(!is_outside_support_cone(b, *mergept, slope));
    }

    SECTION("Points separated by less than critical angle have the lower point as mergept") {
        Vec3f a{0.f, 0.f, 0.f}, b = {-0.5f, -0.5f, -1.f};
        float slope = PI / 4.f;

        auto mergept = branchingtree::find_merge_pt(a, b, slope);

        REQUIRE(bool(mergept));
        REQUIRE((*mergept - b).norm() < EPSILON);
    }
}

TEST_CASE("BranchingSupports::ElevatedSupportsDoNotPierceModel", "[SLASupportGeneration][Branching]") {

    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 10.;
    supportcfg.tree_type = sla::SupportTreeType::Branching;

    for (auto fname : SUPPORT_TEST_MODELS)
        test_support_model_collision(fname, supportcfg);
}

TEST_CASE("BranchingSupports::FloorSupportsDoNotPierceModel", "[SLASupportGeneration][Branching]") {

    sla::SupportTreeConfig supportcfg;
    supportcfg.object_elevation_mm = 0;
    supportcfg.tree_type = sla::SupportTreeType::Branching;

    for (auto fname : SUPPORT_TEST_MODELS)
        test_support_model_collision(fname, supportcfg);
}

TEST_CASE("InitializedRasterShouldBeNONEmpty", "[SLARasterOutput]") {
    // Default Prusa SL1 display parameters
    sla::Resolution res{2560, 1440};
    sla::PixelDim   pixdim{120. / res.width_px, 68. / res.height_px};
    
    sla::RasterGrayscaleAAGammaPower raster(res, pixdim, {}, 1.);
    REQUIRE(raster.resolution().width_px == res.width_px);
    REQUIRE(raster.resolution().height_px == res.height_px);
    REQUIRE(raster.pixel_dimensions().w_mm == Approx(pixdim.w_mm));
    REQUIRE(raster.pixel_dimensions().h_mm == Approx(pixdim.h_mm));
}

TEST_CASE("MirroringShouldBeCorrect", "[SLARasterOutput]") {
    sla::RasterBase::TMirroring mirrorings[] = {sla::RasterBase::NoMirror,
                                                sla::RasterBase::MirrorX,
                                                sla::RasterBase::MirrorY,
                                                sla::RasterBase::MirrorXY};

    sla::RasterBase::Orientation orientations[] =
        {sla::RasterBase::roLandscape, sla::RasterBase::roPortrait};
    
    for (auto orientation : orientations)
        for (auto &mirror : mirrorings)
            check_raster_transformations(orientation, mirror);
}


TEST_CASE("RasterizedPolygonAreaShouldMatch", "[SLARasterOutput]") {
    double disp_w = 120., disp_h = 68.;
    sla::Resolution res{2560, 1440};
    sla::PixelDim pixdim{disp_w / res.width_px, disp_h / res.height_px};
    
    double gamma = 1.;
    sla::RasterGrayscaleAAGammaPower raster(res, pixdim, {}, gamma);
    auto bb = BoundingBox({0, 0}, {scaled(disp_w), scaled(disp_h)});
    
    ExPolygon poly = square_with_hole(10.);
    poly.translate(bb.center().x(), bb.center().y());
    raster.draw(poly);
    
    double a = poly.area() / (scaled<double>(1.) * scaled(1.));
    double ra = raster_white_area(raster);
    double diff = std::abs(a - ra);
    
    REQUIRE(diff <= predict_error(poly, pixdim));
    
    raster.clear();
    poly = square_with_hole(60.);
    poly.translate(bb.center().x(), bb.center().y());
    raster.draw(poly);
    
    a = poly.area() / (scaled<double>(1.) * scaled(1.));
    ra = raster_white_area(raster);
    diff = std::abs(a - ra);
    
    REQUIRE(diff <= predict_error(poly, pixdim));
    
    sla::RasterGrayscaleAA raster0(res, pixdim, {}, [](double) { return 0.; });
    REQUIRE(raster_pxsum(raster0) == 0);
    
    raster0.draw(poly);
    ra = raster_white_area(raster);
    REQUIRE(raster_pxsum(raster0) == 0);
}


TEST_CASE("halfcone test", "[halfcone]") {
    sla::DiffBridge br{Vec3d{1., 1., 1.}, Vec3d{10., 10., 10.}, 0.25, 0.5};

    indexed_triangle_set m = sla::get_mesh(br, 45);

    its_merge_vertices(m);
    its_write_obj(m, "Halfcone.obj");
}

TEST_CASE("Test concurrency")
{
    std::vector<double> vals = grid(0., 100., 10.);

    double ref = std::accumulate(vals.begin(), vals.end(), 0.);

    double s = execution::accumulate(ex_tbb, vals.begin(), vals.end(), 0.);

    REQUIRE(s == Approx(ref));
}
