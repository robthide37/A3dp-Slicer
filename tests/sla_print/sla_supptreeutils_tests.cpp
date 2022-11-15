#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include "libslic3r/Execution/ExecutionSeq.hpp"
#include "libslic3r/SLA/SupportTreeUtils.hpp"

TEST_CASE("Avoid disk below junction", "[suptreeutils]")
{
    // In this test there will be a disk mesh with some radius, centered at
    // (0, 0, 0) and above the disk, a junction from which the support pillar
    // should be routed. The algorithm needs to find an avoidance route.

    using namespace Slic3r;

    constexpr double FromRadius = .5;
    constexpr double EndRadius  = 1.;
    constexpr double CylRadius  = 4.;
    constexpr double CylHeight  = 1.;

    sla::SupportTreeConfig cfg;

    indexed_triangle_set disk = its_make_cylinder(CylRadius, CylHeight);

    // 2.5 * CyRadius height should be enough to be able to insert a bridge
    // with 45 degree tilt above the disk.
    sla::Junction j{Vec3d{0., 0., 2.5 * CylRadius}, FromRadius};

    sla::SupportableMesh sm{disk, sla::SupportPoints{}, cfg};

    sla::GroundConnection conn =
            sla::deepsearch_ground_connection(ex_seq, sm, j, EndRadius, sla::DOWN);

#ifndef NDEBUG

    sla::SupportTreeBuilder builder;

    if (!conn)
        builder.add_junction(j);

    sla::build_ground_connection(builder, sm, conn);

    its_merge(disk, builder.merged_mesh());

    its_write_stl_ascii("output_disk.stl", "disk", disk);
#endif

    REQUIRE(bool(conn));

    // The route should include the source and one avoidance junction.
    REQUIRE(conn.path.size() == 2);

    // Check if the radius increases with each node
    REQUIRE(conn.path.front().r < conn.path.back().r);
    REQUIRE(conn.path.back().r < conn.pillar_base->r_top);

    // The end radius and the pillar base's upper radius should match
    REQUIRE(conn.pillar_base->r_top == Approx(EndRadius));

    // Check if the avoidance junction is indeed outside of the disk barrier's
    // edge.
    auto p = conn.path.back().pos;
    double pR = std::sqrt(p.x() * p.x()) + std::sqrt(p.y() * p.y());
    REQUIRE(pR + FromRadius > CylRadius);
}

TEST_CASE("Avoid disk below junction - Zero elevation", "[suptreeutils]")
{
    // In this test there will be a disk mesh with some radius, centered at
    // (0, 0, 0) and above the disk, a junction from which the support pillar
    // should be routed. The algorithm needs to find an avoidance route.

    using namespace Slic3r;

    constexpr double FromRadius = .5;
    constexpr double EndRadius  = 1.;
    constexpr double CylRadius  = 4.;
    constexpr double CylHeight  = 1.;

    sla::SupportTreeConfig cfg;
    cfg.object_elevation_mm = 0.;

    indexed_triangle_set disk = its_make_cylinder(CylRadius, CylHeight);

    // 2.5 * CyRadius height should be enough to be able to insert a bridge
    // with 45 degree tilt above the disk.
    sla::Junction j{Vec3d{0., 0., 2.5 * CylRadius}, FromRadius};

    sla::SupportableMesh sm{disk, sla::SupportPoints{}, cfg};

    sla::GroundConnection conn =
            sla::deepsearch_ground_connection(ex_seq, sm, j, EndRadius, sla::DOWN);

#ifndef NDEBUG

    sla::SupportTreeBuilder builder;

    if (!conn)
        builder.add_junction(j);

    sla::build_ground_connection(builder, sm, conn);

    its_merge(disk, builder.merged_mesh());

    its_write_stl_ascii("output_disk_ze.stl", "disk", disk);
#endif

    REQUIRE(bool(conn));

    // The route should include the source and one avoidance junction.
    REQUIRE(conn.path.size() == 2);

    // Check if the radius increases with each node
    REQUIRE(conn.path.front().r < conn.path.back().r);
    REQUIRE(conn.path.back().r < conn.pillar_base->r_top);

    // The end radius and the pillar base's upper radius should match
    REQUIRE(conn.pillar_base->r_top == Approx(EndRadius));

    // Check if the avoidance junction is indeed outside of the disk barrier's
    // edge.
    auto p = conn.path.back().pos;
    double pR = std::sqrt(p.x() * p.x()) + std::sqrt(p.y() * p.y());
    REQUIRE(pR + FromRadius > CylRadius);
}

TEST_CASE("Avoid disk below junction with barrier on the side", "[suptreeutils]")
{
    // In this test there will be a disk mesh with some radius, centered at
    // (0, 0, 0) and above the disk, a junction from which the support pillar
    // should be routed. The algorithm needs to find an avoidance route.

    using namespace Slic3r;

    constexpr double FromRadius = .5;
    constexpr double EndRadius  = 1.;
    constexpr double CylRadius  = 4.;
    constexpr double CylHeight  = 1.;
    constexpr double JElevX     = 2.5;

    sla::SupportTreeConfig cfg;

    indexed_triangle_set disk = its_make_cylinder(CylRadius, CylHeight);
    indexed_triangle_set wall = its_make_cube(1., 2 * CylRadius, JElevX * CylRadius);
    its_translate(wall, Vec3f{float(FromRadius), -float(CylRadius), 0.f});
    its_merge(disk, wall);

    // 2.5 * CyRadius height should be enough to be able to insert a bridge
    // with 45 degree tilt above the disk.
    sla::Junction j{Vec3d{0., 0., JElevX * CylRadius}, FromRadius};

    sla::SupportableMesh sm{disk, sla::SupportPoints{}, cfg};

    sla::GroundConnection conn =
            sla::deepsearch_ground_connection(ex_seq, sm, j, EndRadius, sla::DOWN);

#ifndef NDEBUG

    sla::SupportTreeBuilder builder;

    if (!conn)
        builder.add_junction(j);

    sla::build_ground_connection(builder, sm, conn);

    its_merge(disk, builder.merged_mesh());

    its_write_stl_ascii("output_disk_wall.stl", "disk_wall", disk);
#endif

    REQUIRE(bool(conn));

    // The route should include the source and one avoidance junction.
    REQUIRE(conn.path.size() == 2);

    // Check if the radius increases with each node
    REQUIRE(conn.path.front().r < conn.path.back().r);
    REQUIRE(conn.path.back().r < conn.pillar_base->r_top);

    // The end radius end the pillar base's upper radius should match
    REQUIRE(conn.pillar_base->r_top == Approx(EndRadius));

    // Check if the avoidance junction is indeed outside of the disk barrier's
    // edge.
    auto p = conn.path.back().pos;
    double pR = std::sqrt(p.x() * p.x()) + std::sqrt(p.y() * p.y());
    REQUIRE(pR + FromRadius > CylRadius);
}
