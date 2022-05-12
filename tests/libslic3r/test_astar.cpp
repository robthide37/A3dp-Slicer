#include <catch2/catch.hpp>

#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/AStar.hpp"
#include "libslic3r/Execution/ExecutionSeq.hpp"
#include "libslic3r/PointGrid.hpp"

using namespace Slic3r;

struct PointGridTracer {
    using Node = size_t;
    const PointGrid<float> &grid;
    size_t final;

    PointGridTracer(const PointGrid<float> &g, size_t goal) :
        grid{g}, final{goal} {}

    template<class Fn>
    void foreach_reachable(size_t from, Fn &&fn) const
    {
        Vec3i from_crd = grid.get_coord(from);
        REQUIRE(grid.get_idx(from_crd) == from);

        if (size_t i = grid.get_idx(from_crd + Vec3i{ 1,  0,  0}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0,  1,  0}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0,  0,  1}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 1,  1,  0}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0,  1,  1}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 1,  1,  1}); i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{-1,  0,  0}); from_crd.x() > 0 && i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0, -1,  0}); from_crd.y() > 0 && i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0,  0, -1}); from_crd.z() > 0 && i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{-1, -1,  0}); from_crd.x() > 0 && from_crd.y() > 0 && i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{ 0, -1, -1}); from_crd.y() > 0 && from_crd.z() && i < grid.point_count()) fn(i);
        if (size_t i = grid.get_idx(from_crd + Vec3i{-1, -1, -1}); from_crd.x() > 0 && from_crd.y() > 0 && from_crd.z() && i < grid.point_count()) fn(i);

    }

    float distance(size_t a, size_t b) const
    {
        return (grid.get(a) - grid.get(b)).squaredNorm();
    }

    float goal_heuristic(size_t n) const
    {
        return n == final ? -1.f : (grid.get(n) - grid.get(final)).squaredNorm();
    }

    size_t unique_id(size_t n) const { return n; }
};

TEST_CASE("astar algorithm test over 3D point grid", "[AStar]") {
    auto vol = BoundingBox3Base<Vec3f>{{0.f, 0.f, 0.f}, {1.f, 1.f, 1.f}};

    auto pgrid = point_grid(ex_seq, vol, {0.1f, 0.1f, 0.1f});

    size_t target = pgrid.point_count() - 1;

    std::cout << "Tracing route to " << pgrid.get_coord(target).transpose() << std::endl;
    PointGridTracer pgt{pgrid, pgrid.point_count() - 1};
    std::vector<size_t> out;
    bool found = astar::search_route(pgt, size_t(0), std::back_inserter(out));

    std::cout << "Route taken: ";
    for (size_t i : out) {
        std::cout << "(" << pgrid.get_coord(i).transpose() << ") ";
    }
    std::cout << std::endl;

    REQUIRE(found);
}
