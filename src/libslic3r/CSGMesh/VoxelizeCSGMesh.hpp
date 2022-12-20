#ifndef VOXELIZECSGMESH_HPP
#define VOXELIZECSGMESH_HPP

#include <functional>
#include <stack>

#include "CSGMesh.hpp"
#include "libslic3r/OpenVDBUtils.hpp"
#include "libslic3r/Execution/ExecutionTBB.hpp"

namespace Slic3r { namespace csg {

using VoxelizeParams = MeshToGridParams;

// This method can be overriden when a specific CSGPart type supports caching
// of the voxel grid
template<class CSGPartT>
VoxelGridPtr get_voxelgrid(const CSGPartT &csgpart, VoxelizeParams params)
{
    const indexed_triangle_set *its = csg::get_mesh(csgpart);
    VoxelGridPtr ret;

    params.trafo(params.trafo() * csg::get_transform(csgpart));

    if (its)
        ret = mesh_to_grid(*its, params);

    return ret;
}

template<class It>
VoxelGridPtr voxelize_csgmesh(const Range<It>      &csgrange,
                              const VoxelizeParams &params = {})
{
    VoxelGridPtr ret;

    std::vector<VoxelGridPtr> grids (csgrange.size());

    execution::for_each(ex_tbb, size_t(0), csgrange.size(), [&](size_t csgidx) {
        if (params.statusfn() && params.statusfn()(-1))
            return;

        auto it = csgrange.begin();
        std::advance(it, csgidx);
        auto &csgpart = *it;
        grids[csgidx] = get_voxelgrid(csgpart, params);
    }, execution::max_concurrency(ex_tbb));

    size_t csgidx = 0;
    struct Frame { CSGType op = CSGType::Union; VoxelGridPtr grid; };
    std::stack opstack{std::vector<Frame>{}};

    if (!csgrange.empty() && csg::get_stack_operation(*csgrange.begin()) != CSGStackOp::Push)
        opstack.push({});

    for (auto &csgpart : csgrange) {
        if (params.statusfn() && params.statusfn()(-1))
            break;

        auto &partgrid = grids[csgidx++];

        auto op = get_operation(csgpart);

        if (get_stack_operation(csgpart) == CSGStackOp::Push) {
            opstack.push({op, nullptr});
            op = CSGType::Union;
        }

        Frame *top = &opstack.top();

        if (!top->grid && op == CSGType::Union) {
            top->grid = std::move(partgrid);
        } else if (top->grid && partgrid) {
            switch (get_operation(csgpart)) {
            case CSGType::Union:
                grid_union(*(top->grid), *partgrid);
                break;
            case CSGType::Difference:
                grid_difference(*(top->grid), *partgrid);
                break;
            case CSGType::Intersection:
                grid_intersection(*(top->grid), *partgrid);
                break;
            }
        }

        if (get_stack_operation(csgpart) == CSGStackOp::Pop) {
            VoxelGridPtr popgrid = std::move(top->grid);
            auto popop = opstack.top().op;
            opstack.pop();
            VoxelGridPtr &grid = opstack.top().grid;

            switch (popop) {
            case CSGType::Union:
                grid_union(*grid, *popgrid);
                break;
            case CSGType::Difference:
                grid_difference(*grid, *popgrid);
                break;
            case CSGType::Intersection:
                grid_intersection(*grid, *popgrid);
                break;
            }
        }
    }

    ret = std::move(opstack.top().grid);

    return ret;
}

}} // namespace Slic3r::csg

#endif // VOXELIZECSGMESH_HPP
