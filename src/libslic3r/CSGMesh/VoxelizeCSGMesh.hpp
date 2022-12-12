#ifndef VOXELIZECSGMESH_HPP
#define VOXELIZECSGMESH_HPP

#include <functional>

#include "CSGMesh.hpp"
#include "libslic3r/OpenVDBUtils.hpp"

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

    for (auto &csgpart : csgrange) {
        if (params.statusfn() && params.statusfn()(-1))
            break;

        VoxelGridPtr partgrid = get_voxelgrid(csgpart, params);

        if (!ret && get_operation(csgpart) == CSGType::Union) {
            ret = std::move(partgrid);
        } else if (ret) {
            switch (get_operation(csgpart)) {
            case CSGType::Union:
                grid_union(*ret, *partgrid);
                break;
            case CSGType::Difference:
                grid_difference(*ret, *partgrid);
                break;
            case CSGType::Intersection:
                grid_intersection(*ret, *partgrid);
                break;
            }
        }
    }

    return ret;
}

}} // namespace Slic3r::csg

#endif // VOXELIZECSGMESH_HPP
