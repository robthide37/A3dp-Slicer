#ifndef VOXELIZECSGMESH_HPP
#define VOXELIZECSGMESH_HPP

#include "CSGMesh.hpp"
#include "libslic3r/OpenVDBUtils.hpp"

namespace Slic3r { namespace csg {

class VoxelizeParams {
    float m_voxel_scale        = 1.f;
    float m_exterior_bandwidth = 3.0f;
    float m_interior_bandwidth = 3.0f;

public:
    float voxel_scale() const noexcept { return m_voxel_scale; }
    float exterior_bandwidth() const noexcept { return m_exterior_bandwidth; }
    float interior_bandwidth() const noexcept { return m_interior_bandwidth; }

    VoxelizeParams &voxel_scale(float val) noexcept
    {
        m_voxel_scale = val;
        return *this;
    }
    VoxelizeParams &exterior_bandwidth(float val) noexcept
    {
        m_exterior_bandwidth = val;
        return *this;
    }
    VoxelizeParams &interior_bandwidth(float val) noexcept
    {
        m_interior_bandwidth = val;
        return *this;
    }
};

// This method can be overriden when a specific CSGPart type supports caching
// of the voxel grid
template<class CSGPartT>
VoxelGridPtr get_voxelgrid(const CSGPartT &csgpart, const VoxelizeParams &params)
{
    const indexed_triangle_set *its = csg::get_mesh(csgpart);
    VoxelGridPtr ret;

    if (its)
        ret = mesh_to_grid(*csg::get_mesh(csgpart),
                           csg::get_transform(csgpart),
                           params.voxel_scale(),
                           params.exterior_bandwidth(),
                           params.interior_bandwidth());

    return ret;
}

template<class It>
VoxelGridPtr voxelize_csgmesh(const Range<It>      &csgrange,
                              const VoxelizeParams &params = {})
{
    VoxelGridPtr ret;

    for (auto &csgpart : csgrange) {
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
