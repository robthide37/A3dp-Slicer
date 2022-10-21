#ifndef SLICECSGMESH_HPP
#define SLICECSGMESH_HPP

#include "CSGMesh.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/ClipperUtils.hpp"

namespace Slic3r { namespace csg {

template<class ItCSG>
std::vector<ExPolygons> slice_csgmesh_ex(
    const Range<ItCSG>          &csg,
    const std::vector<float>    &slicegrid,
    const MeshSlicingParamsEx   &params,
    const std::function<void()> &throw_on_cancel = [] {})
{
    std::vector<ExPolygons> ret(slicegrid.size());

    MeshSlicingParamsEx params_cpy = params;
    auto trafo = params.trafo;
    for (const auto &m : csg) {
        const indexed_triangle_set *its = csg::get_mesh(m);
        if (!its)
            continue;

        params_cpy.trafo = trafo * csg::get_transform(m).template cast<double>();
        std::vector<ExPolygons> slices = slice_mesh_ex(*its,
                                                       slicegrid, params_cpy,
                                                       throw_on_cancel);

        assert(slices.size() == slicegrid.size());

        for (size_t i = 0; i < slicegrid.size(); ++i) {
            switch(get_operation(m)) {
            case CSGType::Union:
                for (ExPolygon &expoly : slices[i])
                    ret[i].emplace_back(std::move(expoly));

                ret[i] = union_ex(ret[i]);
                break;
            case CSGType::Difference:
                ret[i] = diff_ex(ret[i], slices[i]);
                break;
            case CSGType::Intersection:
                ret[i] = intersection_ex(ret[i], slices[i]);
                break;
            }
        }

        for (ExPolygons &slice : ret) {
            auto it = std::remove_if(slice.begin(), slice.end(), [](const ExPolygon &p){
                return p.area() < double(SCALED_EPSILON) * double(SCALED_EPSILON);
            });

            // Hopefully, ExPolygons are moved, not copied to new positions
            // and that is cheap for expolygons
            slice.erase(it, slice.end());
        }
    }

    return ret;
}

}} // namespace Slic3r::csg

#endif // SLICECSGMESH_HPP
