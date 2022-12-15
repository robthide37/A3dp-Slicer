#ifndef SLICECSGMESH_HPP
#define SLICECSGMESH_HPP

#include "CSGMesh.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Execution/ExecutionTBB.hpp"

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
    auto nonempty_indices = reserve_vector<size_t>(slicegrid.size());

    for (const auto &m : csg) {
        const indexed_triangle_set *its = csg::get_mesh(m);
        if (!its)
            continue;

        params_cpy.trafo = trafo * csg::get_transform(m).template cast<double>();
        std::vector<ExPolygons> slices = slice_mesh_ex(*its,
                                                       slicegrid, params_cpy,
                                                       throw_on_cancel);

        assert(slices.size() == slicegrid.size());

        nonempty_indices.clear();
        for (size_t i = 0; i < slicegrid.size(); ++i) {
            if (get_operation(m) == CSGType::Intersection || !slices[i].empty())
                nonempty_indices.emplace_back(i);
        }

        auto mergefn = [&m, &slices, &ret](size_t i){
            switch(get_operation(m)) {
            case CSGType::Union:
                for (ExPolygon &expoly : slices[i])
                    ret[i].emplace_back(std::move(expoly));

                break;
            case CSGType::Difference:
                ret[i] = diff_ex(ret[i], slices[i]);
                break;
            case CSGType::Intersection:
                ret[i] = intersection_ex(ret[i], slices[i]);
                break;
            }
        };

        execution::for_each(ex_tbb,
                            nonempty_indices.begin(), nonempty_indices.end(),
                            mergefn,
                            execution::max_concurrency(ex_tbb));
    }

    execution::for_each(ex_tbb, ret.begin(), ret.end(), [](ExPolygons &slice) {
        auto it = std::remove_if(slice.begin(), slice.end(), [](const ExPolygon &p){
            return p.area() < double(SCALED_EPSILON) * double(SCALED_EPSILON);
        });

        // Hopefully, ExPolygons are moved, not copied to new positions
        // and that is cheap for expolygons
        slice.erase(it, slice.end());
        slice = union_ex(slice);
    }, execution::max_concurrency(ex_tbb));

    return ret;
}

}} // namespace Slic3r::csg

#endif // SLICECSGMESH_HPP
