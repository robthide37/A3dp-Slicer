#ifndef SLICECSGMESH_HPP
#define SLICECSGMESH_HPP

#include "CSGMesh.hpp"

#include <stack>

#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Execution/ExecutionTBB.hpp"

namespace Slic3r { namespace csg {

template<class ItCSG>
std::vector<ExPolygons> slice_csgmesh_ex(
    const Range<ItCSG>          &csgrange,
    const std::vector<float>    &slicegrid,
    const MeshSlicingParamsEx   &params,
    const std::function<void()> &throw_on_cancel = [] {})
{
    struct Frame { CSGType op; std::vector<ExPolygons> slices; };

    std::stack opstack{std::vector<Frame>{}};

    MeshSlicingParamsEx params_cpy = params;
    auto trafo = params.trafo;
    auto nonempty_indices = reserve_vector<size_t>(slicegrid.size());

    if (!csgrange.empty() && csg::get_stack_operation(*csgrange.begin()) != CSGStackOp::Push)
        opstack.push({CSGType::Union, std::vector<ExPolygons>(slicegrid.size())});

    for (const auto &csgpart : csgrange) {
        const indexed_triangle_set *its = csg::get_mesh(csgpart);

        auto op = get_operation(csgpart);

        if (get_stack_operation(csgpart) == CSGStackOp::Push) {
            opstack.push({op, std::vector<ExPolygons>(slicegrid.size())});
            op = CSGType::Union;
        }

        Frame *top = &opstack.top();

        if (its) {
            params_cpy.trafo = trafo * csg::get_transform(csgpart).template cast<double>();
            std::vector<ExPolygons> slices = slice_mesh_ex(*its,
                                                           slicegrid, params_cpy,
                                                           throw_on_cancel);

            assert(slices.size() == slicegrid.size());

            nonempty_indices.clear();
            for (size_t i = 0; i < slicegrid.size(); ++i) {
                if (op == CSGType::Intersection || !slices[i].empty())
                    nonempty_indices.emplace_back(i);
            }

            auto mergefn = [&csgpart, &slices, &top](size_t i){
                switch(get_operation(csgpart)) {
                case CSGType::Union:
                    for (ExPolygon &expoly : slices[i])
                        top->slices[i].emplace_back(std::move(expoly));

                    break;
                case CSGType::Difference:
                    top->slices[i] = diff_ex(top->slices[i], slices[i]);
                    break;
                case CSGType::Intersection:
                    top->slices[i] = intersection_ex(top->slices[i], slices[i]);
                    break;
                }
            };

            execution::for_each(ex_tbb,
                                nonempty_indices.begin(), nonempty_indices.end(),
                                mergefn,
                                execution::max_concurrency(ex_tbb));
        }

        if (get_stack_operation(csgpart) == CSGStackOp::Pop) {
            std::vector<ExPolygons> popslices = std::move(top->slices);
            auto popop = opstack.top().op;
            opstack.pop();
            std::vector<ExPolygons> &prev_slices = opstack.top().slices;

            nonempty_indices.clear();
            for (size_t i = 0; i < slicegrid.size(); ++i) {
                if (popop == CSGType::Intersection || !popslices[i].empty())
                    nonempty_indices.emplace_back(i);
            }

            auto mergefn2 = [&popslices, &prev_slices, popop](size_t i){
                switch(popop) {
                case CSGType::Union:
                    for (ExPolygon &expoly : popslices[i])
                        prev_slices[i].emplace_back(std::move(expoly));

                    break;
                case CSGType::Difference:
                    prev_slices[i] = diff_ex(prev_slices[i], popslices[i]);
                    break;
                case CSGType::Intersection:
                    prev_slices[i] = intersection_ex(prev_slices[i], popslices[i]);
                    break;
                }
            };

            execution::for_each(ex_tbb,
                                nonempty_indices.begin(), nonempty_indices.end(),
                                mergefn2,
                                execution::max_concurrency(ex_tbb));
        }
    }

    std::vector<ExPolygons> ret = std::move(opstack.top().slices);

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
