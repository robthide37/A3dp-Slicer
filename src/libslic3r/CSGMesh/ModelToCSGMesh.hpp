#ifndef MODELTOCSGMESH_HPP
#define MODELTOCSGMESH_HPP

#include "CSGMesh.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/SLA/Hollowing.hpp"
#include "libslic3r/MeshSplitImpl.hpp"

namespace Slic3r { namespace csg {

enum ModelParts {
    mpartsPositive = 1,
    mpartsNegative = 2,
    mpartsDrillHoles = 4,
    mpartsDoSplits = 8
};

template<class OutIt>
void model_to_csgmesh(const ModelObject &mo,
                      const Transform3d &trafo,
                      OutIt              out,
                      // values of ModelParts ORed
                      int                parts_to_include = mpartsPositive
                      )
{
    bool do_positives  = parts_to_include & mpartsPositive;
    bool do_negatives  = parts_to_include & mpartsNegative;
    bool do_drillholes = parts_to_include & mpartsDrillHoles;
    bool do_splits     = parts_to_include & mpartsDoSplits;

    for (const ModelVolume *vol : mo.volumes) {
        if (vol && vol->mesh_ptr() &&
            ((do_positives && vol->is_model_part()) ||
             (do_negatives && vol->is_negative_volume()))) {

            if (do_splits && its_is_splittable(vol->mesh().its)) {
                CSGPart part_begin{{}, vol->is_model_part() ? CSGType::Union : CSGType::Difference};
                part_begin.stack_operation = CSGStackOp::Push;
                *out = std::move(part_begin);
                ++out;

                its_split(vol->mesh().its, SplitOutputFn{[&out, &vol, &trafo](indexed_triangle_set &&its) {
                              CSGPart part{std::make_unique<indexed_triangle_set>(std::move(its)),
                                       CSGType::Union,
                                       (trafo * vol->get_matrix()).cast<float>()};

                              *out = std::move(part);
                              ++out;
                          }});

                CSGPart part_end{{}};
                part_end.stack_operation = CSGStackOp::Pop;
                *out = std::move(part_end);
                ++out;
            } else {
                CSGPart part{&(vol->mesh().its),
                             vol->is_model_part() ? CSGType::Union : CSGType::Difference,
                             (trafo * vol->get_matrix()).cast<float>()};

                *out = std::move(part);
                ++out;
            }
        }
    }

    if (do_drillholes) {
        sla::DrainHoles drainholes = sla::transformed_drainhole_points(mo, trafo);

        for (const sla::DrainHole &dhole : drainholes) {
            CSGPart part{std::make_unique<const indexed_triangle_set>(
                             dhole.to_mesh()),
                         CSGType::Difference};

            *out = std::move(part);
            ++out;
        }
    }
}

}} // namespace Slic3r::csg

#endif // MODELTOCSGMESH_HPP
