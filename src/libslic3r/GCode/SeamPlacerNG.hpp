#ifndef libslic3r_SeamPlacerNG_hpp_
#define libslic3r_SeamPlacerNG_hpp_

#include <optional>
#include <vector>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/KDTreeIndirect.hpp"

namespace Slic3r {

class PrintObject;
class ExtrusionLoop;
class Print;
class Layer;

namespace EdgeGrid {
class Grid;
}

namespace SeamPlacerImpl {

struct SeamCandidate
{
    SeamCandidate(const Vec3d &pos) : position(pos), visibility(0) {}
    // TODO is there some equivalent for Vec with coordf_t type?
    Vec3d  position;
    size_t visibility;
};

struct KDTreeCoordinateFunctor
{
    KDTreeCoordinateFunctor(std::vector<SeamCandidate> *seam_candidates)
        : seam_candidates(seam_candidates)
    {}
    std::vector<SeamCandidate> *seam_candidates;
    float                       operator()(size_t index, size_t dim) const
    {
        return seam_candidates->operator[](index).position[dim];
    }
};
} // namespace SeamPlacerImpl

class SeamPlacer
{
public:
    std::vector<
        KDTreeIndirect<3, coordf_t, SeamPlacerImpl::KDTreeCoordinateFunctor>>
        seam_candidates_trees;

    void init(const Print &print);

    void plan_perimeters(const std::vector<const ExtrusionEntity *> perimeters,
                         const Layer                               &layer,
                         SeamPosition          seam_position,
                         Point                 last_pos,
                         coordf_t              nozzle_dmr,
                         const PrintObject    *po,
                         const EdgeGrid::Grid *lower_layer_edge_grid)
    {}

    void place_seam(ExtrusionLoop        &loop,
                    const Point          &last_pos,
                    bool                  external_first,
                    double                nozzle_diameter,
                    const EdgeGrid::Grid *lower_layer_edge_grid)
    {
        Point seam = last_pos;
        if (!loop.split_at_vertex(seam))
            // The point is not in the original loop.
            // Insert it.
            loop.split_at(seam, true);
    }
};

} // namespace Slic3r

#endif // libslic3r_SeamPlacerNG_hpp_
