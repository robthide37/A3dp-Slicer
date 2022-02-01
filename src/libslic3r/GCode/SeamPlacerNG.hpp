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

struct SeamCandidate {
    explicit SeamCandidate(const Vec3d &pos) :
            position(pos), visibility(0.0), polygon_start(false) {
    }
    SeamCandidate(const Vec3d &pos, bool polygon_start) :
            position(pos), visibility(0.0), polygon_start(polygon_start) {
    }
    Vec3d position;
    float visibility;
    bool polygon_start;
};

struct HitInfo {
    Vec3d m_position;
    Vec3d m_surface_normal;
};


struct KDTreeCoordinateFunctor {
    KDTreeCoordinateFunctor(std::vector<SeamCandidate> *seam_candidates) :
            seam_candidates(seam_candidates) {
    }
    std::vector<SeamCandidate> *seam_candidates;
    float operator()(size_t index, size_t dim) const {
        return seam_candidates->operator[](index).position[dim];
    }
};

struct HitInfoCoordinateFunctor {
    HitInfoCoordinateFunctor(std::vector<HitInfo> *hit_points) :
        m_hit_points(hit_points) {
    }
    std::vector<HitInfo> *m_hit_points;
    float operator()(size_t index, size_t dim) const {
        return m_hit_points->operator[](index).m_position[dim];
    }
};
} // namespace SeamPlacerImpl

class SeamPlacer {
    using PointTree =
    KDTreeIndirect<3, coordf_t, SeamPlacerImpl::KDTreeCoordinateFunctor>;
    const size_t ray_count_per_object = 100000;

public:
    std::unordered_map<const PrintObject*, PointTree> m_perimeter_points_trees_per_object;
    std::unordered_map<const PrintObject*, std::vector<SeamPlacerImpl::SeamCandidate>> m_perimeter_points_per_object;

    void init(const Print &print);

    void place_seam(const PrintObject *po, ExtrusionLoop &loop, coordf_t unscaled_z, const Point &last_pos,
            bool external_first,
            double nozzle_diameter, const EdgeGrid::Grid *lower_layer_edge_grid);
};

} // namespace Slic3r

#endif // libslic3r_SeamPlacerNG_hpp_
