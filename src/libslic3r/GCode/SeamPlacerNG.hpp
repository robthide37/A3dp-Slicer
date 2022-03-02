#ifndef libslic3r_SeamPlacerNG_hpp_
#define libslic3r_SeamPlacerNG_hpp_

#include <optional>
#include <vector>
#include <memory>
#include <atomic>

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

struct GlobalModelInfo;

enum class EnforcedBlockedSeamPoint {
    Blocked = 0,
    Neutral = 1,
    Enforced = 2,
};

// struct representing single perimeter loop
struct Perimeter {
    size_t start_index;
    size_t end_index; //inclusive!
    size_t seam_index;

    // During alignment, a final position may be stored here. In that case, aligned is set to true.
    // Note that final seam position is not limited to points of the perimeter loop. In theory it can be any position
    bool aligned = false;
    Vec3f final_seam_position;
};

//Struct over which all processing of perimeters is done. For each perimeter point, its respective candidate is created,
// then all the needed attributes are computed and finally, for each perimeter one point is chosen as seam.
// This seam position can be than further aligned
struct SeamCandidate {
    SeamCandidate(const Vec3f &pos, std::shared_ptr<Perimeter> perimeter,
            float local_ccw_angle,
            EnforcedBlockedSeamPoint type) :
            position(pos), perimeter(perimeter), visibility(0.0f), overhang(0.0f), local_ccw_angle(
                    local_ccw_angle), type(type) {
    }
    const Vec3f position;
    // pointer to Perimter loop of this point. It is shared across all points of the loop
    const std::shared_ptr<Perimeter> perimeter;
    float visibility;
    float overhang;
    float local_ccw_angle;
    EnforcedBlockedSeamPoint type;
};

// struct to represent hits of the mesh during occulision raycasting.
struct HitInfo {
    Vec3f position;
    Vec3f surface_normal;
};

struct SeamCandidateCoordinateFunctor {
    SeamCandidateCoordinateFunctor(std::vector<SeamCandidate> *seam_candidates) :
            seam_candidates(seam_candidates) {
    }
    std::vector<SeamCandidate> *seam_candidates;
    float operator()(size_t index, size_t dim) const {
        return seam_candidates->operator[](index).position[dim];
    }
};

struct HitInfoCoordinateFunctor {
    HitInfoCoordinateFunctor(std::vector<HitInfo> *hit_points) :
            hit_points(hit_points) {
    }
    std::vector<HitInfo> *hit_points;
    float operator()(size_t index, size_t dim) const {
        return hit_points->operator[](index).position[dim];
    }
};
} // namespace SeamPlacerImpl

class SeamPlacer {
public:
    using SeamCandidatesTree =
    KDTreeIndirect<3, float, SeamPlacerImpl::SeamCandidateCoordinateFunctor>;
    // Rough estimates of hits of the mesh during raycasting per surface circle defined by considered_area_radius
    static constexpr float expected_hits_per_area = 400.0f;
    // area considered when computing number of rays and then gathering visiblity info from the hits
    static constexpr float considered_area_radius = 4.0f;
    // quadric error limit of quadric decimation function used on the mesh before raycasting
    static constexpr float raycasting_decimation_target_error = 0.1f;

    // cosine sampling power represents how prefered are forward directions when raycasting from given spot
    // in this case, forward direction means towards the center of the mesh
    static constexpr float cosine_hemisphere_sampling_power = 4.0f;

    // arm length used during angles computation
    static constexpr float polygon_local_angles_arm_distance = 1.0f;

    // If enforcer or blocker is closer to the seam candidate than this limit, the seam candidate is set to Blocer or Enforcer
    static constexpr float enforcer_blocker_sqr_distance_tolerance = 0.1f;

    // When searching for seam clusters for alignment:
    // seam_align_tolerable_dist - if seam is closer to the previous seam position projected to the current layer than this value,
    //it belongs automaticaly to the cluster
    static constexpr float seam_align_tolerable_dist = 1.0f;
    // if the seam of the current layer is too far away, and the closest seam candidate is not very good, layer is skipped.
    // this param limits the number of allowed skips
    static constexpr size_t seam_align_tolerable_skips = 4;
    // minimum number of seams needed in cluster to make alignemnt happen
    static constexpr size_t seam_align_minimum_string_seams = 4;

    //The following data structures hold all perimeter points for all PrintObject. The structure is as follows:
    // Map of PrintObjects (PO) -> vector of layers of PO -> vector of perimeter points of the given layer
    std::unordered_map<const PrintObject*, std::vector<std::vector<SeamPlacerImpl::SeamCandidate>>> m_perimeter_points_per_object;
    // Map of PrintObjects (PO) -> vector of layers of PO -> unique_ptr to KD tree of all points of the given layer
    std::unordered_map<const PrintObject*, std::vector<std::unique_ptr<SeamCandidatesTree>>> m_perimeter_points_trees_per_object;

    void init(const Print &print);

    void place_seam(const Layer *layer, ExtrusionLoop &loop, bool external_first);

private:
    void gather_seam_candidates(const PrintObject *po, const SeamPlacerImpl::GlobalModelInfo &global_model_info);
    void calculate_candidates_visibility(const PrintObject *po,
            const SeamPlacerImpl::GlobalModelInfo &global_model_info);
    void calculate_overhangs(const PrintObject *po);
    template<typename Comparator>
    void align_seam_points(const PrintObject *po, const Comparator &comparator);
    template<typename Comparator>
    bool find_next_seam_in_string(const PrintObject *po, Vec3f &last_point_pos,
            size_t layer_idx, const Comparator &comparator,
            std::vector<std::pair<size_t, size_t>> &seam_strings,
            std::vector<std::pair<size_t, size_t>> &outliers);
};

} // namespace Slic3r

#endif // libslic3r_SeamPlacerNG_hpp_
