#ifndef libslic3r_SeamPlacer_hpp_
#define libslic3r_SeamPlacer_hpp_

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
struct SeamComparator;

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
    float flow_width;

    // During alignment, a final position may be stored here. In that case, finalized is set to true.
    // Note that final seam position is not limited to points of the perimeter loop. In theory it can be any position
    // Random position also uses this flexibility to set final seam point position
    bool finalized = false;
    Vec3f final_seam_position;
};

//Struct over which all processing of perimeters is done. For each perimeter point, its respective candidate is created,
// then all the needed attributes are computed and finally, for each perimeter one point is chosen as seam.
// This seam position can be than further aligned
struct SeamCandidate {
    SeamCandidate(const Vec3f &pos, std::shared_ptr<Perimeter> perimeter,
            float local_ccw_angle,
            EnforcedBlockedSeamPoint type) :
            position(pos), perimeter(perimeter), visibility(0.0f), overhang(0.0f), embedded_distance(0.0f), local_ccw_angle(
                    local_ccw_angle), type(type), central_enforcer(false) {
    }
    const Vec3f position;
    // pointer to Perimter loop of this point. It is shared across all points of the loop
    const std::shared_ptr<Perimeter> perimeter;
    float visibility;
    float overhang;
    // distance inside the merged layer regions, for detecting perimter points which are hidden indside the print (e.g. multimaterial join)
    // Negative sign means inside the print, comes from EdgeGrid structure
    float embedded_distance;
    float local_ccw_angle;
    EnforcedBlockedSeamPoint type;
    bool central_enforcer; //marks this candidate as central point of enforced segment on the perimeter - important for alignment
};

struct FaceVisibilityInfo {
    float visibility;
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
} // namespace SeamPlacerImpl

class SeamPlacer {
public:
    using SeamCandidatesTree =
    KDTreeIndirect<3, float, SeamPlacerImpl::SeamCandidateCoordinateFunctor>;
    static constexpr size_t raycasting_decimation_target_triangle_count = 10000;
    static constexpr float raycasting_subdivision_target_length = 2.0f;
    //square of number of rays per triangle
    static constexpr size_t sqr_rays_per_triangle = 7;

    // arm length used during angles computation
    static constexpr float polygon_local_angles_arm_distance = 0.5f;

    // increases angle importance at the cost of deacreasing visibility info importance. must be > 0
    static constexpr float additional_angle_importance = 0.3f;

    // If enforcer or blocker is closer to the seam candidate than this limit, the seam candidate is set to Blocker or Enforcer
    static constexpr float enforcer_blocker_distance_tolerance = 0.3f;
    // For long polygon sides, if they are close to the custom seam drawings, they are oversampled with this step size
    static constexpr float enforcer_blocker_oversampling_distance = 0.1f;

    // When searching for seam clusters for alignment:
    // following value describes, how much worse score can point have and still be picked into seam cluster instead of original seam point on the same layer
    static constexpr float seam_align_score_tolerance = 0.5f;
    // seam_align_tolerable_dist - if next layer closes point is too far away, break string
    static constexpr float seam_align_tolerable_dist = 1.0f;
    // if the seam of the current layer is too far away, and the closest seam candidate is not very good, layer is skipped.
    // this param limits the number of allowed skips
    static constexpr size_t seam_align_tolerable_skips = 4;
    // minimum number of seams needed in cluster to make alignemnt happen
    static constexpr size_t seam_align_minimum_string_seams = 6;
    // points covered by spline; determines number of splines for the given string
    static constexpr size_t seam_align_seams_per_spline = 30;

    //The following data structures hold all perimeter points for all PrintObject. The structure is as follows:
    // Map of PrintObjects (PO) -> vector of layers of PO -> vector of perimeter points of the given layer
    std::unordered_map<const PrintObject*, std::vector<std::vector<SeamPlacerImpl::SeamCandidate>>> m_perimeter_points_per_object;
    // Map of PrintObjects (PO) -> vector of layers of PO -> unique_ptr to KD tree of all points of the given layer
    std::unordered_map<const PrintObject*, std::vector<std::unique_ptr<SeamCandidatesTree>>> m_perimeter_points_trees_per_object;

    void init(const Print &print);

    void place_seam(const Layer *layer, ExtrusionLoop &loop, bool external_first, const Point &last_pos) const;

private:
    void gather_seam_candidates(const PrintObject *po, const SeamPlacerImpl::GlobalModelInfo &global_model_info,
            const SeamPosition configured_seam_preference);
    void calculate_candidates_visibility(const PrintObject *po,
            const SeamPlacerImpl::GlobalModelInfo &global_model_info);
    void calculate_overhangs_and_layer_embedding(const PrintObject *po);
    void align_seam_points(const PrintObject *po, const SeamPlacerImpl::SeamComparator &comparator);
    bool find_next_seam_in_layer(const PrintObject *po,
            std::pair<size_t, size_t> &last_point_indexes,
            size_t layer_idx, const SeamPlacerImpl::SeamComparator &comparator,
            std::vector<std::pair<size_t, size_t>> &seam_string);
};

} // namespace Slic3r

#endif // libslic3r_SeamPlacer_hpp_
