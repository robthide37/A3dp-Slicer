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

struct Perimeter {
    size_t start_index;
    size_t end_index;
    size_t seam_index;

    bool aligned = false;
    Vec3f final_seam_position;
};

struct SeamCandidate {
    SeamCandidate(const Vec3f &pos, std::shared_ptr<Perimeter> perimeter,
            float local_ccw_angle,
            EnforcedBlockedSeamPoint type) :
            position(pos), perimeter(perimeter), visibility(0.0f), overhang(0.0f), higher_layer_overhang(
                    0.0f), local_ccw_angle(
                    local_ccw_angle), type(type) {
    }
    const Vec3f position;
    const std::shared_ptr<Perimeter> perimeter;
    float visibility;
    float overhang;
    float higher_layer_overhang; // represents how much is the position covered by the upper layer, useful for local visibility
    float local_ccw_angle;
    EnforcedBlockedSeamPoint type;
};

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
    static constexpr float expected_hits_per_area = 400.0f;
    static constexpr float considered_area_radius = 1.6f;

    static constexpr float cosine_hemisphere_sampling_power = 8.0f;

    static constexpr float polygon_local_angles_arm_distance = 0.6f;

    static constexpr float enforcer_blocker_sqr_distance_tolerance = 0.2f;

    static constexpr float seam_align_tolerable_dist = 2.0f;
    static constexpr size_t seam_align_tolerable_skips = 10;
    static constexpr size_t seam_align_minimum_string_seams = 4;

    //perimeter points per object per layer idx, and their corresponding KD trees
    std::unordered_map<const PrintObject*, std::vector<std::vector<SeamPlacerImpl::SeamCandidate>>> m_perimeter_points_per_object;
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
