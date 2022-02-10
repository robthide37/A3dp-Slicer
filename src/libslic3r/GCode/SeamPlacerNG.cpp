#include "SeamPlacerNG.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <random>
#include <algorithm>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/SVG.hpp"
#include "libslic3r/Layer.hpp"

// TODO remove
#include <boost/nowide/cstdio.hpp>

namespace Slic3r {

namespace SeamPlacerImpl {

/// Coordinate frame
class Frame {
public:
    Frame() {
        mX = Vec3d(1, 0, 0);
        mY = Vec3d(0, 1, 0);
        mZ = Vec3d(0, 0, 1);
    }

    Frame(const Vec3d &x, const Vec3d &y, const Vec3d &z) :
            mX(x), mY(y), mZ(z) {
    }

    void set_from_z(const Vec3d &z) {
        mZ = z.normalized();
        Vec3d tmpZ = mZ;
        Vec3d tmpX = (std::abs(tmpZ.x()) > 0.99f) ? Vec3d(0, 1, 0) : Vec3d(1, 0, 0);
        mY = (tmpZ.cross(tmpX)).normalized();
        mX = mY.cross(tmpZ);
    }

    Vec3d to_world(const Vec3d &a) const {
        return a.x() * mX + a.y() * mY + a.z() * mZ;
    }

    Vec3d to_local(const Vec3d &a) const {
        return Vec3d(mX.dot(a), mY.dot(a), mZ.dot(a));
    }

    const Vec3d& binormal() const {
        return mX;
    }

    const Vec3d& tangent() const {
        return mY;
    }

    const Vec3d& normal() const {
        return mZ;
    }

private:
    Vec3d mX, mY, mZ;
};

Vec3d sample_sphere_uniform(const Vec2f &samples) {
    float term1 = 2.0f * M_PIf32 * samples.x();
    float term2 = 2.0f * sqrt(samples.y() - samples.y() * samples.y());
    return {cos(term1) * term2, sin(term1) * term2,
        1.0f - 2.0f * samples.y()};
}

Vec3d sample_power_cosine_hemisphere(const Vec2f &samples, float power) {
    float term1 = 2.f * M_PIf32 * samples.x();
    float term2 = pow(samples.y(), 1.f / (power + 1.f));
    float term3 = sqrt(1.f - term2 * term2);

    return Vec3d(cos(term1) * term3, sin(term1) * term3, term2);
}

std::vector<HitInfo> raycast_visibility(size_t ray_count,
        const AABBTreeIndirect::Tree<3, float> &raycasting_tree,
        const indexed_triangle_set &triangles) {
    auto bbox = raycasting_tree.node(0).bbox;
    Vec3d vision_sphere_center = bbox.center().cast<double>();
    float vision_sphere_raidus = (bbox.sizes().maxCoeff() * 0.55); // 0.5 (half) covers whole object,
                                                                   // 0.05 added to avoid corner cases

    // Prepare random samples per ray
//    std::random_device rnd_device;
    std::mt19937 mersenne_engine { 12345 };
    std::uniform_real_distribution<float> dist { 0, 1 };

    auto gen = [&dist, &mersenne_engine]() {
        return Vec2f(dist(mersenne_engine), dist(mersenne_engine));
    };

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: generate random samples: start";
    std::vector<Vec2f> global_dir_random_samples(ray_count);
    generate(begin(global_dir_random_samples), end(global_dir_random_samples), gen);
    std::vector<Vec2f> local_dir_random_samples(ray_count);
    generate(begin(local_dir_random_samples), end(local_dir_random_samples), gen);

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: generate random samples: end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: raycast visibility for " << ray_count << " rays: start";
    // raycast visibility
    std::vector<HitInfo> hit_points = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, ray_count),
            std::vector<HitInfo> { },
            [&](tbb::blocked_range<size_t> r, std::vector<HitInfo> init) {
                for (size_t index = r.begin(); index < r.end(); ++index) {
                    Vec3d global_ray_dir = sample_sphere_uniform(global_dir_random_samples[index]);
                    Vec3d ray_origin = (vision_sphere_center - global_ray_dir * vision_sphere_raidus);
                    Vec3d local_dir = sample_power_cosine_hemisphere(local_dir_random_samples[index], SeamPlacer::cosine_hemisphere_sampling_power);

                    Frame f;
                    f.set_from_z(global_ray_dir);
                    Vec3d final_ray_dir = (f.to_world(local_dir));

                    igl::Hit hitpoint;
                    // FIXME: This AABBTTreeIndirect query will not compile for float ray origin and
                    // direction for some reason
                    auto hit = AABBTreeIndirect::intersect_ray_first_hit(triangles.vertices,
                            triangles.indices, raycasting_tree, ray_origin, final_ray_dir, hitpoint);

                    if (hit) {
                        auto face = triangles.indices[hitpoint.id];
                        auto edge1 = triangles.vertices[face[1]] - triangles.vertices[face[0]];
                        auto edge2 = triangles.vertices[face[2]] - triangles.vertices[face[0]];

                        Vec3d hit_pos = (triangles.vertices[face[0]] + edge1 * hitpoint.u + edge2 * hitpoint.v).cast<
                                double>();
                        Vec3d surface_normal = edge1.cross(edge2).cast<double>().normalized();

                        init.push_back(HitInfo { hit_pos, surface_normal });
                    }
                }
                return init;
            },
            [](std::vector<HitInfo> left, const std::vector<HitInfo>& right) {
                left.insert(left.end(), right.begin(), right.end());
                return left;
            }
            );

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: raycast visibility for " << ray_count << " rays: end";

    return hit_points;
}

std::vector<float> calculate_polygon_angles_at_vertices(const Polygon &polygon, const std::vector<float> &lengths,
        float min_arm_length)
        {
    assert(polygon.points.size() + 1 == lengths.size());
    if (min_arm_length > 0.25f * lengths.back())
        min_arm_length = 0.25f * lengths.back();

    // Find the initial prev / next point span.
    size_t idx_prev = polygon.points.size();
    size_t idx_curr = 0;
    size_t idx_next = 1;
    while (idx_prev > idx_curr && lengths.back() - lengths[idx_prev] < min_arm_length)
        --idx_prev;
    while (idx_next < idx_prev && lengths[idx_next] < min_arm_length)
        ++idx_next;

    std::vector<float> angles(polygon.points.size(), 0.f);
    for (; idx_curr < polygon.points.size(); ++idx_curr) {
        // Move idx_prev up until the distance between idx_prev and idx_curr is lower than min_arm_length.
        if (idx_prev >= idx_curr) {
            while (idx_prev < polygon.points.size()
                    && lengths.back() - lengths[idx_prev] + lengths[idx_curr] > min_arm_length)
                ++idx_prev;
            if (idx_prev == polygon.points.size())
                idx_prev = 0;
        }
        while (idx_prev < idx_curr && lengths[idx_curr] - lengths[idx_prev] > min_arm_length)
            ++idx_prev;
        // Move idx_prev one step back.
        if (idx_prev == 0)
            idx_prev = polygon.points.size() - 1;
        else
            --idx_prev;
        // Move idx_next up until the distance between idx_curr and idx_next is greater than min_arm_length.
        if (idx_curr <= idx_next) {
            while (idx_next < polygon.points.size() && lengths[idx_next] - lengths[idx_curr] < min_arm_length)
                ++idx_next;
            if (idx_next == polygon.points.size())
                idx_next = 0;
        }
        while (idx_next < idx_curr && lengths.back() - lengths[idx_curr] + lengths[idx_next] < min_arm_length)
            ++idx_next;
        // Calculate angle between idx_prev, idx_curr, idx_next.
        const Point &p0 = polygon.points[idx_prev];
        const Point &p1 = polygon.points[idx_curr];
        const Point &p2 = polygon.points[idx_next];
        const Point v1 = p1 - p0;
        const Point v2 = p2 - p1;
        int64_t dot = int64_t(v1(0)) * int64_t(v2(0)) + int64_t(v1(1)) * int64_t(v2(1));
        int64_t cross = int64_t(v1(0)) * int64_t(v2(1)) - int64_t(v1(1)) * int64_t(v2(0));
        float angle = float(atan2(double(cross), double(dot)));
        angles[idx_curr] = angle;
    }

    return angles;
}

struct GlobalModelInfo {
    std::vector<HitInfo> geometry_raycast_hits;
    KDTreeIndirect<3, coordf_t, HitInfoCoordinateFunctor> raycast_hits_tree;
    indexed_triangle_set enforcers;
    indexed_triangle_set blockers;
    AABBTreeIndirect::Tree<3, float> enforcers_tree;
    AABBTreeIndirect::Tree<3, float> blockers_tree;

    GlobalModelInfo() :
            raycast_hits_tree(HitInfoCoordinateFunctor { &geometry_raycast_hits }) {
    }

    double enforcer_distance_check(const Vec3d &position) const {
        size_t hit_idx_out;
        Vec3d closest_vec3d;
        return AABBTreeIndirect::squared_distance_to_indexed_triangle_set(enforcers.vertices, enforcers.indices,
                enforcers_tree, position, hit_idx_out, closest_vec3d);
    }

    double blocker_distance_check(const Vec3d &position) const {
        size_t hit_idx_out;
        Vec3d closest_vec3d;
        return AABBTreeIndirect::squared_distance_to_indexed_triangle_set(blockers.vertices, blockers.indices,
                blockers_tree, position, hit_idx_out, closest_vec3d);
    }

    double calculate_point_visibility(const Vec3d &position, double max_distance) const {
        auto nearby_points = find_nearby_points(raycast_hits_tree, position, max_distance);
        double visibility = 0;
        for (const auto &hit_point_index : nearby_points) {
            double distance =
                    (position - geometry_raycast_hits[hit_point_index].m_position).norm();
            visibility += max_distance - distance; // The further away from the perimeter point,
            // the less representative ray hit is
        }
        return visibility;

    }
}
;

void process_perimeter_polygon(const Polygon &orig_polygon, coordf_t z_coord, std::vector<SeamCandidate> &result_vec,
        const GlobalModelInfo &global_model_info) {
    Polygon polygon = orig_polygon;
    polygon.make_counter_clockwise();
    std::vector<float> lengths = polygon.parameter_by_length();
    std::vector<float> angles = calculate_polygon_angles_at_vertices(polygon, lengths,
            SeamPlacer::polygon_angles_arm_distance);

    Vec3d last_enforcer_checked_point { 0, 0, -1 };
    double enforcer_dist_sqr = global_model_info.enforcer_distance_check(last_enforcer_checked_point);
    Vec3d last_blocker_checked_point { 0, 0, -1 };
    double blocker_dist_sqr = global_model_info.blocker_distance_check(last_blocker_checked_point);

    for (size_t index = 0; index < polygon.size(); ++index) {
        Vec2d unscaled_p = unscale(polygon[index]);
        Vec3d unscaled_position = Vec3d { unscaled_p.x(), unscaled_p.y(), z_coord };
        EnforcedBlockedSeamPoint type = EnforcedBlockedSeamPoint::NONE;

        float ccw_angle = angles[index];

        if (enforcer_dist_sqr >= 0) { // if enforcer dist < 0, it means there are no enforcers, skip
            //if there is enforcer, any other enforcer cannot be in a sphere defined by last check point and enforcer distance
            // so as long as we are at least enforcer_blocker_distance_tolerance deep in that area, and the enforcer distance is greater than
            // enforcer_blocker_distance_tolerance, we are fine.
            if (enforcer_dist_sqr < SeamPlacer::enforcer_blocker_sqr_distance_tolerance
                    ||
                    (last_enforcer_checked_point - unscaled_position).squaredNorm()
                            >= enforcer_dist_sqr - 2 * SeamPlacer::enforcer_blocker_sqr_distance_tolerance) {
                //do check
                enforcer_dist_sqr = global_model_info.enforcer_distance_check(unscaled_position);
                last_enforcer_checked_point = unscaled_position;
                if (enforcer_dist_sqr < SeamPlacer::enforcer_blocker_sqr_distance_tolerance) {
                    type = EnforcedBlockedSeamPoint::ENFORCED;
                }
            }
        }
        //same for blockers
        if (blocker_dist_sqr >= 0) {
            if (blocker_dist_sqr < SeamPlacer::enforcer_blocker_sqr_distance_tolerance
                    ||
                    (last_blocker_checked_point - unscaled_position).squaredNorm()
                            >= blocker_dist_sqr - 2 * SeamPlacer::enforcer_blocker_sqr_distance_tolerance) {
                blocker_dist_sqr = global_model_info.blocker_distance_check(unscaled_position);
                last_blocker_checked_point = unscaled_position;
                if (blocker_dist_sqr < SeamPlacer::enforcer_blocker_sqr_distance_tolerance) {
                    type = EnforcedBlockedSeamPoint::BLOCKED;
                }
            }
        }

        result_vec.emplace_back(unscaled_position, polygon.size() - index - 1, ccw_angle, type);
    }
}

std::pair<size_t, size_t> find_previous_and_next_perimeter_point(const std::vector<SeamCandidate> &perimeter_points,
        size_t index) {
    const SeamCandidate &current = perimeter_points[index];

    size_t prev = index > 0 ? index - 1 : index;
    size_t next = index + 1 < perimeter_points.size() ? index + 1 : index;

    //NOTE: dont forget that m_polygon_index_reverse are reversed indexes, so 0 is last point
    if (current.m_polygon_index_reverse == 0) {
        // next is at the start of loop
        //find start
        size_t start = index;
        while (start > 0 && perimeter_points[start - 1].m_polygon_index_reverse != 0) {
            start--;
        }
        next = start;
    }

    if (index > 1 && perimeter_points[index - 1].m_polygon_index_reverse == 0) {
        //prev is at the end of loop
        prev = index + perimeter_points[index].m_polygon_index_reverse;
    }

    return {prev,next};
}

float calculate_overhang(const SeamCandidate &point, const SeamCandidate &under_a, const SeamCandidate &under_b,
        const SeamCandidate &under_c) {
    auto p = Vec2d { point.m_position.x(), point.m_position.y() };
    auto a = Vec2d { under_a.m_position.x(), under_a.m_position.y() };
    auto b = Vec2d { under_b.m_position.x(), under_b.m_position.y() };
    auto c = Vec2d { under_c.m_position.x(), under_c.m_position.y() };

    auto oriented_line_dist = [](const Vec2d a, const Vec2d b, const Vec2d p) {
        return -((b.x() - a.x()) * (a.y() - p.y()) - (a.x() - p.x()) * (b.y() - a.y())) / (a - b).norm();
    };

    auto dist_ab = oriented_line_dist(a, b, p);
    auto dist_bc = oriented_line_dist(b, c, p);

    if (under_b.m_ccw_angle > 0 && dist_ab > 0 && dist_bc > 0) { //convex shape, p is inside
        return 0;
    }

    if (under_b.m_ccw_angle < 0 && (dist_ab < 0 || dist_bc < 0)) { //concave shape, p is inside
        return 0;
    }

    return Vec2d((p - b).norm(), std::min(abs(dist_ab), abs(dist_bc))).norm();

}

template<typename CompareFunc>
void pick_seam_point(std::vector<SeamCandidate> &perimeter_points, size_t start_index, size_t end_index,
        const CompareFunc &is_first_better) {
    size_t seam_index = start_index;
    for (size_t index = start_index + 1; index <= end_index; ++index) {
        if (is_first_better(perimeter_points[index], perimeter_points[seam_index])) {
            seam_index = index;
        }
    }

    for (size_t index = start_index; index <= end_index; ++index) {
        perimeter_points[index].m_seam_index = seam_index;
        perimeter_points[index].m_nearby_seam_points.get()->store(0, std::memory_order_relaxed);
    }
}

void gather_global_model_info(GlobalModelInfo &result, const PrintObject *po) {
    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: build AABB tree for raycasting and gather occlusion info: start";
// Build AABB tree for raycasting
    auto obj_transform = po->trafo_centered();
    auto triangle_set = po->model_object()->raw_indexed_triangle_set();
    its_transform(triangle_set, obj_transform);

    auto raycasting_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(triangle_set.vertices,
            triangle_set.indices);

    result.geometry_raycast_hits = raycast_visibility(SeamPlacer::ray_count_per_object, raycasting_tree,
            triangle_set);
    result.raycast_hits_tree.build(result.geometry_raycast_hits.size());

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: build AABB tree for raycasting and gather occlusion info: end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: build AABB trees for raycasting enforcers/blockers: start";

    for (const ModelVolume *mv : po->model_object()->volumes) {
        if (mv->is_model_part()) {
            its_merge(result.enforcers, mv->seam_facets.get_facets(*mv, EnforcerBlockerType::ENFORCER));
            its_merge(result.blockers, mv->seam_facets.get_facets(*mv, EnforcerBlockerType::BLOCKER));
        }
    }
    its_transform(result.enforcers, obj_transform);
    its_transform(result.blockers, obj_transform);

    result.enforcers_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(result.enforcers.vertices,
            result.enforcers.indices);
    result.blockers_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(result.blockers.vertices,
            result.blockers.indices);

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: build AABB trees for raycasting enforcers/blockers: end";
}

struct DefaultSeamComparator {
    //is a better?
    bool operator()(const SeamCandidate &a, const SeamCandidate &b) const {
        if (a.m_type > b.m_type) {
            return true;
        }
        if (b.m_type > a.m_type) {
            return false;
        }

//        if (a.m_overhang > 0.2 && b.m_overhang < a.m_overhang) {
//            return false;
//        }
//
//        if (b.m_ccw_angle < -float(0.3 * PI) && a.m_ccw_angle > -float(0.3 * PI)){
//            return false;
//        }

        if (*b.m_nearby_seam_points > *a.m_nearby_seam_points) {
            return false;
        }

        if (b.m_visibility < 1.2*a.m_visibility) {
            return false;
        }


        return true;
    }
}
;

} // namespace SeamPlacerImpl

void SeamPlacer::init(const Print &print) {
using namespace SeamPlacerImpl;
m_perimeter_points_trees_per_object.clear();
m_perimeter_points_per_object.clear();

for (const PrintObject *po : print.objects()) {

    GlobalModelInfo global_model_info { };
    gather_global_model_info(global_model_info, po);

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather and build KD trees with seam candidates: start";

    m_perimeter_points_per_object.emplace(po, po->layer_count());
    m_perimeter_points_trees_per_object.emplace(po, po->layer_count());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, po->layers().size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    std::vector<SeamCandidate> &layer_candidates = m_perimeter_points_per_object[po][layer_idx];
                    const auto layer = po->get_layer(layer_idx);
                    auto unscaled_z = layer->slice_z;
                    for (const LayerRegion *layer_region : layer->regions()) {
                        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                            Polygons polygons;
                            if (ex_entity->is_collection()) { //collection of inner, outer, and overhang perimeters
                                for (const ExtrusionEntity *perimeter :
                static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter) {
                        perimeter->polygons_covered_by_width(polygons, 0);
                    }
                }
            } else {
                polygons = ex_entity->polygons_covered_by_width();
            }
            for (const auto &poly : polygons) {
                process_perimeter_polygon(poly, unscaled_z, layer_candidates,
                        global_model_info);
            }
        }
    }
    auto functor = SeamCandidateCoordinateFunctor { &layer_candidates };
    m_perimeter_points_trees_per_object[po][layer_idx] = (std::make_unique<SeamCandidatesTree>(
            functor, layer_candidates.size()));

}
}   );

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather and build KD tree with seam candidates: end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather visibility data into perimeter points : start";

    tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    for (auto &perimeter_point : m_perimeter_points_per_object[po][layer_idx]) {
                        perimeter_point.m_visibility = global_model_info.calculate_point_visibility(
                                perimeter_point.m_position, considered_hits_distance);
                    }
                }
            });

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather visibility data into perimeter points : end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: compute overhangs : start";

    tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    for (SeamCandidate &perimeter_point : m_perimeter_points_per_object[po][layer_idx]) {
                        if (layer_idx > 0) {
                            size_t closest_supporter = find_closest_point(
                                    *m_perimeter_points_trees_per_object[po][layer_idx - 1],
                                    perimeter_point.m_position);
                            const SeamCandidate &supporter_point =
                                    m_perimeter_points_per_object[po][layer_idx - 1][closest_supporter];

                            auto prev_next = find_previous_and_next_perimeter_point(m_perimeter_points_per_object[po][layer_idx-1], closest_supporter);
                            const SeamCandidate &prev_point =
                                    m_perimeter_points_per_object[po][layer_idx - 1][prev_next.first];
                            const SeamCandidate &next_point =
                                    m_perimeter_points_per_object[po][layer_idx - 1][prev_next.second];

                            perimeter_point.m_overhang = calculate_overhang(perimeter_point, prev_point,
                                    supporter_point, next_point);

                        }
                    }
                }
            });

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: compute overhangs : end";

    for (size_t iteration = 0; iteration < seam_align_iterations; ++iteration) {
        if (iteration > 0) { //skip this in first iteration, no seam has been picked yet
            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: distribute seam positions to other layers : start";

            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
                    [&](tbb::blocked_range<size_t> r) {
                        for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                            std::vector<SeamCandidate> &layer_perimeter_points =
                                    m_perimeter_points_per_object[po][layer_idx];
                            size_t current = 0;
                            while (current < layer_perimeter_points.size()) {
                                auto seam_position =
                                        layer_perimeter_points[layer_perimeter_points[current].m_seam_index].m_position;

                                size_t other_layer_idx_start = std::max(
                                        (int) layer_idx - (int) seam_align_layer_dist, 0);
                                size_t other_layer_idx_end = std::min(layer_idx + seam_align_layer_dist,
                                        m_perimeter_points_per_object[po].size() - 1);

                                for (size_t other_layer_idx = other_layer_idx_start;
                                        other_layer_idx <= other_layer_idx_end; ++other_layer_idx) {

                                    size_t closest_point_idx = find_closest_point(
                                            *m_perimeter_points_trees_per_object[po][other_layer_idx],
                                            seam_position);

                                    m_perimeter_points_per_object[po][other_layer_idx][closest_point_idx].m_nearby_seam_points->fetch_add(
                    1, std::memory_order_relaxed);

        }

                                current += layer_perimeter_points[current].m_polygon_index_reverse + 1;
                            }
                        }
                    });

            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: distribute seam positions to other layers : end";
        }

        BOOST_LOG_TRIVIAL(debug)
                        << "SeamPlacer: find seam for each perimeter polygon and store its position in each member of the polygon : start";

        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                        std::vector<SeamCandidate> &layer_perimeter_points =
                                m_perimeter_points_per_object[po][layer_idx];
                        size_t current = 0;
                        while (current < layer_perimeter_points.size()) {
                            pick_seam_point(layer_perimeter_points, current,
                                    current + layer_perimeter_points[current].m_polygon_index_reverse,
                                    DefaultSeamComparator{});
                            current += layer_perimeter_points[current].m_polygon_index_reverse + 1;
                        }
                    }
                });

        BOOST_LOG_TRIVIAL(debug)
                        << "SeamPlacer: find seam for each perimeter polygon and store its position in each member of the polygon : end";

    }
}
}

void SeamPlacer::place_seam(const PrintObject *po, ExtrusionLoop &loop, coordf_t unscaled_z, int layer_index,
    bool external_first) {
assert(m_perimeter_points_trees_per_object.find(po) != nullptr);
assert(m_perimeter_points_per_object.find(po) != nullptr);

assert(layer_index >= 0);
const auto &perimeter_points_tree = *m_perimeter_points_trees_per_object[po][layer_index];
const auto &perimeter_points = m_perimeter_points_per_object[po][layer_index];

const Point &fp = loop.first_point();

//This is backup check, so that slicer does not crash if something weird is going on
if (perimeter_points.empty()) {
    BOOST_LOG_TRIVIAL(error)
                << "SeamPlacer: Trying to place seam for index which does not contain any outer or overhang perimeter points, maybe new perimeter type option?";
    loop.split_at(fp, true);
} else {
    auto unscaled_p = unscale(fp);
    auto closest_perimeter_point_index = find_closest_point(perimeter_points_tree,
            Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z });
    size_t perimeter_seam_index = perimeter_points[closest_perimeter_point_index].m_seam_index;
    Vec3d seam_position = perimeter_points[perimeter_seam_index].m_position;

    Point seam_point = scaled(Vec2d { seam_position.x(), seam_position.y() });

    if (!loop.split_at_vertex(seam_point))
        // The point is not in the original loop.
        // Insert it.
        loop.split_at(seam_point, true);
}
}

#ifdef DEBUG_FILES
    Slic3r::CNumericLocalesSetter locales_setter;
    FILE *fp = boost::nowide::fopen("perimeters.obj", "w");
    if (fp == nullptr) {
        BOOST_LOG_TRIVIAL(error)
        << "Couldn't open " << "perimeters.obj" << " for writing";
    }

    for (size_t i = 0; i < perimeter_points.size(); ++i)
    fprintf(fp, "v %f %f %f %f\n", perimeter_points[i].m_position[0], perimeter_points[i].m_position[1],
            perimeter_points[i].m_position[2], perimeter_points[i].m_visibility);
    fclose(fp);
#endif

//TODO disable, only debug code
#ifdef DEBUG_FILES
        its_write_obj(triangles, "triangles.obj");

        Slic3r::CNumericLocalesSetter locales_setter;
        FILE *fp = boost::nowide::fopen("hits.obj", "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Couldn't open " << "hits.obj" << " for writing";
        }

        for (size_t i = 0; i < hit_points.size(); ++i)
            fprintf(fp, "v %f %f %f \n", hit_points[i].m_position[0], hit_points[i].m_position[1],
                    hit_points[i].m_position[2]);
        fclose(fp);
    #endif

} // namespace Slic3r
