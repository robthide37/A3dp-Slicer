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

#define DEBUG_FILES

#ifdef DEBUG_FILES
#include "Subdivide.hpp"
#include <boost/nowide/cstdio.hpp>
#endif

namespace Slic3r {

namespace SeamPlacerImpl {

/// Coordinate frame
class Frame {
public:
    Frame() {
        mX = Vec3f(1, 0, 0);
        mY = Vec3f(0, 1, 0);
        mZ = Vec3f(0, 0, 1);
    }

    Frame(const Vec3f &x, const Vec3f &y, const Vec3f &z) :
            mX(x), mY(y), mZ(z) {
    }

    void set_from_z(const Vec3f &z) {
        mZ = z.normalized();
        Vec3f tmpZ = mZ;
        Vec3f tmpX = (std::abs(tmpZ.x()) > 0.99f) ? Vec3f(0, 1, 0) : Vec3f(1, 0, 0);
        mY = (tmpZ.cross(tmpX)).normalized();
        mX = mY.cross(tmpZ);
    }

    Vec3f to_world(const Vec3f &a) const {
        return a.x() * mX + a.y() * mY + a.z() * mZ;
    }

    Vec3f to_local(const Vec3f &a) const {
        return Vec3f(mX.dot(a), mY.dot(a), mZ.dot(a));
    }

    const Vec3f& binormal() const {
        return mX;
    }

    const Vec3f& tangent() const {
        return mY;
    }

    const Vec3f& normal() const {
        return mZ;
    }

private:
    Vec3f mX, mY, mZ;
};

Vec3f sample_sphere_uniform(const Vec2f &samples) {
    float term1 = 2.0f * float(PI) * samples.x();
    float term2 = 2.0f * sqrt(samples.y() - samples.y() * samples.y());
    return {cos(term1) * term2, sin(term1) * term2,
        1.0f - 2.0f * samples.y()};
}

Vec3f sample_power_cosine_hemisphere(const Vec2f &samples, float power) {
    float term1 = 2.f * float(PI) * samples.x();
    float term2 = pow(samples.y(), 1.f / (power + 1.f));
    float term3 = sqrt(1.f - term2 * term2);

    return Vec3f(cos(term1) * term3, sin(term1) * term3, term2);
}

std::vector<HitInfo> raycast_visibility(const AABBTreeIndirect::Tree<3, float> &raycasting_tree,
        const indexed_triangle_set &triangles) {

    auto bbox = raycasting_tree.node(0).bbox;
    Vec3f vision_sphere_center = bbox.center().cast<float>();
    Vec3f side_sizes = bbox.sizes().cast<float>();
    float vision_sphere_raidus = (side_sizes.norm() * 0.55); // 0.5 (half) covers whole object,
    // 0.05 added to avoid corner cases
    // very rough approximation of object surface area from its bounding sphere
    float approx_area = 4 * PI * vision_sphere_raidus * vision_sphere_raidus;
    float local_considered_area = PI * SeamPlacer::considered_area_radius * SeamPlacer::considered_area_radius;
    size_t ray_count = SeamPlacer::expected_hits_per_area * approx_area / local_considered_area;

    // Prepare random samples per ray
    // std::random_device rnd_device;
    // use fixed seed, we can backtrack potential issues easier
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
                    Vec3f global_ray_dir = sample_sphere_uniform(global_dir_random_samples[index]);
                    Vec3f ray_origin = (vision_sphere_center - global_ray_dir * vision_sphere_raidus);
                    Vec3f local_dir = sample_power_cosine_hemisphere(local_dir_random_samples[index], SeamPlacer::cosine_hemisphere_sampling_power);

                    Frame f;
                    f.set_from_z(global_ray_dir);
                    Vec3f final_ray_dir = (f.to_world(local_dir));

                    igl::Hit hitpoint;
                    // FIXME: This AABBTTreeIndirect query will not compile for float ray origin and
                    // direction.
                    Vec3d ray_origin_d = ray_origin.cast<double>();
                    Vec3d final_ray_dir_d = final_ray_dir.cast<double>();
                    bool hit = AABBTreeIndirect::intersect_ray_first_hit(triangles.vertices,
                            triangles.indices, raycasting_tree, ray_origin_d, final_ray_dir_d, hitpoint);

                    if (hit) {
                        auto face = triangles.indices[hitpoint.id];
                        auto edge1 = triangles.vertices[face[1]] - triangles.vertices[face[0]];
                        auto edge2 = triangles.vertices[face[2]] - triangles.vertices[face[0]];

                        Vec3f hit_pos = (triangles.vertices[face[0]] + edge1 * hitpoint.u + edge2 * hitpoint.v);
                        Vec3f surface_normal = edge1.cross(edge2).normalized();

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
        float angle = float(atan2(float(cross), float(dot)));
        angles[idx_curr] = angle;
    }

    return angles;
}

struct GlobalModelInfo {
    std::vector<HitInfo> geometry_raycast_hits;
    KDTreeIndirect<3, float, HitInfoCoordinateFunctor> raycast_hits_tree;
    indexed_triangle_set enforcers;
    indexed_triangle_set blockers;
    AABBTreeIndirect::Tree<3, float> enforcers_tree;
    AABBTreeIndirect::Tree<3, float> blockers_tree;

    GlobalModelInfo() :
            raycast_hits_tree(HitInfoCoordinateFunctor { &geometry_raycast_hits }) {
    }

    bool is_enforced(const Vec3f &position) const {
        float radius = SeamPlacer::enforcer_blocker_sqr_distance_tolerance;
        return AABBTreeIndirect::is_any_triangle_in_radius(enforcers.vertices, enforcers.indices,
                enforcers_tree, position, radius);
    }

    bool is_blocked(const Vec3f &position) const {
        float radius = SeamPlacer::enforcer_blocker_sqr_distance_tolerance;
        return AABBTreeIndirect::is_any_triangle_in_radius(blockers.vertices, blockers.indices,
                blockers_tree, position, radius);
    }

    float calculate_point_visibility(const Vec3f &position) const {
        auto nearby_points = find_nearby_points(raycast_hits_tree, position, SeamPlacer::considered_area_radius);
        float visibility = 0;
        for (const auto &hit_point_index : nearby_points) {
            float distance =
                    (position - geometry_raycast_hits[hit_point_index].position).norm();
            visibility += SeamPlacer::considered_area_radius - distance; // The further away from the perimeter point,
            // the less representative ray hit is
        }
        return visibility;

    }

#ifdef DEBUG_FILES
    void debug_export(const indexed_triangle_set &obj_mesh, const char *file_name) const {
        indexed_triangle_set divided_mesh = subdivide(obj_mesh, SeamPlacer::considered_area_radius);
        Slic3r::CNumericLocalesSetter locales_setter;
        {
            FILE *fp = boost::nowide::fopen(file_name, "w");
            if (fp == nullptr) {
                BOOST_LOG_TRIVIAL(error)
                << "stl_write_obj: Couldn't open " << file_name << " for writing";
                return;
            }

            const auto vis_to_rgb = [](float normalized_visibility) {
                float ratio = 2 * normalized_visibility;
                float blue = std::max(0.0f, 1.0f - ratio);
                float red = std::max(0.0f, ratio - 1.0f);
                float green = std::max(0.0f, 1.0f - blue - red);
                return Vec3f { red, blue, green };
            };

            for (size_t i = 0; i < divided_mesh.vertices.size(); ++i) {
                float visibility = calculate_point_visibility(divided_mesh.vertices[i]);
                float normalized = visibility / SeamPlacer::expected_hits_per_area / 2.0;
                Vec3f color = vis_to_rgb(normalized);
                fprintf(fp, "v %f %f %f  %f %f %f\n",
                        divided_mesh.vertices[i](0), divided_mesh.vertices[i](1), divided_mesh.vertices[i](2),
                        color(0), color(1), color(2)
                                );
            }
            for (size_t i = 0; i < divided_mesh.indices.size(); ++i)
                fprintf(fp, "f %d %d %d\n", divided_mesh.indices[i][0] + 1, divided_mesh.indices[i][1] + 1,
                        divided_mesh.indices[i][2] + 1);
            fclose(fp);
        }
        {
            auto fname = std::string("hits_").append(file_name);
            FILE *fp = boost::nowide::fopen(fname.c_str(), "w");
            if (fp == nullptr) {
                BOOST_LOG_TRIVIAL(error)
                << "Couldn't open " << fname << " for writing";
            }

            for (size_t i = 0; i < geometry_raycast_hits.size(); ++i)
                fprintf(fp, "v %f %f %f \n", geometry_raycast_hits[i].position[0], geometry_raycast_hits[i].position[1],
                        geometry_raycast_hits[i].position[2]);
            fclose(fp);
        }
    }
#endif

}
;

Polygons extract_perimeter_polygons(const Layer *layer) {
    Polygons polygons;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            if (ex_entity->is_collection()) { //collection of inner, outer, and overhang perimeters
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter) {
                        Points p;
                        perimeter->collect_points(p);
                        polygons.emplace_back(p);
                    }
                }
                if (polygons.empty()) {
                    Points p;
                    ex_entity->collect_points(p);
                    polygons.emplace_back(p);
                }
            } else {
                Points p;
                ex_entity->collect_points(p);
                polygons.emplace_back(p);
            }
        }
    }

    if (polygons.empty()) {
        polygons.emplace_back(std::vector { Point { 0, 0 } });
    }

    return polygons;
}

void process_perimeter_polygon(const Polygon &orig_polygon, float z_coord, std::vector<SeamCandidate> &result_vec,
        const GlobalModelInfo &global_model_info) {
    Polygon polygon = orig_polygon;
    bool was_clockwise = polygon.make_counter_clockwise();

    std::vector<float> lengths = polygon.parameter_by_length();
    std::vector<float> angles = calculate_polygon_angles_at_vertices(polygon, lengths,
            SeamPlacer::polygon_angles_arm_distance);
    std::shared_ptr<Perimeter> perimeter = std::make_shared<Perimeter>();

    perimeter->start_index = result_vec.size();
    perimeter->end_index = result_vec.size() + polygon.size() - 1;

    for (size_t index = 0; index < polygon.size(); ++index) {
        Vec2f unscaled_p = unscale(polygon[index]).cast<float>();
        Vec3f unscaled_position = Vec3f { unscaled_p.x(), unscaled_p.y(), z_coord };
        EnforcedBlockedSeamPoint type = EnforcedBlockedSeamPoint::Neutral;

        if (global_model_info.is_enforced(unscaled_position)) {
            type = EnforcedBlockedSeamPoint::Enforced;
        }

        if (global_model_info.is_blocked(unscaled_position)) {
            type = EnforcedBlockedSeamPoint::Blocked;
        }

        float ccw_angle = was_clockwise ? -angles[index] : angles[index];

        result_vec.emplace_back(unscaled_position, perimeter, ccw_angle, type);
    }
}

std::pair<size_t, size_t> find_previous_and_next_perimeter_point(const std::vector<SeamCandidate> &perimeter_points,
        size_t index) {
    const SeamCandidate &current = perimeter_points[index];
    int prev = index - 1;
    int next = index + 1;

    if (index == current.perimeter->start_index) {
        prev = current.perimeter->end_index;
    }

    if (index == current.perimeter->end_index) {
        next = current.perimeter->start_index;
    }

    assert(prev >= 0);
    assert(next >= 0);
    return {size_t(prev),size_t(next)};
}

//NOTE: only rough esitmation of overhang distance
float calculate_overhang(const SeamCandidate &point, const SeamCandidate &under_a, const SeamCandidate &under_b,
        const SeamCandidate &under_c) {
    auto p = Vec2d { point.position.x(), point.position.y() };
    auto a = Vec2d { under_a.position.x(), under_a.position.y() };
    auto b = Vec2d { under_b.position.x(), under_b.position.y() };
    auto c = Vec2d { under_c.position.x(), under_c.position.y() };

    auto oriented_line_dist = [](const Vec2d a, const Vec2d b, const Vec2d p) {
        return -((b.x() - a.x()) * (a.y() - p.y()) - (a.x() - p.x()) * (b.y() - a.y())) / (a - b).norm();
    };

    auto dist_ab = oriented_line_dist(a, b, p);
    auto dist_bc = oriented_line_dist(b, c, p);

    if (under_b.ccw_angle > 0 && dist_ab > 0 && dist_bc > 0) { //convex shape, p is inside
        return 0;
    }

    if (under_b.ccw_angle < 0 && (dist_ab < 0 || dist_bc < 0)) { //concave shape, p is inside
        return 0;
    }

    return (p - b).norm();
}

template<typename CompareFunc>
void pick_seam_point(std::vector<SeamCandidate> &perimeter_points, size_t start_index,
        const CompareFunc &is_first_better) {
    size_t end_index = perimeter_points[start_index].perimeter->end_index;
    size_t seam_index = start_index;
    for (size_t index = start_index + 1; index <= end_index; ++index) {
        if (is_first_better(perimeter_points[index], perimeter_points[seam_index])) {
            seam_index = index;
        }
    }

    perimeter_points[start_index].perimeter->seam_index = seam_index;
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

    result.geometry_raycast_hits = raycast_visibility(raycasting_tree, triangle_set);
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

#ifdef DEBUG_FILES
    auto filename = "visiblity_of_" + std::to_string(po->id().id) + ".obj";
    result.debug_export(triangle_set, filename.c_str());
#endif
}

struct DefaultSeamComparator {
    //is A better?
    bool operator()(const SeamCandidate &a, const SeamCandidate &b) const {
        // Blockers/Enforcers discrimination, top priority
        if (a.type > b.type) {
            return true;
        }
        if (b.type > a.type) {
            return false;
        }

        //avoid overhangs
        if (a.overhang > 0.5f && b.overhang < a.overhang) {
            return false;
        }

        auto angle_score = [](float ccw_angle) {
            if (ccw_angle > 0) {
                float normalized = (ccw_angle / float(PI)) * 0.9f;
                return normalized;
            } else {
                float normalized = (-ccw_angle / float(PI)) * 1.1f;
                return normalized;
            }
        };
        float angle_weight = 2.0f;

        auto vis_score = [](float visibility) {
            return (1.0f - visibility / SeamPlacer::expected_hits_per_area);
        };
        float vis_weight = 1.2f;

        float score_a = angle_score(a.ccw_angle) * angle_weight +
                vis_score(a.visibility) * vis_weight;
        float score_b = angle_score(b.ccw_angle) * angle_weight +
                vis_score(b.visibility) * vis_weight;

        if (score_a > score_b)
            return true;
        else
            return false;
    }
}
;

} // namespace SeamPlacerImpl

void SeamPlacer::gather_seam_candidates(const PrintObject *po,
    const SeamPlacerImpl::GlobalModelInfo &global_model_info) {
using namespace SeamPlacerImpl;

m_perimeter_points_per_object.emplace(po, po->layer_count());
m_perimeter_points_trees_per_object.emplace(po, po->layer_count());

tbb::parallel_for(tbb::blocked_range<size_t>(0, po->layers().size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                std::vector<SeamCandidate> &layer_candidates = m_perimeter_points_per_object[po][layer_idx];
                const Layer *layer = po->get_layer(layer_idx);
                auto unscaled_z = layer->slice_z;
                Polygons polygons = extract_perimeter_polygons(layer);
                for (const auto &poly : polygons) {
                    process_perimeter_polygon(poly, unscaled_z, layer_candidates,
                            global_model_info);
                }
                auto functor = SeamCandidateCoordinateFunctor { &layer_candidates };
                m_perimeter_points_trees_per_object[po][layer_idx] = (std::make_unique<SeamCandidatesTree>(
                        functor, layer_candidates.size()));
            }
        }
);
}

void SeamPlacer::calculate_candidates_visibility(const PrintObject *po,
    const SeamPlacerImpl::GlobalModelInfo &global_model_info) {
using namespace SeamPlacerImpl;

tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                for (auto &perimeter_point : m_perimeter_points_per_object[po][layer_idx]) {
                    perimeter_point.visibility = global_model_info.calculate_point_visibility(perimeter_point.position);
                }
            }
        });
}

void SeamPlacer::calculate_overhangs(const PrintObject *po) {
using namespace SeamPlacerImpl;

tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                for (SeamCandidate &perimeter_point : m_perimeter_points_per_object[po][layer_idx]) {
                    if (layer_idx > 0) {
                        size_t closest_supporter = find_closest_point(
                                *m_perimeter_points_trees_per_object[po][layer_idx - 1],
                                perimeter_point.position);
                        const SeamCandidate &supporter_point =
                                m_perimeter_points_per_object[po][layer_idx - 1][closest_supporter];

                        auto prev_next = find_previous_and_next_perimeter_point(m_perimeter_points_per_object[po][layer_idx-1], closest_supporter);
                        const SeamCandidate &prev_point =
                                m_perimeter_points_per_object[po][layer_idx - 1][prev_next.first];
                        const SeamCandidate &next_point =
                                m_perimeter_points_per_object[po][layer_idx - 1][prev_next.second];

                        perimeter_point.overhang = calculate_overhang(perimeter_point, prev_point,
                                supporter_point, next_point);

                    }
                }
            }
        });
    }

void SeamPlacer::distribute_seam_positions_for_alignment(const PrintObject *po) {
    using namespace SeamPlacerImpl;

}

void SeamPlacer::init(const Print &print) {
    using namespace SeamPlacerImpl;
    m_perimeter_points_trees_per_object.clear();
    m_perimeter_points_per_object.clear();

    for (const PrintObject *po : print.objects()) {

        GlobalModelInfo global_model_info { };
        gather_global_model_info(global_model_info, po);

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: gather_seam_candidates: start";
        gather_seam_candidates(po, global_model_info);
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: gather_seam_candidates: end";

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_candidates_visibility : start";
        calculate_candidates_visibility(po, global_model_info);
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_candidates_visibility : end";

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_overhangs : start";
        calculate_overhangs(po);
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_overhangs : end";

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: distribute_seam_positions_for_alignment, pick_seams : start";

        for (size_t iteration = 0; iteration < seam_align_iterations; ++iteration) {

            if (iteration > 0) { //skip this in first iteration, no seam has been picked yet; nothing to distribute
                distribute_seam_positions_for_alignment(po);
            }

            //pick seam point
            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
                    [&](tbb::blocked_range<size_t> r) {
                        for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                            std::vector<SeamCandidate> &layer_perimeter_points =
                                    m_perimeter_points_per_object[po][layer_idx];
                            size_t current = 0;
                            while (current < layer_perimeter_points.size()) {
                                //NOTE: pick seam point function also resets the m_nearby_seam_points count on all passed points
                                pick_seam_point(layer_perimeter_points, current, DefaultSeamComparator { });
                                current = layer_perimeter_points[current].perimeter->end_index + 1;
                            }
                        }
                    });
        }
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: distribute_seam_positions_for_alignment, pick_seams : end";
    }
}

void SeamPlacer::place_seam(const Layer *layer, ExtrusionLoop &loop, bool external_first) {
    const PrintObject *po = layer->object();
    //NOTE this is necessary, since layer->id() is quite unreliable
    size_t layer_index = std::max(0, int(layer->id()) - int(po->slicing_parameters().raft_layers()));
    double unscaled_z = layer->slice_z;

    const auto &perimeter_points_tree = *m_perimeter_points_trees_per_object[po][layer_index];
    const auto &perimeter_points = m_perimeter_points_per_object[po][layer_index];

    const Point &fp = loop.first_point();

    Vec2f unscaled_p = unscale(fp).cast<float>();
    size_t closest_perimeter_point_index = find_closest_point(perimeter_points_tree,
            Vec3f { unscaled_p.x(), unscaled_p.y(), float(unscaled_z) });
    size_t perimeter_seam_index = perimeter_points[closest_perimeter_point_index].perimeter->seam_index;
    Vec3f seam_position = perimeter_points[perimeter_seam_index].position;

    Point seam_point = scaled(Vec2d { seam_position.x(), seam_position.y() });

    if (!loop.split_at_vertex(seam_point))
// The point is not in the original loop.
// Insert it.
        loop.split_at(seam_point, true);
}

// Disabled debug code, can be used to export debug data into obj files (e.g. point cloud of visibility hits)
#if 0
    #include <boost/nowide/cstdio.hpp>
    Slic3r::CNumericLocalesSetter locales_setter;
    FILE *fp = boost::nowide::fopen("perimeters.obj", "w");
    if (fp == nullptr) {
        BOOST_LOG_TRIVIAL(error)
        << "Couldn't open " << "perimeters.obj" << " for writing";
    }

    for (size_t i = 0; i < perimeter_points.size(); ++i)
    fprintf(fp, "v %f %f %f %f\n", perimeter_points[i].position[0], perimeter_points[i].position[1],
            perimeter_points[i].position[2], perimeter_points[i].visibility);
    fclose(fp);
#endif

#if 0
        its_write_obj(triangles, "triangles.obj");

        Slic3r::CNumericLocalesSetter locales_setter;
        FILE *fp = boost::nowide::fopen("hits.obj", "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Couldn't open " << "hits.obj" << " for writing";
        }

        for (size_t i = 0; i < hit_points.size(); ++i)
            fprintf(fp, "v %f %f %f \n", hit_points[i].position[0], hit_points[i].position[1],
                    hit_points[i].position[2]);
        fclose(fp);
    #endif

} // namespace Slic3r
