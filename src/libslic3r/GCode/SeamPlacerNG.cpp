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
                        Vec3f surface_normal = its_face_normal(triangles, hitpoint.id);

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
        float min_arm_length) {
    std::vector<float> result(polygon.size());

    if (polygon.size() == 1) {
        result[0] = 0.0f;
    }

    auto make_idx_circular = [&](int index) {
        while (index < 0) {
            index += polygon.size();
        }
        return index % polygon.size();
    };

    int idx_prev = 0;
    int idx_curr = 0;
    int idx_next = 0;

    float distance_to_prev = 0;
    float distance_to_next = 0;

    //push idx_prev far enough back as initialization
    while (distance_to_prev < min_arm_length) {
        idx_prev = make_idx_circular(idx_prev - 1);
        distance_to_prev += lengths[idx_prev];
    }

    for (size_t _i = 0; _i < polygon.size(); ++_i) {
        // pull idx_prev to current as much as possible, while respecting the min_arm_length
        while (distance_to_prev - lengths[idx_prev] > min_arm_length) {
            distance_to_prev -= lengths[idx_prev];
            idx_prev = make_idx_circular(idx_prev + 1);
        }

        //push idx_next forward as far as needed
        while (distance_to_next < min_arm_length) {
            distance_to_next += lengths[idx_next];
            idx_next = make_idx_circular(idx_next + 1);
        }

        // Calculate angle between idx_prev, idx_curr, idx_next.
        const Point &p0 = polygon.points[idx_prev];
        const Point &p1 = polygon.points[idx_curr];
        const Point &p2 = polygon.points[idx_next];
        const Point v1 = p1 - p0;
        const Point v2 = p2 - p1;
        int64_t dot = int64_t(v1(0)) * int64_t(v2(0)) + int64_t(v1(1)) * int64_t(v2(1));
        int64_t cross = int64_t(v1(0)) * int64_t(v2(1)) - int64_t(v1(1)) * int64_t(v2(0));
        float angle = float(atan2(float(cross), float(dot)));
        result[idx_curr] = angle;

        // increase idx_curr by one
        float curr_distance = lengths[idx_curr];
        idx_curr++;
        distance_to_prev += curr_distance;
        distance_to_next -= curr_distance;
    }

    return result;
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
        size_t closest_point_index = find_closest_point(raycast_hits_tree, position);
        if (closest_point_index == raycast_hits_tree.npos
                ||
                (position - geometry_raycast_hits[closest_point_index].position).norm()
                        > SeamPlacer::seam_align_tolerable_dist) {
            return 0;
        }
        auto nearby_points = find_nearby_points(raycast_hits_tree, position, SeamPlacer::considered_area_radius);
        Vec3f local_normal = geometry_raycast_hits[closest_point_index].surface_normal;

        float visibility = 0;
        for (const auto &hit_point_index : nearby_points) {
            // The further away from the perimeter point,
            // the less representative ray hit is
            float distance =
                    (position - geometry_raycast_hits[hit_point_index].position).norm();
            visibility += (SeamPlacer::considered_area_radius - distance) *
                    std::max(0.0f, local_normal.dot(geometry_raycast_hits[hit_point_index].surface_normal));
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
                float normalized = visibility
                        / (SeamPlacer::expected_hits_per_area * SeamPlacer::considered_area_radius);
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

            for (size_t i = 0; i < geometry_raycast_hits.size(); ++i) {
                Vec3f surface_normal = (geometry_raycast_hits[i].surface_normal + Vec3f(1.0, 1.0, 1.0)) / 2.0;
                fprintf(fp, "v %f %f %f %f %f %f \n", geometry_raycast_hits[i].position[0],
                        geometry_raycast_hits[i].position[1],
                        geometry_raycast_hits[i].position[2], surface_normal[0], surface_normal[1],
                        surface_normal[2]);
            }
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
    if (orig_polygon.size() == 0) {
        return;
    }

    Polygon polygon = orig_polygon;
    bool was_clockwise = polygon.make_counter_clockwise();

    std::vector<float> lengths { };
    for (size_t point_idx = 0; point_idx < polygon.size() - 1; ++point_idx) {
        lengths.push_back(std::max((unscale(polygon[point_idx]) - unscale(polygon[point_idx + 1])).norm(), 0.01));
    }
    lengths.push_back(std::max((unscale(polygon[0]) - unscale(polygon[polygon.size() - 1])).norm(), 0.01));

    std::vector<float> local_angles = calculate_polygon_angles_at_vertices(polygon, lengths,
            SeamPlacer::polygon_local_angles_arm_distance);
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

        float local_ccw_angle = was_clockwise ? -local_angles[index] : local_angles[index];

        result_vec.emplace_back(unscaled_position, perimeter, local_ccw_angle, type);
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
// value represents distance from edge, positive is overhang, negative is inside shape
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

    if (under_b.local_ccw_angle > 0 && dist_ab > 0 && dist_bc > 0) { //convex shape, p is inside
        return -((p - b).norm() + dist_ab + dist_bc) / 3.0;
    }

    if (under_b.local_ccw_angle < 0 && (dist_ab < 0 || dist_bc < 0)) { //concave shape, p is inside
        return -((p - b).norm() + dist_ab + dist_bc) / 3.0;
    }

    return ((p - b).norm() + dist_ab + dist_bc) / 3.0;
}

template<typename Comparator>
void pick_seam_point(std::vector<SeamCandidate> &perimeter_points, size_t start_index,
        const Comparator &comparator) {
    size_t end_index = perimeter_points[start_index].perimeter->end_index;
    size_t seam_index = start_index;
    for (size_t index = start_index + 1; index <= end_index; ++index) {
        if (comparator.is_first_better(perimeter_points[index], perimeter_points[seam_index])) {
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
    static constexpr float angle_clusters[] { -1.0, 0.4 * PI, 0.5 * PI, 0.6
            * PI, 0.7 * PI, 0.8 * PI, 0.9 * PI };

    const float get_angle_category(float ccw_angle) const {
        float concave_bonus = ccw_angle < 0 ? 0.1 * PI : 0;
        float abs_angle = abs(ccw_angle) + concave_bonus;
        auto category = std::find_if_not(std::begin(angle_clusters), std::end(angle_clusters),
                [&](float category_limit) {
                    return abs_angle > category_limit;
                });
        category--;
        return *category;
    }

    bool is_first_better(const SeamCandidate &a, const SeamCandidate &b) const {
        // Blockers/Enforcers discrimination, top priority
        if (a.type > b.type) {
            return true;
        }
        if (b.type > a.type) {
            return false;
        }

        //avoid overhangs
        if (a.overhang > 0.3f && b.overhang < a.overhang) {
            return false;
        }

        { //local angles
            float a_local_category = get_angle_category(a.local_ccw_angle);
            float b_local_category = get_angle_category(b.local_ccw_angle);

            if (a_local_category > b_local_category) {
                return true;
            }
            if (a_local_category < b_local_category) {
                return false;
            }
        }

        return a.visibility < b.visibility;
    }

    bool is_first_not_much_worse(const SeamCandidate &a, const SeamCandidate &b) const {
        // Blockers/Enforcers discrimination, top priority
        if (a.type > b.type) {
            return true;
        }
        if (b.type > a.type) {
            return false;
        }

        //avoid overhangs
        if (a.overhang > 0.3f && b.overhang < a.overhang) {
            return false;
        }

        { //local angles
            float a_local_category = get_angle_category(a.local_ccw_angle) + 0.2 * PI; //give a slight bonus
            float b_local_category = get_angle_category(b.local_ccw_angle);

            if (a_local_category > b_local_category) {
                return true;
            }
            if (a_local_category < b_local_category) {
                return false;
            }
        }

        return a.visibility < b.visibility * 1.5;
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
                    m_perimeter_points_trees_per_object[po][layer_idx] = std::make_unique<SeamCandidatesTree>(
                            functor, layer_candidates.size());
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
                        perimeter_point.visibility = global_model_info.calculate_point_visibility(
                                perimeter_point.position);
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
                        const auto calculate_layer_overhang = [&](size_t other_layer_idx) {
                            size_t closest_supporter = find_closest_point(
                                    *m_perimeter_points_trees_per_object[po][other_layer_idx],
                                    perimeter_point.position);
                            const SeamCandidate &supporter_point =
                                    m_perimeter_points_per_object[po][other_layer_idx][closest_supporter];

                            auto prev_next = find_previous_and_next_perimeter_point(m_perimeter_points_per_object[po][other_layer_idx], closest_supporter);
                            const SeamCandidate &prev_point =
                                    m_perimeter_points_per_object[po][other_layer_idx][prev_next.first];
                            const SeamCandidate &next_point =
                                    m_perimeter_points_per_object[po][other_layer_idx][prev_next.second];

                            return calculate_overhang(perimeter_point, prev_point,
                                    supporter_point, next_point);
                        };

                        if (layer_idx > 0) { //calculate overhang
                            perimeter_point.overhang = calculate_layer_overhang(layer_idx-1);
                        }
                        if (layer_idx < m_perimeter_points_per_object[po].size() - 1) { //calculate higher_layer_overhang
                            perimeter_point.higher_layer_overhang = calculate_layer_overhang(layer_idx+1);
                        }
                    }
                }
            });
        }

// sadly cannot be const because map access operator[] is not const, since it can create new object
template<typename Comparator>
bool SeamPlacer::find_next_seam_in_string(const PrintObject *po, Vec3f &last_point_pos,
        size_t layer_idx, const Comparator &comparator,
        std::vector<std::pair<size_t, size_t>> &seam_string,
        std::vector<std::pair<size_t, size_t>> &potential_string_seams) {
    using namespace SeamPlacerImpl;

    Vec3f projected_position { last_point_pos.x(), last_point_pos.y(), float(
            po->get_layer(layer_idx)->slice_z) };
    //find closest point in next layer
    size_t closest_point_index = find_closest_point(
            *m_perimeter_points_trees_per_object[po][layer_idx], projected_position);

    SeamCandidate &closest_point = m_perimeter_points_per_object[po][layer_idx][closest_point_index];

    if (closest_point.perimeter->aligned) { //already aligned, skip
        return false;
    }

    //from the closest point, deduce index of seam in the next layer
    SeamCandidate &next_layer_seam =
            m_perimeter_points_per_object[po][layer_idx][closest_point.perimeter->seam_index];

    if ((next_layer_seam.position - projected_position).norm()
            < SeamPlacer::seam_align_tolerable_dist) { //seam point is within limits, put in the close_by_points list
        seam_string.emplace_back(layer_idx, closest_point.perimeter->seam_index);
        last_point_pos = next_layer_seam.position;
        return true;
    } else if ((closest_point.position - projected_position).norm()
            < SeamPlacer::seam_align_tolerable_dist
            && comparator.is_first_not_much_worse(closest_point, next_layer_seam)) {
        //seam point is far, but if the close point is not much worse, do not count it as a skip and add it to potential_string_seams
        potential_string_seams.emplace_back(layer_idx, closest_point_index);
        last_point_pos = closest_point.position;
        return true;
    } else {
        return false;
    }

}

//https://towardsdatascience.com/least-square-polynomial-fitting-using-c-eigen-package-c0673728bd01
template<typename Comparator>
void SeamPlacer::align_seam_points(const PrintObject *po, const Comparator &comparator) {
    using namespace SeamPlacerImpl;

    for (size_t layer_idx = 0; layer_idx < m_perimeter_points_per_object[po].size(); ++layer_idx) {
        std::vector<SeamCandidate> &layer_perimeter_points =
                m_perimeter_points_per_object[po][layer_idx];
        size_t current_point_index = 0;
        while (current_point_index < layer_perimeter_points.size()) {
            if (layer_perimeter_points[current_point_index].perimeter->aligned) {
                //skip
            } else {
                int skips = SeamPlacer::seam_align_tolerable_skips;
                int next_layer = layer_idx + 1;
                Vec3f last_point_pos = layer_perimeter_points[current_point_index].position;

                std::vector<std::pair<size_t, size_t>> seam_string;
                std::vector<std::pair<size_t, size_t>> potential_string_seams;

                //find close by points and outliers; there is a budget of skips allowed
                while (skips >= 0 && next_layer < int(m_perimeter_points_per_object[po].size())) {
                    if (find_next_seam_in_string(po, last_point_pos, next_layer, comparator, seam_string,
                            potential_string_seams)) {
                        //String added, last_point_pos updated, nothing to be done
                    } else {
                        // Layer skipped, reduce number of available skips
                        skips--;
                    }
                    next_layer++;
                }

                if (seam_string.size() >= seam_align_minimum_string_seams) { //string long enough to be worth aligning
                    //do additional check in back direction
                    next_layer = layer_idx - 1;
                    skips = SeamPlacer::seam_align_tolerable_skips;
                    while (skips >= 0 && next_layer >= 0) {
                        if (find_next_seam_in_string(po, last_point_pos, next_layer, comparator,
                                seam_string,
                                potential_string_seams)) {
                            //String added, last_point_pos updated, nothing to be done
                        } else {
                            // Layer skipped, reduce number of available skips
                            skips--;
                        }
                        next_layer--;
                    }

                    // all string seams and potential string seams gathered, now do the alignment
                    seam_string.insert(seam_string.end(), potential_string_seams.begin(), potential_string_seams.end());
                    std::sort(seam_string.begin(), seam_string.end(),
                            [](const std::pair<size_t, size_t> &left, const std::pair<size_t, size_t> &right) {
                                return left.first < right.first;
                            });

                    //https://en.wikipedia.org/wiki/Exponential_smoothing
                    //inititalization
                    float smoothing_factor = SeamPlacer::seam_align_strength;
                    std::pair<size_t, size_t> init = seam_string[0];
                    Vec2f prev_pos_xy = m_perimeter_points_per_object[po][init.first][init.second].position.head<2>();
                    for (const auto &pair : seam_string) {
                        Vec3f current_pos = m_perimeter_points_per_object[po][pair.first][pair.second].position;
                        float current_height = current_pos.z();
                        Vec2f current_pos_xy = current_pos.head<2>();
                        current_pos_xy = smoothing_factor * prev_pos_xy + (1.0 - smoothing_factor) * current_pos_xy;

                        Perimeter *perimeter =
                                m_perimeter_points_per_object[po][pair.first][pair.second].perimeter.get();
                        perimeter->final_seam_position =
                                Vec3f { current_pos_xy.x(), current_pos_xy.y(), current_height };
                        perimeter->aligned = true;
                        prev_pos_xy = current_pos_xy;
                    }

                } // else string is not long enough, so dont do anything
            }
            current_point_index = layer_perimeter_points[current_point_index].perimeter->end_index + 1;
        }
    }

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
        << "SeamPlacer: pick_seam_point : start";
        //pick seam point
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                        std::vector<SeamCandidate> &layer_perimeter_points =
                                m_perimeter_points_per_object[po][layer_idx];
                        size_t current = 0;
                        while (current < layer_perimeter_points.size()) {
                            pick_seam_point(layer_perimeter_points, current, DefaultSeamComparator { });
                            current = layer_perimeter_points[current].perimeter->end_index + 1;
                        }
                    }
                });
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: pick_seam_point : end";

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: align_seam_points : start";
        align_seam_points(po, DefaultSeamComparator { });
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: align_seam_points : end";

    }
}

void SeamPlacer::place_seam(const Layer *layer, ExtrusionLoop &loop, bool external_first) {
    using namespace SeamPlacerImpl;
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
    const Perimeter *perimeter = perimeter_points[closest_perimeter_point_index].perimeter.get();
    Vec3f seam_position = perimeter_points[perimeter->seam_index].position;
    if (perimeter->aligned) {
        seam_position = perimeter->final_seam_position;
    }

    Point seam_point = scaled(Vec2d { seam_position.x(), seam_position.y() });

    if (!loop.split_at_vertex(seam_point))
// The point is not in the original loop.
// Insert it.
        loop.split_at(seam_point, true);
}

} // namespace Slic3r
