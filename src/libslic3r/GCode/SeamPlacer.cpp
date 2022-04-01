#include "SeamPlacer.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <random>
#include <algorithm>
#include <queue>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/QuadricEdgeCollapse.hpp"
#include "libslic3r/Subdivide.hpp"

#include "libslic3r/Geometry/Curves.hpp"

//#define DEBUG_FILES

#ifdef DEBUG_FILES
#include <boost/nowide/cstdio.hpp>
#include <SVG.hpp>
#endif

namespace Slic3r {

namespace SeamPlacerImpl {

template<typename T> int sgn(T val) {
    return int(T(0) < val) - int(val < T(0));
}

// base function: ((e^(((1)/(x^(2)+1)))-1)/(e-1))
// checkout e.g. here: https://www.geogebra.org/calculator
float gauss(float value, float mean_x_coord, float mean_value, float falloff_speed) {
    float shifted = value - mean_x_coord;
    float denominator = falloff_speed * shifted * shifted + 1.0f;
    float exponent = 1.0f / denominator;
    return mean_value * (std::exp(exponent) - 1.0f) / (std::exp(1.0f) - 1.0f);
}

Vec3f value_to_rgbf(float minimum, float maximum, float value) {
    float ratio = 2.0f * (value - minimum) / (maximum - minimum);
    float b = std::max(0.0f, (1.0f - ratio));
    float r = std::max(0.0f, (ratio - 1.0f));
    float g = 1.0f - b - r;
    return Vec3f { r, g, b };
}

Vec3i value_rgbi(float minimum, float maximum, float value) {
    return (value_to_rgbf(minimum, maximum, value) * 255).cast<int>();
}

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

Vec3f sample_hemisphere_uniform(const Vec2f &samples) {
    float term1 = 2.0f * float(PI) * samples.x();
    float term2 = 2.0f * sqrt(samples.y() - samples.y() * samples.y());
    return {cos(term1) * term2, sin(term1) * term2,
        abs(1.0f - 2.0f * samples.y())};
}

Vec3f sample_power_cosine_hemisphere(const Vec2f &samples, float power) {
    float term1 = 2.f * float(PI) * samples.x();
    float term2 = pow(samples.y(), 1.f / (power + 1.f));
    float term3 = sqrt(1.f - term2 * term2);

    return Vec3f(cos(term1) * term3, sin(term1) * term3, term2);
}

std::vector<FaceVisibilityInfo> raycast_visibility(const AABBTreeIndirect::Tree<3, float> &raycasting_tree,
        const indexed_triangle_set &triangles, size_t negative_volumes_start_index) {
    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: raycast visibility for " << triangles.indices.size() << " triangles: start";

    //prepare uniform samples of a hemisphere
    float step_size = 1.0f / SeamPlacer::sqr_rays_per_triangle;
    std::vector<Vec3f> precomputed_sample_directions(
            SeamPlacer::sqr_rays_per_triangle * SeamPlacer::sqr_rays_per_triangle);
    for (size_t x_idx = 0; x_idx < SeamPlacer::sqr_rays_per_triangle; ++x_idx) {
        float sample_x = x_idx * step_size + step_size / 2.0;
        for (size_t y_idx = 0; y_idx < SeamPlacer::sqr_rays_per_triangle; ++y_idx) {
            size_t dir_index = x_idx * SeamPlacer::sqr_rays_per_triangle + y_idx;
            float sample_y = y_idx * step_size + step_size / 2.0;
            precomputed_sample_directions[dir_index] = sample_hemisphere_uniform( { sample_x, sample_y });
        }
    }

    bool model_contains_negative_parts = negative_volumes_start_index < triangles.indices.size();

    std::vector<FaceVisibilityInfo> result(triangles.indices.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, result.size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t face_index = r.begin(); face_index < r.end(); ++face_index) {
                    FaceVisibilityInfo &dest = result[face_index];
                    dest.visibility = 1.0f;
                    constexpr float decrease = 1.0f / (SeamPlacer::sqr_rays_per_triangle * SeamPlacer::sqr_rays_per_triangle);

                    Vec3i face = triangles.indices[face_index];
                    Vec3f A = triangles.vertices[face.x()];
                    Vec3f B = triangles.vertices[face.y()];
                    Vec3f C = triangles.vertices[face.z()];
                    Vec3f center = (A + B + C) / 3.0f;
                    Vec3f normal = ((B - A).cross(C - B)).normalized();
                    // apply the local direction via Frame struct - the local_dir is with respect to +Z being forward
                    Frame f;
                    f.set_from_z(normal);

                    for (const auto &dir : precomputed_sample_directions) {
                        Vec3f final_ray_dir = (f.to_world(dir));
                        if (!model_contains_negative_parts) {
                            igl::Hit hitpoint;
                            // FIXME: This AABBTTreeIndirect query will not compile for float ray origin and
                            // direction.
                            Vec3d final_ray_dir_d = final_ray_dir.cast<double>();
                            Vec3d ray_origin_d = (center + normal * 0.1).cast<double>(); // start above surface.
                            bool hit = AABBTreeIndirect::intersect_ray_first_hit(triangles.vertices,
                                    triangles.indices, raycasting_tree, ray_origin_d, final_ray_dir_d, hitpoint);
                            if (hit) {
                                dest.visibility -= decrease;
                            }
                        } else {
                            std::vector<igl::Hit> hits;
                            int in_negative = 0;
                            if (face_index >= negative_volumes_start_index) { // if casting from negative volume face, invert direction
                                final_ray_dir = -1.0 * final_ray_dir;
                                normal = -normal;
                            }
                            Vec3d final_ray_dir_d = final_ray_dir.cast<double>();
                            Vec3d ray_origin_d = (center + normal * 0.1).cast<double>(); // start above surface.
                            bool some_hit = AABBTreeIndirect::intersect_ray_all_hits(triangles.vertices,
                                    triangles.indices, raycasting_tree,
                                    ray_origin_d, final_ray_dir_d, hits);
                            if (some_hit) {
                                // NOTE: iterating in reverse, from the last hit for one simple reason: We know the state of the ray at that point;
                                //  It cannot be inside model, and it cannot be inside negative volume
                                for (int hit_index = int(hits.size()) - 1; hit_index >= 0; --hit_index) {
                                    if (hits[hit_index].id >= int(negative_volumes_start_index)) { //negative volume hit
                                        Vec3f normal = its_face_normal(triangles, hits[hit_index].id);
                                        in_negative += sgn(normal.dot(final_ray_dir)); // if volume face aligns with ray dir, we are leaving negative space
                                        // which in reverse hit analysis means, that we are entering negative space :) and vice versa
                                    } else if (in_negative <= 0) { //object hit and we are not in negative space
                                        dest.visibility -= decrease;
                                        break;
                                    }
                                }
                            }

                        }
                    }
                }
            });

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: raycast visibility for " << triangles.indices.size() << " triangles: end";

    return result;
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

// structure to store global information about the model - occlusion hits, enforcers, blockers
struct GlobalModelInfo {
    indexed_triangle_set model;
    AABBTreeIndirect::Tree<3, float> model_tree;
    std::vector<FaceVisibilityInfo> visiblity_info;
    indexed_triangle_set enforcers;
    indexed_triangle_set blockers;
    AABBTreeIndirect::Tree<3, float> enforcers_tree;
    AABBTreeIndirect::Tree<3, float> blockers_tree;

    bool is_enforced(const Vec3f &position, float radius) const {
        if (enforcers.empty()) {
            return false;
        }
        float radius_sqr = radius * radius;
        return AABBTreeIndirect::is_any_triangle_in_radius(enforcers.vertices, enforcers.indices,
                enforcers_tree, position, radius_sqr);
    }

    bool is_blocked(const Vec3f &position, float radius) const {
        if (blockers.empty()) {
            return false;
        }
        float radius_sqr = radius * radius;
        return AABBTreeIndirect::is_any_triangle_in_radius(blockers.vertices, blockers.indices,
                blockers_tree, position, radius_sqr);
    }

    float calculate_point_visibility(const Vec3f &position) const {
        size_t hit_idx;
        Vec3f hit_point;
        if (AABBTreeIndirect::squared_distance_to_indexed_triangle_set(model.vertices, model.indices, model_tree,
                position, hit_idx, hit_point) >= 0) {
            return visiblity_info[hit_idx].visibility;
        } else {
            return 0.0f;
        }

    }

#ifdef DEBUG_FILES
    void debug_export(const indexed_triangle_set &obj_mesh, const char *file_name) const {
        indexed_triangle_set divided_mesh = obj_mesh;
        Slic3r::CNumericLocalesSetter locales_setter;

        FILE *fp = boost::nowide::fopen(file_name, "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "stl_write_obj: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < divided_mesh.vertices.size(); ++i) {
            float visibility = calculate_point_visibility(divided_mesh.vertices[i]);
            Vec3f color = value_to_rgbf(0.0f, 1.0f,
                    visibility);
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
#endif
}
;

//Extract perimeter polygons of the given layer
Polygons extract_perimeter_polygons(const Layer *layer, const SeamPosition configured_seam_preference,
        std::vector<const LayerRegion*> &corresponding_regions_out) {
    Polygons polygons;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            if (ex_entity->is_collection()) { //collection of inner, outer, and overhang perimeters
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter
                            || (perimeter->role() == ExtrusionRole::erPerimeter
                                    && configured_seam_preference == spRandom)) { //for random seam alignment, extract all perimeters
                        Points p;
                        perimeter->collect_points(p);
                        polygons.emplace_back(p);
                        corresponding_regions_out.push_back(layer_region);
                    }
                }
                if (polygons.empty()) {
                    Points p;
                    ex_entity->collect_points(p);
                    polygons.emplace_back(p);
                    corresponding_regions_out.push_back(layer_region);
                }
            } else {
                Points p;
                ex_entity->collect_points(p);
                polygons.emplace_back(p);
                corresponding_regions_out.push_back(layer_region);
            }
        }
    }

    if (polygons.empty()) { // If there are no perimeter polygons for whatever reason (disabled perimeters .. ) insert dummy point
        // it is easier than checking everywhere if the layer is not emtpy, no seam will be placed to this layer anyway
        polygons.emplace_back(std::vector { Point { 0, 0 } });
        corresponding_regions_out.push_back(nullptr);
    }

    return polygons;
}

// Insert SeamCandidates created from perimeter polygons in to the result vector.
// Compute its type (Enfrocer,Blocker), angle, and position
//each SeamCandidate also contains pointer to shared Perimeter structure representing the polygon
// if Custom Seam modifiers are present, oversamples the polygon if necessary to better fit user intentions
void process_perimeter_polygon(const Polygon &orig_polygon, float z_coord, const LayerRegion *region,
        const GlobalModelInfo &global_model_info, std::vector<SeamCandidate> &result_vec) {
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

    std::queue<Vec3f> orig_polygon_points { };
    for (size_t index = 0; index < polygon.size(); ++index) {
        Vec2f unscaled_p = unscale(polygon[index]).cast<float>();
        orig_polygon_points.emplace(unscaled_p.x(), unscaled_p.y(), z_coord);
    }
    Vec3f first = orig_polygon_points.front();
    std::queue<Vec3f> oversampled_points { };
    size_t orig_angle_index = 0;
    perimeter->start_index = result_vec.size();
    perimeter->flow_width = region != nullptr ? region->flow(FlowRole::frExternalPerimeter).width() : 0.0f;
    bool some_point_enforced = false;
    while (!orig_polygon_points.empty() || !oversampled_points.empty()) {
        EnforcedBlockedSeamPoint type = EnforcedBlockedSeamPoint::Neutral;
        Vec3f position;
        float local_ccw_angle = 0;
        bool orig_point = false;
        if (!oversampled_points.empty()) {
            position = oversampled_points.front();
            oversampled_points.pop();
        } else {
            position = orig_polygon_points.front();
            orig_polygon_points.pop();
            local_ccw_angle = was_clockwise ? -local_angles[orig_angle_index] : local_angles[orig_angle_index];
            orig_angle_index++;
            orig_point = true;
        }

        if (global_model_info.is_enforced(position, SeamPlacer::enforcer_blocker_distance_tolerance)) {
            type = EnforcedBlockedSeamPoint::Enforced;
            some_point_enforced = true;
        }

        if (global_model_info.is_blocked(position, SeamPlacer::enforcer_blocker_distance_tolerance)) {
            type = EnforcedBlockedSeamPoint::Blocked;
        }

        if (orig_point) {
            Vec3f pos_of_next = orig_polygon_points.empty() ? first : orig_polygon_points.front();
            float distance_to_next = (position - pos_of_next).norm();
            if (global_model_info.is_enforced(position, distance_to_next)
                    || global_model_info.is_blocked(position, distance_to_next)) {
                Vec3f vec_to_next = (pos_of_next - position).normalized();
                float step_size = SeamPlacer::enforcer_blocker_oversampling_distance;
                float step = step_size;
                while (step < distance_to_next) {
                    oversampled_points.push(position + vec_to_next * step);
                    step += step_size;
                }
            }
        }

        result_vec.emplace_back(position, perimeter, local_ccw_angle, type);
    }

    perimeter->end_index = result_vec.size() - 1;

    // We will find first patch of enforced points (patch: continous section of enforced points) and select the middle
    //      point, which will have priority during alignemnt
    // If there are multiple enforced patches in the perimeter, others are ignored
    if (some_point_enforced) {
        size_t first_enforced_idx = perimeter->start_index;
        while (first_enforced_idx <= perimeter->end_index
                && result_vec[first_enforced_idx].type != EnforcedBlockedSeamPoint::Enforced) {
            first_enforced_idx++;
        }
        size_t last_enforced_idx = first_enforced_idx;
        while (last_enforced_idx < perimeter->end_index
                && result_vec[last_enforced_idx + 1].type == EnforcedBlockedSeamPoint::Enforced) {
            last_enforced_idx++;
        }
        size_t central_idx = (first_enforced_idx + last_enforced_idx) / 2;
        result_vec[central_idx].central_enforcer = true;
    }

}

// Get index of previous and next perimeter point of the layer. Because SeamCandidates of all polygons of the given layer
// are sequentially stored in the vector, each perimeter contains info about start and end index. These vales are used to
// deduce index of previous and next neigbour in the corresponding perimeter.
std::pair<size_t, size_t> find_previous_and_next_perimeter_point(const std::vector<SeamCandidate> &perimeter_points,
        size_t point_index) {
    const SeamCandidate &current = perimeter_points[point_index];
    int prev = point_index - 1; //for majority of points, it is true that neighbours lie behind and in front of them in the vector
    int next = point_index + 1;

    if (point_index == current.perimeter->start_index) {
        // if point_index is equal to start, it means that the previous neighbour is at the end
        prev = current.perimeter->end_index;
    }

    if (point_index == current.perimeter->end_index) {
        // if point_index is equal to end, than next neighbour is at the start
        next = current.perimeter->start_index;
    }

    assert(prev >= 0);
    assert(next >= 0);
    return {size_t(prev),size_t(next)};
}

// Computes all global model info - transforms object, performs raycasting,
// stores enforces and blockers
void compute_global_occlusion(GlobalModelInfo &result, const PrintObject *po) {
    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather occlusion meshes: start";
    auto obj_transform = po->trafo_centered();
    indexed_triangle_set triangle_set;
    indexed_triangle_set negative_volumes_set;
    //add all parts
    for (const ModelVolume *model_volume : po->model_object()->volumes) {
        if (model_volume->type() == ModelVolumeType::MODEL_PART
                || model_volume->type() == ModelVolumeType::NEGATIVE_VOLUME) {
            auto model_transformation = model_volume->get_matrix();
            indexed_triangle_set model_its = model_volume->mesh().its;
            its_transform(model_its, model_transformation);
            if (model_volume->type() == ModelVolumeType::MODEL_PART) {
                its_merge(triangle_set, model_its);
            } else {
                its_merge(negative_volumes_set, model_its);
            }
        }
    }
    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: gather occlusion meshes: end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: simplify occlusion meshes: start";

    //simplify raycasting mesh
    its_quadric_edge_collapse(triangle_set, SeamPlacer::raycasting_decimation_target_triangle_count, nullptr, nullptr,
            nullptr);
    triangle_set = its_subdivide(triangle_set, SeamPlacer::raycasting_subdivision_target_length);

    //simplify negative volumes
    its_quadric_edge_collapse(negative_volumes_set, SeamPlacer::raycasting_decimation_target_triangle_count, nullptr,
            nullptr,
            nullptr);
    negative_volumes_set = its_subdivide(negative_volumes_set, SeamPlacer::raycasting_subdivision_target_length);

    size_t negative_volumes_start_index = triangle_set.indices.size();
    its_merge(triangle_set, negative_volumes_set);
    its_transform(triangle_set, obj_transform);

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: simplify occlusion meshes: end";

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer:build AABB tree: start";
    auto raycasting_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(triangle_set.vertices,
            triangle_set.indices);

    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer:build AABB tree: end";
    result.model = triangle_set;
    result.model_tree = raycasting_tree;
    result.visiblity_info = raycast_visibility(raycasting_tree, triangle_set, negative_volumes_start_index);

#ifdef DEBUG_FILES
    auto filename = debug_out_path(("visiblity_of_" + std::to_string(po->id().id) + ".obj").c_str());
    result.debug_export(triangle_set, filename.c_str());
#endif
}

void gather_enforcers_blockers(GlobalModelInfo &result, const PrintObject *po) {
    BOOST_LOG_TRIVIAL(debug)
    << "SeamPlacer: build AABB trees for raycasting enforcers/blockers: start";

    auto obj_transform = po->trafo();

    for (const ModelVolume *mv : po->model_object()->volumes) {
        if (mv->is_seam_painted()) {
            auto model_transformation = mv->get_matrix();

            indexed_triangle_set enforcers = mv->seam_facets.get_facets(*mv, EnforcerBlockerType::ENFORCER);
            its_transform(enforcers, model_transformation);
            its_merge(result.enforcers, enforcers);

            indexed_triangle_set blockers = mv->seam_facets.get_facets(*mv, EnforcerBlockerType::BLOCKER);
            its_transform(blockers, model_transformation);
            its_merge(result.blockers, blockers);
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

struct SeamComparator {
    SeamPosition setup;

    SeamComparator(SeamPosition setup) :
            setup(setup) {
    }

    float compute_angle_penalty(float ccw_angle) const {
        // This function is used:
        // ((ℯ^(((1)/(x^(2)*3+1)))-1)/(ℯ-1))*1+((1)/(2+ℯ^(-x)))
        // looks scary, but it is gaussian combined with sigmoid,
        // so that concave points have much smaller penalty over convex ones

        return gauss(ccw_angle, 0.0f, 1.0f, 3.0f) +
                1.0f / (2 + std::exp(-ccw_angle)); // sigmoid, which heavily favourizes concave angles
    }

    // Standard comparator, must respect the requirements of comparators (e.g. give same result on same inputs) for sorting usage
    // should return if a is better seamCandidate than b
    bool is_first_better(const SeamCandidate &a, const SeamCandidate &b, const Vec2f &preffered_location = Vec2f { 0.0f,
            0.0f }) const {
        if (setup == SeamPosition::spAligned && a.central_enforcer != b.central_enforcer) {
            if (a.central_enforcer) {
                return true;
            }
            if (b.central_enforcer) {
                return false;
            }
        }

        // Blockers/Enforcers discrimination, top priority
        if (a.type > b.type) {
            return true;
        }
        if (b.type > a.type) {
            return false;
        }

        //avoid overhangs
        if (a.overhang > 0.0f && b.overhang < a.overhang) {
            return false;
        }

        // prefer hidden points (more than 1 mm inside)
        if (a.embedded_distance < -1.0f && b.embedded_distance > -1.0f) {
            return true;
        }
        if (b.embedded_distance < -1.0f && a.embedded_distance > -1.0f) {
            return false;
        }

        if (setup == SeamPosition::spRear) {
            return a.position.y() > b.position.y();
        }

        float distance_penalty_a = 1.0f;
        float distance_penalty_b = 1.0f;
        if (setup == spNearest) {
            distance_penalty_a = 1.1f - gauss((a.position.head<2>() - preffered_location).norm(), 0.0f, 1.0f, 0.005f);
            distance_penalty_b = 1.1f - gauss((a.position.head<2>() - preffered_location).norm(), 0.0f, 1.0f, 0.005f);
        }

        //ranges:          [0 - 1]                                              (0 - 1.3]                               [0.1 - 1.1)
        float penalty_a = (a.visibility + SeamPlacer::additional_angle_importance)
                * compute_angle_penalty(a.local_ccw_angle)
                * distance_penalty_a;
        float penalty_b = (b.visibility + SeamPlacer::additional_angle_importance)
                * compute_angle_penalty(b.local_ccw_angle)
                * distance_penalty_b;

        return penalty_a < penalty_b;
    }

    // Comparator used during alignment. If there is close potential aligned point, it is comapred to the current
    // sema point of the perimeter, to find out if the aligned point is not much worse than the current seam
    bool is_first_not_much_worse(const SeamCandidate &a, const SeamCandidate &b) const {
        // Blockers/Enforcers discrimination, top priority
        if (setup == SeamPosition::spAligned && a.central_enforcer != b.central_enforcer) {
            if (a.central_enforcer) {
                return true;
            }
            if (b.central_enforcer) {
                return false;
            }
        }

        if (a.type == EnforcedBlockedSeamPoint::Enforced) {
            return true;
        }

        if (a.type == EnforcedBlockedSeamPoint::Blocked) {
            return false;
        }

        if (a.type > b.type) {
            return true;
        }
        if (b.type > a.type) {
            return false;
        }

        //avoid overhangs
        if (a.overhang > 0.0f && b.overhang < a.overhang) {
            return false;
        }

        // prefer hidden points (more than 1 mm inside)
        if (a.embedded_distance < -1.0f && b.embedded_distance > -1.0f) {
            return true;
        }
        if (b.embedded_distance < -1.0f && a.embedded_distance > -1.0f) {
            return false;
        }

        if (setup == SeamPosition::spRandom) {
            return true;
        }

        if (setup == SeamPosition::spRear) {
            return a.position.y() > b.position.y();
        }

        //ranges:          [0 - 1]                                          (0 - 1.3]                  ;
        float penalty_a = (a.visibility + SeamPlacer::additional_angle_importance)
                * compute_angle_penalty(a.local_ccw_angle);
        float penalty_b = (b.visibility + SeamPlacer::additional_angle_importance)
                * compute_angle_penalty(b.local_ccw_angle);

        return penalty_a <= penalty_b || std::abs(penalty_a - penalty_b) < SeamPlacer::seam_align_score_tolerance;
    }

    //always nonzero, positive
    float get_penalty(const SeamCandidate &a) const {
        if (setup == SeamPosition::spRear) {
            return a.position.y();
        }

        return (a.visibility + SeamPlacer::additional_angle_importance) * compute_angle_penalty(a.local_ccw_angle);
    }
}
;

#ifdef DEBUG_FILES
void debug_export_points(const std::vector<std::vector<SeamPlacerImpl::SeamCandidate>> &object_perimter_points,
        const BoundingBox &bounding_box, std::string object_name, const SeamComparator &comparator) {
    for (size_t layer_idx = 0; layer_idx < object_perimter_points.size(); ++layer_idx) {
        std::string angles_file_name = debug_out_path(
                (object_name + "_angles_" + std::to_string(layer_idx) + ".svg").c_str());
        SVG angles_svg {angles_file_name, bounding_box};
        float min_vis = 0;
        float max_vis = min_vis;

        float min_weight = std::numeric_limits<float>::min();
        float max_weight = min_weight;

        for (const SeamCandidate &point : object_perimter_points[layer_idx]) {
            Vec3i color = value_rgbi(-PI, PI, point.local_ccw_angle);
            std::string fill = "rgb(" + std::to_string(color.x()) + "," + std::to_string(color.y()) + ","
            + std::to_string(color.z()) + ")";
            angles_svg.draw(scaled(Vec2f(point.position.head<2>())), fill);
            min_vis = std::min(min_vis, point.visibility);
            max_vis = std::max(max_vis, point.visibility);

            min_weight = std::min(min_weight, -comparator.get_penalty(point));
            max_weight = std::max(max_weight, -comparator.get_penalty(point));

        }

        std::string visiblity_file_name = debug_out_path(
                (object_name + "_visibility_" + std::to_string(layer_idx) + ".svg").c_str());
        SVG visibility_svg {visiblity_file_name, bounding_box};
        std::string weights_file_name = debug_out_path(
                (object_name + "_weight_" + std::to_string(layer_idx) + ".svg").c_str());
        SVG weight_svg {weights_file_name, bounding_box};
        std::string overhangs_file_name = debug_out_path(
                (object_name + "_overhang_" + std::to_string(layer_idx) + ".svg").c_str());
        SVG overhangs_svg {overhangs_file_name, bounding_box};

        for (const SeamCandidate &point : object_perimter_points[layer_idx]) {
            Vec3i color = value_rgbi(min_vis, max_vis, point.visibility);
            std::string visibility_fill = "rgb(" + std::to_string(color.x()) + "," + std::to_string(color.y()) + ","
            + std::to_string(color.z()) + ")";
            visibility_svg.draw(scaled(Vec2f(point.position.head<2>())), visibility_fill);

            Vec3i weight_color = value_rgbi(min_weight, max_weight, -comparator.get_penalty(point));
            std::string weight_fill = "rgb(" + std::to_string(weight_color.x()) + "," + std::to_string(weight_color.y())
            + ","
            + std::to_string(weight_color.z()) + ")";
            weight_svg.draw(scaled(Vec2f(point.position.head<2>())), weight_fill);

            Vec3i overhang_color = value_rgbi(-0.5, 0.5, std::clamp(point.overhang, -0.5f, 0.5f));
            std::string overhang_fill = "rgb(" + std::to_string(overhang_color.x()) + ","
            + std::to_string(overhang_color.y())
            + ","
            + std::to_string(overhang_color.z()) + ")";
            overhangs_svg.draw(scaled(Vec2f(point.position.head<2>())), overhang_fill);
        }
    }
}
#endif

// Pick best seam point based on the given comparator
void pick_seam_point(std::vector<SeamCandidate> &perimeter_points, size_t start_index,
        const SeamComparator &comparator) {
    size_t end_index = perimeter_points[start_index].perimeter->end_index;

    size_t seam_index = start_index;
    for (size_t index = start_index; index <= end_index; ++index) {
        if (comparator.is_first_better(perimeter_points[index], perimeter_points[seam_index])) {
            seam_index = index;
        }
    }
    perimeter_points[start_index].perimeter->seam_index = seam_index;
}

size_t pick_nearest_seam_point_index(const std::vector<SeamCandidate> &perimeter_points, size_t start_index,
        const Vec2f &preffered_location) {
    size_t end_index = perimeter_points[start_index].perimeter->end_index;
    SeamComparator comparator { spNearest };

    size_t seam_index = start_index;
    for (size_t index = start_index; index <= end_index; ++index) {
        if (comparator.is_first_better(perimeter_points[index], perimeter_points[seam_index], preffered_location)) {
            seam_index = index;
        }
    }
    return seam_index;
}

// picks random seam point uniformly, respecting enforcers blockers and overhang avoidance.
void pick_random_seam_point(std::vector<SeamCandidate> &perimeter_points, size_t start_index) {
    SeamComparator comparator { spRandom };

    // algorithm keeps a list of viable points and their lengths. If it finds a point
    // that is much better than the viable_example_index (e.g. better type, no overhang; see is_first_not_much_worse)
    // then it throws away stored lists and starts from start
    // in the end, the list should contain points with same type (Enforced > Neutral > Blocked) and also only those which are not
    // big overhang.
    size_t viable_example_index = start_index;
    size_t end_index = perimeter_points[start_index].perimeter->end_index;
    std::vector<size_t> viable_indices;
    std::vector<float> viable_edges_lengths;
    std::vector<Vec3f> viable_edges;

    for (size_t index = start_index; index <= end_index; ++index) {
        if (comparator.is_first_not_much_worse(perimeter_points[index], perimeter_points[viable_example_index]) &&
                comparator.is_first_not_much_worse(perimeter_points[viable_example_index], perimeter_points[index])) {
            // index ok, push info into respective vectors
            Vec3f edge_to_next;
            if (index == end_index) {
                edge_to_next = (perimeter_points[start_index].position - perimeter_points[index].position);
            } else
            {
                edge_to_next = (perimeter_points[index + 1].position - perimeter_points[index].position);
            }
            float dist_to_next = edge_to_next.norm();
            viable_indices.push_back(index);
            viable_edges_lengths.push_back(dist_to_next);
            viable_edges.push_back(edge_to_next);
        } else if (comparator.is_first_not_much_worse(perimeter_points[viable_example_index],
                perimeter_points[index])) {
            // index is worse then viable_example_index, skip this point
        } else {
            // index is better than viable example index, update example, clear gathered info, start again
            // clear up all gathered info, start from scratch, update example index
            viable_example_index = index;
            viable_indices.clear();
            viable_edges_lengths.clear();
            viable_edges.clear();

            Vec3f edge_to_next;
            if (index == end_index) {
                edge_to_next = (perimeter_points[start_index].position - perimeter_points[index].position);
            } else {
                edge_to_next = (perimeter_points[index + 1].position - perimeter_points[index].position);
            }
            float dist_to_next = edge_to_next.norm();
            viable_indices.push_back(index);
            viable_edges_lengths.push_back(dist_to_next);
            viable_edges.push_back(edge_to_next);
        }
    }

    // now pick random point from the stored options
    float len_sum = std::accumulate(viable_edges_lengths.begin(), viable_edges_lengths.end(), 0.0f);
    float picked_len = len_sum * (rand() / (float(RAND_MAX) + 1));

    size_t point_idx = 0;
    while (picked_len - viable_edges_lengths[point_idx] > 0) {
        picked_len = picked_len - viable_edges_lengths[point_idx];
        point_idx++;
    }

    Perimeter *perimeter = perimeter_points[start_index].perimeter.get();
    perimeter->seam_index = viable_indices[point_idx];
    perimeter->final_seam_position = perimeter_points[perimeter->seam_index].position
            + viable_edges[point_idx].normalized() * picked_len;
    perimeter->finalized = true;

}

struct EdgeGridWrapper {
    explicit EdgeGridWrapper(ExPolygons ex_polys) :
            ex_polys(ex_polys) {

        grid.create(this->ex_polys, distance_field_resolution);
        grid.calculate_sdf();
    }
    const coord_t distance_field_resolution = coord_t(scale_(1.) + 0.5);
    EdgeGrid::Grid grid;
    ExPolygons ex_polys;
}
;

EdgeGridWrapper compute_layer_merged_edge_grid(const Layer *layer) {
    static const float eps = float(scale_(layer->object()->config().slice_closing_radius.value));
    // merge with offset
    ExPolygons merged = layer->merged(eps);
    // ofsset back
    ExPolygons layer_outline = offset_ex(merged, -eps);
    return EdgeGridWrapper(layer_outline);
}

} // namespace SeamPlacerImpl

// Parallel process and extract each perimeter polygon of the given print object.
// Gather SeamCandidates of each layer into vector and build KDtree over them
// Store results in the SeamPlacer varaibles m_perimeter_points_per_object and m_perimeter_points_trees_per_object
void SeamPlacer::gather_seam_candidates(const PrintObject *po,
        const SeamPlacerImpl::GlobalModelInfo &global_model_info, const SeamPosition configured_seam_preference) {
    using namespace SeamPlacerImpl;

    m_perimeter_points_per_object.emplace(po, po->layer_count());
    m_perimeter_points_trees_per_object.emplace(po, po->layer_count());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, po->layers().size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    std::vector<SeamCandidate> &layer_candidates =
                            m_perimeter_points_per_object[po][layer_idx];
                    const Layer *layer = po->get_layer(layer_idx);
                    auto unscaled_z = layer->slice_z;
                    std::vector<const LayerRegion*> regions;
                    //NOTE corresponding region ptr may be null, if the layer has zero perimeters
                    Polygons polygons = extract_perimeter_polygons(layer, configured_seam_preference, regions);
                    for (size_t poly_index = 0; poly_index < polygons.size(); ++poly_index) {
                        process_perimeter_polygon(polygons[poly_index], unscaled_z,
                                regions[poly_index], global_model_info, layer_candidates);
                    }
                    auto functor = SeamCandidateCoordinateFunctor { &layer_candidates };
                    m_perimeter_points_trees_per_object[po][layer_idx] =
                            std::make_unique<SeamCandidatesTree>(functor, layer_candidates.size());
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

void SeamPlacer::calculate_overhangs_and_layer_embedding(const PrintObject *po) {
    using namespace SeamPlacerImpl;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, m_perimeter_points_per_object[po].size()),
            [&](tbb::blocked_range<size_t> r) {
                std::unique_ptr<EdgeGridWrapper> prev_layer_grid;
                if (r.begin() > 0) { // previous layer exists
                    prev_layer_grid = std::make_unique<EdgeGridWrapper>(
                            compute_layer_merged_edge_grid(po->layers()[r.begin() - 1]));
                }

                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    bool layer_has_multiple_loops =
                            m_perimeter_points_per_object[po][layer_idx][0].perimeter->end_index
                                    < m_perimeter_points_per_object[po][layer_idx].size() - 1;
                    std::unique_ptr<EdgeGridWrapper> current_layer_grid = std::make_unique<EdgeGridWrapper>(
                            compute_layer_merged_edge_grid(po->layers()[layer_idx]));

                    for (SeamCandidate &perimeter_point : m_perimeter_points_per_object[po][layer_idx]) {
                        Point point = Point::new_scale(Vec2f { perimeter_point.position.head<2>() });
                        if (prev_layer_grid.get() != nullptr) {
                            coordf_t overhang_dist;
                            prev_layer_grid->grid.signed_distance(point, scaled(perimeter_point.perimeter->flow_width),
                                    overhang_dist);
                            perimeter_point.overhang =
                                    unscale<float>(overhang_dist) - perimeter_point.perimeter->flow_width;
                        }

                        if (layer_has_multiple_loops) { // search for embedded perimeter points (points hidden inside the print ,e.g. multimaterial join, best position for seam)
                            coordf_t layer_embedded_distance;
                            current_layer_grid->grid.signed_distance(point, scaled(1.0f),
                                    layer_embedded_distance);
                            perimeter_point.embedded_distance = unscale<float>(layer_embedded_distance);
                        }
                    }

                    prev_layer_grid.swap(current_layer_grid);
                }
            }
            );
        }

// Estimates, if there is good seam point in the layer_idx which is close to last_point_pos
// uses comparator.is_first_not_much_worse method to compare current seam with the closest point
// (if current seam is too far away )
// If the current chosen stream is close enough, it is stored in seam_string. returns true and updates last_point_pos
// If the closest point is good enough to replace current chosen seam, it is stored in potential_string_seams, returns true and updates last_point_pos
// Otherwise does nothing, returns false
// sadly cannot be const because map access operator[] is not const, since it can create new object
bool SeamPlacer::find_next_seam_in_layer(const PrintObject *po,
        std::pair<size_t, size_t> &last_point_indexes,
        size_t layer_idx, const SeamPlacerImpl::SeamComparator &comparator,
        std::vector<std::pair<size_t, size_t>> &seam_string) {
    using namespace SeamPlacerImpl;

    const SeamCandidate &last_point =
            m_perimeter_points_per_object[po][last_point_indexes.first][last_point_indexes.second];

    Vec3f projected_position { last_point.position.x(), last_point.position.y(), float(
            po->get_layer(layer_idx)->slice_z) };
    //find closest point in next layer
    size_t closest_point_index = find_closest_point(
            *m_perimeter_points_trees_per_object[po][layer_idx], projected_position);

    SeamCandidate &closest_point = m_perimeter_points_per_object[po][layer_idx][closest_point_index];

    if (closest_point.perimeter->finalized) { //already finalized, skip
        return false;
    }

    //from the closest point, deduce index of seam in the next layer
    SeamCandidate &next_layer_seam =
            m_perimeter_points_per_object[po][layer_idx][closest_point.perimeter->seam_index];

    if (next_layer_seam.central_enforcer
            && (next_layer_seam.position - projected_position).norm() < 3 * SeamPlacer::seam_align_tolerable_dist) {
        seam_string.push_back( { layer_idx, closest_point.perimeter->seam_index });
        last_point_indexes = std::pair<size_t, size_t> { layer_idx, closest_point.perimeter->seam_index };
        return true;
    }

    auto are_similar = [&](const SeamCandidate &a, const SeamCandidate &b) {
        return comparator.is_first_not_much_worse(a, b) && comparator.is_first_not_much_worse(b, a);
    };

    if ((closest_point.position - projected_position).norm() < SeamPlacer::seam_align_tolerable_dist
            && comparator.is_first_not_much_worse(closest_point, next_layer_seam)
            && are_similar(last_point, closest_point)) {
        seam_string.push_back( { layer_idx, closest_point_index });
        last_point_indexes = std::pair<size_t, size_t> { layer_idx, closest_point_index };
        return true;
    } else {
        return false;
    }

}

// clusters already chosen seam points into strings across multiple layers, and then
// aligns the strings via polynomial fit
// Does not change the positions of the SeamCandidates themselves, instead stores
// the new aligned position into the shared Perimeter structure of each perimeter
// Note that this position does not necesarilly lay on the perimeter.
void SeamPlacer::align_seam_points(const PrintObject *po, const SeamPlacerImpl::SeamComparator &comparator) {
    using namespace SeamPlacerImpl;

    // Prepares Debug files for writing.
#ifdef DEBUG_FILES
    Slic3r::CNumericLocalesSetter locales_setter;
    auto clusters_f = debug_out_path(("seam_clusters_of_" + std::to_string(po->id().id) + ".obj").c_str());
    FILE *clusters = boost::nowide::fopen(clusters_f.c_str(), "w");
    if (clusters == nullptr) {
        BOOST_LOG_TRIVIAL(error)
        << "stl_write_obj: Couldn't open " << clusters_f << " for writing";
        return;
    }
    auto aligned_f = debug_out_path(("aligned_clusters_of_" + std::to_string(po->id().id) + ".obj").c_str());
    FILE *aligns = boost::nowide::fopen(aligned_f.c_str(), "w");
    if (aligns == nullptr) {
        BOOST_LOG_TRIVIAL(error)
        << "stl_write_obj: Couldn't open " << clusters_f << " for writing";
        return;
    }
#endif

    //gather vector of all seams on the print_object - pair of layer_index and seam__index within that layer
    std::vector<std::pair<size_t, size_t>> seams;
    for (size_t layer_idx = 0; layer_idx < m_perimeter_points_per_object[po].size(); ++layer_idx) {
        std::vector<SeamCandidate> &layer_perimeter_points =
                m_perimeter_points_per_object[po][layer_idx];
        size_t current_point_index = 0;
        while (current_point_index < layer_perimeter_points.size()) {
            seams.emplace_back(layer_idx, layer_perimeter_points[current_point_index].perimeter->seam_index);
            current_point_index = layer_perimeter_points[current_point_index].perimeter->end_index + 1;
        }
    }

    //sort them before alignment. Alignment is sensitive to initializaion, this gives it better chance to choose something nice
    std::sort(seams.begin(), seams.end(),
            [&](const std::pair<size_t, size_t> &left, const std::pair<size_t, size_t> &right) {
                return comparator.is_first_better(m_perimeter_points_per_object[po][left.first][left.second],
                        m_perimeter_points_per_object[po][right.first][right.second]);
            }
    );

    //align the seam points - start with the best, and check if they are aligned, if yes, skip, else start alignment
    for (const std::pair<size_t, size_t> &seam : seams) {
        size_t layer_idx = seam.first;
        size_t seam_index = seam.second;
        std::vector<SeamCandidate> &layer_perimeter_points =
                m_perimeter_points_per_object[po][layer_idx];
        if (layer_perimeter_points[seam_index].perimeter->finalized) {
            // This perimeter is already aligned, skip seam
            continue;
        } else {

            //initialize searching for seam string - cluster of nearby seams on previous and next layers
            int skips = SeamPlacer::seam_align_tolerable_skips / 2;
            int next_layer = layer_idx + 1;
            std::pair<size_t, size_t> last_point_indexes = std::pair<size_t, size_t>(layer_idx, seam_index);

            std::vector<std::pair<size_t, size_t>> seam_string { std::pair<size_t, size_t>(layer_idx, seam_index) };

            //find seams or potential seams in forward direction; there is a budget of skips allowed
            while (skips >= 0 && next_layer < int(m_perimeter_points_per_object[po].size())) {
                if (find_next_seam_in_layer(po, last_point_indexes, next_layer, comparator, seam_string)) {
                    //String added, last_point_pos updated, nothing to be done
                } else {
                    // Layer skipped, reduce number of available skips
                    skips--;
                }
                next_layer++;
            }

            //do additional check in back direction
            next_layer = layer_idx - 1;
            skips = SeamPlacer::seam_align_tolerable_skips / 2;
            last_point_indexes = std::pair<size_t, size_t>(layer_idx, seam_index);
            while (skips >= 0 && next_layer >= 0) {
                if (find_next_seam_in_layer(po, last_point_indexes, next_layer, comparator, seam_string)) {
                    //String added, last_point_pos updated, nothing to be done
                } else {
                    // Layer skipped, reduce number of available skips
                    skips--;
                }
                next_layer--;
            }

            if (seam_string.size() < seam_align_minimum_string_seams) {
                //string NOT long enough to be worth aligning, skip
                continue;
            }

            // String is long engouh, all string seams and potential string seams gathered, now do the alignment
            //sort by layer index
            std::sort(seam_string.begin(), seam_string.end(),
                    [](const std::pair<size_t, size_t> &left, const std::pair<size_t, size_t> &right) {
                        return left.first < right.first;
                    });

            // gather all positions of seams and their weights (weights are derived as negative penalty, they are made positive in next step)
            std::vector<Vec2f> observations(seam_string.size());
            std::vector<float> observation_points(seam_string.size());
            std::vector<float> weights(seam_string.size());

            //init min_weight by the first point
            float min_weight = -comparator.get_penalty(
                    m_perimeter_points_per_object[po][seam_string[0].first][seam_string[0].second]);

            //gather points positions and weights, update min_weight in each step
            for (size_t index = 0; index < seam_string.size(); ++index) {
                Vec3f pos =
                        m_perimeter_points_per_object[po][seam_string[index].first][seam_string[index].second].position;
                observations[index] = pos.head<2>();
                observation_points[index] = pos.z();
                weights[index] = -comparator.get_penalty(
                        m_perimeter_points_per_object[po][seam_string[index].first][seam_string[index].second]);
                min_weight = std::min(min_weight, weights[index]);
            }

            //makes all weights positive
            for (float &w : weights) {
                w = w - min_weight + 0.01;
            }

            // Curve Fitting
            size_t number_of_segments = std::max(size_t(1),
                    size_t(observations.size() / SeamPlacer::seam_align_seams_per_segment));
            auto curve = Geometry::fit_cubic_bspline(observations, observation_points, weights, number_of_segments);

            // Do alignment - compute fitted point for each point in the string from its Z coord, and store the position into
            // Perimeter structure of the point; also set flag aligned to true
            for (const auto &pair : seam_string) {
                Vec3f current_pos = m_perimeter_points_per_object[po][pair.first][pair.second].position;
                Vec2f fitted_pos = curve.get_fitted_value(current_pos.z());

                Perimeter *perimeter =
                        m_perimeter_points_per_object[po][pair.first][pair.second].perimeter.get();
                perimeter->final_seam_position = Vec3f { fitted_pos.x(), fitted_pos.y(), current_pos.z() };
                perimeter->finalized = true;
            }

#ifdef DEBUG_FILES
            auto randf = []() {
                return float(rand()) / float(RAND_MAX);
            };
            Vec3f color { randf(), randf(), randf() };
            for (size_t i = 0; i < seam_string.size(); ++i) {
                auto orig_seam = m_perimeter_points_per_object[po][seam_string[i].first][seam_string[i].second];
                fprintf(clusters, "v %f %f %f %f %f %f \n", orig_seam.position[0],
                        orig_seam.position[1],
                        orig_seam.position[2], color[0], color[1],
                        color[2]);
            }

            color = Vec3f { randf(), randf(), randf() };
            for (size_t i = 0; i < seam_string.size(); ++i) {
                Perimeter *perimeter =
                        m_perimeter_points_per_object[po][seam_string[i].first][seam_string[i].second].perimeter.get();
                fprintf(aligns, "v %f %f %f %f %f %f \n", perimeter->final_seam_position[0],
                        perimeter->final_seam_position[1],
                        perimeter->final_seam_position[2], color[0], color[1],
                        color[2]);
            }
#endif
        }
    }

#ifdef DEBUG_FILES
    fclose(clusters);
    fclose(aligns);
#endif

}

void SeamPlacer::init(const Print &print) {
    using namespace SeamPlacerImpl;
    m_perimeter_points_trees_per_object.clear();
    m_perimeter_points_per_object.clear();

    for (const PrintObject *po : print.objects()) {

        SeamPosition configured_seam_preference = po->config().seam_position.value;
        SeamComparator comparator { configured_seam_preference };

        GlobalModelInfo global_model_info { };
        gather_enforcers_blockers(global_model_info, po);

        if (configured_seam_preference == spAligned || configured_seam_preference == spNearest) {
            compute_global_occlusion(global_model_info, po);
        }

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: gather_seam_candidates: start";
        gather_seam_candidates(po, global_model_info, configured_seam_preference);
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: gather_seam_candidates: end";

        if (configured_seam_preference == spAligned || configured_seam_preference == spNearest) {
            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: calculate_candidates_visibility : start";
            calculate_candidates_visibility(po, global_model_info);
            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: calculate_candidates_visibility : end";
        }

        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_overhangs and layer embdedding : start";
        calculate_overhangs_and_layer_embedding(po);
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: calculate_overhangs and layer embdedding: end";

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
                            if (configured_seam_preference == spRandom) {
                                pick_random_seam_point(layer_perimeter_points, current);
                            } else {
                                pick_seam_point(layer_perimeter_points, current, comparator);
                            }
                            current = layer_perimeter_points[current].perimeter->end_index + 1;
                        }
                    }
                });
        BOOST_LOG_TRIVIAL(debug)
        << "SeamPlacer: pick_seam_point : end";

        if (configured_seam_preference == spAligned) {
            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: align_seam_points : start";
            align_seam_points(po, comparator);
            BOOST_LOG_TRIVIAL(debug)
            << "SeamPlacer: align_seam_points : end";
        }

#ifdef DEBUG_FILES
        debug_export_points(m_perimeter_points_per_object[po], po->bounding_box(), std::to_string(po->id().id),
                comparator);
#endif
    }
}

void SeamPlacer::place_seam(const Layer *layer, ExtrusionLoop &loop, bool external_first,
        const Point &last_pos) const {
    using namespace SeamPlacerImpl;
    const PrintObject *po = layer->object();
//NOTE this is necessary, since layer->id() is quite unreliable
    size_t layer_index = std::max(0, int(layer->id()) - int(po->slicing_parameters().raft_layers()));
    double unscaled_z = layer->slice_z;

    const auto &perimeter_points_tree = *m_perimeter_points_trees_per_object.find(po)->second[layer_index];
    const auto &perimeter_points = m_perimeter_points_per_object.find(po)->second[layer_index];

    const Point &fp = loop.first_point();

    Vec2f unscaled_p = unscale(fp).cast<float>();
    size_t closest_perimeter_point_index = find_closest_point(perimeter_points_tree,
            Vec3f { unscaled_p.x(), unscaled_p.y(), float(unscaled_z) });
    const Perimeter *perimeter = perimeter_points[closest_perimeter_point_index].perimeter.get();

    size_t seam_index;
    if (po->config().seam_position == spNearest) {
        seam_index = pick_nearest_seam_point_index(perimeter_points, perimeter->start_index,
                unscale(last_pos).cast<float>());
    } else {
        seam_index = perimeter->seam_index;
    }

    Vec3f seam_position = perimeter_points[seam_index].position;
    if (perimeter->finalized) {
        seam_position = perimeter->final_seam_position;
    }
    Point seam_point = scaled(Vec2d { seam_position.x(), seam_position.y() });

    if (!loop.split_at_vertex(seam_point))
// The point is not in the original loop.
// Insert it.
        loop.split_at(seam_point, true);
}

} // namespace Slic3r

