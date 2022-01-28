#include "SeamPlacerNG.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <random>
#include <algorithm>
#include <atomic>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/SVG.hpp"
#include "libslic3r/Layer.hpp"

//TODO remove
#include <boost/nowide/cstdio.hpp>

namespace Slic3r {

namespace SeamPlacerImpl {

/// Coordinate frame
class Frame
{
public:
    Frame()
    {
        mX = Vec3d(1, 0, 0);
        mY = Vec3d(0, 1, 0);
        mZ = Vec3d(0, 0, 1);
    }

    Frame(const Vec3d &x, const Vec3d &y, const Vec3d &z)
        : mX(x), mY(y), mZ(z)
    {}

    void set_from_z(const Vec3d &z)
    {
        mZ         = z.normalized();
        Vec3d tmpZ = mZ;
        Vec3d tmpX = (std::abs(tmpZ.x()) > 0.99f) ? Vec3d(0, 1, 0) :
                                                    Vec3d(1, 0, 0);
        mY         = (tmpZ.cross(tmpX)).normalized();
        mX         = mY.cross(tmpZ);
    }

    Vec3d to_world(const Vec3d &a) const
    {
        return a.x() * mX + a.y() * mY + a.z() * mZ;
    }

    Vec3d to_local(const Vec3d &a) const
    {
        return Vec3d(mX.dot(a), mY.dot(a), mZ.dot(a));
    }

    const Vec3d &binormal() const { return mX; }

    const Vec3d &tangent() const { return mY; }

    const Vec3d &normal() const { return mZ; }

private:
    Vec3d mX, mY, mZ;
};

Vec3d sample_sphere_uniform(const Vec2f &samples)
{
    float term_one = 2.0f * M_PIf32 * samples.x();
    float term_two = 2.0f * sqrt(samples.y() - samples.y() * samples.y());
    return {cos(term_one) * term_two, sin(term_one) * term_two,
            1.0f - 2.0f * samples.y()};
}

Vec3d sample_power_cosine_hemisphere(const Vec2f &samples, float power)
{
    const float term1 = 2.f * M_PIf32 * samples.x();
    const float term2 = pow(samples.y(), 1.f / (power + 1.f));
    const float term3 = sqrt(1.f - term2 * term2);

    return Vec3d(cos(term1) * term3, sin(term1) * term3, term2);
}

void raycast_visibility(
    size_t                                  ray_count,
    const AABBTreeIndirect::Tree<3, float> &raycasting_tree,
    const indexed_triangle_set             &triangles,
    const KDTreeIndirect<3, coordf_t, KDTreeCoordinateFunctor>
                               &perimeter_points_tree,
    std::vector<SeamCandidate> &perimeter_points)
{
    auto  bbox                 = raycasting_tree.node(0).bbox;
    Vec3d vision_sphere_center = bbox.center().cast<double>();
    float vision_sphere_raidus = (bbox.sizes().maxCoeff() *
                                  0.55); // 0.5 (half) covers whole object,
                                         // 0.05 added to avoid corner cases

    // Prepare random samples per ray
    std::random_device                    rnd_device;
    std::mt19937                          mersenne_engine{rnd_device()};
    std::uniform_real_distribution<float> dist{0, 1};

    auto gen = [&dist, &mersenne_engine]() {
        return Vec2f(dist(mersenne_engine), dist(mersenne_engine));
    };

    BOOST_LOG_TRIVIAL(debug) << "PM: generate random samples: start";
    std::vector<Vec2f> global_dir_random_samples(ray_count);
    generate(begin(global_dir_random_samples), end(global_dir_random_samples),
             gen);
    std::vector<Vec2f> local_dir_random_samples(ray_count);
    generate(begin(local_dir_random_samples), end(local_dir_random_samples),
             gen);

    BOOST_LOG_TRIVIAL(debug) << "PM: generate random samples: end";

    std::vector<std::atomic<size_t>> visibility_counters(
        perimeter_points.size());

    BOOST_LOG_TRIVIAL(debug)
        << "PM: raycast visibility for " << ray_count << " rays: start";
    // raycast visibility
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, ray_count),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t index = r.begin(); index < r.end(); ++index) {
                Vec3d global_ray_dir = sample_sphere_uniform(
                    global_dir_random_samples[index]);
                Vec3d ray_origin = (vision_sphere_center +
                                    global_ray_dir * vision_sphere_raidus);
                Vec3d local_dir  = sample_power_cosine_hemisphere(
                     local_dir_random_samples[index], 1.0);

                Frame f;
                f.set_from_z(global_ray_dir);
                Vec3d final_ray_dir = (f.to_world(local_dir));

                igl::Hit hitpoint;
                // FIXME: This query will not compile for float ray origin and
                // direction for some reason
                auto hit = AABBTreeIndirect::intersect_ray_first_hit(
                    triangles.vertices, triangles.indices, raycasting_tree,
                    ray_origin, final_ray_dir, hitpoint);

                if (hit) {
                    auto face  = triangles.indices[hitpoint.id];
                    auto edge1 = triangles.vertices[face[1]] -
                                 triangles.vertices[face[0]];
                    auto edge2 = triangles.vertices[face[2]] -
                                 triangles.vertices[face[0]];

                    Vec3d hit_pos = (triangles.vertices[face[0]] +
                                     edge1 * hitpoint.u + edge2 * hitpoint.v)
                                        .cast<double>();

                    auto perimeter_point_index =
                        find_closest_point(perimeter_points_tree, hit_pos);

                    visibility_counters[perimeter_point_index]
                        .fetch_add(1, std::memory_order_relaxed);
                }
            }
        });

    BOOST_LOG_TRIVIAL(debug)
        << "PM: raycast visibility for " << ray_count << " rays: end";

    BOOST_LOG_TRIVIAL(debug)
        << "PM: assign visibility to perimenter points : start";

    tbb::parallel_for(tbb::blocked_range<size_t>(0, perimeter_points.size()),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t index = r.begin(); index < r.end();
                               ++index) {
                              perimeter_points[index].visibility =
                                  visibility_counters[index];
                          }
                      });
    BOOST_LOG_TRIVIAL(debug)
        << "PM: assign visibility to perimenter points : end";

    its_write_obj(triangles, "triangles.obj");

    Slic3r::CNumericLocalesSetter locales_setter;
    FILE *fp = boost::nowide::fopen("perimeter.obj", "w");
    if (fp == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Couldn't open "
                                 << "perimeter.obj"
                                 << " for writing";
    }

    for (size_t i = 0; i < perimeter_points.size(); ++i)
        fprintf(fp, "v %f %f %f %zu\n", perimeter_points[i].position[0],
                perimeter_points[i].position[1],
                perimeter_points[i].position[2],
                perimeter_points[i].visibility);
    fclose(fp);
}

} // namespace SeamPlacerImpl

void SeamPlacer::init(const Print &print)
{
    using namespace SeamPlacerImpl;
    seam_candidates_trees.clear();

    for (const PrintObject *po : print.objects()) {
        BOOST_LOG_TRIVIAL(debug)
            << "PM: build AABB tree for raycasting: start";
        // Build AABB tree for raycasting
        auto obj_transform = po->trafo_centered();
        auto triangle_set  = po->model_object()->raw_indexed_triangle_set();
        its_transform(triangle_set, obj_transform);

        auto raycasting_tree =
            AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
                triangle_set.vertices, triangle_set.indices);

        BOOST_LOG_TRIVIAL(debug) << "PM: build AABB tree for raycasting: end";

        BOOST_LOG_TRIVIAL(debug)
            << "PM: gather and build KD tree with seam candidates: start";
        // gather seam candidates (perimeter points)
        auto seam_candidates = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0, po->layers().size()),
            std::vector<SeamCandidate>{},
            [&](tbb::blocked_range<size_t> r,
                std::vector<SeamCandidate> init) {
                for (size_t index = r.begin(); index < r.end(); ++index) {
                    const auto layer      = po->layers()[index];
                    auto       unscaled_z = layer->slice_z;
                    for (Points points : layer->lslices) {
                        for (Point point : points) {
                            auto unscaled_p        = unscale(point);
                            auto unscaled_position = Vec3d{unscaled_p.x(),
                                                           unscaled_p.y(),
                                                           unscaled_z};
                            init.emplace_back(unscaled_position);
                        }
                    }
                }
                return init;
            },
            [](std::vector<SeamCandidate> prev,
               std::vector<SeamCandidate> next) {
                prev.insert(prev.end(), next.begin(), next.end());
                return prev;
            });
        // Build KD tree with seam candidates
        auto functor = KDTreeCoordinateFunctor{&seam_candidates};
        auto perimeter_points_tree = KDTreeIndirect<
            3, coordf_t, KDTreeCoordinateFunctor>{functor,
                                                  seam_candidates.size()};

        BOOST_LOG_TRIVIAL(debug)
            << "PM: gather and build KD tree with seam candidates: end";

        raycast_visibility(1000000, raycasting_tree, triangle_set,
                           perimeter_points_tree, seam_candidates);
    }
}
} // namespace Slic3r
