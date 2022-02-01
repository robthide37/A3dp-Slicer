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
    float term_one = 2.0f * M_PIf32 * samples.x();
    float term_two = 2.0f * sqrt(samples.y() - samples.y() * samples.y());
    return {cos(term_one) * term_two, sin(term_one) * term_two,
        1.0f - 2.0f * samples.y()};
}

Vec3d sample_power_cosine_hemisphere(const Vec2f &samples, float power) {
    const float term1 = 2.f * M_PIf32 * samples.x();
    const float term2 = pow(samples.y(), 1.f / (power + 1.f));
    const float term3 = sqrt(1.f - term2 * term2);

    return Vec3d(cos(term1) * term3, sin(term1) * term3, term2);
}

void resolve_geometry_hit(const igl::Hit &hitpoint,
        const indexed_triangle_set &triangles,
        const KDTreeIndirect<3, coordf_t, KDTreeCoordinateFunctor> &perimeter_points_tree,
        const std::vector<SeamCandidate> &perimeter_points,
        std::vector<std::atomic<size_t>> &visibility_counters) {
    auto face = triangles.indices[hitpoint.id];
    auto edge1 = triangles.vertices[face[1]] - triangles.vertices[face[0]];
    auto edge2 = triangles.vertices[face[2]] - triangles.vertices[face[0]];

    Vec3d hit_pos = (triangles.vertices[face[0]] + edge1 * hitpoint.u + edge2 * hitpoint.v).cast<double>();
    auto perimeter_point_index = find_closest_point(perimeter_points_tree, hit_pos);
    auto perimeter_point_pos = perimeter_points[perimeter_point_index].position;
    auto dist_squared = (perimeter_point_pos - hit_pos).squaredNorm();

    visibility_counters[perimeter_point_index].fetch_add(1.0 / dist_squared, std::memory_order_relaxed);
}

std::vector<HitInfo> raycast_visibility(size_t ray_count,
        const AABBTreeIndirect::Tree<3, float> &raycasting_tree,
        const indexed_triangle_set &triangles) {
    auto bbox = raycasting_tree.node(0).bbox;
    Vec3d vision_sphere_center = bbox.center().cast<double>();
    float vision_sphere_raidus = (bbox.sizes().maxCoeff() * 0.55); // 0.5 (half) covers whole object,
                                                                   // 0.05 added to avoid corner cases

    // Prepare random samples per ray
    std::random_device rnd_device;
    std::mt19937 mersenne_engine { rnd_device() };
    std::uniform_real_distribution<float> dist { 0, 1 };

    auto gen = [&dist, &mersenne_engine]() {
        return Vec2f(dist(mersenne_engine), dist(mersenne_engine));
    };

    BOOST_LOG_TRIVIAL(debug)
    << "PM: generate random samples: start";
    std::vector<Vec2f> global_dir_random_samples(ray_count);
    generate(begin(global_dir_random_samples), end(global_dir_random_samples), gen);
    std::vector<Vec2f> local_dir_random_samples(ray_count);
    generate(begin(local_dir_random_samples), end(local_dir_random_samples), gen);

    BOOST_LOG_TRIVIAL(debug)
    << "PM: generate random samples: end";

    BOOST_LOG_TRIVIAL(debug)
    << "PM: raycast visibility for " << ray_count << " rays: start";
    // raycast visibility
    std::vector<HitInfo> hit_points = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, ray_count),
            std::vector<HitInfo> { },
            [&](tbb::blocked_range<size_t> r, std::vector<HitInfo> init) {
                for (size_t index = r.begin(); index < r.end(); ++index) {
                    Vec3d global_ray_dir = sample_sphere_uniform(global_dir_random_samples[index]);
                    Vec3d ray_origin = (vision_sphere_center - global_ray_dir * vision_sphere_raidus);
                    Vec3d local_dir = sample_power_cosine_hemisphere(local_dir_random_samples[index], 1.5);

                    Frame f;
                    f.set_from_z(global_ray_dir);
                    Vec3d final_ray_dir = (f.to_world(local_dir));

                    igl::Hit hitpoint;
                    // FIXME: This query will not compile for float ray origin and
                    // direction for some reason
                    auto hit = AABBTreeIndirect::intersect_ray_first_hit(triangles.vertices,
                            triangles.indices, raycasting_tree, ray_origin, final_ray_dir, hitpoint);

                    if (hit) {
                        auto face = triangles.indices[hitpoint.id];
                        auto edge1 = triangles.vertices[face[1]] - triangles.vertices[face[0]];
                        auto edge2 = triangles.vertices[face[2]] - triangles.vertices[face[0]];

                        Vec3d hit_pos = (triangles.vertices[face[0]] + edge1 * hitpoint.u + edge2 * hitpoint.v).cast<double>();
                        Vec3d surface_normal = edge1.cross(edge2).cast<double>();

                        init.push_back(HitInfo { hit_pos, surface_normal });
                    }
                }
                return init;
            },
            [](std::vector<HitInfo> left, std::vector<HitInfo> right) {
                left.insert(left.end(), right.begin(), right.end());
                return left;
            }
            );

    BOOST_LOG_TRIVIAL(debug)
    << "PM: raycast visibility for " << ray_count << " rays: end";

//TODO disable, only debug code
//#ifdef 0
    its_write_obj(triangles, "triangles.obj");

    Slic3r::CNumericLocalesSetter locales_setter;
    FILE *fp = boost::nowide::fopen("perimeter.obj", "w");
    if (fp == nullptr) {
        BOOST_LOG_TRIVIAL(error)
        << "Couldn't open " << "perimeter.obj" << " for writing";
    }

    for (size_t i = 0; i < hit_points.size(); ++i)
        fprintf(fp, "v %f %f %f \n", hit_points[i].m_position[0], hit_points[i].m_position[1],
                hit_points[i].m_position[2]);
    fclose(fp);
//#endif

    return hit_points;
}

} // namespace SeamPlacerImpl

void SeamPlacer::init(const Print &print) {
    using namespace SeamPlacerImpl;
    m_perimeter_points_trees_per_object.clear();
    m_perimeter_points_per_object.clear();

    for (const PrintObject *po : print.objects()) {
        BOOST_LOG_TRIVIAL(debug)
        << "PM: build AABB tree for raycasting: start";
        // Build AABB tree for raycasting
        auto obj_transform = po->trafo_centered();
        auto triangle_set = po->model_object()->raw_indexed_triangle_set();
        its_transform(triangle_set, obj_transform);

        auto raycasting_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(triangle_set.vertices,
                triangle_set.indices);

        BOOST_LOG_TRIVIAL(debug)
        << "PM: build AABB tree for raycasting: end";

        std::vector<HitInfo> hit_points = raycast_visibility(ray_count_per_object, raycasting_tree, triangle_set);
        HitInfoCoordinateFunctor hit_points_functor { &hit_points };
        KDTreeIndirect<3, coordf_t, HitInfoCoordinateFunctor> hit_points_tree { hit_points_functor, hit_points.size() };

        BOOST_LOG_TRIVIAL(debug)
        << "PM: gather and build KD tree with seam candidates: start";
        // gather seam candidates (perimeter points)
        m_perimeter_points_per_object[po] = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, po->layers().size()),
                std::vector<SeamCandidate> { }, [&](tbb::blocked_range<size_t> r, std::vector<SeamCandidate> init) {
                    for (size_t index = r.begin(); index < r.end(); ++index) {
                        const auto layer = po->layers()[index];
                        auto unscaled_z = layer->slice_z;
                        for (const ExPolygon &expoly : layer->lslices) {
                            // Contour - insert first point marked as Polygon
                            // start, then insert rest sequentially.
                            {
                                auto unscaled_p = unscale(expoly.contour[0]);
                                init.emplace_back(Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z }, true);
                            }
                            for (size_t index = 1; index < expoly.contour.size(); ++index) {
                                auto unscaled_p = unscale(expoly.contour[index]);
                                init.emplace_back(Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z });
                            }

                            for (const Polygon &hole : expoly.holes) {
                                // Perform the same for each hole
                                {
                                    auto unscaled_p = unscale(hole[0]);
                                    init.emplace_back(Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z }, true);
                                }
                                for (size_t index = 1; index < hole.size(); ++index) {
                                    auto unscaled_p = unscale(hole[index]);
                                    init.emplace_back(Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z });
                                }
                            }
                        }
                    }
                    return init;
                },
                [](std::vector<SeamCandidate> prev, std::vector<SeamCandidate> next) {
                    prev.insert(prev.end(), next.begin(), next.end());
                    return prev;
                });
        auto &perimeter_points = m_perimeter_points_per_object[po];
        // Build KD tree with seam candidates
        auto functor = KDTreeCoordinateFunctor { &perimeter_points };
        m_perimeter_points_trees_per_object.emplace(std::piecewise_construct, std::forward_as_tuple(po),
                std::forward_as_tuple(functor, m_perimeter_points_per_object[po].size()));
        SeamPlacer::PointTree &perimeter_points_tree = m_perimeter_points_trees_per_object.find(po)->second;

        BOOST_LOG_TRIVIAL(debug)
        << "PM: gather and build KD tree with seam candidates: end";
    }
}

void SeamPlacer::place_seam(const PrintObject *po, ExtrusionLoop &loop, coordf_t unscaled_z, const Point &last_pos,
        bool external_first,
        double nozzle_diameter, const EdgeGrid::Grid *lower_layer_edge_grid) {
    assert(m_perimeter_points_trees_per_object.find(po) != nullptr);
    assert(m_perimeter_points_per_object.find(po) != nullptr);
    const auto &perimeter_points_tree = m_perimeter_points_trees_per_object.find(po)->second;
    const auto &perimeter_points = m_perimeter_points_per_object.find(po)->second;

    Points loop_points { };
    loop.collect_points(loop_points);

    // vector of pairs: first-> index into perimeter points, second-> index into loop points
    auto closest_perimeter_point_indices = std::vector<std::pair<size_t, size_t>>(loop_points.size());
    for (size_t p_index = 0; p_index < loop_points.size(); ++p_index) {
        auto unscaled_p = unscale(loop_points[p_index]);
        closest_perimeter_point_indices.emplace_back(find_closest_point(perimeter_points_tree,
                Vec3d { unscaled_p.x(), unscaled_p.y(), unscaled_z }), p_index);
    }

    std::sort(closest_perimeter_point_indices.begin(), closest_perimeter_point_indices.end(),
            [&](const std::pair<size_t, size_t> left, const std::pair<size_t, size_t> right) {
                return perimeter_points[left.first].visibility < perimeter_points[right.first].visibility;
            });

    loop.split_at_vertex(loop_points[closest_perimeter_point_indices[0].second]);
//
//    Point seam = last_pos;
//    if (!loop.split_at_vertex(seam))
//        // The point is not in the original loop.
//        // Insert it.
//        loop.split_at(seam, true);
}

} // namespace Slic3r
