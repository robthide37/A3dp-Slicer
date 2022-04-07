#ifndef SRC_LIBSLIC3R_TRIANGLESELECTORWRAPPER_HPP_
#define SRC_LIBSLIC3R_TRIANGLESELECTORWRAPPER_HPP_

#include <memory>

#include "TriangleSelector.hpp"
#include "AABBTreeIndirect.hpp"

namespace Slic3r {

class TriangleSelectorWrapper {
public:
    const TriangleMesh &mesh;
    TriangleSelector selector;
    AABBTreeIndirect::Tree<3, float> triangles_tree;

    TriangleSelectorWrapper(const TriangleMesh &mesh) :
            mesh(mesh), selector(mesh), triangles_tree(AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(mesh.its.vertices, mesh.its.indices)) {

    }

    void enforce_spot(const Vec3f &point, float radius) {
        size_t hit_face_index;
        Vec3f hit_point;
        auto dist = AABBTreeIndirect::squared_distance_to_indexed_triangle_set(mesh.its.vertices, mesh.its.indices,
                triangles_tree,
                point, hit_face_index, hit_point);
        if (dist < 0 || dist > radius)
            return;

        std::unique_ptr<TriangleSelector::Cursor> cursor = std::make_unique<TriangleSelector::Sphere>(point, point,
                radius, Transform3d::Identity(), TriangleSelector::ClippingPlane { });

        selector.select_patch(hit_face_index, std::move(cursor), EnforcerBlockerType::ENFORCER, Transform3d::Identity(), true,
                0.0f);
    }

};

}

#endif /* SRC_LIBSLIC3R_TRIANGLESELECTORWRAPPER_HPP_ */
