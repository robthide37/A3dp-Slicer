#include "Model.hpp"
#include "TriangleSelectorWrapper.hpp"
#include <memory>

namespace Slic3r {

TriangleSelectorWrapper::TriangleSelectorWrapper(const TriangleMesh &mesh) :
        mesh(mesh), selector(mesh), triangles_tree(
                AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(mesh.its.vertices, mesh.its.indices)) {

}

void TriangleSelectorWrapper::enforce_spot(const Vec3f &point, const Vec3f &origin, float radius) {
    std::vector<igl::Hit> hits;
    Vec3f dir = (point - origin).normalized();
    if (AABBTreeIndirect::intersect_ray_all_hits(mesh.its.vertices, mesh.its.indices, triangles_tree,
            Vec3d(origin.cast<double>()),
            Vec3d(dir.cast<double>()),
            hits)) {
        for (int hit_idx = hits.size() - 1; hit_idx >= 0; --hit_idx) {
            const igl::Hit &hit = hits[hit_idx];
            Vec3f pos = origin + dir * hit.t;
            Vec3f face_normal = its_face_normal(mesh.its, hit.id);
            if (point.z() + radius > pos.z() && face_normal.dot(dir) < 0) {
                std::unique_ptr<TriangleSelector::Cursor> cursor = std::make_unique<TriangleSelector::Sphere>(
                        pos, origin, radius, Transform3d::Identity(), TriangleSelector::ClippingPlane { });
                selector.select_patch(hit.id, std::move(cursor), EnforcerBlockerType::ENFORCER, Transform3d::Identity(),
                        true, 0.0f);
                break;
            }
        }
    }
}

}
