#include "CameraUtils.hpp"
#include <igl/project.h> // projecting points

#include "slic3r/GUI/3DScene.hpp" // GLVolume
#include "libslic3r/Geometry/ConvexHull.hpp"

using namespace Slic3r;
using namespace GUI;

Points CameraUtils::project(const Camera &            camera,
                            const std::vector<Vec3d> &points)
{
    Vec4i viewport(camera.get_viewport().data());

    // Convert our std::vector to Eigen dynamic matrix.
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::DontAlign>
        pts(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i)
        pts.block<1, 3>(i, 0) = points[i];

    // Get the projections.
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::DontAlign> projections;
    igl::project(pts, camera.get_view_matrix().matrix(),
                    camera.get_projection_matrix().matrix(), viewport, projections);

    Points result;
    result.reserve(points.size());
    int window_height = viewport[3];

    // convert to points --> loss precision
    for (int i = 0; i < projections.rows(); ++i) {
        double x = projections(i, 0);
        double y = projections(i, 1);
        // opposit direction o Y
        result.emplace_back(x, window_height - y);
    }
    return result;
}

Point CameraUtils::project(const Camera &camera, const Vec3d &point)
{
    // IMPROVE: do it faster when you need it (inspire in project multi point)
    return project(camera, std::vector{point}).front();
}

Slic3r::Polygon CameraUtils::create_hull2d(const Camera &  camera,
                                   const GLVolume &volume)
{
    std::vector<Vec3d>  vertices;
    const TriangleMesh *hull = volume.convex_hull();
    if (hull != nullptr) {
        const indexed_triangle_set &its = hull->its;        
        vertices.reserve(its.vertices.size());
        // cast vector
        for (const Vec3f &vertex : its.vertices)
            vertices.emplace_back(vertex.cast<double>());
    } else {
        // Negative volume doesn't have convex hull so use bounding box
        auto bb = volume.bounding_box();
        Vec3d &min = bb.min;
        Vec3d &max = bb.max;
        vertices   = {min,
                    Vec3d(min.x(), min.y(), max.z()),
                    Vec3d(min.x(), max.y(), min.z()),
                    Vec3d(min.x(), max.y(), max.z()),
                    Vec3d(max.x(), min.y(), min.z()),
                    Vec3d(max.x(), min.y(), max.z()),
                    Vec3d(max.x(), max.y(), min.z()),
                    max};
    }

    const Transform3d &trafoMat =
        volume.get_instance_transformation().get_matrix() *
        volume.get_volume_transformation().get_matrix();
    for (Vec3d &vertex : vertices)
        vertex = trafoMat * vertex.cast<double>();

    Points vertices_2d = project(camera, vertices);
    return Geometry::convex_hull(vertices_2d);
}

#include <igl/unproject.h>
Vec3d CameraUtils::create_ray(const Camera &camera, const Vec2d &coor) {
    Matrix4d modelview  = camera.get_view_matrix().matrix();
    Matrix4d projection = camera.get_projection_matrix().matrix();
    Vec4i    viewport(camera.get_viewport().data());

    Vec3d scene_point(coor.x(), viewport[3] - coor.y(), 0.);
    Vec3d unprojected_point;
    igl::unproject(scene_point, modelview, projection, viewport, unprojected_point);

    Vec3d p0 = camera.get_position();
    Vec3d dir = unprojected_point - p0;
    dir.normalize();
    return dir;
}

Vec2d CameraUtils::get_z0_position(const Camera &camera, const Vec2d & coor)
{
    Vec3d dir = CameraUtils::create_ray(camera, coor);
    Vec3d p0  = camera.get_position();
    // find position of ray cross plane(z = 0)
    double t  = p0.z() / dir.z();
    Vec3d p = p0 - t * dir;
    return Vec2d(p.x(), p.y());
}
