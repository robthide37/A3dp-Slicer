#include "MeshUtils.hpp"

#include "libslic3r/Tesselate.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Model.hpp"

#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/Camera.hpp"

#include <GL/glew.h>

#include <igl/unproject.h>

#include <cstdint>


namespace Slic3r {
namespace GUI {

void MeshClipper::set_behaviour(bool fill_cut, double contour_width)
{
    if (fill_cut != m_fill_cut || contour_width != m_contour_width)
        m_result.reset();
    m_fill_cut = fill_cut;
    m_contour_width = contour_width;
}



void MeshClipper::set_plane(const ClippingPlane& plane)
{
    if (m_plane != plane) {
        m_plane = plane;
        m_result.reset();
    }
}


void MeshClipper::set_limiting_plane(const ClippingPlane& plane)
{
    if (m_limiting_plane != plane) {
        m_limiting_plane = plane;
        m_result.reset();
    }
}



void MeshClipper::set_mesh(const TriangleMesh& mesh)
{
    if (m_mesh != &mesh) {
        m_mesh = &mesh;
        m_result.reset();
    }
}

void MeshClipper::set_negative_mesh(const TriangleMesh& mesh)
{
    if (m_negative_mesh != &mesh) {
        m_negative_mesh = &mesh;
        m_result.reset();
    }
}



void MeshClipper::set_transformation(const Geometry::Transformation& trafo)
{
    if (! m_trafo.get_matrix().isApprox(trafo.get_matrix())) {
        m_trafo = trafo;
        m_result.reset();
    }
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void MeshClipper::render_cut(const ColorRGBA& color)
#else
void MeshClipper::render_cut()
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
    if (! m_result)
        recalculate_triangles();
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* curr_shader = wxGetApp().get_current_shader();
    if (curr_shader != nullptr)
        curr_shader->stop_using();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();
        const Camera& camera = wxGetApp().plater()->get_camera();
        shader->set_uniform("view_model_matrix", camera.get_view_matrix());
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        for (CutIsland& isl : m_result->cut_islands) {
            isl.model.set_color(isl.disabled ? ColorRGBA(1.f, 0.f, 0.f, 1.f) : color);
            isl.model.render();
        }
        shader->stop_using();
    }

    if (curr_shader != nullptr)
        curr_shader->start_using();
#else
    if (m_vertex_array.has_VBOs())
        m_vertex_array.render();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}


#if ENABLE_LEGACY_OPENGL_REMOVAL
void MeshClipper::render_contour(const ColorRGBA& color)
#else
void MeshClipper::render_contour()
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
    if (! m_result)
        recalculate_triangles();
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* curr_shader = wxGetApp().get_current_shader();
    if (curr_shader != nullptr)
        curr_shader->stop_using();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();
        const Camera& camera = wxGetApp().plater()->get_camera();
        shader->set_uniform("view_model_matrix", camera.get_view_matrix());
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        for (CutIsland& isl : m_result->cut_islands) {
            isl.model_expanded.set_color(color);
            isl.model_expanded.render();
        }
        shader->stop_using();
    }

    if (curr_shader != nullptr)
        curr_shader->start_using();
#else
    if (m_vertex_array_expanded.has_VBOs())
        m_vertex_array_expanded.render();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

bool MeshClipper::contains(Vec3d point)
{
    if (!m_result)
        recalculate_triangles();

    for (CutIsland& isl : m_result->cut_islands) {
        BoundingBoxf3 bb = isl.model_expanded.get_bounding_box();

        // instead of using of standard bb.contains(point)
        // because of precision (Note, that model_expanded is pretranslate(0.003 * normal.normalized()))
        constexpr double pres = 0.01;
        bool ret = (point.x() > bb.min.x() || is_approx(point.x(), bb.min.x(), pres)) && (point.x() < bb.max.x() || is_approx(point.x(), bb.max.x(), pres))
                && (point.y() > bb.min.y() || is_approx(point.y(), bb.min.y(), pres)) && (point.y() < bb.max.y() || is_approx(point.y(), bb.max.y(), pres))
                && (point.z() > bb.min.z() || is_approx(point.z(), bb.min.z(), pres)) && (point.z() < bb.max.z() || is_approx(point.z(), bb.max.z(), pres));
        if (ret) {
            // when we detected, that model_expanded's bb contains a point, then check if its polygon contains this point 
            Vec3d point_inv = m_result->trafo.inverse() * point;
            Point pt = Point(scale_(point_inv.x()), scale_(point_inv.y()));
            if (isl.expoly.contains(pt))
                return true;
        }
    }
    return false;
}

bool MeshClipper::has_valid_contour()
{
    if (!m_result)
        recalculate_triangles();

    for (CutIsland& isl : m_result->cut_islands)
        if (isl.model_expanded.get_bounding_box().defined)
            return true;

    return false;
}


void MeshClipper::pass_mouse_click(const Vec3d& point_in)
{
    if (! m_result || m_result->cut_islands.empty())
        return;
    Vec3d point = m_result->trafo.inverse() * point_in;
    Point pt_2d = Point::new_scale(Vec2d(point.x(), point.y()));

    for (CutIsland& isl : m_result->cut_islands) {
        if (isl.expoly_bb.contains(pt_2d) && isl.expoly.contains(pt_2d))
            isl.disabled = ! isl.disabled;
    }
}

void MeshClipper::recalculate_triangles()
{
    m_result = ClipResult();

    auto plane_mesh = Eigen::Hyperplane<double, 3>(m_plane.get_normal(), -m_plane.distance(Vec3d::Zero())).transform(m_trafo.get_matrix().inverse());
    const Vec3d up = plane_mesh.normal();
    const float height_mesh = -plane_mesh.offset();

    // Now do the cutting
    MeshSlicingParams slicing_params;
    slicing_params.trafo.rotate(Eigen::Quaternion<double, Eigen::DontAlign>::FromTwoVectors(up, Vec3d::UnitZ()));

    ExPolygons expolys = union_ex(slice_mesh(m_mesh->its, height_mesh, slicing_params));

    if (m_negative_mesh && !m_negative_mesh->empty()) {
        const ExPolygons neg_expolys = union_ex(slice_mesh(m_negative_mesh->its, height_mesh, slicing_params));
        expolys = diff_ex(expolys, neg_expolys);
    }

    // Triangulate and rotate the cut into world coords:
    Eigen::Quaterniond q;
    q.setFromTwoVectors(Vec3d::UnitZ(), up);
    Transform3d tr = Transform3d::Identity();
    tr.rotate(q);
    tr = m_trafo.get_matrix() * tr;

    m_result->trafo = tr;

    if (m_limiting_plane != ClippingPlane::ClipsNothing())
    {
        // Now remove whatever ended up below the limiting plane (e.g. sinking objects).
        // First transform the limiting plane from world to mesh coords.
        // Note that inverse of tr transforms the plane from world to horizontal.
        const Vec3d normal_old = m_limiting_plane.get_normal().normalized();
        const Vec3d normal_new = (tr.matrix().block<3,3>(0,0).transpose() * normal_old).normalized();

        // normal_new should now be the plane normal in mesh coords. To find the offset,
        // transform a point and set offset so it belongs to the transformed plane.
        Vec3d pt = Vec3d::Zero();
        const double plane_offset = m_limiting_plane.get_data()[3];
        if (std::abs(normal_old.z()) > 0.5) // normal is normalized, at least one of the coords if larger than sqrt(3)/3 = 0.57
            pt.z() = - plane_offset / normal_old.z();
        else if (std::abs(normal_old.y()) > 0.5)
            pt.y() = - plane_offset / normal_old.y();
        else
            pt.x() = - plane_offset / normal_old.x();
        pt = tr.inverse() * pt;
        const double offset = -(normal_new.dot(pt));

        if (std::abs(normal_old.dot(m_plane.get_normal().normalized())) > 0.99) {
            // The cuts are parallel, show all or nothing.
            if (normal_old.dot(m_plane.get_normal().normalized()) < 0.0 && offset < height_mesh)
                expolys.clear();
        } else {
            // The cut is a horizontal plane defined by z=height_mesh.
            // ax+by+e=0 is the line of intersection with the limiting plane.
            // Normalized so a^2 + b^2 = 1.
            const double len = std::hypot(normal_new.x(), normal_new.y());
            if (len == 0.)
                return;
            const double a = normal_new.x() / len;
            const double b = normal_new.y() / len;
            const double e = (normal_new.z() * height_mesh + offset) / len;

            // We need a half-plane to limit the cut. Get angle of the intersecting line.
            double angle = (b != 0.0) ? std::atan(-a / b) : ((a < 0.0) ? -0.5 * M_PI : 0.5 * M_PI);
            if (b > 0) // select correct half-plane
                angle += M_PI;

            // We'll take a big rectangle above x-axis and rotate and translate
            // it so it lies on our line. This will be the figure to subtract
            // from the cut. The coordinates must not overflow after the transform,
            // make the rectangle a bit smaller.
            const coord_t size = (std::numeric_limits<coord_t>::max()/2 - scale_(std::max(std::abs(e*a), std::abs(e*b)))) / 4;
            Polygons ep {Polygon({Point(-size, 0), Point(size, 0), Point(size, 2*size), Point(-size, 2*size)})};
            ep.front().rotate(angle);
            ep.front().translate(scale_(-e * a), scale_(-e * b));
            expolys = diff_ex(expolys, ep);
        }
    }

    tr.pretranslate(0.001 * m_plane.get_normal().normalized()); // to avoid z-fighting
    Transform3d tr2 = tr;
    tr2.pretranslate(0.002 * m_plane.get_normal().normalized());


#if ENABLE_LEGACY_OPENGL_REMOVAL
    std::vector<Vec2f> triangles2d;

    for (const ExPolygon& exp : expolys) {
        triangles2d.clear();

        m_result->cut_islands.push_back(CutIsland());
        CutIsland& isl = m_result->cut_islands.back();

        if (m_fill_cut) {
            triangles2d = triangulate_expolygon_2f(exp, m_trafo.get_matrix().matrix().determinant() < 0.);
            GLModel::Geometry init_data;
            init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3 };
            init_data.reserve_vertices(triangles2d.size());
            init_data.reserve_indices(triangles2d.size());

            // vertices + indices
            for (auto it = triangles2d.cbegin(); it != triangles2d.cend(); it = it + 3) {
                init_data.add_vertex((Vec3f)(tr * Vec3d((*(it + 0)).x(), (*(it + 0)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                init_data.add_vertex((Vec3f)(tr * Vec3d((*(it + 1)).x(), (*(it + 1)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                init_data.add_vertex((Vec3f)(tr * Vec3d((*(it + 2)).x(), (*(it + 2)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                const size_t idx = it - triangles2d.cbegin();
                init_data.add_triangle((unsigned int)idx, (unsigned int)idx + 1, (unsigned int)idx + 2);
            }

            if (!init_data.is_empty())
                isl.model.init_from(std::move(init_data));
        }

        if (m_contour_width != 0. && ! exp.contour.empty()) {
            triangles2d.clear();

            // The contours must not scale with the object. Check the scale factor
            // in the respective directions, create a scaled copy of the ExPolygon
            // offset it and then unscale the result again.

            Transform3d t = tr;
            t.translation() = Vec3d::Zero();
            double scale_x = (t * Vec3d::UnitX()).norm();
            double scale_y = (t * Vec3d::UnitY()).norm();

            // To prevent overflow after scaling, downscale the input if needed:
            double extra_scale = 1.;
            int32_t limit = int32_t(std::min(std::numeric_limits<coord_t>::max() / (2. * scale_x), std::numeric_limits<coord_t>::max() / (2. * scale_y)));
            int32_t max_coord = 0;
            for (const Point& pt : exp.contour)
                max_coord = std::max(max_coord, std::max(std::abs(pt.x()), std::abs(pt.y())));
            if (max_coord + m_contour_width >= limit)
                extra_scale = 0.9 * double(limit) / max_coord;            

            ExPolygon exp_copy = exp;
            if (extra_scale != 1.)
                exp_copy.scale(extra_scale);
            exp_copy.scale(scale_x, scale_y);

            ExPolygons expolys_exp = offset_ex(exp_copy, scale_(m_contour_width));
            expolys_exp = diff_ex(expolys_exp, ExPolygons({exp_copy}));

            for (ExPolygon& e : expolys_exp) {
                e.scale(1./scale_x, 1./scale_y);
                if (extra_scale != 1.)
                    e.scale(1./extra_scale);
            }


            triangles2d = triangulate_expolygons_2f(expolys_exp, m_trafo.get_matrix().matrix().determinant() < 0.);
            GLModel::Geometry init_data = GLModel::Geometry();
            init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3 };
            init_data.reserve_vertices(triangles2d.size());
            init_data.reserve_indices(triangles2d.size());

            // vertices + indices
            for (auto it = triangles2d.cbegin(); it != triangles2d.cend(); it = it + 3) {
                init_data.add_vertex((Vec3f)(tr2 * Vec3d((*(it + 0)).x(), (*(it + 0)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                init_data.add_vertex((Vec3f)(tr2 * Vec3d((*(it + 1)).x(), (*(it + 1)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                init_data.add_vertex((Vec3f)(tr2 * Vec3d((*(it + 2)).x(), (*(it + 2)).y(), height_mesh)).cast<float>(), (Vec3f)up.cast<float>());
                const size_t idx = it - triangles2d.cbegin();
                init_data.add_triangle((unsigned short)idx, (unsigned short)idx + 1, (unsigned short)idx + 2);
            }

            if (!init_data.is_empty())
                isl.model_expanded.init_from(std::move(init_data));
        }

        isl.expoly = std::move(exp);
        isl.expoly_bb = get_extents(exp);
    }
#else
    #error NOT IMPLEMENTED
#endif // ENABLE_LEGACY_OPENGL_REMOVAL



#if ENABLE_LEGACY_OPENGL_REMOVAL
#else
    #error NOT IMPLEMENTED
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}


Vec3f MeshRaycaster::get_triangle_normal(size_t facet_idx) const
{
    return m_normals[facet_idx];
}

#if ENABLE_RAYCAST_PICKING
void MeshRaycaster::line_from_mouse_pos(const Vec2d& mouse_pos, const Transform3d& trafo, const Camera& camera,
                                        Vec3d& point, Vec3d& direction)
#else
void MeshRaycaster::line_from_mouse_pos(const Vec2d& mouse_pos, const Transform3d& trafo, const Camera& camera,
                                        Vec3d& point, Vec3d& direction)
#endif // ENABLE_RAYCAST_PICKING
{
    Matrix4d modelview = camera.get_view_matrix().matrix();
    Matrix4d projection= camera.get_projection_matrix().matrix();
    Vec4i viewport(camera.get_viewport().data());

    Vec3d pt1;
    Vec3d pt2;
    igl::unproject(Vec3d(mouse_pos.x(), viewport[3] - mouse_pos.y(), 0.),
                   modelview, projection, viewport, pt1);
    igl::unproject(Vec3d(mouse_pos.x(), viewport[3] - mouse_pos.y(), 1.),
                   modelview, projection, viewport, pt2);

    Transform3d inv = trafo.inverse();
    pt1 = inv * pt1;
    pt2 = inv * pt2;

    point = pt1;
    direction = pt2-pt1;
}


bool MeshRaycaster::unproject_on_mesh(const Vec2d& mouse_pos, const Transform3d& trafo, const Camera& camera,
                                      Vec3f& position, Vec3f& normal, const ClippingPlane* clipping_plane,
                                      size_t* facet_idx, bool* was_clipping_plane_hit) const
{
    if (was_clipping_plane_hit)
        *was_clipping_plane_hit = false;

    Vec3d point;
    Vec3d direction;
    line_from_mouse_pos(mouse_pos, trafo, camera, point, direction);

    std::vector<AABBMesh::hit_result> hits = m_emesh.query_ray_hits(point, direction);

    if (hits.empty())
        return false; // no intersection found

    unsigned i = 0;

    // Remove points that are obscured or cut by the clipping plane.
    // Also, remove anything below the bed (sinking objects).
    for (i=0; i<hits.size(); ++i) {
        Vec3d transformed_hit = trafo * hits[i].position();
        if (transformed_hit.z() >= SINKING_Z_THRESHOLD &&
            (! clipping_plane || ! clipping_plane->is_point_clipped(transformed_hit)))
            break;
    }

    if (i==hits.size()) {
        // All hits are clipped.
        return false;
    }
    if  ((hits.size()-i) % 2 != 0) {
        // There is an odd number of unclipped hits - meaning the nearest must be from inside the mesh.
        // In that case, calculate intersection with the clipping place.
        if (clipping_plane && was_clipping_plane_hit) {
            direction = direction + point;
            point = trafo * point; // transform to world coords
            direction = trafo * direction - point;

            Vec3d normal = -clipping_plane->get_normal().cast<double>();
            double den = normal.dot(direction);
            if (den != 0.) {
                double t = (-clipping_plane->get_offset() - normal.dot(point))/den;
                position = (point + t * direction).cast<float>();
                *was_clipping_plane_hit = true;
            }
        }
        return false;
    }

    // Now stuff the points in the provided vector and calculate normals if asked about them:
    position = hits[i].position().cast<float>();
    normal = hits[i].normal().cast<float>();

    if (facet_idx)
        *facet_idx = hits[i].face();

    return true;
}

bool MeshRaycaster::unproject_on_mesh(const Vec2d& mouse_pos, const Transform3d& trafo, const Camera& camera,
                                      Vec3d& position, Vec3d& normal) const
{
    Vec3d point;
    Vec3d direction;
    line_from_mouse_pos(mouse_pos, trafo, camera, point, direction);

    std::vector<AABBMesh::hit_result> hits = m_emesh.query_ray_hits(point, direction);

    if (hits.empty())
        return false; // no intersection found

    // Now stuff the points in the provided vector and calculate normals if asked about them:
    position = hits[0].position();
    normal = hits[0].normal();

    return true;
}

bool MeshRaycaster::is_valid_intersection(Vec3d point, Vec3d direction, const Transform3d& trafo) const 
{
    point = trafo.inverse() * point;

    std::vector<AABBMesh::hit_result> hits      = m_emesh.query_ray_hits(point, direction);
    std::vector<AABBMesh::hit_result> neg_hits  = m_emesh.query_ray_hits(point, -direction);

    return !hits.empty() && !neg_hits.empty();
}


std::vector<unsigned> MeshRaycaster::get_unobscured_idxs(const Geometry::Transformation& trafo, const Camera& camera, const std::vector<Vec3f>& points,
                                                       const ClippingPlane* clipping_plane) const
{
    std::vector<unsigned> out;

#if ENABLE_WORLD_COORDINATE
    const Transform3d instance_matrix_no_translation_no_scaling = trafo.get_rotation_matrix();
#else
    const Transform3d& instance_matrix_no_translation_no_scaling = trafo.get_matrix(true,false,true);
#endif // ENABLE_WORLD_COORDINATE
    Vec3d direction_to_camera = -camera.get_dir_forward();
    Vec3d direction_to_camera_mesh = (instance_matrix_no_translation_no_scaling.inverse() * direction_to_camera).normalized().eval();
    direction_to_camera_mesh = direction_to_camera_mesh.cwiseProduct(trafo.get_scaling_factor());
    const Transform3d inverse_trafo = trafo.get_matrix().inverse();

    for (size_t i=0; i<points.size(); ++i) {
        const Vec3f& pt = points[i];
        if (clipping_plane && clipping_plane->is_point_clipped(pt.cast<double>()))
            continue;

        bool is_obscured = false;
        // Cast a ray in the direction of the camera and look for intersection with the mesh:
        std::vector<AABBMesh::hit_result> hits;
        // Offset the start of the ray by EPSILON to account for numerical inaccuracies.
        hits = m_emesh.query_ray_hits((inverse_trafo * pt.cast<double>() + direction_to_camera_mesh * EPSILON),
                                      direction_to_camera_mesh);

        if (! hits.empty()) {
            // If the closest hit facet normal points in the same direction as the ray,
            // we are looking through the mesh and should therefore discard the point:
            if (hits.front().normal().dot(direction_to_camera_mesh.cast<double>()) > 0)
                is_obscured = true;

            // Eradicate all hits that the caller wants to ignore
            for (unsigned j=0; j<hits.size(); ++j) {
                if (clipping_plane && clipping_plane->is_point_clipped(trafo.get_matrix() * hits[j].position())) {
                    hits.erase(hits.begin()+j);
                    --j;
                }
            }

            // FIXME: the intersection could in theory be behind the camera, but as of now we only have camera direction.
            // Also, the threshold is in mesh coordinates, not in actual dimensions.
            if (! hits.empty())
                is_obscured = true;
        }
        if (! is_obscured)
            out.push_back(i);
    }
    return out;
}

#if ENABLE_RAYCAST_PICKING
bool MeshRaycaster::closest_hit(const Vec2d& mouse_pos, const Transform3d& trafo, const Camera& camera,
    Vec3f& position, Vec3f& normal, const ClippingPlane* clipping_plane, size_t* facet_idx) const
{
    Vec3d point;
    Vec3d direction;
    line_from_mouse_pos(mouse_pos, trafo, camera, point, direction);

    const std::vector<AABBMesh::hit_result> hits = m_emesh.query_ray_hits(point, direction.normalized());

    if (hits.empty())
        return false; // no intersection found

    size_t hit_id = 0;
    if (clipping_plane != nullptr) {
        while (hit_id < hits.size() && clipping_plane->is_point_clipped(trafo * hits[hit_id].position())) {
            ++hit_id;
        }
    }

    if (hit_id == hits.size())
        return false; // all points are obscured or cut by the clipping plane.

    const AABBMesh::hit_result& hit = hits[hit_id];

    position = hit.position().cast<float>();
    normal = hit.normal().cast<float>();

    if (facet_idx != nullptr)
        *facet_idx = hit.face();

    return true;
}
#endif // ENABLE_RAYCAST_PICKING

Vec3f MeshRaycaster::get_closest_point(const Vec3f& point, Vec3f* normal) const
{
    int idx = 0;
    Vec3d closest_point;
    Vec3d pointd = point.cast<double>();
    m_emesh.squared_distance(pointd, idx, closest_point);
    if (normal)
        // TODO: consider: get_normal(m_emesh, pointd).cast<float>();
        *normal = m_normals[idx];

    return closest_point.cast<float>();
}

int MeshRaycaster::get_closest_facet(const Vec3f &point) const
{
    int   facet_idx = 0;
    Vec3d closest_point;
    m_emesh.squared_distance(point.cast<double>(), facet_idx, closest_point);
    return facet_idx;
}

} // namespace GUI
} // namespace Slic3r
