// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Geometry/ConvexHull.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/SurfaceMesh.hpp"

#include <numeric>

#include <GL/glew.h>

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA DEFAULT_PLANE_COLOR       = { 0.9f, 0.9f, 0.9f, 0.5f };
static const Slic3r::ColorRGBA DEFAULT_HOVER_PLANE_COLOR = { 0.9f, 0.9f, 0.9f, 0.75f };

GLGizmoMeasure::GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{}



bool GLGizmoMeasure::on_mouse(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Moving()) {
        // only for sure 
        m_mouse_left_down = false;
        return false;
    }
    if (mouse_event.LeftDown()) {
        if (m_hover_id != -1) {
            m_mouse_left_down = true;
            Selection &selection = m_parent.get_selection();
            if (selection.is_single_full_instance()) {
                // Rotate the object so the normal points downward:
                selection.flattening_rotate(m_planes[m_hover_id].normal);
                m_parent.do_rotate(L("Gizmo-Place on Face"));
            }
            return true;
        }

        // fix: prevent restart gizmo when reselect object
        // take responsibility for left up
        if (m_parent.get_first_hover_volume_idx() >= 0) m_mouse_left_down = true;
        
    } else if (mouse_event.LeftUp()) {
        if (m_mouse_left_down) {
            // responsible for mouse left up after selecting plane
            m_mouse_left_down = false;
            return true;
        }
    } else if (mouse_event.Leaving()) {
        m_mouse_left_down = false;
    }
    return false;
}



void GLGizmoMeasure::data_changed()
{
    const Selection &  selection    = m_parent.get_selection();
    const ModelObject *model_object = nullptr;
    if (selection.is_single_full_instance() ||
        selection.is_from_single_object() ) {        
        model_object = selection.get_model()->objects[selection.get_object_idx()];
    }    
    set_flattening_data(model_object);
}



bool GLGizmoMeasure::on_init()
{
    // FIXME m_shortcut_key = WXK_CONTROL_F;
    return true;
}



void GLGizmoMeasure::on_set_state()
{
}



CommonGizmosDataID GLGizmoMeasure::on_get_requirements() const
{
    return CommonGizmosDataID::SelectionInfo;
}



std::string GLGizmoMeasure::on_get_name() const
{
    return _u8L("Measure");
}



bool GLGizmoMeasure::on_is_activable() const
{
    // This is assumed in GLCanvas3D::do_rotate, do not change this
    // without updating that function too.
    return m_parent.get_selection().is_single_full_instance();
}



void GLGizmoMeasure::on_render()
{
    const Selection& selection = m_parent.get_selection();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;
    
    shader->start_using();

    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glEnable(GL_BLEND));

    if (selection.is_single_full_instance()) {
        const Transform3d& m = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        if (this->is_plane_update_necessary())
            update_planes();
        for (int i = 0; i < (int)m_planes.size(); ++i) {
            m_planes[i].vbo.set_color(i == m_hover_id ? DEFAULT_HOVER_PLANE_COLOR : DEFAULT_PLANE_COLOR);
            m_planes[i].vbo.render();
        }
    }

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}





#if ! ENABLE_LEGACY_OPENGL_REMOVAL
    #error NOT IMPLEMENTED
#endif
#if ! ENABLE_GL_SHADERS_ATTRIBUTES
    #error NOT IMPLEMENTED
#endif





void GLGizmoMeasure::on_render_for_picking()
{
    const Selection& selection = m_parent.get_selection();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    shader->start_using();

    glsafe(::glDisable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_BLEND));

    if (selection.is_single_full_instance() && !wxGetKeyState(WXK_CONTROL)) {
        const Transform3d& m = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        if (this->is_plane_update_necessary())
            update_planes();
        for (int i = 0; i < (int)m_planes.size(); ++i) {
            m_planes[i].vbo.set_color(picking_color_component(i));
            m_planes[i].vbo.render();
        }
    }

    glsafe(::glEnable(GL_CULL_FACE));

    shader->stop_using();
}



void GLGizmoMeasure::set_flattening_data(const ModelObject* model_object)
{
    if (model_object != m_old_model_object) {
        m_planes.clear();
        m_planes_valid = false;
    }
}



void GLGizmoMeasure::update_planes()
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    TriangleMesh ch;
    for (const ModelVolume* vol : mo->volumes) {
        if (vol->type() != ModelVolumeType::MODEL_PART)
            continue;
        TriangleMesh vol_ch = vol->mesh();
        vol_ch.transform(vol->get_matrix());
        ch.merge(vol_ch);
    }
    m_planes.clear();
    const Transform3d& inst_matrix = mo->instances.front()->get_matrix();

    // Now we'll go through all the facets and append Points of facets sharing the same normal.
    // This part is still performed in mesh coordinate system.
    const int                num_of_facets  = ch.facets_count();
    std::vector<size_t>      face_to_plane(num_of_facets, 0);
    const std::vector<Vec3f> face_normals   = its_face_normals(ch.its);
    const std::vector<Vec3i> face_neighbors = its_face_neighbors(ch.its);
    std::vector<int>         facet_queue(num_of_facets, 0);
    std::vector<bool>        facet_visited(num_of_facets, false);
    int                      facet_queue_cnt = 0;
    const stl_normal*        normal_ptr      = nullptr;
    int facet_idx = 0;

    auto is_same_normal = [](const stl_normal& a, const stl_normal& b) -> bool {
        return (std::abs(a(0) - b(0)) < 0.001 && std::abs(a(1) - b(1)) < 0.001 && std::abs(a(2) - b(2)) < 0.001);
    };

    while (1) {
        // Find next unvisited triangle:
        for (; facet_idx < num_of_facets; ++ facet_idx)
            if (!facet_visited[facet_idx]) {
                facet_queue[facet_queue_cnt ++] = facet_idx;
                facet_visited[facet_idx] = true;
                normal_ptr = &face_normals[facet_idx];
                face_to_plane[facet_idx] = m_planes.size();
                m_planes.emplace_back();                
                break;
            }
        if (facet_idx == num_of_facets)
            break; // Everything was visited already

        while (facet_queue_cnt > 0) {
            int facet_idx = facet_queue[-- facet_queue_cnt];
            const stl_normal& this_normal = face_normals[facet_idx];
            if (is_same_normal(this_normal, *normal_ptr)) {
                const Vec3i& face = ch.its.indices[facet_idx];
                for (int j=0; j<3; ++j)
                    m_planes.back().vertices.emplace_back(ch.its.vertices[face[j]].cast<double>());

                facet_visited[facet_idx] = true;
                face_to_plane[facet_idx] = m_planes.size() - 1;
                for (int j = 0; j < 3; ++ j)
                    if (int neighbor_idx = face_neighbors[facet_idx][j]; neighbor_idx >= 0 && ! facet_visited[neighbor_idx])
                        facet_queue[facet_queue_cnt ++] = neighbor_idx;
            }
        }
        m_planes.back().normal = normal_ptr->cast<double>();

        Pointf3s& verts = m_planes.back().vertices;
        // Now we'll transform all the points into world coordinates, so that the areas, angles and distances
        // make real sense.
        verts = transform(verts, inst_matrix);
    }

    // Let's prepare transformation of the normal vector from mesh to instance coordinates.
    Geometry::Transformation t(inst_matrix);
    Vec3d scaling = t.get_scaling_factor();
    t.set_scaling_factor(Vec3d(1./scaling(0), 1./scaling(1), 1./scaling(2)));

    // Now we'll go through all the polygons, transform the points into xy plane to process them:
    for (unsigned int polygon_id=0; polygon_id < m_planes.size(); ++polygon_id) {
        Pointf3s& polygon = m_planes[polygon_id].vertices;
        const Vec3d& normal = m_planes[polygon_id].normal;

        // transform the normal according to the instance matrix:
        Vec3d normal_transformed = t.get_matrix() * normal;

        // We are going to rotate about z and y to flatten the plane
        Eigen::Quaterniond q;
        Transform3d m = Transform3d::Identity();
        m.matrix().block(0, 0, 3, 3) = q.setFromTwoVectors(normal_transformed, Vec3d::UnitZ()).toRotationMatrix();
        polygon = transform(polygon, m);

        // Now to remove the inner points. We'll misuse Geometry::convex_hull for that, but since
        // it works in fixed point representation, we will rescale the polygon to avoid overflows.
        // And yes, it is a nasty thing to do. Whoever has time is free to refactor.
        Vec3d bb_size = BoundingBoxf3(polygon).size();
        float sf = std::min(1./bb_size(0), 1./bb_size(1));
        Transform3d tr = Geometry::assemble_transform(Vec3d::Zero(), Vec3d::Zero(), Vec3d(sf, sf, 1.f));
        polygon = transform(polygon, tr);
        polygon = Slic3r::Geometry::convex_hull(polygon);
        polygon = transform(polygon, tr.inverse());

        // We will shrink the polygon a little bit so it does not touch the object edges:
        Vec3d centroid = std::accumulate(polygon.begin(), polygon.end(), Vec3d(0.0, 0.0, 0.0));
        centroid /= (double)polygon.size();
        for (auto& vertex : polygon)
            vertex = 0.95f*vertex + 0.05f*centroid;

        // Raise a bit above the object surface to avoid flickering:
        for (auto& b : polygon)
            b(2) += 0.1f;

        // Transform back to 3D (and also back to mesh coordinates)
        polygon = transform(polygon, inst_matrix.inverse() * m.inverse());
    }

    // We'll sort the planes by area and only keep the 254 largest ones (because of the picking pass limitations):
    std::sort(m_planes.rbegin(), m_planes.rend(), [](const PlaneData& a, const PlaneData& b) { return a.area < b.area; });
    m_planes.resize(std::min((int)m_planes.size(), 254));

    // Planes are finished - let's save what we calculated it from:
    m_volumes_matrices.clear();
    m_volumes_types.clear();
    for (const ModelVolume* vol : mo->volumes) {
        m_volumes_matrices.push_back(vol->get_matrix());
        m_volumes_types.push_back(vol->type());
    }
    m_first_instance_scale = mo->instances.front()->get_scaling_factor();
    m_first_instance_mirror = mo->instances.front()->get_mirror();
    m_old_model_object = mo;

    // And finally create respective VBOs. The polygon is convex with
    // the vertices in order, so triangulation is trivial.
    for (auto& plane : m_planes) {
        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::TriangleFan, GLModel::Geometry::EVertexLayout::P3N3 };
        init_data.reserve_vertices(plane.vertices.size());
        init_data.reserve_indices(plane.vertices.size());
        // vertices + indices
        for (size_t i = 0; i < plane.vertices.size(); ++i) {
            init_data.add_vertex((Vec3f)plane.vertices[i].cast<float>(), (Vec3f)plane.normal.cast<float>());
            init_data.add_index((unsigned int)i);
        }
        plane.vbo.init_from(std::move(init_data));

        // FIXME: vertices should really be local, they need not
        // persist now when we use VBOs
        plane.vertices.clear();
        plane.vertices.shrink_to_fit();
    }

    m_planes_valid = true;
}



bool GLGizmoMeasure::is_plane_update_necessary() const
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    if (m_state != On || ! mo || mo->instances.empty())
        return false;

    if (! m_planes_valid || mo != m_old_model_object
     || mo->volumes.size() != m_volumes_matrices.size())
        return true;

    // We want to recalculate when the scale changes - some planes could (dis)appear.
    if (! mo->instances.front()->get_scaling_factor().isApprox(m_first_instance_scale)
     || ! mo->instances.front()->get_mirror().isApprox(m_first_instance_mirror))
        return true;

    for (unsigned int i=0; i < mo->volumes.size(); ++i)
        if (! mo->volumes[i]->get_matrix().isApprox(m_volumes_matrices[i])
         || mo->volumes[i]->type() != m_volumes_types[i])
            return true;

    return false;
}

} // namespace GUI
} // namespace Slic3r
