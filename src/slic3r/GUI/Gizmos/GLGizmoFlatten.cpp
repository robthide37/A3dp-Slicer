// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoFlatten.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/GUI_App.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_GL_SHADERS_ATTRIBUTES
#include "slic3r/GUI/Plater.hpp"
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

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




// TESTING:
static Halfedge_index hi;
static bool hi_initialized = false;
static std::unique_ptr<SurfaceMesh> sm_ptr;
static Vertex_index src;
static Vertex_index tgt;






GLGizmoFlatten::GLGizmoFlatten(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{
    indexed_triangle_set a = its_make_cone(0.05, .2);
    its_rotate_x(a, M_PI);
    its_translate(a, stl_vertex(0., 0., .8));
    indexed_triangle_set b = its_make_cylinder(.02, 0.8);
    its_merge(a, b);
    arrow.init_from(a);
}

bool GLGizmoFlatten::on_mouse(const wxMouseEvent &mouse_event)
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
                hi = sm_ptr->halfedge(Face_index(m_hover_id));
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

void GLGizmoFlatten::data_changed()
{
    const Selection &  selection    = m_parent.get_selection();
    const ModelObject *model_object = nullptr;
    if (selection.is_single_full_instance() ||
        selection.is_from_single_object() ) {        
        model_object = selection.get_model()->objects[selection.get_object_idx()];
    }    
    set_flattening_data(model_object);
}

bool GLGizmoFlatten::on_init()
{
    m_shortcut_key = WXK_CONTROL_F;
    return true;
}

void GLGizmoFlatten::on_set_state()
{
}

CommonGizmosDataID GLGizmoFlatten::on_get_requirements() const
{
    return CommonGizmosDataID::SelectionInfo;
}

std::string GLGizmoFlatten::on_get_name() const
{
    return _u8L("Place on face");
}

bool GLGizmoFlatten::on_is_activable() const
{
    // This is assumed in GLCanvas3D::do_rotate, do not change this
    // without updating that function too.
    return m_parent.get_selection().is_single_full_instance();
}

void GLGizmoFlatten::on_render()
{
    const Selection& selection = m_parent.get_selection();

#if ENABLE_LEGACY_OPENGL_REMOVAL



    if (! hi_initialized) {
        const indexed_triangle_set& its = m_c->selection_info()->model_object()->volumes.front()->mesh().its;
        sm_ptr.reset(new SurfaceMesh(its));
        hi = sm_ptr->halfedge(Face_index(0));
        hi_initialized = true;
    }
    SurfaceMesh& sm = *sm_ptr;



    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;
    
    shader->start_using();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glEnable(GL_BLEND));

    if (selection.is_single_full_instance()) {
        const Transform3d& m = selection.get_first_volume()->get_instance_transformation().get_matrix();
#if ENABLE_GL_SHADERS_ATTRIBUTES
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_first_volume()->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#else
        glsafe(::glPushMatrix());
        glsafe(::glTranslatef(0.f, 0.f, selection.get_first_volume()->get_sla_shift_z()));
        glsafe(::glMultMatrixd(m.data()));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        if (this->is_plane_update_necessary())
            update_planes();
        for (int i = 0; i < (int)m_planes.size(); ++i) {
#if ENABLE_LEGACY_OPENGL_REMOVAL
        int cur_face = hi.is_invalid() ? 1000000 : sm.face(hi);
        for (int i=0; i < m_planes.size(); ++i) {
            m_planes[i].vbo.set_color(i == m_hover_id ? DEFAULT_HOVER_PLANE_COLOR : DEFAULT_PLANE_COLOR);
            if (i == cur_face)
                m_planes[i].vbo.set_color(i == m_hover_id ? ColorRGBA(.5f, 0.f, 0.f, 1.f) : ColorRGBA(1.f, 0.f, 0.f, 1.f));
             m_planes[i].vbo.render();
#else
            glsafe(::glColor4fv(i == m_hover_id ? DEFAULT_HOVER_PLANE_COLOR.data() : DEFAULT_PLANE_COLOR.data()));
            if (m_planes[i].vbo.has_VBOs())
                m_planes[i].vbo.render();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
        }
#if !ENABLE_GL_SHADERS_ATTRIBUTES
        glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
        }
    }




    /////////////////
    ////////////////
    //////////////////
    
    auto draw_arrow = [&](const Vec3d& from, const Vec3d& to) -> void {
        Vec3d desired_pos = from;
        Vec3d desired_dir = to - from;
        double desired_len = desired_dir.norm();
        desired_dir.normalize();

        Transform3d m = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix();
        m.translate(desired_pos);
        Eigen::Quaterniond q;
        Transform3d rot = Transform3d::Identity();
        rot.matrix().block(0, 0, 3, 3) = q.setFromTwoVectors(Vec3d::UnitZ(), desired_dir).toRotationMatrix();
        Transform3d sc = Transform3d::Identity();
        sc.scale(desired_len);
        m = m*sc*rot;

        const Camera& camera = wxGetApp().plater()->get_camera();
        Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z() * Vec3d::UnitZ()) * m;
        
        shader->set_uniform("view_model_matrix", view_model_matrix);
        arrow.render();
    };

    m_imgui->begin(std::string("DEBUG"));
    bool invalid = hi.is_invalid();
    if (invalid) {
        if (m_imgui->button(std::string("HALFEDGE INVALID (Click to reset)")))
            hi = sm.halfedge(Face_index(0));
    } else {
        m_imgui->text(sm.is_border(hi) ? "BORDER HALFEDGE !" : "Halfedge is not border");
        m_imgui->text((std::string("Face: ") + std::to_string(int(hi.face()))).c_str());
        m_imgui->text(std::string("Target degree:" + std::to_string(sm.degree(sm.target(hi)))));
        m_imgui->text(std::string("Face degree:" + std::to_string(sm.degree(sm.face(hi)))));
    }
    m_imgui->disabled_begin(invalid);
    if (m_imgui->button(std::string("next")))
        hi = sm.next(hi);
    if (m_imgui->button(std::string("prev")))
        hi = sm.prev(hi);
    if (m_imgui->button(std::string("opposite")))
        hi = sm.opposite(hi);
    if (m_imgui->button(std::string("next_around_target")))
        hi = sm.next_around_target(hi);
    if (m_imgui->button(std::string("prev_around_target")))
        hi = sm.prev_around_target(hi);
    if (m_imgui->button(std::string("next_around_source")))
        hi = sm.next_around_source(hi);
    if (m_imgui->button(std::string("prev_around_source")))
        hi = sm.prev_around_source(hi);
    if (m_imgui->button(std::string("remember one")))
        src = sm.target(hi);
    if (m_imgui->button(std::string("switch to halfedge"))) {
        tgt = sm.target(hi);
        hi = sm.halfedge(src, tgt);
    }
        
    if (invalid)
        m_imgui->disabled_end();
    m_imgui->end();

    if (! hi.is_invalid()) {
        Vec3d a = sm.point(sm.source(hi)).cast<double>();
        Vec3d b = sm.point(sm.target(hi)).cast<double>();
        draw_arrow(a, b);
    }
    

    /////////////////
    ////////////////
    //////////////////



    


    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

#if ENABLE_LEGACY_OPENGL_REMOVAL
    shader->stop_using();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

void GLGizmoFlatten::on_render_for_picking()
{
    const Selection& selection = m_parent.get_selection();

#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    shader->start_using();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    glsafe(::glDisable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_BLEND));

    if (selection.is_single_full_instance() && !wxGetKeyState(WXK_CONTROL)) {
        const Transform3d& m = selection.get_first_volume()->get_instance_transformation().get_matrix();
#if ENABLE_GL_SHADERS_ATTRIBUTES
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_first_volume()->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#else
        glsafe(::glPushMatrix());
        glsafe(::glTranslatef(0.f, 0.f, selection.get_first_volume()->get_sla_shift_z()));
        glsafe(::glMultMatrixd(m.data()));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        if (this->is_plane_update_necessary())
            update_planes();
        for (int i = 0; i < (int)m_planes.size(); ++i) {
#if ENABLE_LEGACY_OPENGL_REMOVAL
            m_planes[i].vbo.set_color(picking_color_component(i));
#else
            glsafe(::glColor4fv(picking_color_component(i).data()));
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
            m_planes[i].vbo.render();
        }
#if !ENABLE_GL_SHADERS_ATTRIBUTES
        glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
    }

    glsafe(::glEnable(GL_CULL_FACE));

#if ENABLE_LEGACY_OPENGL_REMOVAL
    shader->stop_using();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

void GLGizmoFlatten::set_flattening_data(const ModelObject* model_object)
{
    if (model_object != m_old_model_object) {
        m_planes.clear();
        m_planes_valid = false;
    }
}

void GLGizmoFlatten::update_planes()
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    TriangleMesh ch;
    for (const ModelVolume* vol : mo->volumes) {
        if (vol->type() != ModelVolumeType::MODEL_PART)
            continue;
        TriangleMesh vol_ch = vol->mesh(); //vol->get_convex_hull();
        vol_ch.transform(vol->get_matrix());
        ch.merge(vol_ch);
    }
    //ch = ch.convex_hull_3d();
    m_planes.clear();
#if ENABLE_WORLD_COORDINATE
    const Transform3d inst_matrix = mo->instances.front()->get_matrix_no_offset();
#else
    const Transform3d& inst_matrix = mo->instances.front()->get_matrix(true);
#endif // ENABLE_WORLD_COORDINATE

    // Following constants are used for discarding too small polygons.
    const float minimal_area = 5.f; // in square mm (world coordinates)
    const float minimal_side = 1.f; // mm

    // Now we'll go through all the facets and append Points of facets sharing the same normal.
    // This part is still performed in mesh coordinate system.
    const int                num_of_facets  = ch.facets_count();
    const std::vector<Vec3f> face_normals   = its_face_normals(ch.its);
    const std::vector<Vec3i> face_neighbors = its_face_neighbors(ch.its);
    std::vector<int>         facet_queue(num_of_facets, 0);
    std::vector<bool>        facet_visited(num_of_facets, false);
    int                      facet_queue_cnt = 0;
    const stl_normal*        normal_ptr      = nullptr;
    
    for (size_t i=0; i<ch.its.indices.size(); ++i) {
        const Vec3i& face = ch.its.indices[i];
        m_planes.emplace_back();
        for (int j=0; j<3; ++j)
            m_planes.back().vertices.emplace_back(ch.its.vertices[face[j]].cast<double>());
        m_planes.back().normal = face_normals[i].cast<double>();        

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
            vertex = 0.9f*vertex + 0.1f*centroid;

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
#if ENABLE_LEGACY_OPENGL_REMOVAL
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
#else
        plane.vbo.reserve(plane.vertices.size());
        for (const auto& vert : plane.vertices)
            plane.vbo.push_geometry(vert, plane.normal);
        for (size_t i=1; i<plane.vertices.size()-1; ++i)
            plane.vbo.push_triangle(0, i, i+1); // triangle fan
        plane.vbo.finalize_geometry(true);
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
        // FIXME: vertices should really be local, they need not
        // persist now when we use VBOs
        plane.vertices.clear();
        plane.vertices.shrink_to_fit();
    }

    m_planes_valid = true;
}


bool GLGizmoFlatten::is_plane_update_necessary() const
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
