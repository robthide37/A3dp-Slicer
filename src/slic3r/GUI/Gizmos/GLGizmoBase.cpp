#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"

#include <GL/glew.h>

#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/Plater.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

// TODO: Display tooltips quicker on Linux

namespace Slic3r {
namespace GUI {

const float GLGizmoBase::Grabber::SizeFactor = 0.05f;
const float GLGizmoBase::Grabber::MinHalfSize = 1.5f;
const float GLGizmoBase::Grabber::DraggingScaleFactor = 1.25f;

#if ENABLE_RAYCAST_PICKING
PickingModel GLGizmoBase::Grabber::s_cube;
PickingModel GLGizmoBase::Grabber::s_cone;
#else
GLModel GLGizmoBase::Grabber::s_cube;
GLModel GLGizmoBase::Grabber::s_cone;
#endif // ENABLE_RAYCAST_PICKING

GLGizmoBase::Grabber::~Grabber()
{
#if ENABLE_RAYCAST_PICKING
    if (s_cube.model.is_initialized())
        s_cube.model.reset();

    if (s_cone.model.is_initialized())
        s_cone.model.reset();
#else
    if (s_cube.is_initialized())
        s_cube.reset();

    if (s_cone.is_initialized())
        s_cone.reset();
#endif // ENABLE_RAYCAST_PICKING
}

float GLGizmoBase::Grabber::get_half_size(float size) const
{
    return std::max(size * SizeFactor, MinHalfSize);
}

float GLGizmoBase::Grabber::get_dragging_half_size(float size) const
{
    return get_half_size(size) * DraggingScaleFactor;
}

#if ENABLE_RAYCAST_PICKING
void GLGizmoBase::Grabber::register_raycasters_for_picking(int id)
{
    picking_id = id;
    // registration will happen on next call to render()
}

void GLGizmoBase::Grabber::unregister_raycasters_for_picking()
{
    wxGetApp().plater()->canvas3D()->remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, picking_id);
    picking_id = -1;
    raycasters = { nullptr };
}
#endif // ENABLE_RAYCAST_PICKING

#if ENABLE_RAYCAST_PICKING
void GLGizmoBase::Grabber::render(float size, const ColorRGBA& render_color)
#else
void GLGizmoBase::Grabber::render(float size, const ColorRGBA& render_color, bool picking)
#endif // ENABLE_RAYCAST_PICKING
{
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_current_shader();
    if (shader == nullptr)
        return;
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_RAYCAST_PICKING
    if (!s_cube.model.is_initialized()) {
#else
    if (!s_cube.is_initialized()) {
#endif // ENABLE_RAYCAST_PICKING
        // This cannot be done in constructor, OpenGL is not yet
        // initialized at that point (on Linux at least).
        indexed_triangle_set its = its_make_cube(1.0, 1.0, 1.0);
        its_translate(its, -0.5f * Vec3f::Ones());
#if ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_RAYCAST_PICKING
        s_cube.model.init_from(its);
        s_cube.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(its)));
#else
        s_cube.init_from(its);
#endif // ENABLE_RAYCAST_PICKING
#else
        s_cube.init_from(its, BoundingBoxf3{ { -0.5, -0.5, -0.5 }, { 0.5, 0.5, 0.5 } });
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    }

#if ENABLE_RAYCAST_PICKING
    if (!s_cone.model.is_initialized()) {
        indexed_triangle_set its = its_make_cone(0.375, 1.5, double(PI) / 18.0);
        s_cone.model.init_from(its);
        s_cone.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(its)));
    }
#else
    if (!s_cone.is_initialized())
        s_cone.init_from(its_make_cone(0.375, 1.5, double(PI) / 18.0));
#endif // ENABLE_RAYCAST_PICKING

    const float half_size = dragging ? get_dragging_half_size(size) : get_half_size(size);

#if ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_RAYCAST_PICKING
    s_cube.model.set_color(render_color);
    s_cone.model.set_color(render_color);
#else
    s_cube.set_color(render_color);
    s_cone.set_color(render_color);
#endif // ENABLE_RAYCAST_PICKING

    const Camera& camera = wxGetApp().plater()->get_camera();
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#if ENABLE_RAYCAST_PICKING
    const Transform3d& view_matrix = camera.get_view_matrix();
    const Matrix3d view_matrix_no_offset = view_matrix.matrix().block(0, 0, 3, 3);
    std::vector<Transform3d> elements_matrices(GRABBER_ELEMENTS_MAX_COUNT, Transform3d::Identity());
    elements_matrices[0] = matrix * Geometry::assemble_transform(center, angles, 2.0 * half_size * Vec3d::Ones());
    Transform3d view_model_matrix = view_matrix * elements_matrices[0];
#else
    const Transform3d& view_matrix = camera.get_view_matrix();
    const Transform3d model_matrix = matrix * Geometry::assemble_transform(center, angles, 2.0 * half_size * Vec3d::Ones());
    const Transform3d view_model_matrix = view_matrix * model_matrix;
#endif // ENABLE_RAYCAST_PICKING

    shader->set_uniform("view_model_matrix", view_model_matrix);
#if ENABLE_RAYCAST_PICKING
    Matrix3d view_normal_matrix = view_matrix_no_offset * elements_matrices[0].matrix().block(0, 0, 3, 3).inverse().transpose();
#else
    const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
#endif // ENABLE_RAYCAST_PICKING
    shader->set_uniform("view_normal_matrix", view_normal_matrix);
#else
    s_cube.set_color(-1, render_color);
    s_cone.set_color(-1, render_color);
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(center.x(), center.y(), center.z()));
    glsafe(::glRotated(Geometry::rad2deg(angles.z()), 0.0, 0.0, 1.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.y()), 0.0, 1.0, 0.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.x()), 1.0, 0.0, 0.0));
    glsafe(::glScaled(2.0 * half_size, 2.0 * half_size, 2.0 * half_size));
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_RAYCAST_PICKING
    s_cube.model.render();
#else
    s_cube.render();
#endif // ENABLE_RAYCAST_PICKING

    auto render_extension = [&view_matrix, &view_matrix_no_offset, shader](const Transform3d& matrix) {
        const Transform3d view_model_matrix = view_matrix * matrix;
        shader->set_uniform("view_model_matrix", view_model_matrix);
        const Matrix3d view_normal_matrix = view_matrix_no_offset * matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
        shader->set_uniform("view_normal_matrix", view_normal_matrix);
        s_cone.model.render();
    };

#if ENABLE_LEGACY_OPENGL_REMOVAL
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosX)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[1] = elements_matrices[0] * Geometry::assemble_transform(Vec3d::UnitX(), Vec3d(0.0, 0.5 * double(PI), 0.0));
        render_extension(elements_matrices[1]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(Vec3d::UnitX(), Vec3d(0.0, 0.5 * double(PI), 0.0)));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegX)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[2] = elements_matrices[0] * Geometry::assemble_transform(-Vec3d::UnitX(), Vec3d(0.0, -0.5 * double(PI), 0.0));
        render_extension(elements_matrices[2]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(-Vec3d::UnitX(), Vec3d(0.0, -0.5 * double(PI), 0.0)));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosY)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[3] = elements_matrices[0] * Geometry::assemble_transform(Vec3d::UnitY(), Vec3d(-0.5 * double(PI), 0.0, 0.0));
        render_extension(elements_matrices[3]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(Vec3d::UnitY(), Vec3d(-0.5 * double(PI), 0.0, 0.0)));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegY)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[4] = elements_matrices[0] * Geometry::assemble_transform(-Vec3d::UnitY(), Vec3d(0.5 * double(PI), 0.0, 0.0));
        render_extension(elements_matrices[4]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(-Vec3d::UnitY(), Vec3d(0.5 * double(PI), 0.0, 0.0)));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosZ)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[5] = elements_matrices[0] * Geometry::assemble_transform(Vec3d::UnitZ());
        render_extension(elements_matrices[5]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(Vec3d::UnitZ()));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegZ)) != 0) {
#if ENABLE_RAYCAST_PICKING
        elements_matrices[6] = elements_matrices[0] * Geometry::assemble_transform(-Vec3d::UnitZ(), Vec3d(double(PI), 0.0, 0.0));
        render_extension(elements_matrices[6]);
#else
        shader->set_uniform("view_model_matrix", view_model_matrix * Geometry::assemble_transform(-Vec3d::UnitZ(), Vec3d(double(PI), 0.0, 0.0)));
        s_cone.render();
#endif // ENABLE_RAYCAST_PICKING
    }
#else
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosX)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(1.0, 0.0, 0.0));
        glsafe(::glRotated(0.5 * Geometry::rad2deg(double(PI)), 0.0, 1.0, 0.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegX)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(-1.0, 0.0, 0.0));
        glsafe(::glRotated(-0.5 * Geometry::rad2deg(double(PI)), 0.0, 1.0, 0.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosY)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(0.0, 1.0, 0.0));
        glsafe(::glRotated(-0.5 * Geometry::rad2deg(double(PI)), 1.0, 0.0, 0.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegY)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(0.0, -1.0, 0.0));
        glsafe(::glRotated(0.5 * Geometry::rad2deg(double(PI)), 1.0, 0.0, 0.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosZ)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(0.0, 0.0, 1.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
    if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegZ)) != 0) {
        glsafe(::glPushMatrix());
        glsafe(::glTranslated(0.0, 0.0, -1.0));
        glsafe(::glRotated(Geometry::rad2deg(double(PI)), 1.0, 0.0, 0.0));
        s_cone.render();
        glsafe(::glPopMatrix());
    }
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if !ENABLE_LEGACY_OPENGL_REMOVAL
    glsafe(::glPopMatrix());
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_RAYCAST_PICKING
    if (raycasters[0] == nullptr) {
        GLCanvas3D& canvas = *wxGetApp().plater()->canvas3D();
        raycasters[0] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cube.mesh_raycaster, elements_matrices[0]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosX)) != 0)
            raycasters[1] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[1]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegX)) != 0)
            raycasters[2] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[2]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosY)) != 0)
            raycasters[3] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[3]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegY)) != 0)
            raycasters[4] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[4]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::PosZ)) != 0)
            raycasters[5] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[5]);
        if ((int(extensions) & int(GLGizmoBase::EGrabberExtension::NegZ)) != 0)
            raycasters[6] = canvas.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, picking_id, *s_cone.mesh_raycaster, elements_matrices[6]);
    }
    else {
        for (size_t i = 0; i < GRABBER_ELEMENTS_MAX_COUNT; ++i) {
            if (raycasters[i] != nullptr)
                raycasters[i]->set_transform(elements_matrices[i]);
        }
    }
#endif // ENABLE_RAYCAST_PICKING
}

GLGizmoBase::GLGizmoBase(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : m_parent(parent)
    , m_icon_filename(icon_filename)
    , m_sprite_id(sprite_id)
    , m_imgui(wxGetApp().imgui())
{
}

void GLGizmoBase::set_hover_id(int id)
{
    // do not change hover id during dragging
    assert(!m_dragging);

    // allow empty grabbers when not using grabbers but use hover_id - flatten, rotate
    if (!m_grabbers.empty() && id >= (int) m_grabbers.size())
        return;
    
    m_hover_id = id;
    on_set_hover_id();    
}

bool GLGizmoBase::update_items_state()
{
    bool res = m_dirty;
    m_dirty  = false;
    return res;
}

#if ENABLE_RAYCAST_PICKING
void GLGizmoBase::register_grabbers_for_picking()
{
    for (size_t i = 0; i < m_grabbers.size(); ++i) {
        m_grabbers[i].register_raycasters_for_picking((m_group_id >= 0) ? m_group_id : i);
    }
}

void GLGizmoBase::unregister_grabbers_for_picking()
{
    for (size_t i = 0; i < m_grabbers.size(); ++i) {
        m_grabbers[i].unregister_raycasters_for_picking();
    }
}
#else
ColorRGBA GLGizmoBase::picking_color_component(unsigned int id) const
{
    id = BASE_ID - id;
    if (m_group_id > -1)
        id -= m_group_id;

    return picking_decode(id);
}
#endif // ENABLE_RAYCAST_PICKING

void GLGizmoBase::render_grabbers(const BoundingBoxf3& box) const
{
    render_grabbers((float)((box.size().x() + box.size().y() + box.size().z()) / 3.0));
}

void GLGizmoBase::render_grabbers(float size) const
{
    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;
    shader->start_using();
    shader->set_uniform("emission_factor", 0.1f);
    for (int i = 0; i < (int)m_grabbers.size(); ++i) {
        if (m_grabbers[i].enabled)
            m_grabbers[i].render(m_hover_id == i, size);
    }
    shader->stop_using();
}

#if !ENABLE_RAYCAST_PICKING
void GLGizmoBase::render_grabbers_for_picking(const BoundingBoxf3& box) const
{
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
        const float mean_size = float((box.size().x() + box.size().y() + box.size().z()) / 3.0);

        for (unsigned int i = 0; i < (unsigned int)m_grabbers.size(); ++i) {
            if (m_grabbers[i].enabled) {
                m_grabbers[i].color = picking_color_component(i);
                m_grabbers[i].render_for_picking(mean_size);
            }
        }
#if ENABLE_LEGACY_OPENGL_REMOVAL
        shader->stop_using();
    }
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}
#endif // !ENABLE_RAYCAST_PICKING

// help function to process grabbers
// call start_dragging, stop_dragging, on_dragging
bool GLGizmoBase::use_grabbers(const wxMouseEvent &mouse_event) {
    bool is_dragging_finished = false;
    if (mouse_event.Moving()) { 
        // it should not happen but for sure
        assert(!m_dragging);
        if (m_dragging) is_dragging_finished = true;
        else return false; 
    } 

    if (mouse_event.LeftDown()) {
        Selection &selection = m_parent.get_selection();        
        if (!selection.is_empty() && m_hover_id != -1 && 
            (m_grabbers.empty() || m_hover_id < static_cast<int>(m_grabbers.size()))) {
            selection.setup_cache();

            m_dragging = true;
            for (auto &grabber : m_grabbers) grabber.dragging = false;
            if (!m_grabbers.empty() && m_hover_id < int(m_grabbers.size()))
                m_grabbers[m_hover_id].dragging = true;            
            
            // prevent change of hover_id during dragging
            m_parent.set_mouse_as_dragging();
            on_start_dragging();

            // Let the plater know that the dragging started
            m_parent.post_event(SimpleEvent(EVT_GLCANVAS_MOUSE_DRAGGING_STARTED));
            m_parent.set_as_dirty();
            return true;
        }
    } else if (m_dragging) {
        // when mouse cursor leave window than finish actual dragging operation
        bool is_leaving = mouse_event.Leaving();
        if (mouse_event.Dragging()) {
            m_parent.set_mouse_as_dragging();
            Point      mouse_coord(mouse_event.GetX(), mouse_event.GetY());
            auto       ray = m_parent.mouse_ray(mouse_coord);
            UpdateData data(ray, mouse_coord);

            on_dragging(data);

            wxGetApp().obj_manipul()->set_dirty();
            m_parent.set_as_dirty();
            return true;
        }
        else if (mouse_event.LeftUp() || is_leaving || is_dragging_finished) {
#if ENABLE_WORLD_COORDINATE
            do_stop_dragging(is_leaving);
#else
            for (auto &grabber : m_grabbers) grabber.dragging = false;
            m_dragging = false;

            // NOTE: This should be part of GLCanvas3D
            // Reset hover_id when leave window
            if (is_leaving) m_parent.mouse_up_cleanup();

            on_stop_dragging();

            // There is prediction that after draggign, data are changed
            // Data are updated twice also by canvas3D::reload_scene.
            // Should be fixed.
            m_parent.get_gizmos_manager().update_data(); 

            wxGetApp().obj_manipul()->set_dirty();

            // Let the plater know that the dragging finished, so a delayed
            // refresh of the scene with the background processing data should
            // be performed.
            m_parent.post_event(SimpleEvent(EVT_GLCANVAS_MOUSE_DRAGGING_FINISHED));
            // updates camera target constraints
            m_parent.refresh_camera_scene_box();
#endif // ENABLE_WORLD_COORDINATE
            return true;
        }
    }
    return false;
}

#if ENABLE_WORLD_COORDINATE
void GLGizmoBase::do_stop_dragging(bool perform_mouse_cleanup)
{
    for (auto& grabber : m_grabbers) grabber.dragging = false;
    m_dragging = false;

    // NOTE: This should be part of GLCanvas3D
    // Reset hover_id when leave window
    if (perform_mouse_cleanup) m_parent.mouse_up_cleanup();

    on_stop_dragging();

    // There is prediction that after draggign, data are changed
    // Data are updated twice also by canvas3D::reload_scene.
    // Should be fixed.
    m_parent.get_gizmos_manager().update_data();

    wxGetApp().obj_manipul()->set_dirty();

    // Let the plater know that the dragging finished, so a delayed
    // refresh of the scene with the background processing data should
    // be performed.
    m_parent.post_event(SimpleEvent(EVT_GLCANVAS_MOUSE_DRAGGING_FINISHED));
    // updates camera target constraints
    m_parent.refresh_camera_scene_box();
}
#endif // ENABLE_WORLD_COORDINATE

std::string GLGizmoBase::format(float value, unsigned int decimals) const
{
    return Slic3r::string_printf("%.*f", decimals, value);
}

void GLGizmoBase::set_dirty() {
    m_dirty = true;
}

void GLGizmoBase::render_input_window(float x, float y, float bottom_limit)
{
    on_render_input_window(x, y, bottom_limit);
    if (m_first_input_window_render) {
        // for some reason, the imgui dialogs are not shown on screen in the 1st frame where they are rendered, but show up only with the 2nd rendered frame
        // so, we forces another frame rendering the first time the imgui window is shown
        m_parent.set_as_dirty();
        m_first_input_window_render = false;
    }
}



std::string GLGizmoBase::get_name(bool include_shortcut) const
{
    int key = get_shortcut_key();
    std::string out = on_get_name();
    if (include_shortcut && key >= WXK_CONTROL_A && key <= WXK_CONTROL_Z)
        out += std::string(" [") + char(int('A') + key - int(WXK_CONTROL_A)) + "]";
    return out;
}


} // namespace GUI
} // namespace Slic3r
