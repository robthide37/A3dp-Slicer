// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp"
#if ENABLE_WORLD_COORDINATE
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#endif // ENABLE_WORLD_COORDINATE

#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/Jobs/RotoptimizeJob.hpp"

#include "libslic3r/PresetBundle.hpp"

#include <GL/glew.h>

namespace Slic3r {
namespace GUI {


const float GLGizmoRotate::Offset = 5.0f;
const unsigned int GLGizmoRotate::AngleResolution = 64;
const unsigned int GLGizmoRotate::ScaleStepsCount = 72;
const float GLGizmoRotate::ScaleStepRad = 2.0f * float(PI) / GLGizmoRotate::ScaleStepsCount;
const unsigned int GLGizmoRotate::ScaleLongEvery = 2;
const float GLGizmoRotate::ScaleLongTooth = 0.1f; // in percent of radius
const unsigned int GLGizmoRotate::SnapRegionsCount = 8;
const float GLGizmoRotate::GrabberOffset = 0.15f; // in percent of radius

GLGizmoRotate::GLGizmoRotate(GLCanvas3D& parent, GLGizmoRotate::Axis axis)
    : GLGizmoBase(parent, "", -1)
    , m_axis(axis)
    , m_angle(0.0)
    , m_center(0.0, 0.0, 0.0)
    , m_radius(0.0f)
    , m_snap_coarse_in_radius(0.0f)
    , m_snap_coarse_out_radius(0.0f)
    , m_snap_fine_in_radius(0.0f)
    , m_snap_fine_out_radius(0.0f)
    , m_drag_color(DEFAULT_DRAG_COLOR)
    , m_highlight_color(DEFAULT_HIGHLIGHT_COLOR)
{
    m_group_id = static_cast<int>(axis);
}

void GLGizmoRotate::set_highlight_color(const ColorRGBA &color)
{
    m_highlight_color = color;
}

void GLGizmoRotate::set_angle(double angle)
{
    if (std::abs(angle - 2.0 * double(PI)) < EPSILON)
        angle = 0.0;

    m_angle = angle;
}

std::string GLGizmoRotate::get_tooltip() const
{
    std::string axis;
    switch (m_axis)
    {
    case X: { axis = "X"; break; }
    case Y: { axis = "Y"; break; }
    case Z: { axis = "Z"; break; }
    }
    return (m_hover_id == 0 || m_grabbers.front().dragging) ? axis + ": " + format(float(Geometry::rad2deg(m_angle)), 4) : "";
}

bool GLGizmoRotate::on_mouse(const wxMouseEvent &mouse_event)
{
    return use_grabbers(mouse_event);
}

void GLGizmoRotate::dragging(const UpdateData &data) { on_dragging(data); }

void GLGizmoRotate::start_dragging()
{
    m_grabbers[0].dragging = true;
    on_start_dragging();
}

void GLGizmoRotate::stop_dragging()
{
    m_grabbers[0].dragging = false;
    on_stop_dragging();
}

void GLGizmoRotate::enable_grabber() { m_grabbers[0].enabled = true; }
void GLGizmoRotate::disable_grabber() { m_grabbers[0].enabled = false; }

bool GLGizmoRotate::on_init()
{
    m_grabbers.push_back(Grabber());
#if ENABLE_GIZMO_GRABBER_REFACTOR
    m_grabbers.back().extensions = (GLGizmoBase::EGrabberExtension)(int(GLGizmoBase::EGrabberExtension::PosY) | int(GLGizmoBase::EGrabberExtension::NegY));
#endif // ENABLE_GIZMO_GRABBER_REFACTOR
    return true;
}

void GLGizmoRotate::on_start_dragging()
{
#if ENABLE_WORLD_COORDINATE
    init_data_from_selection(m_parent.get_selection());
#else
    const BoundingBoxf3& box = m_parent.get_selection().get_bounding_box();
    m_center = box.center();
    m_radius = Offset + box.radius();
    m_snap_coarse_in_radius = m_radius / 3.0f;
    m_snap_coarse_out_radius = 2.0f * m_snap_coarse_in_radius;
    m_snap_fine_in_radius = m_radius;
    m_snap_fine_out_radius = m_snap_fine_in_radius + m_radius * ScaleLongTooth;
#endif // ENABLE_WORLD_COORDINATE
}

void GLGizmoRotate::on_dragging(const UpdateData &data)
{
    const Vec2d mouse_pos = to_2d(mouse_position_in_local_plane(data.mouse_ray, m_parent.get_selection()));

    const Vec2d orig_dir = Vec2d::UnitX();
    const Vec2d new_dir = mouse_pos.normalized();

    double theta = ::acos(std::clamp(new_dir.dot(orig_dir), -1.0, 1.0));
    if (cross2(orig_dir, new_dir) < 0.0)
        theta = 2.0 * (double)PI - theta;

    const double len = mouse_pos.norm();

    // snap to coarse snap region
    if (m_snap_coarse_in_radius <= len && len <= m_snap_coarse_out_radius) {
        const double step = 2.0 * double(PI) / double(SnapRegionsCount);
        theta = step * std::round(theta / step);
    }
    else {
        // snap to fine snap region (scale)
        if (m_snap_fine_in_radius <= len && len <= m_snap_fine_out_radius) {
            const double step = 2.0 * double(PI) / double(ScaleStepsCount);
            theta = step * std::round(theta / step);
        }
    }

    if (theta == 2.0 * double(PI))
        theta = 0.0;

    m_angle = theta;
}

void GLGizmoRotate::on_render()
{
    if (!m_grabbers.front().enabled)
        return;

#if !ENABLE_GIZMO_GRABBER_REFACTOR
    if (!m_cone.is_initialized())
        m_cone.init_from(its_make_cone(1.0, 1.0, double(PI) / 12.0));
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR

    const Selection& selection = m_parent.get_selection();
#if !ENABLE_WORLD_COORDINATE
    const BoundingBoxf3& box = selection.get_bounding_box();
#endif // !ENABLE_WORLD_COORDINATE

    if (m_hover_id != 0 && !m_grabbers.front().dragging) {
#if ENABLE_WORLD_COORDINATE
        init_data_from_selection(selection);
#else
        m_center = box.center();
        m_radius = Offset + box.radius();
        m_snap_coarse_in_radius = m_radius / 3.0f;
        m_snap_coarse_out_radius = 2.0f * m_snap_coarse_in_radius;
        m_snap_fine_in_radius = m_radius;
        m_snap_fine_out_radius = m_radius * (1.0f + ScaleLongTooth);
#endif // ENABLE_WORLD_COORDINATE
    }

    const double grabber_radius = (double)m_radius * (1.0 + (double)GrabberOffset);
    m_grabbers.front().center = Vec3d(::cos(m_angle) * grabber_radius, ::sin(m_angle) * grabber_radius, 0.0);
    m_grabbers.front().angles.z() = m_angle;

    glsafe(::glEnable(GL_DEPTH_TEST));

#if ENABLE_GL_SHADERS_ATTRIBUTES
    m_grabbers.front().matrix = local_transform(selection);
#else
    glsafe(::glPushMatrix());
    transform_to_local(selection);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    glsafe(::glLineWidth((m_hover_id != -1) ? 2.0f : 1.5f));
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();

#if ENABLE_GL_SHADERS_ATTRIBUTES
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() * m_grabbers.front().matrix;
        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

        const bool radius_changed = std::abs(m_old_radius - m_radius) > EPSILON;
        m_old_radius = m_radius;

        ColorRGBA color((m_hover_id != -1) ? m_drag_color : m_highlight_color);
        render_circle(color, radius_changed);
        if (m_hover_id != -1) {
            const bool hover_radius_changed = std::abs(m_old_hover_radius - m_radius) > EPSILON;
            m_old_hover_radius = m_radius;

            render_scale(color, hover_radius_changed);
            render_snap_radii(color, hover_radius_changed);
            render_reference_radius(color, hover_radius_changed);
            render_angle_arc(m_highlight_color, hover_radius_changed);
        }

        render_grabber_connection(color, radius_changed);
        shader->stop_using();
    }
#else
    glsafe(::glColor4fv((m_hover_id != -1) ? m_drag_color.data() : m_highlight_color.data()));

    render_circle();

    if (m_hover_id != -1) {
        render_scale();
        render_snap_radii();
        render_reference_radius();
    }

    glsafe(::glColor4fv(m_highlight_color.data()));

    if (m_hover_id != -1)
        render_angle();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_WORLD_COORDINATE
    render_grabber(m_bounding_box);
#if !ENABLE_GIZMO_GRABBER_REFACTOR
    render_grabber_extension(m_bounding_box, false);
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR
#else
    render_grabber(box);
#if !ENABLE_GIZMO_GRABBER_REFACTOR
    render_grabber_extension(box, false);
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR
#endif // ENABLE_WORLD_COORDINATE

#if !ENABLE_GL_SHADERS_ATTRIBUTES
    glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
}

void GLGizmoRotate::on_render_for_picking()
{
    const Selection& selection = m_parent.get_selection();

    glsafe(::glDisable(GL_DEPTH_TEST));

#if ENABLE_GL_SHADERS_ATTRIBUTES
    m_grabbers.front().matrix = local_transform(selection);
#else
    glsafe(::glPushMatrix());
    transform_to_local(selection);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

#if ENABLE_WORLD_COORDINATE
    render_grabbers_for_picking(m_bounding_box);
#if !ENABLE_GIZMO_GRABBER_REFACTOR
    render_grabber_extension(m_bounding_box, true);
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR
#else
    const BoundingBoxf3& box = selection.get_bounding_box();
    render_grabbers_for_picking(box);
#if !ENABLE_GIZMO_GRABBER_REFACTOR
    render_grabber_extension(box, true);
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR
#endif // ENABLE_WORLD_COORDINATE

#if !ENABLE_GL_SHADERS_ATTRIBUTES
    glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
}

#if ENABLE_WORLD_COORDINATE
void GLGizmoRotate::init_data_from_selection(const Selection& selection)
{
    ECoordinatesType coordinates_type;
    if (selection.is_wipe_tower())
        coordinates_type = ECoordinatesType::Local;
    else
        coordinates_type = wxGetApp().obj_manipul()->get_coordinates_type();
    if (coordinates_type == ECoordinatesType::World) {
        m_bounding_box = selection.get_bounding_box();
        m_center = m_bounding_box.center();
    }
    else if (coordinates_type == ECoordinatesType::Local && selection.is_single_volume_or_modifier()) {
        const GLVolume& v = *selection.get_first_volume();
        m_bounding_box = v.transformed_convex_hull_bounding_box(
            v.get_instance_transformation().get_scaling_factor_matrix() * v.get_volume_transformation().get_scaling_factor_matrix());
        m_center = v.world_matrix() * m_bounding_box.center();
    }
    else {
        m_bounding_box.reset();
        const Selection::IndicesList& ids = selection.get_volume_idxs();
        for (unsigned int id : ids) {
            const GLVolume& v = *selection.get_volume(id);
            m_bounding_box.merge(v.transformed_convex_hull_bounding_box(v.get_volume_transformation().get_matrix()));
        }
        const Geometry::Transformation inst_trafo = selection.get_first_volume()->get_instance_transformation();
        m_bounding_box = m_bounding_box.transformed(inst_trafo.get_scaling_factor_matrix());
        m_center = inst_trafo.get_matrix_no_scaling_factor() * m_bounding_box.center();
    }

    m_radius = Offset + m_bounding_box.radius();
    m_snap_coarse_in_radius = m_radius / 3.0f;
    m_snap_coarse_out_radius = 2.0f * m_snap_coarse_in_radius;
    m_snap_fine_in_radius = m_radius;
    m_snap_fine_out_radius = m_snap_fine_in_radius + m_radius * ScaleLongTooth;

    if (coordinates_type == ECoordinatesType::World)
        m_orient_matrix = Transform3d::Identity();
    else if (coordinates_type == ECoordinatesType::Local && (selection.is_wipe_tower() || selection.is_single_volume_or_modifier())) {
        const GLVolume& v = *selection.get_first_volume();
        m_orient_matrix = v.get_instance_transformation().get_rotation_matrix() * v.get_volume_transformation().get_rotation_matrix();
    }
    else {
        const GLVolume& v = *selection.get_first_volume();
        m_orient_matrix = v.get_instance_transformation().get_rotation_matrix();
    }
}
#endif // ENABLE_WORLD_COORDINATE

void GLGizmoRotate3D::on_render_input_window(float x, float y, float bottom_limit)
{
    if (wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() != ptSLA)
        return;

    RotoptimzeWindow popup{m_imgui, m_rotoptimizewin_state, {x, y, bottom_limit}};
}

void GLGizmoRotate3D::load_rotoptimize_state()
{
    std::string accuracy_str =
        wxGetApp().app_config->get("sla_auto_rotate", "accuracy");

    std::string method_str =
        wxGetApp().app_config->get("sla_auto_rotate", "method_id");

    if (!accuracy_str.empty()) {
        float accuracy = std::stof(accuracy_str);
        accuracy = std::max(0.f, std::min(accuracy, 1.f));

        m_rotoptimizewin_state.accuracy = accuracy;
    }

    if (!method_str.empty()) {
        int method_id = std::stoi(method_str);
        if (method_id < int(RotoptimizeJob::get_methods_count()))
            m_rotoptimizewin_state.method_id = method_id;
    }
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_circle(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_circle() const
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
#if ENABLE_LEGACY_OPENGL_REMOVAL
    if (!m_circle.is_initialized() || radius_changed) {
        m_circle.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::LineLoop, GLModel::Geometry::EVertexLayout::P3 };
        init_data.reserve_vertices(ScaleStepsCount);
        init_data.reserve_indices(ScaleStepsCount);

        // vertices + indices
        for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * ScaleStepRad);
            init_data.add_vertex(Vec3f(::cos(angle) * m_radius, ::sin(angle) * m_radius, 0.0f));
            init_data.add_index(i);
        }

        m_circle.init_from(std::move(init_data));
    }

    m_circle.set_color(color);
    m_circle.render();
#else
    ::glBegin(GL_LINE_LOOP);
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        const float angle = float(i) * ScaleStepRad;
        const float x = ::cos(angle) * m_radius;
        const float y = ::sin(angle) * m_radius;
        const float z = 0.0f;
        ::glVertex3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    }
    glsafe(::glEnd());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_scale(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_scale() const
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
    const float out_radius_long = m_snap_fine_out_radius;
    const float out_radius_short = m_radius * (1.0f + 0.5f * ScaleLongTooth);

#if ENABLE_LEGACY_OPENGL_REMOVAL
    if (!m_scale.is_initialized() || radius_changed) {
        m_scale.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.reserve_vertices(2 * ScaleStepsCount);
        init_data.reserve_indices(2 * ScaleStepsCount);

        // vertices + indices
        for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * ScaleStepRad);
            const float cosa = ::cos(angle);
            const float sina = ::sin(angle);
            const float in_x = cosa * m_radius;
            const float in_y = sina * m_radius;
            const float out_x = (i % ScaleLongEvery == 0) ? cosa * out_radius_long : cosa * out_radius_short;
            const float out_y = (i % ScaleLongEvery == 0) ? sina * out_radius_long : sina * out_radius_short;

            // vertices
            init_data.add_vertex(Vec3f(in_x, in_y, 0.0f));
            init_data.add_vertex(Vec3f(out_x, out_y, 0.0f));

            // indices
            init_data.add_line(i * 2, i * 2 + 1);
        }

        m_scale.init_from(std::move(init_data));
    }

    m_scale.set_color(color);
    m_scale.render();
#else
    ::glBegin(GL_LINES);
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        const float angle = (float)i * ScaleStepRad;
        const float cosa = ::cos(angle);
        const float sina = ::sin(angle);
        const float in_x = cosa * m_radius;
        const float in_y = sina * m_radius;
        const float in_z = 0.0f;
        const float out_x = (i % ScaleLongEvery == 0) ? cosa * out_radius_long : cosa * out_radius_short;
        const float out_y = (i % ScaleLongEvery == 0) ? sina * out_radius_long : sina * out_radius_short;
        const float out_z = 0.0f;
        ::glVertex3f((GLfloat)in_x, (GLfloat)in_y, (GLfloat)in_z);
        ::glVertex3f((GLfloat)out_x, (GLfloat)out_y, (GLfloat)out_z);
    }
    glsafe(::glEnd());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_snap_radii(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_snap_radii() const
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
    const float step = 2.0f * float(PI) / float(SnapRegionsCount);
    const float in_radius = m_radius / 3.0f;
    const float out_radius = 2.0f * in_radius;

#if ENABLE_LEGACY_OPENGL_REMOVAL
    if (!m_snap_radii.is_initialized() || radius_changed) {
        m_snap_radii.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.reserve_vertices(2 * ScaleStepsCount);
        init_data.reserve_indices(2 * ScaleStepsCount);

        // vertices + indices
        for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * step);
            const float cosa = ::cos(angle);
            const float sina = ::sin(angle);
            const float in_x = cosa * in_radius;
            const float in_y = sina * in_radius;
            const float out_x = cosa * out_radius;
            const float out_y = sina * out_radius;

            // vertices
            init_data.add_vertex(Vec3f(in_x, in_y, 0.0f));
            init_data.add_vertex(Vec3f(out_x, out_y, 0.0f));

            // indices
            init_data.add_line(i * 2, i * 2 + 1);
        }

        m_snap_radii.init_from(std::move(init_data));
    }

    m_snap_radii.set_color(color);
    m_snap_radii.render();
#else
    ::glBegin(GL_LINES);
    for (unsigned int i = 0; i < SnapRegionsCount; ++i) {
        const float angle = (float)i * step;
        const float cosa = ::cos(angle);
        const float sina = ::sin(angle);
        const float in_x = cosa * in_radius;
        const float in_y = sina * in_radius;
        const float in_z = 0.0f;
        const float out_x = cosa * out_radius;
        const float out_y = sina * out_radius;
        const float out_z = 0.0f;
        ::glVertex3f((GLfloat)in_x, (GLfloat)in_y, (GLfloat)in_z);
        ::glVertex3f((GLfloat)out_x, (GLfloat)out_y, (GLfloat)out_z);
    }
    glsafe(::glEnd());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_reference_radius(const ColorRGBA& color, bool radius_changed)
{
    if (!m_reference_radius.is_initialized() || radius_changed) {
        m_reference_radius.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex(Vec3f(m_radius * (1.0f + GrabberOffset), 0.0f, 0.0f));

        // indices
        init_data.add_line(0, 1);

        m_reference_radius.init_from(std::move(init_data));
    }

    m_reference_radius.set_color(color);
    m_reference_radius.render();
}
#else
void GLGizmoRotate::render_reference_radius() const
{
    ::glBegin(GL_LINES);
    ::glVertex3f(0.0f, 0.0f, 0.0f);
    ::glVertex3f((GLfloat)(m_radius * (1.0f + GrabberOffset)), 0.0f, 0.0f);
    glsafe(::glEnd());
}
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_angle_arc(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_angle() const
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
    const float step_angle = float(m_angle) / float(AngleResolution);
    const float ex_radius = m_radius * (1.0f + GrabberOffset);

#if ENABLE_LEGACY_OPENGL_REMOVAL
    const bool angle_changed = std::abs(m_old_angle - m_angle) > EPSILON;
    m_old_angle = m_angle;

    if (!m_angle_arc.is_initialized() || radius_changed || angle_changed) {
        m_angle_arc.reset();
        if (m_angle > 0.0f) {
            GLModel::Geometry init_data;
            init_data.format = { GLModel::Geometry::EPrimitiveType::LineStrip, GLModel::Geometry::EVertexLayout::P3 };
            init_data.reserve_vertices(1 + AngleResolution);
            init_data.reserve_indices(1 + AngleResolution);

            // vertices + indices
            for (unsigned int i = 0; i <= AngleResolution; ++i) {
                const float angle = float(i) * step_angle;
                init_data.add_vertex(Vec3f(::cos(angle) * ex_radius, ::sin(angle) * ex_radius, 0.0f));
                init_data.add_index(i);
            }

            m_angle_arc.init_from(std::move(init_data));
        }
    }

    m_angle_arc.set_color(color);
    m_angle_arc.render();
#else
    ::glBegin(GL_LINE_STRIP);
    for (unsigned int i = 0; i <= AngleResolution; ++i) {
        const float angle = float(i) * step_angle;
        const float x = ::cos(angle) * ex_radius;
        const float y = ::sin(angle) * ex_radius;
        const float z = 0.0f;
        ::glVertex3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    }
    glsafe(::glEnd());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void GLGizmoRotate::render_grabber_connection(const ColorRGBA& color, bool radius_changed)
{
    if (!m_grabber_connection.model.is_initialized() || radius_changed || !m_grabber_connection.old_center.isApprox(m_grabbers.front().center)) {
        m_grabber_connection.model.reset();
        m_grabber_connection.old_center = m_grabbers.front().center;

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex((Vec3f)m_grabbers.front().center.cast<float>());

        // indices
        init_data.add_line(0, 1);

        m_grabber_connection.model.init_from(std::move(init_data));
    }

    m_grabber_connection.model.set_color(color);
    m_grabber_connection.model.render();
}
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

void GLGizmoRotate::render_grabber(const BoundingBoxf3& box)
{
#if !ENABLE_LEGACY_OPENGL_REMOVAL
    const double grabber_radius = double(m_radius) * (1.0 + double(GrabberOffset));
    m_grabbers[0].center = Vec3d(::cos(m_angle) * grabber_radius, ::sin(m_angle) * grabber_radius, 0.0);
    m_grabbers[0].angles.z() = m_angle;

    glsafe(::glColor4fv((m_hover_id != -1) ? m_drag_color.data() : m_highlight_color.data()));

    ::glBegin(GL_LINES);
    ::glVertex3f(0.0f, 0.0f, 0.0f);
    ::glVertex3dv(m_grabbers[0].center.data());
    glsafe(::glEnd());
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

    m_grabbers.front().color = m_highlight_color;
    render_grabbers(box);
}

#if !ENABLE_GIZMO_GRABBER_REFACTOR
void GLGizmoRotate::render_grabber_extension(const BoundingBoxf3& box, bool picking)
{
    const float mean_size = float((box.size().x() + box.size().y() + box.size().z()) / 3.0);
    const double size = m_dragging ? double(m_grabbers.front().get_dragging_half_size(mean_size)) : double(m_grabbers.front().get_half_size(mean_size));

#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader(picking ? "flat" : "gouraud_light");
    if (shader == nullptr)
        return;

    m_cone.set_color((!picking && m_hover_id != -1) ? complementary(m_grabbers.front().color) : m_grabbers.front().color);

    shader->start_using();
    shader->set_uniform("emission_factor", 0.1f);
#else
    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    m_cone.set_color(-1, (!picking && m_hover_id != -1) ? complementary(m_grabbers.front().color) : m_grabbers.front().color);
    if (!picking) {
        shader->start_using();
        shader->set_uniform("emission_factor", 0.1f);
    }
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    const Vec3d& center = m_grabbers.front().center;

#if ENABLE_GL_SHADERS_ATTRIBUTES
    const Camera& camera = wxGetApp().plater()->get_camera();
    const Transform3d& view_matrix = camera.get_view_matrix();
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());

    Transform3d view_model_matrix = view_matrix * m_grabbers.front().matrix *
        Geometry::assemble_transform(center, Vec3d(0.5 * PI, 0.0, m_angle)) *
        Geometry::assemble_transform(2.0 * size * Vec3d::UnitZ(), Vec3d::Zero(), Vec3d(0.75 * size, 0.75 * size, 3.0 * size));

    shader->set_uniform("view_model_matrix", view_model_matrix);
    shader->set_uniform("normal_matrix", (Matrix3d)view_model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
#else
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(center.x(), center.y(), center.z()));
    glsafe(::glRotated(Geometry::rad2deg(m_angle), 0.0, 0.0, 1.0));
    glsafe(::glRotated(90.0, 1.0, 0.0, 0.0));
    glsafe(::glTranslated(0.0, 0.0, 2.0 * size));
    glsafe(::glScaled(0.75 * size, 0.75 * size, 3.0 * size));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
    m_cone.render();
#if ENABLE_GL_SHADERS_ATTRIBUTES
    view_model_matrix = view_matrix * m_grabbers.front().matrix *
        Geometry::assemble_transform(center, Vec3d(-0.5 * PI, 0.0, m_angle)) *
        Geometry::assemble_transform(2.0 * size * Vec3d::UnitZ(), Vec3d::Zero(), Vec3d(0.75 * size, 0.75 * size, 3.0 * size));

    shader->set_uniform("view_model_matrix", view_model_matrix);
    shader->set_uniform("normal_matrix", (Matrix3d)view_model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
#else
    glsafe(::glPopMatrix());
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(center.x(), center.y(), center.z()));
    glsafe(::glRotated(Geometry::rad2deg(m_angle), 0.0, 0.0, 1.0));
    glsafe(::glRotated(-90.0, 1.0, 0.0, 0.0));
    glsafe(::glTranslated(0.0, 0.0, 2.0 * size));
    glsafe(::glScaled(0.75 * size, 0.75 * size, 3.0 * size));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
    m_cone.render();
#if !ENABLE_GL_SHADERS_ATTRIBUTES
    glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES

#if !ENABLE_LEGACY_OPENGL_REMOVAL
    if (! picking)
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL
        shader->stop_using();
}
#endif // !ENABLE_GIZMO_GRABBER_REFACTOR

#if ENABLE_GL_SHADERS_ATTRIBUTES
Transform3d GLGizmoRotate::local_transform(const Selection& selection) const
{
    Transform3d ret;

    switch (m_axis)
    {
    case X:
    {
#if ENABLE_WORLD_COORDINATE
        ret = Geometry::rotation_transform(0.5 * PI * Vec3d::UnitY()) * Geometry::rotation_transform(-0.5 * PI * Vec3d::UnitZ());
#else
        ret = Geometry::assemble_transform(Vec3d::Zero(), 0.5 * PI * Vec3d::UnitY()) * Geometry::assemble_transform(Vec3d::Zero(), -0.5 * PI * Vec3d::UnitZ());
#endif // ENABLE_WORLD_COORDINATE
        break;
    }
    case Y:
    {
#if ENABLE_WORLD_COORDINATE
        ret = Geometry::rotation_transform(-0.5 * PI * Vec3d::UnitZ()) * Geometry::rotation_transform(-0.5 * PI * Vec3d::UnitY());
#else
        ret = Geometry::assemble_transform(Vec3d::Zero(), -0.5 * PI * Vec3d::UnitZ()) * Geometry::assemble_transform(Vec3d::Zero(), -0.5 * PI * Vec3d::UnitY());
#endif // ENABLE_WORLD_COORDINATE
        break;
    }
    default:
    case Z:
    {
        ret = Transform3d::Identity();
        break;
    }
    }

#if ENABLE_WORLD_COORDINATE
    return Geometry::translation_transform(m_center) * m_orient_matrix * ret;
#else
    if (selection.is_single_volume() || selection.is_single_modifier() || selection.requires_local_axes())
        ret = selection.get_first_volume()->get_instance_transformation().get_matrix(true, false, true, true) * ret;

    return Geometry::assemble_transform(m_center) * ret;
#endif // ENABLE_WORLD_COORDINATE
}
#else
void GLGizmoRotate::transform_to_local(const Selection& selection) const
{
    glsafe(::glTranslated(m_center.x(), m_center.y(), m_center.z()));

#if ENABLE_WORLD_COORDINATE
    glsafe(::glMultMatrixd(m_orient_matrix.data()));
#else
    if (selection.is_single_volume() || selection.is_single_modifier() || selection.requires_local_axes()) {
        const Transform3d orient_matrix = selection.get_first_volume()->get_instance_transformation().get_matrix(true, false, true, true);
        glsafe(::glMultMatrixd(orient_matrix.data()));
    }
#endif // ENABLE_WORLD_COORDINATE

    switch (m_axis)
    {
    case X:
    {
        glsafe(::glRotatef(90.0f, 0.0f, 1.0f, 0.0f));
        glsafe(::glRotatef(-90.0f, 0.0f, 0.0f, 1.0f));
        break;
    }
    case Y:
    {
        glsafe(::glRotatef(-90.0f, 0.0f, 0.0f, 1.0f));
        glsafe(::glRotatef(-90.0f, 0.0f, 1.0f, 0.0f));
        break;
    }
    default:
    case Z:
    {
        // no rotation
        break;
    }
    }
}
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

Vec3d GLGizmoRotate::mouse_position_in_local_plane(const Linef3& mouse_ray, const Selection& selection) const
{
    double half_pi = 0.5 * double(PI);

    Transform3d m = Transform3d::Identity();

    switch (m_axis)
    {
    case X:
    {
        m.rotate(Eigen::AngleAxisd(half_pi, Vec3d::UnitZ()));
        m.rotate(Eigen::AngleAxisd(-half_pi, Vec3d::UnitY()));
        break;
    }
    case Y:
    {
        m.rotate(Eigen::AngleAxisd(half_pi, Vec3d::UnitY()));
        m.rotate(Eigen::AngleAxisd(half_pi, Vec3d::UnitZ()));
        break;
    }
    default:
    case Z:
    {
        // no rotation applied
        break;
    }
    }

#if ENABLE_WORLD_COORDINATE
    m = m * m_orient_matrix.inverse();
#else
    if (selection.is_single_volume() || selection.is_single_modifier() || selection.requires_local_axes())
        m = m * selection.get_first_volume()->get_instance_transformation().get_matrix(true, false, true, true).inverse();
#endif // ENABLE_WORLD_COORDINATE

    m.translate(-m_center);

    return transform(mouse_ray, m).intersect_plane(0.0);
}

GLGizmoRotate3D::GLGizmoRotate3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_gizmos({ 
        GLGizmoRotate(parent, GLGizmoRotate::X), 
        GLGizmoRotate(parent, GLGizmoRotate::Y),
        GLGizmoRotate(parent, GLGizmoRotate::Z) })
{
    load_rotoptimize_state();
}

bool GLGizmoRotate3D::on_mouse(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Dragging() && m_dragging) {
        // Apply new temporary rotations
#if ENABLE_WORLD_COORDINATE
        TransformationType transformation_type;
        if (m_parent.get_selection().is_wipe_tower())
            transformation_type = TransformationType::Instance_Relative_Joint;
        else {
            switch (wxGetApp().obj_manipul()->get_coordinates_type())
            {
            default:
            case ECoordinatesType::World:    { transformation_type = TransformationType::World_Relative_Joint; break; }
            case ECoordinatesType::Instance: { transformation_type = TransformationType::Instance_Relative_Joint; break; }
            case ECoordinatesType::Local:    { transformation_type = TransformationType::Local_Relative_Joint; break; }
            }
        }
#else
        TransformationType transformation_type(TransformationType::World_Relative_Joint);
#endif // ENABLE_WORLD_COORDINATE
        if (mouse_event.AltDown())
            transformation_type.set_independent();
        m_parent.get_selection().rotate(get_rotation(), transformation_type);
    }
    return use_grabbers(mouse_event);
}

void GLGizmoRotate3D::data_changed() {
    if (m_parent.get_selection().is_wipe_tower()) {
#if !ENABLE_WORLD_COORDINATE
        const DynamicPrintConfig& config = wxGetApp().preset_bundle->prints.get_edited_preset().config;
        const float wipe_tower_rotation_angle =
            dynamic_cast<const ConfigOptionFloat*>(
                config.option("wipe_tower_rotation_angle"))->value;
        set_rotation(Vec3d(0., 0., (M_PI / 180.) * wipe_tower_rotation_angle));
#endif // !ENABLE_WORLD_COORDINATE
        m_gizmos[0].disable_grabber();
        m_gizmos[1].disable_grabber();
    }
    else {
#if !ENABLE_WORLD_COORDINATE
        set_rotation(Vec3d::Zero());
#endif // !ENABLE_WORLD_COORDINATE
        m_gizmos[0].enable_grabber();
        m_gizmos[1].enable_grabber();
    }
#if ENABLE_WORLD_COORDINATE
    set_rotation(Vec3d::Zero());
#endif // ENABLE_WORLD_COORDINATE
}

bool GLGizmoRotate3D::on_init()
{
    for (GLGizmoRotate& g : m_gizmos) 
        if (!g.init()) return false;

    for (unsigned int i = 0; i < 3; ++i)
        m_gizmos[i].set_highlight_color(AXES_COLOR[i]);

    m_shortcut_key = WXK_CONTROL_R;

    return true;
}

std::string GLGizmoRotate3D::on_get_name() const
{
    return _u8L("Rotate");
}

bool GLGizmoRotate3D::on_is_activable() const
{
    return !m_parent.get_selection().is_empty();
}

void GLGizmoRotate3D::on_start_dragging()
{
    assert(0 <= m_hover_id && m_hover_id < 3);
    m_gizmos[m_hover_id].start_dragging();
}

void GLGizmoRotate3D::on_stop_dragging()
{
    assert(0 <= m_hover_id && m_hover_id < 3);
    m_parent.do_rotate(L("Gizmo-Rotate"));
    m_gizmos[m_hover_id].stop_dragging();
}

void GLGizmoRotate3D::on_dragging(const UpdateData &data)
{
    assert(0 <= m_hover_id && m_hover_id < 3);
    m_gizmos[m_hover_id].dragging(data);
}

void GLGizmoRotate3D::on_render()
{
    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    if (m_hover_id == -1 || m_hover_id == 0)
        m_gizmos[X].render();

    if (m_hover_id == -1 || m_hover_id == 1)
        m_gizmos[Y].render();

    if (m_hover_id == -1 || m_hover_id == 2)
        m_gizmos[Z].render();
}

GLGizmoRotate3D::RotoptimzeWindow::RotoptimzeWindow(ImGuiWrapper *   imgui,
                                                    State &          state,
                                                    const Alignment &alignment)
    : m_imgui{imgui}
{
    imgui->begin(_L("Optimize orientation"), ImGuiWindowFlags_NoMove |
                                     ImGuiWindowFlags_AlwaysAutoResize |
                                     ImGuiWindowFlags_NoCollapse);

    // adjust window position to avoid overlap the view toolbar
    float win_h = ImGui::GetWindowHeight();
    float x = alignment.x, y = alignment.y;
    y = std::min(y, alignment.bottom_limit - win_h);
    ImGui::SetWindowPos(ImVec2(x, y), ImGuiCond_Always);

    float max_text_w = 0.;
    auto padding = ImGui::GetStyle().FramePadding;
    padding.x *= 2.f;
    padding.y *= 2.f;

    for (size_t i = 0; i < RotoptimizeJob::get_methods_count(); ++i) {
        float w =
            ImGui::CalcTextSize(RotoptimizeJob::get_method_name(i).c_str()).x +
            padding.x + ImGui::GetFrameHeight();
        max_text_w = std::max(w, max_text_w);
    }

    ImGui::PushItemWidth(max_text_w);

    if (ImGui::BeginCombo("", RotoptimizeJob::get_method_name(state.method_id).c_str())) {
        for (size_t i = 0; i < RotoptimizeJob::get_methods_count(); ++i) {
            if (ImGui::Selectable(RotoptimizeJob::get_method_name(i).c_str())) {
                state.method_id = i;
                wxGetApp().app_config->set("sla_auto_rotate",
                                           "method_id",
                                           std::to_string(state.method_id));
            }

            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", RotoptimizeJob::get_method_description(i).c_str());
        }

        ImGui::EndCombo();
    }

    ImVec2 sz = ImGui::GetItemRectSize();

    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", RotoptimizeJob::get_method_description(state.method_id).c_str());

    ImGui::Separator();

    auto btn_txt = _L("Apply");
    auto btn_txt_sz = ImGui::CalcTextSize(btn_txt.c_str());
    ImVec2 button_sz = {btn_txt_sz.x + padding.x, btn_txt_sz.y + padding.y};
    ImGui::SetCursorPosX(padding.x + sz.x - button_sz.x);

    if (!wxGetApp().plater()->get_ui_job_worker().is_idle())
        imgui->disabled_begin(true);

    if ( imgui->button(btn_txt) ) {
        replace_job(wxGetApp().plater()->get_ui_job_worker(),
                    std::make_unique<RotoptimizeJob>());
    }

    imgui->disabled_end();
}

GLGizmoRotate3D::RotoptimzeWindow::~RotoptimzeWindow()
{
    m_imgui->end();
}

} // namespace GUI
} // namespace Slic3r
