// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp"

#include <GL/glew.h>

#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "libslic3r/PresetBundle.hpp"

#include "slic3r/GUI/Jobs/RotoptimizeJob.hpp"

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
{}

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

bool GLGizmoRotate::on_init()
{
    m_grabbers.push_back(Grabber());
    return true;
}

void GLGizmoRotate::on_start_dragging()
{
    const BoundingBoxf3& box = m_parent.get_selection().get_bounding_box();
    m_center = box.center();
    m_radius = Offset + box.radius();
    m_snap_coarse_in_radius = m_radius / 3.0f;
    m_snap_coarse_out_radius = 2.0f * m_snap_coarse_in_radius;
    m_snap_fine_in_radius = m_radius;
    m_snap_fine_out_radius = m_snap_fine_in_radius + m_radius * ScaleLongTooth;
}

void GLGizmoRotate::on_update(const UpdateData& data)
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

    if (!m_cone.is_initialized())
        m_cone.init_from(its_make_cone(1.0, 1.0, double(PI) / 12.0));

    const Selection& selection = m_parent.get_selection();
    const BoundingBoxf3& box = selection.get_bounding_box();

    if (m_hover_id != 0 && !m_grabbers.front().dragging) {
        m_center = box.center();
        m_radius = Offset + box.radius();
        m_snap_coarse_in_radius = m_radius / 3.0f;
        m_snap_coarse_out_radius = 2.0f * m_snap_coarse_in_radius;
        m_snap_fine_in_radius = m_radius;
        m_snap_fine_out_radius = m_radius * (1.0f + ScaleLongTooth);
    }

    const double grabber_radius = (double)m_radius * (1.0 + (double)GrabberOffset);
    m_grabbers.front().center = Vec3d(::cos(m_angle) * grabber_radius, ::sin(m_angle) * grabber_radius, 0.0);
    m_grabbers.front().angles.z() = m_angle;

    glsafe(::glEnable(GL_DEPTH_TEST));

    glsafe(::glPushMatrix());
    transform_to_local(selection);

    glsafe(::glLineWidth((m_hover_id != -1) ? 2.0f : 1.5f));
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();

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
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    render_grabber(box);
    render_grabber_extension(box, false);

    glsafe(::glPopMatrix());
}

void GLGizmoRotate::on_render_for_picking()
{
    const Selection& selection = m_parent.get_selection();

    glsafe(::glDisable(GL_DEPTH_TEST));

    glsafe(::glPushMatrix());

    transform_to_local(selection);

    const BoundingBoxf3& box = selection.get_bounding_box();
    render_grabbers_for_picking(box);
    render_grabber_extension(box, true);

    glsafe(::glPopMatrix());
}

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

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_circle(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_circle() const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (!m_circle.is_initialized() || radius_changed) {
        m_circle.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::LineLoop, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(ScaleStepsCount);
        init_data.reserve_indices(ScaleStepsCount);

        // vertices + indices
        for (unsigned short i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * ScaleStepRad);
            init_data.add_vertex(Vec3f(::cos(angle) * m_radius, ::sin(angle) * m_radius, 0.0f));
            init_data.add_ushort_index(i);
        }

        m_circle.init_from(std::move(init_data));
    }

    m_circle.set_color(color);
    m_circle.render();
#else
    ::glBegin(GL_LINE_LOOP);
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        float angle = (float)i * ScaleStepRad;
        float x = ::cos(angle) * m_radius;
        float y = ::sin(angle) * m_radius;
        float z = 0.0f;
        ::glVertex3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    }
    glsafe(::glEnd());
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_scale(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_scale() const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
    const float out_radius_long = m_snap_fine_out_radius;
    const float out_radius_short = m_radius * (1.0f + 0.5f * ScaleLongTooth);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (!m_scale.is_initialized() || radius_changed) {
        m_scale.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(2 * ScaleStepsCount);
        init_data.reserve_indices(2 * ScaleStepsCount);

        // vertices + indices
        for (unsigned short i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * ScaleStepRad);
            const float cosa = ::cos(angle);
            const float sina = ::sin(angle);
            const float in_x = cosa * m_radius;
            const float in_y = sina * m_radius;
            const float out_x = (i % ScaleLongEvery == 0) ? cosa * out_radius_long : cosa * out_radius_short;
            const float out_y = (i % ScaleLongEvery == 0) ? sina * out_radius_long : sina * out_radius_short;

            init_data.add_vertex(Vec3f(in_x, in_y, 0.0f));
            init_data.add_vertex(Vec3f(out_x, out_y, 0.0f));
            init_data.add_ushort_index(i * 2);
            init_data.add_ushort_index(i * 2 + 1);
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
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_snap_radii(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_snap_radii() const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
    const float step = 2.0f * float(PI) / float(SnapRegionsCount);
    const float in_radius = m_radius / 3.0f;
    const float out_radius = 2.0f * in_radius;

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (!m_snap_radii.is_initialized() || radius_changed) {
        m_snap_radii.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(2 * ScaleStepsCount);
        init_data.reserve_indices(2 * ScaleStepsCount);

        // vertices + indices
        for (unsigned short i = 0; i < ScaleStepsCount; ++i) {
            const float angle = float(i * step);
            const float cosa = ::cos(angle);
            const float sina = ::sin(angle);
            const float in_x = cosa * in_radius;
            const float in_y = sina * in_radius;
            const float out_x = cosa * out_radius;
            const float out_y = sina * out_radius;

            init_data.add_vertex(Vec3f(in_x, in_y, 0.0f));
            init_data.add_vertex(Vec3f(out_x, out_y, 0.0f));
            init_data.add_ushort_index(i * 2);
            init_data.add_ushort_index(i * 2 + 1);
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
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_reference_radius(const ColorRGBA& color, bool radius_changed)
{
    if (!m_reference_radius.is_initialized() || radius_changed) {
        m_reference_radius.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex(Vec3f(m_radius * (1.0f + GrabberOffset), 0.0f, 0.0f));

        // indices
        init_data.add_ushort_line(0, 1);

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
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_angle_arc(const ColorRGBA& color, bool radius_changed)
#else
void GLGizmoRotate::render_angle() const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
    const float step_angle = float(m_angle) / float(AngleResolution);
    const float ex_radius = m_radius * (1.0f + GrabberOffset);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (!m_angle_arc.is_initialized() || radius_changed) {
        m_angle_arc.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::LineStrip, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(1 + AngleResolution);
        init_data.reserve_indices(1 + AngleResolution);

        // vertices + indices
        for (unsigned short i = 0; i <= AngleResolution; ++i) {
            const float angle = float(i) * step_angle;
            init_data.add_vertex(Vec3f(::cos(angle) * ex_radius, ::sin(angle) * ex_radius, 0.0f));
            init_data.add_ushort_index(i);
        }

        m_angle_arc.init_from(std::move(init_data));
    }

    m_angle_arc.set_color(color);
    m_angle_arc.render();
#else
    ::glBegin(GL_LINE_STRIP);
    for (unsigned int i = 0; i <= AngleResolution; ++i) {
        float angle = (float)i * step_angle;
        float x = ::cos(angle) * ex_radius;
        float y = ::sin(angle) * ex_radius;
        float z = 0.0f;
        ::glVertex3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    }
    glsafe(::glEnd());
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLGizmoRotate::render_grabber_connection(const ColorRGBA& color, bool radius_changed)
{
    if (!m_grabber_connection.model.is_initialized() || radius_changed || !m_grabber_connection.old_center.isApprox(m_grabbers.front().center)) {
        m_grabber_connection.model.reset();
        m_grabber_connection.old_center = m_grabbers.front().center;

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex((Vec3f)m_grabbers.front().center.cast<float>());

        // indices
        init_data.add_ushort_line(0, 1);

        m_grabber_connection.model.init_from(std::move(init_data));
    }

    m_grabber_connection.model.set_color(color);
    m_grabber_connection.model.render();
}
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

void GLGizmoRotate::render_grabber(const BoundingBoxf3& box)
{
#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    const double grabber_radius = double(m_radius) * (1.0 + double(GrabberOffset));
    m_grabbers[0].center = Vec3d(::cos(m_angle) * grabber_radius, ::sin(m_angle) * grabber_radius, 0.0);
    m_grabbers[0].angles.z() = m_angle;

    glsafe(::glColor4fv((m_hover_id != -1) ? m_drag_color.data() : m_highlight_color.data()));

    ::glBegin(GL_LINES);
    ::glVertex3f(0.0f, 0.0f, 0.0f);
    ::glVertex3dv(m_grabbers[0].center.data());
    glsafe(::glEnd());
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

    m_grabbers.front().color = m_highlight_color;
    render_grabbers(box);
}

void GLGizmoRotate::render_grabber_extension(const BoundingBoxf3& box, bool picking)
{
    const float mean_size = float((box.size().x() + box.size().y() + box.size().z()) / 3.0);
    const double size = m_dragging ? double(m_grabbers.front().get_dragging_half_size(mean_size)) : double(m_grabbers.front().get_half_size(mean_size));

#if ENABLE_GLBEGIN_GLEND_REMOVAL
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
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    const Vec3d& center = m_grabbers.front().center;

    glsafe(::glPushMatrix());
    glsafe(::glTranslated(center.x(), center.y(), center.z()));
    glsafe(::glRotated(Geometry::rad2deg(m_angle), 0.0, 0.0, 1.0));
    glsafe(::glRotated(90.0, 1.0, 0.0, 0.0));
    glsafe(::glTranslated(0.0, 0.0, 2.0 * size));
    glsafe(::glScaled(0.75 * size, 0.75 * size, 3.0 * size));
    m_cone.render();
    glsafe(::glPopMatrix());
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(center.x(), center.y(), center.z()));
    glsafe(::glRotated(Geometry::rad2deg(m_angle), 0.0, 0.0, 1.0));
    glsafe(::glRotated(-90.0, 1.0, 0.0, 0.0));
    glsafe(::glTranslated(0.0, 0.0, 2.0 * size));
    glsafe(::glScaled(0.75 * size, 0.75 * size, 3.0 * size));
    m_cone.render();
    glsafe(::glPopMatrix());

#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    if (! picking)
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL
        shader->stop_using();
}

void GLGizmoRotate::transform_to_local(const Selection& selection) const
{
    glsafe(::glTranslated(m_center.x(), m_center.y(), m_center.z()));

    if (selection.is_single_volume() || selection.is_single_modifier() || selection.requires_local_axes()) {
        const Transform3d orient_matrix = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix(true, false, true, true);
        glsafe(::glMultMatrixd(orient_matrix.data()));
    }

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

    if (selection.is_single_volume() || selection.is_single_modifier() || selection.requires_local_axes())
        m = m * selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix(true, false, true, true).inverse();

    m.translate(-m_center);

    return transform(mouse_ray, m).intersect_plane(0.0);
}

GLGizmoRotate3D::GLGizmoRotate3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_gizmos({ GLGizmoRotate(parent, GLGizmoRotate::X), GLGizmoRotate(parent, GLGizmoRotate::Y), GLGizmoRotate(parent, GLGizmoRotate::Z) })
{
    for (unsigned int i = 0; i < 3; ++i) {
        m_gizmos[i].set_group_id(i);
    }

    load_rotoptimize_state();
}

bool GLGizmoRotate3D::on_init()
{
    for (GLGizmoRotate& g : m_gizmos) {
        if (!g.init())
            return false;
    }

    for (unsigned int i = 0; i < 3; ++i) {
        m_gizmos[i].set_highlight_color(AXES_COLOR[i]);
    }

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
    if (0 <= m_hover_id && m_hover_id < 3)
        m_gizmos[m_hover_id].start_dragging();
}

void GLGizmoRotate3D::on_stop_dragging()
{
    if (0 <= m_hover_id && m_hover_id < 3)
        m_gizmos[m_hover_id].stop_dragging();
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
