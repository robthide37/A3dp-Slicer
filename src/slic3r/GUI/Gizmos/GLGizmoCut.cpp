// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoCut.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"

#include <GL/glew.h>

#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/stattext.h>
#include <wx/sizer.h>

#include <algorithm>

#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "libslic3r/AppConfig.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"

namespace Slic3r {
namespace GUI {

static const double Margin = 20.0;
static const ColorRGBA GRABBER_COLOR = ColorRGBA::ORANGE();

#define use_grabber_extension 1

GLGizmoCut3D::GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_connector_type (CutConnectorType::Plug)
    , m_connector_style (size_t(CutConnectorStyle::Prizm))
    , m_connector_shape_id (size_t(CutConnectorShape::Hexagon))
    , m_rotation_gizmo(GLGizmoRotate3D(parent, "", -1))
    , m_rotation_matrix(Slic3r::Matrix3d::Identity())
{
    m_rotation_gizmo.use_only_grabbers();
    m_group_id = 3;
    m_connectors_group_id = 4;

    m_modes = { _u8L("Planar"),  _u8L("By Line"),_u8L("Grid")
//              , _u8L("Radial"), _u8L("Modular")
    };

    m_connector_modes = { _u8L("Auto"), _u8L("Manual") };
    m_connector_types = { _u8L("Plug"), _u8L("Dowel") };

    m_connector_styles = { _u8L("Prizm"), _u8L("Frustum")
//              , _u8L("Claw")
    };

    m_connector_shapes = { _u8L("Triangle"), _u8L("Square"), _u8L("Hexagon"), _u8L("Circle")
//              , _u8L("D-shape")
    };

    m_axis_names = { "X", "Y", "Z" };

    update_connector_shape();
}

std::string GLGizmoCut3D::get_tooltip() const
{
    std::string tooltip = m_rotation_gizmo.get_tooltip();
    if (tooltip.empty() &&
        (m_hover_id == m_group_id || m_grabbers[0].dragging)) {
        double koef = m_imperial_units ? ObjectManipulation::mm_to_in : 1.0;
        std::string unit_str = " " + (m_imperial_units ? _u8L("inch") : _u8L("mm"));
        const BoundingBoxf3 tbb = transformed_bounding_box();
        if (tbb.max.z() >= 0.0) {
            double top = (tbb.min.z() <= 0.0 ? tbb.max.z() : tbb.size().z()) * koef;
            tooltip += format(top, 2) + " " + unit_str + " (" + _u8L("Top part") + ")";
            if (tbb.min.z() <= 0.0)
                tooltip += "\n";
        }
        if (tbb.min.z() <= 0.0) {
            double bottom = (tbb.max.z() <= 0.0 ? tbb.size().z() : (tbb.min.z() * (-1))) * koef;
            tooltip += format(bottom, 2) + " " + unit_str + " (" + _u8L("Bottom part") + ")";
        }
        return tooltip;
    }

    return tooltip;
}

bool GLGizmoCut3D::on_mouse(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Moving()) 
        return false; 

    if (m_rotation_gizmo.on_mouse(mouse_event)) {
        update_clipper();
        return true;
    }

    if (use_grabbers(mouse_event)) 
        return true;

    Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
    Vec2d mouse_pos = mouse_coord.cast<double>();

    static bool pending_right_up = false;
    if (mouse_event.LeftDown()) {
        bool grabber_contains_mouse = (get_hover_id() != -1);
        bool control_down = mouse_event.CmdDown();
        if ((!control_down || grabber_contains_mouse) &&
            gizmo_event(SLAGizmoEventType::LeftDown, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), false))
            return true;
    }
    else if (mouse_event.Dragging()) {
        bool control_down = mouse_event.CmdDown();
        if (m_parent.get_move_volume_id() != -1) {
            // don't allow dragging objects with the Sla gizmo on
            return true;
        }
        else if (!control_down &&
            gizmo_event(SLAGizmoEventType::Dragging, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), false)) {
            // the gizmo got the event and took some action, no need to do
            // anything more here
            m_parent.set_as_dirty();
            return true;
        }
        else if (control_down && (mouse_event.LeftIsDown() || mouse_event.RightIsDown())) {
            // CTRL has been pressed while already dragging -> stop current action
            if (mouse_event.LeftIsDown())
                gizmo_event(SLAGizmoEventType::LeftUp, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), true);
            else if (mouse_event.RightIsDown())
                pending_right_up = false;
        }
    }
    else if (mouse_event.LeftUp() && !m_parent.is_mouse_dragging()) {
        // in case SLA/FDM gizmo is selected, we just pass the LeftUp event
        // and stop processing - neither object moving or selecting is
        // suppressed in that case
        gizmo_event(SLAGizmoEventType::LeftUp, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), mouse_event.CmdDown());
        return true;
    }
    else if (mouse_event.RightDown()) {
        if (m_parent.get_selection().get_object_idx() != -1 &&
            gizmo_event(SLAGizmoEventType::RightDown, mouse_pos, false, false, false)) {
            // we need to set the following right up as processed to avoid showing
            // the context menu if the user release the mouse over the object
            pending_right_up = true;
            // event was taken care of by the SlaSupports gizmo
            return true;
        }
    }
    else if (pending_right_up && mouse_event.RightUp()) {
        pending_right_up = false;
        return true;
    }
    return false;
}

void GLGizmoCut3D::shift_cut_z(double delta)
{
    Vec3d new_cut_center = m_plane_center;
    new_cut_center[Z] += delta;
    set_center_pos(new_cut_center);
}

void GLGizmoCut3D::rotate_vec3d_around_center(Vec3d& vec, const Vec3d& angles, const Vec3d& center)
{
    if (m_rotations != angles) {
        m_rotation_matrix =   Eigen::AngleAxisd(angles[Z], Vec3d::UnitZ())
                            * Eigen::AngleAxisd(angles[Y], Vec3d::UnitY())
                            * Eigen::AngleAxisd(angles[X], Vec3d::UnitX());
        m_rotations = angles;
    }

    vec -= center;
    vec = m_rotation_matrix * vec;
    vec += center;
}

void GLGizmoCut3D::put_connetors_on_cut_plane(const Vec3d& cp_normal, double cp_offset)
{
    ModelObject* mo = m_c->selection_info()->model_object();
    if (CutConnectors& connectors = mo->cut_connectors; !connectors.empty()) {
        const float sla_shift        = m_c->selection_info()->get_sla_shift();
        const Vec3d& instance_offset = mo->instances[m_c->selection_info()->get_active_instance()]->get_offset();

        for (auto& connector : connectors) {
            // convert connetor pos to the world coordinates
            Vec3d pos = connector.pos + instance_offset;
            pos[Z] += sla_shift;
            // scalar distance from point to plane along the normal
            double distance = -cp_normal.dot(pos) + cp_offset;
            // move connector
            connector.pos += distance * cp_normal;
        }
    }
}

void GLGizmoCut3D::update_clipper()
{
    const Vec3d& angles = m_rotation_gizmo.get_rotation();
    BoundingBoxf3 box = bounding_box();

    double radius = box.radius();

    Vec3d beg, end = beg = m_plane_center;
    beg[Z] = box.center().z() - radius;//box.min.z();
    end[Z] = box.center().z() + radius;//box.max.z();

    rotate_vec3d_around_center(beg, angles, m_plane_center);
    rotate_vec3d_around_center(end, angles, m_plane_center);

    double dist = (m_plane_center - beg).norm();

    // calculate normal and offset for clipping plane
    Vec3d normal = end - beg;
    dist = std::clamp(dist, 0.0001, normal.norm());
    normal.normalize();
    const double offset = normal.dot(beg) + dist;

    m_c->object_clipper()->set_range_and_pos(normal, offset, dist);

    put_connetors_on_cut_plane(normal, offset);
}

void GLGizmoCut3D::update_clipper_on_render()
{
    update_clipper();
    suppress_update_clipper_on_render = true;
}

void GLGizmoCut3D::set_center(const Vec3d& center)
{
    set_center_pos(center);
    m_rotation_gizmo.set_center(m_plane_center);
    update_clipper();
}

bool GLGizmoCut3D::render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(label);
    ImGui::SameLine(m_label_width);
    ImGui::PushItemWidth(m_control_width);
    
    size_t selection_out = selection_idx;
    // It is necessary to use BeginGroup(). Otherwise, when using SameLine() is called, then other items will be drawn inside the combobox.
    ImGui::BeginGroup();
    ImVec2 combo_pos = ImGui::GetCursorScreenPos();
    if (ImGui::BeginCombo(("##"+label).c_str(), "")) {
        for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx) {
            ImGui::PushID(int(line_idx));
            ImVec2 start_position = ImGui::GetCursorScreenPos();

            if (ImGui::Selectable("", line_idx == selection_idx))
                selection_out = line_idx;

            ImGui::SameLine();
            ImGui::Text("%s", lines[line_idx].c_str());
            ImGui::PopID();
        }

        ImGui::EndCombo();
    }

    ImVec2      backup_pos = ImGui::GetCursorScreenPos();
    ImGuiStyle& style = ImGui::GetStyle();

    ImGui::SetCursorScreenPos(ImVec2(combo_pos.x + style.FramePadding.x, combo_pos.y + style.FramePadding.y));
    ImGui::Text("%s", lines[selection_out].c_str());
    ImGui::SetCursorScreenPos(backup_pos);
    ImGui::EndGroup();

    bool is_changed = selection_idx != selection_out;
    selection_idx = selection_out;

    return is_changed;
}

bool GLGizmoCut3D::render_double_input(const std::string& label, double& value_in)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(label);
    ImGui::SameLine(m_label_width);
    ImGui::PushItemWidth(m_control_width);

    double value = value_in;
    if (m_imperial_units)
        value *= ObjectManipulation::mm_to_in;
    double old_val = value;
    ImGui::InputDouble(("##" + label).c_str(), &value, 0.0f, 0.0f, "%.2f", ImGuiInputTextFlags_CharsDecimal);

    ImGui::SameLine();
    m_imgui->text(m_imperial_units ? _L("in") : _L("mm"));

    value_in = value * (m_imperial_units ? ObjectManipulation::in_to_mm : 1.0);
    return old_val != value;
}

void GLGizmoCut3D::render_move_center_input(int axis)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(m_axis_names[axis]+":");
    ImGui::SameLine();
    ImGui::PushItemWidth(0.3*m_control_width);

    Vec3d move = m_plane_center;
    double in_val, value = in_val = move[axis];
    if (m_imperial_units)
        value *= ObjectManipulation::mm_to_in;
    ImGui::InputDouble(("##move_" + m_axis_names[axis]).c_str(), &value, 0.0f, 0.0f, "%.2f", ImGuiInputTextFlags_CharsDecimal);
    ImGui::SameLine();

    double val = value * (m_imperial_units ? ObjectManipulation::in_to_mm : 1.0);

    if (in_val != val) {
        move[axis] = val;
        set_center(move);
    }
}

void GLGizmoCut3D::render_rotation_input(int axis)
{
    m_imgui->text(m_axis_names[axis] + ":");
    ImGui::SameLine();

    Vec3d rotation = m_rotation_gizmo.get_rotation();
    double value = Geometry::rad2deg(rotation[axis]);
    if (value > 360)
        value -= 360;

    ImGui::PushItemWidth(0.3*m_control_width);
    ImGui::InputDouble(("##rotate_" + m_axis_names[axis]).c_str(), &value, 0.0f, 0.0f, "%.2f", ImGuiInputTextFlags_CharsDecimal);
    ImGui::SameLine();

    if (double val = Geometry::deg2rad(value); val != rotation[axis]) {
        rotation[axis] = val;
        m_rotation_gizmo.set_rotation(rotation);
        update_clipper();
    }
}

void GLGizmoCut3D::render_connect_type_radio_button(CutConnectorType type)
{
    ImGui::SameLine(type == CutConnectorType::Plug ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_types[size_t(type)], m_connector_type == type)) {
        m_connector_type = type;
        update_connector_shape();
    }
}

void GLGizmoCut3D::render_connect_mode_radio_button(CutConnectorMode mode)
{
    ImGui::SameLine(mode == CutConnectorMode::Auto ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_modes[int(mode)], m_connector_mode == mode))
        m_connector_mode = mode;
}

bool GLGizmoCut3D::render_revert_button(const std::string& label_id)
{
    const ImGuiStyle& style = ImGui::GetStyle();

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, { 1, style.ItemSpacing.y });
    ImGui::SameLine(m_label_width);

    ImGui::PushStyleColor(ImGuiCol_Button, { 0.25f, 0.25f, 0.25f, 0.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, { 0.4f, 0.4f, 0.4f, 1.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, { 0.4f, 0.4f, 0.4f, 1.0f });

    std::string label;
    label += ImGui::RevertButton;
    bool revert = ImGui::Button((label + "##" + label_id).c_str());

    ImGui::PopStyleColor(3);

    if (ImGui::IsItemHovered())
        m_imgui->tooltip(into_u8(_L("Revert")).c_str(), ImGui::GetFontSize() * 20.0f);

    ImGui::PopStyleVar();

    ImGui::SameLine();
    return revert;
}

void GLGizmoCut3D::render_cut_plane()
{
#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_CULL_FACE));
    glsafe(::glEnable(GL_BLEND));
    glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

    shader->start_using();
#if ENABLE_GL_SHADERS_ATTRIBUTES
    const Camera& camera = wxGetApp().plater()->get_camera();
    const Transform3d view_model_matrix = camera.get_view_matrix() * Geometry::assemble_transform(
        m_plane_center,
        m_rotation_gizmo.get_rotation(),
        Vec3d::Ones(),
        Vec3d::Ones()
    );
    shader->set_uniform("view_model_matrix", view_model_matrix);
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#else
    const Vec3d& angles = m_rotation_gizmo.get_rotation();
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(m_plane_center.x(), m_plane_center.y(), m_plane_center.z()));
    glsafe(::glRotated(Geometry::rad2deg(angles.z()), 0.0, 0.0, 1.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.y()), 0.0, 1.0, 0.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.x()), 1.0, 0.0, 0.0));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    if (!m_plane.is_initialized()) {
        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3 };
        init_data.color = { 0.8f, 0.8f, 0.8f, 0.5f };
        init_data.reserve_vertices(4);
        init_data.reserve_indices(6);

        const BoundingBoxf3 bb = bounding_box();
        const float min_x = bb.min.x() - Margin - m_plane_center.x();
        const float max_x = bb.max.x() + Margin - m_plane_center.x();
        const float min_y = bb.min.y() - Margin - m_plane_center.y();
        const float max_y = bb.max.y() + Margin - m_plane_center.y();

        // vertices
        init_data.add_vertex(Vec3f(min_x, min_y, 0.0));
        init_data.add_vertex(Vec3f(max_x, min_y, 0.0));
        init_data.add_vertex(Vec3f(max_x, max_y, 0.0));
        init_data.add_vertex(Vec3f(min_x, max_y, 0.0));

        // indices
        init_data.add_triangle(0, 1, 2);
        init_data.add_triangle(2, 3, 0);

        m_plane.init_from(std::move(init_data));
    }

    m_plane.render();
#if !ENABLE_GL_SHADERS_ATTRIBUTES
    glsafe(::glPopMatrix());
#endif //!ENABLE_GL_SHADERS_ATTRIBUTES
#else
    // Draw the cutting plane
    ::glBegin(GL_QUADS);
    ::glColor4fv(PLANE_COLOR.data());
    ::glVertex3f(min_x, min_y, plane_center.z());
    ::glVertex3f(max_x, min_y, plane_center.z());
    ::glVertex3f(max_x, max_y, plane_center.z());
    ::glVertex3f(min_x, max_y, plane_center.z());
    glsafe(::glEnd());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}

void GLGizmoCut3D::render_cut_center_graber()
{
    const Vec3d& angles = m_rotation_gizmo.get_rotation();
    m_grabbers[0].angles = angles;
    m_grabbers[0].color  = GRABBER_COLOR;

#if use_grabber_extension
    // UI experiments with grabber

    m_grabbers[0].center = m_plane_center;

    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (!shader)
        return;

    auto color = m_hover_id == m_group_id ? complementary(GRABBER_COLOR) : GRABBER_COLOR;
    m_sphere.set_color(color);
    m_cone.set_color(color);

    const Camera& camera = wxGetApp().plater()->get_camera();
    const Grabber& graber = m_grabbers.front();
    const Vec3d& center = graber.center;

    const BoundingBoxf3 box = bounding_box();

    const double mean_size = float((box.size().x() + box.size().y() + box.size().z()) / 6.0);
    const double size = m_dragging ? double(graber.get_dragging_half_size(mean_size)) : double(graber.get_half_size(mean_size));

    const Vec3d scale = Vec3d(0.75 * size, 0.75 * size, 1.8 * size);
    const Vec3d offset = 1.25 * size * Vec3d::UnitZ();

    shader->start_using();
    shader->set_uniform("emission_factor", 0.1f);
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());

    const Transform3d view_matrix = camera.get_view_matrix() * graber.matrix *
        Geometry::assemble_transform(center, angles);

    Transform3d view_model_matrix = view_matrix * Geometry::assemble_transform(Vec3d::Zero(), Vec3d::Zero(), size * Vec3d::Ones());

    shader->set_uniform("view_model_matrix", view_model_matrix);
    shader->set_uniform("normal_matrix", (Matrix3d)view_model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
    m_sphere.render();

    const BoundingBoxf3 tbb = transformed_bounding_box();
    if (tbb.max.z() >= 0.0) {
        view_model_matrix = view_matrix * Geometry::assemble_transform(offset, Vec3d::Zero(), scale);

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("normal_matrix", (Matrix3d)view_model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
        m_cone.render();
    }

    if (tbb.min.z() <= 0.0) {
        view_model_matrix = view_matrix * Geometry::assemble_transform(-offset, PI * Vec3d::UnitX(), scale);

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("normal_matrix", (Matrix3d)view_model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
        m_cone.render();
    }

    shader->stop_using();

#else 
    const BoundingBoxf3 box = bounding_box();

    Vec3d grabber_center = m_plane_center;
    grabber_center[Z] += float(box.radius()/2.0); // Margin
    rotate_vec3d_around_center(grabber_center, angles, m_plane_center);

    m_grabbers[0].center = grabber_center;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
    glsafe(::glLineWidth(m_hover_id == m_group_id ? 2.0f : 1.5f));

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();
#if ENABLE_GL_SHADERS_ATTRIBUTES
        const Camera& camera = wxGetApp().plater()->get_camera();
        shader->set_uniform("view_model_matrix", camera.get_view_matrix());
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        m_grabber_connection.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.color = ColorRGBA::YELLOW();
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex((Vec3f)m_plane_center.cast<float>());
        init_data.add_vertex((Vec3f)m_grabbers[0].center.cast<float>());

        // indices
        init_data.add_line(0, 1);

        m_grabber_connection.init_from(std::move(init_data));

        m_grabber_connection.render();

        shader->stop_using();
    }

#if !ENABLE_LEGACY_OPENGL_REMOVAL
    glsafe(::glColor3f(1.0, 1.0, 0.0));
    ::glBegin(GL_LINES);
    ::glVertex3dv(plane_center.data());
    ::glVertex3dv(m_grabbers[0].center.data());
    glsafe(::glEnd());
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

    shader = wxGetApp().get_shader("gouraud_light");
    if (shader != nullptr) {
        shader->start_using();
        shader->set_uniform("emission_factor", 0.1f);

        m_grabbers[0].render(m_hover_id == m_group_id, float((box.size().x() + box.size().y() + box.size().z()) / 3.0));

        shader->stop_using();
    }
#endif
}

bool GLGizmoCut3D::on_init()
{
    m_grabbers.emplace_back();
    m_shortcut_key = WXK_CONTROL_C;

    if(!m_rotation_gizmo.init())
        return false;

    return true;
}

std::string GLGizmoCut3D::on_get_name() const
{
    return _u8L("Cut");
}

void GLGizmoCut3D::on_set_state()
{
    if (get_state() == On)
        update_bb();

    m_rotation_gizmo.set_center(m_plane_center);
    m_rotation_gizmo.set_state(m_state);

    suppress_update_clipper_on_render = m_state != On;
}

void GLGizmoCut3D::on_set_hover_id() 
{
    m_rotation_gizmo.set_hover_id(m_hover_id < m_group_id ? m_hover_id: -1);
}

bool GLGizmoCut3D::on_is_activable() const
{
    // This is assumed in GLCanvas3D::do_rotate, do not change this
    // without updating that function too.
    return m_parent.get_selection().is_single_full_instance();
}

void GLGizmoCut3D::on_dragging(const UpdateData& data)
{
    if (m_hover_id < m_group_id)
        return;

    CutConnectors& connectors = m_c->selection_info()->model_object()->cut_connectors;

    if (m_hover_id == m_group_id) {
#if use_grabber_extension
        Vec3d starting_box_center = m_plane_center;
        starting_box_center[Z] -= 1.0; // some Margin
        rotate_vec3d_around_center(starting_box_center, m_rotation_gizmo.get_rotation(), m_plane_center);
#else
        const Vec3d& starting_box_center = m_plane_center;
#endif
        const Vec3d& starting_drag_position = m_grabbers[0].center;
        double projection = 0.0;

        Vec3d starting_vec = starting_drag_position - starting_box_center;
        if (starting_vec.norm() != 0.0) {
            Vec3d mouse_dir = data.mouse_ray.unit_vector();
            // finds the intersection of the mouse ray with the plane parallel to the camera viewport and passing throught the starting position
            // use ray-plane intersection see i.e. https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection algebric form
            // in our case plane normal and ray direction are the same (orthogonal view)
            // when moving to perspective camera the negative z unit axis of the camera needs to be transformed in world space and used as plane normal
            Vec3d inters = data.mouse_ray.a + (starting_drag_position - data.mouse_ray.a).dot(mouse_dir) / mouse_dir.squaredNorm() * mouse_dir;
            // vector from the starting position to the found intersection
            Vec3d inters_vec = inters - starting_drag_position;

            starting_vec.normalize();
            // finds projection of the vector along the staring direction
            projection = inters_vec.dot(starting_vec);
        }
        if (wxGetKeyState(WXK_SHIFT))
            projection = m_snap_step * (double)std::round(projection / m_snap_step);

        Vec3d shift = starting_vec * projection;

        // move  cut plane center
        set_center(m_plane_center + shift);
    }

    else if (m_hover_id > m_group_id && m_connector_mode == CutConnectorMode::Manual)
    {
        std::pair<Vec3d, Vec3d> pos_and_normal;
        if (!unproject_on_cut_plane(data.mouse_pos.cast<double>(), pos_and_normal))
            return;
        connectors[m_hover_id - m_connectors_group_id].pos    = pos_and_normal.first;
    }
}

void GLGizmoCut3D::set_center_pos(const Vec3d& center_pos)
{
    m_plane_center = center_pos;

    // !!! ysFIXME add smart clamp calculation
    // Clamp the center position of the cut plane to the object's bounding box
    //m_plane_center = Vec3d(std::clamp(center_pos.x(), m_min_pos.x(), m_max_pos.x()),
    //                       std::clamp(center_pos.y(), m_min_pos.y(), m_max_pos.y()),
    //                       std::clamp(center_pos.z(), m_min_pos.z(), m_max_pos.z()));

    m_center_offset = m_plane_center - m_bb_center;
}

BoundingBoxf3 GLGizmoCut3D::bounding_box() const
{
    BoundingBoxf3 ret;
    const Selection& selection = m_parent.get_selection();
    const Selection::IndicesList& idxs = selection.get_volume_idxs();
    for (unsigned int i : idxs) {
        const GLVolume* volume = selection.get_volume(i);
        // respect just to the solid parts for FFF and ignore pad and supports for SLA
        if (!volume->is_modifier && !volume->is_sla_pad() && !volume->is_sla_support())
            ret.merge(volume->transformed_convex_hull_bounding_box());
    }
    return ret;
}

BoundingBoxf3 GLGizmoCut3D::transformed_bounding_box() const
{
    // #ysFIXME !!!
    BoundingBoxf3 ret;
    const Selection& selection = m_parent.get_selection();
    const Selection::IndicesList& idxs = selection.get_volume_idxs();

    const int instance_idx = selection.get_instance_idx();
    const int object_idx = selection.get_object_idx();
    if (instance_idx < 0 || object_idx < 0)
        return ret;

    const Vec3d& instance_offset = wxGetApp().plater()->model().objects[object_idx]->instances[instance_idx]->get_offset();

    Vec3d cut_center_offset = m_plane_center - instance_offset;
    cut_center_offset[Z] -= selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z();

    const Vec3d& rotation = m_rotation_gizmo.get_rotation();
    const auto move  = Geometry::assemble_transform(-cut_center_offset, Vec3d::Zero(), Vec3d::Ones(), Vec3d::Ones() );
    const auto rot_z = Geometry::assemble_transform(Vec3d::Zero(), Vec3d(0, 0, -rotation.z()), Vec3d::Ones(), Vec3d::Ones());
    const auto rot_y = Geometry::assemble_transform(Vec3d::Zero(), Vec3d(0, -rotation.y(), 0), Vec3d::Ones(), Vec3d::Ones());
    const auto rot_x = Geometry::assemble_transform(Vec3d::Zero(), Vec3d(-rotation.x(), 0, 0), Vec3d::Ones(), Vec3d::Ones());

    const auto cut_matrix = rot_x * rot_y * rot_z * move;

    for (unsigned int i : idxs) {
        const GLVolume* volume = selection.get_volume(i);
        // respect just to the solid parts for FFF and ignore pad and supports for SLA
        if (!volume->is_modifier && !volume->is_sla_pad() && !volume->is_sla_support()) {

            const auto instance_matrix = Geometry::assemble_transform(
                Vec3d::Zero(),  // don't apply offset
                volume->get_instance_rotation().cwiseProduct(Vec3d(1.0, 1.0, 1.0)),
                volume->get_instance_scaling_factor(),
                volume->get_instance_mirror()
            );

            ret.merge(volume->transformed_convex_hull_bounding_box(cut_matrix * instance_matrix * volume->get_volume_transformation().get_matrix()));
        }
    }
    return ret;
}

bool GLGizmoCut3D::update_bb()
{
    const BoundingBoxf3 box = bounding_box();
    if (m_max_pos != box.max && m_min_pos != box.min) {
        m_max_pos = box.max;
        m_min_pos = box.min;
        m_bb_center = box.center();
        set_center_pos(m_bb_center + m_center_offset);

        m_plane.reset();
        m_cone.reset();
        m_sphere.reset();
        if (CommonGizmosDataObjects::SelectionInfo* selection = m_c->selection_info()) {
            m_selected.clear();
            m_selected.resize(selection->model_object()->cut_connectors.size(), false);
        }

        return true;
    }
    return false;
}

void GLGizmoCut3D::on_render()
{
    if (update_bb()) {
        m_rotation_gizmo.set_center(m_plane_center);
        update_clipper_on_render();
    }

    if (!m_cone.is_initialized())
        m_cone.init_from(its_make_cone(1.0, 1.0, double(PI) / 12.0));
    if (!m_sphere.is_initialized())
        m_sphere.init_from(its_make_sphere(1.0, double(PI) / 12.0));

    render_connectors(false);

    m_c->object_clipper()->render_cut();

    if (!m_hide_cut_plane) {
        render_cut_center_graber();
        render_cut_plane();
        if (m_mode == size_t(CutMode::cutPlanar)) {
            if (m_hover_id < m_group_id)
                m_rotation_gizmo.render();
        }
    }

    if (!suppress_update_clipper_on_render)
        update_clipper_on_render();
}

void GLGizmoCut3D::on_render_for_picking()
{
    m_rotation_gizmo.render_for_picking();
    render_grabbers_for_picking(m_parent.get_selection().get_bounding_box());

    render_connectors(true);
}

void GLGizmoCut3D::on_render_input_window(float x, float y, float bottom_limit)
{
    static float last_y = 0.0f;
    static float last_h = 0.0f;

    m_imgui->begin(_L("Cut"), ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

    m_imperial_units    = wxGetApp().app_config->get("use_inches") == "1";
    m_label_width       = m_imgui->get_style_scaling() * 100.0f;
    m_control_width     = m_imgui->get_style_scaling() * 150.0f;

    // adjust window position to avoid overlap the view toolbar
    const float win_h = ImGui::GetWindowHeight();
    y = std::min(y, bottom_limit - win_h);
    ImGui::SetWindowPos(ImVec2(x, y), ImGuiCond_Always);
    if (last_h != win_h || last_y != y) {
        // ask canvas for another frame to render the window in the correct position
        m_imgui->set_requires_extra_frame();
        if (last_h != win_h)
            last_h = win_h;
        if (last_y != y)
            last_y = y;
    }

    render_combo(_u8L("Mode"), m_modes, m_mode);

    bool revert_rotation{ false };
    bool revert_move{ false };

    CutConnectors& connectors = m_c->selection_info()->model_object()->cut_connectors;

    if (m_mode <= size_t(CutMode::cutByLine)) {
        ImGui::Separator();

        if (m_mode == size_t(CutMode::cutPlanar)) {
            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Move center"));

            m_imgui->disabled_begin(m_plane_center == bounding_box().center());
            revert_move = render_revert_button("move");
            m_imgui->disabled_end();

            for (Axis axis : {X, Y, Z})
                render_move_center_input(axis);
            m_imgui->text(m_imperial_units ? _L("in") : _L("mm")); 

            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Rotation"));

            m_imgui->disabled_begin(m_rotation_gizmo.get_rotation() == Vec3d::Zero());
            revert_rotation = render_revert_button("rotation");
            m_imgui->disabled_end();

            for (Axis axis : {X, Y, Z})
                render_rotation_input(axis);
            m_imgui->text(_L("°"));

            ImGui::Separator();

            double koef = m_imperial_units ? ObjectManipulation::mm_to_in : 1.0;
            wxString unit_str = " " + (m_imperial_units ? _L("in") : _L("mm"));

            const BoundingBoxf3 tbb = transformed_bounding_box();
            Vec3d tbb_sz = tbb.size();
            wxString size = "X: " + double_to_string(tbb.size().x() * koef,2) + unit_str +
                         ",  Y: " + double_to_string(tbb.size().y() * koef,2) + unit_str +
                         ",  Z: " + double_to_string(tbb.size().z() * koef,2) + unit_str ;

            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Build size"));
            ImGui::SameLine(m_label_width);
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, size);
#if 0
            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("DistToTop"));
            ImGui::SameLine(m_label_width);
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, double_to_string(tbb.max.z() * koef, 2));

            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("DistToBottom"));
            ImGui::SameLine(m_label_width);
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, double_to_string(tbb.min.z() * koef, 2));
#endif
        }
        else {
            ImGui::AlignTextToFramePadding();
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Connect some two points of object to cteate a cut plane"));
        }
    }

    m_imgui->disabled_begin(!m_keep_lower || !m_keep_upper);
    // Connectors section
    ImGui::Separator();

    ImGui::AlignTextToFramePadding();
    m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Connectors"));

    m_imgui->disabled_begin(connectors.empty());
    ImGui::SameLine(m_label_width);
    if (m_imgui->button("  " + _L("Reset") + "  ##connectors"))
        reset_connectors();
    m_imgui->disabled_end();

    m_imgui->text(_L("Mode"));
    render_connect_mode_radio_button(CutConnectorMode::Auto);
    render_connect_mode_radio_button(CutConnectorMode::Manual);

    m_imgui->text(_L("Type"));
    render_connect_type_radio_button(CutConnectorType::Plug);
    render_connect_type_radio_button(CutConnectorType::Dowel);

    if (render_combo(_u8L("Style"), m_connector_styles, m_connector_style))
        update_connector_shape();
    if (render_combo(_u8L("Shape"), m_connector_shapes, m_connector_shape_id))
        update_connector_shape();

    if (render_double_input(_u8L("Depth ratio"), m_connector_depth_ratio))
        for (auto& connector : connectors)
            connector.height = float(m_connector_depth_ratio);
    if (render_double_input(_u8L("Size"), m_connector_size))
        for (auto& connector : connectors)
            connector.radius = float(m_connector_size * 0.5);

    m_imgui->disabled_end();

    if (m_mode <= size_t(CutMode::cutByLine)) {
        ImGui::Separator();

        ImGui::AlignTextToFramePadding();
        m_imgui->text(_L("After cut") + ": ");
        bool keep = true;

        ImGui::SameLine(m_label_width);
        m_imgui->text(_L("Upper part"));
        ImGui::SameLine(2 * m_label_width);

        m_imgui->disabled_begin(!connectors.empty());
        m_imgui->checkbox(_L("Keep") + "##upper", connectors.empty() ? m_keep_upper : keep);
        m_imgui->disabled_end();
        ImGui::SameLine();
        m_imgui->disabled_begin(!m_keep_upper);
        m_imgui->checkbox(_L("Flip") + "##upper", m_rotate_upper);
        m_imgui->disabled_end();

        m_imgui->text("");
        ImGui::SameLine(m_label_width);
        m_imgui->text(_L("Lower part"));
        ImGui::SameLine(2 * m_label_width);

        m_imgui->disabled_begin(!connectors.empty());
        m_imgui->checkbox(_L("Keep") + "##lower", connectors.empty() ? m_keep_lower : keep);
        m_imgui->disabled_end();
        ImGui::SameLine();
        m_imgui->disabled_begin(!m_keep_lower);
        m_imgui->checkbox(_L("Flip") + "##lower", m_rotate_lower);
        m_imgui->disabled_end();
    }

    ImGui::Separator();

    m_imgui->disabled_begin((!m_keep_upper && !m_keep_lower) || !can_perform_cut());
    const bool cut_clicked = m_imgui->button(_L("Perform cut"));
    m_imgui->disabled_end();

    ImGui::Separator();

    m_imgui->checkbox(_L("Hide cut plane and grabbers"), m_hide_cut_plane);

    ////////
    static bool hide_clipped = true;
    static bool fill_cut = true;
    static float contour_width = 0.2f;
    if (m_imgui->checkbox("hide_clipped", hide_clipped) && !hide_clipped)
        m_clp_normal = m_c->object_clipper()->get_clipping_plane()->get_normal();
    m_imgui->checkbox("fill_cut", fill_cut);
    m_imgui->slider_float("contour_width", &contour_width, 0.f, 3.f);
    m_c->object_clipper()->set_behavior(hide_clipped, fill_cut, contour_width);
    ////////

    m_imgui->end();

    if (cut_clicked && (m_keep_upper || m_keep_lower))
        perform_cut(m_parent.get_selection());

    if (revert_move)
        set_center(bounding_box().center());
    if (revert_rotation) {
        m_rotation_gizmo.set_rotation(Vec3d::Zero());
        update_clipper();
    }
}

// get volume transformation regarding to the "border". Border is related from the siae of connectors
Transform3d GLGizmoCut3D::get_volume_transformation(const ModelVolume* volume) const
{
    bool is_prizm_dowel = m_connector_type == CutConnectorType::Dowel && m_connector_style == size_t(CutConnectorStyle::Prizm);
    const Transform3d connector_trafo = Geometry::assemble_transform(
        is_prizm_dowel ? Vec3d(0.0, 0.0, -m_connector_depth_ratio) : Vec3d::Zero(),
        m_rotation_gizmo.get_rotation(),
        Vec3d(0.5*m_connector_size, 0.5*m_connector_size, is_prizm_dowel ? 2 * m_connector_depth_ratio : m_connector_depth_ratio),
        Vec3d::Ones());
    const Vec3d connector_bb = m_connector_mesh.transformed_bounding_box(connector_trafo).size();

    const Vec3d bb = volume->mesh().bounding_box().size();

    // calculate an unused border - part of the the volume, where we can't put connectors
    const Vec3d border_scale(connector_bb.x() / bb.x(), connector_bb.y() / bb.y(), connector_bb.z() / bb.z());

    const Transform3d vol_matrix = volume->get_matrix();
    const Vec3d vol_trans = vol_matrix.translation();
    // offset of the volume will be changed after scaling, so calculate the needed offset and set it to a volume_trafo
    const Vec3d offset(vol_trans.x() * border_scale.x(), vol_trans.y() * border_scale.y(), vol_trans.z() * border_scale.z());

    // scale and translate volume to suppress to put connectors too close to the border
    return Geometry::assemble_transform(offset, Vec3d::Zero(), Vec3d::Ones() - border_scale, Vec3d::Ones()) * vol_matrix;
}

void GLGizmoCut3D::render_connectors(bool picking)
{
    m_has_invalid_connector = false;

    if (m_connector_mode == CutConnectorMode::Auto)
        return;

    const ModelObject* mo = m_c->selection_info()->model_object();
    const CutConnectors& connectors = mo->cut_connectors;
    if (connectors.size() != m_selected.size())
        return;

#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLShaderProgram* shader = picking ? wxGetApp().get_shader("flat") : wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    shader->start_using();

    ScopeGuard guard([shader]() { shader->stop_using(); });
#else
    GLShaderProgram* shader = picking ? nullptr : wxGetApp().get_shader("gouraud_light");
    if (shader)
        shader->start_using();
    ScopeGuard guard([shader]() { if (shader) shader->stop_using(); });
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_GL_SHADERS_ATTRIBUTES
    const Camera& camera = wxGetApp().plater()->get_camera();
#else
    glsafe(::glPushMatrix());
    glsafe(::glTranslated(0.0, 0.0, m_c->selection_info()->get_sla_shift()));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    ColorRGBA render_color;

    const ModelInstance* mi = mo->instances[m_c->selection_info()->get_active_instance()];
    const Vec3d& instance_offset = mi->get_offset();
    const float sla_shift        = m_c->selection_info()->get_sla_shift();

    const ClippingPlane* cp = m_c->object_clipper()->get_clipping_plane();
    const Vec3d& normal = cp && cp->is_active() ? cp->get_normal() : m_clp_normal;

    const Transform3d instance_trafo = Geometry::assemble_transform(Vec3d(0.0, 0.0, sla_shift), Vec3d::Zero(), Vec3d::Ones(), Vec3d::Ones()) * mi->get_transformation().get_matrix();

    for (size_t i = 0; i < connectors.size(); ++i) {
        const CutConnector& connector = connectors[i];
        const bool& point_selected = m_selected[i];

        double height = connector.height;
        // recalculate connector position to world position
        Vec3d pos = connector.pos + instance_offset;
        if (m_connector_type == CutConnectorType::Dowel &&
            m_connector_style == size_t(CutConnectorStyle::Prizm)) {
            pos -= height * normal;
            height *= 2;
        }
        pos[Z] += sla_shift;

        // First decide about the color of the point.
        if (picking)
            render_color = picking_decode(BASE_ID - i - m_connectors_group_id);
        else {
            if (size_t(m_hover_id- m_connectors_group_id) == i)
                render_color = ColorRGBA::CYAN();
            else { // neither hover nor picking
                int mesh_id = -1;
                for (const ModelVolume* mv : mo->volumes) {
                    ++mesh_id;
                    if (!mv->is_model_part())
                        continue;

                    const Transform3d volume_trafo = get_volume_transformation(mv);

                    if (m_c->raycaster()->raycasters()[mesh_id]->is_valid_intersection(pos, -normal, instance_trafo * volume_trafo)) {
                        render_color = ColorRGBA(1.0f, 1.0f, 1.0f, 0.5f);
                        break;
                    }
                    render_color = ColorRGBA(1.0f, 0.3f, 0.3f, 0.5f);
                }
                if (!m_has_invalid_connector && render_color == ColorRGBA(1.0f, 0.3f, 0.3f, 0.5f))
                    m_has_invalid_connector = true;
            }
        }

#if ENABLE_GL_SHADERS_ATTRIBUTES
        m_connector_shape.set_color(render_color);

        const Transform3d view_model_matrix = camera.get_view_matrix() * Geometry::assemble_transform(
            Vec3d(pos.x(), pos.y(), pos.z()),
            m_rotation_gizmo.get_rotation(),
            Vec3d(connector.radius, connector.radius, height),
            Vec3d::Ones()
        );
        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
#else
        const_cast<GLModel*>(&m_connector_shape)->set_color(-1, render_color);

        glsafe(::glPushMatrix());
        glsafe(::glTranslatef(pos.x(), pos.y(), pos.z()));

        const Vec3d& angles = m_rotation_gizmo.get_rotation();
        glsafe(::glRotated(Geometry::rad2deg(angles.z()), 0.0, 0.0, 1.0));
        glsafe(::glRotated(Geometry::rad2deg(angles.y()), 0.0, 1.0, 0.0));
        glsafe(::glRotated(Geometry::rad2deg(angles.x()), 1.0, 0.0, 0.0));

        glsafe(::glTranslated(0., 0., -0.5*connector.height));
        glsafe(::glScaled(connector.radius, connector.radius, height));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

        m_connector_shape.render();

#if !ENABLE_GL_SHADERS_ATTRIBUTES
        glsafe(::glPopMatrix());
#endif //!ENABLE_GL_SHADERS_ATTRIBUTES
    }

#if !ENABLE_GL_SHADERS_ATTRIBUTES
    glsafe(::glPopMatrix());
#endif //!ENABLE_GL_SHADERS_ATTRIBUTES
}

bool GLGizmoCut3D::can_perform_cut() const
{
    if (m_has_invalid_connector)
        return false;

    BoundingBoxf3 box   = bounding_box();
    double dist = (m_plane_center - box.center()).norm();
    if (dist > box.radius())
        return false;

    return true;
}

void GLGizmoCut3D::perform_cut(const Selection& selection)
{
    const int instance_idx = selection.get_instance_idx();
    const int object_idx = selection.get_object_idx();

    wxCHECK_RET(instance_idx >= 0 && object_idx >= 0, "GLGizmoCut: Invalid object selection");

    // m_cut_z is the distance from the bed. Subtract possible SLA elevation.
    const GLVolume* first_glvolume = selection.get_volume(*selection.get_volume_idxs().begin());
    const double object_cut_z = m_plane_center.z() - first_glvolume->get_sla_shift_z();

    const Vec3d& instance_offset = wxGetApp().plater()->model().objects[object_idx]->instances[instance_idx]->get_offset();

    Vec3d cut_center_offset = m_plane_center - instance_offset;
    cut_center_offset[Z] -= first_glvolume->get_sla_shift_z();

    bool create_dowels_as_separate_object = false;
    if (0.0 < object_cut_z && can_perform_cut()) {
        ModelObject* mo = wxGetApp().plater()->model().objects[object_idx];
        const bool has_connectors = !mo->cut_connectors.empty();
        // update connectors pos as offset of its center before cut performing
        if (has_connectors && m_connector_mode == CutConnectorMode::Manual) {
            for (CutConnector& connector : mo->cut_connectors) {
                connector.rotation = m_rotation_gizmo.get_rotation();

                if (m_connector_type == CutConnectorType::Dowel) {
                    if (m_connector_style == size_t(CutConnectorStyle::Prizm))
                        connector.height *= 2;
                }
                else {
                    // culculate shift of the connector center regarding to the position on the cut plane
                    Vec3d norm = m_grabbers[0].center - m_plane_center;
                    norm.normalize();
                    Vec3d shift = norm * (0.5 * connector.height);
                    connector.pos += shift;
                }
            }
            mo->apply_cut_connectors(_u8L("Connector"), CutConnectorAttributes(CutConnectorType(m_connector_type), CutConnectorStyle(m_connector_style), CutConnectorShape(m_connector_shape_id)));
            if (m_connector_type == CutConnectorType::Dowel)
                create_dowels_as_separate_object = true;
        }

        wxGetApp().plater()->cut(object_idx, instance_idx, cut_center_offset, m_rotation_gizmo.get_rotation(),
            only_if(has_connectors ? true : m_keep_upper, ModelObjectCutAttribute::KeepUpper) |
            only_if(has_connectors ? true : m_keep_lower, ModelObjectCutAttribute::KeepLower) |
            only_if(m_rotate_upper, ModelObjectCutAttribute::FlipUpper) | 
            only_if(m_rotate_lower, ModelObjectCutAttribute::FlipLower) | 
            only_if(create_dowels_as_separate_object, ModelObjectCutAttribute::CreateDowels));
        m_selected.clear();
    }
    else {
        // the object is SLA-elevated and the plane is under it.
    }
}



// Unprojects the mouse position on the mesh and saves hit point and normal of the facet into pos_and_normal
// Return false if no intersection was found, true otherwise.
bool GLGizmoCut3D::unproject_on_cut_plane(const Vec2d& mouse_position, std::pair<Vec3d, Vec3d>& pos_and_normal)
{
    const float sla_shift = m_c->selection_info()->get_sla_shift();

    const ModelObject* mo = m_c->selection_info()->model_object();
    const ModelInstance* mi = mo->instances[m_c->selection_info()->get_active_instance()];
    const Transform3d    instance_trafo = sla_shift > 0.0 ? 
        Geometry::assemble_transform(Vec3d(0.0, 0.0, sla_shift), Vec3d::Zero(), Vec3d::Ones(), Vec3d::Ones()) * mi->get_transformation().get_matrix() : mi->get_transformation().get_matrix();
    const Camera& camera = wxGetApp().plater()->get_camera();

    int mesh_id = -1;
    for (const ModelVolume* mv : mo->volumes) {
        ++mesh_id;
        if (!mv->is_model_part())
            continue;
        Vec3f hit;
        Vec3f normal;
        bool clipping_plane_was_hit = false;

        const Transform3d volume_trafo = get_volume_transformation(mv);

        m_c->raycaster()->raycasters()[mesh_id]->unproject_on_mesh(mouse_position, instance_trafo * volume_trafo,
            camera, hit, normal, m_c->object_clipper()->get_clipping_plane(),
            nullptr, &clipping_plane_was_hit);
        if (clipping_plane_was_hit) {
            // recalculate hit to object's local position
            Vec3d hit_d = hit.cast<double>();
            hit_d -= mi->get_offset();
            hit_d[Z] -= sla_shift;

            // Return both the point and the facet normal.
            pos_and_normal = std::make_pair(hit_d, normal.cast<double>());
            return true;
        }
    }
    return false;
}

void GLGizmoCut3D::reset_connectors()
{
    m_c->selection_info()->model_object()->cut_connectors.clear();
    update_model_object();
    m_selected.clear();
}

void GLGizmoCut3D::update_connector_shape()
{
    if (m_connector_shape.is_initialized())
        m_connector_shape.reset();

    const indexed_triangle_set its = ModelObject::get_connector_mesh({ m_connector_type, CutConnectorStyle(m_connector_style), CutConnectorShape(m_connector_shape_id) });
    m_connector_shape.init_from(its);

    m_connector_mesh.clear();
    m_connector_mesh = TriangleMesh(its);

}

void GLGizmoCut3D::update_model_object() const
{
    const ModelObjectPtrs& mos = wxGetApp().model().objects;
    ModelObject* mo = m_c->selection_info()->model_object();
    wxGetApp().obj_list()->update_info_items(std::find(mos.begin(), mos.end(), mo) - mos.begin());

    m_parent.post_event(SimpleEvent(EVT_GLCANVAS_SCHEDULE_BACKGROUND_PROCESS));
}

bool GLGizmoCut3D::gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down)
{
    if (is_dragging() || m_connector_mode == CutConnectorMode::Auto || (!m_keep_upper || !m_keep_lower))
        return false;

    CutConnectors& connectors = m_c->selection_info()->model_object()->cut_connectors;
    const Camera& camera = wxGetApp().plater()->get_camera();

    int mesh_id = -1;

    // left down without selection rectangle - place connector on the cut plane:
    if (action == SLAGizmoEventType::LeftDown && /*!m_selection_rectangle.is_dragging() && */!shift_down) {
        // If any point is in hover state, this should initiate its move - return control back to GLCanvas:
        if (m_hover_id != -1)
            return false;

        // If there is some selection, don't add new point and deselect everything instead.
        if (m_selection_empty) {
            std::pair<Vec3d, Vec3d> pos_and_normal;
            if (unproject_on_cut_plane(mouse_position.cast<double>(), pos_and_normal)) {
                const Vec3d& hit = pos_and_normal.first;
                const Vec3d& normal = pos_and_normal.second;
                // The clipping plane was clicked, hit containts coordinates of the hit in world coords.
                std::cout << hit.x() << "\t" << hit.y() << "\t" << hit.z() << std::endl;
                Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Add connector"));

                connectors.emplace_back(hit, m_rotation_gizmo.get_rotation(), float(m_connector_size * 0.5), float(m_connector_depth_ratio));
                update_model_object();
                m_selected.push_back(false);
                assert(m_selected.size() == connectors.size());
                m_parent.set_as_dirty();
                m_wait_for_up_event = true;

                return true;
            }
            return false;
        }
        return true;
    }
    else if (action == SLAGizmoEventType::RightDown && !shift_down) {
        // If any point is in hover state, this should initiate its move - return control back to GLCanvas:
        if (m_hover_id < m_connectors_group_id)
            return false;

        Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Delete connector"));

        size_t connector_id = m_hover_id - m_connectors_group_id;
        connectors.erase(connectors.begin() + connector_id);
        update_model_object();
        m_selected.erase(m_selected.begin() + connector_id);
        assert(m_selected.size() == connectors.size());
        m_parent.set_as_dirty();

        return true;
    }
    return false;
}

CommonGizmosDataID GLGizmoCut3D::on_get_requirements() const {
    return CommonGizmosDataID(
                int(CommonGizmosDataID::SelectionInfo)
              | int(CommonGizmosDataID::InstancesHider)
              | int(CommonGizmosDataID::ObjectClipper)
              | int(CommonGizmosDataID::Raycaster));
}

} // namespace GUI
} // namespace Slic3r
