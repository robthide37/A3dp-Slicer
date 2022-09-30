// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoCut.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"

#include <GL/glew.h>

#include <algorithm>

#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/Utils/UndoRedo.hpp"
#include "libslic3r/AppConfig.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"

#include "imgui/imgui_internal.h"

namespace Slic3r {
namespace GUI {

static const ColorRGBA GRABBER_COLOR = ColorRGBA::YELLOW();

// connector colors
static const ColorRGBA PLAG_COLOR           = ColorRGBA::YELLOW();
static const ColorRGBA DOWEL_COLOR          = ColorRGBA::DARK_YELLOW();
static const ColorRGBA HOVERED_PLAG_COLOR   = ColorRGBA::CYAN();
static const ColorRGBA HOVERED_DOWEL_COLOR  = ColorRGBA(0.0f, 0.5f, 0.5f, 1.0f);
static const ColorRGBA SELECTED_PLAG_COLOR  = ColorRGBA::GRAY();
static const ColorRGBA SELECTED_DOWEL_COLOR = ColorRGBA::DARK_GRAY();
static const ColorRGBA CONNECTOR_DEF_COLOR  = ColorRGBA(1.0f, 1.0f, 1.0f, 0.5f);

const unsigned int AngleResolution = 64;
const unsigned int ScaleStepsCount = 72;
const float ScaleStepRad = 2.0f * float(PI) / ScaleStepsCount;
const unsigned int ScaleLongEvery = 2;
const float ScaleLongTooth = 0.1f; // in percent of radius
const unsigned int SnapRegionsCount = 8;

using namespace Geometry;

// Generates mesh for a line
static GLModel::Geometry its_make_line(Vec3f beg_pos, Vec3f end_pos)
{
    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(2);
    init_data.reserve_indices(2);

    // vertices
    init_data.add_vertex(beg_pos);
    init_data.add_vertex(end_pos);

    // indices
    init_data.add_line(0, 1);
    return init_data;
}

// Generates mesh for a square plane
static GLModel::Geometry its_make_square_plane(float radius)
{
    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(4);
    init_data.reserve_indices(6);

    // vertices
    init_data.add_vertex(Vec3f(-radius, -radius, 0.0));
    init_data.add_vertex(Vec3f(radius , -radius, 0.0));
    init_data.add_vertex(Vec3f(radius , radius , 0.0));
    init_data.add_vertex(Vec3f(-radius, radius , 0.0));

    // indices
    init_data.add_triangle(0, 1, 2);
    init_data.add_triangle(2, 3, 0);
    return init_data;
}

//! -- #ysFIXME those functions bodies are ported from GizmoRotation
// Generates mesh for a circle 
static void init_from_circle(GLModel& model, double radius)
{
    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::LineLoop, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(ScaleStepsCount);
    init_data.reserve_indices(ScaleStepsCount);

    // vertices + indices
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        const float angle = float(i * ScaleStepRad);
        init_data.add_vertex(Vec3f(::cos(angle) * float(radius), ::sin(angle) * float(radius), 0.0f));
        init_data.add_index(i);
    }

    model.init_from(std::move(init_data));
    model.set_color(ColorRGBA::WHITE());
}

// Generates mesh for a scale
static void init_from_scale(GLModel& model, double radius)
{
    const float out_radius_long  = float(radius) * (1.0f + ScaleLongTooth);
    const float out_radius_short = float(radius) * (1.0f + 0.5f * ScaleLongTooth);

    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(2 * ScaleStepsCount);
    init_data.reserve_indices(2 * ScaleStepsCount);

    // vertices + indices
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        const float angle = float(i * ScaleStepRad);
        const float cosa = ::cos(angle);
        const float sina = ::sin(angle);
        const float in_x = cosa * float(radius);
        const float in_y = sina * float(radius);
        const float out_x = (i % ScaleLongEvery == 0) ? cosa * out_radius_long : cosa * out_radius_short;
        const float out_y = (i % ScaleLongEvery == 0) ? sina * out_radius_long : sina * out_radius_short;

        // vertices
        init_data.add_vertex(Vec3f(in_x, in_y, 0.0f));
        init_data.add_vertex(Vec3f(out_x, out_y, 0.0f));

        // indices
        init_data.add_line(i * 2, i * 2 + 1);
    }

    model.init_from(std::move(init_data));
    model.set_color(ColorRGBA::WHITE());
}

// Generates mesh for a snap_radii
static void init_from_snap_radii(GLModel& model, double radius)
{
    const float step = 2.0f * float(PI) / float(SnapRegionsCount);
    const float in_radius = float(radius) / 3.0f;
    const float out_radius = 2.0f * in_radius;

    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(2 * ScaleStepsCount);
    init_data.reserve_indices(2 * ScaleStepsCount);

    // vertices + indices
    for (unsigned int i = 0; i < ScaleStepsCount; ++i) {
        const float angle = float(i) * step;
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

    model.init_from(std::move(init_data));
    model.set_color(ColorRGBA::WHITE());
}

// Generates mesh for a angle_arc
static void init_from_angle_arc(GLModel& model, double angle, double radius)
{
    model.reset();

    const float step_angle = float(angle) / float(AngleResolution);
    const float ex_radius = float(radius);

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

    model.init_from(std::move(init_data));
}

//! --

GLGizmoCut3D::GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_connectors_group_id (3)
    , m_connector_type (CutConnectorType::Plug)
    , m_connector_style (size_t(CutConnectorStyle::Prizm))
    , m_connector_shape_id (size_t(CutConnectorShape::Circle))
{
    m_modes = { _u8L("Planar")//, _u8L("Grid")
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
    std::string tooltip;
    if (m_hover_id == Z) {
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
    if (tooltip.empty() && (m_hover_id == X || m_hover_id == Y)) {
        std::string axis = m_hover_id == X ? "X" : "Y";
        return axis + ": " + format(float(rad2deg(m_angle)), 1) + _u8L("°");
    }

    return tooltip;
}

bool GLGizmoCut3D::on_mouse(const wxMouseEvent &mouse_event)
{
    Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
    Vec2d mouse_pos = mouse_coord.cast<double>();

    if (mouse_event.ShiftDown() && mouse_event.LeftDown())
        return gizmo_event(SLAGizmoEventType::LeftDown, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), mouse_event.CmdDown());
    if (cut_line_processing()) {
        if (mouse_event.ShiftDown()) {
            if (mouse_event.Moving()|| mouse_event.Dragging())
                return gizmo_event(SLAGizmoEventType::Moving, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), mouse_event.CmdDown());
            if (mouse_event.LeftUp())
                return gizmo_event(SLAGizmoEventType::LeftUp, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), mouse_event.CmdDown());
        }
        discard_cut_line_processing();
    }
    else if (mouse_event.Moving())
        return false;

    if (use_grabbers(mouse_event)) {
        if (m_hover_id >= m_connectors_group_id) {
            if (mouse_event.LeftDown() && !mouse_event.CmdDown()&& !mouse_event.AltDown())
                unselect_all_connectors();
            if (mouse_event.LeftUp() && !mouse_event.ShiftDown())
                gizmo_event(SLAGizmoEventType::LeftUp, mouse_pos, mouse_event.ShiftDown(), mouse_event.AltDown(), mouse_event.CmdDown());
        }
        return true;
    }

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
    set_center(new_cut_center);
}

void GLGizmoCut3D::rotate_vec3d_around_plane_center(Vec3d&vec)
{
    vec = Transformation( assemble_transform(m_plane_center) * m_rotation_m * assemble_transform(-m_plane_center)).get_matrix() * vec;
}

void GLGizmoCut3D::put_connectors_on_cut_plane(const Vec3d& cp_normal, double cp_offset)
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
    BoundingBoxf3 box = bounding_box();

    Vec3d beg, end = beg = m_plane_center;
    beg[Z] = box.center().z() - m_radius;
    end[Z] = box.center().z() + m_radius;

    rotate_vec3d_around_plane_center(beg);
    rotate_vec3d_around_plane_center(end);

    double dist = (m_plane_center - beg).norm();

    // calculate normal and offset for clipping plane
    Vec3d normal = end - beg;
    if (normal == Vec3d::Zero())
        return;
    dist = std::clamp(dist, 0.0001, normal.norm());
    normal.normalize();
    const double offset = normal.dot(beg) + dist;

    m_c->object_clipper()->set_range_and_pos(normal, offset, dist);

    put_connectors_on_cut_plane(normal, offset);

    if (m_raycasters.empty())
        on_register_raycasters_for_picking();
    else
        update_raycasters_for_picking_transform();
}

void GLGizmoCut3D::update_clipper_on_render()
{
    update_clipper();
    force_update_clipper_on_render = false;
}

void GLGizmoCut3D::set_center(const Vec3d& center)
{
    set_center_pos(center);
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

    if (is_changed)
        update_connector_shape();

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

bool GLGizmoCut3D::render_slider_double_input(const std::string& label, double& value_in, int& tolerance_in)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(label);
    ImGui::SameLine(m_label_width);
    ImGui::PushItemWidth(m_control_width * 0.85f);

    float value = (float)value_in;
    if (m_imperial_units)
        value *= float(ObjectManipulation::mm_to_in);
    float old_val = value;

    const BoundingBoxf3 bbox = bounding_box();
    float mean_size = float((bbox.size().x() + bbox.size().y() + bbox.size().z()) / 9.0);
    float min_size = 1.f;
    if (m_imperial_units) {
        mean_size *= float(ObjectManipulation::mm_to_in);
        min_size  *= float(ObjectManipulation::mm_to_in);
    }
    std::string format = m_imperial_units ? "%.4f  " + _u8L("in") : "%.2f  " + _u8L("mm");

    m_imgui->slider_float(("##" + label).c_str(), &value, min_size, mean_size, format.c_str());
    value_in = (double)(value) * (m_imperial_units ? ObjectManipulation::in_to_mm : 1.0);

    ImGui::SameLine(m_label_width + m_control_width + 3);
    ImGui::PushItemWidth(m_control_width * 0.3f);

    float old_tolerance, tolerance = old_tolerance = (float)tolerance_in;
    m_imgui->slider_float(("##tolerance_" + label).c_str(), &tolerance, 1.f, 20.f, "%.f %%", 1.f, true, _L("Tolerance"));
    tolerance_in = (int)tolerance;

    return old_val != value || old_tolerance != tolerance;
}

void GLGizmoCut3D::render_move_center_input(int axis)
{
    m_imgui->text(m_axis_names[axis]+":");
    ImGui::SameLine();
    ImGui::PushItemWidth(0.3f*m_control_width);

    Vec3d move = m_plane_center;
    double in_val, value = in_val = move[axis];
    if (m_imperial_units)
        value *= ObjectManipulation::mm_to_in;
    ImGui::InputDouble(("##move_" + m_axis_names[axis]).c_str(), &value, 0.0, 0.0, "%.2f", ImGuiInputTextFlags_CharsDecimal);
    ImGui::SameLine();

    double val = value * (m_imperial_units ? ObjectManipulation::in_to_mm : 1.0);

    if (in_val != val) {
        move[axis] = val;
        set_center(move);
    }
}

bool GLGizmoCut3D::render_connect_type_radio_button(CutConnectorType type)
{
    ImGui::SameLine(type == CutConnectorType::Plug ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_types[size_t(type)], m_connector_type == type)) {
        m_connector_type = type;
        update_connector_shape();
        return true;
    }
    return false;
}

void GLGizmoCut3D::render_connect_mode_radio_button(CutConnectorMode mode)
{
    ImGui::SameLine(mode == CutConnectorMode::Auto ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_modes[int(mode)], m_connector_mode == mode))
        m_connector_mode = mode;
}

bool GLGizmoCut3D::render_reset_button(const std::string& label_id, const std::string& tooltip) const
{
    const ImGuiStyle& style = ImGui::GetStyle();

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, { 1, style.ItemSpacing.y });

    ImGui::PushStyleColor(ImGuiCol_Button, { 0.25f, 0.25f, 0.25f, 0.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, { 0.4f, 0.4f, 0.4f, 1.0f });
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, { 0.4f, 0.4f, 0.4f, 1.0f });

    std::string btn_label;
    btn_label += ImGui::RevertButton;
    const bool revert = ImGui::Button((btn_label +"##" + label_id).c_str());

    ImGui::PopStyleColor(3);

    if (ImGui::IsItemHovered())
        m_imgui->tooltip(tooltip.c_str(), ImGui::GetFontSize() * 20.0f);

    ImGui::PopStyleVar();

    return revert;
}

void GLGizmoCut3D::render_cut_plane()
{
    if (cut_line_processing())
        return;

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_CULL_FACE));
    glsafe(::glEnable(GL_BLEND));
    glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

    shader->start_using();

    const Camera& camera = wxGetApp().plater()->get_camera();
    const Transform3d view_model_matrix = camera.get_view_matrix() * translation_transform(m_plane_center) * m_rotation_m;

    shader->set_uniform("view_model_matrix", view_model_matrix);
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());

    if (can_perform_cut())
        m_plane.set_color({ 0.8f, 0.8f, 0.8f, 0.5f });
    else
        m_plane.set_color({ 1.0f, 0.8f, 0.8f, 0.5f });
    m_plane.render();

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}

static float get_grabber_mean_size(const BoundingBoxf3& bb)
{
    return float((bb.size().x() + bb.size().y() + bb.size().z()) / 3.0);
}

void GLGizmoCut3D::render_model(GLModel& model, const ColorRGBA& color, Transform3d view_model_matrix)
{
    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader) {
        shader->start_using();

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", wxGetApp().plater()->get_camera().get_projection_matrix());

        model.set_color(color);
        model.render();

        shader->stop_using();
    }
}

void GLGizmoCut3D::render_line(GLModel& line_model, const ColorRGBA& color, Transform3d view_model_matrix, float width)
{
    GLShaderProgram* shader = OpenGLManager::get_gl_info().is_core_profile() ? wxGetApp().get_shader("dashed_thick_lines") : wxGetApp().get_shader("flat");
    if (shader) {
        shader->start_using();

        const Camera& camera = wxGetApp().plater()->get_camera();
        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", wxGetApp().plater()->get_camera().get_projection_matrix());
        shader->set_uniform("width", width);

        line_model.set_color(color);
        line_model.render();

        shader->stop_using();
    }
}

void GLGizmoCut3D::render_rotation_snapping(Axis axis, const ColorRGBA& color)
{
    GLShaderProgram* line_shader = OpenGLManager::get_gl_info().is_core_profile() ? wxGetApp().get_shader("dashed_thick_lines") : wxGetApp().get_shader("flat");
    if (!line_shader)
        return;

    const Camera& camera = wxGetApp().plater()->get_camera();
    Transform3d view_model_matrix = camera.get_view_matrix() * translation_transform(m_plane_center) * m_start_dragging_m;

    if (axis == X)
        view_model_matrix = view_model_matrix * rotation_transform(0.5 * PI * Vec3d::UnitY()) * rotation_transform(-PI * Vec3d::UnitZ());
    else
        view_model_matrix = view_model_matrix * rotation_transform(-0.5 * PI * Vec3d::UnitZ()) * rotation_transform(-0.5 * PI * Vec3d::UnitY());

    line_shader->start_using();
    line_shader->set_uniform("projection_matrix", camera.get_projection_matrix());
    line_shader->set_uniform("view_model_matrix", view_model_matrix);
    line_shader->set_uniform("width", 0.25f);

    m_circle.render();
    m_scale.render();
    m_snap_radii.render();
    m_reference_radius.render();
    if (m_dragging) {
        line_shader->set_uniform("width", 1.5f);
        m_angle_arc.set_color(color);
        m_angle_arc.render();
    }

    line_shader->stop_using();
}

void GLGizmoCut3D::render_grabber_connection(const ColorRGBA& color, Transform3d view_matrix)
{
    const Transform3d line_view_matrix = view_matrix * scale_transform(Vec3d(1.0, 1.0, m_grabber_connection_len));

    render_line(m_grabber_connection, color, line_view_matrix, 0.2f);
};

void GLGizmoCut3D::render_cut_plane_grabbers()
{
    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    ColorRGBA color = m_hover_id == Z ? complementary(GRABBER_COLOR) : GRABBER_COLOR;

    const Transform3d view_matrix = wxGetApp().plater()->get_camera().get_view_matrix() * translation_transform(m_plane_center) * m_rotation_m;

    const Grabber& grabber = m_grabbers.front();
    const float mean_size = get_grabber_mean_size(bounding_box());

    double size = m_dragging && m_hover_id == Z ? double(grabber.get_dragging_half_size(mean_size)) : double(grabber.get_half_size(mean_size));

    Vec3d cone_scale = Vec3d(0.75 * size, 0.75 * size, 1.8 * size);
    Vec3d offset = 1.25 * size * Vec3d::UnitZ();

    // render Z grabber

    if ((!m_dragging && m_hover_id < 0))
        render_grabber_connection(color, view_matrix);
    render_model(m_sphere.model, color, view_matrix * scale_transform(size));

    if (!m_dragging && m_hover_id < 0 || m_hover_id == Z)
    {
        const BoundingBoxf3 tbb = transformed_bounding_box();
        if (tbb.min.z() <= 0.0)
            render_model(m_cone.model, color, view_matrix * assemble_transform(-offset, PI * Vec3d::UnitX(), cone_scale));

        if (tbb.max.z() >= 0.0)
            render_model(m_cone.model, color, view_matrix * assemble_transform(offset, Vec3d::Zero(), cone_scale));
    }

    // render top sphere for X/Y grabbers

    if ((!m_dragging && m_hover_id < 0) || m_hover_id == X || m_hover_id == Y)
    {
        size = m_dragging ? double(grabber.get_dragging_half_size(mean_size)) : double(grabber.get_half_size(mean_size));
        color = m_hover_id == Y ? complementary(ColorRGBA::GREEN()) :
                m_hover_id == X ? complementary(ColorRGBA::RED())   : ColorRGBA::GRAY();
        render_model(m_sphere.model, color, view_matrix * assemble_transform(m_grabber_connection_len * Vec3d::UnitZ(), Vec3d::Zero(), size * Vec3d::Ones()));
    }

    // render X grabber

    if ((!m_dragging && m_hover_id < 0) || m_hover_id == X)
    {
        size = m_dragging && m_hover_id == X ? double(grabber.get_dragging_half_size(mean_size)) : double(grabber.get_half_size(mean_size));
        cone_scale = Vec3d(0.75 * size, 0.75 * size, 1.8 * size);
        color = m_hover_id == X ? complementary(ColorRGBA::RED()) : ColorRGBA::RED();

        if (m_hover_id == X) {
            render_grabber_connection(color, view_matrix);
            render_rotation_snapping(X, color);
        }

        offset = Vec3d(0.0, 1.25 * size, m_grabber_connection_len);
        render_model(m_cone.model, color, view_matrix * assemble_transform(offset, -0.5 * PI * Vec3d::UnitX(), cone_scale));
        offset = Vec3d(0.0, -1.25 * size, m_grabber_connection_len);
        render_model(m_cone.model, color, view_matrix * assemble_transform(offset, 0.5 * PI * Vec3d::UnitX(), cone_scale));
    }

    // render Y grabber

    if ((!m_dragging && m_hover_id < 0) || m_hover_id == Y)
    {
        size = m_dragging && m_hover_id == Y ? double(grabber.get_dragging_half_size(mean_size)) : double(grabber.get_half_size(mean_size));
        cone_scale = Vec3d(0.75 * size, 0.75 * size, 1.8 * size);
        color = m_hover_id == Y ? complementary(ColorRGBA::GREEN()) : ColorRGBA::GREEN();

        if (m_hover_id == Y) {
            render_grabber_connection(color, view_matrix);
            render_rotation_snapping(Y, color);
        }

        offset = Vec3d(1.25 * size, 0.0, m_grabber_connection_len);
        render_model(m_cone.model, color, view_matrix * assemble_transform(offset, 0.5 * PI * Vec3d::UnitY(), cone_scale));
        offset = Vec3d(-1.25 * size, 0.0, m_grabber_connection_len);
        render_model(m_cone.model, color, view_matrix * assemble_transform(offset, -0.5 * PI * Vec3d::UnitY(), cone_scale));
    }
}

void GLGizmoCut3D::render_cut_line()
{
    if (!cut_line_processing() || m_line_end == Vec3d::Zero())
        return;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    m_cut_line.reset();
    m_cut_line.init_from(its_make_line((Vec3f)m_line_beg.cast<float>(), (Vec3f)m_line_end.cast<float>()));

    render_line(m_cut_line, GRABBER_COLOR, wxGetApp().plater()->get_camera().get_view_matrix(), 0.25f);
}

bool GLGizmoCut3D::on_init()
{
    m_grabbers.emplace_back();
    m_shortcut_key = WXK_CONTROL_C;

    // initiate info shortcuts
    const wxString ctrl = GUI::shortkey_ctrl_prefix();
    const wxString alt  = GUI::shortkey_alt_prefix();

    m_shortcuts.push_back(std::make_pair(_L("Left click"),        _L("Add connector")));
    m_shortcuts.push_back(std::make_pair(_L("Right click"),       _L("Remove connector")));
    m_shortcuts.push_back(std::make_pair(_L("Drag"),              _L("Move connector")));
    m_shortcuts.push_back(std::make_pair(ctrl + _L("Left click"), _L("Add connector to selection")));
    m_shortcuts.push_back(std::make_pair(alt + _L("Left click"),  _L("Remove connector from selection")));
    m_shortcuts.push_back(std::make_pair(ctrl + "A",              _L("Select all connectors")));

    return true;
}

void GLGizmoCut3D::on_load(cereal::BinaryInputArchive& ar)
{
    ar( m_keep_upper, m_keep_lower, m_rotate_lower, m_rotate_upper, m_hide_cut_plane, m_mode, m_connectors_editing,//m_selected,
 //       m_connector_depth_ratio, m_connector_size, m_connector_mode, m_connector_type, m_connector_style, m_connector_shape_id,
        m_ar_plane_center, m_rotation_m);

    set_center_pos(m_ar_plane_center, true);

    force_update_clipper_on_render = true;

    m_parent.request_extra_frame();
}

void GLGizmoCut3D::on_save(cereal::BinaryOutputArchive& ar) const
{ 
    ar( m_keep_upper, m_keep_lower, m_rotate_lower, m_rotate_upper, m_hide_cut_plane, m_mode, m_connectors_editing,//m_selected,
 //       m_connector_depth_ratio, m_connector_size, m_connector_mode, m_connector_type, m_connector_style, m_connector_shape_id,
        m_ar_plane_center, m_start_dragging_m);
}

std::string GLGizmoCut3D::on_get_name() const
{
    return _u8L("Cut");
}

void GLGizmoCut3D::on_set_state()
{
    if (m_state == On) {
        update_bb();

        // initiate archived values
        m_ar_plane_center   = m_plane_center;
        m_start_dragging_m  = m_rotation_m;

        m_parent.request_extra_frame();
    }
    else
        m_c->object_clipper()->release();

	force_update_clipper_on_render = m_state == On;
}

void GLGizmoCut3D::on_register_raycasters_for_picking()
{
    assert(m_raycasters.empty());
    set_volumes_picking_state(false);

    init_picking_models();

    if (m_connectors_editing) {
        if (CommonGizmosDataObjects::SelectionInfo* si = m_c->selection_info()) {
            const CutConnectors& connectors = si->model_object()->cut_connectors;
            for (int i = 0; i < int(connectors.size()); ++i)
                m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, i + m_connectors_group_id, *(m_shapes[connectors[i].attribs]).mesh_raycaster, Transform3d::Identity()));
        }
    }
    else {
        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, X, *m_cone.mesh_raycaster, Transform3d::Identity()));
        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, X, *m_cone.mesh_raycaster, Transform3d::Identity()));

        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, Y, *m_cone.mesh_raycaster, Transform3d::Identity()));
        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, Y, *m_cone.mesh_raycaster, Transform3d::Identity()));

        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, Z, *m_sphere.mesh_raycaster, Transform3d::Identity()));
        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, Z, *m_cone.mesh_raycaster, Transform3d::Identity()));
        m_raycasters.emplace_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, Z, *m_cone.mesh_raycaster, Transform3d::Identity()));
    }

    update_raycasters_for_picking_transform();
}

void GLGizmoCut3D::on_unregister_raycasters_for_picking()
{
    m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo);
    m_raycasters.clear();
    set_volumes_picking_state(true);
}

void GLGizmoCut3D::update_raycasters_for_picking()
{
    on_unregister_raycasters_for_picking();
    on_register_raycasters_for_picking();
}

void GLGizmoCut3D::set_volumes_picking_state(bool state)
{
    std::vector<std::shared_ptr<SceneRaycasterItem>>* raycasters = m_parent.get_raycasters_for_picking(SceneRaycaster::EType::Volume);
    if (raycasters != nullptr) {
        const Selection& selection = m_parent.get_selection();
        const Selection::IndicesList ids = selection.get_volume_idxs();
        for (unsigned int id : ids) {
            const GLVolume* v = selection.get_volume(id);
            auto it = std::find_if(raycasters->begin(), raycasters->end(), [v](std::shared_ptr<SceneRaycasterItem> item) { return item->get_raycaster() == v->mesh_raycaster.get(); });
            if (it != raycasters->end())
                (*it)->set_active(state);
        }
    }
}

void GLGizmoCut3D::update_raycasters_for_picking_transform()
{
    if (m_connectors_editing) {
        CommonGizmosDataObjects::SelectionInfo* si = m_c->selection_info();
        if (!si) 
            return;
        const ModelObject* mo = si->model_object();
        const CutConnectors& connectors = mo->cut_connectors;
        if (connectors.empty())
            return;
        auto inst_id = m_c->selection_info()->get_active_instance();
        if (inst_id < 0)
            return;

        const Vec3d& instance_offset = mo->instances[inst_id]->get_offset();
        const double sla_shift = double(m_c->selection_info()->get_sla_shift());

        const ClippingPlane* cp = m_c->object_clipper()->get_clipping_plane();
        const Vec3d& normal = cp && cp->is_active() ? cp->get_normal() : m_clp_normal;

        for (size_t i = 0; i < connectors.size(); ++i) {
            const CutConnector& connector = connectors[i];

            float height = connector.height;
            // recalculate connector position to world position
            Vec3d pos = connector.pos + instance_offset;
            if (connector.attribs.type == CutConnectorType::Dowel &&
                connector.attribs.style == CutConnectorStyle::Prizm) {
                pos -= height * normal;
                height *= 2;
            }
            pos[Z] += sla_shift;

            const Transform3d scale_trafo = scale_transform(Vec3f(connector.radius, connector.radius, height).cast<double>());
            m_raycasters[i]->set_transform(translation_transform(pos) * m_rotation_m * scale_trafo);
        }
    }
    else {
        const Transform3d trafo = translation_transform(m_plane_center) * m_rotation_m;

        const BoundingBoxf3 box = bounding_box();
        const float mean_size = get_grabber_mean_size(box);

        double size = double(m_grabbers.front().get_half_size(mean_size));
        Vec3d scale = Vec3d(0.75 * size, 0.75 * size, 1.8 * size);

        Vec3d offset = Vec3d(0.0, 1.25 * size, m_grabber_connection_len);
        m_raycasters[0]->set_transform(trafo * assemble_transform(offset, -0.5 * PI * Vec3d::UnitX(), scale));
        offset = Vec3d(0.0, -1.25 * size, m_grabber_connection_len);
        m_raycasters[1]->set_transform(trafo * assemble_transform(offset, 0.5 * PI * Vec3d::UnitX(), scale));

        offset = Vec3d(1.25 * size, 0.0, m_grabber_connection_len);
        m_raycasters[2]->set_transform(trafo * assemble_transform(offset, 0.5 * PI * Vec3d::UnitY(), scale));
        offset = Vec3d(-1.25 * size, 0.0, m_grabber_connection_len);
        m_raycasters[3]->set_transform(trafo * assemble_transform(offset, -0.5 * PI * Vec3d::UnitY(), scale));

        offset = 1.25 * size * Vec3d::UnitZ();
        m_raycasters[4]->set_transform(trafo * scale_transform(size));
        m_raycasters[5]->set_transform(trafo * assemble_transform(-offset, PI * Vec3d::UnitX(), scale));
        m_raycasters[6]->set_transform(trafo * assemble_transform(offset, Vec3d::Zero(), scale));
    }
}

void GLGizmoCut3D::on_set_hover_id() 
{
}

bool GLGizmoCut3D::on_is_activable() const
{
    // This is assumed in GLCanvas3D::do_rotate, do not change this
    // without updating that function too.
    return m_parent.get_selection().is_single_full_instance();
}

Vec3d GLGizmoCut3D::mouse_position_in_local_plane(Axis axis, const Linef3& mouse_ray) const
{
    double half_pi = 0.5 * PI;

    Transform3d m = Transform3d::Identity();

    switch (axis)
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

    m = m * m_start_dragging_m.inverse();
    m.translate(-m_plane_center);

    return transform(mouse_ray, m).intersect_plane(0.0);
}

void GLGizmoCut3D::dragging_grabber_z(const GLGizmoBase::UpdateData &data)
{
    Vec3d starting_box_center = m_plane_center - Vec3d::UnitZ(); // some Margin
    rotate_vec3d_around_plane_center(starting_box_center);

    const Vec3d&starting_drag_position = m_plane_center;
    double      projection             = 0.0;

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

    const Vec3d shift = starting_vec * projection;

    // move  cut plane center
    set_center(m_plane_center + shift);
}

void GLGizmoCut3D::dragging_grabber_xy(const GLGizmoBase::UpdateData &data)
{
    const Vec2d mouse_pos = to_2d(mouse_position_in_local_plane((Axis)m_hover_id, data.mouse_ray));

    const Vec2d orig_dir = Vec2d::UnitX();
    const Vec2d new_dir  = mouse_pos.normalized();

    const double two_pi = 2.0 * PI;

    double theta = ::acos(std::clamp(new_dir.dot(orig_dir), -1.0, 1.0));
    if (cross2(orig_dir, new_dir) < 0.0)
        theta = two_pi - theta;

    const double len = mouse_pos.norm();
    // snap to coarse snap region
    if (m_snap_coarse_in_radius <= len && len <= m_snap_coarse_out_radius) {
        const double step = two_pi / double(SnapRegionsCount);
        theta             = step * std::round(theta / step);
    }
    // snap to fine snap region (scale)
    else if (m_snap_fine_in_radius <= len && len <= m_snap_fine_out_radius) {
        const double step = two_pi / double(ScaleStepsCount);
        theta             = step * std::round(theta / step);
    }

    if (is_approx(theta, two_pi))
        theta = 0.0;
    if (m_hover_id == X)
        theta += 0.5 * PI;

    Vec3d rotation = Vec3d::Zero();
    rotation[m_hover_id] = theta;
    m_rotation_m         = m_start_dragging_m * rotation_transform(rotation);

    m_angle = theta;
    while (m_angle > two_pi)
        m_angle -= two_pi;
    if (m_angle < 0.0)
        m_angle += two_pi;

    update_clipper();
}

void GLGizmoCut3D::dragging_connector(const GLGizmoBase::UpdateData &data)
{
    CutConnectors&          connectors = m_c->selection_info()->model_object()->cut_connectors;
    std::pair<Vec3d, Vec3d> pos_and_normal;
    Vec3d                   pos_world;

    if (unproject_on_cut_plane(data.mouse_pos.cast<double>(), pos_and_normal, pos_world)) {
        connectors[m_hover_id - m_connectors_group_id].pos = pos_and_normal.first;
        update_raycasters_for_picking_transform();
    }
}

void GLGizmoCut3D::on_dragging(const UpdateData& data)
{
    if (m_hover_id < 0)
        return;
    if (m_hover_id == Z)
        dragging_grabber_z(data);
    else if (m_hover_id == X || m_hover_id == Y)
        dragging_grabber_xy(data);
    else if (m_hover_id >= m_connectors_group_id && m_connector_mode == CutConnectorMode::Manual)
        dragging_connector(data);
}

void GLGizmoCut3D::on_start_dragging()
{
    m_angle = 0.0;
    if (m_hover_id >= m_connectors_group_id && m_connector_mode == CutConnectorMode::Manual)
        Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Move connector"), UndoRedo::SnapshotType::GizmoAction);

    if (m_hover_id == X || m_hover_id == Y)
        m_start_dragging_m = m_rotation_m;
}

void GLGizmoCut3D::on_stop_dragging()
{
    if (m_hover_id == X || m_hover_id == Y) {
        m_angle_arc.reset();
        m_angle = 0.0;
        Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Rotate cut plane"), UndoRedo::SnapshotType::GizmoAction);
        m_start_dragging_m = m_rotation_m;
    }
    else if (m_hover_id == Z) {
        Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Move cut plane"), UndoRedo::SnapshotType::GizmoAction);
        m_ar_plane_center = m_plane_center;
    }
}

void GLGizmoCut3D::set_center_pos(const Vec3d& center_pos, bool force/* = false*/)
{
    if (force || transformed_bounding_box(true).contains(center_pos)) {
        m_plane_center = center_pos;
        m_center_offset = m_plane_center - m_bb_center;
    }
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

BoundingBoxf3 GLGizmoCut3D::transformed_bounding_box(bool revert_move /*= false*/) const
{
    // #ysFIXME !!!
    BoundingBoxf3 ret;

    const ModelObject* mo = m_c->selection_info()->model_object();
    if (!mo)
        return ret;
    const int instance_idx = m_c->selection_info()->get_active_instance();
    if (instance_idx < 0)
        return ret;
    const ModelInstance* mi = mo->instances[instance_idx];

    const Vec3d& instance_offset = mi->get_offset();
    Vec3d cut_center_offset = m_plane_center - instance_offset;
    cut_center_offset[Z] -= m_c->selection_info()->get_sla_shift();

    const auto move  = assemble_transform(-cut_center_offset);
    const auto move2 = assemble_transform(m_plane_center);

    const auto cut_matrix = (revert_move ? move2 : Transform3d::Identity()) * m_rotation_m.inverse() * move;

    const Selection& selection = m_parent.get_selection();
    const Selection::IndicesList& idxs = selection.get_volume_idxs();
    for (unsigned int i : idxs) {
        const GLVolume* volume = selection.get_volume(i);
        // respect just to the solid parts for FFF and ignore pad and supports for SLA
        if (!volume->is_modifier && !volume->is_sla_pad() && !volume->is_sla_support()) {

            const auto instance_matrix = assemble_transform(
                Vec3d::Zero(),  // don't apply offset
                volume->get_instance_rotation().cwiseProduct(Vec3d(1.0, 1.0, 1.0)),
                volume->get_instance_scaling_factor(),
                volume->get_instance_mirror()
            );

            auto volume_trafo = instance_matrix * volume->get_volume_transformation().get_matrix();

            ret.merge(volume->transformed_convex_hull_bounding_box(cut_matrix * volume_trafo));
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
        if (box.contains(m_center_offset))
            set_center_pos(m_bb_center + m_center_offset, true);
        else
            set_center_pos(m_bb_center, true);

        m_radius = box.radius();
        m_grabber_connection_len = 0.75 * m_radius;// std::min<double>(0.75 * m_radius, 35.0);
        m_grabber_radius = m_grabber_connection_len * 0.85;

        m_snap_coarse_in_radius   = m_grabber_radius / 3.0;
        m_snap_coarse_out_radius  = m_snap_coarse_in_radius * 2.;
        m_snap_fine_in_radius     = m_grabber_connection_len * 0.85;
        m_snap_fine_out_radius    = m_grabber_connection_len * 1.15;

        m_plane.reset();
        m_cone.reset();
        m_sphere.reset();
        m_grabber_connection.reset();
        m_circle.reset();
        m_scale.reset();
        m_snap_radii.reset();
        m_reference_radius.reset();

        on_unregister_raycasters_for_picking();

        if (CommonGizmosDataObjects::SelectionInfo* selection = m_c->selection_info()) {
            m_selected.clear();
            m_selected.resize(selection->model_object()->cut_connectors.size(), false);
        }

        return true;
    }
    return false;
}

void GLGizmoCut3D::init_picking_models()
{
    if (!m_cone.model.is_initialized()) {
        indexed_triangle_set its = its_make_cone(1.0, 1.0, PI / 12.0);
        m_cone.model.init_from(its);
        m_cone.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(its)));
    }
    if (!m_sphere.model.is_initialized()) {
        indexed_triangle_set its = its_make_sphere(1.0, PI / 12.0);
        m_sphere.model.init_from(its);
        m_sphere.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(its)));
    }
}

void GLGizmoCut3D::init_rendering_items()
{
    if (!m_grabber_connection.is_initialized())
        m_grabber_connection.init_from(its_make_line(Vec3f::Zero(), Vec3f::UnitZ()));
    if (!m_circle.is_initialized())
        init_from_circle(m_circle, m_grabber_radius);
    if (!m_scale.is_initialized())
        init_from_scale(m_scale, m_grabber_radius);
    if (!m_snap_radii.is_initialized())
        init_from_snap_radii(m_snap_radii, m_grabber_radius);
    if (!m_reference_radius.is_initialized()) {
        m_reference_radius.init_from(its_make_line(Vec3f::Zero(), m_grabber_connection_len * Vec3f::UnitX()));
        m_reference_radius.set_color(ColorRGBA::WHITE());
    }
    if (!m_angle_arc.is_initialized() || m_angle != 0.0)
        init_from_angle_arc(m_angle_arc, m_angle, m_grabber_connection_len);

    if (!m_plane.is_initialized() && !m_hide_cut_plane && !m_connectors_editing)
        m_plane.init_from(its_make_square_plane(float(m_radius)));
}

void GLGizmoCut3D::render_clipper_cut()
{
    if (! m_connectors_editing)
        ::glDisable(GL_DEPTH_TEST);
    m_c->object_clipper()->render_cut();
    if (! m_connectors_editing)
        ::glEnable(GL_DEPTH_TEST);
}

void GLGizmoCut3D::on_render()
{
    if (update_bb() || force_update_clipper_on_render) {
        update_clipper_on_render();
        m_c->object_clipper()->set_behavior(m_connectors_editing, m_connectors_editing, 0.4);
    }

    init_picking_models();

    init_rendering_items();

    render_connectors();

    render_clipper_cut();

    if (!m_hide_cut_plane && !m_connectors_editing) {
        render_cut_plane();
        render_cut_plane_grabbers();
    }

    render_cut_line();
}

void GLGizmoCut3D::render_debug_input_window()
{
    m_imgui->begin(wxString("DEBUG"));

    static bool  hide_clipped  = false;
    static bool  fill_cut      = false;
    static float contour_width = 0.4f;

    m_imgui->checkbox(_L("Hide cut plane and grabbers"), m_hide_cut_plane);
    if (m_imgui->checkbox("hide_clipped", hide_clipped) && !hide_clipped)
        m_clp_normal = m_c->object_clipper()->get_clipping_plane()->get_normal();
    m_imgui->checkbox("fill_cut", fill_cut);
    m_imgui->slider_float("contour_width", &contour_width, 0.f, 3.f);
    if (auto oc = m_c->object_clipper())
        oc->set_behavior(hide_clipped || m_connectors_editing, fill_cut || m_connectors_editing, double(contour_width));

    m_imgui->end();
}

void GLGizmoCut3D::adjust_window_position(float x, float y, float bottom_limit)
{
    static float last_y = 0.0f;
    static float last_h = 0.0f;

    const float win_h = ImGui::GetWindowHeight();
    y                 = std::min(y, bottom_limit - win_h);

    ImGui::SetWindowPos(ImVec2(x, y), ImGuiCond_Always);

    if (last_h != win_h || last_y != y) {
        // ask canvas for another frame to render the window in the correct position
        m_imgui->set_requires_extra_frame();
        if (last_h != win_h)
            last_h = win_h;
        if (last_y != y)
            last_y = y;
    }
}

void GLGizmoCut3D::unselect_all_connectors()
{
    std::fill(m_selected.begin(), m_selected.end(), false);
    m_selected_count = 0;
}

void GLGizmoCut3D::select_all_connectors()
{
    std::fill(m_selected.begin(), m_selected.end(), true);
    m_selected_count = int(m_selected.size());
}

void GLGizmoCut3D::render_shortcuts()
{
    if (m_imgui->button("? " + (m_show_shortcuts ? wxString(ImGui::CollapseBtn) : wxString(ImGui::ExpandBtn))))
        m_show_shortcuts = !m_show_shortcuts;

    if (m_show_shortcuts)
        for (const auto&shortcut : m_shortcuts ){
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, shortcut.first);
            ImGui::SameLine(m_label_width);
            m_imgui->text(shortcut.second);
        }
}

void GLGizmoCut3D::apply_selected_connectors(std::function<void(size_t idx)> apply_fn)
{
    for (size_t idx = 0; idx < m_selected.size(); idx++)
        if (m_selected[idx])
            apply_fn(idx);
}

void GLGizmoCut3D::render_connectors_input_window(CutConnectors &connectors)
{
    // add shortcuts panel
    render_shortcuts();

    // Connectors section

    ImGui::Separator();

    // WIP : Auto : Need to implement
    // m_imgui->text(_L("Mode"));
    // render_connect_mode_radio_button(CutConnectorMode::Auto);
    // render_connect_mode_radio_button(CutConnectorMode::Manual);

    ImGui::AlignTextToFramePadding();
    m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Connectors"));

    m_imgui->disabled_begin(connectors.empty());
    ImGui::SameLine(m_label_width);
    if (render_reset_button("connectors", _u8L("Remove connectors")))
        reset_connectors();
    m_imgui->disabled_end();

    m_imgui->text(_L("Type"));
    bool type_changed = render_connect_type_radio_button(CutConnectorType::Plug);
    type_changed     |= render_connect_type_radio_button(CutConnectorType::Dowel);
    if (type_changed)
        apply_selected_connectors([this, &connectors] (size_t idx) { connectors[idx].attribs.type = CutConnectorType(m_connector_type); });

    if (render_combo(_u8L("Style"), m_connector_styles, m_connector_style))
        apply_selected_connectors([this, &connectors](size_t idx) { connectors[idx].attribs.style = CutConnectorStyle(m_connector_style); });

    if (render_combo(_u8L("Shape"), m_connector_shapes, m_connector_shape_id))
        apply_selected_connectors([this, &connectors](size_t idx) { connectors[idx].attribs.shape = CutConnectorShape(m_connector_shape_id); });

    if (render_slider_double_input(_u8L("Depth ratio"), m_connector_depth_ratio, m_connector_depth_ratio_tolerance))
        apply_selected_connectors([this, &connectors](size_t idx) { 
            connectors[idx].height           = float(m_connector_depth_ratio);
            connectors[idx].height_tolerance = 0.01f * float(m_connector_depth_ratio_tolerance);
        });

    if (render_slider_double_input(_u8L("Size"), m_connector_size, m_connector_size_tolerance))
        apply_selected_connectors([this, &connectors](size_t idx) { 
            connectors[idx].radius           = float(m_connector_size * 0.5);
            connectors[idx].radius_tolerance = 0.01f * float(m_connector_size_tolerance);
        });

    ImGui::Separator();

    if (m_imgui->button(_L("Confirm connectors"))) {
        m_clp_normal         = m_c->object_clipper()->get_clipping_plane()->get_normal();
        unselect_all_connectors();
        set_connectors_editing(false);
    }
}

void GLGizmoCut3D::render_build_size()
{
    double              koef     = m_imperial_units ? ObjectManipulation::mm_to_in : 1.0;
    wxString            unit_str = " " + (m_imperial_units ? _L("in") : _L("mm"));
    const BoundingBoxf3 tbb      = transformed_bounding_box();
            
    Vec3d    tbb_sz = tbb.size();
    wxString size   =   "X: " + double_to_string(tbb_sz.x() * koef, 2) + unit_str +
                     ",  Y: " + double_to_string(tbb_sz.y() * koef, 2) + unit_str +
                     ",  Z: " + double_to_string(tbb_sz.z() * koef, 2) + unit_str;

    ImGui::AlignTextToFramePadding();
    m_imgui->text(_L("Build size"));
    ImGui::SameLine(m_label_width);
    m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, size);
}

void GLGizmoCut3D::reset_cut_plane()
{
    set_center(bounding_box().center());
    m_rotation_m = Transform3d::Identity();
    m_angle_arc.reset();
    update_clipper();
}

void GLGizmoCut3D::invalidate_cut_plane()
{
    m_rotation_m    = Transform3d::Identity();
    m_plane_center  = Vec3d::Zero();
    m_min_pos       = Vec3d::Zero();
    m_max_pos       = Vec3d::Zero();
    m_bb_center     = Vec3d::Zero();
    m_center_offset = Vec3d::Zero();
}

void GLGizmoCut3D::set_connectors_editing(bool connectors_editing)
{
    m_connectors_editing = connectors_editing;
    update_raycasters_for_picking();

    m_parent.request_extra_frame();
}

void GLGizmoCut3D::render_cut_plane_input_window(CutConnectors &connectors)
{
    // WIP : cut plane mode
    // render_combo(_u8L("Mode"), m_modes, m_mode);

    if (m_mode == size_t(CutMode::cutPlanar)) {
        ImGui::AlignTextToFramePadding();
        m_imgui->text(wxString(ImGui::InfoMarkerSmall));
        ImGui::SameLine();
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Hold SHIFT key and connect some two points of an object to cut by line"));

        ImGui::Separator();

        render_build_size();

        ImGui::AlignTextToFramePadding();
        m_imgui->text(_L("Cut position: "));
        ImGui::SameLine(m_label_width);
        render_move_center_input(Z);
        ImGui::SameLine();
        if (render_reset_button("cut_plane", _u8L("Reset cutting plane")))
            reset_cut_plane();

        if (wxGetApp().plater()->printer_technology() == ptFFF) {
            m_imgui->disabled_begin(!m_keep_upper || !m_keep_lower);
                if (m_imgui->button(_L("Add/Edit connectors")))
                    set_connectors_editing(true);
            m_imgui->disabled_end();
        }

        ImGui::Separator();

        auto render_part_action_line = [this, connectors](const wxString& label, const wxString& suffix, bool& keep_part, bool& place_on_cut_part, bool& rotate_part) {
            bool keep = true;
            ImGui::AlignTextToFramePadding();
            m_imgui->text(label);

            ImGui::SameLine(m_label_width);

            m_imgui->disabled_begin(!connectors.empty());
                m_imgui->checkbox(_L("Keep") + suffix, connectors.empty() ? keep_part : keep);
            m_imgui->disabled_end();

            ImGui::SameLine(2 * m_label_width);

            m_imgui->disabled_begin(!keep_part);
                if (m_imgui->checkbox(_L("Place on cut") + suffix, place_on_cut_part))
                    rotate_part = false;
                ImGui::SameLine();
                if (m_imgui->checkbox(_L("Flip") + suffix, rotate_part))
                    place_on_cut_part = false;
            m_imgui->disabled_end();
        };

        m_imgui->text(_L("After cut") + ": ");
        render_part_action_line( _L("Upper part"), "##upper", m_keep_upper, m_place_on_cut_upper, m_rotate_upper);
        render_part_action_line( _L("Lower part"), "##lower", m_keep_lower, m_place_on_cut_lower, m_rotate_lower);
    }

    ImGui::Separator();

    m_imgui->disabled_begin(!can_perform_cut());
        if(m_imgui->button(_L("Perform cut")))
            perform_cut(m_parent.get_selection());
    m_imgui->disabled_end();
}

void GLGizmoCut3D::init_input_window_data(CutConnectors &connectors)
{
    m_imperial_units = wxGetApp().app_config->get("use_inches") == "1";
    m_label_width    = m_imgui->get_style_scaling() * 100.0f;
    m_control_width  = m_imgui->get_style_scaling() * 150.0f;

    if (m_selected_count == 1)
        for (size_t idx = 0; idx < m_selected.size(); idx++)
            if (m_selected[idx]) {
                auto&connector                    = connectors[idx];
                m_connector_depth_ratio           = connector.height;
                m_connector_depth_ratio_tolerance = 100 * connector.height_tolerance;
                m_connector_size                  = 2. * connector.radius;
                m_connector_size_tolerance        = 100 * connector.radius_tolerance;
                m_connector_type                  = connector.attribs.type;
                m_connector_style                 = size_t(connector.attribs.style);
                m_connector_shape_id              = size_t(connector.attribs.shape);

                break;
            }
}

void GLGizmoCut3D::render_input_window_warning() const
{
    if (wxGetApp().plater()->printer_technology() == ptFFF && m_has_invalid_connector)
        m_imgui->text(wxString(ImGui::WarningMarkerSmall) + _L("Invalid connectors detected."));
    if (!m_keep_upper && !m_keep_lower)
        m_imgui->text(wxString(ImGui::WarningMarkerSmall) + _L("Invalid state. \nNo one part is selected for keep after cut"));
}

void GLGizmoCut3D::on_render_input_window(float x, float y, float bottom_limit)
{
    m_imgui->begin(get_name(), ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

    // adjust window position to avoid overlap the view toolbar
    adjust_window_position(x, y, bottom_limit);

    CutConnectors& connectors = m_c->selection_info()->model_object()->cut_connectors;

    init_input_window_data(connectors);

    if (m_connectors_editing) // connectors mode
        render_connectors_input_window(connectors); 
    else
        render_cut_plane_input_window(connectors);

    render_input_window_warning();

    m_imgui->end();

    render_debug_input_window();
}

// get volume transformation regarding to the "border". Border is related from the size of connectors
Transform3d GLGizmoCut3D::get_volume_transformation(const ModelVolume* volume) const
{
    bool is_prizm_dowel = m_connector_type == CutConnectorType::Dowel && m_connector_style == size_t(CutConnectorStyle::Prizm);
    const Transform3d connector_trafo = assemble_transform(
        is_prizm_dowel ? Vec3d(0.0, 0.0, -m_connector_depth_ratio) : Vec3d::Zero(),
        Transformation(m_rotation_m).get_rotation(),
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
    return assemble_transform(offset, Vec3d::Zero(), Vec3d::Ones() - border_scale, Vec3d::Ones()) * vol_matrix;
}

void GLGizmoCut3D::render_connectors()
{
    ::glEnable(GL_DEPTH_TEST);

    if (cut_line_processing() || m_connector_mode == CutConnectorMode::Auto || !m_c->selection_info())
        return;

    const ModelObject* mo = m_c->selection_info()->model_object();
    auto inst_id = m_c->selection_info()->get_active_instance();
    if (inst_id < 0)
        return;
    const CutConnectors& connectors = mo->cut_connectors;
    if (connectors.size() != m_selected.size()) {
        // #ysFIXME
        m_selected.clear();
        m_selected.resize(connectors.size(), false);
    }

    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    shader->start_using();

    ScopeGuard guard([shader]() { shader->stop_using(); });

    const Camera& camera = wxGetApp().plater()->get_camera();

    ColorRGBA render_color;

    const ModelInstance* mi = mo->instances[inst_id];
    const Vec3d& instance_offset = mi->get_offset();
    const double sla_shift       = double(m_c->selection_info()->get_sla_shift());

    const ClippingPlane* cp = m_c->object_clipper()->get_clipping_plane();
    const Vec3d& normal = cp && cp->is_active() ? cp->get_normal() : m_clp_normal;

    const Transform3d instance_trafo = assemble_transform(Vec3d(0.0, 0.0, sla_shift)) * mi->get_transformation().get_matrix();

    m_has_invalid_connector = false;

    for (size_t i = 0; i < connectors.size(); ++i) {
        const CutConnector& connector = connectors[i];

        float height = connector.height;
        // recalculate connector position to world position
        Vec3d pos = connector.pos + instance_offset;
        if (connector.attribs.type  == CutConnectorType::Dowel &&
            connector.attribs.style == CutConnectorStyle::Prizm) {
            pos -= height * normal;
            height *= 2;
        }
        pos[Z] += sla_shift;

        // First decide about the color of the point.
        if (!m_connectors_editing)
            render_color = CONNECTOR_DEF_COLOR;
        else if (size_t(m_hover_id - m_connectors_group_id) == i)
            render_color = connector.attribs.type == CutConnectorType::Dowel ? HOVERED_DOWEL_COLOR : HOVERED_PLAG_COLOR;
        else if (m_selected[i])
            render_color = connector.attribs.type == CutConnectorType::Dowel ? SELECTED_DOWEL_COLOR : SELECTED_PLAG_COLOR;
        else // neither hover nor picking
            render_color = connector.attribs.type == CutConnectorType::Dowel ? DOWEL_COLOR          : PLAG_COLOR;
        // ! #ysFIXME rework get_volume_transformation
        if (0) { // else { // neither hover nor picking
            int mesh_id = -1;
            for (const ModelVolume* mv : mo->volumes) {
                ++mesh_id;
                if (!mv->is_model_part())
                    continue;

                const Transform3d volume_trafo = get_volume_transformation(mv);

                if (m_c->raycaster()->raycasters()[mesh_id]->is_valid_intersection(pos, -normal, instance_trafo * volume_trafo)) {
                    render_color = m_connectors_editing ? ColorRGBA(1.0f, 1.0f, 1.0f, 0.5f) : /*ColorRGBA(0.5f, 0.5f, 0.5f, 1.f)*/ColorRGBA(1.0f, 0.3f, 0.3f, 0.5f);
                    break;
                }
                render_color = ColorRGBA(1.0f, 0.3f, 0.3f, 0.5f);
                m_has_invalid_connector = true;
            }
        }

        m_shapes[connector.attribs].model.set_color(render_color);

        const Transform3d scale_trafo = scale_transform(Vec3f(connector.radius, connector.radius, height).cast<double>());
        const Transform3d view_model_matrix = camera.get_view_matrix() * translation_transform(pos) * m_rotation_m * scale_trafo;
        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());

        m_shapes[connector.attribs].model.render();
    }
}

bool GLGizmoCut3D::can_perform_cut() const
{
    if (m_has_invalid_connector || (!m_keep_upper && !m_keep_lower) || m_connectors_editing)
        return false;

    const BoundingBoxf3 tbb = transformed_bounding_box(true);
    return tbb.contains(m_plane_center);
}

void GLGizmoCut3D::apply_connectors_in_model(ModelObject* mo, const bool has_connectors, bool &create_dowels_as_separate_object)
{
    if (has_connectors && m_connector_mode == CutConnectorMode::Manual) {
        m_selected.clear();

        for (CutConnector&connector : mo->cut_connectors) {
            connector.rotation_m = m_rotation_m;

            if (connector.attribs.type == CutConnectorType::Dowel) {
                if (connector.attribs.style == CutConnectorStyle::Prizm)
                    connector.height *= 2;
                create_dowels_as_separate_object = true;
            }
            else {
                // culculate shift of the connector center regarding to the position on the cut plane
                Vec3d shifted_center = m_plane_center + Vec3d::UnitZ();
                rotate_vec3d_around_plane_center(shifted_center);
                Vec3d norm = (shifted_center - m_plane_center).normalized();
                connector.pos += norm * 0.5 * double(connector.height);
            }
        }
        mo->apply_cut_connectors(_u8L("Connector"));
    }
}

void GLGizmoCut3D::perform_cut(const Selection& selection)
{
    const int instance_idx = selection.get_instance_idx();
    const int object_idx = selection.get_object_idx();

    wxCHECK_RET(instance_idx >= 0 && object_idx >= 0, "GLGizmoCut: Invalid object selection");

    Plater* plater = wxGetApp().plater();
    ModelObject* mo = plater->model().objects[object_idx];
    if (!mo)
        return;

    // m_cut_z is the distance from the bed. Subtract possible SLA elevation.
    const double sla_shift_z    = selection.get_first_volume()->get_sla_shift_z();
    const double object_cut_z   = m_plane_center.z() - sla_shift_z;

    const Vec3d instance_offset = mo->instances[instance_idx]->get_offset();
    Vec3d cut_center_offset     = m_plane_center - instance_offset;
    cut_center_offset[Z] -= sla_shift_z;

    if (0.0 < object_cut_z && can_perform_cut()) {
        bool create_dowels_as_separate_object = false;
        const bool has_connectors = !mo->cut_connectors.empty();
        {
            Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Cut by Plane"));
            // update connectors pos as offset of its center before cut performing
            apply_connectors_in_model(mo, has_connectors, create_dowels_as_separate_object);
        }

        plater->cut(object_idx, instance_idx, assemble_transform(cut_center_offset) * m_rotation_m,
                    only_if(has_connectors ? true : m_keep_upper, ModelObjectCutAttribute::KeepUpper) |
                    only_if(has_connectors ? true : m_keep_lower, ModelObjectCutAttribute::KeepLower) |
                    only_if(m_place_on_cut_upper, ModelObjectCutAttribute::PlaceOnCutUpper) | 
                    only_if(m_place_on_cut_lower, ModelObjectCutAttribute::PlaceOnCutLower) | 
                    only_if(m_rotate_upper, ModelObjectCutAttribute::FlipUpper) | 
                    only_if(m_rotate_lower, ModelObjectCutAttribute::FlipLower) | 
                    only_if(create_dowels_as_separate_object, ModelObjectCutAttribute::CreateDowels));
    }
    else {
        // the object is SLA-elevated and the plane is under it.
    }

    invalidate_cut_plane();
}



// Unprojects the mouse position on the mesh and saves hit point and normal of the facet into pos_and_normal
// Return false if no intersection was found, true otherwise.
bool GLGizmoCut3D::unproject_on_cut_plane(const Vec2d& mouse_position, std::pair<Vec3d, Vec3d>& pos_and_normal, Vec3d& pos_world)
{
    const float sla_shift = m_c->selection_info()->get_sla_shift();

    const ModelObject* mo = m_c->selection_info()->model_object();
    const ModelInstance* mi = mo->instances[m_c->selection_info()->get_active_instance()];
    const Transform3d    instance_trafo = sla_shift > 0.0 ? 
        assemble_transform(Vec3d(0.0, 0.0, sla_shift)) * mi->get_transformation().get_matrix() : mi->get_transformation().get_matrix();
    const Camera& camera = wxGetApp().plater()->get_camera();

    int mesh_id = -1;
    for (const ModelVolume* mv : mo->volumes) {
        ++mesh_id;
        if (!mv->is_model_part())
            continue;
        Vec3f normal;
        Vec3f hit;
        bool clipping_plane_was_hit = false;

//        const Transform3d volume_trafo = get_volume_transformation(mv);
        const Transform3d volume_trafo = mv->get_transformation().get_matrix();

        m_c->raycaster()->raycasters()[mesh_id]->unproject_on_mesh(mouse_position, instance_trafo * volume_trafo,
            camera, hit, normal, m_c->object_clipper()->get_clipping_plane(true),
            nullptr, &clipping_plane_was_hit);
        if (clipping_plane_was_hit) {
            // recalculate hit to object's local position
            Vec3d hit_d = hit.cast<double>();
            hit_d -= mi->get_offset();
            hit_d[Z] -= sla_shift;

            // Return both the point and the facet normal.
            pos_and_normal = std::make_pair(hit_d, normal.cast<double>());
            pos_world = hit.cast<double>();
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
    CutConnectorAttributes attribs = { m_connector_type, CutConnectorStyle(m_connector_style), CutConnectorShape(m_connector_shape_id) };
    if (m_shapes.find(attribs) == m_shapes.end()) {
        const indexed_triangle_set its = ModelObject::get_connector_mesh(attribs);
        m_shapes[attribs].model.init_from(its);
        m_shapes[attribs].mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(its)));;
    }

    const indexed_triangle_set its = ModelObject::get_connector_mesh(attribs);
    m_connector_mesh.clear();
    m_connector_mesh = TriangleMesh(its);
}

void GLGizmoCut3D::update_model_object()
{
    const ModelObjectPtrs& mos = wxGetApp().model().objects;
    ModelObject* mo = m_c->selection_info()->model_object();
    wxGetApp().obj_list()->update_info_items(std::find(mos.begin(), mos.end(), mo) - mos.begin());

    m_parent.post_event(SimpleEvent(EVT_GLCANVAS_SCHEDULE_BACKGROUND_PROCESS));

    update_raycasters_for_picking();
}

bool GLGizmoCut3D::cut_line_processing() const
{
    return m_line_beg != Vec3d::Zero();
}

void GLGizmoCut3D::discard_cut_line_processing()
{
    m_line_beg = m_line_end = Vec3d::Zero();
}

bool GLGizmoCut3D::process_cut_line(SLAGizmoEventType action, const Vec2d& mouse_position)
{
    const float sla_shift = m_c->selection_info()->get_sla_shift();
    const ModelObject* mo = m_c->selection_info()->model_object();
    const ModelInstance* mi = mo->instances[m_c->selection_info()->get_active_instance()];
    Transform3d inst_trafo = sla_shift > 0.f ?
        assemble_transform(Vec3d(0.0, 0.0, sla_shift)) * mi->get_transformation().get_matrix() :
        mi->get_transformation().get_matrix();

    const Camera& camera = wxGetApp().plater()->get_camera();

    Vec3d pt;
    Vec3d dir;
    MeshRaycaster::line_from_mouse_pos(mouse_position, Transform3d::Identity(), camera, pt, dir);
    dir.normalize();
    pt += dir; // Move the pt along dir so it is not clipped.

    if (action == SLAGizmoEventType::LeftDown && !cut_line_processing()) {
        m_line_beg = pt;
        m_line_end = pt;
        on_unregister_raycasters_for_picking();
        return true;
    }

    if (cut_line_processing()) {
        m_line_end = pt;
        if (action == SLAGizmoEventType::LeftDown || action == SLAGizmoEventType::LeftUp) {
            Vec3d line_dir = m_line_end - m_line_beg;
            if (line_dir.norm() < 3.0)
                return true;
            Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Cut by line"), UndoRedo::SnapshotType::GizmoAction);

            Vec3d cross_dir = line_dir.cross(dir).normalized();
            Eigen::Quaterniond q;
            Transform3d m = Transform3d::Identity();
            m.matrix().block(0, 0, 3, 3) = q.setFromTwoVectors(Vec3d::UnitZ(), cross_dir).toRotationMatrix();

            m_rotation_m = m;
            m_angle_arc.reset();

            set_center(m_plane_center + cross_dir * (cross_dir.dot(pt - m_plane_center)));

            discard_cut_line_processing();
        }
        else if (action == SLAGizmoEventType::Moving)
            this->set_dirty();
        return true;
    }
    return false;
}

bool GLGizmoCut3D::add_connector(CutConnectors& connectors, const Vec2d& mouse_position)
{
    std::pair<Vec3d, Vec3d> pos_and_normal;
    Vec3d pos_world;
    if (unproject_on_cut_plane(mouse_position.cast<double>(), pos_and_normal, pos_world)) {
        const Vec3d& hit = pos_and_normal.first;

        if (m_connectors_editing) {

            Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Add connector"), UndoRedo::SnapshotType::GizmoAction);

            connectors.emplace_back(hit, m_rotation_m,
                                    float(m_connector_size) * 0.5f, float(m_connector_depth_ratio),
                                    float(m_connector_size_tolerance) * 0.01f, float(m_connector_depth_ratio_tolerance) * 0.01f,
                                    CutConnectorAttributes( CutConnectorType(m_connector_type),
                                                            CutConnectorStyle(m_connector_style),
                                                            CutConnectorShape(m_connector_shape_id)));
            unselect_all_connectors();
            m_selected.push_back(true);
            assert(m_selected.size() == connectors.size());
            update_model_object();
            m_parent.set_as_dirty();
        }
        else {
            // Following would inform the clipper about the mouse click, so it can
            // toggle the respective contour as disabled.
            //m_c->object_clipper()->pass_mouse_click(pos_world);
        }

        return true;
    }
    return false;
}

bool GLGizmoCut3D::delete_selected_connectors(CutConnectors& connectors)
{
    if (connectors.empty())
        return false;

    Plater::TakeSnapshot snapshot(wxGetApp().plater(), _L("Delete connector"), UndoRedo::SnapshotType::GizmoAction);

    // remove  connectors
    for (int i = int(connectors.size()) - 1; i >= 0; i--)
        if (m_selected[i])
            connectors.erase(connectors.begin() + i);
    // remove selections
    m_selected.erase(std::remove_if(m_selected.begin(), m_selected.end(), [](const auto& selected) {
        return selected; }), m_selected.end());

    assert(m_selected.size() == connectors.size());
    update_model_object();
    m_parent.set_as_dirty();
    return true;
}

void GLGizmoCut3D::select_connector(int idx, bool select)
{
    m_selected[idx] = select;
    if (select)
        ++m_selected_count;
    else
        --m_selected_count;
}

bool GLGizmoCut3D::is_selection_changed(bool alt_down, bool control_down)
{
    if (m_hover_id >= m_connectors_group_id) {
        if (alt_down)
            select_connector(m_hover_id - m_connectors_group_id, false);
        else {
            if (!control_down)
                unselect_all_connectors();
            select_connector(m_hover_id - m_connectors_group_id, true);
        }
        return true;
    }
    return false;
}

bool GLGizmoCut3D::gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down)
{
    if (is_dragging() || m_connector_mode == CutConnectorMode::Auto || (!m_keep_upper || !m_keep_lower))
        return false;

    if ( m_hover_id < 0 && shift_down &&  ! m_connectors_editing &&
        (action == SLAGizmoEventType::LeftDown || action == SLAGizmoEventType::LeftUp || action == SLAGizmoEventType::Moving) )
        return process_cut_line(action, mouse_position);

    CutConnectors& connectors = m_c->selection_info()->model_object()->cut_connectors;

    if (action == SLAGizmoEventType::LeftDown && !shift_down) {
        // If there is no selection and no hovering, add new point
        if (m_hover_id == -1 && !control_down && !alt_down)
            if (!add_connector(connectors, mouse_position))
                unselect_all_connectors();
        return true;
    }
    if (!m_connectors_editing)
        return false;

    if (action == SLAGizmoEventType::LeftUp && !shift_down)
        return is_selection_changed(alt_down, control_down);
    
    if (action == SLAGizmoEventType::RightDown && !shift_down) {
        // If any point is in hover state, this should initiate its move - return control back to GLCanvas:
        if (m_hover_id < m_connectors_group_id)
            return false;
        unselect_all_connectors();
        select_connector(m_hover_id - m_connectors_group_id, true);
        return delete_selected_connectors(connectors);
    }
    
    if (action == SLAGizmoEventType::Delete)
        return delete_selected_connectors(connectors);

    if (action == SLAGizmoEventType::SelectAll) {
        select_all_connectors();
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
