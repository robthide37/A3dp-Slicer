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

const double GLGizmoCenterMove::Margin = 20.0;

GLGizmoCenterMove::GLGizmoCenterMove(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoMove3D(parent, "", -1)
{
}

void GLGizmoCenterMove::set_center_pos(const Vec3d& centre_pos)
{
    // Clamp the center position of the cut plane to the object's bounding box
    set_center(Vec3d(std::clamp(centre_pos.x(), m_min_pos.x(), m_max_pos.x()),
        std::clamp(centre_pos.y(), m_min_pos.y(), m_max_pos.y()),
        std::clamp(centre_pos.z(), m_min_pos.z(), m_max_pos.z())));
}

std::string GLGizmoCenterMove::get_tooltip() const
{
    double koef = wxGetApp().app_config->get("use_inches") == "1" ? ObjectManipulation::mm_to_in : 1.0;

    const Vec3d& center_pos = get_center();

    if (m_hover_id == 0 || m_grabbers[0].dragging)
        return "X: " + format(center_pos.x() * koef, 2);
    else if (m_hover_id == 1 || m_grabbers[1].dragging)
        return "Y: " + format(center_pos.y() * koef, 2);
    else if (m_hover_id == 2 || m_grabbers[2].dragging)
        return "Z: " + format(center_pos.z() * koef, 2);
    else
        return "";
}

void GLGizmoCenterMove::on_set_state()
{
    // Reset internal variables on gizmo activation, if bounding box was changed
    if (get_state() == On) {
        const BoundingBoxf3 box = bounding_box();
        if (m_max_pos != box.max && m_min_pos != box.min) {
            m_max_pos = box.max;
            m_min_pos = box.min;
            set_center_pos(box.center());
        }
    }
}

void GLGizmoCenterMove::on_update(const UpdateData& data)
{
    GLGizmoMove3D::on_update(data);
    set_center_pos(get_center());
}

BoundingBoxf3 GLGizmoCenterMove::bounding_box() const
{
    BoundingBoxf3 ret;
    const Selection& selection = m_parent.get_selection();
    const Selection::IndicesList& idxs = selection.get_volume_idxs();
    for (unsigned int i : idxs) {
        const GLVolume* volume = selection.get_volume(i);
        if (!volume->is_modifier)
            ret.merge(volume->transformed_convex_hull_bounding_box());
    }
    return ret;
}



GLGizmoCut3D::GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_rotation_gizmo(GLGizmoRotate3D(parent, "", -1))
    , m_move_gizmo(GLGizmoCenterMove(parent, "", -1))
{
    m_move_gizmo.set_group_id(3);

    m_modes = { _u8L("Planar"),  _u8L("By Line"),_u8L("Grid")
//              , _u8L("Radial"), _u8L("Modular")
    };

    m_connector_modes = { _u8L("Auto"), _u8L("Manual") };
    m_connector_types = { _u8L("Plug"), _u8L("Dowel") };

    m_connector_styles = { _u8L("Prizm"), _u8L("Frustrum")
//              , _u8L("Claw")
    };

    m_connector_shapes = { _u8L("Triangle"), _u8L("Square"), _u8L("Circle"), _u8L("Hexagon")
//              , _u8L("D-shape")
    };

    m_axis_names = { "X", "Y", "Z" };
}

std::string GLGizmoCut3D::get_tooltip() const
{
    std::string tooltip = m_rotation_gizmo.get_tooltip();
    if (tooltip.empty())
        tooltip = m_move_gizmo.get_tooltip();
    return tooltip;
}

void GLGizmoCut3D::shift_cut_z(double delta)
{
    Vec3d new_cut_center = m_move_gizmo.get_center();
    new_cut_center[Z] += delta;
    set_center(new_cut_center);
}

void GLGizmoCut3D::set_center(const Vec3d& center)
{
    m_move_gizmo.set_center_pos(center);
    m_rotation_gizmo.set_center(m_move_gizmo.get_center());
}

void GLGizmoCut3D::render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx)
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

    selection_idx = selection_out;
}

void GLGizmoCut3D::render_double_input(const std::string& label, double& value_in)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(label);
    ImGui::SameLine(m_label_width);
    ImGui::PushItemWidth(m_control_width);

    double value = value_in;
    if (m_imperial_units)
        value *= ObjectManipulation::mm_to_in;
    ImGui::InputDouble(("##" + label).c_str(), &value, 0.0f, 0.0f, "%.2f", ImGuiInputTextFlags_CharsDecimal);

    ImGui::SameLine();
    m_imgui->text(m_imperial_units ? _L("in") : _L("mm"));

    value_in = value * (m_imperial_units ? ObjectManipulation::in_to_mm : 1.0);
}

void GLGizmoCut3D::render_move_center_input(int axis)
{
    ImGui::AlignTextToFramePadding();
    m_imgui->text(m_axis_names[axis]+":");
    ImGui::SameLine();
    ImGui::PushItemWidth(0.3*m_control_width);

    Vec3d move = m_move_gizmo.get_center();
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
    double value = rotation[axis] * (180. / M_PI);
    if (value > 360)
        value -= 360;

    ImGui::PushItemWidth(0.3*m_control_width);
    ImGui::InputDouble(("##rotate_" + m_axis_names[axis]).c_str(), &value, 0.0f, 0.0f, "%.2f", ImGuiInputTextFlags_CharsDecimal);
    ImGui::SameLine();

    rotation[axis] = (M_PI / 180.) * value;
    m_rotation_gizmo.set_rotation(rotation);
}

void GLGizmoCut3D::render_connect_type_radio_button(ConnectorType type)
{
    ImGui::SameLine(type == ConnectorType::Plug ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_types[int(type)], m_connector_type == type))
        m_connector_type = type;
}

void GLGizmoCut3D::render_connect_mode_radio_button(ConnectorMode mode)
{
    ImGui::SameLine(mode == ConnectorMode::Auto ? m_label_width : 2*m_label_width);
    ImGui::PushItemWidth(m_control_width);
    if (m_imgui->radio_button(m_connector_modes[int(mode)], m_connector_mode == mode))
        m_connector_mode = mode;
}

void GLGizmoCut3D::render_cut_plane()
{
    const BoundingBoxf3 box = m_move_gizmo.bounding_box();
    Vec3d plane_center = m_move_gizmo.get_center();// == Vec3d::Zero() ? box.center() : m_move_gizmo.get_center();

    const float min_x = box.min.x() - GLGizmoCenterMove::Margin - plane_center.x();
    const float max_x = box.max.x() + GLGizmoCenterMove::Margin - plane_center.x();
    const float min_y = box.min.y() - GLGizmoCenterMove::Margin - plane_center.y();
    const float max_y = box.max.y() + GLGizmoCenterMove::Margin - plane_center.y();
    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_CULL_FACE));
    glsafe(::glEnable(GL_BLEND));
    glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;
    shader->start_using();

    Vec3d angles = m_rotation_gizmo.get_rotation();

    glsafe(::glPushMatrix());
    glsafe(::glTranslated(plane_center.x(), plane_center.y(), plane_center.z()));
    glsafe(::glRotated(Geometry::rad2deg(angles.z()), 0.0, 0.0, 1.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.y()), 0.0, 1.0, 0.0));
    glsafe(::glRotated(Geometry::rad2deg(angles.x()), 1.0, 0.0, 0.0));

    if (!m_plane.is_initialized()) {
        m_plane.reset();

        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
        init_data.color = { 0.8f, 0.8f, 0.8f, 0.5f };
        init_data.reserve_vertices(4);
        init_data.reserve_indices(6);

        // vertices
        init_data.add_vertex(Vec3f(min_x, min_y, 0.0));
        init_data.add_vertex(Vec3f(max_x, min_y, 0.0));
        init_data.add_vertex(Vec3f(max_x, max_y, 0.0));
        init_data.add_vertex(Vec3f(min_x, max_y, 0.0));

        // indices
        init_data.add_ushort_triangle(0, 1, 2);
        init_data.add_ushort_triangle(2, 3, 0);

        m_plane.init_from(std::move(init_data));
    }

    m_plane.render();
    glsafe(::glPopMatrix());
#else
    // Draw the cutting plane
    ::glBegin(GL_QUADS);
    ::glColor4fv(PLANE_COLOR.data());
    ::glVertex3f(min_x, min_y, plane_center.z());
    ::glVertex3f(max_x, min_y, plane_center.z());
    ::glVertex3f(max_x, max_y, plane_center.z());
    ::glVertex3f(min_x, max_y, plane_center.z());
    glsafe(::glEnd());
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}

bool GLGizmoCut3D::on_init()
{
    if(!m_rotation_gizmo.init())
        return false;
    if(!m_move_gizmo.init())
        return false;

    m_shortcut_key = WXK_CONTROL_C;
    return true;
}

std::string GLGizmoCut3D::on_get_name() const
{
    return _u8L("Cut");
}

void GLGizmoCut3D::on_set_state() 
{
    m_move_gizmo.set_state(m_state);
    m_rotation_gizmo.set_center(m_move_gizmo.get_center());
    m_rotation_gizmo.set_state(m_state);
}

void GLGizmoCut3D::on_set_hover_id() 
{
    int move_group_id       = m_move_gizmo.get_group_id();
    m_rotation_gizmo.   set_hover_id((m_hover_id < move_group_id)  ? m_hover_id                 : -1);
    m_move_gizmo.       set_hover_id((m_hover_id >= move_group_id) ? m_hover_id - move_group_id : -1);
}

void GLGizmoCut3D::on_enable_grabber(unsigned int id) 
{
    m_rotation_gizmo.enable_grabber(id);
    m_move_gizmo.enable_grabber(id- m_move_gizmo.get_group_id());
}

void GLGizmoCut3D::on_disable_grabber(unsigned int id) 
{
    m_rotation_gizmo.disable_grabber(id);
    m_move_gizmo.disable_grabber(id- m_move_gizmo.get_group_id());
}

bool GLGizmoCut3D::on_is_activable() const
{
    return m_move_gizmo.is_activable();
}

void GLGizmoCut3D::on_start_dragging()
{
    m_rotation_gizmo.start_dragging();
    m_move_gizmo.start_dragging();
}

void GLGizmoCut3D::on_stop_dragging()
{
    m_rotation_gizmo.stop_dragging();
    m_rotation_gizmo.set_center(m_move_gizmo.get_center());
    m_move_gizmo.stop_dragging();
}

void GLGizmoCut3D::on_update(const UpdateData& data)
{
    m_move_gizmo.update(data);
    m_rotation_gizmo.update(data);
}

void GLGizmoCut3D::on_render()
{
    render_cut_plane();
    if (m_mode == CutMode::cutPlanar) {
        int move_group_id = m_move_gizmo.get_group_id();
        if (m_hover_id < move_group_id)
            m_rotation_gizmo.render();
        if (m_hover_id == -1 || m_hover_id >= move_group_id)
            m_move_gizmo.render();
    }
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

    if (m_mode <= CutMode::cutByLine) {
        ImGui::Separator();

        if (m_mode == CutMode::cutPlanar) {
            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Move center"));
            ImGui::SameLine(m_label_width);
            for (Axis axis : {X, Y, Z})
                render_move_center_input(axis);
            m_imgui->text(m_imperial_units ? _L("in") : _L("mm"));

            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Rotation"));
            ImGui::SameLine(m_label_width);
            for (Axis axis : {X, Y, Z})
                render_rotation_input(axis);
            m_imgui->text(_L("°"));
        }
        else {
            ImGui::AlignTextToFramePadding();
            ImGui::AlignTextToFramePadding();
            ImGui::AlignTextToFramePadding();
            ImGui::PushTextWrapPos(ImGui::GetCursorPos().x + 3*m_control_width);
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Connect some two points of object to cteate a cut plane"));
            ImGui::PopTextWrapPos();
            ImGui::AlignTextToFramePadding();
            ImGui::AlignTextToFramePadding();
        }

        ImGui::AlignTextToFramePadding();
        m_imgui->text(_L("After cut"));
        ImGui::SameLine(m_label_width);
        m_imgui->checkbox(_L("Keep upper part"), m_keep_upper);
        m_imgui->text("");
        ImGui::SameLine(m_label_width);
        m_imgui->checkbox(_L("Keep lower part"), m_keep_lower);
        m_imgui->text("");
        ImGui::SameLine(m_label_width);
        m_imgui->disabled_begin(!m_keep_lower);
        m_imgui->checkbox(_L("Rotate lower part upwards"), m_rotate_lower);
        m_imgui->disabled_end();
    }

    m_imgui->disabled_begin(!m_keep_lower || !m_keep_upper);
    // Connectors section
    ImGui::Separator();

    ImGui::AlignTextToFramePadding();
    m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Connectors"));

    m_imgui->text(_L("Mode"));
    render_connect_mode_radio_button(ConnectorMode::Auto);
    render_connect_mode_radio_button(ConnectorMode::Manual);

    m_imgui->text(_L("Type"));
    render_connect_type_radio_button(ConnectorType::Plug);
    render_connect_type_radio_button(ConnectorType::Dowel);

    render_combo(_u8L("Style"), m_connector_styles, m_connector_style);
    render_combo(_u8L("Shape"), m_connector_shapes, m_connector_shape);

    render_double_input(_u8L("Depth ratio"), m_connector_depth_ratio);
    render_double_input(_u8L("Size"), m_connector_size);

    m_imgui->disabled_end();

    ImGui::Separator();

    m_imgui->disabled_begin((!m_keep_upper && !m_keep_lower) || !can_perform_cut());
    const bool cut_clicked = m_imgui->button(_L("Perform cut"));
    m_imgui->disabled_end();

    m_imgui->end();

    if (cut_clicked && (m_keep_upper || m_keep_lower))
        perform_cut(m_parent.get_selection());
}

bool GLGizmoCut3D::can_perform_cut() const
{    
    return true;
}

void GLGizmoCut3D::perform_cut(const Selection& selection)
{
    const int instance_idx = selection.get_instance_idx();
    const int object_idx = selection.get_object_idx();

    wxCHECK_RET(instance_idx >= 0 && object_idx >= 0, "GLGizmoCut: Invalid object selection");

    // m_cut_z is the distance from the bed. Subtract possible SLA elevation.
    const GLVolume* first_glvolume = selection.get_volume(*selection.get_volume_idxs().begin());
    const double object_cut_z = m_move_gizmo.get_center().z() - first_glvolume->get_sla_shift_z();

    if (0.0 < object_cut_z && can_perform_cut())
        wxGetApp().plater()->cut(object_idx, instance_idx, object_cut_z,
            only_if(m_keep_upper, ModelObjectCutAttribute::KeepUpper) |
            only_if(m_keep_lower, ModelObjectCutAttribute::KeepLower) |
            only_if(m_rotate_lower, ModelObjectCutAttribute::FlipLower));
    else {
        // the object is SLA-elevated and the plane is under it.
    }
}

} // namespace GUI
} // namespace Slic3r
