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

GLGizmoCut3D::GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_rotation_gizmo(GLGizmoRotate3D(parent, "", -1))
    , m_rotation_matrix(  Eigen::AngleAxisd(0.0, Vec3d::UnitZ())
                        * Eigen::AngleAxisd(0.0, Vec3d::UnitY())
                        * Eigen::AngleAxisd(0.0, Vec3d::UnitX()))
{
    m_rotation_gizmo.use_only_grabbers();
    m_group_id = 3;

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
    if (tooltip.empty()) {
        double koef = wxGetApp().app_config->get("use_inches") == "1" ? ObjectManipulation::mm_to_in : 1.0;
        if (m_hover_id == m_group_id || m_grabbers[0].dragging)
            return "X: " + format(m_plane_center.x() * koef, 2) + "; " +//"\n" +
                   "Y: " + format(m_plane_center.y() * koef, 2) + "; " +//"\n" +
                   "Z: " + format(m_plane_center.z() * koef, 2);
    }

    return tooltip;
}

bool GLGizmoCut3D::on_mouse(const wxMouseEvent &mouse_event)
{
    if (m_rotation_gizmo.on_mouse(mouse_event)) {
        update_clipper();
        return true;
    }
    return use_grabbers(mouse_event);
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

    m_c->object_clipper()->set_range_and_pos(beg, end, dist);
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
    m_c->object_clipper()->render_cut();

    if (m_hide_cut_plane)
        return;

    const BoundingBoxf3 box = bounding_box();

    const float min_x = box.min.x() - Margin - m_plane_center.x();
    const float max_x = box.max.x() + Margin - m_plane_center.x();
    const float min_y = box.min.y() - Margin - m_plane_center.y();
    const float max_y = box.max.y() + Margin - m_plane_center.y();

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_CULL_FACE));
    glsafe(::glEnable(GL_BLEND));
    glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;
    shader->start_using();

    const Vec3d& angles = m_rotation_gizmo.get_rotation();

    glsafe(::glPushMatrix());
    glsafe(::glTranslated(m_plane_center.x(), m_plane_center.y(), m_plane_center.z()));
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

void GLGizmoCut3D::render_cut_center_graber()
{
    const Vec3d& angles = m_rotation_gizmo.get_rotation();
    const BoundingBoxf3 box = bounding_box();

    Vec3d grabber_center = m_plane_center;
    grabber_center[Z] += 10; // Margin

    rotate_vec3d_around_center(grabber_center, angles, m_plane_center);

    m_grabbers[0].center = grabber_center;
    m_grabbers[0].angles = angles;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
    glsafe(::glLineWidth(m_hover_id == m_group_id ? 2.0f : 1.5f));

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();
        //       if (!m_grabber_connection.is_initialized() || z_changed)
        {
            m_grabber_connection.reset();

            GLModel::Geometry init_data;
            init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
            init_data.color = ColorRGBA::YELLOW();
            init_data.reserve_vertices(2);
            init_data.reserve_indices(2);

            // vertices
            init_data.add_vertex((Vec3f)m_plane_center.cast<float>());
            init_data.add_vertex((Vec3f)m_grabbers[0].center.cast<float>());

            // indices
            init_data.add_ushort_line(0, 1);

            m_grabber_connection.init_from(std::move(init_data));
        }
        m_grabber_connection.render();

        shader->stop_using();
    }

    shader = wxGetApp().get_shader("gouraud_light");
    if (shader != nullptr) {
        shader->start_using();
        shader->set_uniform("emission_factor", 0.1f);

        m_grabbers[0].color = GRABBER_COLOR;
        m_grabbers[0].render(m_hover_id == m_group_id, float((box.size().x() + box.size().y() + box.size().z()) / 3.0));

        shader->stop_using();
    }
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
    assert(m_hover_id == m_group_id);

    const Vec3d & starting_box_center = m_plane_center;
    const Vec3d & starting_drag_position = m_grabbers[0].center;
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

    set_center(starting_box_center + starting_vec * projection);
}

void GLGizmoCut3D::set_center_pos(const Vec3d& center_pos)
{
    m_plane_center = center_pos;

    // !!! ysFIXME add smart clamp calculation
    // Clamp the center position of the cut plane to the object's bounding box
    //m_plane_center = Vec3d(std::clamp(center_pos.x(), m_min_pos.x(), m_max_pos.x()),
    //                       std::clamp(center_pos.y(), m_min_pos.y(), m_max_pos.y()),
    //                       std::clamp(center_pos.z(), m_min_pos.z(), m_max_pos.z())));

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

bool GLGizmoCut3D::update_bb()
{
    const BoundingBoxf3 box = bounding_box();
    if (m_max_pos != box.max && m_min_pos != box.min) {
        m_max_pos = box.max;
        m_min_pos = box.min;
        m_bb_center = box.center();
        set_center_pos(m_bb_center + m_center_offset);
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

    render_cut_plane();
    render_cut_center_graber();
    if (m_mode == CutMode::cutPlanar) {
        if (m_hover_id < m_group_id)
            m_rotation_gizmo.render();
    }

    if (!suppress_update_clipper_on_render)
        update_clipper_on_render();
}

void GLGizmoCut3D::on_render_for_picking()
{
    m_rotation_gizmo.render_for_picking();
    render_grabbers_for_picking(m_parent.get_selection().get_bounding_box());
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

    if (m_mode <= CutMode::cutByLine) {
        ImGui::Separator();

        if (m_mode == CutMode::cutPlanar) {
            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Move center"));
            revert_move = render_revert_button("move");
            for (Axis axis : {X, Y, Z})
                render_move_center_input(axis);
            m_imgui->text(m_imperial_units ? _L("in") : _L("mm")); 

            ImGui::AlignTextToFramePadding();
            m_imgui->text(_L("Rotation"));
            revert_rotation = render_revert_button("rotation");
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

    ImGui::Separator();

    m_imgui->checkbox("hide_cut_plane", m_hide_cut_plane);

    ////////
    static bool hide_clipped = true;
    static bool fill_cut = true;
    static float contour_width = 0.2f;
    m_imgui->checkbox("hide_clipped", hide_clipped);
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

bool GLGizmoCut3D::can_perform_cut() const
{
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

    Vec3d instance_offset = wxGetApp().plater()->model().objects[object_idx]->instances[instance_idx]->get_offset();

    Vec3d cut_center_offset = m_plane_center - instance_offset;
    cut_center_offset[Z] -= first_glvolume->get_sla_shift_z();

    if (0.0 < object_cut_z && can_perform_cut())
        wxGetApp().plater()->cut(object_idx, instance_idx, cut_center_offset, m_rotation_gizmo.get_rotation(),
            only_if(m_keep_upper, ModelObjectCutAttribute::KeepUpper) |
            only_if(m_keep_lower, ModelObjectCutAttribute::KeepLower) |
            only_if(m_rotate_lower, ModelObjectCutAttribute::FlipLower));
    else {
        // the object is SLA-elevated and the plane is under it.
    }
}

bool GLGizmoCut3D::gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down)
{
    if (is_dragging())
        return false;
    const ModelObject   *mo                           = m_c->selection_info()->model_object();
    const ModelInstance *mi                           = mo->instances[m_c->selection_info()->get_active_instance()];
    const Transform3d    instance_trafo               = mi->get_transformation().get_matrix();
    const Transform3d    instance_trafo_not_translate = mi->get_transformation().get_matrix(true);
    const Camera& camera = wxGetApp().plater()->get_camera();

    int mesh_id = -1;
    for (const ModelVolume *mv : mo->volumes) {
        ++mesh_id;
        if (! mv->is_model_part())
            continue;
        Vec3f hit;
        Vec3f normal;
        bool clipping_plane_was_hit = false;
        m_c->raycaster()->raycasters()[mesh_id]->unproject_on_mesh(mouse_position, instance_trafo * mv->get_matrix(),
                        camera, hit, normal, m_c->object_clipper()->get_clipping_plane(),
                        nullptr, &clipping_plane_was_hit);
        if (clipping_plane_was_hit) {
            // The clipping plane was clicked, hit containts coordinates of the hit in world coords.
            std::cout << hit.x() << "\t" << hit.y() << "\t" << hit.z() << std::endl;
            return true;
        }
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
