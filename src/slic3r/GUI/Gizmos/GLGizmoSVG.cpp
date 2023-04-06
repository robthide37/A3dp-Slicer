#include "GLGizmoSVG.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/MainFrame.hpp" // to update title when add text
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/MsgDialog.hpp"
#include "slic3r/GUI/format.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "slic3r/GUI/Jobs/EmbossJob.hpp"
#include "slic3r/Utils/UndoRedo.hpp"

#include "libslic3r/Point.hpp"      
#include "libslic3r/SVG.hpp"      // debug store
#include "libslic3r/Geometry.hpp" // covex hull 2d
#include "libslic3r/Timer.hpp" // covex hull 2d

#include "libslic3r/NSVGUtils.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/ClipperUtils.hpp" // union_ex

#include "imgui/imgui_stdlib.h" // using std::string for inputs
#include "nanosvg/nanosvg.h"    // load SVG file

#include <wx/display.h> // detection of change DPI
#include <boost/log/trivial.hpp>

#include <GL/glew.h>
#include <chrono> // measure enumeration of fonts

using namespace Slic3r;
using namespace Slic3r::Emboss;
using namespace Slic3r::GUI;
using namespace Slic3r::GUI::Emboss;

namespace priv {

} // namespace priv

GLGizmoSVG::GLGizmoSVG(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -3)
    , m_volume(nullptr)
    , m_rotate_gizmo(parent, GLGizmoRotate::Axis::Z) // grab id = 2 (Z axis)
    , m_job_cancel(nullptr)
{
    m_rotate_gizmo.set_group_id(0);
    m_rotate_gizmo.set_force_local_coordinate(true);
}

// Private functions to create emboss volume
namespace{

// Variable keep limits for variables
const struct Limits
{
    MinMax<float> emboss{0.01f, 1e4f}; // in mm
    // distance text object from surface
    MinMax<float> angle{-180.f, 180.f}; // in degrees
} limits;

/// <summary>
/// Open file dialog with svg files
/// </summary>
/// <returns>File path to svg</returns>
std::string choose_svg_file();

/// <summary>
/// Let user to choose file with (S)calable (V)ector (G)raphics - SVG.
/// Than let select contour
/// </summary>
/// <returns>EmbossShape to create</returns>
EmbossShape select_shape();

/// <summary>
/// Create new embos data
/// </summary>
/// <param name="cancel">Cancel for previous job</param>
/// <returns>Base data for emboss SVG</returns>
DataBasePtr create_emboss_data_base(std::shared_ptr<std::atomic<bool>> &cancel);

/// <summary>
/// Create symbol '?' as default shape
/// without source file
/// with size 2cm
/// </summary>
/// <returns>Default shape to emboss</returns>
ExPolygons default_shape();

/// <summary>
/// Separate file name from file path.
/// String after last delimiter and before last point 
/// </summary>
/// <param name="file_path">path return by file dialog</param>
/// <returns>File name without directory path</returns>
std::string get_file_name(const std::string &file_path);

/// <summary>
/// Create volume name from shape information
/// </summary>
/// <param name="shape">File path</param>
/// <returns>Name for volume</returns>
std::string volume_name(const EmbossShape& shape);

/// <summary>
/// Create input for volume creation
/// </summary>
/// <param name="canvas">parent of gizmo</param>
/// <param name="raycaster">Keep scene</param>
/// <param name="volume_type">Type of volume to be created</param>
/// <returns>Params</returns>
CreateVolumeParams create_input(GLCanvas3D &canvas, RaycastManager &raycaster, ModelVolumeType volume_type);
} // namespace 


bool GLGizmoSVG::create_volume(ModelVolumeType volume_type, const Vec2d &mouse_pos)
{
    CreateVolumeParams input = create_input(m_parent, m_raycast_manager, volume_type);
    DataBasePtr base = create_emboss_data_base(m_job_cancel);
    return start_create_volume(input, std::move(base), mouse_pos);
}

bool GLGizmoSVG::create_volume(ModelVolumeType volume_type) 
{
    CreateVolumeParams input = create_input(m_parent, m_raycast_manager, volume_type);
    DataBasePtr base = create_emboss_data_base(m_job_cancel);
    return start_create_volume_without_position(input, std::move(base));
}

bool GLGizmoSVG::is_svg(const ModelVolume &volume) {
    return volume.emboss_shape.has_value();
}

bool GLGizmoSVG::is_svg_object(const ModelVolume &volume) {
    if (!volume.emboss_shape.has_value()) return false;
    if (volume.type() != ModelVolumeType::MODEL_PART) return false;
    for (const ModelVolume *v : volume.get_object()->volumes) {
        if (v->id() == volume.id()) continue;
        if (v->type() == ModelVolumeType::MODEL_PART) return false;
    }
    return true;
}

bool GLGizmoSVG::on_mouse_for_rotation(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Moving()) return false;

    bool used = use_grabbers(mouse_event);
    if (!m_dragging) return used;

    if (mouse_event.Dragging()) {
        auto &angle_opt = m_volume->text_configuration->style.prop.angle;
        if (!m_rotate_start_angle.has_value())
            m_rotate_start_angle = angle_opt.has_value() ? *angle_opt : 0.f;        
        double angle = m_rotate_gizmo.get_angle();
        angle -= PI / 2; // Grabber is upward

        // temporary rotation
        const TransformationType transformation_type = m_parent.get_selection().is_single_text() ?
          TransformationType::Local_Relative_Joint : TransformationType::World_Relative_Joint;
        m_parent.get_selection().rotate(Vec3d(0., 0., angle), transformation_type);

        angle += *m_rotate_start_angle;
        // move to range <-M_PI, M_PI>
        Geometry::to_range_pi_pi(angle);
        // propagate angle into property
        angle_opt = static_cast<float>(angle);

        // do not store zero
        if (is_approx(*angle_opt, 0.f))
            angle_opt.reset();
    }
    return used;
}

bool GLGizmoSVG::on_mouse_for_translate(const wxMouseEvent &mouse_event)
{
    // exist selected volume?
    if (m_volume == nullptr)
        return false;

    std::optional<double> up_limit;
    if (m_keep_up)
        up_limit = Slic3r::GUI::up_limit;
    const Camera &camera = wxGetApp().plater()->get_camera();

    bool was_dragging = m_surface_drag.has_value();
    bool res = on_mouse_surface_drag(mouse_event, camera, m_surface_drag, m_parent, m_raycast_manager, up_limit);
    bool is_dragging  = m_surface_drag.has_value();

    // End with surface dragging?
    if (was_dragging && !is_dragging) {
        // Update surface by new position
        if (m_volume->text_configuration->style.prop.use_surface)
            process();

        // Show correct value of height & depth inside of inputs
        calculate_scale();
    }

    // Start with dragging
    else if (!was_dragging && is_dragging) {
        // Cancel job to prevent interuption of dragging (duplicit result)
        if (m_job_cancel != nullptr)
            m_job_cancel->store(true);
    }

    // during drag
    else if (was_dragging && is_dragging) {
        // update scale of selected volume --> should be approx the same
        calculate_scale();

        // Recalculate angle for GUI
        if (!m_keep_up) {
            const GLVolume *gl_volume = get_selected_gl_volume(m_parent.get_selection());
            assert(gl_volume != nullptr);
            m_angle = calc_up(gl_volume->world_matrix(), Slic3r::GUI::up_limit);
        }
    }
    return res;
}

bool GLGizmoSVG::on_mouse(const wxMouseEvent &mouse_event)
{
    // not selected volume
    if (m_volume == nullptr ||
        get_model_volume(m_volume_id, m_parent.get_selection().get_model()->objects) == nullptr ||
        !m_volume->text_configuration.has_value()) return false;

    if (on_mouse_for_rotation(mouse_event)) return true;
    if (on_mouse_for_translate(mouse_event)) return true;

    return false;
}

bool GLGizmoSVG::wants_enter_leave_snapshots() const { return true; }
std::string GLGizmoSVG::get_gizmo_entering_text() const { return _u8L("Enter SVG gizmo"); }
std::string GLGizmoSVG::get_gizmo_leaving_text() const { return _u8L("Leave SVG gizmo"); }
std::string GLGizmoSVG::get_action_snapshot_name() const { return _u8L("SVG actions"); }

bool GLGizmoSVG::on_init()
{
    m_rotate_gizmo.init();
    ColorRGBA gray_color(.6f, .6f, .6f, .3f);
    m_rotate_gizmo.set_highlight_color(gray_color);

    // No shortCut
    // m_shortcut_key = WXK_CONTROL_T;

    // Set rotation gizmo upwardrotate
    m_rotate_gizmo.set_angle(PI / 2);
    return true;
}

std::string GLGizmoSVG::on_get_name() const { return _u8L("SVG"); }

void GLGizmoSVG::on_render() {
    // no volume selected
    if (m_volume == nullptr ||
        get_model_volume(m_volume_id, m_parent.get_selection().get_model()->objects) == nullptr)
        return;

    Selection &selection = m_parent.get_selection();
    if (selection.is_empty()) return;

    // prevent get local coordinate system on multi volumes
    if (!selection.is_single_volume_or_modifier() && 
        !selection.is_single_volume_instance()) return;
    bool is_surface_dragging = m_surface_drag.has_value();
    bool is_parent_dragging = m_parent.is_mouse_dragging();
    // Do NOT render rotation grabbers when dragging object
    bool is_rotate_by_grabbers = m_dragging;
    if (is_rotate_by_grabbers || 
        (!is_surface_dragging && !is_parent_dragging)) {
        glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
        m_rotate_gizmo.render();
    }
}

void GLGizmoSVG::on_register_raycasters_for_picking(){
    m_rotate_gizmo.register_raycasters_for_picking();
}
void GLGizmoSVG::on_unregister_raycasters_for_picking(){
    m_rotate_gizmo.unregister_raycasters_for_picking();
}

void GLGizmoSVG::on_render_input_window(float x, float y, float bottom_limit)
{
    set_volume_by_selection();

    // Configuration creation
    double screen_scale = wxDisplay(wxGetApp().plater()).GetScaleFactor();
    float  main_toolbar_height = m_parent.get_main_toolbar_height();
    if (!m_gui_cfg.has_value() || // Exist configuration - first run
        m_gui_cfg->screen_scale != screen_scale || // change of DPI
        m_gui_cfg->main_toolbar_height != main_toolbar_height // change size of view port
        ) {
        // Create cache for gui offsets
        GuiCfg cfg       = create_gui_configuration();
        cfg.screen_scale = screen_scale;
        cfg.main_toolbar_height = main_toolbar_height;
        m_gui_cfg.emplace(std::move(cfg));
        // set position near toolbar
        m_set_window_offset = ImVec2(-1.f, -1.f);
    }

    const ImVec2 &min_window_size = m_gui_cfg->minimal_window_size;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, min_window_size);

    // Draw origin position of text during dragging
    if (m_surface_drag.has_value()) {
        ImVec2 mouse_pos = ImGui::GetMousePos();
        ImVec2 center(
            mouse_pos.x + m_surface_drag->mouse_offset.x(),
            mouse_pos.y + m_surface_drag->mouse_offset.y());
        ImU32 color = ImGui::GetColorU32(
            m_surface_drag->exist_hit ? 
                ImVec4(1.f, 1.f, 1.f, .75f) : // transparent white
                ImVec4(1.f, .3f, .3f, .75f)
        ); // Warning color
        const float radius = 16.f;
        ImGuiWrapper::draw_cross_hair(center, radius, color);
    }

    // check if is set window offset
    if (m_set_window_offset.has_value()) {
        if (m_set_window_offset->y < 0)
            // position near toolbar
            m_set_window_offset = ImVec2(x, std::min(y, bottom_limit - min_window_size.y));
        
        ImGui::SetNextWindowPos(*m_set_window_offset, ImGuiCond_Always);
        m_set_window_offset.reset();
    } 

    bool is_opened = true;
    ImGuiWindowFlags flag = ImGuiWindowFlags_NoCollapse;
    if (ImGui::Begin(on_get_name().c_str(), &is_opened, flag)) {
        // Need to pop var before draw window
        ImGui::PopStyleVar(); // WindowMinSize
        draw_window();
    } else {
        ImGui::PopStyleVar(); // WindowMinSize
    }

    ImGui::End();
    if (!is_opened)
        close();
}

void GLGizmoSVG::on_set_state()
{
    // enable / disable bed from picking
    // Rotation gizmo must work through bed
    m_parent.set_raycaster_gizmos_on_top(GLGizmoBase::m_state == GLGizmoBase::On);

    m_rotate_gizmo.set_state(GLGizmoBase::m_state);

    // Closing gizmo. e.g. selecting another one
    if (GLGizmoBase::m_state == GLGizmoBase::Off) {
        // refuse outgoing during text preview
        if (false) {
            GLGizmoBase::m_state = GLGizmoBase::On;
            auto notification_manager = wxGetApp().plater()->get_notification_manager();
            notification_manager->push_notification(
                NotificationType::CustomNotification,
                NotificationManager::NotificationLevel::RegularNotificationLevel,
                _u8L("ERROR: Wait until ends or Cancel process."));
            return;
        }
        reset_volume();
    } else if (GLGizmoBase::m_state == GLGizmoBase::On) {
        // Try(when exist) set text configuration by volume 
        set_volume_by_selection();
                
        if (!m_gui_cfg.has_value())
            m_set_window_offset = ImGuiWrapper::change_window_position(on_get_name().c_str(), false);
        else
            m_set_window_offset = ImVec2(-1, -1);
        
        // when open by hyperlink it needs to show up
        // or after key 'T' windows doesn't appear
        // m_parent.set_as_dirty();
    }
}

void GLGizmoSVG::data_changed(bool is_serializing) { 
    set_volume_by_selection();
    if (!is_serializing && m_volume == nullptr)
        close();
}

void GLGizmoSVG::on_start_dragging() { m_rotate_gizmo.start_dragging(); }
void GLGizmoSVG::on_stop_dragging()
{
    m_rotate_gizmo.stop_dragging();

    // TODO: when start second rotatiton previous rotation rotate draggers
    // This is fast fix for second try to rotate
    // When fixing, move grabber above text (not on side)
    m_rotate_gizmo.set_angle(PI/2);

    // apply rotation
    m_parent.do_rotate(L("Text-Rotate"));

    m_rotate_start_angle.reset();

    // recalculate for surface cut
    if (m_volume != nullptr && 
        m_volume->emboss_shape.has_value() &&
        m_volume->emboss_shape->projection.use_surface)
        process();
}
void GLGizmoSVG::on_dragging(const UpdateData &data) { m_rotate_gizmo.dragging(data); }

GLGizmoSVG::GuiCfg GLGizmoSVG::create_gui_configuration()
{
    GuiCfg cfg; // initialize by default values;

    
    return cfg;
}

void GLGizmoSVG::set_volume_by_selection()
{
    const Selection &selection = m_parent.get_selection();
    const GLVolume *gl_volume = get_selected_gl_volume(selection);
    if (gl_volume == nullptr)
        return reset_volume();

    const ModelObjectPtrs &objects = selection.get_model()->objects;
    ModelVolume *volume =get_model_volume(*gl_volume, objects);
    if (volume == nullptr)
        return reset_volume();

    // is same volume as actual selected?
    if (volume->id() == m_volume_id)
        return;

    // Do not use focused input value when switch volume(it must swith value)
    if (m_volume != nullptr && 
        m_volume != volume) // when update volume it changed id BUT not pointer
        ImGuiWrapper::left_inputs();

    // is valid svg volume?
    if (!is_svg(*volume)) 
        return reset_volume();

    // cancel previous job
    if (m_job_cancel != nullptr) {
        m_job_cancel->store(true);
        m_job_cancel = nullptr;
    }

    m_volume = volume;
    m_volume_id = volume->id();

    // Calculate current angle of up vector
    m_angle = calc_up(gl_volume->world_matrix(), Slic3r::GUI::up_limit);    

    // calculate scale for height and depth inside of scaled object instance
    calculate_scale();
}

void GLGizmoSVG::reset_volume()
{
    if (m_volume == nullptr)
        return; // already reseted

    m_volume = nullptr;
    m_volume_id.id = 0;
}

void GLGizmoSVG::calculate_scale() {
    Transform3d to_world = m_parent.get_selection().get_first_volume()->world_matrix();
    auto to_world_linear = to_world.linear();
    auto calc = [&to_world_linear](const Vec3d &axe, std::optional<float>& scale)->bool {
        Vec3d  axe_world = to_world_linear * axe;
        double norm_sq   = axe_world.squaredNorm();
        if (is_approx(norm_sq, 1.)) {
            if (scale.has_value())
                scale.reset();
            else
                return false;
        } else {
            scale = sqrt(norm_sq);
        }
        return true;
    };

    calc(Vec3d::UnitY(), m_scale_height);
    calc(Vec3d::UnitZ(), m_scale_depth);
}

bool GLGizmoSVG::process()
{
    // no volume is selected -> selection from right panel
    assert(m_volume != nullptr);
    if (m_volume == nullptr) return false;
        
    //DataUpdate data{create_emboss_data_base(m_text, m_style_manager, m_job_cancel), m_volume->id()};

    //std::unique_ptr<Job> job = nullptr;

    //// check cutting from source mesh
    //bool &use_surface = data.text_configuration.style.prop.use_surface;
    //bool  is_object   = m_volume->get_object()->volumes.size() == 1;
    //if (use_surface && is_object) 
    //    use_surface = false;
    //
    //if (use_surface) {
    //    // Model to cut surface from.
    //    SurfaceVolumeData::ModelSources sources = create_volume_sources(m_volume);
    //    if (sources.empty()) return false;

    //    Transform3d text_tr = m_volume->get_matrix();
    //    auto& fix_3mf = m_volume->text_configuration->fix_3mf_tr;
    //    if (fix_3mf.has_value())
    //        text_tr = text_tr * fix_3mf->inverse();

    //    // when it is new applying of use surface than move origin onto surfaca
    //    if (!m_volume->text_configuration->style.prop.use_surface) {
    //        auto offset = calc_surface_offset(*m_volume, m_raycast_manager, m_parent.get_selection());
    //        if (offset.has_value())
    //            text_tr *= Eigen::Translation<double, 3>(*offset);
    //    }

    //    bool is_outside = m_volume->is_model_part();
    //    // check that there is not unexpected volume type
    //    assert(is_outside || m_volume->is_negative_volume() ||
    //           m_volume->is_modifier());
    //    UpdateSurfaceVolumeData surface_data{std::move(data), {text_tr, is_outside, std::move(sources)}};
    //    job = std::make_unique<UpdateSurfaceVolumeJob>(std::move(surface_data));                  
    //} else {
    //    job = std::make_unique<UpdateJob>(std::move(data));
    //}

    //auto &worker = wxGetApp().plater()->get_ui_job_worker();
    //queue_job(worker, std::move(job));

    //// notification is removed befor object is changed by job
    //remove_notification_not_valid_font();
    return true;
}

void GLGizmoSVG::close()
{
    // close gizmo == open it again
    auto& mng = m_parent.get_gizmos_manager();
    if (mng.get_current_type() == GLGizmosManager::Svg)
        mng.open_gizmo(GLGizmosManager::Svg);
}

void GLGizmoSVG::draw_window()
{
    if (m_volume != nullptr && m_volume->emboss_shape.has_value())
        ImGui::Text("SVG file path is %s", m_volume->emboss_shape->svg_file_path.c_str());

    if (ImGui::Button("change file")) {
        auto data = create_emboss_data_base(m_job_cancel);
        std::string file = choose_svg_file();
    }
    draw_model_type();        
 }

void GLGizmoSVG::draw_model_type()
{
    assert(m_volume != nullptr);
    bool is_last_solid_part = is_svg_object(*m_volume);
    std::string title = _u8L("Type");
    if (is_last_solid_part) {
        ImVec4 color{.5f, .5f, .5f, 1.f};
        m_imgui->text_colored(color, title.c_str());
    } else {
        ImGui::Text("%s", title.c_str());
    }

    std::optional<ModelVolumeType> new_type;
    ModelVolumeType modifier = ModelVolumeType::PARAMETER_MODIFIER;
    ModelVolumeType negative = ModelVolumeType::NEGATIVE_VOLUME;
    ModelVolumeType part = ModelVolumeType::MODEL_PART;
    ModelVolumeType type = m_volume->type();

    if (ImGui::RadioButton(_u8L("Added").c_str(), type == part))
        new_type = part;
    else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Change to object's part.").c_str());
    ImGui::SameLine();

    std::string last_solid_part_hint = _u8L("You can't change a type of the last solid part of the object.");
    if (ImGui::RadioButton(_u8L("Subtracted").c_str(), type == negative))
        new_type = negative;
    else if (ImGui::IsItemHovered()) {
        if (is_last_solid_part)
            ImGui::SetTooltip("%s", last_solid_part_hint.c_str());
        else if (type != negative)
            ImGui::SetTooltip("%s", _u8L("Change to negative volume.").c_str());
    }

    GUI_App &app = wxGetApp();
    Plater  *plater = app.plater();
    // SLA printers do not use modifiers
    bool is_sla = plater->printer_technology() == ptSLA;
    // In simple mode are not modifiers ?!?
    bool is_simple = app.get_mode() == ConfigOptionMode::comSimple;
    if (!is_sla && !is_simple) {
        ImGui::SameLine();
        if (ImGui::RadioButton(_u8L("Modifier").c_str(), type == modifier))
            new_type = modifier;
        else if (ImGui::IsItemHovered()) {
            if (is_last_solid_part)
                ImGui::SetTooltip("%s", last_solid_part_hint.c_str());
            else if (type != modifier)
                ImGui::SetTooltip("%s", _u8L("Change to modifier.").c_str());
        }
    }

    if (m_volume != nullptr && new_type.has_value() && !is_last_solid_part) {
        Plater::TakeSnapshot snapshot(plater, _L("Change Text Type"), UndoRedo::SnapshotType::GizmoAction);
        m_volume->set_type(*new_type);

         // Update volume position when switch from part or into part
        if (m_volume->text_configuration->style.prop.use_surface) {
            // move inside
            bool is_volume_move_inside  = (type == part);
            bool is_volume_move_outside = (*new_type == part);
            if (is_volume_move_inside || is_volume_move_outside) process();
        }

        // inspiration in ObjectList::change_part_type()
        // how to view correct side panel with objects
        ObjectList *obj_list = app.obj_list();
        wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(
            obj_list->get_selected_obj_idx(),
            [volume = m_volume](const ModelVolume *vol) { return vol == volume; });
        if (!sel.IsEmpty()) obj_list->select_item(sel.front());       

        // NOTE: on linux, function reorder_volumes_and_get_selection call GLCanvas3D::reload_scene(refresh_immediately = false)
        // which discard m_volume pointer and set it to nullptr also selection is cleared so gizmo is automaticaly closed
        auto &mng = m_parent.get_gizmos_manager();
        if (mng.get_current_type() != GLGizmosManager::Svg)
            mng.open_gizmo(GLGizmosManager::Svg);
        // TODO: select volume back - Ask @Sasa
    }
}


/////////////
// private namespace implementation
///////////////
namespace {

std::string get_file_name(const std::string &file_path)
{
    if (file_path.empty())
        return file_path;

    size_t pos_last_delimiter = file_path.find_last_of("/\\");
    if (pos_last_delimiter == std::string::npos) {
        // should not happend that in path is not delimiter
        assert(false);
        pos_last_delimiter = 0;
    }

    size_t pos_point = file_path.find_last_of('.');
    if (pos_point == std::string::npos || pos_point < pos_last_delimiter // last point is inside of directory path
    ) {
        // there is no extension
        assert(false);
        pos_point = file_path.size();
    }

    size_t offset = pos_last_delimiter + 1;             // result should not contain last delimiter ( +1 )
    size_t count  = pos_point - pos_last_delimiter - 1; // result should not contain extension point ( -1 )
    return file_path.substr(offset, count);
}

std::string volume_name(const EmbossShape &shape)
{
    std::string file_name = get_file_name(shape.svg_file_path);
    if (!file_name.empty())
        return file_name;
    return "SVG shape";
}

CreateVolumeParams create_input(GLCanvas3D &canvas, RaycastManager& raycaster, ModelVolumeType volume_type)
{
    auto gizmo = static_cast<unsigned char>(GLGizmosManager::Svg);
    const GLVolume *gl_volume = get_first_hovered_gl_volume(canvas);
    Plater *plater = wxGetApp().plater();
    return CreateVolumeParams{canvas, plater->get_camera(), plater->build_volume(),
        plater->get_ui_job_worker(), volume_type, raycaster, gizmo, gl_volume};
}

std::string choose_svg_file()
{
    wxWindow    *parent       = nullptr;
    wxString     message      = _L("Choose SVG file for emboss:");
    wxString     defaultDir   = wxEmptyString;
    wxString     selectedFile = wxEmptyString;
    wxString     wildcard     = file_wildcards(FT_SVG);
    long         style        = wxFD_OPEN | wxFD_FILE_MUST_EXIST;
    wxFileDialog dialog(parent, message, defaultDir, selectedFile, wildcard, style);
    if (dialog.ShowModal() != wxID_OK) {
        BOOST_LOG_TRIVIAL(warning) << "SVG file for emboss was NOT selected.";
        return {};
    }

    wxArrayString input_files;
    dialog.GetPaths(input_files);
    if (input_files.IsEmpty()) {
        BOOST_LOG_TRIVIAL(warning) << "SVG file dialog result is empty.";
        return {};
    }

    if (input_files.size() != 1)
        BOOST_LOG_TRIVIAL(warning) << "SVG file dialog result contain multiple files but only first is used.";

    auto       &input_file = input_files.front();
    std::string path       = std::string(input_file.c_str());
    return path;
}

ExPolygons default_shape()
{
    std::string file  = Slic3r::resources_dir() + "/icons/question.svg";
    NSVGimage  *image = nsvgParseFromFile(file.c_str(), "px", 96.0f);
    ExPolygons  shape = to_expolygons(image);
    nsvgDelete(image);
    return shape;
}

EmbossShape select_shape()
{
    EmbossShape shape;
    shape.projection.depth       = 10.;
    shape.projection.use_surface = false;

    shape.svg_file_path = choose_svg_file();
    if (shape.svg_file_path.empty())
        return {};

    NSVGimage *image = nsvgParseFromFile(shape.svg_file_path.c_str(), "mm", 96.0f);
    ScopeGuard sg([image]() { nsvgDelete(image); });
    shape.shapes = to_expolygons(image);

    // TODO: get scale from file
    shape.scale = 1.;

    // Must contain some shapes !!!
    if (shape.shapes.empty())
        shape.shapes = default_shape();

    return shape;
}

DataBasePtr create_emboss_data_base(std::shared_ptr<std::atomic<bool>> &cancel)
{
    EmbossShape shape = select_shape();

    if (shape.shapes.empty())
        // canceled selection of SVG file
        return nullptr;

    // Cancel previous Job, when it is in process
    // worker.cancel(); --> Use less in this case I want cancel only previous EmbossJob no other jobs
    // Cancel only EmbossUpdateJob no others
    if (cancel != nullptr)
        cancel->store(true);
    // create new shared ptr to cancel new job
    cancel = std::make_shared<std::atomic<bool>>(false);

    std::string name = volume_name(shape);

    return std::make_unique<DataBase>(name, cancel /*copy*/, std::move(shape));
}
} // namespace

// any existing icon filename to not influence GUI
const std::string GLGizmoSVG::M_ICON_FILENAME = "cut.svg";
