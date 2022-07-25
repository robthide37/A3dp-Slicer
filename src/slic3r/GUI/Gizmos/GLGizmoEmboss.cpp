#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp" // to update title when add text
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/MsgDialog.hpp"
#include "slic3r/GUI/format.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "slic3r/GUI/Jobs/EmbossJob.hpp"
#include "slic3r/GUI/Jobs/NotificationProgressIndicator.hpp"
#include "slic3r/Utils/WxFontUtils.hpp"
#include "slic3r/Utils/FontListSerializable.hpp"
#include "slic3r/Utils/UndoRedo.hpp"

// TODO: remove include
#include "libslic3r/SVG.hpp"      // debug store
#include "libslic3r/Geometry.hpp" // covex hull 2d

#include "libslic3r/NSVGUtils.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/ClipperUtils.hpp" // union_ex
#include "libslic3r/AppConfig.hpp"    // store/load font list
#include "libslic3r/Format/OBJ.hpp" // load obj file for default object
#include "libslic3r/BuildVolume.hpp"

#include "imgui/imgui_stdlib.h" // using std::string for inputs
#include "nanosvg/nanosvg.h"    // load SVG file

#include <wx/font.h>
#include <wx/fontutil.h>
#include <wx/fontdlg.h>
#include <wx/fontenum.h>

#include <GL/glew.h>

// uncomment for easier debug
//#define ALLOW_DEBUG_MODE
#ifdef ALLOW_DEBUG_MODE
#define ALLOW_ADD_FONT_BY_FILE
#define ALLOW_ADD_FONT_BY_OS_SELECTOR
#define SHOW_WX_FONT_DESCRIPTOR // OS specific descriptor | file path --> in edit style <tree header>
#define SHOW_FONT_FILE_PROPERTY // ascent, descent, line gap, cache --> in advanced <tree header>
#define SHOW_FONT_COUNT // count of enumerated font --> in font combo box
#define SHOW_CONTAIN_3MF_FIX // when contain fix matrix --> show gray '3mf' next to close button
#define SHOW_OFFSET_DURING_DRAGGING // when drag with text over surface visualize used center
#define SHOW_IMGUI_ATLAS
#define SHOW_ICONS_TEXTURE
#define SHOW_FINE_POSITION
#define SHOW_WX_WEIGHT_INPUT
#define DRAW_PLACE_TO_ADD_TEXT
#define ALLOW_REVERT_ALL_STYLES
#endif // ALLOW_DEBUG_MODE

//#define SHOW_WX_FONT_DESCRIPTOR
//#define SHOW_FONT_FILE_PROPERTY
//#define SHOW_FONT_COUNT
#define SHOW_CONTAIN_3MF_FIX
//#define ALLOW_ADD_FONT_BY_FILE
//#define ALLOW_ADD_FONT_BY_OS_SELECTOR
//#define ALLOW_REVERT_ALL_STYLES

using namespace Slic3r;
using namespace Slic3r::GUI;

// anonymous namespace for unique names
namespace {
template<typename T>
struct MinMax
{
    T min;
    T max;
};
template<typename T>
struct Limit
{
    MinMax<T> gui;
    MinMax<T> values;
};
struct Limits
{
    MinMax<float> emboss{0.01f, 1e4f};
    MinMax<float> size_in_mm{0.1f, 1000.f};
    Limit<float> boldness{{-200.f, 200.f}, {-2e4f, 2e4f}};
    Limit<float> skew{{-1.f, 1.f}, {-100.f, 100.f}};
    MinMax<int>  char_gap{-20000, 20000};
    MinMax<int>  line_gap{-20000, 20000};
    // distance text object from surface
    MinMax<float> angle{-180.f, 180.f}; // in mm

    template<typename T>
    static bool apply(std::optional<T> &val, const MinMax<T> &limit) {
        if (val.has_value())
            return apply<T>(*val, limit);
        return false;
    }
    template<typename T>
    static bool apply(T &val, const MinMax<T> &limit)
    {
        if (val > limit.max) {
            val = limit.max;
            return true;
        }
        if (val < limit.min) {
            val = limit.min;
            return true;
        }
        return false;
    }
};
static const Limits limits;

static bool is_text_empty(const std::string &text){
    return text.empty() ||
           text.find_first_not_of(" \n\t\r") == std::string::npos;
}
} // namespace

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -2)
    , m_volume(nullptr)
    , m_exist_notification(false)
    , m_rotate_gizmo(parent, GLGizmoRotate::Axis::Z) // grab id = 2 (Z axis)
    , m_font_manager(m_imgui->get_glyph_ranges())
    , m_update_job_cancel(nullptr)
    , m_allow_update_rendered_font(false)
{
    m_rotate_gizmo.set_group_id(0);
    // TODO: add suggestion to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web    
}

void GLGizmoEmboss::set_fine_position()
{
    const Selection &selection = m_parent.get_selection();
    const Selection::IndicesList indices   = selection.get_volume_idxs();
    // no selected volume
    if (indices.empty()) return;
    const GLVolume *volume = selection.get_volume(*indices.begin());
    // bad volume selected (e.g. deleted one)
    if (volume == nullptr) return;

    const Camera &camera = wxGetApp().plater()->get_camera();
    Polygon hull = CameraUtils::create_hull2d(camera, *volume);

    const ImVec2 &windows_size = get_minimal_window_size();
    Size          c_size       = m_parent.get_canvas_size();
    ImVec2 canvas_size(c_size.get_width(), c_size.get_height());
    ImVec2 offset = ImGuiWrapper::suggest_location(windows_size, hull, canvas_size);
    m_set_window_offset = offset;
    return;

    Polygon rect({Point(offset.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y + windows_size.y),
                  Point(offset.x, offset.y + windows_size.y)});
    ImGuiWrapper::draw(hull);
    ImGuiWrapper::draw(rect);
}

void GLGizmoEmboss::create_volume(ModelVolumeType volume_type, const Vec2d& mouse_pos)
{
    assert(volume_type == ModelVolumeType::MODEL_PART ||
           volume_type == ModelVolumeType::NEGATIVE_VOLUME ||
           volume_type == ModelVolumeType::PARAMETER_MODIFIER);
    if (!m_gui_cfg.has_value()) initialize();
    set_default_text();

    Vec2d screen_coor = mouse_pos;
    if (mouse_pos.x() < 0 || mouse_pos.y() < 0) {
        // use center of screen
        auto screen_size = m_parent.get_canvas_size();
        screen_coor.x()  = screen_size.get_width() / 2.;
        screen_coor.y()  = screen_size.get_height() / 2.;        
    }

    if (start_volume_creation(volume_type, screen_coor)) return;

    // start creation of new object
    Plater* plater = wxGetApp().plater();
    const Camera &camera = plater->get_camera();
    const Pointfs &bed_shape = plater->build_volume().bed_shape();
    auto &worker = plater->get_ui_job_worker();

    Emboss::FontFileWithCache font;
    if (m_font_manager.is_activ_font())
        font = m_font_manager.get_activ_style().font_file_with_cache;

    EmbossDataBase base{font, create_configuration(), create_volume_name()};
    EmbossDataCreateObject data{std::move(base),
                                screen_coor,
                                camera,
                                bed_shape};
    queue_job(worker, std::make_unique<EmbossCreateObjectJob>(std::move(data)));
}

bool GLGizmoEmboss::on_mouse_for_rotation(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Moving()) return false;

    bool used = use_grabbers(mouse_event);
    if (!m_dragging) return used;

    assert(m_volume != nullptr);
    assert(m_volume->text_configuration.has_value());

    static std::optional<float> start_angle;
    if (mouse_event.Dragging()) {
        auto &angle_opt = m_volume->text_configuration->font_item.prop.angle;
        if (!start_angle.has_value()) {
            start_angle = angle_opt.has_value() ? *angle_opt : 0.f;    
        }
        // temporary rotation
        TransformationType transformation_type(TransformationType::Local);
        double angle = m_rotate_gizmo.get_angle();
        m_parent.get_selection().rotate(Vec3d(0., 0., angle), transformation_type);

        // propagate angle into property
        angle_opt = static_cast<float>(*start_angle + angle);
        // move to range <-M_PI, M_PI>
        if (*angle_opt > M_PI || *angle_opt < -M_PI) {
            int count = static_cast<int>(std::round(*angle_opt / (2 * M_PI)));
            *angle_opt -= static_cast<float>(count * 2 * M_PI);
        }
        // do not store zero
        if (is_approx(*angle_opt, 0.f))
            angle_opt.reset();
        
        // set into activ style
        if (m_font_manager.is_activ_font()) {
            m_font_manager.get_font_prop().angle = angle_opt;
        }
    } else if (mouse_event.LeftUp()) {
        // apply rotation
        m_parent.do_rotate(L("Text-Rotate"));
        start_angle.reset();
    }
    return used;
}

static Vec2d calc_mouse_to_center_text_offset(const Vec2d& mouse, const ModelVolume& mv) {
    const Transform3d &volume_tr   = mv.get_matrix();
    const Camera      &camera      = wxGetApp().plater()->get_camera();
    assert(mv.text_configuration.has_value());

    auto calc_offset = [&mouse, &volume_tr, &camera, &mv]
    (const Transform3d &instrance_tr) -> Vec2d {
        Transform3d to_world = instrance_tr * volume_tr;  

        // Use fix of .3mf loaded tranformation when exist
        if (mv.text_configuration->fix_3mf_tr.has_value())
            to_world = to_world * (*mv.text_configuration->fix_3mf_tr);        
        // zero point of volume in world coordinate system
        Vec3d volume_center = to_world.translation();
        // screen coordinate of volume center
        Vec2i coor = CameraUtils::project(camera, volume_center);
        return coor.cast<double>() - mouse;
    };

    auto object = mv.get_object();
    assert(!object->instances.empty());
    // Speed up for one instance
    if (object->instances.size() == 1) 
        return calc_offset(object->instances.front()->get_matrix());

    Vec2d  nearest_offset;
    double nearest_offset_size = std::numeric_limits<double>::max();
    for (const ModelInstance *instance : object->instances) {        
        Vec2d offset = calc_offset(instance->get_matrix());
        double offset_size = offset.norm();
        if (nearest_offset_size < offset_size) continue;
        nearest_offset_size = offset_size;
        nearest_offset      = offset;
    }
    return nearest_offset;
}

bool GLGizmoEmboss::on_mouse_for_translate(const wxMouseEvent &mouse_event)
{
    // filter events
    if (!(mouse_event.Dragging() && mouse_event.LeftIsDown()) && 
        !mouse_event.LeftUp() &&
        !mouse_event.LeftDown())
        return false;

    // text volume must be selected
    if (m_volume == nullptr) return false;

    // must exist hover object
    int hovered_id = m_parent.get_first_hover_volume_idx();
    if (hovered_id < 0) return false;

    GLVolume *gl_volume = m_parent.get_volumes().volumes[hovered_id];
    const ModelObjectPtrs &objects = wxGetApp().plater()->model().objects;
    ModelVolume *act_model_volume = get_model_volume(gl_volume, objects);

    // hovered object must be actual text volume
    if (m_volume != act_model_volume) return false;

    const ModelVolumePtrs &volumes = m_volume->get_object()->volumes;
    std::vector<size_t> allowed_volumes_id;
    if (volumes.size() > 1) {
        allowed_volumes_id.reserve(volumes.size() - 1);
        for (auto &v : volumes) { 
            if (v->id() == m_volume->id()) continue;
            allowed_volumes_id.emplace_back(v->id().id);
        }
    }

    // wxCoord == int --> wx/types.h
    Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
    Vec2d mouse_pos = mouse_coord.cast<double>();        
    RaycastManager::AllowVolumes condition(std::move(allowed_volumes_id));

    // detect start text dragging
    if (mouse_event.LeftDown()) {
        // initialize raycasters
        // IMPROVE: move to job, for big scene it slows down 
        ModelObject *act_model_object = act_model_volume->get_object();
        m_raycast_manager.actualize(act_model_object, &condition);
        m_dragging_mouse_offset = calc_mouse_to_center_text_offset(mouse_pos, *m_volume);
        return false;
    }

    // Dragging starts out of window
    if (!m_dragging_mouse_offset.has_value()) return false;

    const Camera &camera = wxGetApp().plater()->get_camera();
    Vec2d offseted_mouse = mouse_pos + *m_dragging_mouse_offset;
    auto hit = m_raycast_manager.unproject(offseted_mouse, camera, &condition);
    if (!hit.has_value()) { 
        // there is no hit
        // show common translation of object
        m_parent.toggle_model_objects_visibility(true);
        m_temp_transformation = {};
        return false; 
    }
        
    if (mouse_event.Dragging()) {
        TextConfiguration &tc = *m_volume->text_configuration;
        // hide common dragging of object
        m_parent.toggle_model_objects_visibility(false, m_volume->get_object(), gl_volume->instance_idx(), m_volume);

        // Calculate temporary position
        Transform3d object_trmat = m_raycast_manager.get_transformation(hit->tr_key);
        Transform3d trmat = Emboss::create_transformation_onto_surface(hit->position, hit->normal);
        const FontProp& font_prop = tc.font_item.prop;
        Emboss::apply_transformation(font_prop, trmat);

        // fix baked transformation from .3mf store process
        if (tc.fix_3mf_tr.has_value())
            trmat = trmat * (*tc.fix_3mf_tr);

        // temp is in wolrld coors
        m_temp_transformation = object_trmat * trmat;
    } else if (mouse_event.LeftUp()) {
        // TODO: Disable apply common transformation after draggig
        // Call after is used for apply transformation after common dragging to rewrite it
        Transform3d volume_trmat =
            gl_volume->get_instance_transformation().get_matrix().inverse() *
            *m_temp_transformation;
        wxGetApp().plater()->CallAfter([volume_trmat, mv = m_volume]() {
            mv->set_transformation(volume_trmat);
        });

        m_parent.toggle_model_objects_visibility(true);
        // Apply temporary position
        m_temp_transformation = {};
        m_dragging_mouse_offset = {};

        // Update surface by new position
        if (m_volume->text_configuration->font_item.prop.use_surface) {
            // need actual position
            m_volume->set_transformation(volume_trmat);
            process();
        }
    }
    return false;
}

bool GLGizmoEmboss::on_mouse(const wxMouseEvent &mouse_event)
{
    // not selected volume
    if (m_volume == nullptr) return false;

    // do not process moving event
    if (mouse_event.Moving()) return false;
    if (on_mouse_for_rotation(mouse_event)) return true;
    if (on_mouse_for_translate(mouse_event)) return true;

    return false;
}

bool GLGizmoEmboss::on_init()
{
    m_rotate_gizmo.init();
    ColorRGBA gray_color(.6f, .6f, .6f, .3f);
    m_rotate_gizmo.set_highlight_color(gray_color);

    m_shortcut_key = WXK_CONTROL_T;
    return true;
}

std::string GLGizmoEmboss::on_get_name() const { return _u8L("Emboss"); }

void GLGizmoEmboss::on_render() {
    // no volume selected
    if (m_volume == nullptr) return;
    Selection &selection = m_parent.get_selection();
    if (selection.is_empty()) return;

    if (m_temp_transformation.has_value()) {
        // draw text volume on temporary position
#if ENABLE_LEGACY_OPENGL_REMOVAL
        GLVolume& gl_volume = *selection.get_volume(*selection.get_volume_idxs().begin());
#else
        const GLVolume& gl_volume = *selection.get_volume(*selection.get_volume_idxs().begin());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_GL_SHADERS_ATTRIBUTES
        GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
#else
        glsafe(::glPushMatrix());
        glsafe(::glMultMatrixd(m_temp_transformation->data()));                
        GLShaderProgram *shader = wxGetApp().get_shader("gouraud_light");
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        shader->start_using();

#if ENABLE_GL_SHADERS_ATTRIBUTES
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d matrix = camera.get_view_matrix() * (*m_temp_transformation);
        shader->set_uniform("view_model_matrix", matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        shader->set_uniform("normal_matrix", (Matrix3d)matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

        // dragging object must be selected so draw it with correct color
        //auto color = gl_volume.color;
        //auto color = gl_volume.render_color;
        auto color = GLVolume::SELECTED_COLOR;
        // Set transparent color for NEGATIVE_VOLUME & PARAMETER_MODIFIER
        bool is_transparent = m_volume->type() != ModelVolumeType::MODEL_PART;        
        if (is_transparent) {
            color.a(0.5f);
            glsafe(::glEnable(GL_BLEND));
            glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));
        }
#if !ENABLE_LEGACY_OPENGL_REMOVAL
        shader->set_uniform("uniform_color", color);
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

        glsafe(::glEnable(GL_DEPTH_TEST));
#if ENABLE_LEGACY_OPENGL_REMOVAL
        gl_volume.model.set_color(color);
        gl_volume.model.render();
#else
        gl_volume.indexed_vertex_array.render();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
        glsafe(::glDisable(GL_DEPTH_TEST));

        if (is_transparent) glsafe(::glDisable(GL_BLEND));

        shader->stop_using();
#if !ENABLE_GL_SHADERS_ATTRIBUTES
        glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
    }

    // Do NOT render rotation grabbers when dragging object
    bool is_rotate_by_grabbers = m_dragging;
    if (!m_parent.is_dragging() || is_rotate_by_grabbers) {
        glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
        m_rotate_gizmo.render();
    }
}

void GLGizmoEmboss::on_render_for_picking() {
    m_rotate_gizmo.render_for_picking();
}

#ifdef SHOW_FINE_POSITION
// draw suggested position of window
static void draw_fine_position(const Selection &selection,
                               const Size      &canvas,
                               const ImVec2    &windows_size)
{
    const Selection::IndicesList indices = selection.get_volume_idxs();
    // no selected volume
    if (indices.empty()) return;
    const GLVolume *volume = selection.get_volume(*indices.begin());
    // bad volume selected (e.g. deleted one)
    if (volume == nullptr) return;

    const Camera   &camera = wxGetApp().plater()->get_camera();
    Slic3r::Polygon hull   = CameraUtils::create_hull2d(camera, *volume);
    ImVec2          canvas_size(canvas.get_width(), canvas.get_height());
    ImVec2 offset = ImGuiWrapper::suggest_location(windows_size, hull,
                                                   canvas_size);
    Slic3r::Polygon rect(
        {Point(offset.x, offset.y), Point(offset.x + windows_size.x, offset.y),
         Point(offset.x + windows_size.x, offset.y + windows_size.y),
         Point(offset.x, offset.y + windows_size.y)});
    ImGuiWrapper::draw(hull);
    ImGuiWrapper::draw(rect);
}
#endif // SHOW_FINE_POSITION

#ifdef DRAW_PLACE_TO_ADD_TEXT
static void draw_place_to_add_text()
{
    ImVec2             mp = ImGui::GetMousePos();
    Vec2d              mouse_pos(mp.x, mp.y);
    const Camera      &camera = wxGetApp().plater()->get_camera();
    Vec3d              p1 = CameraUtils::get_z0_position(camera, mouse_pos);
    std::vector<Vec3d> rect3d{p1 + Vec3d(5, 5, 0), p1 + Vec3d(-5, 5, 0),
                              p1 + Vec3d(-5, -5, 0), p1 + Vec3d(5, -5, 0)};
    Points             rect2d = CameraUtils::project(camera, rect3d);
    ImGuiWrapper::draw(Slic3r::Polygon(rect2d));
}
#endif // DRAW_PLACE_TO_ADD_TEXT

#ifdef SHOW_OFFSET_DURING_DRAGGING
static void draw_mouse_offset(const std::optional<Vec2d> &offset)
{
    if (!offset.has_value()) return;
    // debug draw
    auto   draw_list = ImGui::GetOverlayDrawList();
    ImVec2 p1        = ImGui::GetMousePos();
    ImVec2 p2(p1.x + offset->x(), p1.y + offset->y());
    ImU32  color     = ImGui::GetColorU32(ImGuiWrapper::COL_ORANGE_LIGHT);
    float  thickness = 3.f;
    draw_list->AddLine(p1, p2, color, thickness);
}
#endif // SHOW_OFFSET_DURING_DRAGGING

void GLGizmoEmboss::on_render_input_window(float x, float y, float bottom_limit)
{
    initialize();
    check_selection(); 

    // TODO: fix width - showing scroll in first draw of advanced.
    const ImVec2 &min_window_size = get_minimal_window_size();
    ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, min_window_size);

#ifdef SHOW_FINE_POSITION
    draw_fine_position(m_parent.get_selection(), m_parent.get_canvas_size(), min_window_size);
#endif // SHOW_FINE_POSITION
#ifdef DRAW_PLACE_TO_ADD_TEXT
    draw_place_to_add_text();
#endif // DRAW_PLACE_TO_ADD_TEXT
#ifdef SHOW_OFFSET_DURING_DRAGGING
    draw_mouse_offset(m_dragging_mouse_offset);
#endif // SHOW_OFFSET_DURING_DRAGGING

    // check if is set window offset
    if (m_set_window_offset.has_value()) {
        ImGui::SetNextWindowPos(*m_set_window_offset, ImGuiCond_Always);
        m_set_window_offset.reset();
    }

    ImGuiWindowFlags flag = //ImGuiWindowFlags_AlwaysAutoResize 
               //ImGuiWindowFlags_NoResize         
               ImGuiWindowFlags_NoCollapse
        ;
    bool is_open = true;
    if (ImGui::Begin(on_get_name().c_str(), &is_open, flag)) {
        // Need to pop var before draw window
        ImGui::PopStyleVar(); // WindowMinSize
        draw_window();
    } else {
        ImGui::PopStyleVar(); // WindowMinSize
    }
    ImGui::End();

    // close button in header was hit
    if (!is_open) { 
        close();     
    }
}

void GLGizmoEmboss::on_set_state()
{
    // set manipulator to be able to rotate with text
    ObjectManipulation *manipul = wxGetApp().obj_manipul();
    static ECoordinatesType prev_coordinate_type = ECoordinatesType::World;
    if (GLGizmoBase::m_state == GLGizmoBase::Off)
        manipul->set_coordinates_type(prev_coordinate_type); // set previous state
    else if (GLGizmoBase::m_state == GLGizmoBase::On) {
        prev_coordinate_type = manipul->get_coordinates_type();
        manipul->set_coordinates_type(ECoordinatesType::Local);
    }

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
        m_volume = nullptr;

        remove_notification_not_valid_font();
    } else if (GLGizmoBase::m_state == GLGizmoBase::On) {
        if (!m_gui_cfg.has_value()) initialize();

        // to reload fonts from system, when install new one
        wxFontEnumerator::InvalidateCache();

        // Try(when exist) set configuration by volume
        load_configuration(get_selected_volume());

        // change position of just opened emboss window
        set_fine_position();

        // when open by hyperlink it needs to show up
        // or after key 'T' windows doesn't appear
        m_parent.set_as_dirty();
    }
}

void GLGizmoEmboss::on_start_dragging() { m_rotate_gizmo.start_dragging(); }
void GLGizmoEmboss::on_stop_dragging()
{
    m_rotate_gizmo.stop_dragging();

    // TODO: when start second rotatiton previous rotation rotate draggers
    // This is fast fix for second try to rotate
    // When fixing, move grabber above text (not on side)
    m_rotate_gizmo.set_angle(0);
}
void GLGizmoEmboss::on_dragging(const UpdateData &data) { m_rotate_gizmo.dragging(data); }

void GLGizmoEmboss::initialize()
{
    if (m_gui_cfg.has_value()) return;

    GuiCfg cfg; // initialize by default values;    
    float line_height = ImGui::GetTextLineHeight();
    float line_height_with_spacing = ImGui::GetTextLineHeightWithSpacing();
    float space = line_height_with_spacing - line_height;

    cfg.max_font_name_width = ImGui::CalcTextSize("Maximal font name, extended").x;

    cfg.icon_width = std::ceil(line_height);
    // make size pair number
    if (cfg.icon_width % 2 != 0) ++cfg.icon_width;

    cfg.delete_pos_x = cfg.max_font_name_width + space;
    int count_line_of_text = 3;
    cfg.text_size = ImVec2(-FLT_MIN, line_height_with_spacing * count_line_of_text);
    ImVec2 letter_m_size = ImGui::CalcTextSize("M");
    int count_letter_M_in_input = 12;
    cfg.input_width = letter_m_size.x * count_letter_M_in_input;
    GuiCfg::Translations &tr = cfg.translations;
    tr.type  = _u8L("Type");
    tr.style = _u8L("Style");
    float max_style_text_width = std::max(
        ImGui::CalcTextSize(tr.type.c_str()).x,
        ImGui::CalcTextSize(tr.style.c_str()).x);
    cfg.style_offset = max_style_text_width + 3 * space;

    tr.font  = _u8L("Font");
    tr.size  = _u8L("Height");
    tr.depth = _u8L("Depth");
    float max_text_width = std::max({
        ImGui::CalcTextSize(tr.font.c_str()).x,
        ImGui::CalcTextSize(tr.size.c_str()).x,
        ImGui::CalcTextSize(tr.depth.c_str()).x});
    cfg.input_offset = max_text_width + 3 * space + ImGui::GetTreeNodeToLabelSpacing();

    tr.use_surface = _u8L("Use surface");
    tr.char_gap = _u8L("Char gap");
    tr.line_gap = _u8L("Line gap");
    tr.boldness = _u8L("Boldness");
    tr.italic = _u8L("Skew ratio");
    tr.surface_distance = _u8L("Z-move");
    tr.angle = _u8L("Z-rot");
    tr.collection = _u8L("Collection");
    float max_advanced_text_width = std::max({
        ImGui::CalcTextSize(tr.use_surface.c_str()).x,
        ImGui::CalcTextSize(tr.char_gap.c_str()).x,
        ImGui::CalcTextSize(tr.line_gap.c_str()).x,
        ImGui::CalcTextSize(tr.boldness.c_str()).x,
        ImGui::CalcTextSize(tr.italic.c_str()).x,
        ImGui::CalcTextSize(tr.surface_distance.c_str()).x,
        ImGui::CalcTextSize(tr.angle.c_str()).x,
        ImGui::CalcTextSize(tr.collection.c_str()).x });
    cfg.advanced_input_offset = 
        3 * space + 2* ImGui::GetTreeNodeToLabelSpacing() + max_advanced_text_width;

    // calculate window size
    const ImGuiStyle &style = ImGui::GetStyle();
    float window_title = line_height + 2*style.FramePadding.y;
    float input_height = line_height_with_spacing + 2*style.FramePadding.y;
    float tree_header  = line_height_with_spacing;
    float window_height = 
        window_title + // window title
        cfg.text_size.y +  // text field
        input_height * 6 + // type Radios + style selector + font name +
                           // height + depth + close button
        tree_header +      // advance tree
        2 * style.WindowPadding.y;
    float window_width = cfg.input_offset + cfg.input_width + style.WindowPadding.x * 2;
    cfg.minimal_window_size = ImVec2(window_width, window_height);

    // 6 = charGap, LineGap, Bold, italic, surfDist, angle
    // 4 = 1px for fix each edit image of drag float 
    float advance_height = input_height * 7 + 8;
    cfg.minimal_window_size_with_advance =
        ImVec2(cfg.minimal_window_size.x,
               cfg.minimal_window_size.y + advance_height);

    int max_style_image_width = cfg.max_font_name_width -
                                2 * style.FramePadding.x;
    int max_style_image_height = 2 * input_height;
    cfg.max_style_image_size = Vec2i(max_style_image_width, max_style_image_height);
    cfg.face_name_size.y() = line_height_with_spacing;

    m_gui_cfg.emplace(std::move(cfg));

    init_icons();
    
    const AppConfig *app_cfg = wxGetApp().app_config;
    size_t activ_index = -1;
    FontList font_list = load_font_list_from_app_config(app_cfg, activ_index);
    m_font_manager.add_fonts(font_list);
    if (!m_font_manager.load_font(activ_index) &&
        !m_font_manager.load_first_valid_font()) {
        m_font_manager.add_fonts(create_default_font_list());
        // When default fonts are not loadable use packed one
        if (!m_font_manager.load_first_valid_font()) {
            m_font_manager.add_font(create_default_font());
            [[maybe_unused]]
            bool success = m_font_manager.load_font(0);
            assert(success);
        }
    }
    fill_stored_font_items();
    select_stored_font_item();
    set_default_text();
}

FontList GLGizmoEmboss::create_default_font_list()
{
    // https://docs.wxwidgets.org/3.0/classwx_font.html
    // Predefined objects/pointers: wxNullFont, wxNORMAL_FONT, wxSMALL_FONT, wxITALIC_FONT, wxSWISS_FONT
    return {
        WxFontUtils::get_font_item(*wxNORMAL_FONT, _u8L("NORMAL")), // wxSystemSettings::GetFont(wxSYS_DEFAULT_GUI_FONT)
        WxFontUtils::get_font_item(*wxSMALL_FONT, _u8L("SMALL")),  // A font using the wxFONTFAMILY_SWISS family and 2 points smaller than wxNORMAL_FONT.
        WxFontUtils::get_font_item(*wxITALIC_FONT, _u8L("ITALIC")), // A font using the wxFONTFAMILY_ROMAN family and wxFONTSTYLE_ITALIC style and of the same size of wxNORMAL_FONT.
        WxFontUtils::get_font_item(*wxSWISS_FONT, _u8L("SWISS")),  // A font identic to wxNORMAL_FONT except for the family used which is wxFONTFAMILY_SWISS.
        WxFontUtils::get_font_item(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD), _u8L("MODERN"))
        //, WxFontUtils::get_os_font() == wxNORMAL_FONT
    };
}

FontItem GLGizmoEmboss::create_default_font() {
    std::string font_path = Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf";
    return FontItem{"Default font", font_path, FontItem::Type::file_path};
}

void GLGizmoEmboss::set_default_text()
{
    m_text = _u8L("Embossed text");
}

bool GLGizmoEmboss::start_volume_creation(ModelVolumeType volume_type,
                                          const Vec2d    &screen_coor)
{
    Plater* plater = wxGetApp().plater();
    
    const Selection &selection = m_parent.get_selection();
    if (selection.is_empty()) return false;

    int hovered_id_signed = m_parent.get_first_hover_volume_idx();
    if (hovered_id_signed < 0) return false;

    size_t hovered_id = static_cast<size_t>(hovered_id_signed);
    auto &volumes = m_parent.get_volumes().volumes;
    if (hovered_id >= volumes.size()) return false;  
        
    int object_idx_signed = selection.get_object_idx(); 
    if (object_idx_signed < 0) return false;

    size_t object_idx = static_cast<size_t>(object_idx_signed);
    auto &objects = plater->model().objects;
    if (object_idx >= objects.size()) return false;

    ModelObject *obj = objects[object_idx];
    m_raycast_manager.actualize(obj);

    const Camera &camera = plater->get_camera();
    std::optional<RaycastManager::Hit> hit =
        m_raycast_manager.unproject(screen_coor, camera);
        
    // context menu for add text could be open only by right click on an
    // object. After right click, object is selected and object_idx is set
    // also hit must exist. But there is proper behavior when hit doesn't
    // exists. When this assert appear distquish remove of it.
    assert(hit.has_value());
    if (!hit.has_value()) return false;
        
    Transform3d hit_object_trmat = m_raycast_manager.get_transformation(hit->tr_key);
    
    GLVolume *gl_volume = volumes[hovered_id];
    Transform3d hit_instance_trmat = gl_volume->get_instance_transformation().get_matrix();

    // create volume
    auto &worker = plater->get_ui_job_worker(); 
    
    const Emboss::FontFileWithCache& font = m_font_manager.get_activ_style().font_file_with_cache;
    EmbossDataBase base{font, create_configuration(), create_volume_name()};
    EmbossDataCreateVolume data{std::move(base),
                                volume_type,
                                screen_coor,
                                object_idx_signed,
                                camera,
                                *hit,
                                hit_object_trmat,
                                hit_instance_trmat};
    queue_job(worker, std::make_unique<EmbossCreateVolumeJob>(std::move(data)));
    return true;
}

#include "imgui/imgui_internal.h" // to unfocus input --> ClearActiveID
void GLGizmoEmboss::check_selection()
{
    ModelVolume *vol = get_selected_volume();
    // is same volume selected?
    if (vol != nullptr && m_volume == vol) return;

    // for changed volume notification is NOT valid
    remove_notification_not_valid_font();

    // Do not use focused input value when switch volume(it must swith value)
    if (m_volume != nullptr) ImGui::ClearActiveID();

    // is select embossed volume?
    if (load_configuration(vol)) 
        // successfull load volume for editing
        return;
    
    // behave like adding new text
    m_volume = nullptr;
    set_default_text();
}

ModelVolume *GLGizmoEmboss::get_selected_volume()
{
    return get_selected_volume(m_parent.get_selection(),
                               wxGetApp().plater()->model().objects);
}

ModelVolume *GLGizmoEmboss::get_model_volume(const GLVolume *      gl_volume,
                                             const ModelObjectPtrs& objects)
{
    const GLVolume::CompositeID &id = gl_volume->composite_id;

    if (id.object_id < 0 ||
        static_cast<size_t>(id.object_id) >= objects.size())
        return nullptr;
    ModelObject *object = objects[id.object_id];

    if (id.volume_id < 0 ||
        static_cast<size_t>(id.volume_id) >= object->volumes.size())
        return nullptr;
    return object->volumes[id.volume_id];
}

ModelVolume *GLGizmoEmboss::get_selected_volume(const Selection &selection,
                                                const ModelObjectPtrs& objects)
{
    int object_idx = selection.get_object_idx();
    // is more object selected?
    if (object_idx == -1) return nullptr;

    auto volume_idxs = selection.get_volume_idxs();
    // is more volumes selected?
    if (volume_idxs.size() != 1) return nullptr;
    unsigned int                 vol_id_gl = *volume_idxs.begin();
    const GLVolume *             vol_gl    = selection.get_volume(vol_id_gl);
    return get_model_volume(vol_gl, objects);
}

// Run Job on main thread (blocking) - ONLY DEBUG
static inline void execute_job(std::shared_ptr<Job> j)
{
    struct MyCtl : public Job::Ctl
    {
        void update_status(int st, const std::string &msg = "") override{};
        bool was_canceled() const override { return false; }
        std::future<void> call_on_main_thread(std::function<void()> fn) override
        {
            return std::future<void>{};
        }
    } ctl;
    j->process(ctl);
    wxGetApp().plater()->CallAfter([j]() {
        std::exception_ptr e_ptr = nullptr;
        j->finalize(false, e_ptr);
    });
}

bool GLGizmoEmboss::process()
{
    // no volume is selected -> selection from right panel
    if (m_volume == nullptr) return false;

    // without text there is nothing to emboss
    if (m_text.empty()) return false;

    // exist loaded font
    if (!m_font_manager.is_activ_font()) return false;
    Emboss::FontFileWithCache font = m_font_manager.get_activ_style().font_file_with_cache;
    assert(font.has_value());
    if (!font.has_value()) return false;    
    
    // Can't use cancel, because I want cancel only previous EmbossUpdateJob
    // worker.cancel();
    
    // Cancel only EmbossUpdateJob no others
    if (m_update_job_cancel != nullptr)
        m_update_job_cancel->store(true);
    // create new shared ptr to cancel new job
    m_update_job_cancel = std::make_shared<std::atomic<bool> >(false);
    EmbossDataBase   base{font, create_configuration(), create_volume_name()};
    EmbossDataUpdate data{std::move(base), m_volume->id(), m_update_job_cancel};

    std::unique_ptr<Job> job = nullptr;

    // check cutting from source mesh
    const TextConfiguration &tc = data.text_configuration;
    if (tc.font_item.prop.use_surface) {
        // Model to cut surface from.
        UseSurfaceData::ModelSources sources = UseSurfaceData::create_sources(m_volume);
        if (sources.empty()) return false;

        Transform3d text_tr = m_volume->get_matrix();
        auto& fix_3mf = m_volume->text_configuration->fix_3mf_tr;
        if (fix_3mf.has_value())
            text_tr = text_tr * fix_3mf->inverse();

        bool is_outside = m_volume->is_model_part();
        // check that there is not unexpected volume type
        assert(is_outside || m_volume->is_negative_volume() ||
               m_volume->is_modifier());
        UseSurfaceData surface_data{std::move(data), text_tr, is_outside,
                                    std::move(sources)};
        job = std::make_unique<UseSurfaceJob>(std::move(surface_data));                  
    } else {
        job = std::make_unique<EmbossUpdateJob>(std::move(data));
    }

    //*    
    auto &worker = wxGetApp().plater()->get_ui_job_worker();
    queue_job(worker, std::move(job));
    /*/ // Run Job on main thread (blocking) - ONLY DEBUG
    execute_job(std::move(job));
    // */

    // notification is removed befor object is changed by job
    remove_notification_not_valid_font();
    return true;
}

void GLGizmoEmboss::close()
{
    // remove volume when text is empty
    if (m_volume != nullptr && 
        m_volume->text_configuration.has_value() &&
        is_text_empty(m_text)) {
        Plater &p = *wxGetApp().plater();
        if (is_text_object(m_volume)) {
            // delete whole object
            p.remove(m_parent.get_selection().get_object_idx());
        } else {
            // delete text volume
            p.remove_selected();
        }
    }

    // close gizmo == open it again
    auto& mng = m_parent.get_gizmos_manager();
    if (mng.get_current_type() == GLGizmosManager::Emboss)
        mng.open_gizmo(GLGizmosManager::Emboss);
}

void GLGizmoEmboss::fill_stored_font_items()
{
    m_stored_font_items.clear();
    for (const auto &item : m_font_manager.get_styles()) {
        const FontItem &fi = item.font_item;
        // skip file paths + fonts from other OS(loaded from .3mf)
        if (fi.type != WxFontUtils::get_actual_type()) continue;

        assert(m_stored_font_items.find(fi.name) == m_stored_font_items.end());
        m_stored_font_items[fi.name] = fi; // copy
    }
    select_stored_font_item();
}

void GLGizmoEmboss::select_stored_font_item()
{    
    const std::string &name = m_font_manager.get_font_item().name;
    const auto &it = m_stored_font_items.find(name);
    if (it == m_stored_font_items.end()) { 
        m_stored_font_item.reset();
        m_stored_wx_font.reset();
        return;
    }
    m_stored_font_item = it->second;
    m_stored_wx_font = WxFontUtils::load_wxFont(m_stored_font_item->path);
}

void GLGizmoEmboss::draw_window()
{
#ifdef ALLOW_DEBUG_MODE
    if (ImGui::Button("re-process")) process();
    if (ImGui::Button("add svg")) choose_svg_file();
    if (ImGui::Button("use system font")) {
        size_t font_index = m_font_list.size();
        m_font_list.emplace_back(WxFontUtils::get_os_font());
        bool loaded = load_font(font_index);
    }
#endif //  ALLOW_DEBUG_MODE

    bool is_selected_style = m_font_manager.is_activ_font();
    bool exist_font_file   = is_selected_style &&
                           m_font_manager.get_font_file() != nullptr;
    if (!is_selected_style) {
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Warning: No text style is selected."));
    }else if (!exist_font_file) {
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _L("Warning: No font is selected. Select correct one."));
    }
    draw_text_input();
    draw_model_type();
    draw_style_list();
    m_imgui->disabled_begin(!is_selected_style);
    ImGui::TreePush();
    draw_style_edit();
    ImGui::TreePop();
    if (ImGui::TreeNode(_u8L("advanced").c_str())) {
        if (!m_is_advanced_edit_style) {
            set_minimal_window_size(true);
        } else {
            draw_advanced();
        }
        ImGui::TreePop();
    } else if (m_is_advanced_edit_style) 
        set_minimal_window_size(false);
    m_imgui->disabled_end(); // !is_selected_style
       
#ifdef SHOW_WX_FONT_DESCRIPTOR
    if (is_selected_style)
        m_imgui->text_colored(ImGuiWrapper::COL_GREY_DARK, m_font_manager.get_font_item().path);
#endif // SHOW_WX_FONT_DESCRIPTOR
        
    if (ImGui::Button(_u8L("Close").c_str())) close();

    // Option to create text volume when reselecting volumes
    m_imgui->disabled_begin(!exist_font_file);
    if (m_volume == nullptr) {
        ImGui::SameLine();
        if (ImGui::Button(_u8L("Generate object").c_str()))
            create_volume(ModelVolumeType::MODEL_PART);
    }
    m_imgui->disabled_end();

#ifdef SHOW_CONTAIN_3MF_FIX
    if (m_volume!=nullptr &&
        m_volume->text_configuration.has_value() &&
        m_volume->text_configuration->fix_3mf_tr.has_value()) {
        ImGui::SameLine();
        m_imgui->text_colored(ImGuiWrapper::COL_GREY_DARK, ".3mf");
        if (ImGui::IsItemHovered()) {
            Transform3d &fix = *m_volume->text_configuration->fix_3mf_tr;
            std::stringstream ss;
            ss << fix.matrix();            
            std::string filename = (m_volume->source.input_file.empty())? "unknown.3mf" :
                                   m_volume->source.input_file + ".3mf";
            ImGui::SetTooltip("Text configuation contain \n"
                              "Fix Transformation Matrix \n"
                              "%s\n"
                              "loaded from \"%s\" file.",
                              ss.str().c_str(), filename.c_str()
                );
        }
    }
#endif // SHOW_CONTAIN_3MF_FIX
#ifdef SHOW_ICONS_TEXTURE    
    auto &t = m_icons_texture;
    ImGui::Image((void *) t.get_id(), ImVec2(t.get_width(), t.get_height()));
#endif //SHOW_ICONS_TEXTURE
#ifdef SHOW_IMGUI_ATLAS
    auto &atlas = m_font_manager.m_imgui_font_atlas;
    ImGui::Image(atlas.TexID, ImVec2(atlas.TexWidth, atlas.TexHeight));
#endif // SHOW_IMGUI_ATLAS
}

void GLGizmoEmboss::draw_text_input()
{
    auto create_range_text = [&manager = m_font_manager, &text = m_text, &exist_unknown = m_text_contain_unknown_glyph]() {
        auto &font_file = manager.get_font_file();
        const auto  &cn = manager.get_font_prop().collection_number;
        unsigned int font_index = (cn.has_value()) ? *cn : 0;
        return Emboss::create_range_text(text, *font_file, font_index, &exist_unknown);
    };

    static const ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_AllowTabInput | ImGuiInputTextFlags_AutoSelectAll;

    ImFont *imgui_font = m_font_manager.get_imgui_font();
    if (imgui_font == nullptr) {
        // try create new imgui font
        m_font_manager.create_imgui_font(create_range_text());
        imgui_font = m_font_manager.get_imgui_font();
    }
    bool exist_font = imgui_font != nullptr;
    assert(!exist_font || imgui_font->IsLoaded());
    if (exist_font) ImGui::PushFont(imgui_font);

    // flag for extend font ranges if neccessary
    // ranges can't be extend during font is activ(pushed)
    std::string range_text;
    float  window_height  = ImGui::GetWindowHeight();
    float  minimal_height = get_minimal_window_size().y;
    float  extra_height   = window_height - minimal_height;
    ImVec2 text_size(m_gui_cfg->text_size.x,
                     m_gui_cfg->text_size.y + extra_height);
    if (ImGui::InputTextMultiline("##Text", &m_text, text_size, flags)) {
        process();
        range_text = create_range_text();
    }

    if (exist_font) ImGui::PopFont();

    // show warning about incorrectness view of font
    std::string warning;
    std::string tool_tip;
    if (!exist_font) {
        warning = _u8L("Can't write text by selected font.");
        tool_tip = _u8L("Try to choose another font.");
    } else {
        std::string who;
        auto append_warning = [&who, &tool_tip](std::string w, std::string t) {
            if (!w.empty()) {
                if (!who.empty()) who += " & ";
                who += w;
            }
            if (!t.empty()) {
                if (!tool_tip.empty()) tool_tip += "\n";
                tool_tip += t;            
            }
        };
        if (is_text_empty(m_text))
            append_warning(_u8L("Empty"), 
                _u8L("Embossed text can NOT contain only white spaces."));
        if (m_text_contain_unknown_glyph)
            append_warning(_u8L("Bad symbol"), 
                _u8L("Text contain character glyph (represented by '?') unknown by font."));

        const FontProp &prop = m_font_manager.get_font_prop();
        if (prop.skew.has_value())
            append_warning(_u8L("Skew"), 
                _u8L("Unsupported visualization of font skew for text input."));
        if (prop.boldness.has_value())
            append_warning(_u8L("Boldness"),
                _u8L("Unsupported visualization of font boldness for text input."));
        if (prop.line_gap.has_value())
            append_warning(_u8L("Line gap"),
                _u8L("Unsupported visualization of gap between lines inside text input."));
        float imgui_size = FontManager::get_imgui_font_size(prop, *m_font_manager.get_font_file());
        if (imgui_size > FontManager::max_imgui_font_size)           
            append_warning(_u8L("To tall"),
                _u8L("Diminished font height inside text input."));
        if (imgui_size < FontManager::min_imgui_font_size)       
            append_warning(_u8L("To small"),
                _u8L("Enlarged font height inside text input."));
        if (!who.empty())
            warning = GUI::format(_u8L("%1% is NOT shown."), who);
    }

    if (!warning.empty()) {
        if (ImGui::IsItemHovered() && !tool_tip.empty())
            ImGui::SetTooltip("%s", tool_tip.c_str());
        ImVec2 cursor = ImGui::GetCursorPos();
        float width = ImGui::GetContentRegionAvailWidth();
        ImVec2 size = ImGui::CalcTextSize(warning.c_str());
        ImVec2 padding = ImGui::GetStyle().FramePadding;
        ImGui::SetCursorPos(ImVec2(width - size.x + padding.x,
                                   cursor.y - size.y - padding.y));
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, warning);
        ImGui::SetCursorPos(cursor);
    }

    // IMPROVE: only extend not clear
    // Extend font ranges
    // imgui_font has to be unused
    if (!range_text.empty() &&
        !m_imgui->contain_all_glyphs(imgui_font, range_text) ) { 
        m_font_manager.clear_imgui_font(); 
        m_font_manager.create_imgui_font(range_text);
    }
}

//#define DEBUG_NOT_LOADABLE_FONTS

/// <summary>
/// Keep list of loadable OS fonts
/// Filtrate which can be loaded.
/// Sort alphanumerical.
/// </summary>
class MyFontEnumerator : public wxFontEnumerator
{
    wxFontEncoding m_encoding;
public: 
    std::vector<wxString> m_facenames;
    MyFontEnumerator(wxFontEncoding encoding) : m_encoding(encoding) {}

#ifdef DEBUG_NOT_LOADABLE_FONTS
    std::vector<std::string> m_efacenames;
#endif // DEBUG_NOT_LOADABLE_FONTS

protected:
    /// <summary>
    /// Called by wxFontEnumerator::EnumerateFacenames() for each match.
    /// </summary>
    /// <param name="facename">name identificator to load face by wxFont</param>
    /// <returns> True to continue enumeration or false to stop it.</returns>
    virtual bool OnFacename(const wxString& facename) wxOVERRIDE {
        // vertical font start with @, we will filter it out
        if (facename.empty() || facename[0] == '@') return true;
        wxFont wx_font(wxFontInfo().FaceName(facename).Encoding(m_encoding));

        //*
        // Faster chech if wx_font is loadable but not 100%
        // names could contain not loadable font
        if (!WxFontUtils::can_load(wx_font)) {
#ifdef DEBUG_NOT_LOADABLE_FONTS
            m_efacenames.emplace_back(facename.c_str());
#endif // DEBUG_NOT_LOADABLE_FONTS
            return true; // can't load
        }
        /*/
        // Slow copy of font files to try load font
        // After this all files are loadable
        auto font_file = WxFontUtils::create_font_file(wx_font);
        if (font_file == nullptr) {
#ifdef DEBUG_NOT_LOADABLE_FONTS
            m_efacenames.emplace_back(facename.c_str());
#endif // DEBUG_NOT_LOADABLE_FONTS
            return true; // can't create font file
        } // */
        m_facenames.push_back(facename);
        return true;
    }
};

bool GLGizmoEmboss::select_facename(const wxString &facename) {
    if (!wxFontEnumerator::IsValidFacename(facename)) return false;
    // Select font
    const wxFontEncoding &encoding = m_face_names.encoding;
    wxFont wx_font(wxFontInfo().FaceName(facename).Encoding(encoding));
    if (!wx_font.IsOk()) return false;
    if (!m_font_manager.set_wx_font(wx_font)) return false;
    process();
    return true;
}

void GLGizmoEmboss::init_face_names() {
    if (m_face_names.is_init) return;
    m_face_names.is_init      = true;
    wxFontEncoding   encoding = wxFontEncoding::wxFONTENCODING_SYSTEM;
    MyFontEnumerator fontEnumerator(encoding);
    bool             fixed_width_only = false;
    fontEnumerator.EnumerateFacenames(encoding, fixed_width_only);
    m_face_names.encoding = encoding;
    m_face_names.names    = std::move(fontEnumerator.m_facenames);
    std::sort(m_face_names.names.begin(), m_face_names.names.end());
}

#include "slic3r/GUI/Jobs/CreateFontNameImageJob.hpp"
void GLGizmoEmboss::draw_font_list()
{
    // Create of texture for font name
    auto init_texture = [&face_names = m_face_names, size = m_gui_cfg->face_name_size]() {
        // check if already exists
        GLuint &id = face_names.texture_id; 
        if (id != 0) return;
        // create texture for font
        GLenum target = GL_TEXTURE_2D, format = GL_ALPHA,
               type = GL_UNSIGNED_BYTE;
        GLint level = 0, border = 0;
        glsafe(::glGenTextures(1, &id));
        glsafe(::glBindTexture(target, id));
        glsafe(::glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
        glsafe(::glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
        GLint w = size.x(), h = face_names.names.size() * size.y();
        std::vector<unsigned char> data(w * h, {0});
        glsafe(::glTexImage2D(target, level, GL_ALPHA, w, h, border, format, type, data.data()));
        // bind default texture
        GLuint no_texture_id = 0;
        glsafe(::glBindTexture(target, no_texture_id));

        // no one is initialized yet
        face_names.exist_textures = std::vector<bool>(face_names.names.size(), {false});
    };
    auto init_trancated_names = [&face_names = m_face_names, 
        width = m_gui_cfg->face_name_max_width]() { 
        if (!face_names.names_truncated.empty()) return;
        face_names.names_truncated.reserve(face_names.names.size());
        for (const wxString &face_name : face_names.names) {
            std::string name(face_name.ToUTF8().data());
            face_names.names_truncated.emplace_back(ImGuiWrapper::trunc(name, width));
        }
    };

    // Set partial
    wxString actual_face_name;
    if (m_font_manager.is_activ_font()) {
        std::optional<wxFont> &wx_font_opt = m_font_manager.get_wx_font();
        if (wx_font_opt.has_value())
            actual_face_name = wx_font_opt->GetFaceName();
    }
    const char * selected = (!actual_face_name.empty()) ?
        actual_face_name.ToUTF8().data() : " --- ";
    wxString del_facename;
    if (ImGui::BeginCombo("##font_selector", selected)) {
        if (!m_face_names.is_init) init_face_names();
        init_texture();
        if (m_face_names.names_truncated.empty()) init_trancated_names();
        ImTextureID tex_id = (void *) (intptr_t) m_face_names.texture_id;
        for (const wxString &face_name : m_face_names.names) {
            size_t index = &face_name - &m_face_names.names.front();
            ImGui::PushID(index);
            bool is_selected = (actual_face_name == face_name);
            if (ImGui::Selectable(m_face_names.names_truncated[index].c_str(), is_selected)) {
                if (!select_facename(face_name)) {
                    del_facename = face_name;
                    wxMessageBox(GUI::format(
                        _L("Font face \"%1%\" can't be selected."), face_name));
                }
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", face_name.ToUTF8().data());
            if (is_selected) ImGui::SetItemDefaultFocus();
            if (!m_face_names.exist_textures[index] && 
                ImGui::IsItemVisible()) {
                m_face_names.exist_textures[index] = true;
                std::string text = m_text.empty() ? "AaBbCc" : m_text;
                // render text to texture
                FontImageData data{text,
                                   face_name,
                                   m_face_names.encoding,
                                   m_face_names.texture_id,
                                   index,
                                   m_gui_cfg->face_name_size,
                                   &m_allow_update_rendered_font
                };
                auto job = std::make_unique<CreateFontImageJob>(std::move(data));
                auto& worker = wxGetApp().plater()->get_ui_job_worker();
                queue_job(worker, std::move(job));
            }
            ImGui::SameLine(m_gui_cfg->face_name_offset);
            ImVec2 size(m_gui_cfg->face_name_size.x(),
                        m_gui_cfg->face_name_size.y()),
                uv0(0.f, index / (float) m_face_names.names.size()),
                uv1(1.f, (index + 1) / (float) m_face_names.names.size());
            ImGui::Image(tex_id, size, uv0, uv1);
            ImGui::PopID();
        }        
#ifdef SHOW_FONT_COUNT
        ImGui::TextColored(ImGuiWrapper::COL_GREY_LIGHT, "Count %d",
                           static_cast<int>(m_face_names.names.size()));
#endif // SHOW_FONT_COUNT
        ImGui::EndCombo();
        m_allow_update_rendered_font = true;
    } else if (m_allow_update_rendered_font) { 
        // free texture and set id to zero
        m_allow_update_rendered_font = false;
        glsafe(::glDeleteTextures(1, &m_face_names.texture_id));
        m_face_names.texture_id = 0;
    }

    // delete unloadable face name when appear
    if (!del_facename.empty()) {
        // IMPROVE: store list of deleted facename into app.ini
        std::vector<wxString> &f = m_face_names.names;
        f.erase(std::remove(f.begin(), f.end(), del_facename), f.end());
        m_face_names.names_truncated.clear();
    }

#ifdef ALLOW_ADD_FONT_BY_FILE
    ImGui::SameLine();
    // select font file by file browser
    if (draw_button(IconType::open_file)) {
        if (choose_true_type_file()) { 
            process();
        }
    } else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("add file with font(.ttf, .ttc)").c_str());
#endif //  ALLOW_ADD_FONT_BY_FILE

#ifdef ALLOW_ADD_FONT_BY_OS_SELECTOR
    ImGui::SameLine();
    if (draw_button(IconType::system_selector)) {
        if (choose_font_by_wxdialog()) {
            process();
        }
    } else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Open dialog for choose from fonts.").c_str());
#endif //  ALLOW_ADD_FONT_BY_OS_SELECTOR

}

void GLGizmoEmboss::draw_model_type()
{
    std::optional<ModelVolumeType> new_type;
    ModelVolumeType modifier = ModelVolumeType::PARAMETER_MODIFIER;
    ModelVolumeType negative = ModelVolumeType::NEGATIVE_VOLUME;
    ModelVolumeType part = ModelVolumeType::MODEL_PART;
    ModelVolumeType type = (m_volume != nullptr) ? m_volume->type() :
                                                   ModelVolumeType::INVALID;

    bool is_last_solid_part = is_text_object(m_volume);
    ImGui::Text("%s", m_gui_cfg->translations.type.c_str());
    ImGui::SameLine(m_gui_cfg->style_offset);
    if (type == part) { 
        draw_icon(IconType::part, IconState::hovered);
    } else {
        if (draw_button(IconType::part)) new_type = part;
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Click to change text into object part.").c_str());
    }

    ImGui::SameLine();
    if (type == negative) { 
        draw_icon(IconType::negative, IconState::hovered);    
    } else {
        if (draw_button(IconType::negative, is_last_solid_part))
            new_type = negative;        
        if(ImGui::IsItemHovered()){
            if(is_last_solid_part)
                ImGui::SetTooltip("%s", _u8L("You can't change a type of the last solid part of the object.").c_str());
            else if (type != negative)
                ImGui::SetTooltip("%s", _u8L("Click to change part type into negative volume.").c_str());
        }
    }

    ImGui::SameLine();
    if (type == modifier) {
        draw_icon(IconType::modifier, IconState::hovered);  
    } else {
        if(draw_button(IconType::modifier, is_last_solid_part))
            new_type = modifier;    
        if (ImGui::IsItemHovered()) {
            if(is_last_solid_part)
                ImGui::SetTooltip("%s", _u8L("You can't change a type of the last solid part of the object.").c_str());
            else if (type != modifier)
                ImGui::SetTooltip("%s", _u8L("Click to change part type into modifier.").c_str());
        }
    }

    if (m_volume != nullptr && new_type.has_value() && !is_last_solid_part) {
        GUI_App &app    = wxGetApp();
        Plater * plater = app.plater();
        Plater::TakeSnapshot snapshot(plater, _L("Change Text Type"), UndoRedo::SnapshotType::GizmoAction);
        m_volume->set_type(*new_type);

        // inspiration in ObjectList::change_part_type()
        // how to view correct side panel with objects
        ObjectList *obj_list = app.obj_list();
        wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(
            obj_list->get_selected_obj_idx(),
            [volume = m_volume](const ModelVolume *vol) { return vol == volume; });
        if (!sel.IsEmpty()) obj_list->select_item(sel.front());

        // Update volume position when switch from part or into part
        if (m_volume->text_configuration->font_item.prop.use_surface) {
            // move inside
            bool is_volume_move_inside = (type == part);
            bool is_volume_move_outside = (*new_type == part);
            if (is_volume_move_inside || is_volume_move_outside)
                process();
        }
    }
}

void GLGizmoEmboss::draw_rename_style_botton()
{
    bool start_rename = false;
    if (draw_button(IconType::rename))
        start_rename = true;
    else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Rename actual style.").c_str());

    std::string title = _u8L("Rename style");
    const char * popup_id = title.c_str();
    static FontItem *  rename_item;
    static std::string new_name;
    if (start_rename && !ImGui::IsPopupOpen(popup_id)) {
        ImGui::OpenPopup(popup_id);
        rename_item = &m_font_manager.get_font_item();
        new_name    = rename_item->name; // initialize with original copy
    }

    if (ImGui::BeginPopupModal(popup_id, 0, ImGuiWindowFlags_AlwaysAutoResize)) {
        const std::string &original_style_name = rename_item->name;
        std::string text_in_popup =
            GUI::format(_u8L("Rename style(%1%) for embossing text: "), original_style_name);
        ImGui::Text("%s", text_in_popup.c_str());

        bool is_unique = true;
        for (const auto &item : m_font_manager.get_styles()) {
            const FontItem &fi = item.font_item;
            if (&fi == rename_item)
                continue; // could be same as original name
            if (fi.name == new_name) is_unique = false;
        }
        bool allow_change = false;
        if (new_name.empty()) {
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_DARK, _u8L("Name can't be empty."));
        }else if (!is_unique) { 
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_DARK, _u8L("Name has to be unique."));
        } else {
            ImGui::NewLine();
            allow_change = true;
        }

        ImGuiInputTextFlags flags = ImGuiInputTextFlags_EnterReturnsTrue;
        if ((ImGui::InputText("##rename style", &new_name, flags) && allow_change) ||
            m_imgui->button(_L("ok"), ImVec2(0.f, 0.f), allow_change)) {
            rename_item->name = new_name;
            m_font_manager.get_truncated_name().clear();
            ImGui::CloseCurrentPopup();
            select_stored_font_item();
        }
        ImGui::EndPopup();
    }
}

void GLGizmoEmboss::draw_delete_style_button() {
    const std::vector<FontManager::Item> &styles = m_font_manager.get_styles();    
    bool can_delete = styles.size() > 1;
    std::string title = _u8L("Remove style");
    const char * popup_id = title.c_str();
    static size_t next_style_index = std::numeric_limits<size_t>::max();
    if (draw_button(IconType::erase, !can_delete)) {
        while (true) {
            const FontManager::Item &style = m_font_manager.get_activ_style();
            size_t activ_index = &style - &styles.front();
            next_style_index = (activ_index > 0) ? activ_index - 1 :
                                                   activ_index + 1;
            if (next_style_index >= styles.size()) {
                // can't remove last font style
                break;
            }
            // clean unactivable styles
            if (!m_font_manager.load_font(next_style_index)) {
                m_font_manager.erase(next_style_index);
                continue;
            }

            // load back
            m_font_manager.load_font(activ_index);
            ImGui::OpenPopup(popup_id);
            break;
        }
    }
    if (ImGui::IsItemHovered()) {
        const FontManager::Item &style = m_font_manager.get_activ_style();
        std::string style_name = style.font_item.name;
        std::string tooltip = can_delete ?            
            GUI::format(_L("Delete \"%1%\" style."), style_name) :
            GUI::format(_L("Can't delete \"%1%\". It is last style."), style_name) ;
        ImGui::SetTooltip("%s", tooltip.c_str());
    }

    if (ImGui::BeginPopupModal(popup_id)) {
        const FontManager::Item &style = m_font_manager.get_activ_style();
        std::string style_name = style.font_item.name;
        std::string text_in_popup = GUI::format(_u8L("Are you sure,\nthat you want permanently and unrecoverable \nremove style \"%1%\"?"), style_name);
        ImGui::Text("%s", text_in_popup.c_str());
        if (ImGui::Button(_u8L("Yes").c_str())) {
            m_font_manager.load_font(next_style_index);
            size_t activ_index = &style - &styles.front();
            m_font_manager.erase(activ_index);
            store_font_list_to_app_config();
            ImGui::CloseCurrentPopup();
            process();
        }
        ImGui::SameLine();
        if (ImGui::Button(_u8L("No").c_str()))
            ImGui::CloseCurrentPopup();        
        ImGui::EndPopup();
    }
}

void GLGizmoEmboss::draw_undo_style_button(bool is_stored, bool is_changed)
{
    bool can_undo = is_stored && is_changed;
    if (draw_button(IconType::undo, !can_undo)) {
        FontItem &font_item  = m_font_manager.get_font_item();
        bool is_path_changed = font_item.path != m_stored_font_item->path;
            
        // is rotation changed
        auto &angle = font_item.prop.angle;
        const auto &angle_ = m_stored_font_item->prop.angle;
        // TODO: compare with approx
        if (angle.has_value() != angle_.has_value() ||
            (angle.has_value() && !is_approx(*angle, *angle_))) {
            auto &tc = m_volume->text_configuration;
            if (m_volume != nullptr && tc.has_value()) {
                // change actual text configuration
                tc->font_item.prop.angle = angle_;
                float act_angle  = angle_.has_value() ? *angle_ : .0f;
                float prev_angle = angle.has_value() ? *angle : .0f;
                do_rotate(act_angle - prev_angle);
            }
        }
            
        // is distance changed
        auto &distance = font_item.prop.distance;
        const auto &distance_ = m_stored_font_item->prop.distance;
        if (distance.has_value() != distance_.has_value() ||
            (distance.has_value() && !is_approx(*distance, *distance_))) {
            auto &tc = m_volume->text_configuration;
            if (m_volume != nullptr && tc.has_value()) {
                tc->font_item.prop.distance = distance_;
                float act_distance = distance_.has_value() ? *distance_ : .0f;
                float prev_distance = distance.has_value() ? *distance : .0f;
                do_translate(Vec3d::UnitZ() * (act_distance - prev_distance));
            }
        }

        font_item = *m_stored_font_item;
        if (is_path_changed) {
            m_font_manager.get_wx_font() = WxFontUtils::load_wxFont(font_item.path);
            m_font_manager.wx_font_changed();
        }
        process();
    } else if (ImGui::IsItemHovered()) {
        if (can_undo)
            ImGui::SetTooltip("%s", _u8L("Reload stored values of selected style").c_str());
        else if (!is_stored)
            ImGui::SetTooltip("%s", _u8L("Nothing to reload from").c_str());
        else if (!is_changed)
            ImGui::SetTooltip("%s", _u8L("No change to restore from style").c_str());
    }
}

void GLGizmoEmboss::draw_revert_all_styles_button() {
    if (draw_button(IconType::revert_all)) {
        m_font_manager     = FontManager(m_imgui->get_glyph_ranges());
        FontList font_list = create_default_font_list();
        m_font_manager.add_fonts(font_list);
        // What to do when NO one default font is loadable?
        [[maybe_unused]] bool success = m_font_manager.load_first_valid_font();
        assert(success);
        select_stored_font_item();
        process();
    } else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Revert all styles").c_str());
}

void GLGizmoEmboss::draw_style_list() {
    if (!m_font_manager.is_activ_font()) return;
    const float &max_width = m_gui_cfg->max_font_name_width;
    const FontItem &actual_font_item = m_font_manager.get_font_item();
    std::string &trunc_name = m_font_manager.get_truncated_name();
    if (trunc_name.empty()) {
        // generate trunc name
        const std::string &current_name  = actual_font_item.name;
        trunc_name = ImGuiWrapper::trunc(current_name, max_width);
    }
    ImGui::Text("%s", m_gui_cfg->translations.style.c_str());
    ImGui::SameLine(m_gui_cfg->style_offset);
    ImGui::SetNextItemWidth(m_gui_cfg->input_width);
    if (ImGui::BeginCombo("##style_selector", trunc_name.c_str())) {
        m_font_manager.init_style_images(m_gui_cfg->max_style_image_size, m_text);
        m_font_manager.init_trunc_names(max_width);
        std::optional<std::pair<size_t,size_t>> swap_indexes;
        const std::vector<FontManager::Item> &styles = m_font_manager.get_styles();
        for (const auto &item : styles) {
            size_t             index             = &item - &styles.front();
            const FontItem &   fi                = item.font_item;
            const std::string &actual_style_name = fi.name;
            ImGui::PushID(actual_style_name.c_str());
            bool is_selected = (&fi == &actual_font_item);

            ImVec2 select_size(0,m_gui_cfg->max_style_image_size.y()); // 0,0 --> calculate in draw
            const std::optional<FontManager::StyleImage> img = item.image;            
            // allow click delete button
            ImGuiSelectableFlags_ flags = ImGuiSelectableFlags_AllowItemOverlap; 
            if (ImGui::Selectable(item.truncated_name.c_str(), is_selected, flags, select_size)) {
                if (m_font_manager.load_font(index)) { 
                    process();
                    select_stored_font_item();
                }
            } else if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", actual_style_name.c_str());

            // reorder items
            if (ImGui::IsItemActive() && !ImGui::IsItemHovered()) {
                if (ImGui::GetMouseDragDelta(0).y < 0.f) {
                    if (index > 0) 
                        swap_indexes = {index, index - 1};
                } else if ((index + 1) < styles.size())
                    swap_indexes = {index, index + 1};
                if (swap_indexes.has_value()) 
                    ImGui::ResetMouseDragDelta();
            }

            // draw style name
            if (img.has_value()) {
                ImGui::SameLine();
                ImGui::Image(img->texture_id, img->tex_size, img->uv0, img->uv1);
            }

            ImGui::PopID();
        }
        if (swap_indexes.has_value()) 
            m_font_manager.swap(swap_indexes->first,
                                swap_indexes->second);
        ImGui::EndCombo();
    } else {
        // do not keep in memory style images when no combo box open
        m_font_manager.free_style_images();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Style selector").c_str());
    }
        
    ImGui::SameLine();
    draw_rename_style_botton();

    ImGui::SameLine();
    if (draw_button(IconType::duplicate)) m_font_manager.duplicate();
    else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Duplicate style.").c_str());
    
    // Is style changed against stored one
    FontItem &font_item  = m_font_manager.get_font_item();
    bool is_stored  = m_stored_font_item.has_value();
    bool is_changed = (is_stored) ? !(*m_stored_font_item == font_item) : true;
    // TODO: check order of font items in list to allowe save actual order

    ImGui::SameLine();
    if (draw_button(IconType::save, !is_changed)) {
        // save styles to app config
        store_font_list_to_app_config();
    }else if (ImGui::IsItemHovered()) { 
        if (is_changed) {
            ImGui::SetTooltip("%s", _u8L("Save current settings to selected style").c_str());
        } else {
            ImGui::SetTooltip("%s", _u8L("No changes to save into style").c_str());
        }
    }

    ImGui::SameLine();
    draw_undo_style_button(is_stored, is_changed);    
    
#ifdef ALLOW_REVERT_ALL_STYLES
    ImGui::SameLine();
    draw_revert_all_styles_button();
#endif // ALLOW_REVERT_ALL_STYLES

    // delete button
    ImGui::SameLine();
    draw_delete_style_button();
}

bool GLGizmoEmboss::italic_button()
{
    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font(); 
    const std::shared_ptr<const Emboss::FontFile> &font_file = m_font_manager.get_font_file();
    if (!wx_font.has_value() || font_file == nullptr) { 
        draw_icon(IconType::italic, IconState::disabled);
        return false;
    }

    std::optional<float> &skew = m_font_manager.get_font_prop().skew;
    bool is_font_italic = skew.has_value() || WxFontUtils::is_italic(*wx_font);
    if (is_font_italic) {
        // unset italic
        if (draw_clickable(IconType::italic, IconState::hovered,
                           IconType::unitalic, IconState::hovered)) {
            skew.reset();
            if (wx_font->GetStyle() != wxFontStyle::wxFONTSTYLE_NORMAL) {
                wx_font->SetStyle(wxFontStyle::wxFONTSTYLE_NORMAL);
                m_font_manager.wx_font_changed();
            }
            return true;
        }
        if (ImGui::IsItemHovered()) 
            ImGui::SetTooltip("%s", _u8L("Unset italic").c_str());
    } else {
        // set italic
        if (draw_button(IconType::italic)) {
            auto new_ff = WxFontUtils::set_italic(*wx_font, *font_file);
            if (new_ff != nullptr) {
                m_font_manager.wx_font_changed(std::move(new_ff));
            } else {
                // italic font doesn't exist 
                // add skew when wxFont can't set it
                skew = 0.2f;
            }            
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Set italic").c_str());
    }
    return false;
}

bool GLGizmoEmboss::bold_button() {
    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();
    const std::shared_ptr<const Emboss::FontFile> &font_file = m_font_manager.get_font_file();
    if (!wx_font.has_value() || font_file==nullptr) {
        draw_icon(IconType::bold, IconState::disabled);
        return false;
    }
    
    std::optional<float> &boldness = m_font_manager.get_font_prop().boldness;
    bool is_font_bold = boldness.has_value() || WxFontUtils::is_bold(*wx_font);
    if (is_font_bold) {
        // unset bold
        if (draw_clickable(IconType::bold, IconState::hovered,
                           IconType::unbold, IconState::hovered)) {
            boldness.reset();
            if (wx_font->GetWeight() != wxFontWeight::wxFONTWEIGHT_NORMAL) {
                wx_font->SetWeight(wxFontWeight::wxFONTWEIGHT_NORMAL);
                m_font_manager.wx_font_changed();
            }
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Unset bold").c_str());
    } else {
        // set bold
        if (draw_button(IconType::bold)) {
            auto new_ff = WxFontUtils::set_bold(*wx_font, *font_file);
            if (new_ff != nullptr) {
                m_font_manager.wx_font_changed(std::move(new_ff));
            } else {
                // bold font can't be loaded
                // set up boldness
                boldness = 20.f;
                //font_file->cache.empty();
            }
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Set bold").c_str());
    }
    return false;
}

template<typename T, typename Draw>
bool GLGizmoEmboss::revertible(const std::string &name,
                               T                 &value,
                               T                 *default_value,
                               bool               exist_change,
                               const std::string &undo_tooltip,
                               float              undo_offset,
                               Draw               draw)
{
    if (exist_change)
        ImGuiWrapper::text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, name);
    else
        ImGuiWrapper::text(name);
    // render revert changes button
    if (exist_change) {        
        ImGui::SameLine(undo_offset);
        if (draw_button(IconType::undo)) {
            value = *default_value;
            return true;
        } else if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", undo_tooltip.c_str());
    }
    return draw();
}


bool GLGizmoEmboss::rev_input(const std::string  &name,
                              float              &value,
                              float              *default_value,
                              const std::string  &undo_tooltip,
                              float               step,
                              float               step_fast,
                              const char         *format,
                              ImGuiInputTextFlags flags)
{
    // draw offseted input
    auto draw_offseted_input = [&]()->bool{
        float input_offset = m_gui_cfg->input_offset;
        float input_width  = m_gui_cfg->input_width;
        ImGui::SameLine(input_offset);
        ImGui::SetNextItemWidth(input_width);
        return ImGui::InputFloat(("##" + name).c_str(),
            &value, step, step_fast, format, flags);
    };
    float undo_offset = ImGui::GetStyle().FramePadding.x;
    bool exist_change = default_value != nullptr ? (!is_approx(value, *default_value)) : false;
    return revertible(name, value, default_value, exist_change, undo_tooltip, undo_offset, draw_offseted_input);
}

bool GLGizmoEmboss::rev_checkbox(const std::string &name,
                                 bool              &value,
                                 bool              *default_value,
                                 const std::string &undo_tooltip)
{
    // draw offseted input
    auto draw_offseted_input = [&]() -> bool {
        ImGui::SameLine(m_gui_cfg->advanced_input_offset);
        return ImGui::Checkbox(("##" + name).c_str(), &value);
    };
    float undo_offset  = ImGui::GetStyle().FramePadding.x;
    bool  exist_change = default_value != nullptr ?
                             (value != *default_value) :
                             false;
    return revertible(name, value, default_value, exist_change, undo_tooltip,
                      undo_offset, draw_offseted_input);
}

void GLGizmoEmboss::draw_style_edit() {
    const GuiCfg::Translations &tr = m_gui_cfg->translations;

    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();
    FontItem &fi = m_font_manager.get_font_item();
    // TODO: should not be there
    // when actual font not loaded try to load
    if (!wx_font.has_value() && fi.type == WxFontUtils::get_actual_type())
        wx_font = WxFontUtils::load_wxFont(fi.path);

    bool is_font_changed = false;
    if (m_stored_font_item.has_value() && wx_font.has_value() && m_stored_wx_font.has_value()) {
        bool is_font_face_changed = m_stored_wx_font->GetFaceName() != wx_font->GetFaceName();

        const std::optional<float> &skew = m_font_manager.get_font_prop().skew;
        bool is_italic = skew.has_value() || WxFontUtils::is_italic(*wx_font);
        const std::optional<float> &skew_stored = m_stored_font_item->prop.skew;
        bool is_stored_italic = skew_stored.has_value() || WxFontUtils::is_italic(*m_stored_wx_font);
        bool is_italic_changed = is_italic != is_stored_italic;

        const std::optional<float> &boldness = m_font_manager.get_font_prop().boldness;
        bool is_bold = boldness.has_value() || WxFontUtils::is_bold(*wx_font);
        const std::optional<float> &boldness_stored = m_stored_font_item->prop.boldness;
        bool is_stored_bold = boldness_stored.has_value() || WxFontUtils::is_bold(*m_stored_wx_font);
        bool is_bold_changed = is_bold != is_stored_bold;

        bool is_font_style_changed = is_italic_changed || is_bold_changed;

        is_font_changed = is_font_face_changed || is_font_style_changed;
    }

    if (is_font_changed)
        ImGuiWrapper::text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, tr.font);
    else
        ImGuiWrapper::text(tr.font);
    ImGui::SameLine(m_gui_cfg->input_offset);
    ImGui::SetNextItemWidth(m_gui_cfg->input_width);
    draw_font_list();
    ImGui::SameLine();
    bool exist_change = false;
    exist_change |= italic_button();
    ImGui::SameLine();
    exist_change |= bold_button();
    
    if (is_font_changed) {
        ImGui::SameLine(ImGui::GetStyle().FramePadding.x);
        if (draw_button(IconType::undo)) {
            fi.path = m_stored_font_item->path;
            fi.prop.boldness = m_stored_font_item->prop.boldness;
            fi.prop.skew = m_stored_font_item->prop.skew;

            wx_font = WxFontUtils::load_wxFont(fi.path);
            m_font_manager.wx_font_changed();
            exist_change = true;
        } else if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Revert font changes.").c_str());
    }

    if (exist_change) {
        m_font_manager.clear_glyphs_cache();
        process();
    }

    FontProp &font_prop = fi.prop;
    float * def_size = m_stored_font_item.has_value() ? 
        &m_stored_font_item->prop.size_in_mm : nullptr;
    if (rev_input(tr.size, font_prop.size_in_mm, def_size,
                  _u8L("Revert text size."), 0.1f, 1.f, "%.1f mm")) {
        // size can't be zero or negative
        Limits::apply(font_prop.size_in_mm, limits.size_in_mm);

        // only different value need process
        if (m_volume != nullptr &&
            m_volume->text_configuration.has_value() &&
            !is_approx(font_prop.size_in_mm, m_volume->text_configuration->font_item.prop.size_in_mm)) {

            // store font size into path
            if (fi.type == WxFontUtils::get_actual_type()) {
                if (wx_font.has_value()) {
                    wx_font->SetPointSize(
                        static_cast<int>(font_prop.size_in_mm));
                    m_font_manager.wx_font_changed();
                }
            }
            process();
        }
    }

#ifdef SHOW_WX_WEIGHT_INPUT
    if (wx_font.has_value()) {
        ImGui::Text("%s", "weight");
        ImGui::SameLine(m_gui_cfg->input_offset);
        ImGui::SetNextItemWidth(m_gui_cfg->input_width);
        int weight     = wx_font->GetNumericWeight();
        int min_weight = 1, max_weight = 1000;
        if (ImGui::SliderInt("##weight", &weight, min_weight, max_weight)) {
            wx_font->SetNumericWeight(weight);
            m_font_manager.wx_font_changed();
            process();
        }

        wxFont f       = wx_font->Bold();
        bool   disable = f == *wx_font;
        ImGui::SameLine();
        if (draw_button(IconType::bold, disable)) {
            *wx_font = f;
            m_font_manager.wx_font_changed();
            process();
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("wx Make bold").c_str());
    }
#endif // SHOW_WX_WEIGHT_INPUT

    float *def_depth = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.emboss : nullptr;
    if (rev_input(tr.depth, font_prop.emboss, def_depth,
        _u8L("Revert embossed depth."), 0.1f, 0.25, "%.2f mm")) {
        Limits::apply(font_prop.emboss, limits.emboss);
        process();
    }    
}

bool GLGizmoEmboss::rev_slider(const std::string &name,
                               std::optional<int>& value,
                               std::optional<int> *default_value,
                               const std::string &undo_tooltip,
                               int                v_min,
                               int                v_max,
                               const std::string& format,
                               const wxString    &tooltip)
{    
    auto draw_slider_optional_int = [&]() -> bool {
        float slider_offset = m_gui_cfg->advanced_input_offset;
        float slider_width  = m_gui_cfg->input_width;
        ImGui::SameLine(slider_offset);
        ImGui::SetNextItemWidth(slider_width);
        return m_imgui->slider_optional_int( ("##" + name).c_str(), value, 
            v_min, v_max, format.c_str(), 1.f, false, tooltip);
    };
    float undo_offset = ImGui::GetStyle().FramePadding.x;    
    bool exist_change = default_value != nullptr ? (value != *default_value) : false;
    return revertible(name, value, default_value, exist_change, 
        undo_tooltip, undo_offset, draw_slider_optional_int);
}

bool GLGizmoEmboss::rev_slider(const std::string &name,
                               std::optional<float>& value,
                               std::optional<float> *default_value,
                               const std::string &undo_tooltip,
                               float                v_min,
                               float                v_max,
                               const std::string& format,
                               const wxString    &tooltip)
{    
    auto draw_slider_optional_float = [&]() -> bool {
        float slider_offset = m_gui_cfg->advanced_input_offset;
        float slider_width  = m_gui_cfg->input_width;
        ImGui::SameLine(slider_offset);
        ImGui::SetNextItemWidth(slider_width);
        return m_imgui->slider_optional_float(("##" + name).c_str(), value,
            v_min, v_max, format.c_str(), 1.f, false, tooltip);
    };
    float undo_offset = ImGui::GetStyle().FramePadding.x;    
    bool exist_change = (default_value != nullptr)? 
        (!is_approx(value, *default_value)) : false;
    return revertible(name, value, default_value, exist_change,
        undo_tooltip, undo_offset, draw_slider_optional_float);
}

bool GLGizmoEmboss::rev_slider(const std::string &name,
                               float             &value,
                               float             *default_value,
                               const std::string &undo_tooltip,
                               float              v_min,
                               float              v_max,
                               const std::string &format,
                               const wxString    &tooltip)
{    
    auto draw_slider_float = [&]() -> bool {
        float slider_offset = m_gui_cfg->advanced_input_offset;
        float slider_width  = m_gui_cfg->input_width;
        ImGui::SameLine(slider_offset);
        ImGui::SetNextItemWidth(slider_width);
        return m_imgui->slider_float("##" + name, &value, v_min, v_max,
            format.c_str(), 1.f, false, tooltip);
    };
    float undo_offset = ImGui::GetStyle().FramePadding.x;    
    bool exist_change = default_value != nullptr ? 
        (!is_approx(value, *default_value)) : false;
    return revertible(name, value, default_value, exist_change,
        undo_tooltip, undo_offset, draw_slider_float);
}

void GLGizmoEmboss::do_translate(const Vec3d &relative_move)
{
    assert(m_volume != nullptr);
    assert(m_volume->text_configuration.has_value());
    Selection &selection = m_parent.get_selection();
    assert(!selection.is_empty());
    selection.setup_cache();
    selection.translate(relative_move, TransformationType::Local);

    std::string snapshot_name; // empty meand no store undo / redo
    // NOTE: it use L instead of _L macro because prefix _ is appended inside
    // function do_move
    // snapshot_name = L("Set surface distance");
    m_parent.do_move(snapshot_name);
}

void GLGizmoEmboss::do_rotate(float relative_z_angle)
{
    assert(m_volume != nullptr);
    assert(m_volume->text_configuration.has_value());
    Selection &selection = m_parent.get_selection();
    assert(!selection.is_empty());
    selection.setup_cache();
    selection.rotate(Vec3d(0., 0., relative_z_angle), TransformationType::Local);

    std::string snapshot_name; // empty meand no store undo / redo
    // NOTE: it use L instead of _L macro because prefix _ is appended
    // inside function do_move
    // snapshot_name = L("Set text rotation");
    m_parent.do_rotate(snapshot_name);
}

void GLGizmoEmboss::draw_advanced()
{
    const std::shared_ptr<const Emboss::FontFile>& font_file = m_font_manager.get_font_file();
    if (font_file == nullptr) { 
        ImGui::Text("%s", _u8L("Advanced font options could be change only for corect font.\nStart with select correct font.").c_str());
        return;
    }

    FontProp &font_prop = m_font_manager.get_font_item().prop;
    const auto  &cn = m_font_manager.get_font_prop().collection_number;
    unsigned int font_index = (cn.has_value()) ? *cn : 0;
    const auto  &font_info  = font_file->infos[font_index];

#ifdef SHOW_FONT_FILE_PROPERTY
    ImGui::SameLine();
    auto& ff = m_font_manager.get_activ_style().font_file_with_cache;
    int cache_size = ff.has_value()? (int)ff.cache->size() : 0;
    std::string ff_property = 
        "ascent=" + std::to_string(font_info.ascent) +
        ", descent=" + std::to_string(font_info.descent) +
        ", lineGap=" + std::to_string(font_info.linegap) +
        ", unitPerEm=" + std::to_string(font_info.unit_per_em) + 
        ", cache(" + std::to_string(cache_size) + " glyphs)";
    if (font_file->infos.size() > 1) { 
        unsigned int collection = font_prop.collection_number.has_value() ?
            *font_prop.collection_number : 0;
        ff_property += ", collect=" + std::to_string(collection+1) + "/" + std::to_string(font_file->infos.size());
    }
    m_imgui->text_colored(ImGuiWrapper::COL_GREY_DARK, ff_property);
#endif // SHOW_FONT_FILE_PROPERTY

    bool exist_change = false;
    auto &tr = m_gui_cfg->translations;

    bool *def_use_surface = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.use_surface : nullptr;
    if (rev_checkbox(tr.use_surface, font_prop.use_surface, def_use_surface,
                     _u8L("Revert using of model surface."))) {
        if (font_prop.use_surface) { 
            font_prop.distance.reset();
            if (font_prop.emboss < 0.1)
                font_prop.emboss = 1;
        }
        process();
    }

    std::string units = _u8L("font points");
    std::string units_fmt = "%.0f " + units;
    
    // input gap between letters
    auto def_char_gap = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.char_gap : nullptr;

    int half_ascent = font_info.ascent / 2;
    int min_char_gap = -half_ascent, max_char_gap = half_ascent;
    if (rev_slider(tr.char_gap, font_prop.char_gap, def_char_gap, _u8L("Revert gap between letters"), 
        min_char_gap, max_char_gap, units_fmt, _L("Distance between letters"))){
        Limits::apply(font_prop.char_gap, limits.char_gap);
        // char gap is stored inside of imgui font atlas
        m_font_manager.clear_imgui_font();
        exist_change = true;
    }

    // input gap between lines
    auto def_line_gap = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.line_gap : nullptr;
    int  min_line_gap = -half_ascent, max_line_gap = half_ascent;
    if (rev_slider(tr.line_gap, font_prop.line_gap, def_line_gap, _u8L("Revert gap between lines"), 
        min_line_gap, max_line_gap, units_fmt, _L("Distance between lines"))){
        Limits::apply(font_prop.line_gap, limits.line_gap);
        // line gap is planed to be stored inside of imgui font atlas
        m_font_manager.clear_imgui_font();
        exist_change = true;
    }

    // input boldness
    auto def_boldness = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.boldness : nullptr;
    if (rev_slider(tr.boldness, font_prop.boldness, def_boldness, _u8L("Undo boldness"), 
        limits.boldness.gui.min, limits.boldness.gui.max, units_fmt, _L("Tiny / Wide glyphs"))){
        Limits::apply(font_prop.boldness, limits.boldness.values);
        exist_change = true;
    }

    // input italic
    auto def_skew = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.skew : nullptr;
    if (rev_slider(tr.italic, font_prop.skew, def_skew, _u8L("Undo letter's skew"),
        limits.skew.gui.min, limits.skew.gui.max, "%.2f", _L("Italic strength ratio"))){
        Limits::apply(font_prop.skew, limits.skew.values);
        exist_change = true;
    }
    
    // input surface distance
    bool allowe_surface_distance = 
        m_volume != nullptr &&
        m_volume->text_configuration.has_value() &&
        !m_volume->text_configuration->font_item.prop.use_surface &&
        !is_text_object(m_volume);    
    std::optional<float> &distance = font_prop.distance;
    float prev_distance = distance.has_value() ? *distance : .0f,
          min_distance = -2 * font_prop.emboss,
          max_distance = 2 * font_prop.emboss;
    auto def_distance = m_stored_font_item.has_value() ?
        &m_stored_font_item->prop.distance : nullptr;    
    m_imgui->disabled_begin(!allowe_surface_distance);
    if (rev_slider(tr.surface_distance, distance, def_distance, _u8L("Undo translation"), 
        min_distance, max_distance, "%.2f mm", _L("Distance center of text from model surface")) &&
        m_volume != nullptr && m_volume->text_configuration.has_value()){
        m_volume->text_configuration->font_item.prop.distance = font_prop.distance;        
        float act_distance = font_prop.distance.has_value() ? *font_prop.distance : .0f;
        do_translate(Vec3d::UnitZ() * (act_distance - prev_distance));
    }
    m_imgui->disabled_end();

    // slider for Clock-wise angle in degress
    // stored angle is optional CCW and in radians
    std::optional<float> &angle = font_prop.angle;
    float prev_angle = angle.has_value() ? *angle : .0f;
    // Convert stored value to degress
    // minus create clock-wise roation from CCW
    float angle_deg = angle.has_value() ?
        static_cast<float>(-(*angle) * 180 / M_PI) : .0f;
    float def_angle_deg_val = 
        (!m_stored_font_item.has_value() || 
         !m_stored_font_item->prop.angle.has_value()) ?
        0.f : (*m_stored_font_item->prop.angle * -180 / M_PI);
    float* def_angle_deg = m_stored_font_item.has_value() ?
        &def_angle_deg_val : nullptr;
    if (rev_slider(tr.angle, angle_deg, def_angle_deg, _u8L("Undo rotation"), 
        limits.angle.min, limits.angle.max, u8"%.2f °",
                   _L("Rotate text Clock-wise."))) {
        // convert back to radians and CCW
        angle = -angle_deg * M_PI / 180.0;
        if (angle < M_PI || angle > M_PI) {
            int count = std::round(*angle / 2 * M_PI);
            angle     = *angle - count * (2 * M_PI);
        }
        if (is_approx(*angle, 0.f))
            angle.reset();
        
        if (m_volume != nullptr && m_volume->text_configuration.has_value()) {
            m_volume->text_configuration->font_item.prop.angle = angle;
            float act_angle = angle.has_value() ? *angle : .0f;
            do_rotate(act_angle - prev_angle);
        }
    }

    // when more collection add selector
    if (font_file->infos.size() > 1) {
        ImGui::Text("%s", tr.collection.c_str());
        ImGui::SameLine(m_gui_cfg->advanced_input_offset);
        ImGui::SetNextItemWidth(m_gui_cfg->input_width);
        unsigned int selected = font_prop.collection_number.has_value() ?
                               *font_prop.collection_number : 0;
        if (ImGui::BeginCombo("## Font collection", std::to_string(selected).c_str())) {
            for (unsigned int i = 0; i < font_file->infos.size(); ++i) {
                ImGui::PushID(1 << (10 + i));
                bool is_selected = (i == selected);
                if (ImGui::Selectable(std::to_string(i).c_str(), is_selected)) {
                    if (i == 0) font_prop.collection_number.reset();
                    else font_prop.collection_number = i;
                    exist_change = true;
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        } else if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("%s", _u8L("Select from True Type Collection.").c_str());
        }
    }

    if (exist_change) {
        m_font_manager.clear_glyphs_cache();
        process();
    }
#ifdef ALLOW_DEBUG_MODE
    ImGui::Text("family = %s", (font_prop.family.has_value() ?
                                    font_prop.family->c_str() :
                                    " --- "));
    ImGui::Text("face name = %s", (font_prop.face_name.has_value() ?
                                       font_prop.face_name->c_str() :
                                       " --- "));
    ImGui::Text("style = %s",
                (font_prop.style.has_value() ? font_prop.style->c_str() :
                                                 " --- "));
    ImGui::Text("weight = %s", (font_prop.weight.has_value() ?
                                    font_prop.weight->c_str() :
                                    " --- "));

    std::string descriptor = fi.path;
    ImGui::Text("descriptor = %s", descriptor.c_str());
#endif // ALLOW_DEBUG_MODE
}

void GLGizmoEmboss::set_minimal_window_size(bool is_advance_edit_style)
{
    ImVec2 window_size = ImGui::GetWindowSize();
    const ImVec2& min_win_size_prev = get_minimal_window_size();
    //ImVec2 diff(window_size.x - min_win_size_prev.x,
    //            window_size.y - min_win_size_prev.y);
    float diff_y = window_size.y - min_win_size_prev.y;
    m_is_advanced_edit_style = is_advance_edit_style;
    const ImVec2 &min_win_size = get_minimal_window_size();
    ImGui::SetWindowSize(ImVec2(0.f, min_win_size.y + diff_y),
                         ImGuiCond_Always);
}

const ImVec2 &GLGizmoEmboss::get_minimal_window_size() const
{
    return m_is_advanced_edit_style ?
               m_gui_cfg->minimal_window_size_with_advance :
               m_gui_cfg->minimal_window_size;
}

#ifdef ALLOW_ADD_FONT_BY_OS_SELECTOR
bool GLGizmoEmboss::choose_font_by_wxdialog()
{
    wxFontData data;
    data.EnableEffects(false);
    data.RestrictSelection(wxFONTRESTRICT_SCALABLE);
    // set previous selected font
    FontItem &selected_font_item = m_font_manager.get_font_item();
    if (selected_font_item.type == WxFontUtils::get_actual_type()) {
        std::optional<wxFont> selected_font = WxFontUtils::load_wxFont(
            selected_font_item.path);
        if (selected_font.has_value()) data.SetInitialFont(*selected_font);
    }

    wxFontDialog font_dialog(wxGetApp().mainframe, data);
    if (font_dialog.ShowModal() != wxID_OK) return false;

    data                = font_dialog.GetFontData();
    wxFont   wx_font       = data.GetChosenFont();
    size_t   font_index = m_font_manager.get_fonts().size();
    FontItem font_item  = WxFontUtils::get_font_item(wx_font);
    m_font_manager.add_font(font_item);

    // Check that deserialization NOT influence font
    // false - use direct selected wxFont in dialog
    // true - use font item (serialize and deserialize wxFont)
    bool use_deserialized_font = false;

    // Try load and use new added font
    if ((use_deserialized_font && !m_font_manager.load_font(font_index)) ||
        (!use_deserialized_font && !m_font_manager.load_font(font_index, wx_font))) {
        m_font_manager.erase(font_index);
        wxString message = GUI::format_wxstr(
            _L("Font '%1%' can't be used. Please select another."),
            font_item.name);
        wxString      title = _L("Selected font is NOT True-type.");
        MessageDialog not_loaded_font_message(nullptr, message, title, wxOK);
        not_loaded_font_message.ShowModal();
        return choose_font_by_wxdialog();
    }

    // fix dynamic creation of italic font
    const auto& cn = m_font_manager.get_font_prop().collection_number;
    unsigned int font_collection = cn.has_value() ? *cn : 0;
    if (WxFontUtils::is_italic(wx_font) &&
        !Emboss::is_italic(*m_font_manager.get_font_file(), font_collection)) {
        m_font_manager.get_font_item().prop.skew = 0.2;
    }
    return true;
}
#endif // ALLOW_ADD_FONT_BY_OS_SELECTOR

#ifdef ALLOW_ADD_FONT_BY_FILE
bool GLGizmoEmboss::choose_true_type_file()
{
    wxArrayString input_files;
    wxString      fontDir      = wxEmptyString;
    wxString      selectedFile = wxEmptyString;
    wxFileDialog  dialog(nullptr, _L("Choose one or more files (TTF, TTC):"),
                        fontDir, selectedFile, file_wildcards(FT_FONTS),
                        wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (dialog.ShowModal() == wxID_OK) dialog.GetPaths(input_files);
    if (input_files.IsEmpty()) return false;
    size_t index = m_font_manager.get_fonts().size();
    // use first valid font
    for (auto &input_file : input_files) {
        std::string path = std::string(input_file.c_str());
        std::string name = get_file_name(path);
        //make_unique_name(name, m_font_list);
        const FontProp& prop = m_font_manager.get_font_prop();
        FontItem fi{ name, path, FontItem::Type::file_path, prop };
        m_font_manager.add_font(fi);
        // set first valid added font as active
        if (m_font_manager.load_font(index)) return true;
        m_font_manager.erase(index);       
    }
    return false;
}
#endif // ALLOW_ADD_FONT_BY_FILE

bool GLGizmoEmboss::choose_svg_file()
{
    wxArrayString input_files;
    wxString      fontDir      = wxEmptyString;
    wxString      selectedFile = wxEmptyString;
    wxFileDialog  dialog(nullptr, _L("Choose SVG file:"), fontDir,
                        selectedFile, file_wildcards(FT_SVG),
                        wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (dialog.ShowModal() == wxID_OK) dialog.GetPaths(input_files);
    if (input_files.IsEmpty()) return false;
    if (input_files.size() != 1) return false;
    auto &      input_file = input_files.front();
    std::string path       = std::string(input_file.c_str());
    std::string name       = get_file_name(path);

    NSVGimage *image = nsvgParseFromFile(path.c_str(), "mm", 96.0f);
    ExPolygons polys = NSVGUtils::to_ExPolygons(image);
    nsvgDelete(image);

    BoundingBox bb;
    for (const auto &p : polys) bb.merge(p.contour.points);
    const FontProp &fp = m_font_manager.get_font_item().prop;
    float scale   = fp.size_in_mm / std::max(bb.max.x(), bb.max.y());
    auto  project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(fp.emboss / scale), scale);
    indexed_triangle_set its = Emboss::polygons2model(polys, *project);
    return false;
    // test store:
    // for (auto &poly : polys) poly.scale(1e5);
    // SVG svg("converted.svg", BoundingBox(polys.front().contour.points));
    // svg.draw(polys);
    //return add_volume(name, its);
}

TextConfiguration GLGizmoEmboss::create_configuration()
{
    if (!m_font_manager.is_activ_font()) 
        return TextConfiguration{FontItem{}, m_text};
    
    FontItem &fi = m_font_manager.get_font_item();
    // actualize font path
    if (fi.type == WxFontUtils::get_actual_type()) {
        std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();
        if (wx_font.has_value())
            fi.path = WxFontUtils::store_wxFont(*wx_font);
    }
    return TextConfiguration{m_font_manager.get_font_item(), m_text};
}

bool GLGizmoEmboss::load_configuration(ModelVolume *volume)
{
    if (volume == nullptr) return false;
    if (!volume->text_configuration.has_value()) return false;

    TextConfiguration &configuration = *volume->text_configuration;
    FontItem &         c_font_item   = configuration.font_item;

    // try to find font in local font list
    auto is_config = [&c_font_item](const FontManager::Item &font_item) -> bool {
        const FontItem &fi = font_item.font_item;
        return fi.path == c_font_item.path && fi.prop == c_font_item.prop;
    };
    const auto& styles = m_font_manager.get_styles();
    auto it = std::find_if(styles.begin(), styles.end(), is_config);
    bool found_font = it != styles.end();
    size_t font_index;
    if (!found_font) {
        // font is not in list
        // add font to list
        font_index = styles.size();
        m_font_manager.add_font(c_font_item);        
    } else {
        // font is found in list
        font_index = it - styles.begin();
    }

    m_text      = configuration.text;
    m_volume    = volume;

    if (!m_font_manager.load_font(font_index)) {
        // create similar font
        auto wx_font = WxFontUtils::create_wxFont(c_font_item, configuration.font_item.prop);
        if (wx_font.has_value()) {
            // fix not loadable font item
            FontItem &fi = m_font_manager.get_font_item();
            FontItem fi_new = WxFontUtils::get_font_item(*wx_font);
            fi_new.name = fi.name; // use previous name
            fi = fi_new; // rewrite font item
            fi.prop = configuration.font_item.prop;
            if (!m_font_manager.load_font(font_index, *wx_font)) return false;
        } else {
            // can't create similar font use previous
            m_font_manager.erase(font_index);
        }
        create_notification_not_valid_font(configuration);
    }
    return true;
}

void GLGizmoEmboss::create_notification_not_valid_font(
    const TextConfiguration &tc)
{
    // not neccessary, but for sure that old notification doesnt exist
    if (m_exist_notification) remove_notification_not_valid_font();
    m_exist_notification = true;

    auto type = NotificationType::UnknownFont;
    auto level =
        NotificationManager::NotificationLevel::WarningNotificationLevel;

    const auto &fi            = m_font_manager.get_font_item();
    const auto &origin_family = tc.font_item.prop.face_name;
    const auto &actual_family = fi.prop.face_name;

    const std::string &origin_font_name = origin_family.has_value() ?
                                              *origin_family :
                                              tc.font_item.path;
    const std::string &actual_font_name = actual_family.has_value() ?
                                              *actual_family :
                                              fi.name;

    std::string text =
        GUI::format(_L("Can't load exactly same font(\"%1%\"), "
                       "Aplication select similar one(\"%2%\"). "
                       "When you edit text, similar font will be applied."),
                    origin_font_name, actual_font_name);
    auto notification_manager = wxGetApp().plater()->get_notification_manager();
    notification_manager->push_notification(type, level, text);
}

void GLGizmoEmboss::remove_notification_not_valid_font()
{
    if (!m_exist_notification) return;
    m_exist_notification      = false;
    auto type                 = NotificationType::UnknownFont;
    auto notification_manager = wxGetApp().plater()->get_notification_manager();
    notification_manager->close_notification_of_type(type);
}

std::string GLGizmoEmboss::create_volume_name()
{
    const size_t &max_len = m_gui_cfg->max_count_char_in_volume_name;
    return _u8L("Text") + " - " +
           ((m_text.size() > max_len) ?
                (m_text.substr(0, max_len - 3) + " ..") :
                m_text);
}

void GLGizmoEmboss::init_icons()
{
    // icon order has to match the enum IconType
    std::vector<std::string> filenames{
        "edit_button.svg",
        "delete.svg",
        "add_copies.svg", 
        "save.svg", 
        "undo.svg",    
        "make_italic.svg",
        "make_unitalic.svg",
        "make_bold.svg",
        "make_unbold.svg",   
        "search.svg",
        "open.svg",
        "revert_all_.svg",
        "add_text_part.svg",
        "add_text_negative.svg",
        "add_text_modifier.svg"
    };
    assert(filenames.size() == static_cast<size_t>(IconType::_count));
    std::string path = resources_dir() + "/icons/";
    for (std::string &filename : filenames) filename = path + filename;

    // state order has to match the enum IconState
    std::vector<std::pair<int, bool>> states;
    states.push_back(std::make_pair(1, false)); // Activable
    states.push_back(std::make_pair(0, false));  // Hovered
    states.push_back(std::make_pair(2, false)); // Disabled

    bool compress = false;
    bool is_loaded = m_icons_texture.load_from_svg_files_as_sprites_array(
        filenames, states, m_gui_cfg->icon_width, compress);

    if (!is_loaded ||
        (size_t)m_icons_texture.get_width() < (states.size() * m_gui_cfg->icon_width) ||
        (size_t)m_icons_texture.get_height() < (filenames.size() * m_gui_cfg->icon_width)) { 
        // bad load of icons, but all usage of m_icons_texture check that texture is initialized
        assert(false);
        m_icons_texture.reset();
    }
}

void GLGizmoEmboss::draw_icon(IconType icon, IconState state, ImVec2 size)
{
    // canot draw count
    assert(icon != IconType::_count);
    if (icon == IconType::_count) return;

    unsigned int icons_texture_id = m_icons_texture.get_id();
    int          tex_width        = m_icons_texture.get_width();
    int          tex_height       = m_icons_texture.get_height();
    // is icon loaded
    if ((icons_texture_id == 0) || (tex_width <= 1) || (tex_height <= 1)){
        ImGui::Text("▮");
        return;
    }
        
    int icon_width = m_gui_cfg->icon_width;
    ImTextureID tex_id = (void *) (intptr_t) (GLuint) icons_texture_id;
    int start_x = static_cast<unsigned>(state) * (icon_width + 1) + 1,
        start_y = static_cast<unsigned>(icon) * (icon_width + 1) + 1;

    ImVec2 uv0(start_x / (float) tex_width, 
               start_y / (float) tex_height);
    ImVec2 uv1((start_x + icon_width) / (float) tex_width,
               (start_y + icon_width) / (float) tex_height);
        
    if (size.x < 1 || size.y < 1)
        size = ImVec2(m_gui_cfg->icon_width, m_gui_cfg->icon_width);

    ImGui::Image(tex_id, size, uv0, uv1);
}

void GLGizmoEmboss::draw_transparent_icon()
{
    unsigned int icons_texture_id = m_icons_texture.get_id();
    int          tex_width        = m_icons_texture.get_width();
    int          tex_height       = m_icons_texture.get_height();
    // is icon loaded
    if ((icons_texture_id == 0) || (tex_width <= 1) || (tex_height <= 1)) {
        ImGui::Text("▯");
        return;
    }

    ImTextureID tex_id = (void *) (intptr_t) (GLuint) icons_texture_id;
    int    icon_width = m_gui_cfg->icon_width;
    ImVec2 icon_size(icon_width, icon_width);
    ImVec2 pixel_size(1.f / tex_width, 1.f / tex_height);
    // zero pixel is transparent in texture
    ImGui::Image(tex_id, icon_size, ImVec2(0, 0), pixel_size);
}

bool GLGizmoEmboss::draw_clickable(
    IconType icon, IconState state, 
    IconType hover_icon, IconState hover_state)
{
    // check of hover
    float cursor_x = ImGui::GetCursorPosX();
    draw_transparent_icon();
    ImGui::SameLine(cursor_x);

    if (ImGui::IsItemHovered()) {
        // redraw image
        draw_icon(hover_icon, hover_state);
    } else {
        // redraw normal image
        draw_icon(icon, state);
    }
    return ImGui::IsItemClicked();
}

bool GLGizmoEmboss::draw_button(IconType icon, bool disable)
{
    if (disable) {
        draw_icon(icon, IconState::disabled);
        return false;
    }
    return draw_clickable(
        icon, IconState::activable,
        icon, IconState::hovered
    );
}

static std::string APP_CONFIG_ACTIVE_FONT = "activ_font";
FontList GLGizmoEmboss::load_font_list_from_app_config(const AppConfig *cfg, size_t &activ_font_index)
{    
    FontList    result;
    unsigned    index        = 1;
    std::string section_name = FontListSerializable::create_section_name(index++);
    while (cfg->has_section(section_name)) {
        std::optional<FontItem> fi = FontListSerializable::load_font_item(cfg->get_section(section_name));
        if (fi.has_value()) result.emplace_back(*fi);
        section_name = FontListSerializable::create_section_name(index++);
    }
    if (result.empty())
        return create_default_font_list();

    // read selected (activ) font
    if (cfg->has_section(AppConfig::SECTION_FONT)) {
        auto section = cfg->get_section(AppConfig::SECTION_FONT);
        auto it      = section.find(APP_CONFIG_ACTIVE_FONT);
        if (it != section.end()) {
            size_t active_font = static_cast<size_t>(std::atoi(it->second.c_str()));
            if (active_font > 0 && active_font <= result.size()) {
                // stored index start from 1 (readablity in config file)
                activ_font_index = active_font - 1;
            }
        }
    }
    return result;
}

void GLGizmoEmboss::store_font_list_to_app_config()
{
    AppConfig *cfg   = wxGetApp().app_config;
    unsigned   index = 1;
    const auto& styles = m_font_manager.get_styles();

    // store actual font index
    size_t activ_index = &m_font_manager.get_activ_style() - &styles.front();
    cfg->clear_section(AppConfig::SECTION_FONT);
    // activ font first index is +1 to correspond with section name
    std::string activ_font = std::to_string(activ_index + 1);
    cfg->set(AppConfig::SECTION_FONT, APP_CONFIG_ACTIVE_FONT, activ_font);

    for (const auto &item : styles) {
        const FontItem &fi = item.font_item;
        // skip file paths + fonts from other OS(loaded from .3mf)
        if (fi.type != WxFontUtils::get_actual_type()) continue;
        FontListSerializable::store_font_item(*cfg, fi, index++);
    }
    fill_stored_font_items();

    // remove rest of font sections
    std::string section_name = FontListSerializable::create_section_name(index);
    while (cfg->has_section(section_name)) {
        cfg->clear_section(section_name);
        section_name = FontListSerializable::create_section_name(++index);
    }
    cfg->save();
}

bool GLGizmoEmboss::is_text_object(const ModelVolume *text) {
    if (text == nullptr) return false;
    if (!text->text_configuration.has_value()) return false;
    if (text->type() != ModelVolumeType::MODEL_PART) return false;
    for (const ModelVolume *v : text->get_object()->volumes) {
        if (v == text) continue;
        if (v->type() == ModelVolumeType::MODEL_PART) return false;
    }
    return true;
}

//void GLGizmoEmboss::store_font_item_to_app_config() const
//{
//    AppConfig *cfg = wxGetApp().app_config;
//    // index of section start from 1
//    const auto &act_item = m_font_manager.get_activ_style();
//    const FontItem &fi = act_item.font_item;
//
//    //size_t index = &m_font_manager.get_activ_style() -
//    //               &m_font_manager.get_fonts().front();
//    // fix index when, not serialized font is in list
//    size_t index = 0;
//    for (const auto &item : m_font_manager.get_fonts()) {
//        if (fi.type != WxFontUtils::get_actual_type()) continue;
//        if (&item == &act_item) break;
//        ++index;
//    }
//    
//    FontListSerializable::store_font_item(*cfg, fi, index);
//}

std::string GLGizmoEmboss::get_file_name(const std::string &file_path)
{
    size_t pos_last_delimiter = file_path.find_last_of("/\\");
    size_t pos_point          = file_path.find_last_of('.');
    size_t offset             = pos_last_delimiter + 1;
    size_t count              = pos_point - pos_last_delimiter - 1;
    return file_path.substr(offset, count);
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
