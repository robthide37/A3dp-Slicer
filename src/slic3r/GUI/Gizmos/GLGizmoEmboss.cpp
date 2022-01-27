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

// TODO: remove include
#include "libslic3r/SVG.hpp"      // debug store
#include "libslic3r/Geometry.hpp" // covex hull 2d

#include "libslic3r/NSVGUtils.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/ClipperUtils.hpp" // union_ex
#include "libslic3r/AppConfig.hpp"    // store/load font list
#include "libslic3r/MapUtils.hpp"
#include "libslic3r/Format/OBJ.hpp" // load obj file for default object

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
#define SHOW_IMGUI_ATLAS
#endif // ALLOW_DEBUG_MODE

#define ALLOW_ADD_FONT_BY_FILE
#define ALLOW_ADD_FONT_BY_OS_SELECTOR

using namespace Slic3r;
using namespace Slic3r::GUI;

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -2)
    , m_volume(nullptr)
    , m_exist_notification(false)
    , m_is_initialized(false) // initialize on first opening gizmo
    , m_rotate_gizmo(parent, GLGizmoRotate::Axis::Z) // grab id = 2 (Z axis)
    , m_font_manager(m_imgui->get_glyph_ranges())
{
    m_rotate_gizmo.set_group_id(0);
    // TODO: add suggestion to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web    
}

void GLGizmoEmboss::set_fine_position()
{
    const Selection &            selection = m_parent.get_selection();
    const Selection::IndicesList indices   = selection.get_volume_idxs();
    // no selected volume
    if (indices.empty()) return;
    const GLVolume *volume = selection.get_volume(*indices.begin());
    // bad volume selected (e.g. deleted one)
    if (volume == nullptr) return;

    const Camera &camera = wxGetApp().plater()->get_camera();
    Polygon hull = CameraUtils::create_hull2d(camera, *volume);

    const ImVec2 &windows_size = get_minimal_window_size();
    ImVec2 offset = ImGuiWrapper::suggest_location(windows_size, hull);
    m_set_window_offset = offset;
    return;

    Polygon rect({Point(offset.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y + windows_size.y),
                  Point(offset.x, offset.y + windows_size.y)});
    ImGuiWrapper::draw(hull);
    ImGuiWrapper::draw(rect);
}

#ifdef ALLOW_DEBUG_MODE
static void draw_fine_position(const Selection &selection)
{
    const Selection::IndicesList indices = selection.get_volume_idxs();
    // no selected volume
    if (indices.empty()) return;
    const GLVolume *volume = selection.get_volume(*indices.begin());
    // bad volume selected (e.g. deleted one)
    if (volume == nullptr) return;

    const Camera &camera = wxGetApp().plater()->get_camera();
    Slic3r::Polygon hull   = CameraUtils::create_hull2d(camera, *volume);

    ImVec2 windows_size(174, 202);
    ImVec2 offset       = ImGuiWrapper::suggest_location(windows_size, hull);
    Slic3r::Polygon rect(
        {Point(offset.x, offset.y), Point(offset.x + windows_size.x, offset.y),
         Point(offset.x + windows_size.x, offset.y + windows_size.y),
         Point(offset.x, offset.y + windows_size.y)});
    ImGuiWrapper::draw(hull);
    ImGuiWrapper::draw(rect);
}
#endif // ALLOW_DEBUG_MODE

void GLGizmoEmboss::create_volume(ModelVolumeType volume_type, const Vec2d& mouse_pos)
{
    assert(volume_type == ModelVolumeType::MODEL_PART ||
           volume_type == ModelVolumeType::NEGATIVE_VOLUME ||
           volume_type == ModelVolumeType::PARAMETER_MODIFIER);

    if (!m_is_initialized) initialize();
    const Selection &selection = m_parent.get_selection();
    if(selection.is_empty()) return;

    set_default_text();
    
    // By position of cursor create transformation to put text on surface of model
    Transform3d transformation;        
    const ModelObjectPtrs &objects = wxGetApp().plater()->model().objects;
    m_raycast_manager.actualize(objects);
    auto hit = m_raycast_manager.unproject(mouse_pos);
    if (hit.has_value()) {
        transformation = Emboss::create_transformation_onto_surface(hit->position, hit->normal);
    } else {
        // there is no hit with object         
        // TODO: calculate X,Y offset position for lay on platter by mouse position
        transformation = Transform3d::Identity();
    }

    create_emboss_volume(create_mesh(), transformation, create_volume_name(),
                         create_configuration(), volume_type,
                         selection.get_object_idx());
}

bool GLGizmoEmboss::on_mouse_for_rotation(const wxMouseEvent &mouse_event)
{
    if (mouse_event.Dragging()) {
        if (m_dragging) {
            // temporary rotation
            TransformationType transformation_type(
                TransformationType::Local_Relative_Independent);
            Vec3d rotation(0., 0., m_rotate_gizmo.get_angle());
            m_parent.get_selection().rotate(rotation, transformation_type);
        }
    } else if (mouse_event.LeftUp()) {
        if (m_dragging) {
            // apply rotation
            m_parent.do_rotate(L("Text-Rotate"));
        }
    }
    return false;
}

bool GLGizmoEmboss::on_mouse_for_translate(const wxMouseEvent &mouse_event)
{
    // filter events
    if (!mouse_event.Dragging() && 
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

    RaycastManager::SkipVolume skip(m_volume->id().id);
    // detect start text dragging
    if (mouse_event.LeftDown()) {
        // initialize raycasters
        // TODO: move to job, for big scene it slow down
        m_raycast_manager.actualize(objects, &skip);
        return false;
    }

    // wxCoord == int --> wx/types.h
    Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
    Vec2d mouse_pos = mouse_coord.cast<double>();
    auto hit = m_raycast_manager.unproject(mouse_pos, &skip);
    if (!hit.has_value()) { 
        // there is no hit
        // show common translation of object
        m_parent.toggle_model_objects_visibility(true);
        m_temp_transformation = {};
        return false; 
    }
        
    Transform3d object_trmat = m_raycast_manager.get_transformation(hit->tr_key);
    Transform3d trmat = Emboss::create_transformation_onto_surface(hit->position, hit->normal);
    if (mouse_event.Dragging()) {
        // hide common dragging of object
        m_parent.toggle_model_objects_visibility(false, m_volume->get_object(), gl_volume->instance_idx(), m_volume);

        // Show temporary position
        // TODO: store z-rotation and aply after transformation matrix
        m_temp_transformation = object_trmat * trmat;
    } else if (mouse_event.LeftUp()) {

        // TODO: Disable apply common transformation after draggig
        // Call after is used for apply transformation after common dragging to rewrite it
        ModelVolume *mv = m_volume;
        wxGetApp().plater()->CallAfter([trmat, mv]() {
            mv->set_transformation(trmat);
        });

        m_parent.toggle_model_objects_visibility(true);
        // Apply temporary position
        m_temp_transformation = {};     
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
    std::array<float, 4> gray_color = {.6f, .6f, .6f, .3f};
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
        const GLVolume& gl_volume = *selection.get_volume(*selection.get_volume_idxs().begin());
        glsafe(::glPushMatrix());
        glsafe(::glMultMatrixd(m_temp_transformation->data()));                
        GLShaderProgram *shader = wxGetApp().get_shader("gouraud_light");
        shader->start_using();

        // dragging object must be selected so draw it with correct color
        //auto color = gl_volume.color;
        //auto color = gl_volume.render_color;
        auto color = GLVolume::SELECTED_COLOR;
        // Set transparent color for NEGATIVE_VOLUME & PARAMETER_MODIFIER

        bool is_transparent = m_volume->type() != ModelVolumeType::MODEL_PART;        
        if (is_transparent) {
            color[3] = 0.5f;
            glsafe(::glEnable(GL_BLEND));
            glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));
        }
        shader->set_uniform("uniform_color", color);

        glsafe(::glEnable(GL_DEPTH_TEST));
        gl_volume.indexed_vertex_array.render();
        glsafe(::glDisable(GL_DEPTH_TEST));

        if (is_transparent) glsafe(::glDisable(GL_BLEND));

        shader->stop_using();
        glsafe(::glPopMatrix());
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

void GLGizmoEmboss::on_render_input_window(float x, float y, float bottom_limit)
{
    initialize();
    check_selection();

    // TODO: fix width - showing scroll in first draw of advanced.
    const ImVec2 &min_window_size = get_minimal_window_size();
    ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, min_window_size);

#ifdef ALLOW_DEBUG_MODE
    // draw suggested position of window
    draw_fine_position(m_parent.get_selection());
#endif // ALLOW_DEBUG_MODE

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
    ImGui::Begin(on_get_name().c_str(), &is_open, flag);
    draw_window();
    ImGui::End();

    if (!is_open) close();
    ImGui::PopStyleVar(); // WindowMinSize
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
        if (!m_is_initialized) initialize();

        // to reload fonts from system, when install new one
        wxFontEnumerator::InvalidateCache();

        const Selection &selection = m_parent.get_selection();
        bool create_new_object = selection.is_empty();
        // When add Text on empty plate, Create new object with volume
        if (create_new_object) {
            set_default_text();
            create_emboss_object(create_mesh(), create_volume_name(), create_configuration());

            // gizmo will open when successfuly create new object
            GLGizmoBase::m_state = GLGizmoBase::Off;
            return;
        }

        // Try(when exist) set configuration by volume
        load_configuration(get_selected_volume());

        // change position of just opened emboss window
        set_fine_position();

        // when open by hyperlink it needs to show up
        // or after key 'T' windows doesn't appear
        m_parent.reload_scene(true); 
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

void GLGizmoEmboss::initialize()
{
    if (m_is_initialized) return;
    m_is_initialized = true;

    GuiCfg cfg; // initialize by default values;
    
    float line_height = ImGui::GetTextLineHeight();
    float line_height_with_spacing = ImGui::GetTextLineHeightWithSpacing();
    float space = line_height_with_spacing - line_height;

    cfg.max_font_name_width = ImGui::CalcTextSize("Maximal font name").x;
    cfg.icon_width = line_height;
    float icon_width_with_spacing = cfg.icon_width + space;
    float scroll_width = icon_width_with_spacing; // TODO: fix it
    cfg.combo_font_width = cfg.max_font_name_width + space
                                  + icon_width_with_spacing
                                  + scroll_width;
    cfg.delete_pos_x = cfg.max_font_name_width + space;
    int count_line_of_text = 3;
    cfg.text_size = ImVec2(-FLT_MIN, line_height_with_spacing * count_line_of_text);
    ImVec2 letter_m_size = ImGui::CalcTextSize("M");
    int count_letter_M_in_input = 6;
    cfg.advanced_input_width = letter_m_size.x * count_letter_M_in_input;
    GuiCfg::Translations &tr = cfg.translations;
    tr.font                  = _u8L("Font");
    tr.size                  = _u8L("Height");
    tr.depth                 = _u8L("Depth");
    cfg.style_edit_text_width =
        3 * space + ImGui::GetTreeNodeToLabelSpacing() +
        std::max(ImGui::CalcTextSize(tr.font.c_str()).x,
                 std::max(ImGui::CalcTextSize(tr.size.c_str()).x,
                          ImGui::CalcTextSize(tr.depth.c_str()).x));

    // calculate window size
    const ImGuiStyle &style = ImGui::GetStyle();
    float window_title = line_height + 2*style.FramePadding.y;
    float input_height = line_height_with_spacing + 2*style.FramePadding.y;
    float tree_header  = line_height_with_spacing;
    float window_height = 
        window_title + // window title
        cfg.text_size.y +  // text field
        input_height * 3 + // type Radios + style selector + close button
        tree_header +      // Edit style
        2 * style.WindowPadding.y;
    float window_width = cfg.combo_font_width + style.WindowPadding.x * 2;
    cfg.minimal_window_size = ImVec2(window_width, window_height);

    float addition_edit_height        = input_height * 3 + tree_header;
    cfg.minimal_window_size_with_edit = ImVec2(cfg.minimal_window_size.x,
                                               cfg.minimal_window_size.y +
                                                   addition_edit_height);

    float advance_height = input_height * 4;
    cfg.minimal_window_size_with_advance =
        ImVec2(cfg.minimal_window_size_with_edit.x,
               cfg.minimal_window_size_with_edit.y + advance_height);

    cfg.min_style_image_height = line_height_with_spacing;
    cfg.max_style_image_width  = cfg.max_font_name_width -
                                2 * style.FramePadding.x;

    m_gui_cfg.emplace(cfg);

    // TODO: What to do when icon was NOT loaded? Generate them?
    bool success = init_icons();
    assert(success);

    const AppConfig *app_cfg = wxGetApp().app_config;
    FontList font_list = load_font_list_from_app_config(app_cfg);
    m_font_manager.add_fonts(font_list);
    if (!m_font_manager.load_first_valid_font()) {
        FontList font_list = FontListSerializable::create_default_font_list();
        m_font_manager.add_fonts(font_list);
        // TODO: What to do when default fonts are not loadable?
        bool success = m_font_manager.load_first_valid_font();
        assert(success);
    }
    set_default_text();
}

void GLGizmoEmboss::set_default_text()
{
    m_text = _u8L("Embossed text");
}

Slic3r::TriangleMesh GLGizmoEmboss::create_default_mesh()
{
    // When cant load any font use default object loaded from file
    std::string  path = Slic3r::resources_dir() + "/data/embossed_text.stl";
    TriangleMesh triangle_mesh;
    if (!load_obj(path.c_str(), &triangle_mesh)) {
        // when can't load mesh use cube
        return TriangleMesh(its_make_cube(36., 4., 2.5));
    }
    return triangle_mesh;
}

Slic3r::TriangleMesh GLGizmoEmboss::create_mesh()
{
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new embosstext object
    std::shared_ptr<Emboss::FontFile>& font_file = m_font_manager.get_font_file();
    if (font_file == nullptr) return create_default_mesh();
    const FontItem  &fi = m_font_manager.get_font_item();
    TriangleMesh result = create_mesh(m_text.c_str(), *font_file, fi.prop);
    if (result.its.empty()) return create_default_mesh();
    return result;
}

Slic3r::TriangleMesh GLGizmoEmboss::create_mesh(const char *      text,
                                                Emboss::FontFile &font,
                                                const FontProp &  font_prop)
{
    ExPolygons shapes   = Emboss::text2shapes(font, text, font_prop);
    float      scale    = font_prop.size_in_mm / font.ascent;
    float      depth    = font_prop.emboss / scale;
    auto       projectZ = std::make_unique<Emboss::ProjectZ>(depth);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    return TriangleMesh(Emboss::polygons2model(shapes, project));
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

bool GLGizmoEmboss::process()
{
    // no volume is selected -> selection from right panel
    if (m_volume == nullptr) return false;

    // exist loaded font?
    std::shared_ptr<Emboss::FontFile>& font_file = m_font_manager.get_font_file();
    if (font_file == nullptr) return false;
    auto data = std::make_unique<EmbossData>(font_file,
                                             create_configuration(),
                                             create_volume_name(), m_volume);
        
    auto &worker = wxGetApp().plater()->get_ui_job_worker();
    replace_job(worker, std::make_unique<EmbossJob>(std::move(data)));
    
    // notification is removed befor object is changed by job
    remove_notification_not_valid_font();
    return true;
}

void GLGizmoEmboss::close()
{
    // close gizmo == open it again
    m_parent.get_gizmos_manager().open_gizmo(GLGizmosManager::Emboss);
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
    bool exist_font_file = m_font_manager.get_font_file() != nullptr;
    if (!exist_font_file) {
        m_imgui->text_colored(
            ImGuiWrapper::COL_ORANGE_LIGHT,
            _L("Warning: No font is selected. Select correct one."));
    }
    draw_text_input();
    draw_model_type();
    draw_style_list();
    if (ImGui::TreeNode(_u8L("Edit style").c_str())) {
        draw_style_edit();
        ImGui::TreePop();
        if (!m_is_edit_style)
            set_minimal_window_size(true, m_is_advanced_edit_style);
    } else if (m_is_edit_style)
        set_minimal_window_size(false, m_is_advanced_edit_style);

    if (ImGui::Button(_u8L("Close").c_str())) close();

    // Option to create text volume when reselecting volumes
    m_imgui->disabled_begin(!exist_font_file);
    if (m_volume == nullptr) {
        ImGui::SameLine();
        if (ImGui::Button(_u8L("Generate preview").c_str())) { 
            const Selection &s = m_parent.get_selection();
            auto selected_indices = s.get_instance_idxs();
            if (selected_indices.empty()) { 
                create_emboss_object(create_mesh(), create_volume_name(), create_configuration());
            } else {
                create_volume(ModelVolumeType::MODEL_PART);
            }
        }
    }
    m_imgui->disabled_end();

#ifdef SHOW_IMGUI_ATLAS
    auto &atlas = m_font_manager.m_imgui_font_atlas;
    ImGui::Image(atlas.TexID, ImVec2(atlas.TexWidth, atlas.TexHeight));
#endif // SHOW_IMGUI_ATLAS
}

void GLGizmoEmboss::draw_text_input()
{
    static const ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_AllowTabInput | ImGuiInputTextFlags_AutoSelectAll;

    ImFont *imgui_font = m_font_manager.get_imgui_font(m_text);
    bool    exist_font = imgui_font != nullptr && imgui_font->IsLoaded();
    if (exist_font) ImGui::PushFont(imgui_font);

    bool   exist_change   = false;
    float  window_height  = ImGui::GetWindowHeight();
    float  minimal_height = get_minimal_window_size().y;
    float  extra_height   = window_height - minimal_height;
    ImVec2 text_size(m_gui_cfg->text_size.x,
                     m_gui_cfg->text_size.y + extra_height);
    if (ImGui::InputTextMultiline("##Text", &m_text, text_size, flags)) {
        process();
        exist_change = true;
    }

    if (exist_font) ImGui::PopFont();

    // show warning about incorrectness view of font
    // TODO: add char gap and line gap
    std::string warning;
    const FontProp& prop = m_font_manager.get_font_prop();
    if (prop.skew.has_value())
        warning = prop.boldness.has_value() ?
                      _u8L("Italic & Bold is NOT shown") :
                      _u8L("Italic is NOT shown");
    else if (prop.boldness.has_value())
        warning = _u8L("Boldness is NOT shown");

    if (!warning.empty()) {
        ImVec2 cursor = ImGui::GetCursorPos();
        float width = ImGui::GetContentRegionAvailWidth();
        ImVec2 size = ImGui::CalcTextSize(warning.c_str());
        ImVec2 padding = ImGui::GetStyle().FramePadding;
        ImGui::SetCursorPos(ImVec2(width - size.x + padding.x,
                                   cursor.y - size.y - padding.y));
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, warning);
        ImGui::SetCursorPos(cursor);
    }

    // Extend font ranges
    // imgui_font has to be unused
    if (exist_change) m_font_manager.clear_imgui_font();
}

void GLGizmoEmboss::draw_font_list()
{
    class MyFontEnumerator : public wxFontEnumerator
    {
        wxArrayString m_facenames;
        wxFontEncoding m_encoding;
        bool m_fixed_width_only;
        bool m_is_init;
    public: 
        MyFontEnumerator(wxFontEncoding encoding, bool fixed_width_only)
            : m_encoding(encoding)
            , m_fixed_width_only(fixed_width_only)
            , m_is_init(false)
        {}
        const wxArrayString& get_facenames() const{ return m_facenames; }
        bool is_init() const { return m_is_init; }
        bool init() {
            if (m_is_init) return false;
            m_is_init = true;
            if (!wxFontEnumerator::EnumerateFacenames(m_encoding, m_fixed_width_only)) return false;
            if (m_facenames.empty()) return false;
            return true;
        }
    protected:
        virtual bool OnFacename(const wxString& facename) wxOVERRIDE {
            // vertical font start with @, we will filter it out
            if (facename.empty() || facename[0] == '@') return true;
            wxFont wx_font(wxFontInfo().FaceName(facename).Encoding(m_encoding));
            void *addr = WxFontUtils::can_load(wx_font);
            if (addr == nullptr) return true; // can't load
            m_facenames.Add(facename);
            return true;
        }
    };
    wxFontEncoding encoding = wxFontEncoding::wxFONTENCODING_SYSTEM;
    bool fixed_width_only = false;
    static MyFontEnumerator fontEnumerator(encoding, fixed_width_only);

    std::optional<wxFont> &wx_font_opt = m_font_manager.get_wx_font();
    wxString actual_face_name = wx_font_opt.has_value() ?
        wx_font_opt->GetFaceName() : "";
    const char * selected = (!actual_face_name.empty()) ?
        actual_face_name.ToUTF8().data() : " --- ";
    if (ImGui::BeginCombo("##font_selector", selected)) {
        if(!fontEnumerator.is_init()) fontEnumerator.init();
        const wxArrayString &face_names = fontEnumerator.get_facenames();
        //const wxArrayString &face_names = wxFontEnumerator::GetFacenames(encoding, fixed_width_only);
        ImGui::TextColored(ImGuiWrapper::COL_GREY_LIGHT, "count %d", static_cast<int>(face_names.size()));
        for (const wxString &face_name : face_names) {
            size_t index = &face_name - &face_names.front();
            ImGui::PushID(index);
            bool is_selected = (actual_face_name == face_name);
            if (ImGui::Selectable(face_name.ToUTF8().data(), is_selected) &&
                wxFontEnumerator::IsValidFacename(face_name)) {
                // Select font
                wxFont wx_font(wxFontInfo().FaceName(face_name).Encoding(encoding));
                m_font_manager.set_wx_font(wx_font);
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }

#ifdef ALLOW_ADD_FONT_BY_FILE
    ImGui::SameLine();
    // select font file by file browser
    if (draw_button(IconType::open_file)) {
        choose_true_type_file();
    } else if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("add file with font(.ttf, .ttc)").c_str());
#endif //  ALLOW_ADD_FONT_BY_FILE

#ifdef ALLOW_ADD_FONT_BY_OS_SELECTOR
    ImGui::SameLine();
    if (draw_button(IconType::system_selector)) {
        choose_font_by_wxdialog();
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
    bool is_last_solid_part = false;
    if (type == part) {
        is_last_solid_part = true;
        for (const ModelVolume* vol : m_volume->get_object()->volumes) {
            if (vol == m_volume) continue;
            if (vol->type() == part) {
                is_last_solid_part = false;
                break;
            }
        }
    }
    if (ImGui::RadioButton("modifier", type == modifier) && !is_last_solid_part)
        new_type = modifier;
    if(is_last_solid_part && ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("You can't change a type of the last solid part of the object.").c_str());

    ImGui::SameLine();
    if (ImGui::RadioButton("negative", type == negative) && !is_last_solid_part)
        new_type = negative;
    if(is_last_solid_part && ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("You can't change a type of the last solid part of the object.").c_str());

    ImGui::SameLine();
    if (ImGui::RadioButton("part", type == part))
        new_type = part;

    ImGui::SameLine();
    m_imgui->disabled_begin(true);
    ImGui::RadioButton("baked in", false);
    m_imgui->disabled_end();

    if (m_volume != nullptr && new_type.has_value() && !is_last_solid_part) {
        GUI_App &app    = wxGetApp();
        Plater * plater = app.plater();
        plater->take_snapshot(_L("Change Part Type"));         
        m_volume->set_type(*new_type);

        // inspiration in ObjectList::change_part_type()
        // how to view correct side panel with objects
        ObjectList *obj_list = app.obj_list();
        ModelVolume * volume = m_volume;
        wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(
            obj_list->get_selected_obj_idx(),
            [volume](const ModelVolume *vol) { return vol == volume; });
        if (!sel.IsEmpty()) obj_list->select_item(sel.front());

        // TODO: fix way of view - color after change from volume to negative
    }
}

void GLGizmoEmboss::draw_rename_style(bool start_rename)
{
    // rename modal window popup
    const char *       rename_popup_id = "Rename_font";
    static FontItem *  rename_item;
    static std::string new_name;
    if (start_rename && !ImGui::IsPopupOpen(rename_popup_id)) {
        ImGui::OpenPopup(rename_popup_id);
        rename_item = &m_font_manager.get_font_item();
        new_name    = rename_item->name; // initialize with original copy
    }

    if (ImGui::BeginPopupModal(rename_popup_id, 0,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        const std::string &original_font_name = rename_item->name;
        std::string        text_in_popup =
            GUI::format(_u8L("Change font name (%1%): "), original_font_name);
        text_in_popup += "\n" +
                         _u8L("NOTE: Name has to be unique in font list.");
        ImGui::Text("%s", text_in_popup.c_str());
        ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);

        bool is_unique = true;
        for (const auto &item : m_font_manager.get_fonts()) {
            const FontItem &fi = item.font_item;
            if (&fi == rename_item)
                continue; // could be same as original name
            if (fi.name == new_name) is_unique = false;
        }
        bool allow_change = is_unique && !new_name.empty();

        ImGuiInputTextFlags flags = ImGuiInputTextFlags_EnterReturnsTrue;
        if ((ImGui::InputText("##font name", &new_name, flags) && allow_change) ||
            m_imgui->button(_L("ok"), ImVec2(0.f, 0.f), allow_change)) {
            rename_item->name = new_name;
            m_font_manager.get_truncated_name() = "";
            m_font_manager.free_style_images();
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

void GLGizmoEmboss::draw_style_list() {
    const float &max_width = m_gui_cfg->max_font_name_width;
    std::optional<size_t> delete_index;
    const FontItem &actual_font_item = m_font_manager.get_font_item();
    std::string &trunc_name = m_font_manager.get_truncated_name();
    if (trunc_name.empty()) {
        // generate trunc name
        const std::string &current_name  = actual_font_item.name;
        trunc_name = ImGuiWrapper::trunc(current_name, max_width);
    }

    ImGui::Text("%s", _u8L("Style").c_str());
    ImGui::SameLine();
    ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
    if (ImGui::BeginCombo("##style_selector", trunc_name.c_str())) {
        m_font_manager.init_style_images(m_gui_cfg->max_style_image_width);
        const auto &fonts = m_font_manager.get_fonts();
        for (const auto &item : fonts) {
            size_t             index             = &item - &fonts.front();
            const FontItem &   fi                = item.font_item;
            const std::string &actual_style_name = fi.name;
            ImGui::PushID(actual_style_name.c_str());
            std::string name_truncated = 
                ImGuiWrapper::trunc(actual_style_name, max_width);

            bool is_selected = (&fi == &actual_font_item);
            ImGuiSelectableFlags_ flags =
                ImGuiSelectableFlags_AllowItemOverlap; // allow click buttons

            const FontManager::StyleImage &img = *item.image;
            ImVec2 select_size(0.f, std::max(img.tex_size.y, m_gui_cfg->min_style_image_height));
            if (ImGui::Selectable("##style_select", is_selected, flags, select_size)) {
                if (m_font_manager.load_font(index)) process();
            } else if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", actual_style_name.c_str());

            // reorder items
            if (ImGui::IsItemActive() && !ImGui::IsItemHovered()) {
                std::optional<size_t> other_index;
                if (ImGui::GetMouseDragDelta(0).y < 0.f) {
                    if (index > 0) other_index = index - 1;
                } else if ((index + 1) < fonts.size())
                    other_index = index + 1;

                if (other_index.has_value()) {
                    m_font_manager.swap(index, *other_index);
                    ImGui::ResetMouseDragDelta();
                }
            }

            // draw style name
            ImGui::SameLine();
            ImGui::Image(img.texture_id, img.tex_size, img.uv0, img.uv1);

            // delete button
            ImGui::SameLine(m_gui_cfg->delete_pos_x);
            if (draw_button(IconType::erase, is_selected) && !is_selected)
                delete_index = index;
            if (ImGui::IsItemHovered()) {
                std::string tooltip = (is_selected)?
                    GUI::format(_L("Active style \"%1%\" can't be deleted."), actual_style_name) :
                    GUI::format(_L("Delete \"%1%\" style."), actual_style_name);
                ImGui::SetTooltip("%s", tooltip.c_str());
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }

    // delete font item
    if (delete_index.has_value())
        m_font_manager.erase(*delete_index);
    
    ImGui::SameLine();
    bool start_rename = false;
    if (draw_button(IconType::rename)) start_rename = true;
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Rename actual style.").c_str());
    draw_rename_style(start_rename);

    ImGui::SameLine();
    if (draw_button(IconType::duplicate))
        m_font_manager.duplicate();
    
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", _u8L("Duplicate style.").c_str());
    

    // TODO: Is style changed against stored one
    bool is_changed = false;

    ImGui::SameLine();
    if (draw_button(IconType::save, !is_changed)) {
        // TODO: make save style
        store_font_list_to_app_config();
        //store_font_item_to_app_config();
    }
    if (ImGui::IsItemHovered()) 
        if (is_changed)
            ImGui::SetTooltip("%s", _u8L("Save current settings to selected style").c_str());
        else
            ImGui::SetTooltip("%s", _u8L("No changes to save into style").c_str());

    if (is_changed) {
        ImGui::SameLine();
        if (draw_button(IconType::undo)) {
           // TODO: make undo changes in style
        }    
        if (ImGui::IsItemHovered()) 
            ImGui::SetTooltip("%s", _u8L("Reload original value of selected style").c_str());
    }
}

bool GLGizmoEmboss::italic_button()
{
    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font(); 
    if (!wx_font.has_value()) { 
        draw_icon(IconType::italic, IconState::disabled);
        return false;
    }

    std::optional<float> &skew = m_font_manager.get_font_prop().skew;
    bool is_font_italic = skew.has_value() || WxFontUtils::is_italic(*wx_font);
    if (is_font_italic) {
        if (draw_button(IconType::unitalic)) {
            skew.reset();
            wx_font->SetStyle(wxFontStyle::wxFONTSTYLE_NORMAL);
            return true;
        }
        if (ImGui::IsItemHovered()) 
            ImGui::SetTooltip("%s", _u8L("Unset italic").c_str());
    } else {
        if (draw_button(IconType::italic)) {
            std::shared_ptr<Emboss::FontFile> &font_file =
                m_font_manager.get_font_file();
            bool is_set = WxFontUtils::set_italic(*wx_font, font_file);
            // add skew when wxFont can't set it
            if (!is_set) skew = 0.2f;
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Set italic").c_str());
    }
    return false;
}

bool GLGizmoEmboss::bold_button() {
    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();
    if (!wx_font.has_value()) {
        draw_icon(IconType::bold, IconState::disabled);
        return false;
    }
    
    std::optional<float> &boldness = m_font_manager.get_font_prop().boldness;
    bool is_font_bold = boldness.has_value() || WxFontUtils::is_bold(*wx_font);
    if (is_font_bold) {
        if (draw_button(IconType::unbold)) {
            boldness.reset();
            wx_font->SetWeight(wxFontWeight::wxFONTWEIGHT_NORMAL);
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Unset bold").c_str());
    } else {
        if (draw_button(IconType::bold)) {
            std::shared_ptr<Emboss::FontFile> &font_file = 
                m_font_manager.get_font_file();
            bool is_set = WxFontUtils::set_bold(*wx_font, font_file);
            // add boldness when wxFont can't set it
            if (!is_set) boldness = 20.f;
            return true;
        }
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", _u8L("Set bold").c_str());
    }
    return false;
}

void GLGizmoEmboss::draw_style_edit() {
    const GuiCfg::Translations &tr = m_gui_cfg->translations;
    ImGui::Text("%s", tr.font.c_str());
    ImGui::SameLine(m_gui_cfg->style_edit_text_width);
    ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
    draw_font_list();
    ImGui::SameLine();
    bool exist_change = false;
    exist_change |= italic_button();
    ImGui::SameLine();
    exist_change |= bold_button();

    FontItem &fi = m_font_manager.get_font_item();
    std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();

    // TODO: should not be there
    // when actual font not loaded try to load
    if (!wx_font.has_value() && fi.type == WxFontUtils::get_actual_type())
        wx_font = WxFontUtils::load_wxFont(fi.path);

    FontProp &font_prop = fi.prop;

    ImGui::Text("%s", tr.size.c_str());
    ImGui::SameLine(m_gui_cfg->style_edit_text_width);
    ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
    if (ImGui::InputFloat("##line height", &font_prop.size_in_mm, 0.1f, 1.f, "%.1f mm")) {
        // size can't be zero or negative
        if (font_prop.size_in_mm < std::numeric_limits<float>::epsilon()) 
            font_prop.size_in_mm = 10.f;

        // store font size into path
        if (fi.type == WxFontUtils::get_actual_type()) {
            if (wx_font.has_value()) {
                // TODO: check difference wx_font->Scale(float x)
                wx_font->SetPointSize(font_prop.size_in_mm);
                fi.path = WxFontUtils::store_wxFont(*wx_font);
            }
        }
        exist_change = true;
    }
    if (exist_change) {
        m_font_manager.clear_imgui_font();
        m_font_manager.free_style_images();
        process();
    }

    ImGui::Text("%s", tr.depth.c_str());
    ImGui::SameLine(m_gui_cfg->style_edit_text_width);
    ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
    if (ImGui::InputFloat("##size in Z", &font_prop.emboss, 0.1f, 0.25, "%.2f mm")) {
        process();
    }

    if (ImGui::TreeNode(_u8L("advanced").c_str())) {
        draw_advanced();
        ImGui::TreePop();
        if (!m_is_advanced_edit_style)
            set_minimal_window_size(true, true);
    } else if (m_is_advanced_edit_style) 
        set_minimal_window_size(true, false);    
}

void GLGizmoEmboss::draw_advanced()
{
    std::shared_ptr<Emboss::FontFile>& font_file = m_font_manager.get_font_file();
    if (font_file == nullptr) { 
        ImGui::Text("%s", _u8L("Advanced font options could be change only for corect font.\nStart with select correct font.").c_str());
        return;
    }

    FontProp &font_prop = m_font_manager.get_font_item().prop;
    bool   exist_change = false;    
    
    ImGui::SetNextItemWidth(2 * m_gui_cfg->advanced_input_width);
    if (ImGuiWrapper::input_optional_int(_u8L("CharGap[in font points]").c_str(), font_prop.char_gap)) {
        // char gap is stored inside of imgui font atlas
        m_font_manager.clear_imgui_font();
        exist_change = true;
    }

    ImGui::SetNextItemWidth(2*m_gui_cfg->advanced_input_width);
    if (ImGuiWrapper::input_optional_int(_u8L("LineGap [in font points]").c_str(), font_prop.line_gap))
        exist_change = true;

    ImGui::SetNextItemWidth(2 * m_gui_cfg->advanced_input_width);
    if (m_imgui->slider_optional_float(_u8L("Boldness [in font points]").c_str(), font_prop.boldness, -200.f, 200.f, "%.0f", 1.f, false, _L("tiny / wide chars")))
        exist_change = true;
    

    ImGui::SetNextItemWidth(2 * m_gui_cfg->advanced_input_width);
    if (m_imgui->slider_optional_float(_u8L("Italic [Skew ratio]").c_str(), font_prop.skew, -1.f, 1.f, "%.2f", 1.f, false, _L("italic strength")))
        exist_change = true;
    
    // when more collection add selector
    if (font_file->count > 1) {
        ImGui::SetNextItemWidth(m_gui_cfg->advanced_input_width);
        if (ImGui::BeginCombo(_u8L("Font collection").c_str(),
                              std::to_string(font_file->index).c_str())) {
            for (unsigned int i = 0; i < font_file->count; ++i) {
                ImGui::PushID(1 << (10 + i));
                bool is_selected = i == font_file->index;
                if (ImGui::Selectable(std::to_string(i).c_str(), is_selected)) {
                    font_file->index = i;
                    exist_change = true;
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
    }

    if (exist_change) {
        m_font_manager.free_style_images();
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

void GLGizmoEmboss::set_minimal_window_size(bool is_edit_style,
                                            bool is_advance_edit_style)
{
    ImVec2 window_size = ImGui::GetWindowSize();
    const ImVec2& min_win_size_prev = get_minimal_window_size();
    ImVec2 diff(window_size.x - min_win_size_prev.x,
                window_size.y - min_win_size_prev.y);
    float diff_y = ImGui::GetWindowSize().y - get_minimal_window_size().y;
    m_is_edit_style = is_edit_style;
    m_is_advanced_edit_style = is_advance_edit_style;
    const ImVec2 &min_win_size = get_minimal_window_size();
    ImGui::SetWindowSize(ImVec2(0.f, min_win_size.y + diff.y),
                         ImGuiCond_Always);
}

const ImVec2 &GLGizmoEmboss::get_minimal_window_size() const
{
    return (m_is_edit_style) ?
               ((m_is_advanced_edit_style) ?
                    m_gui_cfg->minimal_window_size_with_advance :
                    m_gui_cfg->minimal_window_size_with_edit) :
               m_gui_cfg->minimal_window_size;
}

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
        (!use_deserialized_font && !m_font_manager.load_font(font_index, wx_font)) ||
        !process()) {
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
    if (WxFontUtils::is_italic(wx_font) &&
        !Emboss::is_italic(*m_font_manager.get_font_file())) {
        m_font_manager.get_font_item().prop.skew = 0.2;
    }

    return true;
}

bool GLGizmoEmboss::choose_true_type_file()
{
    wxArrayString input_files;
    wxString      fontDir      = wxEmptyString;
    wxString      selectedFile = wxEmptyString;
    wxFileDialog  dialog(nullptr, _L("Choose one or more files (TTF, TTC):"),
                        fontDir, selectedFile, file_wildcards(FT_FONTS),
                        wxFD_OPEN | wxFD_MULTIPLE | wxFD_FILE_MUST_EXIST);
    if (dialog.ShowModal() == wxID_OK) dialog.GetPaths(input_files);
    if (input_files.IsEmpty()) return false;
    bool font_loaded = false;
    for (auto &input_file : input_files) {
        std::string path = std::string(input_file.c_str());
        std::string name = get_file_name(path);
        //make_unique_name(name, m_font_list);
        FontItem fi(name, path, FontItem::Type::file_path, FontProp());
        m_font_manager.add_font(fi);

        // set first valid added font as active
        if (!font_loaded) {
            size_t index = m_font_manager.get_fonts().size() - 1;
            if (!m_font_manager.load_font(index))
                m_font_manager.erase(index);
            else
                font_loaded = true;
        }
    }
    if (font_loaded) process();
    return font_loaded;
}

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
    FontItem &fi = m_font_manager.get_font_item();
    // actualize font path
    if (fi.type == WxFontUtils::get_actual_type()) {
        std::optional<wxFont> &wx_font = m_font_manager.get_wx_font();
        if (wx_font.has_value())
            fi.path = WxFontUtils::store_wxFont(*wx_font);
    }
    return TextConfiguration(m_font_manager.get_font_item(), m_text);
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
    const auto& fonts = m_font_manager.get_fonts();
    auto it = std::find_if(fonts.begin(), fonts.end(), is_config);
    bool found_font = it != fonts.end();
    size_t font_index;
    if (!found_font) {
        // font is not in list
        // add font to list
        font_index = fonts.size();
        m_font_manager.add_font(c_font_item);        
    } else {
        // font is found in list
        font_index = it - fonts.begin();
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

bool GLGizmoEmboss::init_icons()
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
        "open.svg"
    };
    assert(filenames.size() == static_cast<size_t>(IconType::_count));
    std::string path = resources_dir() + "/icons/";
    for (std::string &filename : filenames) filename = path + filename;

    // state order has to match the enum IconState
    std::vector<std::pair<int, bool>> states;
    states.push_back(std::make_pair(1, false)); // Activable
    states.push_back(std::make_pair(0, true));  // Hovered
    states.push_back(std::make_pair(2, false)); // Disabled

    unsigned int sprite_size_px = std::ceil(m_gui_cfg->icon_width);
    // make size pair number
    if (sprite_size_px % 2 != 0) ++sprite_size_px;
    bool compress = false;
    return m_icons_texture.load_from_svg_files_as_sprites_array(filenames,
                                                                states,
                                                                sprite_size_px,
                                                                compress);
}

void GLGizmoEmboss::draw_icon(IconType icon, IconState state)
{
    // canot draw count
    assert(icon != IconType::_count);
    if (icon == IconType::_count) return; 

    unsigned int icons_texture_id = m_icons_texture.get_id();
    int          tex_width        = m_icons_texture.get_width();
    int          tex_height       = m_icons_texture.get_height();

    // is icon loaded
    if ((icons_texture_id == 0) || (tex_width <= 1) || (tex_height <= 1))
        return;
    ImTextureID tex_id = (void *) (intptr_t) (GLuint) icons_texture_id;

    size_t count_icons  = static_cast<size_t>(IconType::_count);
    size_t count_states = 3; // activable | hovered | disabled
    ImVec2 icon_size(tex_width / count_states, tex_height / count_icons);

    ImVec2 start(static_cast<unsigned>(state) * icon_size.x,
                 static_cast<unsigned>(icon) * icon_size.y);

    ImVec2 uv0(start.x / tex_width, start.y / tex_height);
    ImVec2 uv1((start.x + icon_size.x) / tex_width,
               (start.y + icon_size.y) / tex_height);

    ImGui::Image(tex_id, icon_size, uv0, uv1);
}

bool GLGizmoEmboss::draw_button(IconType icon, bool disable)
{
    if (disable) {
        draw_icon(icon, IconState::disabled);
        return false;
    }

    float cursor_x = ImGui::GetCursorPosX();
    draw_icon(icon, IconState::activable);
    if (ImGui::IsItemClicked()) return true;
    if (ImGui::IsItemHovered()) {
        // redraw image over previous
        ImGui::SameLine(cursor_x);
        draw_icon(icon, IconState::hovered);
    }
    return false;
}

FontList GLGizmoEmboss::load_font_list_from_app_config(const AppConfig *cfg)
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
        return FontListSerializable::create_default_font_list();
    return result;
}

void GLGizmoEmboss::store_font_list_to_app_config() const
{
    AppConfig *cfg   = wxGetApp().app_config;
    unsigned   index = 1;
    for (const auto& item : m_font_manager.get_fonts()) {
        const FontItem &fi = item.font_item;
        // skip file paths + fonts from other OS(loaded from .3mf)
        if (fi.type != WxFontUtils::get_actual_type()) continue;
        FontListSerializable::store_font_item(*cfg, fi, index++);
    }

    // remove rest of font sections
    std::string section_name = FontListSerializable::create_section_name(index);
    while (cfg->has_section(section_name)) {
        cfg->clear_section(section_name);
        section_name = FontListSerializable::create_section_name(++index);
    }
}

void GLGizmoEmboss::store_font_item_to_app_config() const
{
    AppConfig *cfg = wxGetApp().app_config;
    // index of section start from 1
    const auto &act_item = m_font_manager.get_font();
    const FontItem &fi = act_item.font_item;

    //size_t index = &m_font_manager.get_font() -
    //               &m_font_manager.get_fonts().front();
    // fix index when, not serialized font is in list
    size_t index = 0;
    for (const auto &item : m_font_manager.get_fonts()) {
        if (fi.type != WxFontUtils::get_actual_type()) continue;
        if (&item == &act_item) break;
        ++index;
    }
    
    FontListSerializable::store_font_item(*cfg, fi, index);
}

std::string GLGizmoEmboss::get_file_name(const std::string &file_path)
{
    size_t pos_last_delimiter = file_path.find_last_of('\\');
    size_t pos_point          = file_path.find_last_of('.');
    size_t offset             = pos_last_delimiter + 1;
    size_t count              = pos_point - pos_last_delimiter - 1;
    return file_path.substr(offset, count);
}

namespace Slic3r::GUI {

class Priv
{
public:
    Priv() = delete;
    struct EmbossObject
    {
        TriangleMesh      mesh;
        std::string       name;
        TextConfiguration cfg;

        EmbossObject(TriangleMesh &&   mesh,
                     std::string       name,
                     TextConfiguration cfg)
            : mesh(std::move(mesh)), name(name), cfg(cfg)
        {}
    };
    struct EmbossVolume : public EmbossObject
    {
        ModelVolumeType type;
        size_t          object_idx;
        Transform3d     transformation;
        EmbossVolume(TriangleMesh &&   mesh,
                     Transform3d       transformation,
                     std::string       name,
                     TextConfiguration cfg,
                     ModelVolumeType   type,
                     size_t            object_idx)
            : EmbossObject(std::move(mesh), name, cfg)
            , transformation(transformation)
            , type(type)
            , object_idx(object_idx)
        {}
    };
    static void create_emboss_object(EmbossObject &data);
    static void create_emboss_volume(EmbossVolume &data);    
};

} // namespace Slic3r

void GLGizmoEmboss::create_emboss_object(TriangleMesh &&   mesh,
                                         std::string       name,
                                         TextConfiguration cfg)
{
    // Move data to call after is not working
    // data are owen by lambda
    auto data = new Priv::EmbossObject(std::move(mesh), name, cfg);
    wxGetApp().plater()->CallAfter([data]() {
        ScopeGuard sg([data]() { delete data; });
        Priv::create_emboss_object(*data);
    });
}

void GLGizmoEmboss::create_emboss_volume(TriangleMesh &&   mesh,
                                         Transform3d       transformation,
                                         std::string       name,
                                         TextConfiguration cfg,
                                         ModelVolumeType   type,
                                         size_t            object_idx)
{
    // Move data to call after is not working
    // data are owen by lambda
    auto data = new Priv::EmbossVolume(std::move(mesh), transformation, name, cfg, type,
                                       object_idx);
    wxGetApp().plater()->CallAfter([data]() {
        ScopeGuard sg([data]() { delete data; });
        Priv::create_emboss_volume(*data);
    });
}

void Priv::create_emboss_object(EmbossObject &data)
{
    GUI_App &        app      = wxGetApp();
    Plater *         plater   = app.plater();
    ObjectList *     obj_list = app.obj_list();
    GLCanvas3D *     canvas   = plater->canvas3D();
    GLGizmosManager &manager  = canvas->get_gizmos_manager();

    plater->take_snapshot(_L("Add Emboss text object"));
    // Create new object and change selection
    bool center = true;
    obj_list->load_mesh_object(std::move(data.mesh), data.name, center,
                               &data.cfg);

    // new object successfuly added so open gizmo
    assert(manager.get_current_type() != GLGizmosManager::Emboss);
    manager.open_gizmo(GLGizmosManager::Emboss);

    // redraw scene
    canvas->reload_scene(true);

    // Gizmo is not open during time of creation object
    // When cursor move and no one object is selected than Manager::reset_all()
}

void Priv::create_emboss_volume(EmbossVolume &data)
{
    GUI_App &   app      = wxGetApp();
    Plater *    plater   = app.plater();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = plater->canvas3D();

    size_t          object_idx = data.object_idx;
    ModelVolumeType type       = data.type;

    // create new volume inside of object
    Model &model = plater->model();
    if (model.objects.size() <= object_idx) return;
    ModelObject *obj    = model.objects[object_idx];
    ModelVolume *volume = obj->add_volume(std::move(data.mesh));

    // set a default extruder value, since user can't add it manually
    volume->config.set_key_value("extruder", new ConfigOptionInt(0));

    // do not allow model reload from disk
    volume->source.is_from_builtin_objects = true;
    volume->set_type(type);
    volume->name               = data.name;
    volume->text_configuration = data.cfg;
    volume->set_transformation(data.transformation);

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) {
        return vol == volume;
    };
    wxDataViewItemArray sel =
        obj_list->reorder_volumes_and_get_selection((int) object_idx,
                                                    add_to_selection);
    if (!sel.IsEmpty()) obj_list->select_item(sel.front());

    // update printable state on canvas
    if (type == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object(object_idx);

    obj_list->selection_changed();

    // WHY selection_changed set manipulation to world ???
    // so I set it back to local --> RotationGizmo need it
    ObjectManipulation *manipul = wxGetApp().obj_manipul();
    manipul->set_coordinates_type(ECoordinatesType::Local);


    // redraw scene
    canvas->reload_scene(true);
}

void GLGizmoEmboss::update_emboss_volume(TriangleMesh &&   mesh,
                                         const std::string& name,
                                         const TextConfiguration& cfg,
                                         ModelVolume *     volume)
{
    // for sure that some object is created from shape
    if (mesh.its.indices.empty()) return;

    GUI_App &        app      = wxGetApp(); // may be move to input
    Plater *         plater   = app.plater();
    ObjectList *     obj_list = app.obj_list();
    GLCanvas3D *     canvas   = plater->canvas3D();
    GLGizmosManager &manager  = canvas->get_gizmos_manager();

    // Check emboss gizmo is still open
    if (manager.get_current_type() != GLGizmosManager::Emboss) return;

    plater->take_snapshot(_L("Emboss text") + ": " + name);

    // find volume by object id - NOT WORK
    // -> edit text change volume id so could apper not found volume
    // ModelVolume *volume = nullptr;
    // Model &model = plater->model();
    // for (auto obj : model.objects)
    //    for (auto vol : obj->volumes)
    //        if (vol->id() == volume_id) {
    //            volume = vol;
    //            break;
    //        }
    // if (volume == nullptr) return;
    assert(volume != nullptr);

    // update volume
    volume->set_mesh(std::move(mesh));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();
    volume->name               = name;
    volume->text_configuration = cfg;

    // update volume in right panel( volume / object name)
    const Selection &selection = canvas->get_selection();
    const GLVolume * gl_volume = selection.get_volume(
        *selection.get_volume_idxs().begin());
    int object_idx = gl_volume->object_idx();
    int volume_idx = gl_volume->volume_idx();
    obj_list->update_name_in_list(object_idx, volume_idx);

    // update printable state on canvas
    if (volume->type() == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object(
            (size_t) object_idx);

    // redraw scene
    canvas->reload_scene(true);
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
