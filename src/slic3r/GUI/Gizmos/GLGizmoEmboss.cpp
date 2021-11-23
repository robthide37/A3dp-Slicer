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

// TODO: remove include
#include "libslic3r/SVG.hpp"      // debug store
#include "libslic3r/Geometry.hpp" // covex hull 2d

#include "libslic3r/Model.hpp"
#include "libslic3r/ClipperUtils.hpp" // union_ex
#include "libslic3r/AppConfig.hpp"    // store/load font list
#include "libslic3r/MapUtils.hpp"

#include "imgui/imgui_stdlib.h" // using std::string for inputs
#include "nanosvg/nanosvg.h"    // load SVG file

#include <wx/font.h>
#include <wx/fontutil.h>
#include <wx/fontdlg.h>

#include <GL/glew.h>

#ifdef _WIN32
// no specific include - already inside of libslic3r/Emboss.hpp
#elif defined(__APPLE__)
#include <wx/uri.h>
#include <CoreText/CTFont.h>
#include <wx/osx/core/cfdictionary.h>
#elif defined(__linux__)
// Font Config must exist
#include <wx/filename.h>
#include <fontconfig/fontconfig.h>
#endif

// uncomment for easier debug
//#define ALLOW_DEBUG_MODE

namespace Slic3r {
class WxFontUtils
{
public:
    WxFontUtils() = delete;

    // os specific load of wxFont
    static std::optional<Slic3r::Emboss::Font> load_font(const wxFont &font);

    static FontItem::Type get_actual_type();
    static FontItem get_font_item(const wxFont &font);

    // load font used by Operating system as default GUI
    static FontItem get_os_font();
    static std::string get_human_readable_name(const wxFont &font);

    // serialize / deserialize font
    static std::string store_wxFont(const wxFont &font);
    static std::optional<wxFont> load_wxFont(const std::string &font_descriptor);

    // Try to create similar font, loaded from 3mf from different Computer
    static std::optional<wxFont> create_wxFont(const FontItem &fi, const FontProp &fp);
    // update font property by wxFont
    static void update_property(FontProp& font_prop, const wxFont& font);

    // map to convert wxFont type to string and vice versa
    static const std::map<wxFontFamily, std::string> from_family;
    static const std::map<std::string, wxFontFamily> to_family;

    static const std::map<wxFontStyle, std::string> from_style;
    static const std::map<std::string, wxFontStyle>   to_style;

    static const std::map<wxFontWeight, std::string> from_weight;
    static const std::map<std::string, wxFontWeight> to_weight;
};

class NSVGUtils
{
public:
    NSVGUtils() = delete;

    // inspired by nanosvgrast.h function nsvgRasterize->nsvg__flattenShape
    static void flatten_cubic_bez(Polygon &polygon, float tessTol,
        Vec2f p1, Vec2f p2, Vec2f p3,Vec2f p4, int level);
    // convert svg image to ExPolygons
    static ExPolygons to_ExPolygons(NSVGimage *image, float tessTol = 10., int max_level = 10);
};

#ifdef __linux__
class FontConfigHelp
{
public:
    FontConfigHelp() { FcInit(); }
    ~FontConfigHelp() { FcFini(); }
    // inspired by wxpdfdoc - https://github.com/utelle/wxpdfdoc/blob/5bdcdb9953327d06dc50ec312685ccd9bc8400e0/src/pdffontmanager.cpp
    std::string get_font_path(const wxFont &font) {   
        FcConfig *fc = FcInitLoadConfigAndFonts();
        if (fc == nullptr) return "";
        ScopeGuard sg_fc([fc]() { FcConfigDestroy(fc); });

        wxString fontDesc = font.GetNativeFontInfoUserDesc();
        wxString faceName = font.GetFaceName();
        const wxScopedCharBuffer faceNameBuffer = faceName.ToUTF8();
        const char* fontFamily = faceNameBuffer;

        // Check font slant
        int slant = FC_SLANT_ROMAN;
        if (fontDesc.Find(wxS("Oblique")) != wxNOT_FOUND)     slant = FC_SLANT_OBLIQUE;
        else if (fontDesc.Find(wxS("Italic")) != wxNOT_FOUND) slant = FC_SLANT_ITALIC;

        // Check font weight
        int weight = FC_WEIGHT_NORMAL;
        if (fontDesc.Find(wxS("Book")) != wxNOT_FOUND)             weight = FC_WEIGHT_BOOK;
        else if (fontDesc.Find(wxS("Medium")) != wxNOT_FOUND)      weight = FC_WEIGHT_MEDIUM;
#ifdef FC_WEIGHT_ULTRALIGHT
        else if (fontDesc.Find(wxS("Ultra-Light")) != wxNOT_FOUND) weight = FC_WEIGHT_ULTRALIGHT;
#endif
        else if (fontDesc.Find(wxS("Light")) != wxNOT_FOUND)       weight = FC_WEIGHT_LIGHT;
        else if (fontDesc.Find(wxS("Semi-Bold")) != wxNOT_FOUND)   weight = FC_WEIGHT_DEMIBOLD;
#ifdef FC_WEIGHT_ULTRABOLD
        else if (fontDesc.Find(wxS("Ultra-Bold")) != wxNOT_FOUND)  weight = FC_WEIGHT_ULTRABOLD;
#endif
        else if (fontDesc.Find(wxS("Bold")) != wxNOT_FOUND)        weight = FC_WEIGHT_BOLD;
        else if (fontDesc.Find(wxS("Heavy")) != wxNOT_FOUND)       weight = FC_WEIGHT_BLACK;

        // Check font width
        int width = FC_WIDTH_NORMAL;
        if (fontDesc.Find(wxS("Ultra-Condensed")) != wxNOT_FOUND)      width = FC_WIDTH_ULTRACONDENSED;
        else if (fontDesc.Find(wxS("Extra-Condensed")) != wxNOT_FOUND) width = FC_WIDTH_EXTRACONDENSED;
        else if (fontDesc.Find(wxS("Semi-Condensed")) != wxNOT_FOUND)  width = FC_WIDTH_SEMICONDENSED;
        else if (fontDesc.Find(wxS("Condensed")) != wxNOT_FOUND)       width = FC_WIDTH_CONDENSED;
        else if (fontDesc.Find(wxS("Ultra-Expanded")) != wxNOT_FOUND)  width = FC_WIDTH_ULTRAEXPANDED;
        else if (fontDesc.Find(wxS("Extra-Expanded")) != wxNOT_FOUND)  width = FC_WIDTH_EXTRAEXPANDED;
        else if (fontDesc.Find(wxS("Semi-Expanded")) != wxNOT_FOUND)   width = FC_WIDTH_SEMIEXPANDED;
        else if (fontDesc.Find(wxS("Expanded")) != wxNOT_FOUND)        width = FC_WIDTH_EXPANDED;

        FcResult res;
        FcPattern* matchPattern = FcPatternBuild(NULL, FC_FAMILY, FcTypeString, (FcChar8*) fontFamily, NULL);
        ScopeGuard sg_mp([matchPattern]() { FcPatternDestroy(matchPattern); });
        
        FcPatternAddInteger(matchPattern, FC_SLANT, slant);        
        FcPatternAddInteger(matchPattern, FC_WEIGHT, weight);        
        FcPatternAddInteger(matchPattern, FC_WIDTH, width);

        FcConfigSubstitute (NULL, matchPattern, FcMatchPattern);
        FcDefaultSubstitute (matchPattern);

        FcPattern* resultPattern = FcFontMatch (NULL, matchPattern, &res);
        if (resultPattern == nullptr) return "";
        ScopeGuard sg_rp([resultPattern]() { FcPatternDestroy(resultPattern); });
            

        FcChar8 *fileName;
        if (FcPatternGetString(resultPattern, FC_FILE, 0, &fileName) != FcResultMatch) return "";
        wxString fontFileName = wxString::FromUTF8((char *) fileName);

        if (fontFileName.IsEmpty()) return "";

        // find full file path
        wxFileName myFileName(fontFileName);
        if (!myFileName.IsOk()) return "";

        if (myFileName.IsRelative()) {
            // Check whether the file is relative to the current working directory
            if (!(myFileName.MakeAbsolute() && myFileName.FileExists())) {
                return "";
                // File not found, search in given search paths
                //wxString foundFileName = m_searchPaths.FindAbsoluteValidPath(fileName);
                //if (!foundFileName.IsEmpty()) {
                //    myFileName.Assign(foundFileName);
                //}
            }
        }

        if (!myFileName.FileExists() || !myFileName.IsFileReadable()) return "";
            
        // File exists and is accessible
        wxString fullFileName = myFileName.GetFullPath();
        return std::string(fullFileName.c_str());
    }
};

#endif // __linux__

} // namespace Slic3r

using namespace Slic3r;
using namespace Slic3r::GUI;

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -2)
    , m_font_selected(0)
    , m_font(nullptr)
    , m_volume(nullptr)
    , m_exist_notification(false)
    , m_is_initialized(false) // initialize on first opening gizmo
    , m_job(std::make_unique<EmbossJob>())
{
    // TODO: add suggestion to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web  
    
}

GLGizmoEmboss::~GLGizmoEmboss() {}

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

    // TODO: fix width
    ImVec2 windows_size = m_gui_cfg->draw_advanced ?
                              m_gui_cfg->minimal_window_size_with_advance :
                              m_gui_cfg->minimal_window_size;
    ImVec2 offset = ImGuiWrapper::suggest_location(windows_size, hull);
    m_gui_cfg->offset   = offset;
    return;

    Polygon rect({Point(offset.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y),
                  Point(offset.x + windows_size.x, offset.y + windows_size.y),
                  Point(offset.x, offset.y + windows_size.y)});
    ImGuiWrapper::draw(hull);
    ImGuiWrapper::draw(rect);
    // m_gui_cfg->offset = offset;
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

void GLGizmoEmboss::set_volume_type(ModelVolumeType volume_type)
{
    if (m_volume == nullptr) return;
    m_volume->set_type(volume_type);
    m_parent.reload_scene(true);
}

bool GLGizmoEmboss::on_init()
{
    // m_grabbers.emplace_back();
    m_shortcut_key = WXK_CONTROL_T;
    return true;
}

std::string GLGizmoEmboss::on_get_name() const { return _u8L("Emboss"); }

void GLGizmoEmboss::on_render() {}
void GLGizmoEmboss::on_render_for_picking() {}

void GLGizmoEmboss::on_render_input_window(float x, float y, float bottom_limit)
{
    check_selection();


    ImVec2 min_window_size = m_gui_cfg->draw_advanced ?
                                 m_gui_cfg->minimal_window_size_with_advance :
                                 m_gui_cfg->minimal_window_size;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, min_window_size);
    // ImGui::SetNextWindowSize(ImVec2(0, min_window_size.y),
    // ImGuiCond_::ImGuiCond_Always);

#ifdef ALLOW_DEBUG_MODE
    // draw suggested position of window
    draw_fine_position(m_parent.get_selection());
#endif // ALLOW_DEBUG_MODE

    // check if is set window offset
    if (m_gui_cfg->offset.has_value()) {
        ImGui::SetNextWindowPos(*m_gui_cfg->offset, ImGuiCond_Always);
        // clear request on offset
        m_gui_cfg->offset = {};
    }

    int flag = // ImGuiWindowFlags_AlwaysAutoResize |
               // ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoCollapse;
    m_imgui->begin(on_get_name(), flag);

    ImGui::GetWindowSize();
    draw_window();
    m_imgui->end();
    ImGui::PopStyleVar(); // WindowMinSize
}

void GLGizmoEmboss::on_set_state()
{
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
        set_fine_position();

        // When add Text on empty plate, Create new object with volume
        if (m_parent.get_selection().is_empty()) {
            if (!create_default_model_object())
                GLGizmoBase::m_state = GLGizmoBase::Off;
            return;
        }

        // Try set selected volume
        if (!load_configuration(get_selected_volume())) {
            // No volume with text selected, create new one
            set_default_configuration();
            process();
        }

        // when open by hyperlink it needs to show up
        m_parent.reload_scene(true);
    }
}

void GLGizmoEmboss::initialize()
{
    if (m_is_initialized) return;
    m_is_initialized = true;
    m_gui_cfg.emplace(GuiCfg());
    float space = ImGui::GetTextLineHeightWithSpacing() -
                  ImGui::GetTextLineHeight();
    m_gui_cfg->max_font_name_width = ImGui::CalcTextSize("Maximal font name").x;
    m_gui_cfg->icon_width = ImGui::GetTextLineHeight();
    m_gui_cfg->icon_width_with_spacing = m_gui_cfg->icon_width + space;

    float scroll_width          = m_gui_cfg->icon_width_with_spacing; // fix
    m_gui_cfg->combo_font_width = m_gui_cfg->max_font_name_width +
                                  2 * m_gui_cfg->icon_width_with_spacing +
                                  scroll_width;

    m_gui_cfg->rename_pos_x = m_gui_cfg->max_font_name_width + space;
    m_gui_cfg->delete_pos_x = m_gui_cfg->rename_pos_x +
                              m_gui_cfg->icon_width_with_spacing;

    m_gui_cfg->text_size = ImVec2(-FLT_MIN, ImGui::GetTextLineHeight() *
                                                m_gui_cfg->count_line_of_text);

    m_gui_cfg->advanced_input_width = ImGui::CalcTextSize("M").x * 6;

    // TODO: What to do when icon was NOT loaded?
    bool success = init_icons();

    load_font_list();

    m_font_selected = 0;

    bool is_font_loaded = load_font();
    while (!is_font_loaded && !m_font_list.empty()) {
        // can't load so erase it from list
        m_font_list.erase(m_font_list.begin() + m_font_selected);
        m_font_selected = 0; // select first
        is_font_loaded  = load_font();
    }
    set_default_configuration();
}

FontList GLGizmoEmboss::create_default_font_list() {
    return {
        WxFontUtils::get_font_item(wxFont(5, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL))
        , WxFontUtils::get_font_item(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD))
        , WxFontUtils::get_os_font()
    };
}

void GLGizmoEmboss::set_default_configuration()
{
    m_text      = _u8L("Embossed text");
    m_font_prop = FontProp();
    load_font(); // reload actual font - because of font size
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
    if (load_configuration(vol)) {
        // successfull load volume for editing
        return;
    }

    // behave like adding new text
    m_volume = nullptr;
    set_default_configuration();
}

ModelVolume *GLGizmoEmboss::get_selected_volume()
{
    return get_selected_volume(m_parent.get_selection(),
                               wxGetApp().plater()->model().objects);
}

ModelVolume *GLGizmoEmboss::get_selected_volume(const Selection &selection,
                                                const ModelObjectPtrs objects)
{
    int object_idx = selection.get_object_idx();
    // is more object selected?
    if (object_idx == -1) return nullptr;

    auto volume_idxs = selection.get_volume_idxs();
    // is more volumes selected?
    if (volume_idxs.size() != 1) return nullptr;
    unsigned int                 vol_id_gl = *volume_idxs.begin();
    const GLVolume *             vol_gl    = selection.get_volume(vol_id_gl);
    const GLVolume::CompositeID &id        = vol_gl->composite_id;

    if (id.object_id < 0 || static_cast<size_t>(id.object_id) >= objects.size())
        return nullptr;
    ModelObject *object = objects[id.object_id];

    if (id.volume_id < 0 || static_cast<size_t>(id.volume_id) >= object->volumes.size()) 
        return nullptr;
    return object->volumes[id.volume_id];
}

bool GLGizmoEmboss::process()
{
    // exist loaded font?
    if (m_font == nullptr) return false;


    EmbossJob::Data data;
    data.font = m_font;
    data.text_configuration = create_configuration();
    data.volume_name = create_volume_name();
    data.volume_ptr  = m_volume;
    data.object_idx = m_parent.get_selection().get_object_idx();
    m_job->restart(data);

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
    if (m_font == nullptr) {
        ImGui::Text(_u8L("Warning: No font is selected. Select correct one.").c_str());
    }
    draw_font_list();
    draw_text_input();

    bool &advanced = m_gui_cfg->draw_advanced;
    if (ImGui::Checkbox(_u8L("Advance").c_str(), &advanced)) {
        ImVec2 window_size =
            advanced ?
                ImVec2(0, m_gui_cfg->minimal_window_size_with_advance.y) :
                m_gui_cfg->minimal_window_size;
        ImGui::SetWindowSize(window_size, ImGuiCond_Always);
    }
    if (advanced) draw_advanced();

    if (ImGui::Button(_u8L("Close").c_str())) close();

    // Option to create text volume when reselecting volumes
    m_imgui->disabled_begin(m_font == nullptr);
    if (m_volume == nullptr) {
        ImGui::SameLine();
        if (ImGui::Button(_u8L("Generate preview").c_str())) process();
    }
    m_imgui->disabled_end();
}

void GLGizmoEmboss::draw_font_list()
{
    const float &         max_width = m_gui_cfg->max_font_name_width;
    std::optional<size_t> rename_index;
    std::string current_name = imgui_trunc(m_font_list[m_font_selected].name,
                                           max_width);
    ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
    if (ImGui::BeginCombo("##font_selector", current_name.c_str())) {
        // first line
        if (ImGui::Button(_u8L("Choose font").c_str())) {
            choose_font_by_wxdialog();
            store_font_list();
            ImGui::CloseCurrentPopup();
        } else if (ImGui::IsItemHovered())
            ImGui::SetTooltip(
                _u8L("Choose from installed font inside dialog.").c_str());
                
#ifdef ALLOW_DEBUG_MODE
        ImGui::SameLine();
        // select font file by file browser
         if (ImGui::Button(_u8L("Add File").c_str())) {
            choose_true_type_file();
            store_font_list();
            ImGui::CloseCurrentPopup();
        } else if (ImGui::IsItemHovered()) ImGui::SetTooltip(_u8L("add file with font(.ttf, .ttc)").c_str());
#endif //  ALLOW_DEBUG_MODE

        ImGui::Separator();

        for (FontItem &f : m_font_list) {
            ImGui::PushID(f.name.c_str());
            std::string name        = imgui_trunc(f.name, max_width);
            size_t      index       = &f - &m_font_list.front();
            bool        is_selected = index == m_font_selected;
            auto        flags =
                ImGuiSelectableFlags_AllowItemOverlap; // allow clic buttons
            if (ImGui::Selectable(name.c_str(), is_selected, flags)) {
                size_t prev_font_selected = m_font_selected;
                m_font_selected           = index;
                if (!load_font()) {
                    m_font_selected = prev_font_selected;
                } else {
                    process();
                }
            } else if (ImGui::IsItemHovered())
                ImGui::SetTooltip(f.name.c_str());

            // draw buttons rename and delete
            ImGui::SameLine();
            ImGui::SetCursorPosX(m_gui_cfg->rename_pos_x);
            if (draw_button(IconType::rename)) rename_index = index;
            ImGui::SameLine();
            ImGui::SetCursorPosX(m_gui_cfg->delete_pos_x);
            if (draw_button(IconType::erase, is_selected)) {
                m_font_list.erase(m_font_list.begin() + index);
                // fix selected index
                if (index < m_font_selected) --m_font_selected;
                store_font_list();
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }

    // rename modal window popup
    const char *  rename_popup_id = "Rename_font";
    static size_t rename_id;
    if (rename_index.has_value() && !ImGui::IsPopupOpen(rename_popup_id)) {
        ImGui::OpenPopup(rename_popup_id);
        rename_id = *rename_index;
    }
    if (ImGui::BeginPopupModal(rename_popup_id)) {
        FontItem &  fi = m_font_list[rename_id];
        std::string rename_popup =
            GUI::format(_u8L("Change font name (%1%): "), fi.name);
        ImGui::Text(rename_popup.c_str());
        ImGui::SetNextItemWidth(m_gui_cfg->combo_font_width);
        if (ImGui::InputText("##font name", &fi.name,
                             ImGuiInputTextFlags_EnterReturnsTrue) ||
            ImGui::Button("ok")) {
            ImGui::CloseCurrentPopup();
            store_font_list();
        }
        ImGui::EndPopup();
    }
}

void GLGizmoEmboss::draw_text_input()
{
    static const ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_AllowTabInput | ImGuiInputTextFlags_AutoSelectAll;

    ImVector<ImFont *> &fonts      = m_imgui_font_atlas.Fonts;
    ImFont *            imgui_font = fonts.empty() ? nullptr : fonts.front();
    bool exist_font = imgui_font != nullptr && imgui_font->IsLoaded();
    if (exist_font) ImGui::PushFont(imgui_font);

    bool exist_change = false;
    float window_height = ImGui::GetWindowHeight();
    float minimal_height = m_gui_cfg->draw_advanced ?
                               m_gui_cfg->minimal_window_size_with_advance.y :
                               m_gui_cfg->minimal_window_size.y;
    float extra_height   = window_height - minimal_height;
    ImVec2 text_size(m_gui_cfg->text_size.x,
                     m_gui_cfg->text_size.y + extra_height);
    if (ImGui::InputTextMultiline("##Text", &m_text, text_size, flags)) {
        process();
        exist_change = true;
    }

    if (exist_font) ImGui::PopFont();

    // imgui_font has to be unused
    if (exist_change) check_imgui_font_range();
}

void GLGizmoEmboss::draw_advanced()
{
    ImGui::SetNextItemWidth(m_gui_cfg->advanced_input_width);
    if (ImGui::InputFloat(_u8L("Size[in mm]").c_str(),
                          &m_font_prop.size_in_mm)) {
        if (m_font_prop.size_in_mm < 0.1) m_font_prop.size_in_mm = 10;
        // store font size into path
        FontItem &fi = m_font_list[m_font_selected];
        if (fi.type == WxFontUtils::get_actual_type()) {
            std::optional<wxFont> wx_font = WxFontUtils::load_wxFont(fi.path);
            if (wx_font.has_value()) {
                wx_font->SetPointSize(m_font_prop.size_in_mm);
                fi.path = WxFontUtils::store_wxFont(*wx_font);
            }
        }
        load_imgui_font();
        if (m_font != nullptr) m_font->cache.clear();
        process();
    }
    ImGui::SetNextItemWidth(m_gui_cfg->advanced_input_width);
    if (ImGui::InputFloat(_u8L("Emboss[in mm]").c_str(), &m_font_prop.emboss))
        process();
    ImGui::SetNextItemWidth(2*m_gui_cfg->advanced_input_width);
    if (ImGui::InputInt(_u8L("CharGap[in font points]").c_str(),
                        &m_font_prop.char_gap))
        process();
    ImGui::SetNextItemWidth(2*m_gui_cfg->advanced_input_width);
    if (ImGui::InputInt(_u8L("LineGap[in font points]").c_str(),
                        &m_font_prop.line_gap))
        process();

    // when more collection add selector
    if (m_font != nullptr && m_font->count > 1) {
        ImGui::SetNextItemWidth(m_gui_cfg->advanced_input_width);
        if (ImGui::BeginCombo(_u8L("Font collection").c_str(),
                              std::to_string(m_font->index).c_str())) {
            for (unsigned int i = 0; i < m_font->count; ++i) {
                ImGui::PushID(1 << 10 + i);
                if (ImGui::Selectable(std::to_string(i).c_str(),
                                      i == m_font->index)) {
                    m_font->index = i;
                    m_font->cache.clear();
                    process();
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
    }

#ifdef ALLOW_DEBUG_MODE
    std::string descriptor = m_font_list[m_font_selected].path;
    ImGui::Text("family = %s", (m_font_prop.family.has_value() ?
                                    m_font_prop.family->c_str() :
                                    " --- "));
    ImGui::Text("face name = %s", (m_font_prop.face_name.has_value() ?
                                       m_font_prop.face_name->c_str() :
                                       " --- "));
    ImGui::Text("style = %s",
                (m_font_prop.style.has_value() ? m_font_prop.style->c_str() :
                                                 " --- "));
    ImGui::Text("weight = %s", (m_font_prop.weight.has_value() ?
                                    m_font_prop.weight->c_str() :
                                    " --- "));
    ImGui::Text("descriptor = %s", descriptor.c_str());
    ImGui::Image(m_imgui_font_atlas.TexID,
                 ImVec2(m_imgui_font_atlas.TexWidth,
                        m_imgui_font_atlas.TexHeight));
#endif // ALLOW_DEBUG_MODE
}

bool GLGizmoEmboss::create_default_model_object()
{
    set_default_configuration();
    // Is created default model?
    if (process()) return true;

    // can't create object,
    // e.g. selected font don't have letter for "Emboss text"

    // try select another font
    for (size_t font_index = 0; font_index < m_font_list.size();
         ++font_index) {
        if (!load_font(font_index)) continue;
        // Is fixed by change to font from font list?
        if (process()) return true;
    }

    // try add system font and use it
    size_t font_index = m_font_list.size();
    m_font_list.push_back(WxFontUtils::get_os_font());
    if (!load_font(font_index)) {
        // TODO: Solve wrong os font !!!
        return false;
    }
    // Is fixed by change font to os font?
    if (process()) return true;

    // try change default text
    // Do NOT translate it!
    m_text = u8"text";
    // Is fixed by change translation of Emboss text?
    if (process()) return true;

    // Bad os font can't process "text" with os default font
    return false;
}

bool GLGizmoEmboss::load_font(size_t font_index)
{
    std::swap(font_index, m_font_selected);
    bool is_loaded = load_font();
    if (!is_loaded) std::swap(font_index, m_font_selected);
    return is_loaded;
}

bool GLGizmoEmboss::load_font()
{
    if (m_font_selected >= m_font_list.size()) return false;
    FontItem &fi = m_font_list[m_font_selected];
    if (fi.type == FontItem::Type::file_path) {
        // fill font name after load from .3mf
        if (fi.name.empty())
            fi.name = Slic3r::GUI::GLGizmoEmboss::get_file_name(fi.path);
        std::optional<Emboss::Font> font_opt = Emboss::load_font(
            fi.path.c_str());
        if (!font_opt.has_value()) return false;
        // TODO: fix copy of font data
        m_font = std::make_shared<Emboss::Font>(std::move(*font_opt));
        load_imgui_font();
        return true;
    }
    if (fi.type != WxFontUtils::get_actual_type()) return false;
    std::optional<wxFont> wx_font = WxFontUtils::load_wxFont(fi.path);
    if (!wx_font.has_value()) return false;

    // fill font name after load from .3mf
    if (fi.name.empty())
        fi.name = WxFontUtils::get_human_readable_name(*wx_font);
    return load_font(*wx_font);
}

bool GLGizmoEmboss::load_font(const wxFont &font)
{
    auto font_opt = WxFontUtils::load_font(font);
    if (!font_opt.has_value()) return false;
    // TODO: fix copy of font data
    m_font = std::make_shared<Emboss::Font>(std::move(*font_opt));
    WxFontUtils::update_property(m_font_prop, font);
    m_font_prop.emboss = m_font_prop.size_in_mm / 2.f;
    load_imgui_font();
    return true;
}

void GLGizmoEmboss::check_imgui_font_range()
{
    const char *text = m_text.c_str();

    const ImFont *font = m_imgui_font_atlas.Fonts.front();
    if (!font->IsLoaded()) {
        // when create font no one letter in text was inside font
        // check text again
        load_imgui_font();
        return;
    }
    if (font->ConfigData == nullptr) return;
    const ImWchar *ranges       = font->ConfigData->GlyphRanges;
    auto           is_in_ranges = [ranges](unsigned int letter) -> bool {
        for (const ImWchar *range = ranges; range[0] && range[1]; range += 2) {
            ImWchar from = range[0];
            ImWchar to   = range[1];
            if (from <= letter && letter <= to) return true;
            if (letter < to) return false; // ranges should be sorted
        }
        return false;
    };

    bool exist_unknown = false;
    while (*text) {
        unsigned int c     = 0;
        int          c_len = ImTextCharFromUtf8(&c, text, NULL);
        text += c_len;
        if (c_len == 0) break;
        if (!is_in_ranges(c)) {
            exist_unknown = true;
            break;
        }
    }
    if (exist_unknown) load_imgui_font();
}

void GLGizmoEmboss::load_imgui_font()
{
    if (m_font == nullptr) return;

    ImFontGlyphRangesBuilder builder;
    builder.AddRanges(m_imgui->get_glyph_ranges());
    builder.AddText(m_text.c_str());

    m_imgui_font_ranges.clear();

    builder.BuildRanges(&m_imgui_font_ranges);
    int font_size = static_cast<int>(
        std::round(std::abs(m_font_prop.size_in_mm / 0.3528)));
    if (font_size < m_gui_cfg->min_imgui_font_size)
        font_size = m_gui_cfg->min_imgui_font_size;
    if (font_size > m_gui_cfg->max_imgui_font_size)
        font_size = m_gui_cfg->max_imgui_font_size;

    ImFontConfig font_config;
    font_config.FontDataOwnedByAtlas = false;
    m_imgui_font_atlas.Flags |= ImFontAtlasFlags_NoMouseCursors |
                                ImFontAtlasFlags_NoPowerOfTwoHeight;
    m_imgui_font_atlas.Clear();
    m_imgui_font_atlas.AddFontFromMemoryTTF((void *) m_font->buffer.data(),
                                            m_font->buffer.size(), font_size,
                                            &font_config,
                                            m_imgui_font_ranges.Data);

    unsigned char *pixels;
    int            width, height;
    m_imgui_font_atlas.GetTexDataAsAlpha8(&pixels, &width, &height);

    // Upload texture to graphics system
    GLint last_texture;
    glsafe(::glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture));
    ScopeGuard sg([last_texture]() {
        glsafe(::glBindTexture(GL_TEXTURE_2D, last_texture));
    });

    GLuint font_texture;
    glsafe(::glGenTextures(1, &font_texture));
    glsafe(::glBindTexture(GL_TEXTURE_2D, font_texture));
    glsafe(::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
    glsafe(::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
    glsafe(::glPixelStorei(GL_UNPACK_ROW_LENGTH, 0));
    glsafe(::glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, width, height, 0,
                          GL_ALPHA, GL_UNSIGNED_BYTE, pixels));

    // Store our identifier
    m_imgui_font_atlas.TexID = (ImTextureID) (intptr_t) font_texture;
}

bool GLGizmoEmboss::choose_font_by_wxdialog()
{
    wxFontData data;
    data.EnableEffects(false);
    data.RestrictSelection(wxFONTRESTRICT_SCALABLE);
    // set previous selected font
    FontItem &selected_font_item = m_font_list[m_font_selected];
    if (selected_font_item.type == WxFontUtils::get_actual_type()) {
        std::optional<wxFont> selected_font = WxFontUtils::load_wxFont(
            selected_font_item.path);
        if (selected_font.has_value()) data.SetInitialFont(*selected_font);
    }

    wxFontDialog font_dialog(wxGetApp().mainframe, data);
    if (font_dialog.ShowModal() != wxID_OK) return false;

    data                = font_dialog.GetFontData();
    wxFont   font       = data.GetChosenFont();
    size_t   font_index = m_font_list.size();
    FontItem font_item  = WxFontUtils::get_font_item(font);
    m_font_list.emplace_back(font_item);
    FontProp old_font_prop = m_font_prop; // copy

    // Check that deserialization NOT influence font
    // false - use direct selected wxFont in dialog
    // true - use font item (serialize and deserialize wxFont)
    bool use_deserialized_font = false;
    if (!use_deserialized_font) m_font_selected = font_index;
    // Try load and use new added font
    if ((!use_deserialized_font && !load_font(font)) ||
        (use_deserialized_font && !load_font(font_index)) || !process()) {
        // reverse index for font_selected
        std::swap(font_index, m_font_selected); // when not process
        // remove form font list
        m_font_list.pop_back();
        // reverse property
        m_font_prop      = old_font_prop; // when not process
        wxString message = GUI::format_wxstr(
            _L("Font '%1%' can't be used. Please select another."),
            font_item.name);
        wxString      title = _L("Selected font is NOT True-type.");
        MessageDialog not_loaded_font_message(nullptr, message, title, wxOK);
        not_loaded_font_message.ShowModal();
        return choose_font_by_wxdialog();
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
        m_font_list.emplace_back(name, path, FontItem::Type::file_path);

        // set first valid added font as active
        if (!font_loaded) {
            if (!load_font(m_font_list.size() - 1))
                m_font_list.pop_back();
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
    float scale   = m_font_prop.size_in_mm / std::max(bb.max.x(), bb.max.y());
    auto  project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(m_font_prop.emboss / scale), scale);
    indexed_triangle_set its = Emboss::polygons2model(polys, *project);
    // test store:
    // for (auto &poly : polys) poly.scale(1e5);
    // SVG svg("converted.svg", BoundingBox(polys.front().contour.points));
    // svg.draw(polys);
    //return add_volume(name, its);
}

TextConfiguration GLGizmoEmboss::create_configuration()
{
    return TextConfiguration(m_font_list[m_font_selected], m_font_prop,
                             m_text);
}

bool GLGizmoEmboss::load_configuration(ModelVolume *volume)
{
    if (volume == nullptr) return false;
    if (!volume->text_configuration.has_value()) return false;

    TextConfiguration &configuration = *volume->text_configuration;
    FontItem &         c_font_item   = configuration.font_item;

    // try to find font in local font list
    auto is_config = [&c_font_item](const FontItem &font_item) -> bool {
        return font_item.path == c_font_item.path;
    };
    auto it = std::find_if(m_font_list.begin(), m_font_list.end(), is_config);

    size_t prev_font_selected = m_font_selected;

    if (it == m_font_list.end()) {
        // font is not in list
        // add font to list
        m_font_selected = m_font_list.size();
        m_font_list.emplace_back(c_font_item);
    } else {
        // font is found in list
        m_font_selected = it - m_font_list.begin();
    }

    m_font_prop = configuration.font_prop;
    m_text      = configuration.text;
    m_volume    = volume;

    if (!load_font()) {
        // create similar font
        auto wx_font = WxFontUtils::create_wxFont(c_font_item,
                                                  configuration.font_prop);
        if (wx_font.has_value()) {
            // fix not loadable font item
            m_font_list[m_font_selected] = WxFontUtils::get_font_item(
                *wx_font);
            if (!load_font(*wx_font)) return false;
        } else {
            // can't create similar font use previous
            m_font_list.erase(m_font_list.begin() + m_font_selected);
            m_font_selected = prev_font_selected;
            if (!load_font()) return false;
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

    const auto &origin_family = tc.font_prop.face_name;
    const auto &actual_family = m_font_prop.face_name;
    const auto &fi            = m_font_list[m_font_selected];

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
    std::string path = resources_dir() + "/icons/white/";

    // icon order has to match the enum IconType
    std::vector<std::string> filenames = {path + "wrench.svg",
                                          path + "delete.svg"};

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
    unsigned int icons_texture_id = m_icons_texture.get_id();
    int          tex_width        = m_icons_texture.get_width();
    int          tex_height       = m_icons_texture.get_height();

    // is icon loaded
    if ((icons_texture_id == 0) || (tex_width <= 1) || (tex_height <= 1))
        return;
    ImTextureID tex_id = (void *) (intptr_t) (GLuint) icons_texture_id;
    // ImVec2 image_size(tex_width, tex_height);

    size_t count_icons  = 2; // wrench | delete
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
    float line_spacing = ImGui::GetTextLineHeightWithSpacing() -
                         ImGui::GetTextLineHeight();
    float cursor_pos_y = ImGui::GetCursorPosY();
    ImGui::SetCursorPosY(cursor_pos_y - line_spacing / 2);
    ScopeGuard sg([cursor_pos_y]() {
        ImGui::SetCursorPosY(cursor_pos_y);
        ImGui::NewLine();
    });

    if (disable) {
        draw_icon(icon, IconState::disabled);
        if (ImGui::IsItemHovered() && icon == IconType::erase)
            ImGui::SetTooltip(_u8L("Active font can't be removed").c_str());
        return false;
    }

    float cursor_x = ImGui::GetCursorPosX();

    draw_icon(icon, IconState::activable);
    if (ImGui::IsItemClicked()) return true;
    if (ImGui::IsItemHovered()) {
        switch (icon) {
        case IconType::rename:
            ImGui::SetTooltip(_u8L("rename").c_str());
            break;
        case IconType::erase:
            ImGui::SetTooltip(_u8L("delete").c_str());
            break;
        default: break;
        }
        // redraw image over previous
        ImGui::SameLine();
        ImGui::SetCursorPosX(cursor_x);

        draw_icon(icon, IconState::hovered);
        if (ImGui::IsItemClicked()) return true;
    }
    return false;
}

void GLGizmoEmboss::load_font_list()
{
    const AppConfig *cfg          = wxGetApp().app_config;
    unsigned         index        = 1;
    std::string      section_name = get_app_config_font_section(index++);
    while (cfg->has_section(section_name)) {
        std::optional<FontItem> fi = get_font_item(
            cfg->get_section(section_name));
        if (fi.has_value()) m_font_list.emplace_back(*fi);
        section_name = get_app_config_font_section(index++);
    }
    if (m_font_list.empty()) m_font_list = create_default_font_list();
}

void GLGizmoEmboss::store_font_list()
{
    AppConfig *cfg   = wxGetApp().app_config;
    unsigned   index = 1;
    for (const FontItem &fi : m_font_list) {
        // skip file paths + fonts from other OS
        if (fi.type != WxFontUtils::get_actual_type()) continue;
        set_font_item(*cfg, fi, index++);
    }

    // remove rest of font sections
    std::string section_name = get_app_config_font_section(index);
    while (cfg->has_section(section_name)) {
        cfg->clear_section(section_name);
        section_name = get_app_config_font_section(++index);
    }
}

const std::string GLGizmoEmboss::APP_CONFIG_FONT_NAME       = "name";
const std::string GLGizmoEmboss::APP_CONFIG_FONT_DESCRIPTOR = "descriptor";

std::string GLGizmoEmboss::get_app_config_font_section(unsigned index)
{
    return AppConfig::SECTION_FONT + ':' + std::to_string(index);
}

std::optional<FontItem> GLGizmoEmboss::get_font_item(
    const std::map<std::string, std::string> &app_cfg_section)
{
    auto path_it = app_cfg_section.find(APP_CONFIG_FONT_DESCRIPTOR);
    if (path_it == app_cfg_section.end()) return {};
    const std::string &path = path_it->second;

    auto name_it = app_cfg_section.find(APP_CONFIG_FONT_NAME);
    static const std::string default_name = "font_name";
    const std::string &      name = (name_it == app_cfg_section.end()) ?
                                        default_name :
                                        name_it->second;
    return FontItem(name, path, WxFontUtils::get_actual_type());
}

void GLGizmoEmboss::set_font_item(AppConfig &     cfg,
                                  const FontItem &fi,
                                  unsigned        index)
{
    std::string section_name = get_app_config_font_section(index);
    cfg.clear_section(section_name);
    cfg.set(section_name, APP_CONFIG_FONT_NAME, fi.name);
    cfg.set(section_name, APP_CONFIG_FONT_DESCRIPTOR, fi.path);
}

std::string GLGizmoEmboss::get_file_name(const std::string &file_path)
{
    size_t pos_last_delimiter = file_path.find_last_of('\\');
    size_t pos_point          = file_path.find_last_of('.');
    size_t offset             = pos_last_delimiter + 1;
    size_t count              = pos_point - pos_last_delimiter - 1;
    return file_path.substr(offset, count);
}

std::string GLGizmoEmboss::imgui_trunc(const std::string &text, float width)
{
    static const char *tail       = " ..";
    float              tail_width = ImGui::CalcTextSize(tail).x;
    float              text_width = ImGui::CalcTextSize(text.c_str()).x;
    if (text_width < width) return text;
    float    letter_width  = ImGui::CalcTextSize("n").x;
    float    allowed_width = width - tail_width;
    unsigned count_letter  = static_cast<unsigned>(allowed_width /
                                                  letter_width);
    text_width = ImGui::CalcTextSize(text.substr(0, count_letter).c_str()).x;
    if (text_width < allowed_width) {
        // increase letter count
        do {
            ++count_letter;
            text_width =
                ImGui::CalcTextSize(text.substr(0, count_letter).c_str()).x;
        } while (text_width < allowed_width);
        --count_letter;
    } else {
        // decrease letter count
        do {
            --count_letter;
            text_width =
                ImGui::CalcTextSize(text.substr(0, count_letter).c_str()).x;
        } while (text_width > allowed_width);
    }
    return text.substr(0, count_letter) + tail;
}

/// <summary>
/// WxFontUtils - Start definition
/// </summary>

std::optional<Emboss::Font> WxFontUtils::load_font(const wxFont &font)
{
    if (!font.IsOk()) return {};
#ifdef _WIN32
    return Emboss::load_font(font.GetHFONT());
#elif defined(__APPLE__)
    // use file path
    const wxNativeFontInfo *info = font.GetNativeFontInfo();
    if (info == nullptr) return {};
    CTFontDescriptorRef descriptor = info->GetCTFontDescriptor();
    CFURLRef            typeref    = (CFURLRef)
        CTFontDescriptorCopyAttribute(descriptor, kCTFontURLAttribute);
    CFStringRef url = CFURLGetString(typeref);
    if (url == NULL) return {};
    wxString file_uri;
    wxCFTypeRef(url).GetValue(file_uri);
    std::string file_path(wxURI::Unescape(file_uri).c_str());
    size_t      start = std::string("file://").size();
    if (file_path.empty() || file_path.size() <= start) return {};
    file_path = file_path.substr(start, file_path.size() - start);
    return Emboss::load_font(file_path.c_str());
#elif defined(__linux__)
    static FontConfigHelp help;
    std::string           font_path = help.get_font_path(font);
    if (font_path.empty()) return {};
    return Emboss::load_font(font_path.c_str());
#else
    // HERE is place to add implementation for another platform
    // to convert wxFont to font data as windows or font file path as linux
    return {};
#endif
}

FontItem::Type WxFontUtils::get_actual_type()
{
#ifdef _WIN32
    return FontItem::Type::wx_win_font_descr;
#elif defined(__APPLE__)
    return FontItem::Type::wx_mac_font_descr;
#elif defined(__linux__)
    return FontItem::Type::wx_mac_font_descr;
#else
    return FontItem::Type::undefined;
#endif
}

FontItem WxFontUtils::get_font_item(const wxFont &font)
{
    std::string    name     = get_human_readable_name(font);
    std::string    fontDesc = store_wxFont(font);
    FontItem::Type type     = get_actual_type();
    // wxFont         f        = font; // copy
    return FontItem(name, fontDesc, type);
}

FontItem WxFontUtils::get_os_font()
{
    wxSystemFont system_font = wxSYS_DEFAULT_GUI_FONT;
    wxFont       font        = wxSystemSettings::GetFont(system_font);
    FontItem     fi          = get_font_item(font);
    fi.name += std::string(" (" + _u8L("OS default") + ")");
    return get_font_item(font);
}

std::string WxFontUtils::get_human_readable_name(const wxFont &font)
{
    if (!font.IsOk()) return "Font is NOT ok.";
    // Face name is optional in wxFont
    if (!font.GetFaceName().empty()) {
        return std::string(font.GetFaceName().c_str());
    } else {
        return std::string((font.GetFamilyString() + " " +
                            font.GetStyleString() + " " +
                            font.GetWeightString())
                               .c_str());
    }
}

std::string WxFontUtils::store_wxFont(const wxFont &font)
{
    // wxString os = wxPlatformInfo::Get().GetOperatingSystemIdName();
    wxString font_descriptor = font.GetNativeFontInfoDesc();
    return std::string(font_descriptor.c_str());
}

std::optional<wxFont> WxFontUtils::load_wxFont(
    const std::string &font_descriptor)
{
    wxString font_descriptor_wx(font_descriptor);
    wxFont   wx_font(font_descriptor_wx);
    if (!wx_font.IsOk()) return {};
    return wx_font;
}

const std::map<wxFontFamily, std::string> WxFontUtils::from_family(
    {{wxFONTFAMILY_DEFAULT, "default"},
     {wxFONTFAMILY_DECORATIVE, "decorative"},
     {wxFONTFAMILY_ROMAN, "roman"},
     {wxFONTFAMILY_SCRIPT, "script"},
     {wxFONTFAMILY_SWISS, "swiss"},
     {wxFONTFAMILY_MODERN, "modern"},
     {wxFONTFAMILY_TELETYPE, "teletype"},
     {wxFONTFAMILY_MAX, "max"},
     {wxFONTFAMILY_UNKNOWN, "unknown"}});
const std::map<std::string, wxFontFamily> WxFontUtils::to_family =
    MapUtils::create_oposit(WxFontUtils::from_family);

const std::map<wxFontStyle, std::string> WxFontUtils::from_style(
    {{wxFONTSTYLE_ITALIC, "italic"},
     {wxFONTSTYLE_SLANT, "slant"},
     {wxFONTSTYLE_NORMAL, "normal"}});
const std::map<std::string, wxFontStyle> WxFontUtils::to_style =
    MapUtils::create_oposit(WxFontUtils::from_style);

const std::map<wxFontWeight, std::string> WxFontUtils::from_weight(
    {{wxFONTWEIGHT_THIN, "thin"},
     {wxFONTWEIGHT_EXTRALIGHT, "extraLight"},
     {wxFONTWEIGHT_LIGHT, "light"},
     {wxFONTWEIGHT_NORMAL, "normal"},
     {wxFONTWEIGHT_MEDIUM, "medium"},
     {wxFONTWEIGHT_SEMIBOLD, "semibold"},
     {wxFONTWEIGHT_BOLD, "bold"},
     {wxFONTWEIGHT_EXTRABOLD, "extraBold"},
     {wxFONTWEIGHT_HEAVY, "heavy"},
     {wxFONTWEIGHT_EXTRAHEAVY, "extraHeavy"}});
const std::map<std::string, wxFontWeight> WxFontUtils::to_weight =
    MapUtils::create_oposit(WxFontUtils::from_weight);

std::optional<wxFont> WxFontUtils::create_wxFont(const FontItem &fi,
                                                 const FontProp &fp)
{
    double     point_size = static_cast<double>(fp.size_in_mm);
    wxFontInfo info(point_size);
    if (fp.family.has_value()) {
        auto it = to_family.find(*fp.style);
        if (it != to_family.end()) info.Family(it->second);
    }
    if (fp.face_name.has_value()) {
        wxString face_name(*fp.face_name);
        info.FaceName(face_name);
    }
    if (fp.style.has_value()) {
        auto it = to_style.find(*fp.style);
        if (it != to_style.end()) info.Style(it->second);
    }
    if (fp.weight.has_value()) {
        auto it = to_weight.find(*fp.weight);
        if (it != to_weight.end()) info.Weight(it->second);
    }

    // Improve: load descriptor instead of store to font property to 3mf
    // switch (fi.type) {
    // case FontItem::Type::wx_lin_font_descr:
    // case FontItem::Type::wx_win_font_descr:
    // case FontItem::Type::file_path:
    // case FontItem::Type::undefined:
    // default:
    //}

    wxFont font(info);
    if (!font.IsOk()) return {};
    return font;
}

void WxFontUtils::update_property(FontProp &font_prop, const wxFont &font)
{
    // The point size is defined as 1/72 of the Anglo-Saxon inch (25.4 mm): it
    // is approximately 0.0139 inch or 352.8 um. But it is too small, so I
    // decide use point size as mm for emboss
    font_prop.size_in_mm = font.GetPointSize(); // *0.3528f;

    wxString    wx_face_name = font.GetFaceName();
    std::string face_name((const char *) wx_face_name.ToUTF8());
    if (!face_name.empty()) font_prop.face_name = face_name;

    wxFontFamily wx_family = font.GetFamily();
    if (wx_family != wxFONTFAMILY_DEFAULT) {
        auto it = from_family.find(wx_family);
        if (it != from_family.end()) font_prop.family = it->second;
    }

    wxFontStyle wx_style = font.GetStyle();
    if (wx_style != wxFONTSTYLE_NORMAL) {
        auto it = from_style.find(wx_style);
        if (it != from_style.end()) font_prop.style = it->second;
    }

    wxFontWeight wx_weight = font.GetWeight();
    if (wx_weight != wxFONTWEIGHT_NORMAL) {
        auto it = from_weight.find(wx_weight);
        if (it != from_weight.end()) font_prop.weight = it->second;
    }
}

void NSVGUtils::flatten_cubic_bez(Polygon &polygon,
                                  float    tessTol,
                                  Vec2f    p1,
                                  Vec2f    p2,
                                  Vec2f    p3,
                                  Vec2f    p4,
                                  int      level)
{
    Vec2f p12  = (p1 + p2) * 0.5f;
    Vec2f p23  = (p2 + p3) * 0.5f;
    Vec2f p34  = (p3 + p4) * 0.5f;
    Vec2f p123 = (p12 + p23) * 0.5f;

    Vec2f pd  = p4 - p1;
    Vec2f pd2 = p2 - p4;
    float d2  = std::abs(pd2.x() * pd.y() - pd2.y() * pd.x());
    Vec2f pd3 = p3 - p4;
    float d3  = std::abs(pd3.x() * pd.y() - pd3.y() * pd.x());
    float d23 = d2 + d3;

    if ((d23 * d23) < tessTol * (pd.x() * pd.x() + pd.y() * pd.y())) {
        polygon.points.emplace_back(p4.x(), p4.y());
        return;
    }

    --level;
    if (level == 0) return;
    Vec2f p234  = (p23 + p34) * 0.5f;
    Vec2f p1234 = (p123 + p234) * 0.5f;
    flatten_cubic_bez(polygon, tessTol, p1, p12, p123, p1234, level);
    flatten_cubic_bez(polygon, tessTol, p1234, p234, p34, p4, level);
}

ExPolygons NSVGUtils::to_ExPolygons(NSVGimage *image,
                                    float      tessTol,
                                    int        max_level)
{
    Polygons polygons;
    for (NSVGshape *shape = image->shapes; shape != NULL;
         shape            = shape->next) {
        if (!(shape->flags & NSVG_FLAGS_VISIBLE)) continue;
        Slic3r::Polygon polygon;
        if (shape->fill.type != NSVG_PAINT_NONE) {
            for (NSVGpath *path = shape->paths; path != NULL;
                 path           = path->next) {
                // Flatten path
                polygon.points.emplace_back(path->pts[0], path->pts[1]);
                for (size_t i = 0; i < path->npts - 1; i += 3) {
                    float *p = &path->pts[i * 2];
                    Vec2f  p1(p[0], p[1]), p2(p[2], p[3]), p3(p[4], p[5]),
                        p4(p[6], p[7]);
                    flatten_cubic_bez(polygon, tessTol, p1, p2, p3, p4,
                                      max_level);
                }
                if (path->closed) {
                    polygons.push_back(polygon);
                    polygon = Slic3r::Polygon();
                }
            }
        }
        polygons.push_back(polygon);
    }

    // Fix Y axis
    for (Polygon &polygon : polygons)
        for (Point &p : polygon.points) p.y() *= -1;

    return Slic3r::union_ex(polygons);
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
