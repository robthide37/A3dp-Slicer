#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "libslic3r/Model.hpp"

#include <wx/font.h>
#include <wx/fontdlg.h>


namespace Slic3r {
class WxFontUtils
{
public:
    WxFontUtils() = delete;

    // os specific load of wxFont
    static std::optional<Slic3r::Emboss::Font> load_font(const wxFont &font);
    // Must be in gui because of wxWidget
    static std::optional<Slic3r::Emboss::Font> load_font(const Emboss::FontItem &fi);


    static Slic3r::Emboss::FontItem get_font_item(const wxFont &font);

    // load font used by Operating system as default GUI
    static Slic3r::Emboss::FontItem get_os_font();
    static std::string get_human_readable_name(const wxFont &font);

    // serialize / deserialize font
    static std::string store_wxFont(const wxFont &font);
    static wxFont      load_wxFont(const std::string &font_descriptor);
};
} // namespace Slic3r

using namespace Slic3r;
using namespace Slic3r::GUI;

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -2)
    , m_font_list({
        {"NotoSans Regular", Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
        {"NotoSans CJK", Slic3r::resources_dir() + "/fonts/NotoSansCJK-Regular.ttc"}})
    , m_font_selected(0)
    , m_text_size(255)
    , m_text(new char[m_text_size])
    , m_volume(nullptr)
    , m_volume_type(ModelVolumeType::MODEL_PART)
{
    // TODO: suggest to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web
    bool is_font_loaded = load_font();
    add_fonts(Emboss::get_font_list());
    add_fonts({WxFontUtils::get_os_font()});

    if (!is_font_loaded) { 
        // can't load so erase it from list
        m_font_list.erase(m_font_list.begin() + m_font_selected);
        m_font_selected = 0; // select first
        do{
            is_font_loaded = load_font();
            if (!is_font_loaded) m_font_list.erase(m_font_list.begin());
        } while (!is_font_loaded && !m_font_list.empty());
    }

    int index = 0;
    for (char &c : _u8L("Embossed text")) { m_text[index++] = c; }
    m_text[index] = '\0';
}

GLGizmoEmboss::~GLGizmoEmboss() {}

bool GLGizmoEmboss::on_init()
{
    //m_grabbers.emplace_back();
    m_shortcut_key = WXK_CONTROL_Q;
    return true;
}

std::string GLGizmoEmboss::on_get_name() const
{
    return (_L("Emboss")).ToUTF8().data();
}

void GLGizmoEmboss::on_render() {}
void GLGizmoEmboss::on_render_for_picking() {}

void GLGizmoEmboss::on_render_input_window(float x, float y, float bottom_limit)
{
    if (!m_gui_cfg.has_value()) m_gui_cfg.emplace(GuiCfg());

    int flag = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse;
    m_imgui->begin(on_get_name(), flag);
    draw_font_list();

    static std::string fontName;
    if (ImGui::Button(_L("choose font").c_str())) {
        static wxFontData data; // keep last selected font
        wxFontDialog font_dialog((wxWindow*)wxGetApp().mainframe, data);
        font_dialog.SetTitle(_L("Select font for Emboss"));
        if (font_dialog.ShowModal() == wxID_OK) {        
            data        = font_dialog.GetFontData();
            wxFont font = data.GetChosenFont();            
            auto fontOpt = WxFontUtils::load_font(font);
            if (fontOpt.has_value()) {
                Emboss::FontItem fontItem  = WxFontUtils::get_font_item(font);
                m_font_selected           = m_font_list.size();
                add_fonts({fontItem});
                m_font = fontOpt;
                process();            
            }
        }
    }
    if (!fontName.empty()) ImGui::Text(fontName.c_str());

    ImGui::SameLine();
    draw_add_button();

    ImGui::InputFloat("Size[in mm]", &m_font_prop.size_in_mm);
    ImGui::InputFloat("Emboss[in mm]", &m_font_prop.emboss);
    if (ImGui::InputFloat("Flatness", &m_font_prop.flatness))
        if(m_font.has_value()) m_font->cache.clear();
    ImGui::InputInt("CharGap[in font points]", &m_font_prop.char_gap);
    ImGui::InputInt("LineGap[in font points]", &m_font_prop.line_gap);    
    ImGui::InputFloat3("Origin", m_orientation.origin.data());
    //if (ImGui::InputFloat3("Normal", m_normal.data())) m_normal.normalize();    
    //if (ImGui::InputFloat3("Up", m_up.data())) m_up.normalize();
    
    m_imgui->disabled_begin(!m_font.has_value());
    if (ImGui::Button("Preview")) process();
    m_imgui->disabled_end();
    ImVec2 input_size(-FLT_MIN, ImGui::GetTextLineHeight() * 6);
    ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_::ImGuiInputTextFlags_AllowTabInput 
        | ImGuiInputTextFlags_::ImGuiInputTextFlags_AutoSelectAll 
        //| ImGuiInputTextFlags_::ImGuiInputTextFlags_CallbackResize
        //|ImGuiInputTextFlags_::ImGuiInputTextFlags_CtrlEnterForNewLine
        ;

     
    if (ImGui::InputTextMultiline("##Text", m_text.get(), m_text_size, input_size, flags)) {
        process();
    }

    // change text size
    int max_text_size = static_cast<int>(m_text_size);
    if (ImGui::InputInt("max text size", &max_text_size, 8, 64)) {
        if (max_text_size < 4) max_text_size = 4;
        std::unique_ptr<char[]> newData(new char[max_text_size]);
        size_t index = 0;
        while ((index+1) < max_text_size) {
            if (m_text.get()[index] == '\0') break;
            newData.get()[index] = m_text.get()[index];
            ++index;
        }
        newData.get()[index] = '\0';
        m_text = std::move(newData);
        m_text_size = max_text_size;
    }

    // draw 2d triangle in IMGUI
    ImVec2 t0(25, 25);
    ImVec2 t1(150, 5);
    ImVec2 t2(10, 100);
    ImU32  c = ImGui::ColorConvertFloat4ToU32(ImVec4(.0, .8, .2, 1.));
    ImGui::GetOverlayDrawList()->AddTriangleFilled(t0, t1, t2, c);
    
    m_imgui->end();
}


bool GLGizmoEmboss::on_is_activable() const
{
    return !m_parent.get_selection().is_empty();
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
    } else if (GLGizmoBase::m_state == GLGizmoBase::On) {
        // when open by hyperlink it needs to show up
        //request_rerender();
    }
}

// IMPROVE: Do not use gizmo_event - especialy smth with prefix SLA,
// use Bind into wxGLCanvas?
bool GLGizmoEmboss::gizmo_event(SLAGizmoEventType action,
                         const Vec2d &     mouse_position,
                         bool              shift_down,
                         bool              alt_down,
                         bool              control_down)
{
    /* if (action == SLAGizmoEventType::LeftUp) {
        const Camera &       camera    = wxGetApp().plater()->get_camera();
        const Selection &    selection = m_parent.get_selection();
        const ModelObject *  mo = m_c->selection_info()->model_object();
        const ModelInstance *mi = mo->instances[selection.get_instance_idx()];
        const Transform3d &  instance_trafo = mi->get_transformation()
                                                .get_matrix();

        // Precalculate transformations of individual meshes.
        std::vector<Transform3d> trafo_matrices;
        for (const ModelVolume *mv : mo->volumes)
            if (mv->is_model_part())
                trafo_matrices.emplace_back(instance_trafo * mv->get_matrix());

        Vec3f  normal      = Vec3f::Zero();
        Vec3f  hit         = Vec3f::Zero();
        size_t facet       = 0;
        Vec3f  closest_hit = Vec3f::Zero();
        double closest_hit_squared_distance =
            std::numeric_limits<double>::max();
        size_t closest_facet       = 0;
        int    closest_hit_mesh_id = -1;

        // Cast a ray on all meshes, pick the closest hit and save it for the
        // respective mesh
        for (int mesh_id = 0; mesh_id < int(trafo_matrices.size());
             ++mesh_id) {
            if (m_c->raycaster()->raycasters()[mesh_id]->unproject_on_mesh(
                    mouse_position, trafo_matrices[mesh_id], camera, hit,
                    normal, m_c->object_clipper()->get_clipping_plane(),
                    &facet)) {

                // In case this hit is clipped, skip it.

                // Is this hit the closest to the camera so far?
                double hit_squared_distance = (camera.get_position() -
                                               trafo_matrices[mesh_id] *
                                                   hit.cast<double>())
                                                  .squaredNorm();
                if (hit_squared_distance < closest_hit_squared_distance) {
                    closest_hit_squared_distance = hit_squared_distance;
                    closest_facet                = facet;
                    closest_hit_mesh_id          = mesh_id;
                    closest_hit                  = hit;
                }
            }
        }

        // get intersection

    }
    */
    return false;
}

void GLGizmoEmboss::process() {
    if (!m_font.has_value()) return;

    Polygons polygons = Emboss::text2polygons(*m_font, m_text.get(), m_font_prop);
    if (polygons.empty()) return; 
        
    float scale = m_font_prop.size_in_mm / m_font->ascent;
    auto project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(m_font_prop.emboss / scale), scale);
    indexed_triangle_set its = Emboss::polygons2model(polygons, *project);
    if (its.indices.empty()) return;

    // add object
    TriangleMesh tm(std::move(its));
    tm.repair();

    //tm.WriteOBJFile("text_preview.obj");

    const Selection &selection  = m_parent.get_selection();
    int              object_idx = selection.get_object_idx();
    GUI_App &        app        = wxGetApp();
    Plater *         plater     = app.plater();
    ModelObject *    obj        = plater->model().objects[object_idx];
      
    std::string volume_name = create_volume_name();
    plater->take_snapshot(_L("Add") + " " + volume_name);
    
    if (m_volume == nullptr) {
        m_volume = obj->add_volume(std::move(tm), m_volume_type); 
    } else {
        m_volume->set_mesh(std::move(tm));
        m_volume->set_new_unique_id();
        m_volume->translate(-m_volume->source.mesh_offset);

        m_volume->center_geometry_after_creation(true);
        m_volume->calculate_convex_hull();
        m_volume->get_object()->invalidate_bounding_box();
    }
    m_volume->name = volume_name;

    // set a default extruder value, since user can't add it manually
    m_volume->config.set_key_value("extruder", new ConfigOptionInt(0));

    // select new added volume
    ModelVolume *new_volume = m_volume;
    ObjectList *obj_list = app.obj_list();
    obj_list->select_item([new_volume, object_idx, obj_list]() {
        wxDataViewItemArray items = obj_list->reorder_volumes_and_get_selection(
            object_idx, [new_volume](const ModelVolume *volume) { return volume == new_volume; });
        if (items.IsEmpty()) return wxDataViewItem();
        return items.front();
    });

    if (m_volume_type == ModelVolumeType::MODEL_PART)
        // update printable state on canvas
        m_parent.update_instance_printable_state_for_object((size_t) object_idx);

    obj_list->selection_changed();
    m_parent.reload_scene(true);
}

void GLGizmoEmboss::close() {
    // close gizmo == open it again
    GLGizmosManager &gizmos_mgr = m_parent.get_gizmos_manager();
    gizmos_mgr.open_gizmo(GLGizmosManager::EType::Emboss);
}

void GLGizmoEmboss::draw_add_button() {
    if (ImGui::Button(_u8L("Add").c_str())) {
        wxArrayString input_files;
        wxString      fontDir = wxEmptyString;
        wxString      selectedFile = wxEmptyString;
        wxFileDialog  dialog(
            nullptr, _L("Choose one or more files (TTF, TTC):"),
                            fontDir, selectedFile,
                            file_wildcards(FT_FONTS),
                            wxFD_OPEN | wxFD_MULTIPLE | wxFD_FILE_MUST_EXIST);
        if (dialog.ShowModal() == wxID_OK) dialog.GetPaths(input_files);
        if (input_files.IsEmpty()) return;

        Emboss::FontList font_list;
        font_list.reserve(input_files.size());
        for (auto &input_file : input_files) {
            std::string path = std::string(input_file.c_str());
            size_t      pos  = path.find_last_of('\\');
            size_t      pos2  = path.find_last_of('.');
            std::string name = path.substr(pos + 1, pos2 - pos-1);
            font_list.emplace_back(name, path);
        }
        // set last added font as active
        m_font_selected = m_font_list.size() + font_list.size() - 1;
        add_fonts(font_list);
        load_font();
    }
    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::Text("add file with font(.ttf, .ttc)");
        ImGui::EndTooltip();
    }
}

void GLGizmoEmboss::draw_font_list()
{
    auto &current = m_font_list[m_font_selected];
    if (ImGui::BeginCombo("##font_selector", current.name.c_str())) {
        for (const Emboss::FontItem &f : m_font_list) {
            ImGui::PushID((void *) &f.name);
            std::string name =
                (f.name.size() < m_gui_cfg->max_font_name) ?
                    f.name :
                    (f.name.substr(0, m_gui_cfg->max_font_name - 3) + " ..");
            if (ImGui::Selectable(name.c_str(), &f == &current)) {
                size_t prev_font_selected = m_font_selected;
                m_font_selected = &f - &m_font_list.front();
                if (!load_font()) { 
                    m_font_selected = prev_font_selected;
                } else {
                    process();
                }
            }
            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                ImGui::Text((f.name + " " + f.path).c_str());
                ImGui::EndTooltip();
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    if (m_font.has_value() && m_font->count > 1) {
        ImGui::SameLine();
        if (ImGui::BeginCombo("##font_collection_selector",
                              std::to_string(m_font->index).c_str())) {
            for (size_t i = 0; i < m_font->count; ++i) {
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
}

bool GLGizmoEmboss::load_font() { 
    if (m_font_selected >= m_font_list.size()) return false;
    auto font = WxFontUtils::load_font(m_font_list[m_font_selected]);
    if (!font.has_value()) return false;
    m_font = font;
    return true;
}

std::optional<Emboss::Font> WxFontUtils::load_font(const Emboss::FontItem &fi)
{
    switch (fi.type) {
    case Emboss::FontItem::Type::file_path:
        return Emboss::load_font(fi.path.c_str());
    case Emboss::FontItem::Type::wx_font_descr:
        return WxFontUtils::load_font(WxFontUtils::load_wxFont(fi.path));
    }
    return {};
}

std::optional<Slic3r::Emboss::Font> WxFontUtils::load_font(const wxFont &font)
{
    if (!font.IsOk()) return {};
#ifdef _WIN32
    return Slic3r::Emboss::load_font(font.GetHFONT());
#elif __linux__ 
    // use file path 
    return {};
#elif __APPLE__
    const wxNativeFontInfo *info = font.GetNativeFontInfo();
    CTFontDescriptorRef     descriptor = info3->GetCTFontDescriptor();
    CFDictionaryRef attribs = CTFontDescriptorCopyAttributes(descriptor);
    CFStringRef url = (CFStringRef)CTFontDescriptorCopyAttribute(descriptor, kCTFontURLAttribute);
    std::string str(CFStringGetCStringPtr(CFURLGetString(anUrl),kCFStringEncodingUTF8));
    return Emboss::load_font(str);
#endif
}

Slic3r::Emboss::FontItem WxFontUtils::get_font_item(const wxFont &font)
{
    std::string name     = get_human_readable_name(font);
    std::string fontDesc = store_wxFont(font);
    return Emboss::FontItem(name, fontDesc, Emboss::FontItem::Type::wx_font_descr);
}

Slic3r::Emboss::FontItem WxFontUtils::get_os_font()
{
    wxSystemSettings ss;
    wxFont ss_font = ss.GetFont(wxSYS_ANSI_VAR_FONT);
    Emboss::FontItem fi = get_font_item(ss_font);
    fi.name += +" (" + _u8L("OS default") + ")";
    return get_font_item(ss_font);
}

std::string WxFontUtils::get_human_readable_name(const wxFont &font)
{
    if (!font.IsOk()) return "Font is NOT ok.";
    // Face name is optional in wxFont
    if (!font.GetFaceName().empty()) {
        return std::string(font.GetFaceName().c_str());
    } else {
        return std::string((
                font.GetFamilyString() + " " + 
                font.GetStyleString() + " " +
                font.GetWeightString()
            ).c_str());
    }
}

std::string WxFontUtils::store_wxFont(const wxFont &font)
{
    //wxString os = wxPlatformInfo::Get().GetOperatingSystemIdName();    
    wxString font_descriptor = font.GetNativeFontInfoDesc();
    return std::string(font_descriptor.c_str());
}

wxFont WxFontUtils::load_wxFont(const std::string &font_descriptor)
{
    wxString font_descriptor_wx(font_descriptor);
    return wxFont(font_descriptor_wx);
}

void GLGizmoEmboss::sort_fonts() {
    // initialize original index locations
    std::vector<size_t> idx(m_font_list.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
        [this](size_t i1, size_t i2) {
            return m_font_list[i1].name < m_font_list[i2].name;
        });
    
    Emboss::FontList font_list;    
    font_list.reserve(m_font_list.size());
    size_t selected = 0;
    for (const size_t &i : idx) { 
        if (i == m_font_selected) selected = &i - &idx.front();
        font_list.emplace_back(m_font_list[i]); 
    }
    m_font_list = font_list;
    m_font_selected = selected;    
}

void GLGizmoEmboss::add_fonts(const Emboss::FontList &font_list) {
    m_font_list.insert(m_font_list.end(), font_list.begin(), font_list.end());
    sort_fonts();
}

std::string GLGizmoEmboss::create_volume_name()
{
    size_t      max_len = 20;
    std::string text((const char *)m_text.get());
    if (text.size() > max_len) 
        text = text.substr(0, max_len - 3) + " ..";    
    return _u8L("Text") + ": " + text;
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
