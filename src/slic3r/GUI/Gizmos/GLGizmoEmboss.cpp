#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp" // to update title when add text
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"
// TODO: remove include
#include "libslic3r/SVG.hpp" // debug store 

#include "libslic3r/Model.hpp"

#include "imgui/imgui_stdlib.h" // using std::string for inputs
#include "nanosvg/nanosvg.h" // load SVG file

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
    static std::optional<Slic3r::Emboss::Font> load_font(const FontItem &fi);

    static FontItem get_font_item(const wxFont &font);

    // load font used by Operating system as default GUI
    static FontItem get_os_font();
    static std::string get_human_readable_name(const wxFont &font);

    // serialize / deserialize font
    static std::string store_wxFont(const wxFont &font);
    static wxFont      load_wxFont(const std::string &font_descriptor);
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
} // namespace Slic3r

using namespace Slic3r;
using namespace Slic3r::GUI;

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &parent)
    : GLGizmoBase(parent, M_ICON_FILENAME, -2)
    , m_font_selected(0)
    , m_volume(nullptr)
    , m_volume_type(ModelVolumeType::MODEL_PART)
    , m_is_initialized(false) // initialize on first opening gizmo
{
    // TODO: suggest to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web    
}

GLGizmoEmboss::~GLGizmoEmboss() {}

bool GLGizmoEmboss::on_init()
{
    //m_grabbers.emplace_back();
    m_shortcut_key = WXK_CONTROL_T;
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
    check_selection();

    int flag = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse;
    m_imgui->begin(on_get_name(), flag);
    draw_window();

    m_imgui->end(); // 
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
        if (!m_is_initialized) initialize();
        Selection &s = m_parent.get_selection();
        // When add Text on empty plate,
        // Create object with volume
        if (s.is_empty()) {
            set_default_configuration();
            process();
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

    m_font_list = {
        {"NotoSans Regular", Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
        {"NotoSans CJK", Slic3r::resources_dir() + "/fonts/NotoSansCJK-Regular.ttc"}};
    m_font_selected = 0;

    bool is_font_loaded = load_font();
    FontList fl = Emboss::get_font_list();
    m_font_list.insert(m_font_list.end(), fl.begin(), fl.end());
    m_font_list.emplace_back(WxFontUtils::get_os_font());
    while (!is_font_loaded && !m_font_list.empty()) {
        // can't load so erase it from list
        m_font_list.erase(m_font_list.begin() + m_font_selected);
        m_font_selected = 0; // select first    
        is_font_loaded  = load_font();
    }
    sort_fonts();
    set_default_configuration();
}

void GLGizmoEmboss::set_default_configuration() {
    m_text = _u8L("Embossed text");
    m_font_prop = FontProp();
    m_volume_type = ModelVolumeType::MODEL_PART;
    // may be set default font?
}

#include "imgui/imgui_internal.h" // to unfocus input --> ClearActiveID
void GLGizmoEmboss::check_selection()
{
    ModelVolume* vol = get_selected_volume();
    // is same volume selected?
    if (vol!= nullptr && m_volume == vol) return;

    // Do not use focused input value when switch volume(it must swith value)
    if (m_volume != nullptr)
        ImGui::ClearActiveID();

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

    if (id.object_id >= objects.size()) return nullptr;
    ModelObject *object = objects[id.object_id];

    if (id.volume_id >= object->volumes.size()) return nullptr;
    return object->volumes[id.volume_id];
}

// create_text_volume()
bool GLGizmoEmboss::process() {
    if (!m_font.has_value()) return false;

    ExPolygons shapes = Emboss::text2shapes(*m_font, m_text.c_str(), m_font_prop);
    if (shapes.empty()) return false; 
        
    float scale = m_font_prop.size_in_mm / m_font->ascent;
    auto project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(m_font_prop.emboss / scale), scale);
    indexed_triangle_set its = Emboss::polygons2model(shapes, *project);
    if (its.indices.empty()) return false;

    // add object
    TriangleMesh tm(std::move(its));
    tm.repair();

    //tm.WriteOBJFile("text_preview.obj");
         
    GUI_App &app    = wxGetApp();
    Plater * plater = app.plater();
    std::string volume_name = create_volume_name();
    plater->take_snapshot(_L("Add") + " " + volume_name);    
    if (m_volume == nullptr) {
        const Selection &selection = m_parent.get_selection();
        if (selection.is_empty()) {
            // create new object
            app.obj_list()->load_mesh_object(tm, volume_name);
            app.mainframe->update_title();            
            // get new created volume
            m_volume = app.obj_list()->objects()->back()->volumes.front();
            m_volume->text_configuration = create_configuration();

            // load mesh cause close gizmo, soo I open it again
            m_parent.get_gizmos_manager().open_gizmo(GLGizmosManager::EType::Emboss);
            return true;
        } else {
            Model &model = plater->model();
            int object_idx = selection.get_object_idx();
            ModelObject *obj = model.objects[object_idx];
            m_volume = obj->add_volume(std::move(tm), m_volume_type);
        }
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
    m_volume->text_configuration = create_configuration();

    // select new added volume
    ModelObject *mo = m_volume->get_object();
    // Editing object volume change its name
    if (mo->volumes.size() == 1)  mo->name = volume_name;
    ObjectList *obj_list = app.obj_list();
    const ModelObjectPtrs &objs = *obj_list->objects();
    auto item = find(objs.begin(), objs.end(), mo);
    assert(item != objs.end());
    int object_idx = item - objs.begin();
    ModelVolume *new_volume = m_volume; // copy pointer for lambda
    obj_list->select_item([new_volume, object_idx, obj_list]() {
        wxDataViewItemArray items = obj_list->reorder_volumes_and_get_selection(
            object_idx, [new_volume](const ModelVolume *volume) { return volume == new_volume; });
        if (items.IsEmpty()) return wxDataViewItem();
        return items.front();
    });

    if (m_volume->type() == ModelVolumeType::MODEL_PART)
        // update printable state on canvas
        m_parent.update_instance_printable_state_for_object((size_t) object_idx);

    obj_list->selection_changed();
    m_parent.reload_scene(true);
    return true;
}

void GLGizmoEmboss::close() {
    // close gizmo == open it again
    m_parent.get_gizmos_manager().open_gizmo(GLGizmosManager::Emboss);
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

        FontList font_list;
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

void GLGizmoEmboss::draw_window()
{
    if (!m_font.has_value()) {
        ImGui::Text("Warning: No font is selected. Select correct one.");
    }

    draw_font_list();

    if (ImGui::Button(_L("choose font").c_str())) { choose_font_by_dialog(); }

    ImGui::SameLine();
    if (ImGui::Button(_L("use system font").c_str())) {
        wxSystemSettings ss;
        wxFont           f          = ss.GetFont(wxSYS_DEFAULT_GUI_FONT);
        size_t           font_index = m_font_list.size();
        FontItem         fi         = WxFontUtils::get_font_item(f);
        m_font_list.emplace_back(fi);
        bool loaded = load_font(font_index);
    }

    ImGui::SameLine();
    draw_add_button();

    if (ImGui::Button("add svg")) {
        std::string filePath =
            "C:/Users/filip/Downloads/fontawesome-free-5.15.4-web/"
            "fontawesome-free-5.15.4-web/svgs/solid/bicycle.svg";
        filePath         = "C:/Users/filip/Downloads/circle.svg";
        NSVGimage *image = nsvgParseFromFile(filePath.c_str(), "mm", 96.0f);
        ExPolygons polys = NSVGUtils::to_ExPolygons(image);

        for (auto &poly : polys) poly.scale(1e5);
        SVG svg("converted.svg", BoundingBox(polys.front().contour.points));
        svg.draw(polys);

        nsvgDelete(image);
    }

    if (ImGui::InputFloat("Size[in mm]", &m_font_prop.size_in_mm)) {
        if (m_font_prop.size_in_mm < 0.1) m_font_prop.size_in_mm = 10;
        process();
    }
    if (ImGui::InputFloat("Emboss[in mm]", &m_font_prop.emboss)) process();
    if (ImGui::InputFloat("Flatness", &m_font_prop.flatness)) {
        if (m_font.has_value()) m_font->cache.clear();
        process();
    }
    if (ImGui::InputInt("CharGap[in font points]", &m_font_prop.char_gap))
        process();
    if (ImGui::InputInt("LineGap[in font points]", &m_font_prop.line_gap))
        process();

    // ImGui::InputFloat3("Origin", m_orientation.origin.data());
    // if (ImGui::InputFloat3("Normal", m_normal.data())) m_normal.normalize();
    // if (ImGui::InputFloat3("Up", m_up.data())) m_up.normalize();

    ImVec2              input_size(-FLT_MIN, ImGui::GetTextLineHeight() * 6);
    ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_::ImGuiInputTextFlags_AllowTabInput |
        ImGuiInputTextFlags_::ImGuiInputTextFlags_AutoSelectAll
        //| ImGuiInputTextFlags_::ImGuiInputTextFlags_CallbackResize
        //|ImGuiInputTextFlags_::ImGuiInputTextFlags_CtrlEnterForNewLine
        ;

    if (ImGui::InputTextMultiline("##Text", &m_text, input_size, flags))
        process();    

    
    if (ImGui::Button(_L("Close").c_str())) close();

    // Option to create text volume when reselect volumes
    m_imgui->disabled_begin(!m_font.has_value());
    if (m_volume == nullptr) {
        ImGui::SameLine();
        if (ImGui::Button(_L("Generate preview").c_str())) process();
    }
    m_imgui->disabled_end();
}

void GLGizmoEmboss::draw_font_list()
{
    auto &current = m_font_list[m_font_selected];
    if (ImGui::BeginCombo("##font_selector", current.name.c_str())) {
        for (const FontItem &f : m_font_list) {
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

bool GLGizmoEmboss::load_font(size_t font_index)
{
    std::swap(font_index, m_font_selected);
    bool is_loaded = load_font();
    if (!is_loaded) std::swap(font_index, m_font_selected);    
    return is_loaded;
}

bool GLGizmoEmboss::load_font() { 
    if (m_font_selected >= m_font_list.size()) return false;
    auto font = WxFontUtils::load_font(m_font_list[m_font_selected]);
    if (!font.has_value()) return false;
    m_font = font;
    return true;
}

bool GLGizmoEmboss::choose_font_by_dialog() {
    // keep last selected font did not work
    // static wxFontData data; 
    // wxFontDialog      font_dialog((wxWindow *) wxGetApp().mainframe, data);

    wxFontDialog font_dialog(nullptr);
    font_dialog.SetTitle(_L("Select font for Emboss"));
    if (font_dialog.ShowModal() != wxID_OK) return false;
    wxFontData data   = font_dialog.GetFontData();
    wxFont font       = data.GetChosenFont();
    size_t font_index = m_font_list.size();
    m_font_list.emplace_back(WxFontUtils::get_font_item(font));
    if (!load_font(font_index)) { 
        m_font_list.pop_back(); 
        return false;
    }        
    sort_fonts();
    process();
    return true;
}

void GLGizmoEmboss::sort_fonts() {
    // initialize original index locations
    std::vector<size_t> idx(m_font_list.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
        [this](size_t i1, size_t i2) {
            return m_font_list[i1].name < m_font_list[i2].name;
        });
    
    FontList font_list;    
    font_list.reserve(m_font_list.size());
    size_t selected = 0;
    for (const size_t &i : idx) { 
        if (i == m_font_selected) selected = &i - &idx.front();
        font_list.emplace_back(m_font_list[i]); 
    }
    m_font_list = font_list;
    m_font_selected = selected;    
}

void GLGizmoEmboss::add_fonts(const FontList &font_list) {
    m_font_list.insert(m_font_list.end(), font_list.begin(), font_list.end());
    sort_fonts();
}

TextConfiguration GLGizmoEmboss::create_configuration() {
    return TextConfiguration(m_font_list[m_font_selected], m_font_prop, m_text);
}

bool GLGizmoEmboss::load_configuration(ModelVolume *volume)
{
    if (volume == nullptr) return false;
    if (!volume->text_configuration.has_value()) return false;
    const TextConfiguration &configuration = *volume->text_configuration;
    size_t index = m_font_list.size();
    for (const auto &font_item : m_font_list) {        
        if (font_item.type == configuration.font_item.type &&
            font_item.path == configuration.font_item.path) {
            index = &font_item - &m_font_list.front();
        }
    }
    size_t prev_font_selected = m_font_selected;
    // when not in font list add to list
    if (index >= m_font_list.size()) {
        m_font_selected = m_font_list.size();
        add_fonts({configuration.font_item});
    } else {
        m_font_selected = index;
    }
    // When can't load font
    if (!load_font()) {
        // remove bad loadabled font, for correct prev index
        m_font_list.erase(m_font_list.begin() + m_font_selected);
        m_font_selected = prev_font_selected;
        return false;
    }

    m_font_prop = configuration.font_prop;
    m_text = configuration.text;
    m_volume_type = volume->type(); // not neccesary
    m_volume = volume;
    return true;
}

std::string GLGizmoEmboss::create_volume_name()
{
    const size_t max_len = 20;
    return _u8L("Text") + " - " + 
        ((m_text.size() > max_len)? 
        (m_text.substr(0, max_len - 3) + " ..") : m_text);
}

std::optional<Emboss::Font> WxFontUtils::load_font(const FontItem &fi)
{
    switch (fi.type) {
    case FontItem::Type::file_path:
        return Emboss::load_font(fi.path.c_str());
    case FontItem::Type::wx_font_descr:
        return WxFontUtils::load_font(WxFontUtils::load_wxFont(fi.path));
    }
    return {};
}

std::optional<Emboss::Font> WxFontUtils::load_font(const wxFont &font)
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

FontItem WxFontUtils::get_font_item(const wxFont &font)
{
    std::string name     = get_human_readable_name(font);
    std::string fontDesc = store_wxFont(font);
    return FontItem(name, fontDesc, FontItem::Type::wx_font_descr);
}

FontItem WxFontUtils::get_os_font()
{
    wxSystemFont system_font = wxSYS_DEFAULT_GUI_FONT;
    wxFont       font        = wxSystemSettings::GetFont(system_font);
    FontItem     fi          = get_font_item(font);
    fi.name += +" (" + _u8L("OS default") + ")";
    return get_font_item(font);
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

void NSVGUtils::flatten_cubic_bez(Polygon &polygon,
                                  float    tessTol,
                                  Vec2f    p1,
                                  Vec2f    p2,
                                  Vec2f    p3,
                                  Vec2f    p4,
                                  int      level)
{
    Vec2f p12 = (p1 + p2) * 0.5f;
    Vec2f p23 = (p2 + p3) * 0.5f;
    Vec2f p34 = (p3 + p4) * 0.5f;
    Vec2f p123 = (p12 + p23) * 0.5f;

    Vec2f pd = p4 - p1;
    Vec2f pd2 = p2 - p4;
    float d2  = std::abs(pd2.x() * pd.y() - pd2.y() * pd.x());
    Vec2f pd3 = p3 - p4;
    float d3 = std::abs(pd3.x() * pd.y() - pd3.y() * pd.x());
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
                    Vec2f  
                        p1(p[0], p[1]),
                        p2(p[2], p[3]),
                        p3(p[4], p[5]),
                        p4(p[6], p[7]);
                    flatten_cubic_bez(polygon, tessTol, p1, p2, p3, p4, max_level);
                }
                if (path->closed) {
                    polygons.push_back(polygon);
                    polygon = Slic3r::Polygon();
                }
            }
        }
        polygons.push_back(polygon);
    }
    return union_ex(polygons);
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
