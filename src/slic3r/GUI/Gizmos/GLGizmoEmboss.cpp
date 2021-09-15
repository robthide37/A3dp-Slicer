#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "libslic3r/Model.hpp"

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

// took from nanosvgrast.h function nsvgRasterize->nsvg__flattenShape
void flatten_cubic_bez(Slic3r::Polygon &polygon,
                       float            tessTol,
                       float            x1,
                       float            y1,
                       float            x2,
                       float            y2,
                       float            x3,
                       float            y3,
                       float            x4,
                       float            y4,
                       int              level)
{
    float x12, y12, x23, y23, x34, y34, x123, y123, x234, y234, x1234,
        y1234;
    float dx, dy, d2, d3;

    if (level == 0) return;

    x12  = (x1 + x2) * 0.5f;
    y12  = (y1 + y2) * 0.5f;
    x23  = (x2 + x3) * 0.5f;
    y23  = (y2 + y3) * 0.5f;
    x34  = (x3 + x4) * 0.5f;
    y34  = (y3 + y4) * 0.5f;
    x123 = (x12 + x23) * 0.5f;
    y123 = (y12 + y23) * 0.5f;

    dx = x4 - x1;
    dy = y4 - y1;
    d2 = std::abs(((x2 - x4) * dy - (y2 - y4) * dx));
    d3 = std::abs(((x3 - x4) * dy - (y3 - y4) * dx));

    if ((d2 + d3) * (d2 + d3) < tessTol * (dx * dx + dy * dy)) {
        polygon.points.emplace_back(x4, y4);
        return;
    }

    --level;
    if (level == 0) return;
    x234  = (x23 + x34) * 0.5f;
    y234  = (y23 + y34) * 0.5f;
    x1234 = (x123 + x234) * 0.5f;
    y1234 = (y123 + y234) * 0.5f;
    flatten_cubic_bez(polygon, tessTol, x1, y1, x12, y12, x123, y123, x1234, y1234, level);
    flatten_cubic_bez(polygon, tessTol, x1234, y1234, x234, y234, x34, y34, x4, y4, level);
}

Slic3r::ExPolygons to_ExPolygons(NSVGimage *image,
                                 float      tessTol   = 10.,
                                 int        max_level = 10)
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
                    flatten_cubic_bez(polygon, tessTol, p[0], p[1], p[2],
                                      p[3], p[4], p[5], p[6], p[7],
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
    return union_ex(polygons);
}
#include "libslic3r/SVG.hpp"

void GLGizmoEmboss::on_render_input_window(float x, float y, float bottom_limit)
{
    if (!m_gui_cfg.has_value()) m_gui_cfg.emplace(GuiCfg());
    check_selection();

    int flag = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse;
    m_imgui->begin(on_get_name(), flag);
    
    if (!m_font.has_value()) { 
        ImGui::Text("Warning: No font is selected. Select correct one.");
    }

    draw_font_list();

    if (ImGui::Button(_L("choose font").c_str())) {
        choose_font_by_dialog();
    }

    ImGui::SameLine();
    if (ImGui::Button(_L("use system font").c_str())) {
        wxSystemSettings ss;
        wxFont f = ss.GetFont(wxSYS_DEFAULT_GUI_FONT);
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
        NSVGimage *image = nsvgParseFromFile(filePath.c_str(), "mm", 96.0f);
        ExPolygons polys = to_ExPolygons(image);

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
    if (ImGui::InputInt("CharGap[in font points]", &m_font_prop.char_gap)) process();    
    if (ImGui::InputInt("LineGap[in font points]", &m_font_prop.line_gap)) process();
        
    //ImGui::InputFloat3("Origin", m_orientation.origin.data());
    //if (ImGui::InputFloat3("Normal", m_normal.data())) m_normal.normalize();    
    //if (ImGui::InputFloat3("Up", m_up.data())) m_up.normalize();
    
    m_imgui->disabled_begin(!m_font.has_value());

    // create default text
    if (m_volume == nullptr) { 
        if (ImGui::Button("Generate preview")) process();
    }

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
        set_max_text_size(static_cast<size_t>(max_text_size));
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
        if(!set_volume()) set_default_configuration();

        // when open by hyperlink it needs to show up
        m_parent.reload_scene(true);
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

void GLGizmoEmboss::set_default_configuration() {
    set_text(_u8L("Embossed text"));
    m_font_prop = FontProp();
    // may be set default font?
}

void GLGizmoEmboss::check_selection()
{
    ModelVolume* vol = get_selected_volume();
    // is same volume selected?
    if (vol!= nullptr && m_volume == vol) return;

    // Do not use actual edited value when switch volume
    ImGui::SetKeyboardFocusHere(0); 

    // is selected volume embossed?
    if (vol!= nullptr && vol->text_configuration.has_value()) {
        m_volume = vol;
        load_configuration(*vol->text_configuration);
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
    int              object_idx = selection.get_object_idx();
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
    m_volume->text_configuration = create_configuration();

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

void GLGizmoEmboss::set_text(const std::string &text) {
    if (text.size() > m_text_size-1)
        set_max_text_size(text.size() + 1);
    
    int index = 0;
    for (const char &c : text) m_text[index++] = c; 
    m_text[index] = '\0';
}

void GLGizmoEmboss::set_max_text_size(size_t size) {
    if (size < 4) size = 4;
    std::unique_ptr<char[]> newData(new char[size]);
    size_t                  index = 0;
    while ((index + 1) < size) {
        if (m_text.get()[index] == '\0') break;
        newData.get()[index] = m_text.get()[index];
        ++index;
    }
    newData.get()[index] = '\0';
    m_text               = std::move(newData);
    m_text_size          = size;
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
    wxSystemSettings ss;
    wxFont ss_font = ss.GetFont(wxSYS_ANSI_VAR_FONT);
    FontItem fi = get_font_item(ss_font);
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

bool GLGizmoEmboss::set_volume()
{
    ModelVolume *vol = get_selected_volume();
    // Is selected only one volume
    if (vol == nullptr) return false;

    // Is volume created by Emboss?
    if (!vol->text_configuration.has_value()) return false;

    // set selected volume
    m_volume = vol;
    load_configuration(*vol->text_configuration);
    return true;
}

TextConfiguration GLGizmoEmboss::create_configuration() {
    std::string text((const char *) m_text.get());
    return TextConfiguration(m_font_list[m_font_selected], m_font_prop, text);
}

bool GLGizmoEmboss::load_configuration(const TextConfiguration &configuration)
{
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
    set_text(configuration.text);
    return true;
}

std::string GLGizmoEmboss::create_volume_name()
{
    size_t      max_len = 20;
    std::string text((const char *)m_text.get());
    if (text.size() > max_len) 
        text = text.substr(0, max_len - 3) + " ..";    
    return _u8L("Text") + " - " + text;
}

// any existing icon filename to not influence GUI
const std::string GLGizmoEmboss::M_ICON_FILENAME = "cut.svg";
