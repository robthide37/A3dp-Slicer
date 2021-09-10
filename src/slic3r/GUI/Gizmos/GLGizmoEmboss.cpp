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

namespace Slic3r::GUI {

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &       parent,
                             const std::string &icon_filename,
                             unsigned int       sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_font_list({
        {"NotoSans Regular", Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
        {"NotoSans CJK", Slic3r::resources_dir() + "/fonts/NotoSansCJK-Regular.ttc"}})
    , m_font_selected(0)
    , m_text_size(255)
    , m_text(new char[m_text_size])
    , m_scale(0.01f)
    , m_emboss(5.f) 
    , m_volume(nullptr)
{
    // TODO: suggest to use https://fontawesome.com/
    // (copy & paste) unicode symbols from web

    bool is_font_loaded = load_font();
    add_fonts(Emboss::get_font_list());

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
    auto& current = m_font_list[m_font_selected];
    if (ImGui::BeginCombo("##font_selector", current.name.c_str())) {
        for (const Emboss::FontItem &f : m_font_list) {
            ImGui::PushID((void*)&f.name);
            std::string name = (f.name.size() < m_gui_cfg->max_font_name) ?
                f.name : (f.name.substr(0,m_gui_cfg->max_font_name - 3) + " ..");
            if (ImGui::Selectable(name.c_str(), &f == &current)) {
                m_font_selected = &f - &m_font_list.front();
                load_font();
                process();
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
    ImGui::SameLine();
    if (m_font.has_value()) {
        if (ImGui::BeginCombo("##font_collection_selector", std::to_string(m_font->index).c_str())) {
            for (size_t i = 0; i < m_font->count; ++i) {
                ImGui::PushID(1 << 10 + i);
                if (ImGui::Selectable(std::to_string(i).c_str(),
                                      i == m_font->index)) {
                    m_font->index = i;
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
    }

    static std::string fontName;
    if (ImGui::Button(_L("choose font").c_str())) {
        static wxFontData data; // keep last selected font
        wxFontDialog font_dialog(nullptr, data);
        font_dialog.SetTitle(_L("Select font FFF"));
        if (font_dialog.ShowModal() == wxID_OK) {        
            data        = font_dialog.GetFontData();
            wxFont font = data.GetChosenFont();
            fontName    = boost::nowide::narrow(font.GetFaceName());
            m_font = Emboss::load_font(font.GetHFONT());           
            m_font_glyph_cache.clear();
            process();            
        }
    }
    if (!fontName.empty()) ImGui::Text(fontName.c_str());

    ImGui::SameLine();
    draw_add_button();

    ImGui::InputFloat("Scale", &m_scale);
    ImGui::InputFloat("Emboss", &m_emboss);
    if (ImGui::InputFloat("Flatness", &m_font_prop.flatness))
        m_font_glyph_cache.clear();
    ImGui::InputInt("CharGap", &m_font_prop.char_gap);
    ImGui::InputInt("LineGap", &m_font_prop.line_gap);
    
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
                NotificationManager::NotificationLevel::RegularNotification,
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

    Polygons polygons = Emboss::text2polygons(*m_font, m_text.get(), m_font_prop, m_font_glyph_cache);
    if (polygons.empty()) return;

    auto project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(m_emboss/m_scale), m_scale);
    indexed_triangle_set its = Emboss::polygons2model(polygons, *project);
    if (its.indices.empty()) return;

    // add object
    TriangleMesh tm(std::move(its));
    tm.repair();

    //tm.WriteOBJFile("text_preview.obj");

    const Selection &selection  = m_parent.get_selection();
    int              object_idx = selection.get_object_idx();
    ModelObject *    obj = wxGetApp().plater()->model().objects[object_idx];
    
    //if (m_volume != nullptr) {
    //    // TODO: fix index of m_volume
    //    size_t m_volume_index = obj->volumes.size() - 1;
    //    obj->delete_volume(m_volume_index);
    //} 
    //m_volume = obj->add_volume(std::move(tm), ModelVolumeType::MODEL_PART);    
    
    if (m_volume == nullptr) {
        m_volume = obj->add_volume(std::move(tm), ModelVolumeType::MODEL_PART); 
    } else {
        m_volume->set_mesh(std::move(tm));
        m_volume->set_new_unique_id();
        m_volume->translate(-m_volume->source.mesh_offset);

        m_volume->center_geometry_after_creation(true);
        m_volume->calculate_convex_hull();
        m_volume->get_object()->invalidate_bounding_box();
    }
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

bool GLGizmoEmboss::load_font()
{
    auto font_path = m_font_list[m_font_selected].path.c_str();
    m_font         = Emboss::load_font(font_path);
    m_font_glyph_cache.clear();
    return m_font.has_value();
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

} // namespace Slic3r::GUI
