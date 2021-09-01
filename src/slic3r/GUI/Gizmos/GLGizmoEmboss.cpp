#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "libslic3r/Model.hpp"

namespace Slic3r::GUI {

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &       parent,
                             const std::string &icon_filename,
                             unsigned int       sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_fonts({
        {"NotoSans Regular", Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
        {"NotoSans CJK", Slic3r::resources_dir() + "/fonts/NotoSansCJK-Regular.ttc"},
        //{"Font Awsome", Slic3r::resources_dir() + "/fonts/fa-solid-900.ttf"},
        {"Arial", "C:/windows/fonts/arialbd.ttf"}})
    , m_fonts_selected(0)
    , m_text_size(255)
    , m_text(new char[m_text_size])
    , m_scale(0.01f)
    , m_emboss(5.f)
    , m_flatness(2.f)
{
    load_font();
    // TODO: suggest to use https://fontawesome.com/
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
    int flag = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse;
    m_imgui->begin(on_get_name(), flag);

    auto& current = m_fonts[m_fonts_selected];
    if (ImGui::BeginCombo("##font_selector", current.name.c_str())) {
        for (const MyFont &f : m_fonts) {
            ImGui::PushID((void*)&f.name);
            if (ImGui::Selectable(f.name.c_str(), &f == &current)) {
                m_fonts_selected = &f - &m_fonts.front();
                load_font();                
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    m_imgui->disabled_begin(!m_font.has_value() || m_font->count == 1);
    if (ImGui::BeginCombo("##font_collection_selector", std::to_string(m_font->index).c_str())) {
        for (size_t i = 0; i < m_font->count; ++i) {
            ImGui::PushID(1<<10 + i);
            if (ImGui::Selectable(std::to_string(i).c_str(), i == m_font->index)) {
                m_font->index = i;
            }
            ImGui::PopID();   
        }
        ImGui::EndCombo();
    }
    m_imgui->disabled_end();

    ImGui::SameLine();
    draw_add_button();

    ImGui::InputFloat("Scale", &m_scale);
    ImGui::InputFloat("Emboss", &m_emboss);
    ImGui::InputFloat("Flatness", &m_flatness);
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
    int max_text_size = m_text_size;
    if (ImGui::InputInt("max text size", &max_text_size, 8, 64)) {
        if (max_text_size < 4) max_text_size = 4;
        std::unique_ptr<char[]> newData(new char[max_text_size]);
        size_t index = 0;
        while (index < max_text_size-1) {
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

        // refuse outgoing during simlification
        if (false) {
            GLGizmoBase::m_state = GLGizmoBase::On;
            auto notification_manager = wxGetApp().plater()->get_notification_manager();
            notification_manager->push_notification(
                NotificationType::CustomNotification,
                NotificationManager::NotificationLevel::RegularNotification,
                _u8L("ERROR: Wait until ends or Cancel process."));
            return;
        }

    } else if (GLGizmoBase::m_state == GLGizmoBase::On) {
        // when open by hyperlink it needs to show up
        //request_rerender();
    }
}

void GLGizmoEmboss::process() {
    auto project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(m_emboss/m_scale), m_scale);

    Polygons polygons        = Emboss::text2polygons(*m_font, m_text.get(), m_flatness);
    if (polygons.empty()) return;

    indexed_triangle_set its = Emboss::polygons2model(polygons, *project);
    if (its.indices.empty()) return;

    // add object
    TriangleMesh tm(its);
    tm.repair();

    tm.WriteOBJFile("text_preview.obj");

    const Selection &selection  = m_parent.get_selection();
    int              object_idx = selection.get_object_idx();
    ModelObject *    obj = wxGetApp().plater()->model().objects[object_idx];
    ModelVolume *    v   = obj->volumes.front();
    v->set_mesh(tm);
    v->set_new_unique_id();
    obj->invalidate_bounding_box();
    m_parent.reload_scene(true);
}

void GLGizmoEmboss::close() {
    // close gizmo == open it again
    GLGizmosManager &gizmos_mgr = m_parent.get_gizmos_manager();
    gizmos_mgr.open_gizmo(GLGizmosManager::EType::Simplify);
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

        for (auto &input_file : input_files) { 
            std::string name = input_file.AfterLast('\\').c_str();
            std::string path = input_file.c_str();
            m_fonts.emplace_back(name, path);
        }
        // set last added font as active
        m_fonts_selected = m_fonts.size() - 1;
        load_font();
        //load_files(input_files);
    }
    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::Text("add file with font(.ttf, .ttc)");
        ImGui::EndTooltip();
    }
}

void GLGizmoEmboss::load_font()
{
    auto font_path = m_fonts[m_fonts_selected].file_path.c_str();
    m_font         = Emboss::load_font(font_path);
}

} // namespace Slic3r::GUI
