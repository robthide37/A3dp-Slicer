#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/Emboss.hpp"

namespace Slic3r::GUI {

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &       parent,
                             const std::string &icon_filename,
                             unsigned int       sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_fonts({{"NotoSans Regular",
                Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
               {"Arial", "C:/windows/fonts/arialbd.ttf"}})
    , m_selected(1)
    , m_text(_u8L("Embossed text"))
{
    // TODO: suggest to use https://fontawesome.com/
    m_text.reserve(255);
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
    ImGui::Text("welcome in emboss text");
    auto& current = m_fonts[m_selected];
    if (ImGui::BeginCombo("##font_selector", current.name.c_str())) {
        for (const MyFont &f : m_fonts) {
            ImGui::PushID((void*)&f.name);
            if (ImGui::Selectable(f.name.c_str(), &f == &current)) {
                m_selected = &f - &m_fonts.front();
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    ImGui::Button("Add");
    static float scale = 0.01f;
    ImGui::InputFloat("Scale", &scale);
    static float emboss = 100.f;
    ImGui::InputFloat("Emboss", &emboss);

    if(ImGui::Button("Preview")){ 
        
        auto project = std::make_unique<Emboss::ProjectScale>(
            std::make_unique<Emboss::ProjectZ>(emboss)
            ,scale);

        auto font_path  = m_fonts[m_selected].file_path.c_str();
        auto font_opt = Emboss::load_font(font_path);
        if (font_opt.has_value()) {
            Polygons polygons = Emboss::text2polygons(*font_opt, m_text);
            indexed_triangle_set its = Emboss::polygons2model(polygons, *project);
            // add object
            TriangleMesh tm(its);
            tm.repair();

            tm.WriteOBJFile("text_preview.obj");

            const Selection &selection  = m_parent.get_selection();
            int              object_idx = selection.get_object_idx();
            ModelObject *obj = wxGetApp().plater()->model().objects[object_idx];
            ModelVolume *v = obj->volumes.front();
            v->set_mesh(tm);
            v->set_new_unique_id();
            obj->invalidate_bounding_box();
            m_parent.reload_scene(true);
        }        
    }

    ImVec2 input_size(-FLT_MIN, ImGui::GetTextLineHeight() * 6);
    ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_::ImGuiInputTextFlags_AllowTabInput 
        | ImGuiInputTextFlags_::ImGuiInputTextFlags_AutoSelectAll
        //|ImGuiInputTextFlags_::ImGuiInputTextFlags_CtrlEnterForNewLine
        ;
    if (ImGui::InputTextMultiline("##Text", (char *) m_text.c_str(),
        m_text.capacity() + 1, input_size, flags)) {
        
    }

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

void GLGizmoEmboss::close() {
    // close gizmo == open it again
    GLGizmosManager &gizmos_mgr = m_parent.get_gizmos_manager();
    gizmos_mgr.open_gizmo(GLGizmosManager::EType::Simplify);
}
} // namespace Slic3r::GUI
