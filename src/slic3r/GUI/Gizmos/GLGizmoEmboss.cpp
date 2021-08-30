#include "GLGizmoEmboss.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/QuadricEdgeCollapse.hpp"

namespace Slic3r::GUI {

GLGizmoEmboss::GLGizmoEmboss(GLCanvas3D &       parent,
                             const std::string &icon_filename,
                             unsigned int       sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
    , m_fonts({{"NotoSans Regular",
                Slic3r::resources_dir() + "/fonts/NotoSans-Regular.ttf"},
               {"Arial", "C:/windows/fonts/arialbd.ttf"}})
    , m_selected(1)
    , m_text_size(255)
    , m_text(new char[255])
{
    // TODO: suggest to use https://fontawesome.com/
    int index = 0;
    for (char &c : _u8L("Embossed text")) { 
        m_text[index++] = c;
    }
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

    size_t text_size = 255;
    //sizeof(*m_text.get());
    ImVec2 input_size(-FLT_MIN, ImGui::GetTextLineHeight() * 6);
    ImGuiInputTextFlags flags =
        ImGuiInputTextFlags_::ImGuiInputTextFlags_AllowTabInput 
        | ImGuiInputTextFlags_::ImGuiInputTextFlags_AutoSelectAll
        //|ImGuiInputTextFlags_::ImGuiInputTextFlags_CtrlEnterForNewLine
        ;
    if (ImGui::InputTextMultiline("##Text", m_text.get(), text_size, input_size, flags)) {
        
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
