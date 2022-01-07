#include "EmbossJob.hpp"

#include "libslic3r/Model.hpp"

#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Gizmos/GLGizmoEmboss.hpp"

using namespace Slic3r;
using namespace GUI;

void EmbossJob::process(Ctl &ctl) {
    // Changing cursor to busy must be inside main thread 
    // GTK is not thread safe.
    //wxBeginBusyCursor();
    //ScopeGuard sg([]() { wxEndBusyCursor(); });

    // only for sure
    assert(m_input != nullptr);

    // check if exist valid font
    if (m_input->font == nullptr) return;

    const TextConfiguration &cfg  = m_input->text_configuration;
    const std::string &      text = cfg.text;
    // Do NOT process empty string
    if (text.empty()) return;

    const FontProp &prop = cfg.font_item.prop;
    ExPolygons shapes = Emboss::text2shapes(*m_input->font, text.c_str(), prop);

    if (ctl.was_canceled()) return;

    // exist 2d shape made by text ?
    // (no shape means that font hasn't any of text symbols)
    if (shapes.empty()) return;

    float scale    = prop.size_in_mm / m_input->font->ascent;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(prop.emboss / scale);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    m_result = TriangleMesh(Emboss::polygons2model(shapes, project));

    if (ctl.was_canceled()) return;

    // center triangle mesh
    Vec3d shift = m_result.bounding_box().center();
    m_result.translate(-shift.cast<float>());    
}

void EmbossJob::finalize(bool canceled, std::exception_ptr &) {
    if (canceled) return;

    // for sure that some object is created from shape
    if (m_result.its.indices.empty()) return;

    GLGizmoEmboss::update_emboss_volume(std::move(m_result),
                                        m_input->volume_name,
                                        m_input->text_configuration,
                                        m_input->volume);
}
