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

void EmbossJob::process(std::unique_ptr<EmbossData> input)
{
    // Changing cursor to busy
    wxBeginBusyCursor();
    ScopeGuard sg([]() { wxEndBusyCursor(); });

    // only for sure
    assert(input != nullptr);

    // check if exist valid font
    if (input->font == nullptr) return;

    const TextConfiguration &cfg  = input->text_configuration;
    const std::string &      text = cfg.text;
    // Do NOT process empty string
    if (text.empty()) return;

    const FontProp &prop = cfg.font_prop;
    ExPolygons shapes = Emboss::text2shapes(*input->font, text.c_str(), prop);

    if (is_stoping()) return;

    // exist 2d shape made by text ?
    // (no shape means that font hasn't any of text symbols)
    if (shapes.empty()) return;

    float scale    = prop.size_in_mm / input->font->ascent;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(prop.emboss / scale);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    TriangleMesh result(Emboss::polygons2model(shapes, project));

    // for sure that some object is created from shape
    if (result.its.indices.empty()) return;

    // center triangle mesh
    Vec3d shift = result.bounding_box().center();
    result.translate(-shift.cast<float>());

    if (is_stoping()) return;
    GLGizmoEmboss::update_emboss_volume(std::move(result), input->volume_name,
                                        input->text_configuration,
                                        input->volume);
}
