#include "EmbossJob.hpp"

#include "libslic3r/Model.hpp"

#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"

using namespace Slic3r;
using namespace GUI;

namespace Priv {
static void finalize(const EmbossData &input, const indexed_triangle_set &result);
// TODO: May be move to objec list utils
static void select_volume(ModelVolume *volume);
} // namespace Priv

void EmbossJob::process(std::unique_ptr<EmbossData> input, StopCondition is_stop)
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

    if (is_stop()) return;

    // exist 2d shape made by text ?
    // (no shape means that font hasn't any of text symbols)
    if (shapes.empty()) return;

    float scale    = prop.size_in_mm / input->font->ascent;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(prop.emboss / scale);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    auto   its = std::make_unique<indexed_triangle_set>(Emboss::polygons2model(shapes, project));

    if (is_stop()) return;

    // for sure that some object is created from shape
    if (its->indices.empty()) return;

    // Call setting model in UI thread

    // Not work !!! CallAfter use only reference to lambda not lambda itself
    // wxGetApp().plater()->CallAfter(
    //    [its2 = std::move(its), input2 = std::move(input)]() {
    //        Priv::finalize(*input2, *its2);
    //    });

    struct Data
    {
        std::unique_ptr<EmbossData>           input;
        std::unique_ptr<indexed_triangle_set> result;
        Data(std::unique_ptr<EmbossData>           input,
             std::unique_ptr<indexed_triangle_set> result)
            : input(std::move(input)), result(std::move(result))
        {}
    };
    Data *data = new Data(std::move(input), std::move(its));
    // How to proof that call, will be done on exit?
    wxGetApp().plater()->CallAfter([data]() {
        // because of finalize exception delete must be in ScopeGuard
        ScopeGuard sg([data]() { delete data; }); 
        Priv::finalize(*data->input, *data->result);        
    });
}

void Priv::finalize(const EmbossData &input, const indexed_triangle_set &result)
{
    GUI_App &        app    = wxGetApp(); // may be move to input
    Plater *         plater = app.plater();
    GLCanvas3D *     canvas = plater->canvas3D();
    GLGizmosManager &manager= canvas->get_gizmos_manager();

    // TODO: Solve gizmo with empty selection first
    // Check emboss gizmo is still open
    //if (manager.get_current_type() != GLGizmosManager::Emboss) return;

    // it is sad, but there is no move constructor --> copy
    TriangleMesh tm(std::move(result)); 
    
    // center triangle mesh
    Vec3d shift = tm.bounding_box().center();
    tm.translate(-shift.cast<float>());
    const std::string &name = input.volume_name;

    plater->take_snapshot(_L("Emboss text") + ": " + name);
    ModelVolume *volume = input.volume_ptr;
    assert(volume != nullptr);

    // update volume
    volume->set_mesh(std::move(tm));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();
    volume->name               = name;
    volume->text_configuration = input.text_configuration;

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    //select_volume(volume);

    // Job promiss to refresh is not working
    canvas->reload_scene(true);
}

void Priv::select_volume(ModelVolume *volume)
{
    if (volume == nullptr) return;

    GUI_App &   app      = wxGetApp();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = app.plater()->canvas3D();

    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) {
        return vol == volume;
    };

    const Selection &selection = canvas->get_selection();
    int obj_idx = selection.get_object_idx();
    wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(obj_idx, add_to_selection);

    if (!sel.IsEmpty()) obj_list->select_item(sel.front());

    if (volume->type() == ModelVolumeType::MODEL_PART)
        // update printable state on canvas
        canvas->update_instance_printable_state_for_object((size_t)obj_idx);

    obj_list->selection_changed();
}