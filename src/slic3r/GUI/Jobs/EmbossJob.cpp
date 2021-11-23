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

//EmbossJob::EmbossJob(): Job(std::make_shared<NotificationProgressIndicator>(wxGetApp().plater()->get_notification_manager())){}
EmbossJob::EmbossJob() : Job(std::make_shared<EmbossJob::Progress>()) {}

EmbossJob::~EmbossJob() { 
    Job::cancel(); 
    Job::join();
}

void EmbossJob::restart(const Data &data)
{
    if (Job::is_running()) { 
        Job::cancel();
        std::lock_guard lk(m_mutex);
        bool create_restart = !m_data_next.has_value();
        m_data_next = data; // copy

        if (create_restart) {
            wxGetApp().plater()->CallAfter([this]() {
                const Plater *plater = wxGetApp().plater();
                if (plater == nullptr) return;
                const GLCanvas3D *canvas = plater->canvas3D();
                if (canvas == nullptr) return;
                const GLGizmosManager &mng = canvas->get_gizmos_manager();
                // check if simplify is still activ gizmo
                if (mng.get_current_type() != GLGizmosManager::Emboss) return;

                // when emboss is current, this object exist
                Job::join();
                {
                    // Not neccessary to lock mutex because job is stoped, 
                    // but for sure that no other thread try to restart job.
                    // Is not critical and could be there.
                    std::lock_guard lk(m_mutex);
                    assert(m_data_next.has_value());
                    m_data = std::move(m_data_next);
                    // after move optional, data is UNDEFINED
                    m_data_next = {};
                }
                Job::start();
            });
        }
    } else {
        m_data = data; // copy
        Job::start();
    }
}

void EmbossJob::prepare() {}

void EmbossJob::process() {
    // check if data was set
    if (!m_data.has_value()) return;
    // check if exist valid font
    if (m_data->font == nullptr) return;

    // Do NOT process empty string
    const TextConfiguration &cfg  = m_data->text_configuration;
    const std::string &text = cfg.text;
    if (text.empty()) return;

    const FontProp &prop = cfg.font_prop;
    ExPolygons shapes = Emboss::text2shapes(*m_data->font, text.c_str(), prop);
    
    // exist 2d shape made by text ?
    // (no shape means that font hasn't any of text symbols)
    if (shapes.empty()) return;

    float scale   = prop.size_in_mm / m_data->font->ascent;
    auto  project = std::make_unique<Emboss::ProjectScale>(
        std::make_unique<Emboss::ProjectZ>(prop.emboss / scale), scale);

    indexed_triangle_set its = Emboss::polygons2model(shapes, *project);

    // for sure that some object is created from shape
    if (its.indices.empty()) return;
    m_result = its;
}

void EmbossJob::finalize() {
    // check result was created
    if (!m_result.has_value()) return;
    if (m_result->indices.empty()) return;
    // move object data
    TriangleMesh tm(std::move(*m_result));
    m_result = {};

    // center triangle mesh
    Vec3d shift = tm.bounding_box().center();
    tm.translate(-shift.cast<float>());

    GUI_App &   app    = wxGetApp();
    Plater *    plater = app.plater();
    GLCanvas3D *canvas = plater->canvas3D();
    const std::string &name = m_data->volume_name;

    plater->take_snapshot(_L("Emboss text") + ": " + name);
    ModelVolume *volume = m_data->volume_ptr;
    if (volume == nullptr) {
        // decide to add as volume or new object  
        if (m_data->object_idx < 0) {
            // create new object
            app.obj_list()->load_mesh_object(tm, name, true, &m_data->text_configuration);
            app.mainframe->update_title();

            // TODO: find Why ???
            // load mesh cause close gizmo, on windows but not on linux
            // Open gizmo again when it is closed
            GLGizmosManager &mng = canvas->get_gizmos_manager();
            if (mng.get_current_type() != GLGizmosManager::Emboss)
                mng.open_gizmo(GLGizmosManager::Emboss);
            return;
        } else {
            // create new volume inside of object
            Model &model = plater->model();
            if (model.objects.size() <= m_data->object_idx) return;
            ModelObject *obj   = model.objects[m_data->object_idx];
            volume             = obj->add_volume(std::move(tm));
            // set a default extruder value, since user can't add it manually
            volume->config.set_key_value("extruder", new ConfigOptionInt(0));
        }
    } else {
        // update volume
        volume->set_mesh(std::move(tm));
        volume->set_new_unique_id();
        volume->calculate_convex_hull();
        volume->get_object()->invalidate_bounding_box();
    }

    volume->name = name;
    volume->text_configuration = m_data->text_configuration;

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    select_volume(volume);

    // Job promiss to refresh is not working
    canvas->reload_scene(true);
}

void EmbossJob::select_volume(ModelVolume *volume)
{    
    if (volume == nullptr) return;

    ObjectList * obj_list = wxGetApp().obj_list();
    
    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) {
        return vol == volume;
    };
    const Selection &selection = wxGetApp().plater()->canvas3D()->get_selection();
    wxDataViewItemArray sel =
        obj_list->reorder_volumes_and_get_selection(selection.get_object_idx(),
                                                    add_to_selection);

    if (!sel.IsEmpty()) obj_list->select_item(sel.front());
    obj_list->selection_changed();
}