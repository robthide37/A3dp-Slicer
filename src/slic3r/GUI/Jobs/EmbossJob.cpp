#include "EmbossJob.hpp"

#include <libslic3r/Model.hpp>
#include <libslic3r/Format/OBJ.hpp> // load_obj for default mesh

#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/Gizmos/GLGizmoEmboss.hpp"
#include "slic3r/GUI/CameraUtils.hpp"

using namespace Slic3r;
using namespace GUI;

void EmbossUpdateJob::process(Ctl &ctl)
{
    // check if exist valid font
    if (m_input->font_file == nullptr) return;

    const TextConfiguration &cfg  = m_input->text_configuration;
    m_result = EmbossCreateJob::create_mesh(
        cfg.text.c_str(), *m_input->font_file, cfg.font_item.prop, ctl);
    if (m_result.its.empty()) return;    
    if (ctl.was_canceled()) return;

    // center triangle mesh
    Vec3d shift = m_result.bounding_box().center();
    m_result.translate(-shift.cast<float>());    
}

void EmbossUpdateJob::finalize(bool canceled, std::exception_ptr &)
{
    if (canceled) return;

    // for sure that some object is created from shape
    if (m_result.its.indices.empty()) return;

    GUI_App &        app      = wxGetApp(); // may be move to input
    Plater *         plater   = app.plater();
    ObjectList *     obj_list = app.obj_list();
    GLCanvas3D *     canvas   = plater->canvas3D();
    GLGizmosManager &manager  = canvas->get_gizmos_manager();

    // Check emboss gizmo is still open
    if (manager.get_current_type() != GLGizmosManager::Emboss) return;

    plater->take_snapshot(_L("Emboss text") + ": " + m_input->volume_name);

    ModelVolume *volume = m_input->volume;
    // find volume by object id - NOT WORK
    // -> edit text change volume id so could apper not found volume
    // ModelVolume *volume = nullptr;
    // Model &model = plater->model();
    // for (auto obj : model.objects)
    //    for (auto vol : obj->volumes)
    //        if (vol->id() == volume_id) {
    //            volume = vol;
    //            break;
    //        }
    // if (volume == nullptr) return;
    assert(volume != nullptr);

    // update volume
    volume->set_mesh(std::move(m_result));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();
    volume->name               = m_input->volume_name;
    volume->text_configuration = m_input->text_configuration;

    // update volume in right panel( volume / object name)
    const Selection &selection = canvas->get_selection();
    const GLVolume * gl_volume = selection.get_volume(
        *selection.get_volume_idxs().begin());
    int object_idx = gl_volume->object_idx();
    int volume_idx = gl_volume->volume_idx();
    obj_list->update_name_in_list(object_idx, volume_idx);

    // update printable state on canvas
    if (volume->type() == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object(
            (size_t) object_idx);

    // redraw scene
    canvas->reload_scene(true);
}

void EmbossCreateJob::process(Ctl &ctl) {
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new emboss text object
    m_result = (m_input->font_file == nullptr) ?
            create_default_mesh() :
            create_mesh(m_input->text_configuration.text.c_str(),
                        *m_input->font_file,
                        m_input->text_configuration.font_item.prop, ctl);
    if (m_result.its.empty()) m_result = create_default_mesh();
    if (ctl.was_canceled()) return;

    std::optional<RaycastManager::Hit> hit;
    if (m_input->object_idx.has_value()) {
        // By position of cursor create transformation to put text on surface of model
        ModelObject *obj = wxGetApp().plater()->model().objects[*m_input->object_idx];
        m_input->raycast_manager->actualize(obj);
        if (ctl.was_canceled()) return;
        hit = m_input->raycast_manager->unproject(m_input->screen_coor, m_input->camera);

        // context menu for add text could be open only by right click on an
        // object. After right click, object is selected and object_idx is set
        // also hit must exist. But there is proper behavior when hit doesn't
        // exists. When this assert appear distquish remove of it.
        assert(hit.has_value());
        if (!hit.has_value()) m_input->object_idx.reset();
    }
    if (!hit.has_value()) {
        // create new object
        // calculate X,Y offset position for lay on platter in place of
        // mouse click
        Vec2d bed_coor = CameraUtils::get_z0_position(m_input->camera, m_input->screen_coor);

        // check point is on build plate:
        Points  bed_shape_;
        bed_shape_.reserve(m_input->bed_shape.size());
        for (const Vec2d &p : m_input->bed_shape)
            bed_shape_.emplace_back(p.cast<int>());
        Polygon bed(bed_shape_);
        if (!bed.contains(bed_coor.cast<int>()))
            // mouse pose is out of build plate so create object in center of plate
            bed_coor = bed.centroid().cast<double>();

        double z = m_input->text_configuration.font_item.prop.emboss / 2;
        Vec3d  offset(bed_coor.x(), bed_coor.y(), z);
        offset -= m_result.center();
        Transform3d::TranslationType tt(offset.x(), offset.y(), offset.z());
        m_transformation = Transform3d(tt);
    } else {
        //m_transformation = Emboss::create_transformation_onto_surface(hit->position, hit->normal);
        assert(m_input->hit_vol_tr.has_value());
        if (m_input->hit_vol_tr.has_value()) {             
            Transform3d object_trmat = m_input->raycast_manager->get_transformation(hit->tr_key);
            const FontProp &font_prop = m_input->text_configuration.font_item.prop;
            Transform3d surface_trmat = Emboss::create_transformation_onto_surface(hit->position, hit->normal);
            Emboss::apply_transformation(font_prop, surface_trmat);
            m_transformation = m_input->hit_vol_tr->inverse() * object_trmat * surface_trmat;
        }        
    }
}

void EmbossCreateJob::finalize(bool canceled, std::exception_ptr &)
{
    if (canceled) return;
        
    GUI_App &   app      = wxGetApp();
    Plater *    plater   = app.plater();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = plater->canvas3D();

    // decide if create object or volume
    bool create_object = !m_input->object_idx.has_value();
    if (create_object) {
        plater->take_snapshot(_L("Add Emboss text object"));
        // Create new object and change selection
        bool center = false;
        obj_list->load_mesh_object(std::move(m_result), m_input->volume_name, center,
                                   &m_input->text_configuration, &m_transformation);

        // When add new object selection is empty.
        // When cursor move and no one object is selected than Manager::reset_all()
        // So Gizmo could be closed on end of creation object
        GLGizmosManager &manager = canvas->get_gizmos_manager();
        if (manager.get_current_type() != GLGizmosManager::Emboss)
            manager.open_gizmo(GLGizmosManager::Emboss);
    } else {
        // create volume in object
        size_t object_idx = *m_input->object_idx;        
        ModelVolumeType type = m_input->volume_type;

        // create new volume inside of object
        Model &model = plater->model();
        if (model.objects.size() <= object_idx) return;
        ModelObject *obj    = model.objects[object_idx];
        ModelVolume *volume = obj->add_volume(std::move(m_result));

        // set a default extruder value, since user can't add it manually
        volume->config.set_key_value("extruder", new ConfigOptionInt(0));

        // do not allow model reload from disk
        volume->source.is_from_builtin_objects = true;
        volume->set_type(type);
        volume->name               = m_input->volume_name;
        volume->text_configuration = m_input->text_configuration;
        volume->set_transformation(m_transformation);

        // update volume name in object list
        // updata selection after new volume added
        // change name of volume in right panel
        // select only actual volume
        // when new volume is created change selection to this volume
        auto add_to_selection = [volume](const ModelVolume *vol) {
            return vol == volume;
        };
        wxDataViewItemArray sel =
            obj_list->reorder_volumes_and_get_selection((int) object_idx,
                                                        add_to_selection);
        if (!sel.IsEmpty()) obj_list->select_item(sel.front());

        // update printable state on canvas
        if (type == ModelVolumeType::MODEL_PART)
            canvas->update_instance_printable_state_for_object(object_idx);

        obj_list->selection_changed();

        // WHY selection_changed set manipulation to world ???
        // so I set it back to local --> RotationGizmo need it
        ObjectManipulation *manipul = wxGetApp().obj_manipul();
        manipul->set_coordinates_type(ECoordinatesType::Local);
    }
    // redraw scene
    canvas->reload_scene(true);
}

TriangleMesh EmbossCreateJob::create_default_mesh()
{
    // When cant load any font use default object loaded from file
    std::string  path = Slic3r::resources_dir() + "/data/embossed_text.stl";
    TriangleMesh triangle_mesh;
    if (!load_obj(path.c_str(), &triangle_mesh)) {
        // when can't load mesh use cube
        return TriangleMesh(its_make_cube(36., 4., 2.5));
    }
    return triangle_mesh;
}

TriangleMesh EmbossCreateJob::create_mesh(const char *      text,
                                          Emboss::FontFile &font,
                                          const FontProp &  font_prop,
                                          Ctl &             ctl)
{
    ExPolygons shapes = Emboss::text2shapes(font, text, font_prop);
    if (shapes.empty()) return {};
    if (ctl.was_canceled()) return {};

    float scale    = font_prop.size_in_mm / font.unit_per_em;
    float depth    = font_prop.emboss / scale;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(depth);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    if (ctl.was_canceled()) return {};
    return TriangleMesh(Emboss::polygons2model(shapes, project));
}