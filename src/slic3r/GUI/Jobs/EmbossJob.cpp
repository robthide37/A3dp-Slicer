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
#include "slic3r/GUI/format.hpp"
#include "slic3r/Utils/UndoRedo.hpp"

using namespace Slic3r;
using namespace GUI;

// private namespace
namespace priv{
// <summary>
/// Create mesh from text
/// </summary>
/// <param name="text">Text to convert on mesh</param>
/// <param name="font">Define shape of characters.
/// NOTE: Can't be const cache glyphs</param>
/// <param name="font_prop">Property of font</param>
/// <param name="ctl">Control for job, check of cancelation</param>
/// <returns>Triangle mesh model</returns>
static TriangleMesh create_mesh(const char                *text,
                                Emboss::FontFileWithCache &font,
                                const FontProp            &font_prop,
                                GUI::Job::Ctl             &ctl);
/// <summary>
/// Create default mesh for embossed text
/// </summary>
/// <returns>Not empty model(index trinagle set - its)</returns>
static TriangleMesh create_default_mesh();
}


/////////////////
/// Update Volume
EmbossUpdateJob::EmbossUpdateJob(EmbossDataUpdate&& input)
    : m_input(std::move(input))
{}

void EmbossUpdateJob::process(Ctl &ctl)
{
    // check if exist valid font
    if (!m_input.font_file.has_value()) return;

    const TextConfiguration &cfg  = m_input.text_configuration;
    m_result = priv::create_mesh(cfg.text.c_str(), m_input.font_file,
                                 cfg.font_item.prop, ctl);
    if (m_result.its.empty()) return;    
    if (ctl.was_canceled()) return;

    // center triangle mesh
    Vec3d shift = m_result.bounding_box().center();
    m_result.translate(-shift.cast<float>());    
}

void EmbossUpdateJob::finalize(bool canceled, std::exception_ptr &)
{
    if (canceled || *m_input.cancel) return;

    // for sure that some object is created from shape
    if (m_result.its.indices.empty()) return;

    GUI_App &        app      = wxGetApp(); // may be move to input
    Plater *         plater   = app.plater();
    ObjectList *     obj_list = app.obj_list();
    GLCanvas3D *     canvas   = plater->canvas3D();
    GLGizmosManager &manager  = canvas->get_gizmos_manager();

    // Check emboss gizmo is still open
    if (manager.get_current_type() != GLGizmosManager::Emboss) return;

    std::string snap_name = GUI::format(_L("Text: %1%"), m_input.text_configuration.text);
    Plater::TakeSnapshot snapshot(plater, snap_name, UndoRedo::SnapshotType::GizmoAction);

    ModelVolume *volume = nullptr;
    Model &model = plater->model();
    for (auto obj : model.objects)
        for (auto vol : obj->volumes)
            if (vol->id() == m_input.volume_id) {
                volume = vol;
                break;
            }

    // could appear when user delete edited volume
    if (volume == nullptr)
        return;

    // update volume
    volume->set_mesh(std::move(m_result));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();
    volume->name               = m_input.volume_name;
    volume->text_configuration = m_input.text_configuration;

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


/////////////////
/// Create Volume
EmbossCreateVolumeJob::EmbossCreateVolumeJob(EmbossDataCreateVolume &&input)
    : m_input(std::move(input))
{}

void EmbossCreateVolumeJob::process(Ctl &ctl) {
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new emboss text object
    const char *text = m_input.text_configuration.text.c_str();
    FontProp   &prop = m_input.text_configuration.font_item.prop;

    m_result = priv::create_mesh(text, m_input.font_file, prop, ctl);
    if (m_result.its.empty()) m_result = priv::create_default_mesh();

    if (ctl.was_canceled()) return;
        
    // Create new volume inside of object
    const FontProp &font_prop = m_input.text_configuration.font_item.prop;
    Transform3d surface_trmat = Emboss::create_transformation_onto_surface(
            m_input.hit.position, m_input.hit.normal);
    Emboss::apply_transformation(font_prop, surface_trmat);
    m_transformation = m_input.hit_instance_tr.inverse() *
                       m_input.hit_object_tr * surface_trmat;    
}

void EmbossCreateVolumeJob::finalize(bool canceled, std::exception_ptr &) {
    if (canceled) return;

    GUI_App    &app      = wxGetApp();
    Plater     *plater   = app.plater();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = plater->canvas3D();
    Model &model = plater->model();

     // create volume in object
    size_t object_idx = m_input.object_idx;
    assert(model.objects.size() > object_idx);
    if (model.objects.size() <= object_idx) return;

    Plater::TakeSnapshot snapshot(plater, _L("Add Emboss text Volume"));

    ModelObject *obj    = model.objects[object_idx];
    ModelVolumeType type   = m_input.volume_type;
    ModelVolume *volume = obj->add_volume(std::move(m_result), type);

    // set a default extruder value, since user can't add it manually
    volume->config.set_key_value("extruder", new ConfigOptionInt(0));

    // do not allow model reload from disk
    volume->source.is_from_builtin_objects = true;

    volume->name               = m_input.volume_name;
    volume->text_configuration = std::move(m_input.text_configuration);
    volume->set_transformation(m_transformation);

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) {
        return vol == volume;
    };
    wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(
        (int) object_idx, add_to_selection);
    if (!sel.IsEmpty()) obj_list->select_item(sel.front());

    // update printable state on canvas
    if (type == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object(object_idx);

    obj_list->selection_changed();

    // WHY selection_changed set manipulation to world ???
    // so I set it back to local --> RotationGizmo need it
    ObjectManipulation *manipul = wxGetApp().obj_manipul();
    manipul->set_coordinates_type(ECoordinatesType::Local);

    // redraw scene
    canvas->reload_scene(true);
}


/////////////////
/// Create Object
EmbossCreateObjectJob::EmbossCreateObjectJob(EmbossDataCreateObject &&input)
    : m_input(std::move(input))
{}

void EmbossCreateObjectJob::process(Ctl &ctl) 
{
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new emboss text object
    const char *text = m_input.text_configuration.text.c_str();
    FontProp   &prop = m_input.text_configuration.font_item.prop;

    m_result = priv::create_mesh(text, m_input.font_file, prop, ctl);
    if (m_result.its.empty()) m_result = priv::create_default_mesh();

    if (ctl.was_canceled()) return;

    // Create new object
    // calculate X,Y offset position for lay on platter in place of
    // mouse click
    Vec2d bed_coor = CameraUtils::get_z0_position(
        m_input.camera, m_input.screen_coor);

    // check point is on build plate:
    Points bed_shape_;
    bed_shape_.reserve(m_input.bed_shape.size());
    for (const Vec2d &p : m_input.bed_shape)
        bed_shape_.emplace_back(p.cast<int>());
    Polygon bed(bed_shape_);
    if (!bed.contains(bed_coor.cast<int>()))
        // mouse pose is out of build plate so create object in center of plate
        bed_coor = bed.centroid().cast<double>();

    double z = m_input.text_configuration.font_item.prop.emboss / 2;
    Vec3d  offset(bed_coor.x(), bed_coor.y(), z);
    offset -= m_result.center();
    Transform3d::TranslationType tt(offset.x(), offset.y(), offset.z());
    m_transformation = Transform3d(tt);
}

void EmbossCreateObjectJob::finalize(bool canceled, std::exception_ptr &)
{
    if (canceled) return;

    GUI_App    &app      = wxGetApp();
    Plater     *plater   = app.plater();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = plater->canvas3D();

    Plater::TakeSnapshot snapshot(plater, _L("Add Emboss text object"));

    // Create new object and change selection
    bool center = false;
    obj_list->load_mesh_object(std::move(m_result), m_input.volume_name,
                                center, &m_input.text_configuration,
                                &m_transformation);

    // When add new object selection is empty.
    // When cursor move and no one object is selected than
    // Manager::reset_all() So Gizmo could be closed before end of creation object
    GLGizmosManager &manager = canvas->get_gizmos_manager();
    if (manager.get_current_type() != GLGizmosManager::Emboss)
        manager.open_gizmo(GLGizmosManager::Emboss);   

    // redraw scene
    canvas->reload_scene(true);
}


////////////////////////////
/// private namespace implementation
TriangleMesh priv::create_mesh(const char                *text,
                               Emboss::FontFileWithCache &font,
                               const FontProp            &font_prop,
                               GUI::Job::Ctl             &ctl)
{
    assert(font.has_value());
    if (!font.has_value()) return {};

    ExPolygons shapes = Emboss::text2shapes(font, text, font_prop);
    if (shapes.empty()) return {};
    if (ctl.was_canceled()) return {};

    int unit_per_em = font.font_file->unit_per_em;
    float scale    = font_prop.size_in_mm / unit_per_em;
    float depth    = font_prop.emboss / scale;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(depth);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    if (ctl.was_canceled()) return {};
    return TriangleMesh(Emboss::polygons2model(shapes, project));
}

TriangleMesh priv::create_default_mesh()
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