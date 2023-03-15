#include "EmbossJob.hpp"

#include <stdexcept>

#include <libslic3r/Model.hpp>
#include <libslic3r/Format/OBJ.hpp> // load_obj for default mesh
#include <libslic3r/CutSurface.hpp> // use surface cuts
#include <libslic3r/BuildVolume.hpp> // create object

#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Gizmos/GLGizmoEmboss.hpp"
#include "slic3r/GUI/Selection.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "slic3r/GUI/format.hpp"
#include "slic3r/GUI/3DScene.hpp"
#include "slic3r/GUI/Jobs/Worker.hpp" 
#include "slic3r/Utils/UndoRedo.hpp"
#include "slic3r/Utils/RaycastManager.hpp"

using namespace Slic3r;
using namespace Slic3r::Emboss;
using namespace Slic3r::GUI;
using namespace Slic3r::GUI::Emboss;

// Private implementation for create volume and objects jobs
namespace Slic3r::GUI::Emboss {
/// <summary>
/// Hold neccessary data to create ModelVolume in job
/// Volume is created on the surface of existing volume in object.
/// NOTE: EmbossDataBase::font_file doesn't have to be valid !!!
/// </summary>
struct DataCreateVolume
{
    // Hold data about shape
    DataBasePtr base;

    // define embossed volume type
    ModelVolumeType volume_type;

    // parent ModelObject index where to create volume
    ObjectID object_id;

    // new created volume transformation
    std::optional<Transform3d> trmat;

    // Define which gizmo open on the success
    GLGizmosManager::EType gizmo;
};

/// <summary>
/// Create new TextVolume on the surface of ModelObject
/// Should not be stopped
/// NOTE: EmbossDataBase::font_file doesn't have to be valid !!!
/// </summary>
class CreateVolumeJob : public Job
{
    DataCreateVolume m_input;
    TriangleMesh     m_result;

public:
    explicit CreateVolumeJob(DataCreateVolume &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to create ModelObject in job
/// Object is placed on bed under screen coor
/// OR to center of scene when it is out of bed shape
/// </summary>
struct DataCreateObject
{
    // Hold data about shape
    DataBasePtr base;

    // define position on screen where to create object
    Vec2d screen_coor;

    // projection property
    Camera camera;

    // shape of bed in case of create volume on bed
    std::vector<Vec2d> bed_shape;

    // Define which gizmo open on the success
    GLGizmosManager::EType gizmo;
};

/// <summary>
/// Create new TextObject on the platter
/// Should not be stopped
/// </summary>
class CreateObjectJob : public Job
{
    DataCreateObject m_input;
    TriangleMesh     m_result;
    Transform3d      m_transformation;

public:
    explicit CreateObjectJob(DataCreateObject &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to create(cut) volume from surface object in job
/// </summary>
struct CreateSurfaceVolumeData : public SurfaceVolumeData
{
    // Hold data about shape
    DataBasePtr base;

    // define embossed volume type
    ModelVolumeType volume_type;

    // parent ModelObject index where to create volume
    ObjectID object_id;

    // Define which gizmo open on the success
    GLGizmosManager::EType gizmo;
};

/// <summary>
/// Cut surface from object and create cutted volume
/// Should not be stopped
/// </summary>
class CreateSurfaceVolumeJob : public Job
{
    CreateSurfaceVolumeData m_input;
    TriangleMesh            m_result;

public:
    explicit CreateSurfaceVolumeJob(CreateSurfaceVolumeData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};
} // namespace Slic3r::GUI::Emboss

// private namespace
namespace priv{
// create sure that emboss object is bigger than source object [in mm]
constexpr float safe_extension = 1.0f;

/// <summary>
/// Assert check of inputs data
/// </summary>
/// <param name="input"></param>
/// <returns></returns>
static bool check(const DataBase &input, bool check_fontfile = true, bool use_surface = false);
static bool check(GLGizmosManager::EType gizmo);
static bool check(const DataCreateVolume &input, bool is_main_thread = false);
static bool check(const DataCreateObject &input);
static bool check(const DataUpdate &input, bool is_main_thread = false, bool use_surface = false);
static bool check(const CreateSurfaceVolumeData &input, bool is_main_thread = false);
static bool check(const UpdateSurfaceVolumeData &input, bool is_main_thread = false);

// <summary>
/// Try to create mesh from text
/// </summary>
/// <param name="input">Text to convert on mesh
/// + Shape of characters + Property of font</param>
/// <param name="font">Font file with cache
/// NOTE: Cache glyphs is changed</param>
/// <param name="was_canceled">To check if process was canceled</param>
/// <returns>Triangle mesh model</returns>
template<typename Fnc> static TriangleMesh try_create_mesh(DataBase &input, const Fnc& was_canceled);
template<typename Fnc> static TriangleMesh create_mesh(DataBase &input, const Fnc& was_canceled, Job::Ctl &ctl);

/// <summary>
/// Create default mesh for embossed text
/// </summary>
/// <returns>Not empty model(index trinagle set - its)</returns>
static TriangleMesh create_default_mesh();

/// <summary>
/// Must be called on main thread
/// </summary>
/// <param name="mesh">New mesh data</param>
/// <param name="data">Text configuration, ...</param>
/// <param name="mesh">Transformation of volume</param>
static void update_volume(TriangleMesh &&mesh, const DataUpdate &data, Transform3d *tr = nullptr);

/// <summary>
/// Update name in right panel
/// </summary>
/// <param name="obj_list">Right panel data</param>
/// <param name="volume">Volume with just changed name</param>
static void update_name_in_list(const ObjectList &object_list, const ModelVolume &volume);

/// <summary>
/// Add new volume to object
/// </summary>
/// <param name="mesh">triangles of new volume</param>
/// <param name="object_id">Object where to add volume</param>
/// <param name="type">Type of new volume</param>
/// <param name="trmat">Transformation of volume inside of object</param>
/// <param name="data">Text configuration and New VolumeName</param>
/// <param name="gizmo">Gizmo to open</param>
static void create_volume(TriangleMesh &&mesh, const ObjectID& object_id, const ModelVolumeType type, 
    const std::optional<Transform3d>& trmat, const DataBase &data, GLGizmosManager::EType gizmo);

/// <summary>
/// Select Volume from objects
/// </summary>
/// <param name="objects">All objects in scene</param>
/// <param name="volume_id">Identifier of volume in object</param>
/// <returns>Pointer to volume when exist otherwise nullptr</returns>
static ModelVolume *get_volume(ModelObjectPtrs &objects, const ObjectID &volume_id);

/// <summary>
/// Create projection for cut surface from mesh
/// </summary>
/// <param name="tr">Volume transformation in object</param>
/// <param name="shape_scale">Convert shape to milimeters</param>
/// <param name="z_range">Bounding box 3d of model volume for projection ranges</param> 
/// <returns>Orthogonal cut_projection</returns>
static OrthoProject create_projection_for_cut(Transform3d tr, double shape_scale, const std::pair<float, float> &z_range);

/// <summary>
/// Create tranformation for emboss Cutted surface
/// </summary>
/// <param name="is_outside">True .. raise, False .. engrave</param>
/// <param name="emboss">Depth of embossing</param>
/// <param name="tr">Text voliume transformation inside object</param>
/// <param name="cut">Cutted surface from model</param>
/// <returns>Projection</returns>
static OrthoProject3d create_emboss_projection(bool is_outside, float emboss, Transform3d tr, SurfaceCut &cut);

/// <summary>
/// Cut surface into triangle mesh
/// </summary>
/// <param name="base">(can't be const - cache of font)</param>
/// <param name="input2">SurfaceVolume data</param>
/// <param name="was_canceled">Check to interupt execution</param>
/// <returns>Extruded object from cuted surace</returns>
static TriangleMesh cut_surface(/*const*/ DataBase &input1, const SurfaceVolumeData &input2, const std::function<bool()>& was_canceled);

/// <summary>
/// Copied triangles from object to be able create mesh for cut surface from
/// </summary>
/// <param name="volumes">Source object volumes for cut surface from</param>
/// <param name="text_volume_id">Source volume id</param>
/// <returns>Source data for cut surface from</returns>
static SurfaceVolumeData::ModelSources create_sources(const ModelVolumePtrs &volumes, std::optional<size_t> text_volume_id = {});

static void create_message(const std::string &message); // only in finalize
static bool process(std::exception_ptr &eptr);
static bool finalize(bool canceled, std::exception_ptr &eptr, const DataBase &input);

class JobException : public std::runtime_error { 
public: JobException(const char* message):runtime_error(message){}};
static auto was_canceled(Job::Ctl &ctl, DataBase &base){
    return [&ctl, &cancel = base.cancel]() -> bool {
        if (cancel->load())
            return true;
        return ctl.was_canceled();
    };
}

}// namespace priv

void Slic3r::GUI::Emboss::DataBase::write(ModelVolume &volume) const{
    volume.name         = volume_name;
    volume.emboss_shape = shape;
}

/////////////////
/// Create Volume
CreateVolumeJob::CreateVolumeJob(DataCreateVolume &&input): m_input(std::move(input)){ assert(priv::check(m_input, true)); }

void CreateVolumeJob::process(Ctl &ctl) {
    if (!priv::check(m_input)) 
        throw std::runtime_error("Bad input data for EmbossCreateVolumeJob.");

    m_result = priv::create_mesh(*m_input.base, priv::was_canceled(ctl, *m_input.base), ctl);
    // center result
    Vec3f c = m_result.bounding_box().center().cast<float>();
    if (!c.isApprox(Vec3f::Zero())) m_result.translate(-c);
}
void CreateVolumeJob::finalize(bool canceled, std::exception_ptr &eptr) {
    if (!priv::finalize(canceled, eptr, *m_input.base))
        return;
    if (m_result.its.empty()) 
        return priv::create_message(_u8L("Can't create empty volume."));
    priv::create_volume(std::move(m_result), m_input.object_id, m_input.volume_type, m_input.trmat, *m_input.base, m_input.gizmo);
}


/////////////////
/// Create Object
CreateObjectJob::CreateObjectJob(DataCreateObject &&input): m_input(std::move(input)){ assert(priv::check(m_input)); }

void CreateObjectJob::process(Ctl &ctl) 
{
    if (!priv::check(m_input))
        throw std::runtime_error("Bad input data for EmbossCreateObjectJob.");

    if (m_input.base->shape.distance.has_value())
        m_input.base->shape.distance.reset();

    // can't create new object with using surface
    if (m_input.base->shape.use_surface)
        m_input.base->shape.use_surface = false;

    auto was_canceled = priv::was_canceled(ctl, *m_input.base);
    m_result = priv::create_mesh(*m_input.base, was_canceled, ctl);
    if (was_canceled()) return;

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
    Slic3r::Polygon bed(bed_shape_);
    if (!bed.contains(bed_coor.cast<int>()))
        // mouse pose is out of build plate so create object in center of plate
        bed_coor = bed.centroid().cast<double>();

    // TODO: need TextConfiguration refactor to work !!!
    double z = m_input.base->shape.depth / 2;

    Vec3d  offset(bed_coor.x(), bed_coor.y(), z);
    offset -= m_result.center();
    Transform3d::TranslationType tt(offset.x(), offset.y(), offset.z());
    m_transformation = Transform3d(tt);
}

void CreateObjectJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    if (!priv::finalize(canceled, eptr, *m_input.base))
        return;

    // only for sure
    if (m_result.empty()) 
        return priv::create_message(_u8L("Can't create empty object."));

    GUI_App &app    = wxGetApp();
    Plater  *plater = app.plater();
    plater->take_snapshot(_L("Add Emboss text object"));

    Model& model = plater->model();
#ifdef _DEBUG
    check_model_ids_validity(model);
#endif /* _DEBUG */
    {
        // INFO: inspiration for create object is from ObjectList::load_mesh_object()
        ModelObject *new_object = model.add_object();
        new_object->name = m_input.base->volume_name;
        new_object->add_instance(); // each object should have at list one instance

        ModelVolume *new_volume = new_object->add_volume(std::move(m_result));
        // set a default extruder value, since user can't add it manually
        new_volume->config.set_key_value("extruder", new ConfigOptionInt(0));
        // write emboss data into volume
        m_input.base->write(*new_volume);

        // set transformation
        Slic3r::Geometry::Transformation tr(m_transformation);
        new_object->instances.front()->set_transformation(tr);
        new_object->ensure_on_bed();

        // Actualize right panel and set inside of selection
        app.obj_list()->paste_objects_into_list({model.objects.size() - 1});
    }
#ifdef _DEBUG
    check_model_ids_validity(model);
#endif /* _DEBUG */

    // When add new object selection is empty.
    // When cursor move and no one object is selected than
    // Manager::reset_all() So Gizmo could be closed before end of creation object
    GLCanvas3D      *canvas  = plater->canvas3D();
    GLGizmosManager &manager = canvas->get_gizmos_manager();
    if (manager.get_current_type() != m_input.gizmo)
        manager.open_gizmo(m_input.gizmo);

    // redraw scene
    canvas->reload_scene(true);
}

/////////////////
/// Update Volume
UpdateJob::UpdateJob(DataUpdate&& input): m_input(std::move(input)){ assert(priv::check(m_input, true)); }

void UpdateJob::process(Ctl &ctl)
{
    if (!priv::check(m_input))
        throw std::runtime_error("Bad input data for EmbossUpdateJob.");

    auto was_canceled = priv::was_canceled(ctl, *m_input.base);
    m_result = priv::try_create_mesh(*m_input.base, was_canceled);
    if (was_canceled()) return;
    if (m_result.its.empty())
        throw priv::JobException(_u8L("Created text volume is empty. Change text or font.").c_str());

    // center triangle mesh
    Vec3d shift = m_result.bounding_box().center();
    m_result.translate(-shift.cast<float>());    
}

void UpdateJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    if (!priv::finalize(canceled, eptr, *m_input.base))
        return;
    priv::update_volume(std::move(m_result), m_input);    
}

/////////////////
/// Create Surface volume
CreateSurfaceVolumeJob::CreateSurfaceVolumeJob(CreateSurfaceVolumeData &&input) 
    : m_input(std::move(input))
{
    assert(priv::check(m_input, true));
}

void CreateSurfaceVolumeJob::process(Ctl &ctl) {
    if (!priv::check(m_input)) 
        throw std::runtime_error("Bad input data for CreateSurfaceVolumeJob.");
    m_result = priv::cut_surface(*m_input.base, m_input, priv::was_canceled(ctl, *m_input.base));
}

void CreateSurfaceVolumeJob::finalize(bool canceled, std::exception_ptr &eptr) {
    if (!priv::finalize(canceled, eptr, *m_input.base))
        return; 
    priv::create_volume(std::move(m_result), m_input.object_id,
        m_input.volume_type, m_input.transform, *m_input.base, m_input.gizmo);
}

/////////////////
/// Cut Surface
UpdateSurfaceVolumeJob::UpdateSurfaceVolumeJob(UpdateSurfaceVolumeData &&input)
    : m_input(std::move(input))
{
    assert(priv::check(m_input, true));
}

void UpdateSurfaceVolumeJob::process(Ctl &ctl)
{
    if (!priv::check(m_input)) 
        throw std::runtime_error("Bad input data for UseSurfaceJob.");
    m_result = priv::cut_surface(*m_input.base, m_input, priv::was_canceled(ctl, *m_input.base));
}

void UpdateSurfaceVolumeJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    if (!priv::finalize(canceled, eptr, *m_input.base))
        return;

    // when start using surface it is wanted to move text origin on surface of model
    // also when repeteadly move above surface result position should match
    Transform3d *tr = &m_input.transform;
    priv::update_volume(std::move(m_result), m_input, tr);
}

namespace priv {
/// <summary>
/// Check if volume type is possible use for new text volume
/// </summary>
/// <param name="volume_type">Type</param>
/// <returns>True when allowed otherwise false</returns>
static bool is_valid(ModelVolumeType volume_type);

/// <summary>
/// Start job for add new volume to object with given transformation
/// </summary>
/// <param name="worker">Define where to queue the job. e.g. wxGetApp().plater()->get_ui_job_worker()</param>
/// <param name="object">Define where to add</param>
/// <param name="volume_tr">Wanted volume transformation, when not set will be calculated after creation to be near the object</param>
/// <param name="data">Define what to emboss - shape</param>
/// <param name="volume_type">Type of volume: Part, negative, modifier</param>
/// <param name="gizmo">Define which gizmo open on the success</param>
/// <returns>Nullptr when job is sucessfully add to worker otherwise return data to be processed different way</returns>
static bool start_create_volume_job(Worker &worker, const ModelObject &object, const std::optional<Transform3d>& volume_tr, DataBasePtr data, ModelVolumeType volume_type, GLGizmosManager::EType gizmo);

/// <summary>
/// Find volume in selected objects with closest convex hull to screen center.
/// </summary>
/// <param name="selection">Define where to search for closest</param>
/// <param name="screen_center">Canvas center(dependent on camera settings)</param>
/// <param name="objects">Actual objects</param>
/// <param name="closest_center">OUT: coordinate of controid of closest volume</param>
/// <returns>closest volume when exists otherwise nullptr</returns>
static const GLVolume *find_closest(
    const Selection &selection, const Vec2d &screen_center, 
    const Camera &camera, const ModelObjectPtrs &objects, Vec2d *closest_center);

/// <summary>
/// Start job for add object with text into scene
/// </summary>
/// <param name="plater">Contain worker and build shape</param>
/// <param name="emboss_data">Define params of text</param>
/// <param name="coor">Screen coordinat, where to create new object laying on bed</param>
/// <param name="gizmo">Define which gizmo open on the success</param>
/// <returns>True when can add job to worker otherwise FALSE</returns>
static bool start_create_object_job(Plater &plater, DataBasePtr emboss_data, const Vec2d &coor, GLGizmosManager::EType gizmo);

/// <summary>
/// Start job to create volume on the surface of object
/// </summary>
/// <param name="plater">scene RayCasters + Objects + Camera + worker</param>
/// <param name="data">Describe what to emboss</param>
/// <param name="volume_type">Type of new volume</param>
/// <param name="screen_coor">Where to add</param>
/// <param name="gl_volume">Surface point is belonge to</param>
/// <param name="raycaster">For project on surface</param>
/// <param name="gizmo">Define which gizmo open on the success</param>
/// <param name="distance">Distance from surface</param>
/// <param name="angle">Angle around emboss direction</param>
/// <param name="success">True on add job to worker otherwise false</param>
/// <returns>Nullptr when job is sucessfully add to worker otherwise return data to be processed different way</returns>
static DataBasePtr start_create_volume_on_surface_job(Plater                     &plater,
                                                      DataBasePtr                 data,
                                                      ModelVolumeType             volume_type,
                                                      const Vec2d                &screen_coor,
                                                      const GLVolume             &gl_volume,
                                                      RaycastManager           &raycaster,
                                                      GLGizmosManager::EType      gizmo,
                                                      const std::optional<float> &distance,
                                                      const std::optional<float> &angle,
                                                      bool                       &success);

}

namespace Slic3r::GUI::Emboss {

SurfaceVolumeData::ModelSources create_volume_sources(const ModelVolume &text_volume)
{
    const ModelVolumePtrs &volumes = text_volume.get_object()->volumes;
    // no other volume in object
    if (volumes.size() <= 1)
        return {};
    return priv::create_sources(volumes, text_volume.id().id);
}

bool start_create_volume(Plater                     *plater_ptr,
                         DataBasePtr                 data,
                         ModelVolumeType             volume_type,
                         RaycastManager             &raycaster,
                         unsigned char               gizmo,
                         const Vec2d                &mouse_pos,
                         const std::optional<float> &distance,
                         const std::optional<float> &angle)
{
    if (data == nullptr)
        return false;
    if (!priv::is_valid(volume_type))
        return false;

    assert(plater_ptr);
    if (plater_ptr == nullptr)
        return false;
    Plater &plater = *plater_ptr;

    const GLCanvas3D *canvas_ptr = plater.get_current_canvas3D();
    assert(canvas_ptr);
    if (canvas_ptr == nullptr)
        return false;

    auto gizmo_type = static_cast<GLGizmosManager::EType>(gizmo);
    const GLVolume *gl_volume = get_first_hovered_gl_volume(*canvas_ptr);
    if (gl_volume == nullptr)
        // object is not under mouse position soo create object on plater
        return priv::start_create_object_job(plater, std::move(data), mouse_pos, gizmo_type);

    bool success = true;
    DataBasePtr data2 = priv::start_create_volume_on_surface_job(plater, std::move(data), 
        volume_type, mouse_pos, *gl_volume, raycaster, gizmo_type, distance, angle, success);

    // Is successfull created job for add volume on surface?
    if (success)
        return true;

    // not success and consume data
    if (data2 == nullptr)
        return false;
    
    // Can't create on coordinate try to create somewhere
    return start_create_volume_without_position(plater_ptr, std::move(data2), volume_type, raycaster, gizmo, distance, angle);
}

bool start_create_volume_without_position(Plater                     *plater_ptr,
                                          DataBasePtr                 data,
                                          ModelVolumeType             volume_type,
                                          RaycastManager             &raycaster,
                                          unsigned char               gizmo,
                                          const std::optional<float> &distance,
                                          const std::optional<float> &angle)
{
    if (!priv::is_valid(volume_type))
        return false;

    assert(plater_ptr);
    if (plater_ptr == nullptr)
        return false;
    Plater &plater = *plater_ptr;

    GLCanvas3D *canvas_ptr = plater.get_current_canvas3D();
    assert(canvas_ptr);
    if (canvas_ptr == nullptr)
        return false;
    const GLCanvas3D &canvas = *canvas_ptr;

    // select position by camera position and view direction
    const Selection &selection = canvas.get_selection();
    int object_idx = selection.get_object_idx();

    Size s = canvas.get_canvas_size();
    Vec2d screen_center(s.get_width() / 2., s.get_height() / 2.);
    const ModelObjectPtrs &objects = selection.get_model()->objects;

    GLGizmosManager::EType gizmo_type = static_cast<GLGizmosManager::EType>(gizmo);

    // No selected object so create new object
    if (selection.is_empty() || object_idx < 0 || 
        static_cast<size_t>(object_idx) >= objects.size()) 
        // create Object on center of screen
        // when ray throw center of screen not hit bed it create object on center of bed
        return priv::start_create_object_job(plater, std::move(data), screen_center, gizmo_type);

    // create volume inside of selected object
    Vec2d coor;
    const Camera &camera = wxGetApp().plater()->get_camera();
    const GLVolume *gl_volume = priv::find_closest(selection, screen_center, camera, objects, &coor);

    if (gl_volume == nullptr)
        return priv::start_create_object_job(plater, std::move(data), screen_center, gizmo_type);
    
    bool success = true;
    DataBasePtr data2 = priv::start_create_volume_on_surface_job(plater, std::move(data), 
        volume_type, coor, *gl_volume, raycaster, gizmo_type, distance, angle, success);

    // Is successfull created job for add volume on surface?
    if (success)
        return true;

    // not success and consume data
    if (data2 == nullptr)
        return false;

    // In centroid of convex hull is not hit with object. e.g. torid
    // soo create transfomation on border of object

    // there is no point on surface so no use of surface will be applied
    if (data2->shape.use_surface)
        data2->shape.use_surface = false;
        
    Worker &worker = plater.get_ui_job_worker();
    const ModelObject *object = get_model_object(*gl_volume, objects);
    if (object == nullptr)
        return false;

    return priv::start_create_volume_job(worker, *object, {}, std::move(data2), volume_type, gizmo_type);
}

} // namespace Slic3r::GUI::Emboss

////////////////////////////
/// private namespace implementation
bool priv::check(const DataBase &input, bool check_fontfile, bool use_surface)
{
    bool res = true;
    //if (check_fontfile) {
    //    assert(input.font_file.has_value());
    //    res &= input.font_file.has_value();
    //}
    //assert(!input.text_configuration.fix_3mf_tr.has_value());
    //res &= !input.text_configuration.fix_3mf_tr.has_value();
    //assert(!input.text_configuration.text.empty());
    //res &= !input.text_configuration.text.empty();
    assert(!input.volume_name.empty());
    res &= !input.volume_name.empty();
    //assert(input.text_configuration.style.prop.use_surface == use_surface);
    //res &= input.text_configuration.style.prop.use_surface == use_surface;
    return res; 
}

bool priv::check(GLGizmosManager::EType gizmo)
{
    assert(gizmo == GLGizmosManager::Emboss || gizmo == GLGizmosManager::Svg);
    return gizmo == GLGizmosManager::Emboss || gizmo == GLGizmosManager::Svg;
}

bool priv::check(const DataCreateVolume &input, bool is_main_thread) {
    bool check_fontfile = false;
    assert(input.base != nullptr);
    bool res = input.base != nullptr;
    res &= check(*input.base, check_fontfile);
    assert(input.volume_type != ModelVolumeType::INVALID);
    res &= input.volume_type != ModelVolumeType::INVALID;
    assert(input.object_id.id >= 0);
    res &= input.object_id.id >= 0;
    res &= check(input.gizmo);
    return res; 
}
bool priv::check(const DataCreateObject &input) {
    bool check_fontfile = false;
    assert(input.base != nullptr);
    bool res = input.base != nullptr;
    res &= check(*input.base, check_fontfile);
    assert(input.screen_coor.x() >= 0.);
    res &= input.screen_coor.x() >= 0.;
    assert(input.screen_coor.y() >= 0.);
    res &= input.screen_coor.y() >= 0.;
    assert(input.bed_shape.size() >= 3); // at least triangle
    res &= input.bed_shape.size() >= 3;
    res &= check(input.gizmo);
    return res;
}
bool priv::check(const DataUpdate &input, bool is_main_thread, bool use_surface){
    bool check_fontfile = true;
    assert(input.base != nullptr);
    bool res = input.base != nullptr;
    res &= check(*input.base, check_fontfile, use_surface);
    assert(input.volume_id.id >= 0);
    res &= input.volume_id.id >= 0;
    if (is_main_thread)
        assert(get_volume(wxGetApp().model().objects, input.volume_id) != nullptr);
    assert(input.base->cancel != nullptr);
    res &= input.base->cancel != nullptr;
    if (is_main_thread)
        assert(!input.base->cancel->load());
    return res;
}
bool priv::check(const CreateSurfaceVolumeData &input, bool is_main_thread)
{
    bool use_surface = true;
    assert(input.base != nullptr);
    bool res = input.base != nullptr;
    res &= check(*input.base, is_main_thread, use_surface);
    assert(!input.sources.empty());
    res &= !input.sources.empty();
    res &= check(input.gizmo);
    return res;
}
bool priv::check(const UpdateSurfaceVolumeData &input, bool is_main_thread){
    bool use_surface = true;
    assert(input.base != nullptr);
    bool res = input.base != nullptr;
    res &= check(*input.base, is_main_thread, use_surface);
    assert(!input.sources.empty());
    res &= !input.sources.empty();
    return res;
}

template<typename Fnc>
TriangleMesh priv::try_create_mesh(DataBase &base, const Fnc& was_canceled)
{
    const EmbossShape& shape = base.create_shape();
    if (shape.shapes.empty()) return {};  
    double depth = shape.depth / shape.scale;
    auto  projectZ = std::make_unique<ProjectZ>(depth);
    ProjectScale project(std::move(projectZ), shape.scale);
    if (was_canceled()) return {};
    return TriangleMesh(polygons2model(shape.shapes, project)); 
}

template<typename Fnc>
TriangleMesh priv::create_mesh(DataBase &input, const Fnc& was_canceled, Job::Ctl& ctl)
{
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new emboss text object
    TriangleMesh result = try_create_mesh(input, was_canceled);
    if (was_canceled()) return {};

    if (result.its.empty()) {
        result = priv::create_default_mesh();
        if (was_canceled()) return {};
        // only info
        ctl.call_on_main_thread([]() {
            create_message(_u8L("It is used default volume for embossed "
                                "text, try to change text or font to fix it."));
        });
    }

    assert(!result.its.empty());
    return result;
}

TriangleMesh priv::create_default_mesh()
{
    // When cant load any font use default object loaded from file
    std::string  path = Slic3r::resources_dir() + "/data/embossed_text.obj";
    TriangleMesh triangle_mesh;
    if (!load_obj(path.c_str(), &triangle_mesh)) {
        // when can't load mesh use cube
        return TriangleMesh(its_make_cube(36., 4., 2.5));
    }
    return triangle_mesh;
}

void UpdateJob::update_volume(ModelVolume *volume, TriangleMesh &&mesh, const DataBase &base)
{
    // check inputs
    bool is_valid_input = 
        volume != nullptr &&
        !mesh.empty() && 
        !base.volume_name.empty();
    assert(is_valid_input);
    if (!is_valid_input) return;

    // update volume
    volume->set_mesh(std::move(mesh));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();

    // write data from base into volume
    base.write(*volume);
        
    GUI_App &app = wxGetApp(); // may be move to input
    if (volume->name != base.volume_name) {
        volume->name = base.volume_name;
        
        ObjectList *obj_list = app.obj_list();
        if (obj_list != nullptr)
            priv::update_name_in_list(*obj_list, *volume);
    }

    // When text is object.
    // When text positive volume is lowest part of object than modification of text 
    // have to move object on bed.
    if (volume->type() == ModelVolumeType::MODEL_PART)
        volume->get_object()->ensure_on_bed();

    // redraw scene
    GLCanvas3D *canvas = app.plater()->canvas3D();

    bool refresh_immediately = false;
    canvas->reload_scene(refresh_immediately);

    // Change buttons "Export G-code" into "Slice now"
    canvas->post_event(SimpleEvent(EVT_GLCANVAS_SCHEDULE_BACKGROUND_PROCESS));
}

void priv::update_name_in_list(const ObjectList &object_list, const ModelVolume &volume)
{
    const ModelObjectPtrs *objects_ptr = object_list.objects();
    if (objects_ptr == nullptr)
        return;

    const ModelObjectPtrs &objects = *objects_ptr;
    const ModelObject *object = volume.get_object();
    const ObjectID& object_id = object->id();

    // search for index of object
    int object_index = -1;
    for (size_t i = 0; i < objects.size(); ++i) 
        if (objects[i]->id() == object_id) {
            object_index = i;
            break;
        }

    const ModelVolumePtrs volumes = object->volumes;
    const ObjectID& volume_id = volume.id();

    // search for index of volume
    int volume_index = -1;
    for (size_t i = 0; i < volumes.size(); ++i) 
        if (volumes[i]->id() == volume_id) {
            volume_index = i;
            break;
        }

    if (object_index < 0 || volume_index < 0)
        return;

    object_list.update_name_in_list(object_index, volume_index);
}

void priv::update_volume(TriangleMesh &&mesh, const DataUpdate &data, Transform3d* tr)
{
    // for sure that some object will be created
    if (mesh.its.empty()) 
        return create_message("Empty mesh can't be created.");

    Plater *plater = wxGetApp().plater();
    // Check gizmo is still open otherwise job should be canceled
    assert(plater->canvas3D()->get_gizmos_manager().get_current_type() == GLGizmosManager::Emboss || 
           plater->canvas3D()->get_gizmos_manager().get_current_type() == GLGizmosManager::Svg);

    std::string snap_name = GUI::format(_L("Change: %1%"), data.base->volume_name);
    Plater::TakeSnapshot snapshot(plater, snap_name, UndoRedo::SnapshotType::GizmoAction);
    ModelVolume *volume = get_volume(plater->model().objects, data.volume_id);

    // could appear when user delete edited volume
    if (volume == nullptr)
        return;

    if (tr) {
        volume->set_transformation(*tr);
    } else {
        // apply fix matrix made by store to .3mf
        const auto &tc = volume->text_configuration;
        assert(tc.has_value());
        if (tc.has_value() && tc->fix_3mf_tr.has_value())
            volume->set_transformation(volume->get_matrix() * tc->fix_3mf_tr->inverse());
    }

    UpdateJob::update_volume(volume, std::move(mesh), *data.base);
}

void priv::create_volume(TriangleMesh                    &&mesh,
                         const ObjectID                   &object_id,
                         const ModelVolumeType             type,
                         const std::optional<Transform3d> &trmat,
                         const DataBase                   &data,
                         GLGizmosManager::EType            gizmo)
{
    GUI_App         &app      = wxGetApp();
    Plater          *plater   = app.plater();
    ObjectList      *obj_list = app.obj_list();
    GLCanvas3D      *canvas   = plater->canvas3D();
    ModelObjectPtrs &objects  = plater->model().objects;

    ModelObject *obj = nullptr;
    size_t object_idx = 0;
    for (; object_idx < objects.size(); ++object_idx) {
        ModelObject *o = objects[object_idx];
        if (o->id() == object_id) { 
            obj = o;
            break;
        }   
    }

    // Parent object for text volume was propably removed.
    // Assumption: User know what he does, so text volume is no more needed.
    if (obj == nullptr) 
        return priv::create_message(_u8L("Bad object to create volume."));

    if (mesh.its.empty()) 
        return priv::create_message(_u8L("Can't create empty volume."));

    plater->take_snapshot(_L("Add Emboss text Volume"));

    BoundingBoxf3 instance_bb;
    if (!trmat.has_value()) {
        // used for align to instance
        size_t instance_index = 0; // must exist
        instance_bb = obj->instance_bounding_box(instance_index);
    }

    // NOTE: be carefull add volume also center mesh !!!
    // So first add simple shape(convex hull is also calculated)
    ModelVolume *volume = obj->add_volume(make_cube(1., 1., 1.), type);

    // TODO: Refactor to create better way to not set cube at begining
    // Revert mesh centering by set mesh after add cube
    volume->set_mesh(std::move(mesh));
    volume->calculate_convex_hull();

    // set a default extruder value, since user can't add it manually
    volume->config.set_key_value("extruder", new ConfigOptionInt(0));

    // do not allow model reload from disk
    volume->source.is_from_builtin_objects = true;

    volume->name = data.volume_name; // copy

    if (trmat.has_value()) {
        volume->set_transformation(*trmat);    
    } else {
        assert(!data.shape.use_surface);
        // Create transformation for volume near from object(defined by glVolume)
        // Transformation is inspired add generic volumes in ObjectList::load_generic_subobject
        Vec3d volume_size = volume->mesh().bounding_box().size();
        // Translate the new modifier to be pickable: move to the left front corner of the instance's bounding box, lift to print bed.
        Vec3d offset_tr(0, // center of instance - Can't suggest width of text before it will be created
                        -instance_bb.size().y() / 2 - volume_size.y() / 2, // under
                        volume_size.z() / 2 - instance_bb.size().z() / 2); // lay on bed        
        // use same instance as for calculation of instance_bounding_box
        Transform3d tr = obj->instances.front()->get_transformation().get_matrix_no_offset().inverse();
        Transform3d volume_trmat = tr * Eigen::Translation3d(offset_tr);
        volume->set_transformation(volume_trmat);
    }

    data.write(*volume);

    // update printable state on canvas
    if (type == ModelVolumeType::MODEL_PART) {
        volume->get_object()->ensure_on_bed();
        canvas->update_instance_printable_state_for_object(object_idx);
    }

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) { return vol == volume; };
    wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(object_idx, add_to_selection);
    if (!sel.IsEmpty()) obj_list->select_item(sel.front());

    obj_list->selection_changed();

    // Now is valid text volume selected open emboss gizmo
    GLGizmosManager &manager = canvas->get_gizmos_manager();
    if (manager.get_current_type() != GLGizmosManager::Emboss) 
        manager.open_gizmo(GLGizmosManager::Emboss);

    // redraw scene
    canvas->reload_scene(true);
}

ModelVolume *priv::get_volume(ModelObjectPtrs &objects,
                               const ObjectID  &volume_id)
{
    for (ModelObject *obj : objects)
        for (ModelVolume *vol : obj->volumes)
            if (vol->id() == volume_id) return vol;
    return nullptr;
};

OrthoProject priv::create_projection_for_cut(
    Transform3d                    tr,
    double                         shape_scale,
    const std::pair<float, float> &z_range)
{
    double min_z = z_range.first - priv::safe_extension;
    double max_z = z_range.second + priv::safe_extension;
    assert(min_z < max_z);
    // range between min and max value
    double projection_size = max_z - min_z; 
    Matrix3d transformation_for_vector = tr.linear();
    // Projection must be negative value.
    // System of text coordinate
    // X .. from left to right
    // Y .. from bottom to top
    // Z .. from text to eye
    Vec3d untransformed_direction(0., 0., projection_size);
    Vec3d project_direction = transformation_for_vector * untransformed_direction;

    // Projection is in direction from far plane
    tr.translate(Vec3d(0., 0., min_z));
    tr.scale(shape_scale);
    return OrthoProject(tr, project_direction);
}

OrthoProject3d priv::create_emboss_projection(
    bool is_outside, float emboss, Transform3d tr, SurfaceCut &cut)
{
    // Offset of clossed side to model
    const float surface_offset = 0.015f; // [in mm]
    float 
        front_move = (is_outside) ? emboss : surface_offset,
        back_move  = -((is_outside) ? surface_offset : emboss);    
    its_transform(cut, tr.pretranslate(Vec3d(0., 0., front_move)));    
    Vec3d from_front_to_back(0., 0., back_move - front_move);
    return OrthoProject3d(from_front_to_back);
}

// input can't be const - cache of font
TriangleMesh priv::cut_surface(DataBase& base, const SurfaceVolumeData& input2, const std::function<bool()>& was_canceled)
{
    EmbossShape& emboss_shape = base.create_shape();
    ExPolygons& shapes = emboss_shape.shapes; 
    if (shapes.empty())
        throw JobException(_u8L("Font doesn't have any shape for given text.").c_str());

    if (was_canceled()) return {};

    // Define alignment of text - left, right, center, top bottom, ....
    BoundingBox bb  = get_extents(shapes);
    Point       projection_center = bb.center();
    for (ExPolygon &shape : shapes) shape.translate(-projection_center);
    bb.translate(-projection_center);

    const SurfaceVolumeData::ModelSources &sources = input2.sources;
    const SurfaceVolumeData::ModelSource  *biggest = &sources.front();

    size_t biggest_count = 0;
    // convert index from (s)ources to (i)ndexed (t)riangle (s)ets
    std::vector<size_t> s_to_itss(sources.size(), std::numeric_limits<size_t>::max());
    std::vector<indexed_triangle_set>  itss;
    itss.reserve(sources.size());
    for (const SurfaceVolumeData::ModelSource &s : sources) {
        Transform3d mesh_tr_inv       = s.tr.inverse();
        Transform3d cut_projection_tr = mesh_tr_inv * input2.transform;
        std::pair<float, float> z_range{0., 1.};
        OrthoProject cut_projection = create_projection_for_cut(cut_projection_tr, emboss_shape.scale, z_range);
        // copy only part of source model
        indexed_triangle_set its = its_cut_AoI(s.mesh->its, bb, cut_projection);
        if (its.indices.empty()) continue;
        if (biggest_count < its.vertices.size()) {
            biggest_count = its.vertices.size();
            biggest       = &s;
        }
        size_t source_index = &s - &sources.front();
        size_t its_index = itss.size();
        s_to_itss[source_index] = its_index;
        itss.emplace_back(std::move(its));
    }
    if (itss.empty())
        throw JobException(_u8L("There is no volume in projection direction.").c_str());

    Transform3d tr_inv = biggest->tr.inverse();
    Transform3d cut_projection_tr = tr_inv * input2.transform;

    size_t        itss_index = s_to_itss[biggest - &sources.front()];
    BoundingBoxf3 mesh_bb    = bounding_box(itss[itss_index]);
    for (const SurfaceVolumeData::ModelSource &s : sources) {
        size_t itss_index = s_to_itss[&s - &sources.front()];
        if (itss_index == std::numeric_limits<size_t>::max()) continue;
        if (&s == biggest) 
            continue;

        Transform3d tr = s.tr * tr_inv;
        bool fix_reflected = true;
        indexed_triangle_set &its = itss[itss_index];
        its_transform(its, tr, fix_reflected);
        BoundingBoxf3 its_bb = bounding_box(its);
        mesh_bb.merge(its_bb);
    }

    // tr_inv = transformation of mesh inverted
    Transform3d   emboss_tr  = cut_projection_tr.inverse();
    BoundingBoxf3 mesh_bb_tr = mesh_bb.transformed(emboss_tr);
    std::pair<float, float> z_range{mesh_bb_tr.min.z(), mesh_bb_tr.max.z()};
    OrthoProject cut_projection = create_projection_for_cut(cut_projection_tr, emboss_shape.scale, z_range);
    float projection_ratio = (-z_range.first + safe_extension) / (z_range.second - z_range.first + 2 * safe_extension);

    bool is_text_reflected = Slic3r::has_reflection(input2.transform);
    if (is_text_reflected) {
        // revert order of points in expolygons
        // CW --> CCW
        for (ExPolygon &shape : shapes) {
            shape.contour.reverse();
            for (Slic3r::Polygon &hole : shape.holes)
                hole.reverse();
        }
    }

    // Use CGAL to cut surface from triangle mesh
    SurfaceCut cut = cut_surface(shapes, itss, cut_projection, projection_ratio);

    if (is_text_reflected) {
        for (SurfaceCut::Contour &c : cut.contours)
            std::reverse(c.begin(), c.end());
        for (Vec3i &t : cut.indices)
            std::swap(t[0], t[1]);
    }

    if (cut.empty()) throw JobException(_u8L("There is no valid surface for text projection.").c_str());
    if (was_canceled()) return {};

    // !! Projection needs to transform cut    
    OrthoProject3d projection = create_emboss_projection(input2.is_outside, emboss_shape.depth, emboss_tr, cut);
    indexed_triangle_set new_its = cut2model(cut, projection);
    assert(!new_its.empty());
    if (was_canceled()) return {};
    return TriangleMesh(std::move(new_its));
}

SurfaceVolumeData::ModelSources priv::create_sources(const ModelVolumePtrs &volumes, std::optional<size_t> text_volume_id)
{
    SurfaceVolumeData::ModelSources result;
    result.reserve(volumes.size() - 1);
    for (const ModelVolume *v : volumes) {
        if (text_volume_id.has_value() && v->id().id == *text_volume_id)
            continue;
        // skip modifiers and negative volumes, ...
        if (!v->is_model_part())
            continue;
        const TriangleMesh &tm = v->mesh();
        if (tm.empty())
            continue;
        if (tm.its.empty())
            continue;
        result.push_back({v->get_mesh_shared_ptr(), v->get_matrix()});
    }
    return result;
}

bool priv::process(std::exception_ptr &eptr) { 
    if (!eptr) return false;
    try {
        std::rethrow_exception(eptr);
    } catch (priv::JobException &e) {
        create_message(e.what());
        eptr = nullptr;
    }
    return true;
}

bool priv::finalize(bool canceled, std::exception_ptr &eptr, const DataBase &input)
{
    // doesn't care about exception when process was canceled by user
    if (canceled || input.cancel->load()) {
        eptr = nullptr;
        return false;
    }
    return !process(eptr);
}

bool priv::is_valid(ModelVolumeType volume_type)
{
    if (volume_type == ModelVolumeType::MODEL_PART || 
        volume_type == ModelVolumeType::NEGATIVE_VOLUME ||
        volume_type == ModelVolumeType::PARAMETER_MODIFIER )
        return true;

    BOOST_LOG_TRIVIAL(error) << "Can't create embossed volume with this type: " << (int) volume_type;
    return false;
}

bool priv::start_create_volume_job(
    Worker &worker, const ModelObject &object, const std::optional<Transform3d>& volume_tr, DataBasePtr data, ModelVolumeType volume_type, GLGizmosManager::EType gizmo)
{
    bool &use_surface = data->shape.use_surface;
    std::unique_ptr<GUI::Job> job;
    if (use_surface) {
        // Model to cut surface from.
        SurfaceVolumeData::ModelSources sources = create_sources(object.volumes);
        if (sources.empty() || !volume_tr.has_value()) {
            use_surface = false;
        } else {
            bool is_outside = volume_type == ModelVolumeType::MODEL_PART;
            // check that there is not unexpected volume type
            assert(is_outside || volume_type == ModelVolumeType::NEGATIVE_VOLUME || volume_type == ModelVolumeType::PARAMETER_MODIFIER);
            SurfaceVolumeData       sfvd{*volume_tr, is_outside, std::move(sources)};
            CreateSurfaceVolumeData surface_data{std::move(sfvd), std::move(data), volume_type, object.id(), gizmo};
            job = std::make_unique<CreateSurfaceVolumeJob>(std::move(surface_data));
        }
    }
    if (!use_surface) {
        // create volume
        DataCreateVolume create_volume_data{std::move(data), volume_type, object.id(), volume_tr, gizmo};
        job = std::make_unique<CreateVolumeJob>(std::move(create_volume_data));
    }
    return queue_job(worker, std::move(job));
}

const GLVolume * priv::find_closest(
    const Selection &selection, const Vec2d &screen_center, const Camera &camera, const ModelObjectPtrs &objects, Vec2d *closest_center)
{
    assert(closest_center != nullptr);
    const GLVolume               *closest = nullptr;
    const Selection::IndicesList &indices = selection.get_volume_idxs();
    assert(!indices.empty()); // no selected volume
    if (indices.empty())
        return closest;

    double center_sq_distance = std::numeric_limits<double>::max();
    for (unsigned int id : indices) {
        const GLVolume *gl_volume = selection.get_volume(id);
        const ModelVolume *volume = get_model_volume(*gl_volume, objects); 
        if (volume == nullptr || !volume->is_model_part())
            continue;
        Slic3r::Polygon hull        = CameraUtils::create_hull2d(camera, *gl_volume);
        Vec2d           c           = hull.centroid().cast<double>();
        Vec2d           d           = c - screen_center;
        bool            is_bigger_x = std::fabs(d.x()) > std::fabs(d.y());
        if ((is_bigger_x && d.x() * d.x() > center_sq_distance) || (!is_bigger_x && d.y() * d.y() > center_sq_distance))
            continue;

        double distance = d.squaredNorm();
        if (center_sq_distance < distance)
            continue;
        center_sq_distance = distance;

        *closest_center = c;
        closest         = gl_volume;
    }
    return closest;
}

bool priv::start_create_object_job(Plater &plater, DataBasePtr emboss_data, const Vec2d &coor, GLGizmosManager::EType gizmo)
{
    const Camera  &camera    = plater.get_camera();
    const Pointfs &bed_shape = plater.build_volume().bed_shape();

    DataCreateObject data{std::move(emboss_data), coor, camera, bed_shape, gizmo};
    auto job = std::make_unique<CreateObjectJob>(std::move(data));

    Worker &worker = plater.get_ui_job_worker();
    return queue_job(worker, std::move(job));
}

DataBasePtr priv::start_create_volume_on_surface_job(Plater                     &plater,
                                                     DataBasePtr                 data,
                                                     ModelVolumeType             volume_type,
                                                     const Vec2d                &screen_coor,
                                                     const GLVolume             &gl_volume,
                                                     RaycastManager             &raycaster,
                                                     GLGizmosManager::EType      gizmo,
                                                     const std::optional<float> &distance,
                                                     const std::optional<float> &angle,
                                                     bool                       &success)
{
    success = false;
    const ModelObjectPtrs &objects  = plater.model().objects;
    const ModelVolume* volume = get_model_volume(gl_volume, objects);
    const ModelInstance *instance = get_model_instance(gl_volume, objects);
    if (volume == nullptr || instance == nullptr ||
        volume->get_object() == nullptr) {
        // weird situation
        assert(false);
        return data;
    }
    
    const Camera &camera = plater.get_camera();
    GLCanvas3D &canvas = *plater.get_current_canvas3D();

    auto cond = RaycastManager::AllowVolumes({volume->id().id});
    RaycastManager::Meshes meshes = create_meshes(canvas, cond);
    raycaster.actualize(*instance, &cond, &meshes);
    std::optional<RaycastManager::Hit> hit = ray_from_camera(raycaster, screen_coor, camera, &cond);

    // context menu for add text could be open only by right click on an
    // object. After right click, object is selected and object_idx is set
    // also hit must exist. But there is options to add text by object list
    if (!hit.has_value())
        // When model is broken. It could appear that hit miss the object.
        // So add part near by in simmilar manner as right panel do
        return data;

    // Create result volume transformation
    Transform3d surface_trmat = create_transformation_onto_surface(hit->position, hit->normal, Slic3r::GUI::up_limit);

    apply_transformation(angle, distance, surface_trmat);

    Transform3d transform = instance->get_matrix().inverse() * surface_trmat;

    // Try to cast ray into scene and find object for add volume
    Worker &worker = plater.get_ui_job_worker();

    success = priv::start_create_volume_job(worker, *volume->get_object(), transform, std::move(data), volume_type, gizmo);
    return nullptr;
}

#include <wx/msgdlg.h>

void priv::create_message(const std::string &message) {
    wxMessageBox(wxString(message), _L("Issue during embossing the text."),
                 wxOK | wxICON_WARNING);
}
