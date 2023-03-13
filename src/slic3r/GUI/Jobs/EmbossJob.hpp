#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include <atomic>
#include <memory>
#include <string>
#include "libslic3r/Emboss.hpp"
#include "libslic3r/EmbossShape.hpp"
#include "libslic3r/Point.hpp" // Transform3d

#include "slic3r/Utils/RaycastManager.hpp"

#include "slic3r/GUI/Jobs/EmbossJob.hpp" // Emboss::DataBase
#include "slic3r/GUI/Camera.hpp"

#include "Job.hpp"

// forward declarations
namespace Slic3r {
class GLVolume;
class ModelVolume;
class ModelObject;
class TriangleMesh;
typedef std::vector<ModelObject *> ModelObjectPtrs;
typedef std::vector<ModelVolume *> ModelVolumePtrs;
namespace GUI {
class Selection;
class RaycastManager;
class Worker;
}}

namespace Slic3r::GUI::Emboss {

/// <summary>
/// Base data hold data for create emboss shape
/// </summary>
class DataBase
{
public:
    DataBase(std::string volume_name, std::shared_ptr<std::atomic<bool>> cancel) : volume_name(volume_name), cancel(std::move(cancel)) {}
    DataBase(std::string volume_name, std::shared_ptr<std::atomic<bool>> cancel, EmbossShape shape)
        : volume_name(volume_name), cancel(std::move(cancel)), shape(std::move(shape))
    {}
    virtual ~DataBase() {}

    /// <summary>
    /// Create shape
    /// e.g. Text extract glyphs from font
    /// Not 'const' function because it could modify shape
    /// </summary>
    virtual EmbossShape &create_shape() { return shape; };

    /// <summary>
    /// Write data how to reconstruct shape to volume
    /// </summary>
    /// <param name="volume">Data object for store emboss params</param>
    virtual void write(ModelVolume &volume) const
    {
        volume.name         = volume_name;
        volume.emboss_shape = shape;
    };
        
    // new volume name
    std::string volume_name;

    // flag that job is canceled
    // for time after process.
    std::shared_ptr<std::atomic<bool>> cancel;

    // shape to emboss
    EmbossShape shape;
};

/// <summary>
/// Hold neccessary data to create ModelVolume in job
/// Volume is created on the surface of existing volume in object.
/// NOTE: EmbossDataBase::font_file doesn't have to be valid !!!
/// </summary>
struct DataCreateVolume
{
    // Hold data about shape
    std::unique_ptr<DataBase> base;

    // define embossed volume type
    ModelVolumeType volume_type;

    // parent ModelObject index where to create volume
    ObjectID object_id;

    // new created volume transformation
    Transform3d trmat;
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
    CreateVolumeJob(DataCreateVolume&& input);
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
    std::unique_ptr<DataBase> base;

    // define position on screen where to create object
    Vec2d screen_coor;

    // projection property
    Camera camera;

    // shape of bed in case of create volume on bed
    std::vector<Vec2d> bed_shape;
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
    CreateObjectJob(DataCreateObject&& input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to update embossed text object in job
/// </summary>
struct DataUpdate
{
    // Hold data about shape
    std::unique_ptr<DataBase> base;

    // unique identifier of volume to change
    ObjectID volume_id;
};

/// <summary>
/// Update text shape in existing text volume
/// Predict that there is only one runnig(not canceled) instance of it
/// </summary>
class UpdateJob : public Job
{
    DataUpdate   m_input;
    TriangleMesh m_result;

public:
    // move params to private variable
    UpdateJob(DataUpdate &&input);

    /// <summary>
    /// Create new embossed volume by m_input data and store to m_result
    /// </summary>
    /// <param name="ctl">Control containing cancel flag</param>
    void process(Ctl &ctl) override;

    /// <summary>
    /// Update volume - change object_id
    /// </summary>
    /// <param name="canceled">Was process canceled.
    /// NOTE: Be carefull it doesn't care about
    /// time between finished process and started finalize part.</param>
    /// <param name="">unused</param>
    void finalize(bool canceled, std::exception_ptr &eptr) override;

    /// <summary>
    /// Update text volume
    /// </summary>
    /// <param name="volume">Volume to be updated</param>
    /// <param name="mesh">New Triangle mesh for volume</param>
    /// <param name="base">Data to write into volume</param>
    static void update_volume(ModelVolume *volume, TriangleMesh &&mesh, const DataBase &base);
};

struct SurfaceVolumeData
{
    // Transformation of volume inside of object
    Transform3d transform;

    // Define projection move
    // True (raised) .. move outside from surface
    // False (engraved).. move into object
    bool is_outside;

    struct ModelSource
    {
        // source volumes
        std::shared_ptr<const TriangleMesh> mesh;
        // Transformation of volume inside of object
        Transform3d tr;
    };
    using ModelSources = std::vector<ModelSource>;
    ModelSources sources;
};

/// <summary>
/// Hold neccessary data to create(cut) volume from surface object in job
/// </summary>
struct CreateSurfaceVolumeData : public SurfaceVolumeData{     
    // Hold data about shape
    std::unique_ptr<DataBase> base;

    // define embossed volume type
    ModelVolumeType volume_type;

    // parent ModelObject index where to create volume
    ObjectID object_id;
};

/// <summary>
/// Cut surface from object and create cutted volume
/// Should not be stopped
/// </summary>
class CreateSurfaceVolumeJob : public Job
{
    CreateSurfaceVolumeData m_input;
    TriangleMesh           m_result;

public:
    CreateSurfaceVolumeJob(CreateSurfaceVolumeData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to update embossed text object in job
/// </summary>
struct UpdateSurfaceVolumeData : public DataUpdate, public SurfaceVolumeData{};

/// <summary>
/// Update text volume to use surface from object
/// </summary>
class UpdateSurfaceVolumeJob : public Job
{
    UpdateSurfaceVolumeData m_input;
    TriangleMesh            m_result;

public:
    // move params to private variable
    UpdateSurfaceVolumeJob(UpdateSurfaceVolumeData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Copied triangles from object to be able create mesh for cut surface from
/// </summary>
/// <param name="volumes">Source object volumes for cut surface from</param>
/// <param name="text_volume_id">Source volume id</param>
/// <returns>Source data for cut surface from</returns>
SurfaceVolumeData::ModelSources create_sources(const ModelVolumePtrs &volumes, std::optional<size_t> text_volume_id = {});

/// <summary>
/// Copied triangles from object to be able create mesh for cut surface from
/// </summary>
/// <param name="text_volume">Define text in object</param>
/// <returns>Source data for cut surface from</returns>
SurfaceVolumeData::ModelSources create_volume_sources(const ModelVolume *text_volume);

using DataBasePtr = std::unique_ptr<DataBase>;

/// <summary>
/// Start job for add new volume to object with given transformation
/// </summary>
/// <param name="worker">Define where to queue the job. e.g. wxGetApp().plater()->get_ui_job_worker()</param>
/// <param name="object">Define where to add</param>
/// <param name="volume_tr">Wanted volume transformation</param>
/// <param name="data">Define what to emboss - shape</param>
/// <param name="volume_type">Type of volume: Part, negative, modifier</param>
/// <return>True on success otherwise false</return>
bool start_create_volume_job(Worker &worker, const ModelObject &object, const Transform3d volume_tr, DataBasePtr data, ModelVolumeType volume_type);

/// <summary>
/// Start job for add new volume on surface of object defined by screen coor
/// </summary>
/// <param name="screen_coor">Mouse position which define position</param>
/// <param name="volume">Volume to find surface for create</param>
/// <param name="instance">Instance to find surface for create</param>
/// <param name="raycaster">Ability to ray cast to model</param>
/// <param name="canvas">Contain already used scene RayCasters</param>
/// <param name="angle">Initial z move</param>
/// <param name="angle">Initial z rotation</param>
/// <returns>Volume transformation otherwise there is no hit surface by screen coor</returns>
std::optional<Transform3d> create_volume_transformation_on_surface(const Vec2d                &screen_coor,
                                                                   const Camera               &camera,
                                                                   const ModelVolume          &volume,
                                                                   const ModelInstance        &instance,
                                                                   RaycastManager             &raycaster,
                                                                   GLCanvas3D                 &canvas,
                                                                   const std::optional<float> &distance = {},
                                                                   const std::optional<float> &angle    = {});

/// <summary>
/// Create transformation for volume near from object(defined by glVolume)
/// </summary>
/// <param name="gl_volume">Define object</param>
/// <param name="objects">All objects</param>
/// <param name="volume_height">Y Size of embossed volume [mm in instance]</param>
/// <param name="volume_depth">Z size of embossed volume - emboss depth[mm in instance]</param>
/// <returns>Transformation for new created volume</returns>
Transform3d create_volume_transformation(const GLVolume& gl_volume, const ModelObjectPtrs &objects, float volume_height, float volume_depth);

/// <summary>
/// Find volume in selected objects with closest convex hull to screen center.
/// </summary>
/// <param name="selection">Define where to search for closest</param>
/// <param name="screen_center">Canvas center(dependent on camera settings)</param>
/// <param name="objects">Actual objects</param>
/// <param name="closest_center">OUT: coordinate of controid of closest volume</param>
/// <returns>closest volume when exists otherwise nullptr</returns>
const GLVolume *find_closest(
    const Selection &selection, const Vec2d &screen_center, const Camera &camera, const ModelObjectPtrs &objects, Vec2d *closest_center);

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
