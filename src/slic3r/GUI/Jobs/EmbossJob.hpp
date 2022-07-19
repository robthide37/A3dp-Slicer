#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include <atomic>
#include <memory>
#include <string>
#include <libslic3r/Emboss.hpp>
#include <libslic3r/ModelVolumeType.hpp>
#include "slic3r/Utils/RaycastManager.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "Job.hpp"

namespace Slic3r {
class ModelVolume;
class TriangleMesh;
}

namespace Slic3r::GUI {

/// <summary>
/// Base data holder for embossing
/// </summary>
struct EmbossDataBase
{
    // Keep pointer on Data of font (glyph shapes)
    Emboss::FontFileWithCache font_file;
    // font item is not used for create object
    TextConfiguration text_configuration;
    // new volume name created from text
    std::string volume_name;
};

/// <summary>
/// Hold neccessary data to create ModelVolume in job
/// Volume is created on the surface of existing volume in object.
/// NOTE: EmbossDataBase::font_file doesn't have to be valid !!!
/// </summary>
struct EmbossDataCreateVolume : public EmbossDataBase
{
    // define embossed volume type
    ModelVolumeType volume_type;

    // define position on screen where to create object
    Vec2d screen_coor;

    // parent ModelObject index where to create volume
    int object_idx;

    // projection property
    Camera camera;

    // used to find point on surface where to create new object
    RaycastManager::SurfacePoint hit;
    Transform3d                  hit_object_tr;
    Transform3d                  hit_instance_tr;
};

/// <summary>
/// Create new TextVolume on the surface of ModelObject
/// Should not be stopped
/// NOTE: EmbossDataBase::font_file doesn't have to be valid !!!
/// </summary>
class EmbossCreateVolumeJob : public Job
{
    EmbossDataCreateVolume m_input;
    TriangleMesh           m_result;
    Transform3d            m_transformation;

public:
    EmbossCreateVolumeJob(EmbossDataCreateVolume&& input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to create ModelObject in job
/// Object is placed on bed under screen coor
/// OR to center of scene when it is out of bed shape
/// </summary>
struct EmbossDataCreateObject : public EmbossDataBase
{
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
class EmbossCreateObjectJob : public Job
{
    EmbossDataCreateObject m_input;
    TriangleMesh           m_result;
    Transform3d            m_transformation;
public:
    EmbossCreateObjectJob(EmbossDataCreateObject&& input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

/// <summary>
/// Hold neccessary data to update embossed text object in job
/// </summary>
struct EmbossDataUpdate : public EmbossDataBase
{
    // unique identifier of volume to change
    ObjectID volume_id;

    // flag that job is canceled
    // for time after process.
    std::shared_ptr<std::atomic<bool>> cancel;
};

/// <summary>
/// Update text shape in existing text volume
/// Predict that there is only one runnig(not canceled) instance of it
/// </summary>
class EmbossUpdateJob : public Job
{
    EmbossDataUpdate m_input;
    TriangleMesh     m_result;

public:
    // move params to private variable
    EmbossUpdateJob(EmbossDataUpdate &&input);

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
};

/// <summary>
/// Hold neccessary data to update embossed text object in job
/// </summary>
struct UseSurfaceData : public EmbossDataUpdate
{
    // Transformation of text volume inside of object
    Transform3d text_tr;

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

    /// <summary>
    /// Copied triangles from object to be able create mesh for cut surface from
    /// </summary>
    /// <param name="text_volume">Define text in object</param>
    /// <returns>Source data for cut surface from</returns>
    static ModelSources create_sources(const ModelVolume *text_volume);
};

/// <summary>
/// Update text volume to use surface from object
/// </summary>
class UseSurfaceJob : public Job
{
    UseSurfaceData m_input;
    TriangleMesh   m_result;

public:
    // move params to private variable
    UseSurfaceJob(UseSurfaceData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &eptr) override;
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
