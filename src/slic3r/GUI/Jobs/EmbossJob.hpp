#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

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

struct EmbossDataUpdate;
struct EmbossDataCreateVolume;
struct EmbossDataCreateObject;

/// <summary>
/// Update text shape in existing text volume
/// </summary>
class EmbossUpdateJob : public Job
{
    std::unique_ptr<EmbossDataUpdate> m_input;
    TriangleMesh                      m_result;

public:
    EmbossUpdateJob(std::unique_ptr<EmbossDataUpdate> input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};


/// <summary>
/// Create new TextVolume on the surface of ModelObject
/// </summary>
class EmbossCreateVolumeJob : public Job
{
    std::unique_ptr<EmbossDataCreateVolume> m_input;
    TriangleMesh                            m_result;
    Transform3d                             m_transformation;

public:
    EmbossCreateVolumeJob(std::unique_ptr<EmbossDataCreateVolume> input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};

/// <summary>
/// Create new TextObject on the platter
/// </summary>
class EmbossCreateObjectJob : public Job
{
    std::unique_ptr<EmbossDataCreateObject> m_input;
    TriangleMesh                            m_result;
    Transform3d                             m_transformation;

public:
    EmbossCreateObjectJob(std::unique_ptr<EmbossDataCreateObject> input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};

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
    EmbossDataBase(Emboss::FontFileWithCache font_file,
                   TextConfiguration         text_configuration,
                   std::string               volume_name)
        : font_file(font_file)
        , text_configuration(text_configuration)
        , volume_name(volume_name)
    {}
};

/// <summary>
/// Hold neccessary data to update embossed text object in job
/// </summary>
struct EmbossDataUpdate : public EmbossDataBase
{
    // unique identifier of volume to change
    // I can't proove of alive pointer
    ModelVolume *volume;

    // unique identifier of volume to change
    // Change of volume change id, last change could disapear
    // ObjectID     volume_id;
    EmbossDataUpdate(Emboss::FontFileWithCache font_file,
                     TextConfiguration         text_configuration,
                     std::string               volume_name,
                     ModelVolume              *volume)
        : EmbossDataBase(font_file, text_configuration, volume_name)
        , volume(volume)
    {}
};

/// <summary>
/// Hold neccessary data to create ModelVolume in job
/// Volume is created on the surface of existing volume in object.
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

    EmbossDataCreateVolume(Emboss::FontFileWithCache font_file,
                           const TextConfiguration  &text_configuration,
                           const std::string        &volume_name,
                           ModelVolumeType           volume_type,
                           Vec2d                     screen_coor,
                           int                       object_idx,
                           const Camera             &camera,
                           const RaycastManager::SurfacePoint &hit,
                           const Transform3d                  &hit_object_tr,
                           const Transform3d                  &hit_instance_tr)
        : EmbossDataBase(font_file, text_configuration, volume_name)
        , volume_type(volume_type)
        , screen_coor(screen_coor)
        , object_idx(object_idx)
        , camera(camera)
        , hit(hit)
        , hit_object_tr(hit_object_tr)
        , hit_instance_tr(hit_instance_tr)
    {}
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

    EmbossDataCreateObject(Emboss::FontFileWithCache font_file,
                     const TextConfiguration  &text_configuration,
                     const std::string        &volume_name,
                     Vec2d                     screen_coor,
                     const Camera             &camera,
                     const std::vector<Vec2d> &bed_shape)
        : EmbossDataBase(font_file, text_configuration, volume_name)
        , screen_coor(screen_coor)
        , camera(camera)
        , bed_shape(bed_shape)
    {}
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
