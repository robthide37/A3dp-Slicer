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
struct EmbossDataCreate;

class EmbossUpdateJob : public Job
{
    std::unique_ptr<EmbossDataUpdate> m_input;
    TriangleMesh                m_result;
public:
    EmbossUpdateJob(std::unique_ptr<EmbossDataUpdate> input) : m_input(std::move(input)) {}
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};

class EmbossCreateJob : public Job
{
    std::unique_ptr<EmbossDataCreate> m_input;
    TriangleMesh m_result;
    Transform3d m_transformation; 
public:
    EmbossCreateJob(std::unique_ptr<EmbossDataCreate> input): m_input(std::move(input)){}
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;

    // <summary>
    /// Create mesh from text
    /// </summary>
    /// <param name="text">Text to convert on mesh</param>
    /// <param name="font">Define shape of characters.
    /// NOTE: Can't be const cache glyphs</param>
    /// <param name="font_prop">Property of font</param>
    /// <param name="ctl">Control for job, check of cancelation</param>
    /// <returns>Triangle mesh model</returns>
    static TriangleMesh create_mesh(const char *      text,
                                    Emboss::FontFile &font,
                                    const FontProp &  font_prop,
                                    Ctl &             ctl);

private:
    static TriangleMesh create_default_mesh();
};


/// <summary>
/// Base data holder for embossing
/// </summary>
struct EmbossDataBase
{
    // Pointer on Data of font (glyph shapes)
    std::shared_ptr<Emboss::FontFile> font_file;
    // font item is not used for create object
    TextConfiguration text_configuration;
    // new volume name created from text
    std::string volume_name;
    EmbossDataBase(std::shared_ptr<Emboss::FontFile> font_file,
                   TextConfiguration                 text_configuration,
                   std::string                       volume_name)
        : font_file(std::move(font_file))
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
    EmbossDataUpdate(std::shared_ptr<Emboss::FontFile> font_file,
                     TextConfiguration                 text_configuration,
                     std::string                       volume_name,
                     ModelVolume *                     volume)
        : EmbossDataBase(std::move(font_file), text_configuration, volume_name)
        , volume(volume)
    {}
};

/// <summary>
/// Hold neccessary data to create embossed text object in job
/// </summary>
struct EmbossDataCreate: public EmbossDataBase
{
    // define embossed volume type
    ModelVolumeType volume_type;

    // define position on screen where to create object
    Vec2d screen_coor;

    // when exist ModelObject where to create volume
    std::optional<int> object_idx;

    // hitted instance transformation
    std::optional<Transform3d> hit_vol_tr;

    // projection property
    Camera camera;

    // shape of bed in case of create volume on bed
    std::vector<Vec2d> bed_shape;

    // used to find point on surface where to create new object
    RaycastManager *raycast_manager;
    // It is inside of GLGizmoEmboss object,
    // so I hope it will survive

    EmbossDataCreate(std::shared_ptr<Emboss::FontFile> font_file,
                     const TextConfiguration &         text_configuration,
                     const std::string &               volume_name,
                     ModelVolumeType                   volume_type,
                     Vec2d                             screen_coor,
                     std::optional<int>                object_idx,
                     const std::optional<Transform3d>& hit_vol_tr,
                     const Camera&                     camera,
                     const std::vector<Vec2d> &        bed_shape,
                     RaycastManager *                  raycast_manager)
        : EmbossDataBase(std::move(font_file), text_configuration, volume_name)
        , volume_type(volume_type)
        , screen_coor(screen_coor)
        , object_idx(object_idx)
        , hit_vol_tr(hit_vol_tr)
        , camera(camera)
        , bed_shape(bed_shape)
        , raycast_manager(raycast_manager)
    {}
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
