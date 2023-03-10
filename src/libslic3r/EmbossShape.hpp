#ifndef slic3r_EmbossShape_hpp_
#define slic3r_EmbossShape_hpp_

#include <string>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/archives/binary.hpp>
#include "Point.hpp" // Transform3d
#include "ExPolygon.hpp"
#include "ExPolygonSerialize.hpp"

namespace Slic3r {

/// <summary>
/// Contain plane shape information to be able emboss it and edit it
/// </summary>
struct EmbossShape 
{
    // shape defined by integer point consist only by lines not curves
    ExPolygons shapes;

    // scale of shape, multiplier to get 3d point in mm from integer shape
    double scale;
        
    // Emboss depth, Size in local Z direction
    double depth; // [in loacal mm] 
    // NOTE: User should see and modify mainly world size not local

    // Flag that result volume use surface cutted from source objects
    bool use_surface = false;

    // distance from surface point
    // used for move over model surface
    // When not set value is zero and is not stored
    // NOTE: Can't be used together with use_surface
    std::optional<float> distance; // [in mm]

    // !!! Volume stored in .3mf has transformed vertices.
    // (baked transformation into vertices position)
    // Only place for fill this is when load from .3mf
    // This is correction for volume transformation
    // Stored_Transform3d * fix_3mf_tr = Transform3d_before_store_to_3mf
    std::optional<Slic3r::Transform3d> fix_3mf_tr;

    // file(.svg) path to source of shape
    // When empty can't reload from disk
    std::string svg_file_path;
    
    // undo / redo stack recovery
    template<class Archive> void save(Archive &ar) const
    {
        ar(shapes, scale, depth, use_surface, svg_file_path);
        cereal::save(ar, distance);
        cereal::save(ar, fix_3mf_tr);
    }
    template<class Archive> void load(Archive &ar)
    {
        ar(shapes, scale, depth, use_surface, svg_file_path);
        cereal::load(ar, distance);
        cereal::load(ar, fix_3mf_tr);
    }
};

} // namespace Slic3r

#endif // slic3r_EmbossShape_hpp_
