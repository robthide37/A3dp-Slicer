#ifndef slic3r_EmbossShape_hpp_
#define slic3r_EmbossShape_hpp_

#include <string>
#include <optional>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/archives/binary.hpp>
#include "Point.hpp" // Transform3d
#include "ExPolygon.hpp"
#include "ExPolygonSerialize.hpp"

namespace Slic3r {

struct EmbossProjection
{
    // Emboss depth, Size in local Z direction
    double depth = 1.; // [in loacal mm] 
    // NOTE: User should see and modify mainly world size not local

    // Flag that result volume use surface cutted from source objects
    bool use_surface = false;

    // enum class Align {
    //     left,
    //     right,
    //     center,
    //     top_left,
    //     top_right,
    //     top_center,
    //     bottom_left,
    //     bottom_right,
    //     bottom_center
    // };
    //// change pivot of volume
    //// When not set, center is used and is not stored
    // std::optional<Align> align;

    // compare TextStyle
    bool operator==(const EmbossProjection &other) const {
        return depth == other.depth && use_surface == other.use_surface;
    }

    // undo / redo stack recovery
    template<class Archive> void serialize(Archive &ar) { ar(depth, use_surface); }
};

/// <summary>
/// Contain plane shape information to be able emboss it and edit it
/// </summary>
struct EmbossShape 
{
    // shape defined by integer point consist only by lines not curves
    ExPolygons shapes;

    // scale of shape, multiplier to get 3d point in mm from integer shape
    double scale = 1.;

    // Define how to emboss shape
    EmbossProjection projection;

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
        ar(shapes, scale, projection, svg_file_path);
        cereal::save(ar, fix_3mf_tr);
    }
    template<class Archive> void load(Archive &ar)
    {
        ar(shapes, scale, projection, svg_file_path);
        cereal::load(ar, fix_3mf_tr);
    }
};

} // namespace Slic3r

#endif // slic3r_EmbossShape_hpp_
