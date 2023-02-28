#ifndef slic3r_SurfaceDrag_hpp_
#define slic3r_SurfaceDrag_hpp_

#include <optional>
#include "libslic3r/Point.hpp" // Vec2d, Transform3d
#include "slic3r/Utils/RaycastManager.hpp"
#include "wx/event.h" // wxMouseEvent

namespace Slic3r {
class GLVolume;
} // namespace Slic3r

namespace Slic3r::GUI {
class GLCanvas3D;
struct Camera;

// Data for drag&drop over surface with mouse
struct SurfaceDrag
{
    // hold screen coor offset of cursor from object center
    Vec2d mouse_offset;

    // Start dragging text transformations to world
    Transform3d world;

    // Invers transformation of text volume instance
    // Help convert world transformation to instance space
    Transform3d instance_inv;

    // Dragged gl volume
    GLVolume *gl_volume;

    // condition for raycaster
    RaycastManager::AllowVolumes condition;

    // Flag whether coordinate hit some volume
    bool exist_hit = true;
};

bool on_mouse_surface_drag(const wxMouseEvent         &mouse_event,
                           const Camera               &camera,
                           std::optional<SurfaceDrag> &surface_drag,
                           GLCanvas3D                 &canvas,
                           RaycastManager             &raycast_manager);

} // namespace Slic3r::GUI
#endif // slic3r_SurfaceDrag_hpp_