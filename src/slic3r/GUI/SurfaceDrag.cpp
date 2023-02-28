#include "SurfaceDrag.hpp"

#include "libslic3r/Model.hpp" // ModelVolume
#include "GLCanvas3D.hpp"
#include "slic3r/Utils/RaycastManager.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "libslic3r/Emboss.hpp"

namespace Slic3r::GUI {
    
/// <summary>
/// Calculate offset from mouse position to center of text
/// </summary>
/// <param name="screen_coor">Position on screen[in Px] e.g. mouse position</param>
/// <param name="volume">Selected volume(text)</param>
/// <param name="camera">Actual position and view direction of camera</param>
/// <returns>Offset in screen coordinate</returns>
static Vec2d calc_screen_offset_to_volume_center(const Vec2d &screen_coor, const ModelVolume &volume, const Camera &camera)
{
    const Transform3d &volume_tr = volume.get_matrix();
    assert(volume.text_configuration.has_value());

    auto calc_offset = [&screen_coor, &volume_tr, &camera, &volume](const Transform3d &instrance_tr) -> Vec2d {
        Transform3d to_world = instrance_tr * volume_tr;

        // Use fix of .3mf loaded tranformation when exist
        if (volume.text_configuration->fix_3mf_tr.has_value())
            to_world = to_world * (*volume.text_configuration->fix_3mf_tr);
        // zero point of volume in world coordinate system
        Vec3d volume_center = to_world.translation();
        // screen coordinate of volume center
        Vec2i coor = CameraUtils::project(camera, volume_center);
        return coor.cast<double>() - screen_coor;
    };

    auto object = volume.get_object();
    assert(!object->instances.empty());
    // Speed up for one instance
    if (object->instances.size() == 1)
        return calc_offset(object->instances.front()->get_matrix());

    Vec2d  nearest_offset;
    double nearest_offset_size = std::numeric_limits<double>::max();
    for (const ModelInstance *instance : object->instances) {
        Vec2d  offset      = calc_offset(instance->get_matrix());
        double offset_size = offset.norm();
        if (nearest_offset_size < offset_size)
            continue;
        nearest_offset_size = offset_size;
        nearest_offset      = offset;
    }
    return nearest_offset;
}

 // Calculate scale in world
static std::optional<double> calc_scale(const Matrix3d &from, const Matrix3d &to, const Vec3d &dir)
{
    Vec3d  from_dir      = from * dir;
    Vec3d  to_dir        = to * dir;
    double from_scale_sq = from_dir.squaredNorm();
    double to_scale_sq   = to_dir.squaredNorm();
    if (is_approx(from_scale_sq, to_scale_sq, 1e-3))
        return {}; // no scale
    return sqrt(from_scale_sq / to_scale_sq);
}

bool on_mouse_surface_drag(const wxMouseEvent         &mouse_event,
                           const Camera               &camera,
                           std::optional<SurfaceDrag> &surface_drag,
                           GLCanvas3D                 &canvas,
                           RaycastManager             &raycast_manager)
{
    // Fix when leave window during dragging
    // Fix when click right button
    if (surface_drag.has_value() && !mouse_event.Dragging()) {
        // write transformation from UI into model
        canvas.do_move(L("Surface move"));

        // allow moving with object again
        canvas.enable_moving(true);
        surface_drag.reset();

        // only left up is correct
        // otherwise it is fix state and return false
        return mouse_event.LeftUp();
    }

    if (mouse_event.Moving())
        return false;

    // detect start text dragging
    if (mouse_event.LeftDown()) {
        // selected volume
        GLVolume *gl_volume = get_selected_gl_volume(canvas);
        if (gl_volume == nullptr)
            return false;

        // is selected volume closest hovered?
        const GLVolumePtrs &gl_volumes = canvas.get_volumes().volumes;
        int hovered_idx = canvas.get_first_hover_volume_idx();
        if (hovered_idx < 0 || 
            hovered_idx >= gl_volumes.size() || 
            gl_volumes[hovered_idx] != gl_volume)
            return false;

        const ModelObject   *object   = get_model_object(*gl_volume, canvas.get_model()->objects);
        const ModelInstance *instance = get_model_instance(*gl_volume, *object);
        const ModelVolume   *volume   = get_model_volume(*gl_volume, *object);

        assert(object != nullptr && instance != nullptr && volume != nullptr);
        if (object == nullptr || instance == nullptr || volume == nullptr)
            return false;

        const ModelVolumePtrs &volumes = object->volumes;
        std::vector<size_t>    allowed_volumes_id;
        if (volumes.size() > 1) {
            allowed_volumes_id.reserve(volumes.size() - 1);
            for (auto &v : volumes) {
                // skip actual selected object
                if (v->id() == volume->id())
                    continue;
                // drag only above part not modifiers or negative surface
                if (!v->is_model_part())
                    continue;
                allowed_volumes_id.emplace_back(v->id().id);
            }
        }
        RaycastManager::AllowVolumes condition(std::move(allowed_volumes_id));
        RaycastManager::Meshes meshes = create_meshes(canvas, condition);
        // initialize raycasters
        // INFO: It could slows down for big objects
        // (may be move to thread and do not show drag until it finish)
        raycast_manager.actualize(instance, &condition, &meshes);

        // wxCoord == int --> wx/types.h
        Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
        Vec2d mouse_pos    = mouse_coord.cast<double>();
        Vec2d mouse_offset = calc_screen_offset_to_volume_center(mouse_pos, *volume, camera);

        Transform3d volume_tr = gl_volume->get_volume_transformation().get_matrix();

        if (volume->text_configuration.has_value()) {
            const TextConfiguration &tc = *volume->text_configuration;
            // fix baked transformation from .3mf store process
            if (tc.fix_3mf_tr.has_value())
                volume_tr = volume_tr * tc.fix_3mf_tr->inverse();
        }

        Transform3d instance_tr     = instance->get_matrix();
        Transform3d instance_tr_inv = instance_tr.inverse();
        Transform3d world_tr        = instance_tr * volume_tr;
        surface_drag                = SurfaceDrag{mouse_offset, world_tr, instance_tr_inv, gl_volume, condition};

        // disable moving with object by mouse
        canvas.enable_moving(false);
        return true;
    }

    // Dragging starts out of window
    if (!surface_drag.has_value())
        return false;

    if (mouse_event.Dragging()) {
        // wxCoord == int --> wx/types.h
        Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
        Vec2d mouse_pos      = mouse_coord.cast<double>();
        Vec2d offseted_mouse = mouse_pos + surface_drag->mouse_offset;

        std::optional<RaycastManager::Hit> hit = ray_from_camera(
            raycast_manager, offseted_mouse, camera, &surface_drag->condition);

        surface_drag->exist_hit = hit.has_value();
        if (!hit.has_value()) {
            // cross hair need redraw
            canvas.set_as_dirty();
            return true;
        }

        auto world_linear = surface_drag->world.linear();
        // Calculate offset: transformation to wanted position
        {
            // Reset skew of the text Z axis:
            // Project the old Z axis into a new Z axis, which is perpendicular to the old XY plane.
            Vec3d old_z         = world_linear.col(2);
            Vec3d new_z         = world_linear.col(0).cross(world_linear.col(1));
            world_linear.col(2) = new_z * (old_z.dot(new_z) / new_z.squaredNorm());
        }

        Vec3d       text_z_world     = world_linear.col(2); // world_linear * Vec3d::UnitZ()
        auto        z_rotation       = Eigen::Quaternion<double, Eigen::DontAlign>::FromTwoVectors(text_z_world, hit->normal);
        Transform3d world_new        = z_rotation * surface_drag->world;
        auto        world_new_linear = world_new.linear();

        // Fix direction of up vector
        {
            Vec3d z_world = world_new_linear.col(2);
            z_world.normalize();
            Vec3d wanted_up = Emboss::suggest_up(z_world);

            Vec3d y_world    = world_new_linear.col(1);
            auto  y_rotation = Eigen::Quaternion<double, Eigen::DontAlign>::FromTwoVectors(y_world, wanted_up);

            world_new        = y_rotation * world_new;
            world_new_linear = world_new.linear();
        }
        
        // Edit position from right
        Transform3d volume_new{Eigen::Translation<double, 3>(surface_drag->instance_inv * hit->position)};
        volume_new.linear() = surface_drag->instance_inv.linear() * world_new_linear;

        // Check that transformation matrix is valid transformation
        assert(volume_new.matrix()(0, 0) == volume_new.matrix()(0, 0)); // Check valid transformation not a NAN
        if (volume_new.matrix()(0, 0) != volume_new.matrix()(0, 0))
            return true;

        // Check that scale in world did not changed
        assert(!calc_scale(world_linear, world_new_linear, Vec3d::UnitY()).has_value());
        assert(!calc_scale(world_linear, world_new_linear, Vec3d::UnitZ()).has_value());

        const ModelVolume *volume = get_model_volume(*surface_drag->gl_volume, canvas.get_model()->objects);
        if (volume != nullptr && volume->text_configuration.has_value()) {
            const TextConfiguration &tc = *volume->text_configuration;
            // fix baked transformation from .3mf store process
            if (tc.fix_3mf_tr.has_value())
                volume_new = volume_new * (*tc.fix_3mf_tr);

            // apply move in Z direction and rotation by up vector
            Emboss::apply_transformation(tc.style.prop, volume_new);
        }

        // Update transformation for all instances
        for (GLVolume *vol : canvas.get_volumes().volumes) {
            if (vol->object_idx() != surface_drag->gl_volume->object_idx() || vol->volume_idx() != surface_drag->gl_volume->volume_idx())
                continue;
            vol->set_volume_transformation(volume_new);
        }

        canvas.set_as_dirty();
        return true;
    }
    return false;
}

} // namespace Slic3r::GUI