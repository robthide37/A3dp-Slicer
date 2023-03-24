#include "SurfaceDrag.hpp"

#include "libslic3r/Model.hpp" // ModelVolume
#include "GLCanvas3D.hpp"
#include "slic3r/Utils/RaycastManager.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "libslic3r/Emboss.hpp"

namespace{
// Distance of embossed volume from surface to be represented as distance surface
// Maximal distance is also enlarge by size of emboss depth
constexpr Slic3r::MinMax<double> surface_distance_sq{1e-4, 10.}; // [in mm]
}

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
    assert(volume.emboss_shape.has_value());

    auto calc_offset = [&screen_coor, &volume_tr, &camera, &volume](const Transform3d &instrance_tr) -> Vec2d {
        Transform3d to_world = instrance_tr * volume_tr;

        // Use fix of .3mf loaded tranformation when exist        
        if (std::optional<Transform3d> fix = volume.emboss_shape->fix_3mf_tr;
            fix.has_value())
            to_world = to_world * (*fix);
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

 // Calculate scale in world for check in debug
[[maybe_unused]] static std::optional<double> calc_scale(const Matrix3d &from, const Matrix3d &to, const Vec3d &dir)
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
                           RaycastManager             &raycast_manager,
                           std::optional<double>       up_limit)
{
    // Fix when leave window during dragging
    // Fix when click right button
    if (surface_drag.has_value() && !mouse_event.Dragging()) {
        // write transformation from UI into model
        canvas.do_move(L("Surface move"));

        // allow moving with object again
        canvas.enable_moving(true);
        canvas.enable_picking(true);
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
        GLVolume *gl_volume_ptr = get_selected_gl_volume(canvas);
        if (gl_volume_ptr == nullptr)
            return false;
        const GLVolume &gl_volume = *gl_volume_ptr;

        // is selected volume closest hovered?
        const GLVolumePtrs &gl_volumes = canvas.get_volumes().volumes;        
        if (int hovered_idx = canvas.get_first_hover_volume_idx();
            hovered_idx < 0)
            return false;        
        else if (auto hovered_idx_ = static_cast<size_t>(hovered_idx);
            hovered_idx_ >= gl_volumes.size() || 
            gl_volumes[hovered_idx_] != gl_volume_ptr)
            return false;

        const ModelObject *object = get_model_object(gl_volume, canvas.get_model()->objects);
        assert(object != nullptr);
        if (object == nullptr)
            return false;

        const ModelInstance *instance = get_model_instance(gl_volume, *object);
        const ModelVolume   *volume   = get_model_volume(gl_volume, *object);
        assert(instance != nullptr && volume != nullptr);
        if (object == nullptr || instance == nullptr || volume == nullptr)
            return false;

        // allowed drag&drop by canvas for object
        if (volume->is_the_only_one_part())
            return false;

        const ModelVolumePtrs &volumes = object->volumes;
        std::vector<size_t>    allowed_volumes_id;
        if (volumes.size() > 1) {
            allowed_volumes_id.reserve(volumes.size() - 1);
            for (const ModelVolume *v : volumes) {
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
        raycast_manager.actualize(*instance, &condition, &meshes);

        // wxCoord == int --> wx/types.h
        Vec2i mouse_coord(mouse_event.GetX(), mouse_event.GetY());
        Vec2d mouse_pos    = mouse_coord.cast<double>();
        Vec2d mouse_offset = calc_screen_offset_to_volume_center(mouse_pos, *volume, camera);

        Transform3d volume_tr = gl_volume.get_volume_transformation().get_matrix();
        
        // fix baked transformation from .3mf store process
        if (const std::optional<EmbossShape> &es_opt = volume->emboss_shape; 
            es_opt.has_value()) {
            const std::optional<Slic3r::Transform3d> &fix = es_opt->fix_3mf_tr;
            if (fix.has_value())
                volume_tr = volume_tr * fix->inverse();
        }

        Transform3d instance_tr     = instance->get_matrix();
        Transform3d instance_tr_inv = instance_tr.inverse();
        Transform3d world_tr        = instance_tr * volume_tr;
        std::optional<float> start_angle;
        if (up_limit.has_value())
            start_angle = Emboss::calc_up(world_tr, *up_limit);

        std::optional<float> start_distance;        
        if (!volume->emboss_shape->projection.use_surface)
            start_distance = calc_distance(gl_volume, raycast_manager, &condition);
        surface_drag = SurfaceDrag{mouse_offset, world_tr, instance_tr_inv, gl_volume_ptr, condition, start_angle, start_distance};

        // disable moving with object by mouse
        canvas.enable_moving(false);
        canvas.enable_picking(false);
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

        // Fix direction of up vector to zero initial rotation
        if(up_limit.has_value()){
            Vec3d z_world = world_new_linear.col(2);
            z_world.normalize();
            Vec3d wanted_up = Emboss::suggest_up(z_world, *up_limit);

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
        // fix baked transformation from .3mf store process
        if (volume != nullptr && volume->emboss_shape.has_value()) {
            const std::optional<Slic3r::Transform3d> &fix = volume->emboss_shape->fix_3mf_tr;
            if (fix.has_value())
                volume_new = volume_new * (*fix);

            // apply move in Z direction and rotation by up vector
            Emboss::apply_transformation(surface_drag->start_angle, surface_drag->start_distance, volume_new);
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

std::optional<Vec3d> calc_surface_offset(const Selection &selection, RaycastManager &raycast_manager) {
    const GLVolume *gl_volume_ptr = get_selected_gl_volume(selection);
    if (gl_volume_ptr == nullptr)
        return {};
    const GLVolume& gl_volume = *gl_volume_ptr;

    const ModelObjectPtrs &objects = selection.get_model()->objects;
    const ModelVolume* volume = get_model_volume(gl_volume, objects);
    if (volume == nullptr)
        return {};

    const ModelInstance* instance = get_model_instance(gl_volume, objects);
    if (instance == nullptr)
        return {};

    // Move object on surface
    auto cond = RaycastManager::SkipVolume(volume->id().id);
    raycast_manager.actualize(*instance, &cond);

    Transform3d to_world = world_matrix_fixed(gl_volume, selection.get_model()->objects);
    Vec3d point     = to_world * Vec3d::Zero();
    Vec3d direction = to_world.linear() * (-Vec3d::UnitZ());

    // ray in direction of text projection(from volume zero to z-dir)
    std::optional<RaycastManager::Hit> hit_opt = raycast_manager.closest_hit(point, direction, &cond);

    // Try to find closest point when no hit object in emboss direction
    if (!hit_opt.has_value()) {
        std::optional<RaycastManager::ClosePoint> close_point_opt = raycast_manager.closest(point);

        // It should NOT appear. Closest point always exists.
        assert(close_point_opt.has_value());
        if (!close_point_opt.has_value())
            return {};

        // It is no neccesary to move with origin by very small value
        if (close_point_opt->squared_distance < EPSILON)
            return {};

        const RaycastManager::ClosePoint &close_point = *close_point_opt;
        Transform3d hit_tr = raycast_manager.get_transformation(close_point.tr_key);
        Vec3d    hit_world = hit_tr * close_point.point;
        Vec3d offset_world = hit_world - point; // vector in world
        Vec3d offset_volume = to_world.inverse().linear() * offset_world;
        return offset_volume;
    }

    // It is no neccesary to move with origin by very small value
    const RaycastManager::Hit &hit = *hit_opt;
    if (hit.squared_distance < EPSILON)
        return {};
    Transform3d hit_tr = raycast_manager.get_transformation(hit.tr_key);
    Vec3d hit_world    = hit_tr * hit.position;
    Vec3d offset_world = hit_world - point; // vector in world
    // TIP: It should be close to only z move
    Vec3d offset_volume = to_world.inverse().linear() * offset_world;
    return offset_volume;
}

std::optional<float> calc_distance(const GLVolume &gl_volume, RaycastManager &raycaster, GLCanvas3D &canvas)
{
    const ModelObject *object = get_model_object(gl_volume, canvas.get_model()->objects);
    assert(object != nullptr);
    if (object == nullptr)
        return {};

    const ModelInstance *instance = get_model_instance(gl_volume, *object);
    const ModelVolume   *volume   = get_model_volume(gl_volume, *object);
    assert(instance != nullptr && volume != nullptr);
    if (object == nullptr || instance == nullptr || volume == nullptr)
        return {};

    if (volume->is_the_only_one_part())
        return {};

    const ModelVolumePtrs &volumes = object->volumes;
    std::vector<size_t> allowed_volumes_id;
    allowed_volumes_id.reserve(volumes.size() - 1);
    for (const ModelVolume *v : volumes) {
        // skip actual selected object
        if (v->id() == volume->id())
            continue;
        // collect hit only from object parts not modifiers neither negative
        if (!v->is_model_part())
            continue;
        allowed_volumes_id.emplace_back(v->id().id);
    }
    RaycastManager::AllowVolumes condition(std::move(allowed_volumes_id));
    RaycastManager::Meshes meshes = create_meshes(canvas, condition);
    raycaster.actualize(*instance, &condition, &meshes);
    return calc_distance(gl_volume, raycaster, &condition);
}

std::optional<float> calc_distance(const GLVolume &gl_volume, const RaycastManager &raycaster, const RaycastManager::ISkip *condition)
{
    Transform3d w = gl_volume.world_matrix();
    Vec3d p = w.translation();
    const Vec3d& dir = get_z_base(w);
    auto hit_opt = raycaster.closest_hit(p, dir, condition);
    if (!hit_opt.has_value())
        return {};
    const RaycastManager::Hit &hit = *hit_opt;

    // too small distance is calculated as zero distance
    if (hit.squared_distance < ::surface_distance_sq.min)
        return {};

    // check maximal distance
    const BoundingBoxf3& bb = gl_volume.bounding_box();
    double max_squared_distance = std::max(std::pow(2 * bb.size().z(), 2), ::surface_distance_sq.max);
    if (hit.squared_distance > max_squared_distance)
        return {};

    // calculate sign
    float sign = ((hit.position - p).dot(dir) > 0)? 1.f : -1.f;

    // distiguish sign
    return sign * static_cast<float>(sqrt(hit.squared_distance));
}

Transform3d world_matrix_fixed(const GLVolume &gl_volume, const ModelObjectPtrs &objects)
{
    Transform3d res = gl_volume.world_matrix();

    const ModelVolume *mv = get_model_volume(gl_volume, objects);
    if (!mv)
        return res;

    const std::optional<EmbossShape> &es = mv->emboss_shape;
    if (!es.has_value())
        return res;

    const std::optional<Transform3d> &fix = es->fix_3mf_tr;
    if (!fix.has_value())
        return res;

    return res * fix->inverse();
}

Transform3d world_matrix_fixed(const Selection &selection)
{
    const GLVolume *gl_volume = get_selected_gl_volume(selection);
    assert(gl_volume != nullptr);
    if (gl_volume == nullptr)
        return Transform3d::Identity();

    return world_matrix_fixed(*gl_volume, selection.get_model()->objects);
}

bool face_selected_volume_to_camera(const Camera &camera, GLCanvas3D &canvas)
{
    const Vec3d &cam_dir = camera.get_dir_forward();
    Selection   &sel     = canvas.get_selection();
    if (sel.is_empty())
        return false;

    // camera direction transformed into volume coordinate system
    Transform3d to_world   = world_matrix_fixed(sel);
    Vec3d       cam_dir_tr = to_world.inverse().linear() * cam_dir;
    cam_dir_tr.normalize();

    Vec3d emboss_dir(0., 0., -1.);

    // check wether cam_dir is already used
    if (is_approx(cam_dir_tr, emboss_dir))
        return false;

    assert(sel.get_volume_idxs().size() == 1);
    GLVolume *gl_volume = sel.get_volume(*sel.get_volume_idxs().begin());

    Transform3d vol_rot;
    Transform3d vol_tr = gl_volume->get_volume_transformation().get_matrix();
    // check whether cam_dir is opposit to emboss dir
    if (is_approx(cam_dir_tr, -emboss_dir)) {
        // rotate 180 DEG by y
        vol_rot = Eigen::AngleAxis(M_PI_2, Vec3d(0., 1., 0.));
    } else {
        // calc params for rotation
        Vec3d axe = emboss_dir.cross(cam_dir_tr);
        axe.normalize();
        double angle = std::acos(emboss_dir.dot(cam_dir_tr));
        vol_rot      = Eigen::AngleAxis(angle, axe);
    }

    Vec3d       offset     = vol_tr * Vec3d::Zero();
    Vec3d       offset_inv = vol_rot.inverse() * offset;
    Transform3d res        = vol_tr * Eigen::Translation<double, 3>(-offset) * vol_rot * Eigen::Translation<double, 3>(offset_inv);
    // Transform3d res = vol_tr * vol_rot;
    gl_volume->set_volume_transformation(Geometry::Transformation(res));
    get_model_volume(*gl_volume, sel.get_model()->objects)->set_transformation(res);
    return true;
}

void do_local_z_rotate(GLCanvas3D &canvas, double relative_angle)
{
    Selection &selection = canvas.get_selection();

    assert(!selection.is_empty());
    if(selection.is_empty()) return;

    selection.setup_cache();
    TransformationType transformation_type = TransformationType::Local_Relative_Joint;
    selection.rotate(Vec3d(0., 0., relative_angle), transformation_type);

    std::string snapshot_name; // empty meand no store undo / redo
    // NOTE: it use L instead of _L macro because prefix _ is appended
    // inside function do_move
    // snapshot_name = L("Set text rotation");
    canvas.do_rotate(snapshot_name);
}

void do_local_z_move(GLCanvas3D &canvas, double relative_move) {
    
    Selection &selection = canvas.get_selection();
    assert(!selection.is_empty());
    if (selection.is_empty()) return;

    selection.setup_cache();
    Vec3d translate = Vec3d::UnitZ() * relative_move;
    selection.translate(translate, TransformationType::Local);

    std::string snapshot_name; // empty mean no store undo / redo
    // NOTE: it use L instead of _L macro because prefix _ is appended inside
    // function do_move
    // snapshot_name = L("Set surface distance");
    canvas.do_move(snapshot_name);
}

} // namespace Slic3r::GUI