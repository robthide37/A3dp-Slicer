#ifndef slic3r_GLGizmoFlatten_hpp_
#define slic3r_GLGizmoFlatten_hpp_

#include "GLGizmoBase.hpp"
#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/GLModel.hpp"
#if ENABLE_RAYCAST_PICKING
#include "slic3r/GUI/MeshUtils.hpp"
#endif // ENABLE_RAYCAST_PICKING
#else
#include "slic3r/GUI/3DScene.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

namespace Slic3r {

enum class ModelVolumeType : int;


namespace GUI {


class GLGizmoFlatten : public GLGizmoBase
{
// This gizmo does not use grabbers. The m_hover_id relates to polygon managed by the class itself.

private:

    struct PlaneData {
        std::vector<Vec3d> vertices; // should be in fact local in update_planes()
#if ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_RAYCAST_PICKING
        PickingModel vbo;
#else
        GLModel vbo;
#endif // ENABLE_RAYCAST_PICKING
#else
        GLIndexedVertexArray vbo;
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
        Vec3d normal;
        float area;
#if ENABLE_RAYCAST_PICKING
        int picking_id{ -1 };
#endif // ENABLE_RAYCAST_PICKING
    };

    // This holds information to decide whether recalculation is necessary:
    std::vector<Transform3d> m_volumes_matrices;
    std::vector<ModelVolumeType> m_volumes_types;
    Vec3d m_first_instance_scale;
    Vec3d m_first_instance_mirror;

    std::vector<PlaneData> m_planes;
#if ENABLE_RAYCAST_PICKING
    std::vector<std::shared_ptr<SceneRaycasterItem>> m_planes_casters;
#endif // ENABLE_RAYCAST_PICKING
    bool m_mouse_left_down = false; // for detection left_up of this gizmo
    const ModelObject* m_old_model_object = nullptr;
    std::vector<const Transform3d*> instances_matrices;

    void update_planes();
    bool is_plane_update_necessary() const;

public:
    GLGizmoFlatten(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    void set_flattening_data(const ModelObject* model_object);
        
    /// <summary>
    /// Apply rotation on select plane
    /// </summary>
    /// <param name="mouse_event">Keep information about mouse click</param>
    /// <returns>Return True when use the information otherwise False.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    void data_changed() override;
protected:
    bool on_init() override;
    std::string on_get_name() const override;
    bool on_is_activable() const override;
    void on_render() override;
#if ENABLE_RAYCAST_PICKING
    virtual void on_register_raycasters_for_picking() override;
    virtual void on_unregister_raycasters_for_picking() override;
#else
    void on_render_for_picking() override;
#endif // ENABLE_RAYCAST_PICKING
    void on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoFlatten_hpp_
