#ifndef slic3r_GLGizmoMeasure_hpp_
#define slic3r_GLGizmoMeasure_hpp_

#include "GLGizmoBase.hpp"
#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "slic3r/GUI/GLModel.hpp"
#else
#include "slic3r/GUI/3DScene.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL


namespace Slic3r {

enum class ModelVolumeType : int;


namespace GUI {


class GLGizmoMeasure : public GLGizmoBase
{
// This gizmo does not use grabbers. The m_hover_id relates to polygon managed by the class itself.

private:

    int m_currently_shown_plane = 0;

    struct PlaneData {
        std::vector<std::vector<Vec3d>> borders; // should be in fact local in update_planes()
        std::vector<int> facets;
        std::vector<GLModel> vbos;
        Vec3d normal;
        float area;
    };

    // This holds information to decide whether recalculation is necessary:
    std::vector<Transform3d> m_volumes_matrices;
    std::vector<ModelVolumeType> m_volumes_types;
    Vec3d m_first_instance_scale;
    Vec3d m_first_instance_mirror;

    std::vector<PlaneData> m_planes;
    std::vector<size_t> m_face_to_plane;
    bool m_mouse_left_down = false; // for detection left_up of this gizmo
    bool m_planes_valid = false;
    const ModelObject* m_old_model_object = nullptr;
    std::vector<const Transform3d*> instances_matrices;

    void update_planes();
    bool is_plane_update_necessary() const;
    void set_flattening_data(const ModelObject* model_object);

public:
    GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);
        
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
    void on_render_for_picking() override;
    void on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoMeasure_hpp_
