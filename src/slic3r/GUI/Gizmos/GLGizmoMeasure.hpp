#ifndef slic3r_GLGizmoMeasure_hpp_
#define slic3r_GLGizmoMeasure_hpp_

#if ENABLE_MEASURE_GIZMO

#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "slic3r/GUI/GUI_Utils.hpp"
#include "libslic3r/Measure.hpp"

namespace Slic3r {

class ModelVolume;

enum class ModelVolumeType : int;

namespace Measure { class Measuring; }


namespace GUI {

enum class SLAGizmoEventType : unsigned char;

class GLGizmoMeasure : public GLGizmoBase
{
// This gizmo does not use grabbers. The m_hover_id relates to polygon managed by the class itself.

    enum class EMode : unsigned char
    {
        BasicSelection,
        ExtendedSelection
    };

    EMode m_mode{ EMode::BasicSelection };
    std::unique_ptr<Measure::Measuring> m_measuring;

    PickingModel m_sphere;
    PickingModel m_cylinder;
    PickingModel m_circle;
    PickingModel m_plane;

    std::vector<GLModel> m_plane_models_cache;
    std::map<int, std::shared_ptr<SceneRaycasterItem>> m_raycasters;
    std::optional<Measure::SurfaceFeature> m_curr_feature;

    // This holds information to decide whether recalculation is necessary:
    std::vector<Transform3d> m_volumes_matrices;
    std::vector<ModelVolumeType> m_volumes_types;
    Vec3d m_first_instance_scale{ Vec3d::Ones() };
    Vec3d m_first_instance_mirror{ Vec3d::Ones() };

    bool m_mouse_left_down = false; // for detection left_up of this gizmo
    const ModelObject* m_old_model_object = nullptr;
    const ModelVolume* m_old_model_volume = nullptr;

    Vec2d m_mouse_pos{ Vec2d::Zero() };

    KeyAutoRepeatFilter m_ctrl_kar_filter;

    void update_if_needed();

public:
    GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    /// <summary>
    /// Apply rotation on select plane
    /// </summary>
    /// <param name="mouse_event">Keep information about mouse click</param>
    /// <returns>Return True when use the information otherwise False.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    void data_changed() override;

    bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);

protected:
    bool on_init() override;
    std::string on_get_name() const override;
    bool on_is_activable() const override;
    void on_render() override;
    void on_render_for_picking() override;
    void on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;

    virtual void on_render_input_window(float x, float y, float bottom_limit) override;
    virtual void on_register_raycasters_for_picking() override;
    virtual void on_unregister_raycasters_for_picking() override;
};

} // namespace GUI
} // namespace Slic3r

#endif // ENABLE_MEASURE_GIZMO

#endif // slic3r_GLGizmoMeasure_hpp_
