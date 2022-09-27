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
    enum class EMode : unsigned char
    {
        BasicSelection,
        ExtendedSelection
    };

    struct SelectedFeatures
    {
        struct Item
        {
            std::string source;
            std::optional<Measure::SurfaceFeature> feature;

            bool operator == (const Item& other) const {
                if (this->source != other.source) return false;
                return this->feature == other.feature;
            }

            bool operator != (const Item& other) const {
                return !operator == (other);
            }

            void reset() {
                source.clear();
                feature.reset();
            }
        };

        Item first;
        Item second;

        void reset() {
            first.reset();
            second.reset();
        }

        bool operator == (const SelectedFeatures & other) const {
            if (this->first != other.first) return false;
            return this->second == other.second;
        }

        bool operator != (const SelectedFeatures & other) const {
            return !operator == (other);
        }
    };

    EMode m_mode{ EMode::BasicSelection };
    std::unique_ptr<Measure::Measuring> m_measuring;

    PickingModel m_sphere;
    PickingModel m_cylinder;
    PickingModel m_circle;
    PickingModel m_plane;
    struct Dimensioning
    {
        GLModel line;
        GLModel triangle;
        GLModel arc;
    };
    Dimensioning m_dimensioning;

    Transform3d m_volume_matrix{ Transform3d::Identity() };
    std::vector<GLModel> m_plane_models_cache;
    std::map<int, std::shared_ptr<SceneRaycasterItem>> m_raycasters;
    std::optional<Measure::SurfaceFeature> m_curr_feature;
    std::optional<Vec3d> m_curr_point_on_feature_position;
    struct SceneRaycasterState
    {
        std::shared_ptr<SceneRaycasterItem> raycaster{ nullptr };
        bool state{true};

    };
    std::vector<SceneRaycasterState> m_scene_raycasters;

    // These hold information to decide whether recalculation is necessary:
    std::vector<Transform3d> m_volumes_matrices;
    std::vector<ModelVolumeType> m_volumes_types;
    Vec3d m_first_instance_scale{ Vec3d::Ones() };
    Vec3d m_first_instance_mirror{ Vec3d::Ones() };
    float m_last_inv_zoom{ 0.0f };
    std::optional<Measure::SurfaceFeature> m_last_circle;
    int m_last_plane_idx{ -1 };

    bool m_mouse_left_down{ false }; // for detection left_up of this gizmo
    const ModelObject* m_old_model_object{ nullptr };
    const ModelVolume* m_old_model_volume{ nullptr };

    Vec2d m_mouse_pos{ Vec2d::Zero() };

    KeyAutoRepeatFilter m_ctrl_kar_filter;

    SelectedFeatures m_selected_features;

    void update_if_needed();

    void disable_scene_raycasters();
    void restore_scene_raycasters_state();

    void render_dimensioning();

#if ENABLE_MEASURE_GIZMO_DEBUG
    void render_debug_dialog();
#endif // ENABLE_MEASURE_GIZMO_DEBUG

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
