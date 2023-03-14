#ifndef slic3r_GLGizmoSVG_hpp_
#define slic3r_GLGizmoSVG_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/SurfaceDrag.hpp"
#include "slic3r/Utils/RaycastManager.hpp"

#include <optional>
#include <memory>
#include <atomic>

#include "libslic3r/Emboss.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Model.hpp"

#include <imgui/imgui.h>
#include <GL/glew.h>

namespace Slic3r{
    //class AppConfig;
    class GLVolume;
    class ModelVolume;
    enum class ModelVolumeType : int;
}

namespace Slic3r::GUI {
class GLGizmoSVG : public GLGizmoBase
{
public:
    GLGizmoSVG(GLCanvas3D &parent);

    /// <summary>
    /// Create new embossed text volume by type on position of mouse
    /// </summary>
    /// <param name="volume_type">Object part / Negative volume / Modifier</param>
    /// <param name="mouse_pos">Define position of new volume</param>
    /// <returns>True on succesfull start creation otherwise False</returns>
    bool create_volume(ModelVolumeType volume_type, const Vec2d &mouse_pos); // first open file dialog

    /// <summary>
    /// Create new text without given position
    /// </summary>
    /// <param name="volume_type">Object part / Negative volume / Modifier</param>
    /// <returns>True on succesfull start creation otherwise False</returns>
    bool create_volume(ModelVolumeType volume_type); // first open file dialog

    /// <summary>
    /// Check whether volume is object containing only emboss volume
    /// </summary>
    /// <param name="volume">Pointer to volume</param>
    /// <returns>True when object otherwise False</returns>
    static bool is_svg_object(const ModelVolume &volume);

    /// <summary>
    /// Check whether volume has emboss data
    /// </summary>
    /// <param name="volume">Pointer to volume</param>
    /// <returns>True when constain emboss data otherwise False</returns>
    static bool is_svg(const ModelVolume &volume);

protected:
    bool on_init() override;
    std::string on_get_name() const override;
    void on_render() override;
    virtual void on_register_raycasters_for_picking() override;
    virtual void on_unregister_raycasters_for_picking() override;
    void on_render_input_window(float x, float y, float bottom_limit) override;
    bool on_is_activable() const override { return true; }
    bool on_is_selectable() const override { return false; }
    void on_set_state() override;    
    void data_changed() override; // selection changed
    void on_set_hover_id() override{ m_rotate_gizmo.set_hover_id(m_hover_id); }
    void on_enable_grabber(unsigned int id) override { m_rotate_gizmo.enable_grabber(); }
    void on_disable_grabber(unsigned int id) override { m_rotate_gizmo.disable_grabber(); }
    void on_start_dragging() override;
    void on_stop_dragging() override;
    void on_dragging(const UpdateData &data) override;    

    /// <summary>
    /// Rotate by text on dragging rotate grabers
    /// </summary>
    /// <param name="mouse_event">Information about mouse</param>
    /// <returns>Propagete normaly return false.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    bool wants_enter_leave_snapshots() const override { return true; }
    std::string get_gizmo_entering_text() const override { return _u8L("Enter SVG gizmo"); }
    std::string get_gizmo_leaving_text() const override { return _u8L("Leave SVG gizmo"); }
    std::string get_action_snapshot_name() override { return _u8L("SVG actions"); }
private:
    void set_volume_by_selection();
    void reset_volume();

    // create volume from text - main functionality
    bool process();
    void close();
    void draw_window();
    void draw_model_type();

    // process mouse event
    bool on_mouse_for_rotation(const wxMouseEvent &mouse_event);
    bool on_mouse_for_translate(const wxMouseEvent &mouse_event);

    // This configs holds GUI layout size given by translated texts.
    // etc. When language changes, GUI is recreated and this class constructed again,
    // so the change takes effect. (info by GLGizmoFdmSupports.hpp)
    struct GuiCfg
    {
        // Detect invalid config values when change monitor DPI
        double screen_scale;
        float  main_toolbar_height;

        // Zero means it is calculated in init function
        ImVec2 minimal_window_size = ImVec2(0, 0);

        float input_width = 0.f;
        float input_offset = 0.f;

        // Only translations needed for calc GUI size
        struct Translations
        {
            std::string font;
        };
        Translations translations;
    };
    std::optional<const GuiCfg> m_gui_cfg;
    static GuiCfg create_gui_configuration();

    // actual selected only one volume - with emboss data
    ModelVolume *m_volume;

    // When work with undo redo stack there could be situation that 
    // m_volume point to unexisting volume so One need also objectID
    ObjectID m_volume_id;

    // cancel for previous update of volume to cancel finalize part
    std::shared_ptr<std::atomic<bool>> m_job_cancel;

    // Rotation gizmo
    GLGizmoRotate m_rotate_gizmo;
    std::optional<float> m_angle;
    // Value is set only when dragging rotation to calculate actual angle
    std::optional<float> m_rotate_start_angle;

    // TODO: it should be accessible by other gizmo too.
    // May be move to plater?
    RaycastManager m_raycast_manager;
    
    // When true keep up vector otherwise relative rotation
    bool m_keep_up = true;
        
    // setted only when wanted to use - not all the time
    std::optional<ImVec2> m_set_window_offset;

    // Keep data about dragging only during drag&drop
    std::optional<SurfaceDrag> m_surface_drag;

    // For volume on scaled objects
    std::optional<float> m_scale_height;
    std::optional<float> m_scale_depth;
    void calculate_scale();

    // only temporary solution
    static const std::string M_ICON_FILENAME;
};

} // namespace Slic3r::GUI

#endif // slic3r_GLGizmoSVG_hpp_
