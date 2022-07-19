#ifndef slic3r_GLGizmoCut_hpp_
#define slic3r_GLGizmoCut_hpp_

#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "GLGizmoMove.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ObjectID.hpp"

namespace Slic3r {

enum class CutConnectorType : int;
class ModelVolume;

namespace GUI {
class Selection;

enum class SLAGizmoEventType : unsigned char;

class GLGizmoCut3D : public GLGizmoBase
{
    Transform3d                 m_rotation_m{ Transform3d::Identity() };
    double                      m_snap_step{ 1.0 };
    int                         m_connectors_group_id;

    // archived values 
    Vec3d m_ar_plane_center { Vec3d::Zero() };
    Vec3d m_ar_rotations    { Vec3d::Zero() };

    Vec3d m_plane_center{ Vec3d::Zero() };
    // data to check position of the cut palne center on gizmo activation
    Vec3d m_min_pos{ Vec3d::Zero() };
    Vec3d m_max_pos{ Vec3d::Zero() };
    Vec3d m_bb_center{ Vec3d::Zero() };
    Vec3d m_center_offset{ Vec3d::Zero() };

    // values from RotationGizmo
    float m_radius{ 0.0f };
    float m_snap_coarse_in_radius{ 0.0f };
    float m_snap_coarse_out_radius{ 0.0f };
    float m_snap_fine_in_radius{ 0.0f };
    float m_snap_fine_out_radius{ 0.0f };

    GLModel         m_connector_shape;
    TriangleMesh    m_connector_mesh;
    // workaround for using of the clipping plane normal
    Vec3d           m_clp_normal{ Vec3d::Ones() };

    Vec3d           m_line_beg{ Vec3d::Zero() };
    Vec3d           m_line_end{ Vec3d::Zero() };

#if ENABLE_LEGACY_OPENGL_REMOVAL
    GLModel m_plane;
    GLModel m_grabber_connection;
    GLModel m_cut_line;
    GLModel m_cone;
    GLModel m_sphere;
    Vec3d   m_old_center;
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    bool m_keep_upper{ true };
    bool m_keep_lower{ true };
    bool m_place_on_cut_upper{ true };
    bool m_place_on_cut_lower{ true };
    bool m_rotate_upper{ false };
    bool m_rotate_lower{ false };

    bool m_hide_cut_plane{ false };
    bool m_connectors_editing{ false };

    double m_connector_depth_ratio{ 3.0 };
    double m_connector_size{ 2.5 };

    float m_label_width{ 150.0 };
    float m_control_width{ 200.0 };
    bool  m_imperial_units{ false };
    bool  force_update_clipper_on_render{false};

    mutable std::vector<bool> m_selected; // which pins are currently selected
    bool m_selection_empty      = true;
    bool m_wait_for_up_event    = false;

    bool m_has_invalid_connector{ false };

    Matrix3d m_rotation_matrix;
    Vec3d    m_rotations{ Vec3d::Zero() };

    enum class CutMode {
        cutPlanar
        , cutGrig
        //,cutRadial
        //,cutModular
    };

    enum class CutConnectorMode {
        Auto
        , Manual
    };

    std::vector<std::string> m_modes;
    size_t m_mode{ size_t(CutMode::cutPlanar) };

    std::vector<std::string> m_connector_modes;
    CutConnectorMode m_connector_mode{ CutConnectorMode::Manual };

    std::vector<std::string> m_connector_types;
    CutConnectorType m_connector_type;

    std::vector<std::string> m_connector_styles;
    size_t m_connector_style;

    std::vector<std::string> m_connector_shapes;
    size_t m_connector_shape_id;

    std::vector<std::string> m_axis_names;

public:
    GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    std::string get_tooltip() const override;
    bool unproject_on_cut_plane(const Vec2d& mouse_pos, std::pair<Vec3d, Vec3d>& pos_and_normal);
    bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);

    /// <summary>
    /// Drag of plane
    /// </summary>
    /// <param name="mouse_event">Keep information about mouse click</param>
    /// <returns>Return True when use the information otherwise False.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    void shift_cut_z(double delta);
    void rotate_vec3d_around_center(Vec3d& vec, const Vec3d& angles, const Vec3d& center);
    void put_connetors_on_cut_plane(const Vec3d& cp_normal, double cp_offset);
    void update_clipper();
    void update_clipper_on_render();
    void set_connectors_editing() { m_connectors_editing = true; }

    BoundingBoxf3   bounding_box() const;
    BoundingBoxf3   transformed_bounding_box(bool revert_move = false) const;

protected:
    bool on_init() override;
    void on_load(cereal::BinaryInputArchive& ar) override;
    void on_save(cereal::BinaryOutputArchive& ar) const override;
    std::string on_get_name() const override;
    void on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;
    void on_set_hover_id() override;
    bool on_is_activable() const override;
    Vec3d mouse_position_in_local_plane(Axis axis, const Linef3& mouse_ray) const;
    void on_dragging(const UpdateData& data) override;
    void on_start_dragging() override;
    void on_stop_dragging() override;
    void on_render() override;
    void on_render_for_picking() override;
    void on_render_input_window(float x, float y, float bottom_limit) override;

    bool wants_enter_leave_snapshots() const override       { return true; }
    std::string get_gizmo_entering_text() const override    { return _u8L("Entering Cut gizmo"); }
    std::string get_gizmo_leaving_text() const override     { return _u8L("Leaving Cut gizmo"); }
    std::string get_action_snapshot_name() override         { return _u8L("Cut gizmo editing"); }

private:
    void set_center(const Vec3d& center);
    bool render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx);
    bool render_double_input(const std::string& label, double& value_in);
    bool render_slicer_double_input(const std::string& label, double& value_in);
    void render_move_center_input(int axis);
    void render_connect_mode_radio_button(CutConnectorMode mode);
    bool render_revert_button(const std::string& label);
    void render_connect_type_radio_button(CutConnectorType type);
    Transform3d get_volume_transformation(const ModelVolume* volume) const;
    void render_connectors(bool picking);

    bool can_perform_cut() const;
    bool cut_line_processing() const;

    void render_cut_plane();
    void render_cut_center_graber(bool picking = false);
    void render_cut_line();
    void perform_cut(const Selection& selection);
    void set_center_pos(const Vec3d& center_pos, bool force = false);
    bool update_bb();
    void reset_connectors();
    void update_connector_shape();
    void update_model_object() const;
    bool process_cut_line(SLAGizmoEventType action, const Vec2d& mouse_position);
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoCut_hpp_
