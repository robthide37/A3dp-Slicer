#ifndef slic3r_GLGizmoCut_hpp_
#define slic3r_GLGizmoCut_hpp_

#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLSelectionRectangle.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Model.hpp"

namespace Slic3r {

enum class CutConnectorType : int;
class ModelVolume;
struct CutConnectorAttributes;

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

    Vec3d m_plane_center{ Vec3d::Zero() };
    // data to check position of the cut palne center on gizmo activation
    Vec3d m_min_pos{ Vec3d::Zero() };
    Vec3d m_max_pos{ Vec3d::Zero() };
    Vec3d m_bb_center{ Vec3d::Zero() };
    Vec3d m_center_offset{ Vec3d::Zero() };

    // values from RotationGizmo
    double m_radius{ 0.0 };
    double m_grabber_radius{ 0.0 };
    double m_grabber_connection_len{ 0.0 };

    double m_snap_coarse_in_radius{ 0.0 };
    double m_snap_coarse_out_radius{ 0.0 };
    double m_snap_fine_in_radius{ 0.0 };
    double m_snap_fine_out_radius{ 0.0 };

    // dragging angel in hovered axes
    Transform3d m_start_dragging_m{ Transform3d::Identity() };
    double m_angle{ 0.0 };

    TriangleMesh    m_connector_mesh;
    // workaround for using of the clipping plane normal
    Vec3d           m_clp_normal{ Vec3d::Ones() };

    Vec3d           m_line_beg{ Vec3d::Zero() };
    Vec3d           m_line_end{ Vec3d::Zero() };

    Vec2d           m_ldown_mouse_position{ Vec2d::Zero() };

    GLModel m_plane;
    GLModel m_grabber_connection;
    GLModel m_cut_line;

    PickingModel m_sphere;
    PickingModel m_cone;
    std::map<CutConnectorAttributes, PickingModel> m_shapes;
    std::vector<std::shared_ptr<SceneRaycasterItem>> m_raycasters;

    GLModel m_circle;
    GLModel m_scale;
    GLModel m_snap_radii;
    GLModel m_reference_radius;
    GLModel m_angle_arc;

    Vec3d   m_old_center;

    bool m_keep_upper{ true };
    bool m_keep_lower{ true };
    bool m_place_on_cut_upper{ true };
    bool m_place_on_cut_lower{ false };
    bool m_rotate_upper{ false };
    bool m_rotate_lower{ false };

    bool m_hide_cut_plane{ false };
    bool m_connectors_editing{ false };

    float m_connector_depth_ratio{ 3.f };
    float m_connector_size{ 2.5f };

    float m_connector_depth_ratio_tolerance{ 0.1f };
    float m_connector_size_tolerance{ 0.f };

    float m_label_width{ 150.0 };
    float m_control_width{ 200.0 };
    bool  m_imperial_units{ false };
    bool  force_update_clipper_on_render{false};

    mutable std::vector<bool> m_selected; // which pins are currently selected
    int  m_selected_count{ 0 };

    GLSelectionRectangle m_selection_rectangle;

    bool m_has_invalid_connector{ false };

    bool                                        m_show_shortcuts{ false };
    std::vector<std::pair<wxString, wxString>>  m_shortcuts;

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
    bool unproject_on_cut_plane(const Vec2d& mouse_pos, std::pair<Vec3d, Vec3d>& pos_and_normal, Vec3d& pos_world);
    bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);

    bool is_in_editing_mode() const override { return m_connectors_editing; }
    bool is_selection_rectangle_dragging() const override { return m_selection_rectangle.is_dragging(); }

    /// <summary>
    /// Drag of plane
    /// </summary>
    /// <param name="mouse_event">Keep information about mouse click</param>
    /// <returns>Return True when use the information otherwise False.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    void shift_cut_z(double delta);
    void rotate_vec3d_around_plane_center(Vec3d&vec);
    void put_connectors_on_cut_plane(const Vec3d& cp_normal, double cp_offset);
    void update_clipper();
    void update_clipper_on_render();
    void set_connectors_editing() { m_connectors_editing = true; }
    void invalidate_cut_plane();

    BoundingBoxf3   bounding_box() const;
    BoundingBoxf3   transformed_bounding_box(bool revert_move = false) const;

protected:
    bool               on_init() override;
    void               on_load(cereal::BinaryInputArchive&ar) override;
    void               on_save(cereal::BinaryOutputArchive&ar) const override;
    std::string        on_get_name() const override;
    void               on_set_state() override;
    CommonGizmosDataID on_get_requirements() const override;
    void               on_set_hover_id() override;
    bool               on_is_activable() const override;
    Vec3d              mouse_position_in_local_plane(Axis axis, const Linef3&mouse_ray) const;
    void               dragging_grabber_z(const GLGizmoBase::UpdateData &data);
    void               dragging_grabber_xy(const GLGizmoBase::UpdateData &data);
    void               dragging_connector(const GLGizmoBase::UpdateData &data);
    void               on_dragging(const UpdateData&data) override;
    void               on_start_dragging() override;
    void               on_stop_dragging() override;
    void               on_render() override;

    void render_debug_input_window();
    void adjust_window_position(float x, float y, float bottom_limit);
    void unselect_all_connectors();
    void select_all_connectors();
    void render_shortcuts();
    void apply_selected_connectors(std::function<void(size_t idx)> apply_fn);
    void render_connectors_input_window(CutConnectors &connectors);
    void render_build_size();
    void reset_cut_plane();
    void set_connectors_editing(bool connectors_editing);
    void render_cut_plane_input_window(CutConnectors &connectors);
    void init_input_window_data(CutConnectors &connectors);
    void render_input_window_warning() const;
    bool add_connector(CutConnectors&connectors, const Vec2d&mouse_position);
    bool delete_selected_connectors(CutConnectors&connectors);
    void select_connector(int idx, bool select);
    bool is_selection_changed(bool alt_down, bool shift_down);
    void process_selection_rectangle(CutConnectors &connectors);

    virtual void on_register_raycasters_for_picking() override;
    virtual void on_unregister_raycasters_for_picking() override;
    void update_raycasters_for_picking();
    void set_volumes_picking_state(bool state);
    void update_raycasters_for_picking_transform();

    void on_render_input_window(float x, float y, float bottom_limit) override;

    bool wants_enter_leave_snapshots() const override       { return true; }
    std::string get_gizmo_entering_text() const override    { return _u8L("Entering Cut gizmo"); }
    std::string get_gizmo_leaving_text() const override     { return _u8L("Leaving Cut gizmo"); }
    std::string get_action_snapshot_name() override         { return _u8L("Cut gizmo editing"); }

private:
    void set_center(const Vec3d& center);
    bool render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx);
    bool render_double_input(const std::string& label, double& value_in);
    bool render_slider_double_input(const std::string& label, float& value_in, float& tolerance_in);
    void render_move_center_input(int axis);
    void render_connect_mode_radio_button(CutConnectorMode mode);
    bool render_reset_button(const std::string& label_id, const std::string& tooltip) const;
    bool render_connect_type_radio_button(CutConnectorType type);
    Transform3d get_volume_transformation(const ModelVolume* volume) const;
    void render_connectors();

    bool can_perform_cut() const;
    void apply_connectors_in_model(ModelObject* mo, const bool has_connectors, bool &create_dowels_as_separate_object);
    bool cut_line_processing() const;
    void discard_cut_line_processing();

    void render_cut_plane();
    void render_model(GLModel& model, const ColorRGBA& color, Transform3d view_model_matrix);
    void render_line(GLModel& line_model, const ColorRGBA& color, Transform3d view_model_matrix, float width);
    void render_rotation_snapping(Axis axis, const ColorRGBA& color);
    void render_grabber_connection(const ColorRGBA& color, Transform3d view_matrix);
    void render_cut_plane_grabbers();
    void render_cut_line();
    void perform_cut(const Selection&selection);
    void set_center_pos(const Vec3d&center_pos, bool force = false);
    bool update_bb();
    void init_picking_models();
    void init_rendering_items();
    void render_clipper_cut();
    void clear_selection();
    void reset_connectors();
    void init_connector_shapes();
    void update_connector_shape();
    void validate_connector_settings();
    void update_model_object();
    bool process_cut_line(SLAGizmoEventType action, const Vec2d& mouse_position);
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoCut_hpp_
