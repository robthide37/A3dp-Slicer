#ifndef slic3r_GLGizmoCut_hpp_
#define slic3r_GLGizmoCut_hpp_

#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "GLGizmoMove.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ObjectID.hpp"

namespace Slic3r {
namespace GUI {

enum class SLAGizmoEventType : unsigned char;

class GLGizmoCenterMove : public GLGizmoMove3D
{
public:
    static const double Margin;
private:

    Vec3d m_min_pos{ Vec3d::Zero() };
    Vec3d m_max_pos{ Vec3d::Zero() };

public:
    GLGizmoCenterMove(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);
    std::string get_tooltip() const override;

protected:
    virtual void on_set_state() override;
    virtual void on_update(const UpdateData& data) override;

public:
    void            set_center_pos(const Vec3d& center_pos);
    BoundingBoxf3   bounding_box() const;
};


class GLGizmoCut3D : public GLGizmoBase
{
    GLGizmoRotate3D             m_rotation_gizmo;
    GLGizmoCenterMove           m_move_gizmo;

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLModel m_plane;
    float m_old_z{ 0.0f };
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    bool m_keep_upper{ true };
    bool m_keep_lower{ true };
    bool m_rotate_lower{ false };

    double m_connector_depth_ratio{ 1.5 };
    double m_connector_size{ 5.0 };

    float m_label_width{ 150.0 };
    float m_control_width{ 200.0 };
    bool  m_imperial_units{ false };

    enum CutMode {
        cutPlanar
        , cutByLine
        , cutGrig
        //,cutRadial
        //,cutModular
    };

    enum ConnectorMode {
        Auto
        , Manual
    };

    enum ConnectorType {
        Plug
        , Dowel
    };

    enum ConnectorStyle {
        Prizm
        , Frustrum
        //,Claw
    };

    enum ConnectorShape {
        Triangle
        , Square
        , Circle
        , Hexagon
        //,D-shape
    };

    std::vector<std::string> m_modes;
    size_t m_mode{ size_t(cutPlanar) };

    std::vector<std::string> m_connector_modes;
    ConnectorMode m_connector_mode{ Auto };

    std::vector<std::string> m_connector_types;
    ConnectorType m_connector_type{ Plug };

    std::vector<std::string> m_connector_styles;
    size_t m_connector_style{ size_t(Prizm) };

    std::vector<std::string> m_connector_shapes;
    size_t m_connector_shape{ size_t(Hexagon) };

    std::vector<std::string> m_axis_names;

public:
    GLGizmoCut3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    std::string get_tooltip() const override;
    bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);

    void shift_cut_z(double delta);

protected:
    bool on_init() override;
    void on_load(cereal::BinaryInputArchive& ar)  override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    void on_save(cereal::BinaryOutputArchive& ar) const override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    std::string on_get_name() const override;
    void on_set_state() override;
    bool on_is_activable() const override;
    void on_start_dragging() override;
    void on_update(const UpdateData& data) override;
    CommonGizmosDataID on_get_requirements() const override;
    void on_set_hover_id() override;
    void on_enable_grabber(unsigned int id) override;
    void on_disable_grabber(unsigned int id) override;
    bool on_is_activable() const override;
    void on_start_dragging() override;
    void on_stop_dragging() override;
    void on_update(const UpdateData& data) override;
    void on_render() override;
    void on_render_for_picking() override {
        m_rotation_gizmo.render_for_picking();
        m_move_gizmo.render_for_picking();
    }
    void on_render_input_window(float x, float y, float bottom_limit) override;


private:
    void set_center(const Vec3d& center);
    void render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx);
    void render_double_input(const std::string& label, double& value_in);
    void render_move_center_input(int axis);
    void render_rotation_input(int axis);
    void render_connect_mode_radio_button(ConnectorMode mode);
    void render_connect_type_radio_button(ConnectorType type);
    bool can_perform_cut() const;

    void render_cut_plane();
    void perform_cut(const Selection& selection);
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoCut_hpp_
