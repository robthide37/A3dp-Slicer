#ifndef slic3r_GLGizmoCut_hpp_
#define slic3r_GLGizmoCut_hpp_

#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ObjectID.hpp"

namespace Slic3r {
namespace GUI {

class GLGizmoCut : public GLGizmoBase
{
    static const double Offset;
    static const double Margin;

    double m_cut_z{ 0.0 };
    double m_rotation_x{ 0.0 };
    double m_rotation_y{ 0.0 };
    double m_rotation_z{ 0.0 };
    double m_max_z{ 0.0 };
    double m_start_z{ 0.0 };
    Vec3d m_drag_pos;
    Vec3d m_drag_center;
    bool m_keep_upper{ true };
    bool m_keep_lower{ true };
    bool m_rotate_lower{ false };
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLModel m_plane;
    GLModel m_grabber_connection;
    float m_old_z{ 0.0f };
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    double m_connector_depth_ratio{ 1.5 };
    double m_connector_size{ 5.0 };

    float m_label_width{ 150.0 };
    float m_control_width{ 200.0 };
    bool  m_imperial_units{ false };

public:
    enum CutMode {
         cutPlanar = 0
        ,cutGrig
        //,cutRadial
        //,cutModular
    };

    enum ConnectorType {
         Plug = 0
        ,Dowel
    };

    enum ConnectorStyle {
         Prizm = 0
        ,Frustrum
        //,Claw
    };

    enum ConnectorShape {
         Triangle = 0
        ,Square
        ,Circle
        ,Hexagon
//        ,D-shape
    };

private:

    std::vector<std::string> m_modes;
    size_t m_mode{ size_t(cutPlanar) };

    std::vector<std::string> m_connector_types;
    ConnectorType m_connector_type{ Plug };

    std::vector<std::string> m_connector_styles;
    size_t m_connector_style{ size_t(Prizm) };

    std::vector<std::string> m_connector_shapes;
    size_t m_connector_shape{ size_t(Hexagon) };

    struct CutContours
    {
        TriangleMesh mesh;
        GLModel contours;
        double cut_z{ 0.0 };
        Vec3d position{ Vec3d::Zero() };
        Vec3d shift{ Vec3d::Zero() };
        ObjectID object_id;
        int instance_idx{ -1 };
        std::vector<ObjectID> volumes_idxs;
    };

    CutContours m_cut_contours;

public:
    GLGizmoCut(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    double get_cut_z() const { return m_cut_z; }
    void set_cut_z(double cut_z);

    std::string get_tooltip() const override;

protected:
    virtual bool on_init() override;
    virtual void on_load(cereal::BinaryInputArchive& ar)  override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    virtual void on_save(cereal::BinaryOutputArchive& ar) const override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    virtual std::string on_get_name() const override;
    virtual void on_set_state() override;
    virtual bool on_is_activable() const override;
    virtual void on_start_dragging() override;
    virtual void on_update(const UpdateData& data) override;
    virtual void on_render() override;
    virtual void on_render_for_picking() override;
    virtual void on_render_input_window(float x, float y, float bottom_limit) override;

    void render_combo(const std::string& label, const std::vector<std::string>& lines, size_t& selection_idx);
    void render_double_input(const std::string& label, double& value_in);
    void render_rotation_input(const std::string& label, double& value_in);
    void render_radio_button(ConnectorType type);

private:
    void perform_cut(const Selection& selection);
    double calc_projection(const Linef3& mouse_ray) const;
    BoundingBoxf3 bounding_box() const;
    void update_contours();
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoCut_hpp_
