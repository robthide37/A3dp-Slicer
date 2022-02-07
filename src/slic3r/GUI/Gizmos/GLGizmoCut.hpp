#ifndef slic3r_GLGizmoCut_hpp_
#define slic3r_GLGizmoCut_hpp_

#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ObjectID.hpp"

namespace Slic3r {
namespace GUI {

enum class SLAGizmoEventType : unsigned char;

class GLGizmoCut : public GLGizmoBase
{
    static const double Offset;
    static const double Margin;

    double m_cut_z{ 0.0 };
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

    // struct CutContours
    // {
    //     TriangleMesh mesh;
    //     GLModel contours;
    //     double cut_z{ 0.0 };
    //     Vec3d position{ Vec3d::Zero() };
    //     Vec3d shift{ Vec3d::Zero() };
    //     ObjectID object_id;
    //     int instance_idx{ -1 };
    //     std::vector<ObjectID> volumes_idxs;
    // };

    // CutContours m_cut_contours;

public:
    GLGizmoCut(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);

    double get_cut_z() const { return m_cut_z; }
    void set_cut_z(double cut_z);

    std::string get_tooltip() const override;
    bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);

protected:
    bool on_init() override;
    void on_load(cereal::BinaryInputArchive& ar)  override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    void on_save(cereal::BinaryOutputArchive& ar) const override { ar(m_cut_z, m_keep_upper, m_keep_lower, m_rotate_lower); }
    std::string on_get_name() const override;
    void on_set_state() override;
    bool on_is_activable() const override;
    void on_start_dragging() override;
    void on_update(const UpdateData& data) override;
    void on_render() override;
    void on_render_for_picking() override;
    void on_render_input_window(float x, float y, float bottom_limit) override;
    CommonGizmosDataID on_get_requirements() const override;


private:
    void perform_cut(const Selection& selection);
    double calc_projection(const Linef3& mouse_ray) const;
    BoundingBoxf3 bounding_box() const;
    //void update_contours();
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoCut_hpp_
