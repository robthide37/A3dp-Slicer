#ifndef slic3r_GLGizmoMove_hpp_
#define slic3r_GLGizmoMove_hpp_

#include "GLGizmoBase.hpp"


namespace Slic3r {
namespace GUI {

class GLGizmoMove3D : public GLGizmoBase
{
    static const double Offset;

    Vec3d m_displacement{ Vec3d::Zero() };
    double m_snap_step{ 1.0 };
    Vec3d m_starting_drag_position{ Vec3d::Zero() };
    Vec3d m_starting_box_center{ Vec3d::Zero() };
    Vec3d m_starting_box_bottom_center{ Vec3d::Zero() };

    GLModel m_cone;
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    struct GrabberConnection
    {
        GLModel model;
        Vec3d old_center{ Vec3d::Zero() };
    };
    std::array<GrabberConnection, 3> m_grabber_connections;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

public:
    GLGizmoMove3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);
    virtual ~GLGizmoMove3D() = default;

    double get_snap_step(double step) const { return m_snap_step; }
    void set_snap_step(double step) { m_snap_step = step; }

    const Vec3d& get_displacement() const { return m_displacement; }

    std::string get_tooltip() const override;

protected:
    virtual bool on_init() override;
    virtual std::string on_get_name() const override;
    virtual bool on_is_activable() const override;
    virtual void on_start_dragging() override;
    virtual void on_stop_dragging() override;
    virtual void on_update(const UpdateData& data) override;
    virtual void on_render() override;
    virtual void on_render_for_picking() override;

private:
    double calc_projection(const UpdateData& data) const;
    void render_grabber_extension(Axis axis, const BoundingBoxf3& box, bool picking) const;
};



} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoMove_hpp_
