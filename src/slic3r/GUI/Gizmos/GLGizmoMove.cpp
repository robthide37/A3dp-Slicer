// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMove.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"

#include <GL/glew.h>

#include <wx/utils.h> 

namespace Slic3r {
namespace GUI {

const double GLGizmoMove3D::Offset = 10.0;

GLGizmoMove3D::GLGizmoMove3D(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{}

std::string GLGizmoMove3D::get_tooltip() const
{
    const Selection& selection = m_parent.get_selection();
    bool show_position = selection.is_single_full_instance();
    const Vec3d& position = selection.get_bounding_box().center();

    if (m_hover_id == 0 || m_grabbers[0].dragging)
        return "X: " + format(show_position ? position(0) : m_displacement(0), 2);
    else if (m_hover_id == 1 || m_grabbers[1].dragging)
        return "Y: " + format(show_position ? position(1) : m_displacement(1), 2);
    else if (m_hover_id == 2 || m_grabbers[2].dragging)
        return "Z: " + format(show_position ? position(2) : m_displacement(2), 2);
    else
        return "";
}

bool GLGizmoMove3D::on_init()
{
    for (int i = 0; i < 3; ++i) {
        m_grabbers.push_back(Grabber());
    }

    m_shortcut_key = WXK_CONTROL_M;

    return true;
}

std::string GLGizmoMove3D::on_get_name() const
{
    return _u8L("Move");
}

bool GLGizmoMove3D::on_is_activable() const
{
    return !m_parent.get_selection().is_empty();
}

void GLGizmoMove3D::on_start_dragging()
{
    if (m_hover_id != -1) {
        m_displacement = Vec3d::Zero();
        const BoundingBoxf3& box = m_parent.get_selection().get_bounding_box();
        m_starting_drag_position = m_grabbers[m_hover_id].center;
        m_starting_box_center = box.center();
        m_starting_box_bottom_center = box.center();
        m_starting_box_bottom_center(2) = box.min(2);
    }
}

void GLGizmoMove3D::on_stop_dragging()
{
    m_displacement = Vec3d::Zero();
}

void GLGizmoMove3D::on_update(const UpdateData& data)
{
    if (m_hover_id == 0)
        m_displacement.x() = calc_projection(data);
    else if (m_hover_id == 1)
        m_displacement.y() = calc_projection(data);
    else if (m_hover_id == 2)
        m_displacement.z() = calc_projection(data);
}

void GLGizmoMove3D::on_render()
{
    if (!m_cone.is_initialized())
        m_cone.init_from(its_make_cone(1.0, 1.0, double(PI) / 18.0));

    const Selection& selection = m_parent.get_selection();

    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
    glsafe(::glEnable(GL_DEPTH_TEST));

    const BoundingBoxf3& box = selection.get_bounding_box();
    const Vec3d& center = box.center();

    // x axis
    m_grabbers[0].center = { box.max.x() + Offset, center.y(), center.z() };
    m_grabbers[0].color = AXES_COLOR[0];

    // y axis
    m_grabbers[1].center = { center.x(), box.max.y() + Offset, center.z() };
    m_grabbers[1].color = AXES_COLOR[1];

    // z axis
    m_grabbers[2].center = { center.x(), center.y(), box.max.z() + Offset };
    m_grabbers[2].color = AXES_COLOR[2];

    glsafe(::glLineWidth((m_hover_id != -1) ? 2.0f : 1.5f));

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    auto render_grabber_connection = [this, &center](unsigned int id) {
        if (m_grabbers[id].enabled) {
            if (!m_grabber_connections[id].model.is_initialized() || !m_grabber_connections[id].old_center.isApprox(center)) {
                m_grabber_connections[id].old_center = center;
                m_grabber_connections[id].model.reset();

                GLModel::Geometry init_data;
                init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3, GLModel::Geometry::EIndexType::USHORT };
                init_data.color = AXES_COLOR[id];
                init_data.reserve_vertices(2);
                init_data.reserve_indices(2);

                // vertices
                init_data.add_vertex((Vec3f)center.cast<float>());
                init_data.add_vertex((Vec3f)m_grabbers[id].center.cast<float>());

                // indices
                init_data.add_ushort_line(0, 1);

                m_grabber_connections[id].model.init_from(std::move(init_data));
            }

            m_grabber_connections[id].model.render();
        }
    };
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    if (m_hover_id == -1) {
#if ENABLE_GLBEGIN_GLEND_REMOVAL
        GLShaderProgram* shader = wxGetApp().get_shader("flat");
        if (shader != nullptr) {
            shader->start_using();
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

            // draw axes
            for (unsigned int i = 0; i < 3; ++i) {
#if ENABLE_GLBEGIN_GLEND_REMOVAL
                render_grabber_connection(i);
#else
                if (m_grabbers[i].enabled) {
                    glsafe(::glColor4fv(AXES_COLOR[i].data()));
                    ::glBegin(GL_LINES);
                    ::glVertex3dv(center.data());
                    ::glVertex3dv(m_grabbers[i].center.data());
                    glsafe(::glEnd());
                }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
            }

#if ENABLE_GLBEGIN_GLEND_REMOVAL
            shader->stop_using();
        }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

        // draw grabbers
        render_grabbers(box);
        for (unsigned int i = 0; i < 3; ++i) {
            if (m_grabbers[i].enabled)
                render_grabber_extension((Axis)i, box, false);
        }
    }
    else {
        // draw axis
#if ENABLE_GLBEGIN_GLEND_REMOVAL
        GLShaderProgram* shader = wxGetApp().get_shader("flat");
        if (shader != nullptr) {
            shader->start_using();
            render_grabber_connection(m_hover_id);
            shader->stop_using();
        }

        shader = wxGetApp().get_shader("gouraud_light");
#else
        glsafe(::glColor4fv(AXES_COLOR[m_hover_id].data()));
        ::glBegin(GL_LINES);
        ::glVertex3dv(center.data());
        ::glVertex3dv(m_grabbers[m_hover_id].center.data());
        glsafe(::glEnd());

        GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
        if (shader != nullptr) {
            shader->start_using();
            shader->set_uniform("emission_factor", 0.1f);
            // draw grabber
            float mean_size = (float)((box.size().x() + box.size().y() + box.size().z()) / 3.0);
            m_grabbers[m_hover_id].render(true, mean_size);
            shader->stop_using();
        }
        render_grabber_extension((Axis)m_hover_id, box, false);
    }
}

void GLGizmoMove3D::on_render_for_picking()
{
    glsafe(::glDisable(GL_DEPTH_TEST));

    const BoundingBoxf3& box = m_parent.get_selection().get_bounding_box();
    render_grabbers_for_picking(box);
    render_grabber_extension(X, box, true);
    render_grabber_extension(Y, box, true);
    render_grabber_extension(Z, box, true);
}

double GLGizmoMove3D::calc_projection(const UpdateData& data) const
{
    double projection = 0.0;

    Vec3d starting_vec = m_starting_drag_position - m_starting_box_center;
    double len_starting_vec = starting_vec.norm();
    if (len_starting_vec != 0.0) {
        Vec3d mouse_dir = data.mouse_ray.unit_vector();
        // finds the intersection of the mouse ray with the plane parallel to the camera viewport and passing throught the starting position
        // use ray-plane intersection see i.e. https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection algebric form
        // in our case plane normal and ray direction are the same (orthogonal view)
        // when moving to perspective camera the negative z unit axis of the camera needs to be transformed in world space and used as plane normal
        Vec3d inters = data.mouse_ray.a + (m_starting_drag_position - data.mouse_ray.a).dot(mouse_dir) / mouse_dir.squaredNorm() * mouse_dir;
        // vector from the starting position to the found intersection
        Vec3d inters_vec = inters - m_starting_drag_position;

        // finds projection of the vector along the staring direction
        projection = inters_vec.dot(starting_vec.normalized());
    }

    if (wxGetKeyState(WXK_SHIFT))
        projection = m_snap_step * (double)std::round(projection / m_snap_step);

    return projection;
}

void GLGizmoMove3D::render_grabber_extension(Axis axis, const BoundingBoxf3& box, bool picking)
{
    float mean_size = float((box.size().x() + box.size().y() + box.size().z()) / 3.0);
    double size = m_dragging ? double(m_grabbers[axis].get_dragging_half_size(mean_size)) : double(m_grabbers[axis].get_half_size(mean_size));

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLShaderProgram* shader = wxGetApp().get_shader(picking ? "flat" : "gouraud_light");
#else
    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    if (shader == nullptr)
        return;

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    m_cone.set_color((!picking && m_hover_id != -1) ? complementary(m_grabbers[axis].color) : m_grabbers[axis].color);
    shader->start_using();
    shader->set_uniform("emission_factor", 0.1f);
#else
    m_cone.set_color(-1, (!picking && m_hover_id != -1) ? complementary(m_grabbers[axis].color) : m_grabbers[axis].color);
    if (!picking) {
        shader->start_using();
        shader->set_uniform("emission_factor", 0.1f);
    }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    glsafe(::glPushMatrix());
    glsafe(::glTranslated(m_grabbers[axis].center.x(), m_grabbers[axis].center.y(), m_grabbers[axis].center.z()));
    if (axis == X)
        glsafe(::glRotated(90.0, 0.0, 1.0, 0.0));
    else if (axis == Y)
        glsafe(::glRotated(-90.0, 1.0, 0.0, 0.0));

    glsafe(::glTranslated(0.0, 0.0, 2.0 * size));
    glsafe(::glScaled(0.75 * size, 0.75 * size, 3.0 * size));
    m_cone.render();
    glsafe(::glPopMatrix());

#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    if (! picking)
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL
        shader->stop_using();
}

} // namespace GUI
} // namespace Slic3r
