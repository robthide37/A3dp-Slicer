#include "libslic3r/libslic3r.h"

#include "CoordAxes.hpp"
#include "GUI_App.hpp"
#include "3DScene.hpp"
#if ENABLE_LEGACY_OPENGL_REMOVAL
#include "Plater.hpp"
#include "Camera.hpp"
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#include <GL/glew.h>

#if ENABLE_WORLD_COORDINATE

namespace Slic3r {
namespace GUI {

const float CoordAxes::DefaultStemRadius = 0.5f;
const float CoordAxes::DefaultStemLength = 25.0f;
const float CoordAxes::DefaultTipRadius = 2.5f * CoordAxes::DefaultStemRadius;
const float CoordAxes::DefaultTipLength = 5.0f;

#if ENABLE_LEGACY_OPENGL_REMOVAL
void CoordAxes::render(const Transform3d& trafo, float emission_factor)
#else
void CoordAxes::render(float emission_factor)
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
{
#if ENABLE_LEGACY_OPENGL_REMOVAL
    auto render_axis = [this](GLShaderProgram& shader, const Transform3d& transform) {
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d& view_matrix = camera.get_view_matrix();
        const Transform3d matrix = view_matrix * transform;
        shader.set_uniform("view_model_matrix", matrix);
        shader.set_uniform("projection_matrix", camera.get_projection_matrix());
        shader.set_uniform("view_normal_matrix", (Matrix3d)(view_matrix.matrix().block(0, 0, 3, 3) * transform.matrix().block(0, 0, 3, 3).inverse().transpose()));
        m_arrow.render();
#else
    auto render_axis = [this](const Transform3f& transform) {
        glsafe(::glPushMatrix());
        glsafe(::glMultMatrixf(transform.data()));
        m_arrow.render();
        glsafe(::glPopMatrix());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    };

    if (!m_arrow.is_initialized())
        m_arrow.init_from(stilized_arrow(16, m_tip_radius, m_tip_length, m_stem_radius, m_stem_length));

    GLShaderProgram* curr_shader = wxGetApp().get_current_shader();
    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    if (curr_shader != nullptr)
        curr_shader->stop_using();

    shader->start_using();
    shader->set_uniform("emission_factor", emission_factor);

    // x axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::X());
    render_axis(*shader, trafo * Geometry::assemble_transform(m_origin, { 0.0, 0.5 * M_PI, 0.0 }));
#else
    m_arrow.set_color(-1, ColorRGBA::X());
    render_axis(Geometry::assemble_transform(m_origin, { 0.0, 0.5 * M_PI, 0.0 }).cast<float>());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    // y axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Y());
    render_axis(*shader, trafo * Geometry::assemble_transform(m_origin, { -0.5 * M_PI, 0.0, 0.0 }));
#else
    m_arrow.set_color(-1, ColorRGBA::Y());
    render_axis(Geometry::assemble_transform(m_origin, { -0.5 * M_PI, 0.0, 0.0 }).cast<float>());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    // z axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Z());
    render_axis(*shader, trafo * Geometry::assemble_transform(m_origin));
#else
    m_arrow.set_color(-1, ColorRGBA::Z());
    render_axis(Geometry::assemble_transform(m_origin).cast<float>());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    shader->stop_using();
    if (curr_shader != nullptr)
        curr_shader->start_using();
}

} // GUI
} // Slic3r

#endif // ENABLE_WORLD_COORDINATE
