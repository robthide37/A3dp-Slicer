#include "libslic3r/libslic3r.h"

#include "CoordAxes.hpp"
#include "GUI_App.hpp"
#include "3DScene.hpp"

#include <GL/glew.h>

#if ENABLE_WORLD_COORDINATE_SHOW_AXES

namespace Slic3r {
namespace GUI {

const float CoordAxes::DefaultStemRadius = 0.5f;
const float CoordAxes::DefaultStemLength = 25.0f;
const float CoordAxes::DefaultTipRadius = 2.5f * CoordAxes::DefaultStemRadius;
const float CoordAxes::DefaultTipLength = 5.0f;

void CoordAxes::render(float emission_factor)
{
    auto render_axis = [this](const Transform3f& transform) {
        glsafe(::glPushMatrix());
        glsafe(::glMultMatrixf(transform.data()));
        m_arrow.render();
        glsafe(::glPopMatrix());
    };

    if (!m_arrow.is_initialized())
        m_arrow.init_from(stilized_arrow(16, m_tip_radius, m_tip_length, m_stem_radius, m_stem_length));

    GLShaderProgram* curr_shader = wxGetApp().get_current_shader();
    bool shader_differs = (curr_shader == nullptr || curr_shader->get_name() != "gouraud_light");

    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    if (shader_differs) {
        if (curr_shader != nullptr)
            curr_shader->stop_using();
        shader->start_using();
    }
    shader->set_uniform("emission_factor", emission_factor);

    // x axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::X());
#else
    m_arrow.set_color(-1, ColorRGBA::X());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    render_axis(Geometry::assemble_transform(m_origin, { 0.0, 0.5 * M_PI, 0.0 }).cast<float>());

    // y axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Y());
#else
    m_arrow.set_color(-1, ColorRGBA::Y());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    render_axis(Geometry::assemble_transform(m_origin, { -0.5 * M_PI, 0.0, 0.0 }).cast<float>());

    // z axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Z());
#else
    m_arrow.set_color(-1, ColorRGBA::Z());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    render_axis(Geometry::assemble_transform(m_origin).cast<float>());

    if (shader_differs) {
        shader->stop_using();
        if (curr_shader != nullptr)
            curr_shader->start_using();
    }
}

} // GUI
} // Slic3r

#endif // ENABLE_WORLD_COORDINATE_SHOW_AXES
