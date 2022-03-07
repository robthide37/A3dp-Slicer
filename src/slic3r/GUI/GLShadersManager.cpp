#include "libslic3r/libslic3r.h"
#include "libslic3r/Platform.hpp"
#include "GLShadersManager.hpp"
#include "3DScene.hpp"
#include "GUI_App.hpp"

#include <cassert>
#include <algorithm>
#include <string_view>
using namespace std::literals;

#include <GL/glew.h>

namespace Slic3r {

std::pair<bool, std::string> GLShadersManager::init()
{
    std::string error;

    auto append_shader = [this, &error](const std::string& name, const GLShaderProgram::ShaderFilenames& filenames, 
        const std::initializer_list<std::string_view> &defines = {}) {
        m_shaders.push_back(std::make_unique<GLShaderProgram>());
        if (!m_shaders.back()->init_from_files(name, filenames, defines)) {
            error += name + "\n";
            // if any error happens while initializating the shader, we remove it from the list
            m_shaders.pop_back();
            return false;
        }
        return true;
    };

    assert(m_shaders.empty());

    bool valid = true;

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    // basic shader, used to render all what was previously rendered using the immediate mode
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("flat_attr", { "flat_attr.vs", "flat.fs" });
#else
    valid &= append_shader("flat", { "flat.vs", "flat.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // basic shader for textures, used to render textures
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("flat_texture_attr", { "flat_texture_attr.vs", "flat_texture.fs" });
#else
    valid &= append_shader("flat_texture", { "flat_texture.vs", "flat_texture.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // used to render 3D scene background
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("background_attr", { "background_attr.vs", "background.fs" });
#else
    valid &= append_shader("background", { "background.vs", "background.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
#if ENABLE_SHOW_TOOLPATHS_COG
    // used to render toolpaths center of gravity
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("toolpaths_cog_attr", { "toolpaths_cog_attr.vs", "toolpaths_cog.fs" });
#else
    valid &= append_shader("toolpaths_cog", { "toolpaths_cog.vs", "toolpaths_cog.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
#endif // ENABLE_SHOW_TOOLPATHS_COG
    // used to render bed axes and model, selection hints, gcode sequential view marker model, preview shells, options in gcode preview
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("gouraud_light_attr", { "gouraud_light_attr.vs", "gouraud_light.fs" });
#else
    valid &= append_shader("gouraud_light", { "gouraud_light.vs", "gouraud_light.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // used to render printbed
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    valid &= append_shader("printbed_attr", { "printbed_attr.vs", "printbed.fs" });
#else
    valid &= append_shader("printbed", { "printbed.vs", "printbed.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // used to render options in gcode preview
    if (GUI::wxGetApp().is_gl_version_greater_or_equal_to(3, 3)) {
#if ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
        valid &= append_shader("gouraud_light_instanced_attr", { "gouraud_light_instanced_attr.vs", "gouraud_light_instanced.fs" });
#else
        valid &= append_shader("gouraud_light_instanced", { "gouraud_light_instanced.vs", "gouraud_light_instanced.fs" });
#endif // ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    }
#if !ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // used to render extrusion and travel paths as lines in gcode preview
    valid &= append_shader("toolpaths_lines", { "toolpaths_lines.vs", "toolpaths_lines.fs" });
#endif // !ENABLE_GLBEGIN_GLEND_SHADERS_ATTRIBUTES
    // used to render objects in 3d editor
    valid &= append_shader("gouraud", { "gouraud.vs", "gouraud.fs" }
#if ENABLE_ENVIRONMENT_MAP
        , { "ENABLE_ENVIRONMENT_MAP"sv }
#endif // ENABLE_ENVIRONMENT_MAP
        );
    // used to render variable layers heights in 3d editor
    valid &= append_shader("variable_layer_height", { "variable_layer_height.vs", "variable_layer_height.fs" });
    // used to render highlight contour around selected triangles inside the multi-material gizmo
    valid &= append_shader("mm_contour", { "mm_contour.vs", "mm_contour.fs" });
    // Used to render painted triangles inside the multi-material gizmo. Triangle normals are computed inside fragment shader.
    // For Apple's on Arm CPU computed triangle normals inside fragment shader using dFdx and dFdy has the opposite direction.
    // Because of this, objects had darker colors inside the multi-material gizmo.
    // Based on https://stackoverflow.com/a/66206648, the similar behavior was also spotted on some other devices with Arm CPU.
    // Since macOS 12 (Monterey), this issue with the opposite direction on Apple's Arm CPU seems to be fixed, and computed
    // triangle normals inside fragment shader have the right direction.
    if (platform_flavor() == PlatformFlavor::OSXOnArm && wxPlatformInfo::Get().GetOSMajorVersion() < 12)
        valid &= append_shader("mm_gouraud", {"mm_gouraud.vs", "mm_gouraud.fs"}, {"FLIP_TRIANGLE_NORMALS"sv});
    else
        valid &= append_shader("mm_gouraud", {"mm_gouraud.vs", "mm_gouraud.fs"});

    return { valid, error };
}

void GLShadersManager::shutdown()
{
    m_shaders.clear();
}

GLShaderProgram* GLShadersManager::get_shader(const std::string& shader_name)
{
    auto it = std::find_if(m_shaders.begin(), m_shaders.end(), [&shader_name](std::unique_ptr<GLShaderProgram>& p) { return p->get_name() == shader_name; });
    return (it != m_shaders.end()) ? it->get() : nullptr;
}

GLShaderProgram* GLShadersManager::get_current_shader()
{
    GLint id = 0;
    glsafe(::glGetIntegerv(GL_CURRENT_PROGRAM, &id));
    if (id == 0)
        return nullptr;

    auto it = std::find_if(m_shaders.begin(), m_shaders.end(), [id](std::unique_ptr<GLShaderProgram>& p) { return static_cast<GLint>(p->get_id()) == id; });
    return (it != m_shaders.end()) ? it->get() : nullptr;
}

} // namespace Slic3r

