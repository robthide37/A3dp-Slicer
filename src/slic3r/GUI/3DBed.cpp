#include "libslic3r/libslic3r.h"

#include "3DBed.hpp"

#include "libslic3r/Polygon.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Geometry/Circle.hpp"
#include "libslic3r/Tesselate.hpp"
#include "libslic3r/PresetBundle.hpp"

#include "GUI_App.hpp"
#include "GLCanvas3D.hpp"
#if ENABLE_GL_SHADERS_ATTRIBUTES
#include "Plater.hpp"
#include "Camera.hpp"
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

#include <GL/glew.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/log/trivial.hpp>

static const float GROUND_Z = -0.02f;
static const Slic3r::ColorRGBA DEFAULT_MODEL_COLOR             = Slic3r::ColorRGBA::DARK_GRAY();
static const Slic3r::ColorRGBA PICKING_MODEL_COLOR             = Slic3r::ColorRGBA::BLACK();
static const Slic3r::ColorRGBA DEFAULT_SOLID_GRID_COLOR        = { 0.9f, 0.9f, 0.9f, 1.0f };
static const Slic3r::ColorRGBA DEFAULT_TRANSPARENT_GRID_COLOR  = { 0.9f, 0.9f, 0.9f, 0.6f };

namespace Slic3r {
namespace GUI {

#if !ENABLE_LEGACY_OPENGL_REMOVAL
bool GeometryBuffer::set_from_triangles(const std::vector<Vec2f> &triangles, float z)
{
    if (triangles.empty()) {
        m_vertices.clear();
        return false;
    }

    assert(triangles.size() % 3 == 0);
    m_vertices = std::vector<Vertex>(triangles.size(), Vertex());

    Vec2f min = triangles.front();
    Vec2f max = min;

    for (size_t v_count = 0; v_count < triangles.size(); ++ v_count) {
        const Vec2f &p = triangles[v_count];
        Vertex      &v = m_vertices[v_count];
        v.position   = Vec3f(p.x(), p.y(), z);
        v.tex_coords = p;
        min = min.cwiseMin(p).eval();
        max = max.cwiseMax(p).eval();
    }

    Vec2f size = max - min;
    if (size.x() != 0.f && size.y() != 0.f) {
        Vec2f inv_size = size.cwiseInverse();
        inv_size.y() *= -1;
        for (Vertex& v : m_vertices) {
            v.tex_coords -= min;
            v.tex_coords.x() *= inv_size.x();
            v.tex_coords.y() *= inv_size.y();
        }
    }

    return true;
}

bool GeometryBuffer::set_from_lines(const Lines& lines, float z)
{
    m_vertices.clear();

    unsigned int v_size = 2 * (unsigned int)lines.size();
    if (v_size == 0)
        return false;

    m_vertices = std::vector<Vertex>(v_size, Vertex());

    unsigned int v_count = 0;
    for (const Line& l : lines) {
        Vertex& v1 = m_vertices[v_count];
        v1.position[0] = unscale<float>(l.a(0));
        v1.position[1] = unscale<float>(l.a(1));
        v1.position[2] = z;
        ++v_count;

        Vertex& v2 = m_vertices[v_count];
        v2.position[0] = unscale<float>(l.b(0));
        v2.position[1] = unscale<float>(l.b(1));
        v2.position[2] = z;
        ++v_count;
    }

    return true;
}

const float* GeometryBuffer::get_vertices_data() const
{
    return (m_vertices.size() > 0) ? (const float*)m_vertices.data() : nullptr;
}
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

#if !ENABLE_WORLD_COORDINATE
const float Bed3D::Axes::DefaultStemRadius = 0.5f;
const float Bed3D::Axes::DefaultStemLength = 25.0f;
const float Bed3D::Axes::DefaultTipRadius = 2.5f * Bed3D::Axes::DefaultStemRadius;
const float Bed3D::Axes::DefaultTipLength = 5.0f;

void Bed3D::Axes::render()
{
#if ENABLE_GL_SHADERS_ATTRIBUTES
    auto render_axis = [this](GLShaderProgram* shader, const Transform3d& transform) {
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d matrix = camera.get_view_matrix() * transform;
        shader->set_uniform("view_model_matrix", matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        shader->set_uniform("normal_matrix",     (Matrix3d)matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
#else
    auto render_axis = [this](const Transform3f& transform) {
        glsafe(::glPushMatrix());
        glsafe(::glMultMatrixf(transform.data()));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        m_arrow.render();
#if !ENABLE_GL_SHADERS_ATTRIBUTES
        glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
    };

    if (!m_arrow.is_initialized())
        m_arrow.init_from(stilized_arrow(16, DefaultTipRadius, DefaultTipLength, DefaultStemRadius, m_stem_length));

    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;

    glsafe(::glEnable(GL_DEPTH_TEST));

    shader->start_using();
    shader->set_uniform("emission_factor", 0.0f);

    // x axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::X());
#else
    m_arrow.set_color(-1, ColorRGBA::X());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_GL_SHADERS_ATTRIBUTES
    render_axis(shader, Geometry::assemble_transform(m_origin, { 0.0, 0.5 * M_PI, 0.0 }));
#else
    render_axis(Geometry::assemble_transform(m_origin, { 0.0, 0.5 * M_PI, 0.0 }).cast<float>());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    // y axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Y());
#else
    m_arrow.set_color(-1, ColorRGBA::Y());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_GL_SHADERS_ATTRIBUTES
    render_axis(shader, Geometry::assemble_transform(m_origin, { -0.5 * M_PI, 0.0, 0.0 }));
#else
    render_axis(Geometry::assemble_transform(m_origin, { -0.5 * M_PI, 0.0, 0.0 }).cast<float>());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    // z axis
#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_arrow.set_color(ColorRGBA::Z());
#else
    m_arrow.set_color(-1, ColorRGBA::Z());
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
#if ENABLE_GL_SHADERS_ATTRIBUTES
    render_axis(shader, Geometry::assemble_transform(m_origin));
#else
    render_axis(Geometry::assemble_transform(m_origin).cast<float>());
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    shader->stop_using();

    glsafe(::glDisable(GL_DEPTH_TEST));
}
#endif // !ENABLE_WORLD_COORDINATE

bool Bed3D::set_shape(const Pointfs& bed_shape, const double max_print_height, const std::string& custom_texture, const std::string& custom_model, bool force_as_custom)
{
    auto check_texture = [](const std::string& texture) {
        boost::system::error_code ec; // so the exists call does not throw (e.g. after a permission problem)
        return !texture.empty() && (boost::algorithm::iends_with(texture, ".png") || boost::algorithm::iends_with(texture, ".svg")) && boost::filesystem::exists(texture, ec);
    };

    auto check_model = [](const std::string& model) {
        boost::system::error_code ec;
        return !model.empty() && boost::algorithm::iends_with(model, ".stl") && boost::filesystem::exists(model, ec);
    };

    Type type;
    std::string model;
    std::string texture;
    if (force_as_custom)
        type = Type::Custom;
    else {
        auto [new_type, system_model, system_texture] = detect_type(bed_shape);
        type = new_type;
        model = system_model;
        texture = system_texture;
    }

    std::string texture_filename = custom_texture.empty() ? texture : custom_texture;
    if (! texture_filename.empty() && ! check_texture(texture_filename)) {
        BOOST_LOG_TRIVIAL(error) << "Unable to load bed texture: " << texture_filename;
        texture_filename.clear();
    }

    std::string model_filename = custom_model.empty() ? model : custom_model;
    if (! model_filename.empty() && ! check_model(model_filename)) {
        BOOST_LOG_TRIVIAL(error) << "Unable to load bed model: " << model_filename;
        model_filename.clear();
    }

    
    if (m_build_volume.bed_shape() == bed_shape && m_build_volume.max_print_height() == max_print_height && m_type == type && m_texture_filename == texture_filename && m_model_filename == model_filename)
        // No change, no need to update the UI.
        return false;

    m_type = type;
    m_build_volume = BuildVolume { bed_shape, max_print_height };
    m_texture_filename = texture_filename;
    m_model_filename = model_filename;
    m_extended_bounding_box = this->calc_extended_bounding_box();

#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_contour = ExPolygon(Polygon::new_scale(bed_shape));
    m_polygon = offset(m_contour.contour, (float)m_contour.contour.bounding_box().radius() * 1.7f, jtRound, scale_(0.5)).front();

    m_triangles.reset();
    m_gridlines.reset();
    m_contourlines.reset();
#else
    ExPolygon poly{ Polygon::new_scale(bed_shape) };

    calc_triangles(poly);

    const BoundingBox& bed_bbox = poly.contour.bounding_box();
    calc_gridlines(poly, bed_bbox);

    calc_contourlines(poly);

    m_polygon = offset(poly.contour, (float)bed_bbox.radius() * 1.7f, jtRound, scale_(0.5)).front();

    this->release_VBOs();
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
    m_texture.reset();
    m_model.reset();

    // Set the origin and size for rendering the coordinate system axes.
    m_axes.set_origin({ 0.0, 0.0, static_cast<double>(GROUND_Z) });
    m_axes.set_stem_length(0.1f * static_cast<float>(m_build_volume.bounding_volume().max_size()));

    // Let the calee to update the UI.
    return true;
}

bool Bed3D::contains(const Point& point) const
{
    return m_polygon.contains(point);
}

Point Bed3D::point_projection(const Point& point) const
{
    return m_polygon.point_projection(point);
}

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render(GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix, bool bottom, float scale_factor, bool show_axes, bool show_texture)
{
    render_internal(canvas, view_matrix, projection_matrix, bottom, scale_factor, show_axes, show_texture, false);
}

void Bed3D::render_for_picking(GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix, bool bottom, float scale_factor)
{
    render_internal(canvas, view_matrix, projection_matrix, bottom, scale_factor, false, false, true);
}
#else
void Bed3D::render(GLCanvas3D& canvas, bool bottom, float scale_factor, bool show_axes, bool show_texture)
{
    render_internal(canvas, bottom, scale_factor, show_axes, show_texture, false);
}

void Bed3D::render_for_picking(GLCanvas3D& canvas, bool bottom, float scale_factor)
{
    render_internal(canvas, bottom, scale_factor, false, false, true);
}
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_internal(GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix, bool bottom, float scale_factor,
    bool show_axes, bool show_texture, bool picking)
#else
void Bed3D::render_internal(GLCanvas3D& canvas, bool bottom, float scale_factor,
    bool show_axes, bool show_texture, bool picking)
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
{
    m_scale_factor = scale_factor;

    if (show_axes)
        render_axes();

    glsafe(::glEnable(GL_DEPTH_TEST));

#if ENABLE_LEGACY_OPENGL_REMOVAL
    m_model.set_color(picking ? PICKING_MODEL_COLOR : DEFAULT_MODEL_COLOR);
#else
    m_model.set_color(-1, picking ? PICKING_MODEL_COLOR : DEFAULT_MODEL_COLOR);
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

    switch (m_type)
    {
#if ENABLE_GL_SHADERS_ATTRIBUTES
    case Type::System: { render_system(canvas, view_matrix, projection_matrix, bottom, show_texture); break; }
    default:
    case Type::Custom: { render_custom(canvas, view_matrix, projection_matrix, bottom, show_texture, picking); break; }
#else
    case Type::System: { render_system(canvas, bottom, show_texture); break; }
    default:
    case Type::Custom: { render_custom(canvas, bottom, show_texture, picking); break; }
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
    }

    glsafe(::glDisable(GL_DEPTH_TEST));
}

// Calculate an extended bounding box from axes and current model for visualization purposes.
BoundingBoxf3 Bed3D::calc_extended_bounding_box() const
{
    BoundingBoxf3 out { m_build_volume.bounding_volume() };
    const Vec3d size = out.size();
    // ensures that the bounding box is set as defined or the following calls to merge() will not work as intented
    if (size.x() > 0.0 && size.y() > 0.0 && !out.defined)
        out.defined = true;
    // Reset the build volume Z, we don't want to zoom to the top of the build volume if it is empty.
    out.min.z() = 0.0;
    out.max.z() = 0.0;
    // extend to contain axes
    out.merge(m_axes.get_origin() + m_axes.get_total_length() * Vec3d::Ones());
#if ENABLE_WORLD_COORDINATE
    out.merge(out.min + Vec3d(-m_axes.get_tip_radius(), -m_axes.get_tip_radius(), out.max.z()));
#else
    out.merge(out.min + Vec3d(-Axes::DefaultTipRadius, -Axes::DefaultTipRadius, out.max.z()));
#endif // ENABLE_WORLD_COORDINATE
    // extend to contain model, if any
    BoundingBoxf3 model_bb = m_model.get_bounding_box();
    if (model_bb.defined) {
        model_bb.translate(m_model_offset);
        out.merge(model_bb);
    }
    return out;
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void Bed3D::init_triangles()
{
    if (m_triangles.is_initialized())
        return;

    if (m_contour.empty())
        return;

    const std::vector<Vec2f> triangles = triangulate_expolygon_2f(m_contour, NORMALS_UP);
    if (triangles.empty() || triangles.size() % 3 != 0)
        return;

    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3T2 };
    init_data.reserve_vertices(triangles.size());
    init_data.reserve_indices(triangles.size() / 3);

    Vec2f min = triangles.front();
    Vec2f max = min;
    for (const Vec2f& v : triangles) {
        min = min.cwiseMin(v).eval();
        max = max.cwiseMax(v).eval();
    }

    const Vec2f size = max - min;
    if (size.x() <= 0.0f || size.y() <= 0.0f)
        return;

    Vec2f inv_size = size.cwiseInverse();
    inv_size.y() *= -1.0f;

    // vertices + indices
    unsigned int vertices_counter = 0;
    for (const Vec2f& v : triangles) {
        const Vec3f p = { v.x(), v.y(), GROUND_Z };
        init_data.add_vertex(p, (Vec2f)(v - min).cwiseProduct(inv_size).eval());
        ++vertices_counter;
        if (vertices_counter % 3 == 0)
            init_data.add_triangle(vertices_counter - 3, vertices_counter - 2, vertices_counter - 1);
    }

    m_triangles.init_from(std::move(init_data));
}

void Bed3D::init_gridlines()
{
    if (m_gridlines.is_initialized())
        return;

    if (m_contour.empty())
        return;

    const BoundingBox& bed_bbox = m_contour.contour.bounding_box();
    const coord_t step = scale_(10.0);

    Polylines axes_lines;
    for (coord_t x = bed_bbox.min.x(); x <= bed_bbox.max.x(); x += step) {
        Polyline line;
        line.append(Point(x, bed_bbox.min.y()));
        line.append(Point(x, bed_bbox.max.y()));
        axes_lines.push_back(line);
    }
    for (coord_t y = bed_bbox.min.y(); y <= bed_bbox.max.y(); y += step) {
        Polyline line;
        line.append(Point(bed_bbox.min.x(), y));
        line.append(Point(bed_bbox.max.x(), y));
        axes_lines.push_back(line);
    }

    // clip with a slightly grown expolygon because our lines lay on the contours and may get erroneously clipped
    Lines gridlines = to_lines(intersection_pl(axes_lines, offset(m_contour, float(SCALED_EPSILON))));

    // append bed contours
    Lines contour_lines = to_lines(m_contour);
    std::copy(contour_lines.begin(), contour_lines.end(), std::back_inserter(gridlines));

    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(2 * gridlines.size());
    init_data.reserve_indices(2 * gridlines.size());

    for (const Slic3r::Line& l : gridlines) {
        init_data.add_vertex(Vec3f(unscale<float>(l.a.x()), unscale<float>(l.a.y()), GROUND_Z));
        init_data.add_vertex(Vec3f(unscale<float>(l.b.x()), unscale<float>(l.b.y()), GROUND_Z));
        const unsigned int vertices_counter = (unsigned int)init_data.vertices_count();
        init_data.add_line(vertices_counter - 2, vertices_counter - 1);
    }

    m_gridlines.init_from(std::move(init_data));
}
#else
void Bed3D::calc_triangles(const ExPolygon& poly)
{
    if (! m_triangles.set_from_triangles(triangulate_expolygon_2f(poly, NORMALS_UP), GROUND_Z))
        BOOST_LOG_TRIVIAL(error) << "Unable to create bed triangles";
}

void Bed3D::calc_gridlines(const ExPolygon& poly, const BoundingBox& bed_bbox)
{
    Polylines axes_lines;
    for (coord_t x = bed_bbox.min.x(); x <= bed_bbox.max.x(); x += scale_(10.0)) {
        Polyline line;
        line.append(Point(x, bed_bbox.min.y()));
        line.append(Point(x, bed_bbox.max.y()));
        axes_lines.push_back(line);
    }
    for (coord_t y = bed_bbox.min.y(); y <= bed_bbox.max.y(); y += scale_(10.0)) {
        Polyline line;
        line.append(Point(bed_bbox.min.x(), y));
        line.append(Point(bed_bbox.max.x(), y));
        axes_lines.push_back(line);
    }

    // clip with a slightly grown expolygon because our lines lay on the contours and may get erroneously clipped
    Lines gridlines = to_lines(intersection_pl(axes_lines, offset(poly, (float)SCALED_EPSILON)));

    // append bed contours
    Lines contour_lines = to_lines(poly);
    std::copy(contour_lines.begin(), contour_lines.end(), std::back_inserter(gridlines));

    if (!m_gridlines.set_from_lines(gridlines, GROUND_Z))
        BOOST_LOG_TRIVIAL(error) << "Unable to create bed grid lines\n";
}
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if ENABLE_LEGACY_OPENGL_REMOVAL
void Bed3D::init_contourlines()
{
    if (m_contourlines.is_initialized())
        return;

    if (m_contour.empty())
        return;

    const Lines contour_lines = to_lines(m_contour);

    GLModel::Geometry init_data;
    init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
    init_data.reserve_vertices(2 * contour_lines.size());
    init_data.reserve_indices(2 * contour_lines.size());

    for (const Slic3r::Line& l : contour_lines) {
        init_data.add_vertex(Vec3f(unscale<float>(l.a.x()), unscale<float>(l.a.y()), GROUND_Z));
        init_data.add_vertex(Vec3f(unscale<float>(l.b.x()), unscale<float>(l.b.y()), GROUND_Z));
        const unsigned int vertices_counter = (unsigned int)init_data.vertices_count();
        init_data.add_line(vertices_counter - 2, vertices_counter - 1);
    }

    m_contourlines.init_from(std::move(init_data));
    m_contourlines.set_color({ 1.0f, 1.0f, 1.0f, 0.5f });
}
#else
void Bed3D::calc_contourlines(const ExPolygon& poly)
{
    const Lines contour_lines = to_lines(poly);
    if (!m_contourlines.set_from_lines(contour_lines, GROUND_Z))
        BOOST_LOG_TRIVIAL(error) << "Unable to create bed contour lines\n";
}
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

// Try to match the print bed shape with the shape of an active profile. If such a match exists,
// return the print bed model.
std::tuple<Bed3D::Type, std::string, std::string> Bed3D::detect_type(const Pointfs& shape)
{
    auto bundle = wxGetApp().preset_bundle;
    if (bundle != nullptr) {
        const Preset* curr = &bundle->printers.get_selected_preset();
        while (curr != nullptr) {
            if (curr->config.has("bed_shape")) {
                if (shape == dynamic_cast<const ConfigOptionPoints*>(curr->config.option("bed_shape"))->values) {
                    std::string model_filename = PresetUtils::system_printer_bed_model(*curr);
                    std::string texture_filename = PresetUtils::system_printer_bed_texture(*curr);
                    if (!model_filename.empty() && !texture_filename.empty())
                        return { Type::System, model_filename, texture_filename };
                }
            }

            curr = bundle->printers.get_preset_parent(*curr);
        }
    }

    return { Type::Custom, {}, {} };
}

void Bed3D::render_axes()
{
    if (m_build_volume.valid())
#if ENABLE_WORLD_COORDINATE
#if ENABLE_GL_SHADERS_ATTRIBUTES
        m_axes.render(Transform3d::Identity(), 0.25f);
#else
        m_axes.render(0.25f);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
#else
        m_axes.render();
#endif // ENABLE_WORLD_COORDINATE
}

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_system(GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix, bool bottom, bool show_texture)
{
    if (!bottom)
        render_model(view_matrix, projection_matrix);

    if (show_texture)
        render_texture(bottom, canvas, view_matrix, projection_matrix);
    else if (bottom)
        render_contour(view_matrix, projection_matrix);
}
#else
void Bed3D::render_system(GLCanvas3D& canvas, bool bottom, bool show_texture)
{
    if (!bottom)
        render_model();

    if (show_texture)
        render_texture(bottom, canvas);
    else if (bottom)
        render_contour();
}
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_texture(bool bottom, GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix)
#else
void Bed3D::render_texture(bool bottom, GLCanvas3D& canvas)
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
{
    if (m_texture_filename.empty()) {
        m_texture.reset();
#if ENABLE_GL_SHADERS_ATTRIBUTES
        render_default(bottom, false, true, view_matrix, projection_matrix);
#else
        render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        return;
    }

    if (m_texture.get_id() == 0 || m_texture.get_source() != m_texture_filename) {
        m_texture.reset();

        if (boost::algorithm::iends_with(m_texture_filename, ".svg")) {
            // use higher resolution images if graphic card and opengl version allow
            GLint max_tex_size = OpenGLManager::get_gl_info().get_max_tex_size();
            if (m_temp_texture.get_id() == 0 || m_temp_texture.get_source() != m_texture_filename) {
                // generate a temporary lower resolution texture to show while no main texture levels have been compressed
                if (!m_temp_texture.load_from_svg_file(m_texture_filename, false, false, false, max_tex_size / 8)) {
#if ENABLE_GL_SHADERS_ATTRIBUTES
                    render_default(bottom, false, true, view_matrix, projection_matrix);
#else
                    render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
                    return;
                }
                canvas.request_extra_frame();
            }

            // starts generating the main texture, compression will run asynchronously
            if (!m_texture.load_from_svg_file(m_texture_filename, true, true, true, max_tex_size)) {
#if ENABLE_GL_SHADERS_ATTRIBUTES
                render_default(bottom, false, true, view_matrix, projection_matrix);
#else
                render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
                return;
            }
        } 
        else if (boost::algorithm::iends_with(m_texture_filename, ".png")) {
            // generate a temporary lower resolution texture to show while no main texture levels have been compressed
            if (m_temp_texture.get_id() == 0 || m_temp_texture.get_source() != m_texture_filename) {
                if (!m_temp_texture.load_from_file(m_texture_filename, false, GLTexture::None, false)) {
#if ENABLE_GL_SHADERS_ATTRIBUTES
                    render_default(bottom, false, true, view_matrix, projection_matrix);
#else
                    render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
                    return;
                }
                canvas.request_extra_frame();
            }

            // starts generating the main texture, compression will run asynchronously
            if (!m_texture.load_from_file(m_texture_filename, true, GLTexture::MultiThreaded, true)) {
#if ENABLE_GL_SHADERS_ATTRIBUTES
                render_default(bottom, false, true, view_matrix, projection_matrix);
#else
                render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
                return;
            }
        }
        else {
#if ENABLE_GL_SHADERS_ATTRIBUTES
            render_default(bottom, false, true, view_matrix, projection_matrix);
#else
            render_default(bottom, false, true);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
            return;
        }
    }
    else if (m_texture.unsent_compressed_data_available()) {
        // sends to gpu the already available compressed levels of the main texture
        m_texture.send_compressed_data_to_gpu();

        // the temporary texture is not needed anymore, reset it
        if (m_temp_texture.get_id() != 0)
            m_temp_texture.reset();

        canvas.request_extra_frame();
    }

#if ENABLE_LEGACY_OPENGL_REMOVAL
    init_triangles();

    GLShaderProgram* shader = wxGetApp().get_shader("printbed");
    if (shader != nullptr) {
        shader->start_using();
#if ENABLE_GL_SHADERS_ATTRIBUTES
        shader->set_uniform("view_model_matrix", view_matrix);
        shader->set_uniform("projection_matrix", projection_matrix);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        shader->set_uniform("transparent_background", bottom);
        shader->set_uniform("svg_source", boost::algorithm::iends_with(m_texture.get_source(), ".svg"));

        glsafe(::glEnable(GL_DEPTH_TEST));
        if (bottom)
            glsafe(::glDepthMask(GL_FALSE));

        glsafe(::glEnable(GL_BLEND));
        glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

        if (bottom)
            glsafe(::glFrontFace(GL_CW));

        // show the temporary texture while no compressed data is available
        GLuint tex_id = (GLuint)m_temp_texture.get_id();
        if (tex_id == 0)
            tex_id = (GLuint)m_texture.get_id();

        glsafe(::glBindTexture(GL_TEXTURE_2D, tex_id));
        m_triangles.render();
        glsafe(::glBindTexture(GL_TEXTURE_2D, 0));

        if (bottom)
            glsafe(::glFrontFace(GL_CCW));

        glsafe(::glDisable(GL_BLEND));
        if (bottom)
            glsafe(::glDepthMask(GL_TRUE));

        shader->stop_using();
    }
#else
    if (m_triangles.get_vertices_count() > 0) {
        GLShaderProgram* shader = wxGetApp().get_shader("printbed");
        if (shader != nullptr) {
            shader->start_using();
            shader->set_uniform("transparent_background", bottom);
            shader->set_uniform("svg_source", boost::algorithm::iends_with(m_texture.get_source(), ".svg"));

            if (m_vbo_id == 0) {
                glsafe(::glGenBuffers(1, &m_vbo_id));
                glsafe(::glBindBuffer(GL_ARRAY_BUFFER, m_vbo_id));
                glsafe(::glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)m_triangles.get_vertices_data_size(), (const GLvoid*)m_triangles.get_vertices_data(), GL_STATIC_DRAW));
                glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
            }

            glsafe(::glEnable(GL_DEPTH_TEST));
            if (bottom)
                glsafe(::glDepthMask(GL_FALSE));

            glsafe(::glEnable(GL_BLEND));
            glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

            if (bottom)
                glsafe(::glFrontFace(GL_CW));

            const unsigned int stride = m_triangles.get_vertex_data_size();

            // show the temporary texture while no compressed data is available
            GLuint tex_id = (GLuint)m_temp_texture.get_id();
            if (tex_id == 0)
                tex_id = (GLuint)m_texture.get_id();

            glsafe(::glBindTexture(GL_TEXTURE_2D, tex_id));
            glsafe(::glBindBuffer(GL_ARRAY_BUFFER, m_vbo_id));

            glsafe(::glVertexPointer(3, GL_FLOAT, stride, (const void*)(intptr_t)m_triangles.get_position_offset()));
            glsafe(::glEnableClientState(GL_VERTEX_ARRAY));

            glsafe(::glTexCoordPointer(2, GL_FLOAT, stride, (const void*)(intptr_t)m_triangles.get_tex_coords_offset()));
            glsafe(::glEnableClientState(GL_TEXTURE_COORD_ARRAY));

            glsafe(::glDrawArrays(GL_TRIANGLES, 0, (GLsizei)m_triangles.get_vertices_count()));

            glsafe(::glDisableClientState(GL_TEXTURE_COORD_ARRAY));
            glsafe(::glDisableClientState(GL_VERTEX_ARRAY));

            glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
            glsafe(::glBindTexture(GL_TEXTURE_2D, 0));

            if (bottom)
                glsafe(::glFrontFace(GL_CCW));

            glsafe(::glDisable(GL_BLEND));
            if (bottom)
                glsafe(::glDepthMask(GL_TRUE));

            shader->stop_using();
        }
    }
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_model(const Transform3d& view_matrix, const Transform3d& projection_matrix)
#else
void Bed3D::render_model()
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
{
    if (m_model_filename.empty())
        return;

    if (m_model.get_filename() != m_model_filename && m_model.init_from_file(m_model_filename)) {
#if ENABLE_LEGACY_OPENGL_REMOVAL
        m_model.set_color(DEFAULT_MODEL_COLOR);
#else
        m_model.set_color(-1, DEFAULT_MODEL_COLOR);
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

        // move the model so that its origin (0.0, 0.0, 0.0) goes into the bed shape center and a bit down to avoid z-fighting with the texture quad
        m_model_offset = to_3d(m_build_volume.bounding_volume2d().center(), -0.03);

        // update extended bounding box
        m_extended_bounding_box = this->calc_extended_bounding_box();
    }

    if (!m_model.get_filename().empty()) {
        GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
        if (shader != nullptr) {
            shader->start_using();
            shader->set_uniform("emission_factor", 0.0f);
#if ENABLE_GL_SHADERS_ATTRIBUTES
            const Transform3d matrix = view_matrix * Geometry::assemble_transform(m_model_offset);
            shader->set_uniform("view_model_matrix", matrix);
            shader->set_uniform("projection_matrix", projection_matrix);
            shader->set_uniform("normal_matrix", (Matrix3d)matrix.matrix().block(0, 0, 3, 3).inverse().transpose());
#else
            glsafe(::glPushMatrix());
            glsafe(::glTranslated(m_model_offset.x(), m_model_offset.y(), m_model_offset.z()));
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
            m_model.render();
#if !ENABLE_GL_SHADERS_ATTRIBUTES
            glsafe(::glPopMatrix());
#endif // !ENABLE_GL_SHADERS_ATTRIBUTES
            shader->stop_using();
        }
    }
}

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_custom(GLCanvas3D& canvas, const Transform3d& view_matrix, const Transform3d& projection_matrix, bool bottom, bool show_texture, bool picking)
#else
void Bed3D::render_custom(GLCanvas3D& canvas, bool bottom, bool show_texture, bool picking)
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
{
    if (m_texture_filename.empty() && m_model_filename.empty()) {
#if ENABLE_GL_SHADERS_ATTRIBUTES
        render_default(bottom, picking, show_texture, view_matrix, projection_matrix);
#else
        render_default(bottom, picking, show_texture);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
        return;
    }

    if (!bottom)
#if ENABLE_GL_SHADERS_ATTRIBUTES
        render_model(view_matrix, projection_matrix);
#else
        render_model();
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

    if (show_texture)
#if ENABLE_GL_SHADERS_ATTRIBUTES
        render_texture(bottom, canvas, view_matrix, projection_matrix);
    else if (bottom)
        render_contour(view_matrix, projection_matrix);
#else
        render_texture(bottom, canvas);
    else if (bottom)
        render_contour();
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
}

#if ENABLE_GL_SHADERS_ATTRIBUTES
void Bed3D::render_default(bool bottom, bool picking, bool show_texture, const Transform3d& view_matrix, const Transform3d& projection_matrix)
#else
void Bed3D::render_default(bool bottom, bool picking, bool show_texture)
#endif // ENABLE_GL_SHADERS_ATTRIBUTES
{
    m_texture.reset();

#if ENABLE_LEGACY_OPENGL_REMOVAL
    init_gridlines();
    init_triangles();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();

#if ENABLE_GL_SHADERS_ATTRIBUTES
        shader->set_uniform("view_model_matrix", view_matrix);
        shader->set_uniform("projection_matrix", projection_matrix);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

        glsafe(::glEnable(GL_DEPTH_TEST));
        glsafe(::glEnable(GL_BLEND));
        glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

        const bool has_model = !m_model.get_filename().empty();

        if (!has_model && !bottom) {
            // draw background
            glsafe(::glDepthMask(GL_FALSE));
            m_triangles.set_color(picking ? PICKING_MODEL_COLOR : DEFAULT_MODEL_COLOR);
            m_triangles.render();
            glsafe(::glDepthMask(GL_TRUE));
        }

        if (!picking && show_texture) {
            // draw grid
            glsafe(::glLineWidth(1.5f * m_scale_factor));
            m_gridlines.set_color(has_model && !bottom ? DEFAULT_SOLID_GRID_COLOR : DEFAULT_TRANSPARENT_GRID_COLOR);
            m_gridlines.render();
        }
        else if (!show_texture)
            render_contour(view_matrix, projection_matrix);

        glsafe(::glDisable(GL_BLEND));

        shader->stop_using();
    }
#else
    const unsigned int triangles_vcount = m_triangles.get_vertices_count();
    if (triangles_vcount > 0) {
        const bool has_model = !m_model.get_filename().empty();

        glsafe(::glEnable(GL_DEPTH_TEST));
        glsafe(::glEnable(GL_BLEND));
        glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

        glsafe(::glEnableClientState(GL_VERTEX_ARRAY));

        if (!has_model && !bottom) {
            // draw background
            glsafe(::glDepthMask(GL_FALSE));
            glsafe(::glColor4fv(picking ? PICKING_MODEL_COLOR.data() : DEFAULT_MODEL_COLOR.data()));
            glsafe(::glNormal3d(0.0f, 0.0f, 1.0f));
            glsafe(::glVertexPointer(3, GL_FLOAT, m_triangles.get_vertex_data_size(), (GLvoid*)m_triangles.get_vertices_data()));
            glsafe(::glDrawArrays(GL_TRIANGLES, 0, (GLsizei)triangles_vcount));
            glsafe(::glDepthMask(GL_TRUE));
        }

        if (!picking && show_texture) {
            // draw grid
            glsafe(::glLineWidth(1.5f * m_scale_factor));
            glsafe(::glColor4fv(has_model && !bottom ? DEFAULT_SOLID_GRID_COLOR.data() : DEFAULT_TRANSPARENT_GRID_COLOR.data()));
            glsafe(::glVertexPointer(3, GL_FLOAT, m_gridlines.get_vertex_data_size(), (GLvoid*)m_gridlines.get_vertices_data()));
            glsafe(::glDrawArrays(GL_LINES, 0, (GLsizei)m_gridlines.get_vertices_count()));
        }
        else if (!show_texture)
            render_contour();

        glsafe(::glDisableClientState(GL_VERTEX_ARRAY));

        glsafe(::glDisable(GL_BLEND));
    }
#endif // ENABLE_LEGACY_OPENGL_REMOVAL
}

#if ENABLE_LEGACY_OPENGL_REMOVAL
void Bed3D::render_contour(const Transform3d& view_matrix, const Transform3d& projection_matrix)
{
    init_contourlines();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader != nullptr) {
        shader->start_using();

#if ENABLE_GL_SHADERS_ATTRIBUTES
        shader->set_uniform("view_model_matrix", view_matrix);
        shader->set_uniform("projection_matrix", projection_matrix);
#endif // ENABLE_GL_SHADERS_ATTRIBUTES

        glsafe(::glEnable(GL_DEPTH_TEST));
        glsafe(::glEnable(GL_BLEND));
        glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

        // draw contour
        glsafe(::glLineWidth(1.5f * m_scale_factor));
        m_contourlines.render();

        glsafe(::glDisable(GL_BLEND));

        shader->stop_using();
    }
}
#else
void Bed3D::render_contour()
{
    glsafe(::glEnableClientState(GL_VERTEX_ARRAY));
    glsafe(::glLineWidth(1.5f * m_scale_factor));
    glsafe(::glColor4f(1.0f, 1.0f, 1.0f, 0.5f));
    glsafe(::glVertexPointer(3, GL_FLOAT, m_contourlines.get_vertex_data_size(), (GLvoid*)m_contourlines.get_vertices_data()));
    glsafe(::glDrawArrays(GL_LINES, 0, (GLsizei)m_contourlines.get_vertices_count()));
    glsafe(::glDisableClientState(GL_VERTEX_ARRAY));
}
#endif // ENABLE_LEGACY_OPENGL_REMOVAL

#if !ENABLE_LEGACY_OPENGL_REMOVAL
void Bed3D::release_VBOs()
{
    if (m_vbo_id > 0) {
        glsafe(::glDeleteBuffers(1, &m_vbo_id));
        m_vbo_id = 0;
    }
}
#endif // !ENABLE_LEGACY_OPENGL_REMOVAL

} // GUI
} // Slic3r
