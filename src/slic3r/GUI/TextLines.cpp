#include "TextLines.hpp"

#include <GL/glew.h>

#include "libslic3r/Model.hpp"

#include "libslic3r/Emboss.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/Tesselate.hpp"

#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/ExPolygonsIndex.hpp"

#include "slic3r/GUI/Selection.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GLModel.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/3DScene.hpp"

using namespace Slic3r;
using namespace Slic3r::GUI;

namespace {

const Slic3r::Polygon *largest(const Slic3r::Polygons &polygons)
{
    if (polygons.empty())
        return nullptr;
    if (polygons.size() == 1)
        return &polygons.front();

    // compare polygon to find largest
    size_t                 biggest_size = 0;
    const Slic3r::Polygon *result       = nullptr;
    for (const Slic3r::Polygon &polygon : polygons) {
        Point  s    = polygon.bounding_box().size();
        size_t size = s.x() * s.y();
        if (size <= biggest_size)
            continue;
        biggest_size = size;
        result       = &polygon;
    }
    return result;
}

GLModel create_model(const Slic3r::Polygon &polygon, float width_half = 0.5f, ColorRGBA color = ColorRGBA(0.f, 1.f, .2f, 0.5f))
{
    // Improve: Create torus instead of flat path (with model overlaps)

    assert(!polygon.empty());
    if (polygon.empty())
        return {};

    // add a small positive offset to avoid z-fighting
    float                  offset               = static_cast<float>(scale_(0.015f));
    Polygons               polygons_expanded    = expand(polygon, offset);
    const Slic3r::Polygon *polygon_expanded_ptr = largest(polygons_expanded);
    assert(polygon_expanded_ptr != nullptr);
    if (polygon_expanded_ptr == nullptr || polygon_expanded_ptr->empty())
        return {};
    const Slic3r::Polygon &polygon_expanded = *polygon_expanded_ptr;

    // inspired by 3DScene.cpp void GLVolume::SinkingContours::update()
    GLModel::Geometry init_data;
    init_data.format = {GLModel::Geometry::EPrimitiveType::Triangles, GUI::GLModel::Geometry::EVertexLayout::P3};
    init_data.color  = color;

    size_t count = polygon_expanded.size();
    init_data.reserve_vertices(2 * count);
    init_data.reserve_indices(2 * count);

    for (const Point &point : polygon_expanded.points) {
        Vec2f point_d = unscale(point).cast<float>();
        Vec3f vertex(point_d.x(), point_d.y(), width_half);
        init_data.add_vertex(vertex);
        vertex.z() *= -1;
        init_data.add_vertex(vertex);
    }

    unsigned int prev_i = count - 1;
    for (unsigned int i = 0; i < count; ++i) {
        // t .. top
        // b .. bottom
        unsigned int t1 = prev_i * 2;
        unsigned int b1 = t1 + 1;
        unsigned int t2 = i * 2;
        unsigned int b2 = t2 + 1;
        init_data.add_triangle(t1, b1, t2);
        init_data.add_triangle(b2, t2, b1);
        prev_i = i;
    }

    GLModel gl_model;
    gl_model.init_from(std::move(init_data));
    return gl_model;
}

GLModel create_model(const Polygons &polygons, float width_half = 0.5f, ColorRGBA color = ColorRGBA(0.f, 1.f, .2f, 0.5f))
{
    assert(!polygons.empty());

    // add a small positive offset to avoid z-fighting
    float    offset            = static_cast<float>(scale_(0.015f));
    Polygons polygons_expanded = expand(polygons, offset);

    // inspired by 3DScene.cpp void GLVolume::SinkingContours::update()
    GLModel::Geometry init_data;
    init_data.format = {GLModel::Geometry::EPrimitiveType::Triangles, GUI::GLModel::Geometry::EVertexLayout::P3};
    init_data.color  = color;

    size_t count = count_points(polygons);
    init_data.reserve_vertices(2 * count);
    init_data.reserve_indices(2 * count);

    unsigned int vertices_counter = 0;
    for (const Slic3r::Polygon &polygon : polygons_expanded) {
        for (const Point &point : polygon.points) {
            Vec2f point_d = unscale(point).cast<float>();
            Vec3f vertex(point_d.x(), point_d.y(), width_half);
            init_data.add_vertex(vertex);
            vertex.z() *= -1;
            init_data.add_vertex(vertex);
        }

        auto         points_count = static_cast<unsigned int>(polygon.points.size());
        unsigned int prev_i       = points_count - 1;
        for (unsigned int i = 0; i < points_count; i++) {
            // t .. top
            // b .. bottom
            unsigned int t1 = vertices_counter + prev_i * 2;
            unsigned int b1 = t1 + 1;
            unsigned int t2 = vertices_counter + i * 2;
            unsigned int b2 = t2 + 1;
            init_data.add_triangle(t1, b1, t2);
            init_data.add_triangle(b2, t2, b1);
            prev_i = i;
        }
        vertices_counter += 2 * points_count;
    }

    GLModel gl_model;
    gl_model.init_from(std::move(init_data));
    return gl_model;
}

} // namespace

void TextLinesModel::init(const Selection &selection)
{
    const GLVolume *gl_volume_ptr = selection.get_first_volume();
    if (gl_volume_ptr == nullptr)
        return;
    const GLVolume        &gl_volume = *gl_volume_ptr;
    const ModelObjectPtrs &objects   = selection.get_model()->objects;
    const ModelObject     *mo_ptr    = get_model_object(gl_volume, objects);
    if (mo_ptr == nullptr)
        return;
    const ModelObject &mo = *mo_ptr;

    const ModelVolume *mv_ptr = get_model_volume(gl_volume, objects);
    if (mv_ptr == nullptr)
        return;
    const ModelVolume &mv = *mv_ptr;

    const std::optional<TextConfiguration> tc_opt = mv.text_configuration;
    if (!tc_opt.has_value())
        return;

    unsigned count_lines = Emboss::get_count_lines(tc_opt->text);
    if (count_lines == 0)
        return;
    
    // TODO: Calc correct line height by line gap + font file info
    // Be carefull it is not correct !!!
    double line_height = tc_opt->style.prop.size_in_mm;
    double first_line_center = -(count_lines / 2) * line_height - ((count_lines % 2 == 0)? line_height/2. : 0.);
    std::vector<float> line_centers(count_lines);
    for (size_t i = 0; i < count_lines; ++i)
        line_centers[i] = first_line_center + i * line_height;

    const Transform3d &mv_trafo = gl_volume.get_volume_transformation().get_matrix();

    // contour transformation
    auto rot = Eigen::AngleAxis(M_PI_2, Vec3d::UnitX());
    Transform3d c_trafo = mv_trafo * rot;
    Transform3d c_trafo_inv = c_trafo.inverse();

    std::vector<Polygons> line_contours(count_lines);
    for (const ModelVolume *volume : mo.volumes) {
        // only part could be surface for volumes
        if (!volume->is_model_part())
            continue;

        // is selected volume
        if (mv.id() == volume->id())
            continue;

        MeshSlicingParams slicing_params;
        slicing_params.trafo = c_trafo_inv * volume->get_matrix();

        for (size_t i = 0; i < count_lines; ++i) {
            const Polygons polys = Slic3r::slice_mesh(volume->mesh().its, line_centers[i], slicing_params);
            if (polys.empty())
                continue;
            Polygons &contours = line_contours[i];
            contours.insert(contours.end(), polys.begin(), polys.end());
        }
    }

    lines.reserve(count_lines);
    lines.clear();
    // select closest contour
    Vec2d zero(0., 0.);
    for (const Polygons &polygons : line_contours) {
        // Improve: use int values and polygons only
        // Slic3r::Polygons polygons = union_(polygons);
        // std::vector<Slic3r::Line> lines = to_lines(polygons);
        // AABBTreeIndirect::Tree<2, Point> tree;
        // size_t line_idx;
        // Point hit_point;
        // Point::Scalar distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, point, line_idx, hit_point);

        ExPolygons expolygons = union_ex(polygons);
        std::vector<Linef> linesf = to_linesf(expolygons);
        AABBTreeIndirect::Tree2d tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(linesf);

        size_t line_idx;
        Vec2d  hit_point;
        double distance = AABBTreeLines::squared_distance_to_indexed_lines(linesf, tree, zero, line_idx, hit_point);

        // conversion between index of point and expolygon
        ExPolygonsIndices cvt(expolygons);
        ExPolygonsIndex index = cvt.cvt(static_cast<uint32_t>(line_idx));

        const Polygon& polygon = index.is_contour() ? expolygons[index.expolygons_index].contour :
                                    expolygons[index.expolygons_index].holes[index.hole_index()];

        Point hit_point_int = hit_point.cast<Point::coord_type>();
        TextLine tl{polygon, index.point_index, hit_point_int};
        lines.emplace_back(tl);
    }

    ColorRGBA color(.7f, .7f, .7f, .7f); // Gray

    // TODO: create model from all lines
    model = create_model(lines.front().polygon, 0.7f, color);
}

void TextLinesModel::render(const Transform3d &text_world)
{
    if (!model.is_initialized())
        return;

    GUI_App &app = wxGetApp();
    const GLShaderProgram *shader = app.get_shader("flat");
    if (shader == nullptr)
        return;

    const Camera &camera = app.plater()->get_camera();
    auto          rot    = Eigen::AngleAxis(M_PI_2, Vec3d::UnitX());

    shader->start_using();
    shader->set_uniform("view_model_matrix", camera.get_view_matrix() * text_world * rot);
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());

    bool is_depth_test = glIsEnabled(GL_DEPTH_TEST);
    if (!is_depth_test)
        glsafe(::glEnable(GL_DEPTH_TEST));

    bool is_blend = glIsEnabled(GL_BLEND);
    if (!is_blend)
        glsafe(::glEnable(GL_BLEND));
    // glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

    model.render();

    if (!is_depth_test)
        glsafe(::glDisable(GL_DEPTH_TEST));
    if (!is_blend)
        glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}
