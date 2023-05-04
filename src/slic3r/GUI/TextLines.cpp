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
using namespace Slic3r::Emboss;
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

indexed_triangle_set its_create_belt(const Slic3r::Polygon &polygon, float width_half) {
    // Improve: Create torus instead of flat belt path (with model overlaps)
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
    indexed_triangle_set model;
    size_t count = polygon_expanded.size();
    model.vertices.reserve(2 * count);
    model.indices.reserve(2 * count);

    for (const Point &point : polygon_expanded.points) {
        Vec2f point_d = unscale(point).cast<float>();
        Vec3f vertex(point_d.x(), point_d.y(), width_half);
        model.vertices.push_back(vertex);
        vertex.z() *= -1;
        model.vertices.push_back(vertex);
    }

    unsigned int prev_i = count - 1;
    for (unsigned int i = 0; i < count; ++i) {
        // t .. top
        // b .. bottom
        unsigned int t1 = prev_i * 2;
        unsigned int b1 = t1 + 1;
        unsigned int t2 = i * 2;
        unsigned int b2 = t2 + 1;
        model.indices.emplace_back(t1, b1, t2);
        model.indices.emplace_back(b2, t2, b1);
        prev_i = i;
    }
    return model;
}

indexed_triangle_set its_create_torus(const Slic3r::Polygon &polygon, float radius, size_t steps = 20)
{
    assert(!polygon.empty());
    if (polygon.empty())
        return {};

    size_t count = polygon.size();
    if (count < 3)
        return {};

    // convert and scale to float
    std::vector<Vec2f> points_d;
    points_d.reserve(count);
    for (const Point &point : polygon.points)
        points_d.push_back(unscale(point).cast<float>());

    // pre calculate normalized line directions
    auto calc_line_norm = [](const Vec2f &f, const Vec2f &s) -> Vec2f { return  (s - f).normalized(); };    
    std::vector<Vec2f> line_norm(points_d.size());
    for (size_t i = 0; i < count - 1; ++i)
        line_norm[i] = calc_line_norm(points_d[i], points_d[i + 1]);
    line_norm.back() = calc_line_norm(points_d.back(), points_d.front());

    // calculate normals for each point
    auto calc_norm = [](const Vec2f &prev, const Vec2f &next) -> Vec2f {
        Vec2f dir = prev + next;
        return Vec2f(-dir.x(), dir.y());
    };
    std::vector<Vec2f> points_norm(points_d.size());
    points_norm.front() = calc_norm(line_norm.back(), line_norm[1]);
    for (size_t i = 1; i < points_d.size() - 1; ++i)
        points_norm[i] = calc_norm(line_norm[i - 1], line_norm[i + 1]);        
    points_norm.back() = calc_norm(line_norm[points_d.size() - 2], line_norm.front());
    
    // precalculate sinus and cosinus
    double angle_step = 2 * M_PI / steps;
    std::vector<std::pair<double, float>> sin_cos;
    sin_cos.reserve(steps);
    for (size_t s = 0; s < steps; ++s) {
        double angle = s * angle_step;
        sin_cos.emplace_back(
            radius * std::sin(angle), 
            static_cast<float>(radius * std::cos(angle))
        );
    }
    
    // create torus model along polygon path
    indexed_triangle_set model;
    model.vertices.reserve(steps * count);
    model.indices.reserve(2 * steps * count);
    for (size_t i = 0; i < count; ++i) {
        const Vec2f point_d = points_d[i];
        const Vec2f norm = points_norm[i];
        for (const auto &[s, c] : sin_cos) {
            Vec2f xy = s * norm + point_d;
            model.vertices.emplace_back(xy.x(), xy.y(), c);
        }
    }

    unsigned int prev_i = count - 1;
    for (unsigned int i = 0; i < count; ++i) {        
        // TODO: solve <180, =180 and >180 angle
        // to not create self intersection

        // t .. top
        // b .. bottom
        unsigned int prev_t = (prev_i+1) * steps - 1;
        unsigned int t = (i+1) * steps - 1;
        for (size_t s = 0; s < steps; ++s) {
            unsigned int prev_b = prev_i * steps + s;
            unsigned int b = i * steps + s;
            model.indices.emplace_back(prev_t, prev_b, t);
            model.indices.emplace_back(b, t, prev_b);
            prev_t = prev_b;
            t = b;
        }
        prev_i = i;
    }
    return model;
}

// select closest contour for each line
TextLines select_closest_contour(const std::vector<Polygons> &line_contours) {
    TextLines result;
    result.reserve(line_contours.size());
    Vec2d zero(0., 0.);
    for (const Polygons &polygons : line_contours){
        if (polygons.empty()) {
            result.emplace_back();
            continue;
        }
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

        const Slic3r::Polygon& polygon = index.is_contour() ?
            expolygons[index.expolygons_index].contour :
            expolygons[index.expolygons_index].holes[index.hole_index()];

        Point hit_point_int = hit_point.cast<Point::coord_type>();
        TextLine tl{polygon, PolygonPoint{index.point_index, hit_point_int}};
        result.emplace_back(tl);
    }
    return result;
}

inline Eigen::AngleAxis<double> get_rotation() { return Eigen::AngleAxis(-M_PI_2, Vec3d::UnitX()); }

indexed_triangle_set create_its(const TextLines &lines) 
{ 
    const float model_half_width = 0.75; // [in volume mm]
    indexed_triangle_set its;
    // create model from polygons
    for (const TextLine &line : lines) {
        const Slic3r::Polygon &polygon = line.polygon;
        if (polygon.empty()) continue;
        indexed_triangle_set line_its = its_create_belt(polygon, model_half_width); 
        //indexed_triangle_set line_its = its_create_torus(polygon, model_half_width);
        auto transl = Eigen::Translation3d(0., line.y, 0.);
        Transform3d tr = transl * get_rotation();
        its_transform(line_its, tr);
        its_merge(its, line_its);
    }
    return its;
}

GLModel::Geometry create_geometry(const TextLines &lines)
{
    indexed_triangle_set its = create_its(lines);

    GLModel::Geometry geometry;
    geometry.format = {GLModel::Geometry::EPrimitiveType::Triangles, GUI::GLModel::Geometry::EVertexLayout::P3};
    ColorRGBA color(.7f, .7f, .7f, .7f); // Transparent Gray
    geometry.color = color;

    geometry.reserve_vertices(its.vertices.size());
    for (Vec3f vertex : its.vertices)
        geometry.add_vertex(vertex);

    geometry.reserve_indices(its.indices.size() * 3);
    for (Vec3i t : its.indices)
        geometry.add_triangle(t[0], t[1], t[2]);
    return geometry;    
}
} // namespace

void TextLinesModel::init(const Selection &selection, double line_height)
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
    
    double first_line_center = offset + (count_lines / 2) * line_height - ((count_lines % 2 == 0) ? line_height / 2. : 0.);
    std::vector<float> line_centers(count_lines);
    for (size_t i = 0; i < count_lines; ++i)
        line_centers[i] = static_cast<float>(first_line_center - i * line_height);

    const Transform3d &mv_trafo = gl_volume.get_volume_transformation().get_matrix();

    // contour transformation
    Transform3d c_trafo = mv_trafo * get_rotation();
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

    m_lines = select_closest_contour(line_contours);
    assert(m_lines.size() == count_lines);
    assert(line_centers.size() == count_lines);
    for (size_t i = 0; i < count_lines; ++i)
        m_lines[i].y = line_centers[i];

    m_model.reset();
    //*
    m_model.init_from(create_geometry(m_lines));
    /*/
    ColorRGBA color(.7f, .7f, .7f, .7f); // Transparent Gray
    m_model.set_color(color);
    m_model.init_from(create_its(m_lines));
    //*/
}

void TextLinesModel::render(const Transform3d &text_world)
{
    if (!m_model.is_initialized())
        return;

    GUI_App &app = wxGetApp();
    const GLShaderProgram *shader = app.get_shader("flat");
    if (shader == nullptr)
        return;

    const Camera &camera = app.plater()->get_camera();

    shader->start_using();
    shader->set_uniform("view_model_matrix", camera.get_view_matrix() * text_world);
    shader->set_uniform("projection_matrix", camera.get_projection_matrix());

    bool is_depth_test = glIsEnabled(GL_DEPTH_TEST);
    if (!is_depth_test)
        glsafe(::glEnable(GL_DEPTH_TEST));

    bool is_blend = glIsEnabled(GL_BLEND);
    if (!is_blend)
        glsafe(::glEnable(GL_BLEND));
    // glsafe(::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));

    m_model.render();

    if (!is_depth_test)
        glsafe(::glDisable(GL_DEPTH_TEST));
    if (!is_blend)
        glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}

double TextLinesModel::calc_line_height(const Slic3r::Emboss::FontFile &ff, const FontProp &fp)
{
    int line_height = Emboss::get_line_height(ff, fp); // In shape size
    double scale = Emboss::get_shape_scale(fp, ff);
    return line_height * scale;
}
