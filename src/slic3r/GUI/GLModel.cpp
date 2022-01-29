#include "libslic3r/libslic3r.h"
#include "GLModel.hpp"

#include "3DScene.hpp"
#include "GUI_App.hpp"
#include "GLShader.hpp"

#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Polygon.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <GL/glew.h>

namespace Slic3r {
namespace GUI {

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::Geometry::add_vertex(const Vec2f& position)
{
    assert(format.vertex_layout == EVertexLayout::P2);
    vertices.emplace_back(position.x());
    vertices.emplace_back(position.y());
}

void GLModel::Geometry::add_vertex(const Vec2f& position, const Vec2f& tex_coord)
{
    assert(format.vertex_layout == EVertexLayout::P2T2);
    vertices.emplace_back(position.x());
    vertices.emplace_back(position.y());
    vertices.emplace_back(tex_coord.x());
    vertices.emplace_back(tex_coord.y());
}

void GLModel::Geometry::add_vertex(const Vec3f& position)
{
    assert(format.vertex_layout == EVertexLayout::P3);
    vertices.emplace_back(position.x());
    vertices.emplace_back(position.y());
    vertices.emplace_back(position.z());
}

void GLModel::Geometry::add_vertex(const Vec3f& position, const Vec3f& normal)
{
    assert(format.vertex_layout == EVertexLayout::P3N3);
    vertices.emplace_back(position.x());
    vertices.emplace_back(position.y());
    vertices.emplace_back(position.z());
    vertices.emplace_back(normal.x());
    vertices.emplace_back(normal.y());
    vertices.emplace_back(normal.z());
}

void GLModel::Geometry::add_ushort_index(unsigned short id)
{
    if (format.index_type != EIndexType::USHORT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned short), &id, sizeof(unsigned short));
}

void GLModel::Geometry::add_uint_index(unsigned int id)
{
    if (format.index_type != EIndexType::UINT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned int), &id, sizeof(unsigned int));
}

void GLModel::Geometry::add_ushort_line(unsigned short id1, unsigned short id2)
{
    if (format.index_type != EIndexType::USHORT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + 2 * sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - 2 * sizeof(unsigned short), &id1, sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned short), &id2, sizeof(unsigned short));
}

void GLModel::Geometry::add_uint_line(unsigned int id1, unsigned int id2)
{
    if (format.index_type != EIndexType::UINT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + 2 * sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - 2 * sizeof(unsigned int), &id1, sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned int), &id2, sizeof(unsigned int));
}

void GLModel::Geometry::add_ushort_triangle(unsigned short id1, unsigned short id2, unsigned short id3)
{
    if (format.index_type != EIndexType::USHORT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + 3 * sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - 3 * sizeof(unsigned short), &id1, sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - 2 * sizeof(unsigned short), &id2, sizeof(unsigned short));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned short), &id3, sizeof(unsigned short));
}

void GLModel::Geometry::add_uint_triangle(unsigned int id1, unsigned int id2, unsigned int id3)
{
    if (format.index_type != EIndexType::UINT) {
        assert(false);
        return;
    }
    indices.resize(indices.size() + 3 * sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - 3 * sizeof(unsigned int), &id1, sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - 2 * sizeof(unsigned int), &id2, sizeof(unsigned int));
    ::memcpy(indices.data() + indices.size() - sizeof(unsigned int), &id3, sizeof(unsigned int));
}

Vec2f GLModel::Geometry::extract_position_2(size_t id) const
{
    const size_t p_stride = position_stride_floats(format);
    if (p_stride != 2) {
        assert(false);
        return { FLT_MAX, FLT_MAX };
    }

    if (vertices_count() <= id) {
        assert(false);
        return { FLT_MAX, FLT_MAX };
    }

    const float* start = &vertices[id * vertex_stride_floats(format) + position_offset_floats(format)];
    return { *(start + 0), *(start + 1) };
}

Vec3f GLModel::Geometry::extract_position_3(size_t id) const
{
    const size_t p_stride = position_stride_floats(format);
    if (p_stride != 3) {
        assert(false);
        return { FLT_MAX, FLT_MAX, FLT_MAX };
    }

    if (vertices_count() <= id) {
        assert(false);
        return { FLT_MAX, FLT_MAX, FLT_MAX };
    }

    const float* start = &vertices[id * vertex_stride_floats(format) + position_offset_floats(format)];
    return { *(start + 0), *(start + 1), *(start + 2) };
}

Vec3f GLModel::Geometry::extract_normal_3(size_t id) const
{
    const size_t n_stride = normal_stride_floats(format);
    if (n_stride != 3) {
        assert(false);
        return { FLT_MAX, FLT_MAX, FLT_MAX };
    }

    if (vertices_count() <= id) {
        assert(false);
        return { FLT_MAX, FLT_MAX, FLT_MAX };
    }

    const float* start = &vertices[id * vertex_stride_floats(format) + normal_offset_floats(format)];
    return { *(start + 0), *(start + 1), *(start + 2) };
}

Vec2f GLModel::Geometry::extract_tex_coord_2(size_t id) const
{
    const size_t t_stride = tex_coord_stride_floats(format);
    if (t_stride != 2) {
        assert(false);
        return { FLT_MAX, FLT_MAX };
    }

    if (vertices_count() <= id) {
        assert(false);
        return { FLT_MAX, FLT_MAX };
    }

    const float* start = &vertices[id * vertex_stride_floats(format) + tex_coord_offset_floats(format)];
    return { *(start + 0), *(start + 1) };
}

unsigned int GLModel::Geometry::extract_uint_index(size_t id) const
{
    if (format.index_type != EIndexType::UINT) {
        assert(false);
        return -1;
    }

    if (indices_count() <= id) {
        assert(false);
        return -1;
    }

    unsigned int ret = -1;
    ::memcpy(&ret, indices.data() + id * index_stride_bytes(format), sizeof(unsigned int));
    return ret;
}

unsigned short GLModel::Geometry::extract_ushort_index(size_t id) const
{
    if (format.index_type != EIndexType::USHORT) {
        assert(false);
        return -1;
    }

    if (indices_count() <= id) {
        assert(false);
        return -1;
    }

    unsigned short ret = -1;
    ::memcpy(&ret, indices.data() + id * index_stride_bytes(format), sizeof(unsigned short));
    return ret;
}

size_t GLModel::Geometry::vertex_stride_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2:   { return 2; }
    case EVertexLayout::P2T2: { return 4; }
    case EVertexLayout::P3:   { return 3; }
    case EVertexLayout::P3N3: { return 6; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::position_stride_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2:
    case EVertexLayout::P2T2: { return 2; }
    case EVertexLayout::P3:
    case EVertexLayout::P3N3: { return 3; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::position_offset_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2:
    case EVertexLayout::P2T2:
    case EVertexLayout::P3:
    case EVertexLayout::P3N3: { return 0; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::normal_stride_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P3N3: { return 3; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::normal_offset_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P3N3: { return 3; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::tex_coord_stride_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2T2: { return 2; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::tex_coord_offset_floats(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2T2: { return 2; }
    default:                  { assert(false); return 0; }
    };
}

size_t GLModel::Geometry::index_stride_bytes(const Format& format)
{
    switch (format.index_type)
    {
    case EIndexType::UINT:   { return sizeof(unsigned int); }
    case EIndexType::USHORT: { return sizeof(unsigned short); }
    default:                 { assert(false); return 0; }
    };
}

bool GLModel::Geometry::has_position(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2:
    case EVertexLayout::P2T2:
    case EVertexLayout::P3:
    case EVertexLayout::P3N3: { return true; }
    default:                  { assert(false); return false; }
    };
}

bool GLModel::Geometry::has_normal(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2:
    case EVertexLayout::P2T2:
    case EVertexLayout::P3:   { return false; }
    case EVertexLayout::P3N3: { return true; }
    default:                  { assert(false); return false; }
    };
}

bool GLModel::Geometry::has_tex_coord(const Format& format)
{
    switch (format.vertex_layout)
    {
    case EVertexLayout::P2T2: { return true; }
    case EVertexLayout::P2:
    case EVertexLayout::P3:
    case EVertexLayout::P3N3: { return false; }
    default:                  { assert(false); return false; }
    };
}
#else
size_t GLModel::Geometry::vertices_count() const
{
    size_t ret = 0;
    for (const Entity& entity : entities) {
        ret += entity.positions.size();
    }
    return ret;
}

size_t GLModel::Geometry::indices_count() const
{
    size_t ret = 0;
    for (const Entity& entity : entities) {
        ret += entity.indices.size();
    }
    return ret;
}
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::init_from(Geometry&& data)
#else
void GLModel::init_from(const Geometry& data)
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (is_initialized()) {
        // call reset() if you want to reuse this model
        assert(false);
        return;
    }

    if (data.vertices.empty() || data.indices.empty()) {
        assert(false);
        return;
    }

    m_render_data.geometry = std::move(data);

    // update bounding box
    for (size_t i = 0; i < vertices_count(); ++i) {
        const size_t position_stride = Geometry::position_stride_floats(data.format);
        if (position_stride == 3)
            m_bounding_box.merge(m_render_data.geometry.extract_position_3(i).cast<double>());
        else if (position_stride == 2) {
            const Vec2f position = m_render_data.geometry.extract_position_2(i);
            m_bounding_box.merge(Vec3f(position.x(), position.y(), 0.0f).cast<double>());
        }
    }
#else
    if (!m_render_data.empty()) // call reset() if you want to reuse this model
        return;

    for (const Geometry::Entity& entity : data.entities) {
        if (entity.positions.empty() || entity.indices.empty())
            continue;

        assert(entity.normals.empty() || entity.normals.size() == entity.positions.size());

        RenderData rdata;
        rdata.type = entity.type;
        rdata.color = entity.color;

        // vertices/normals data
        std::vector<float> vertices(6 * entity.positions.size());
        for (size_t i = 0; i < entity.positions.size(); ++i) {
            const size_t offset = i * 6;
            ::memcpy(static_cast<void*>(&vertices[offset]), static_cast<const void*>(entity.positions[i].data()), 3 * sizeof(float));
            if (!entity.normals.empty())
                ::memcpy(static_cast<void*>(&vertices[3 + offset]), static_cast<const void*>(entity.normals[i].data()), 3 * sizeof(float));
        }

        // indices data
        std::vector<unsigned int> indices = entity.indices;

        rdata.indices_count = static_cast<unsigned int>(indices.size());

        // update bounding box
        for (size_t i = 0; i < entity.positions.size(); ++i) {
            m_bounding_box.merge(entity.positions[i].cast<double>());
        }

        send_to_gpu(rdata, vertices, indices);
        m_render_data.emplace_back(rdata);
    }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::init_from(const TriangleMesh& mesh)
{
    init_from(mesh.its);
}

void GLModel::init_from(const indexed_triangle_set& its)
#else
void GLModel::init_from(const indexed_triangle_set& its, const BoundingBoxf3 &bbox)
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (is_initialized()) {
        // call reset() if you want to reuse this model
        assert(false);
        return;
    }

    if (its.vertices.empty() || its.indices.empty()){
        assert(false);
        return;
    }

    Geometry& data = m_render_data.geometry;
    data.format = { Geometry::EPrimitiveType::Triangles, Geometry::EVertexLayout::P3N3, Geometry::EIndexType::UINT };
    data.vertices.reserve(3 * its.indices.size() * Geometry::vertex_stride_floats(data.format));
    data.indices.reserve(3 * its.indices.size() * Geometry::index_stride_bytes(data.format));

    // vertices + indices
    unsigned int vertices_counter = 0;
    for (uint32_t i = 0; i < its.indices.size(); ++i) {
        stl_triangle_vertex_indices face = its.indices[i];
        stl_vertex                  vertex[3] = { its.vertices[face[0]], its.vertices[face[1]], its.vertices[face[2]] };
        stl_vertex                  n = face_normal_normalized(vertex);
        for (size_t j = 0; j < 3; ++j) {
            data.add_vertex(vertex[j], n);
        }
        vertices_counter += 3;
        data.add_uint_triangle(vertices_counter - 3, vertices_counter - 2, vertices_counter - 1);
    }

    // update bounding box
    for (size_t i = 0; i < vertices_count(); ++i) {
        m_bounding_box.merge(m_render_data.geometry.extract_position_3(i).cast<double>());
    }
#else
    if (!m_render_data.empty()) // call reset() if you want to reuse this model
        return;

    RenderData data;
    data.type = EPrimitiveType::Triangles;

    std::vector<float> vertices = std::vector<float>(18 * its.indices.size());
    std::vector<unsigned int> indices = std::vector<unsigned int>(3 * its.indices.size());

    unsigned int vertices_count = 0;
    for (uint32_t i = 0; i < its.indices.size(); ++i) {
        stl_triangle_vertex_indices face      = its.indices[i];
        stl_vertex                  vertex[3] = { its.vertices[face[0]], its.vertices[face[1]], its.vertices[face[2]] };
        stl_vertex                  n         = face_normal_normalized(vertex);
        for (size_t j = 0; j < 3; ++ j) {
            size_t offset = i * 18 + j * 6;
            ::memcpy(static_cast<void*>(&vertices[offset]), static_cast<const void*>(vertex[j].data()), 3 * sizeof(float));
            ::memcpy(static_cast<void*>(&vertices[3 + offset]), static_cast<const void*>(n.data()), 3 * sizeof(float));
        }
        for (size_t j = 0; j < 3; ++j)
            indices[i * 3 + j] = vertices_count + j;
        vertices_count += 3;
    }

    data.indices_count = static_cast<unsigned int>(indices.size());
    m_bounding_box = bbox;

    send_to_gpu(data, vertices, indices);
    m_render_data.emplace_back(data);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if !ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::init_from(const indexed_triangle_set& its)
{
    init_from(its, bounding_box(its));
}
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

void GLModel::init_from(const Polygons& polygons, float z)
{
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (is_initialized()) {
        // call reset() if you want to reuse this model
        assert(false);
        return;
    }

    if (polygons.empty()) {
        assert(false);
        return;
    }

    Geometry& data = m_render_data.geometry;
    data.format = { Geometry::EPrimitiveType::Lines, Geometry::EVertexLayout::P3, Geometry::EIndexType::UINT };

    size_t segments_count = 0;
    for (const Polygon& polygon : polygons) {
        segments_count += polygon.points.size();
    }

    data.vertices.reserve(2 * segments_count * Geometry::vertex_stride_floats(data.format));
    data.indices.reserve(2 * segments_count * Geometry::index_stride_bytes(data.format));

    // vertices + indices
    unsigned int vertices_counter = 0;
    for (const Polygon& poly : polygons) {
        for (size_t i = 0; i < poly.points.size(); ++i) {
            const Point& p0 = poly.points[i];
            const Point& p1 = (i == poly.points.size() - 1) ? poly.points.front() : poly.points[i + 1];
            data.add_vertex(Vec3f(unscale<float>(p0.x()), unscale<float>(p0.y()), z));
            data.add_vertex(Vec3f(unscale<float>(p1.x()), unscale<float>(p1.y()), z));
            vertices_counter += 2;
            data.add_uint_line(vertices_counter - 2, vertices_counter - 1);
        }
    }

    // update bounding box
    for (size_t i = 0; i < vertices_count(); ++i) {
        m_bounding_box.merge(m_render_data.geometry.extract_position_3(i).cast<double>());
    }
#else
    auto append_polygon = [](const Polygon& polygon, float z, GUI::GLModel::Geometry& data) {
        if (!polygon.empty()) {
            GUI::GLModel::Geometry::Entity entity;
            entity.type = GUI::GLModel::EPrimitiveType::LineLoop;
            // contour
            entity.positions.reserve(polygon.size() + 1);
            entity.indices.reserve(polygon.size() + 1);
            unsigned int id = 0;
            for (const Point& p : polygon) {
                Vec3f position = unscale(p.x(), p.y(), 0.0).cast<float>();
                position.z() = z;
                entity.positions.emplace_back(position);
                entity.indices.emplace_back(id++);
            }
            data.entities.emplace_back(entity);
        }
    };

    Geometry init_data;
    for (const Polygon& polygon : polygons) {
        append_polygon(polygon, z, init_data);
    }
    init_from(init_data);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

bool GLModel::init_from_file(const std::string& filename)
{
    if (!boost::filesystem::exists(filename))
        return false;

    if (!boost::algorithm::iends_with(filename, ".stl"))
        return false;

    Model model;
    try {
        model = Model::read_from_file(filename);
    }
    catch (std::exception&) {
        return false;
    }

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    init_from(model.mesh());
#else
    const TriangleMesh& mesh = model.mesh();
    init_from(mesh.its, mesh.bounding_box());
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    m_filename = filename;

    return true;
}

#if !ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::set_color(int entity_id, const ColorRGBA& color)
{
    for (size_t i = 0; i < m_render_data.size(); ++i) {
        if (entity_id == -1 || static_cast<int>(i) == entity_id)
            m_render_data[i].color = color;
    }
}

ColorRGBA GLModel::get_color(size_t entity_id) const
{
    if (entity_id < 0 || entity_id >= m_render_data.size()) return ColorRGBA{};
    return m_render_data[entity_id].color;
}
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

void GLModel::reset()
{
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    // release gpu memory
    if (m_render_data.ibo_id > 0) {
        glsafe(::glDeleteBuffers(1, &m_render_data.ibo_id));
        m_render_data.ibo_id = 0;
    }
    if (m_render_data.vbo_id > 0) {
        glsafe(::glDeleteBuffers(1, &m_render_data.vbo_id));
        m_render_data.vbo_id = 0;
    }

    m_render_data.vertices_count = 0;
    m_render_data.indices_count  = 0;
    m_render_data.geometry.vertices = std::vector<float>();
    m_render_data.geometry.indices  = std::vector<unsigned char>();
#else
    for (RenderData& data : m_render_data) {
        // release gpu memory
        if (data.ibo_id > 0)
            glsafe(::glDeleteBuffers(1, &data.ibo_id));
        if (data.vbo_id > 0)
            glsafe(::glDeleteBuffers(1, &data.vbo_id));
    }

    m_render_data.clear();
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    m_bounding_box = BoundingBoxf3();
    m_filename = std::string();
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
static GLenum get_primitive_mode(const GLModel::Geometry::Format& format)
{
    switch (format.type)
    {
    case GLModel::Geometry::EPrimitiveType::Points:        { return GL_POINTS; }
    default:
    case GLModel::Geometry::EPrimitiveType::Triangles:     { return GL_TRIANGLES; }
    case GLModel::Geometry::EPrimitiveType::TriangleStrip: { return GL_TRIANGLE_STRIP; }
    case GLModel::Geometry::EPrimitiveType::TriangleFan:   { return GL_TRIANGLE_FAN; }
    case GLModel::Geometry::EPrimitiveType::Lines:         { return GL_LINES; }
    case GLModel::Geometry::EPrimitiveType::LineStrip:     { return GL_LINE_STRIP; }
    case GLModel::Geometry::EPrimitiveType::LineLoop:      { return GL_LINE_LOOP; }
    }
}

static GLenum get_index_type(const GLModel::Geometry::Format& format)
{
    switch (format.index_type)
    {
    default:
    case GLModel::Geometry::EIndexType::UINT:   { return GL_UNSIGNED_INT; }
    case GLModel::Geometry::EIndexType::USHORT: { return GL_UNSIGNED_SHORT; }
    }
}

void GLModel::render()
#else
void GLModel::render() const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
    GLShaderProgram* shader = wxGetApp().get_current_shader();

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (shader == nullptr)
        return;

    // sends data to gpu if not done yet
    if (m_render_data.vbo_id == 0 || m_render_data.ibo_id == 0) {
        if (m_render_data.geometry.vertices_count() > 0 && m_render_data.geometry.indices_count() > 0 && !send_to_gpu())
            return;
    }

    const Geometry& data = m_render_data.geometry;

    GLenum mode = get_primitive_mode(data.format);
    GLenum index_type = get_index_type(data.format);

    const size_t vertex_stride_bytes = Geometry::vertex_stride_bytes(data.format);
    const bool position  = Geometry::has_position(data.format);
    const bool normal    = Geometry::has_normal(data.format);
    const bool tex_coord = Geometry::has_tex_coord(data.format);

    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, m_render_data.vbo_id));

    if (position) {
        glsafe(::glVertexPointer(Geometry::position_stride_floats(data.format), GL_FLOAT, vertex_stride_bytes, (const void*)Geometry::position_offset_bytes(data.format)));
        glsafe(::glEnableClientState(GL_VERTEX_ARRAY));
    }
    if (normal) {
        glsafe(::glNormalPointer(GL_FLOAT, vertex_stride_bytes, (const void*)Geometry::normal_offset_bytes(data.format)));
        glsafe(::glEnableClientState(GL_NORMAL_ARRAY));
    }
    if (tex_coord) {
        glsafe(::glTexCoordPointer(Geometry::tex_coord_stride_floats(data.format), GL_FLOAT, vertex_stride_bytes, (const void*)Geometry::tex_coord_offset_bytes(data.format)));
        glsafe(::glEnableClientState(GL_TEXTURE_COORD_ARRAY));
    }

    shader->set_uniform("uniform_color", data.color);

    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_render_data.ibo_id));
    glsafe(::glDrawElements(mode, indices_count(), index_type, nullptr));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

    if (tex_coord)
        glsafe(::glDisableClientState(GL_TEXTURE_COORD_ARRAY));
    if (normal)
        glsafe(::glDisableClientState(GL_NORMAL_ARRAY));
    if (position)
        glsafe(::glDisableClientState(GL_VERTEX_ARRAY));

    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
#else
    for (const RenderData& data : m_render_data) {
        if (data.vbo_id == 0 || data.ibo_id == 0)
            continue;

        GLenum mode;
        switch (data.type)
        {
        default:
        case EPrimitiveType::Triangles: { mode = GL_TRIANGLES; break; }
        case EPrimitiveType::Lines:     { mode = GL_LINES; break; }
        case EPrimitiveType::LineStrip: { mode = GL_LINE_STRIP; break; }
        case EPrimitiveType::LineLoop:  { mode = GL_LINE_LOOP; break; }
        }

        glsafe(::glBindBuffer(GL_ARRAY_BUFFER, data.vbo_id));
        glsafe(::glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), (const void*)0));
        glsafe(::glNormalPointer(GL_FLOAT, 6 * sizeof(float), (const void*)(3 * sizeof(float))));

        glsafe(::glEnableClientState(GL_VERTEX_ARRAY));
        glsafe(::glEnableClientState(GL_NORMAL_ARRAY));

        if (shader != nullptr)
            shader->set_uniform("uniform_color", data.color);
        else
            glsafe(::glColor4fv(data.color.data()));

        glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data.ibo_id));
        glsafe(::glDrawElements(mode, static_cast<GLsizei>(data.indices_count), GL_UNSIGNED_INT, (const void*)0));
        glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

        glsafe(::glDisableClientState(GL_NORMAL_ARRAY));
        glsafe(::glDisableClientState(GL_VERTEX_ARRAY));

        glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
    }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
void GLModel::render_instanced(unsigned int instances_vbo, unsigned int instances_count)
#else
void GLModel::render_instanced(unsigned int instances_vbo, unsigned int instances_count) const
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
{
    if (instances_vbo == 0)
        return;

    GLShaderProgram* shader = wxGetApp().get_current_shader();
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    if (shader == nullptr || !boost::algorithm::iends_with(shader->get_name(), "_instanced"))
        return;

    // vertex attributes
    GLint position_id = shader->get_attrib_location("v_position");
    GLint normal_id   = shader->get_attrib_location("v_normal");
    if (position_id == -1 || normal_id == -1)
        return;

    // instance attributes
    GLint offset_id = shader->get_attrib_location("i_offset");
    GLint scales_id = shader->get_attrib_location("i_scales");
    if (offset_id == -1 || scales_id == -1)
        return;

    if (m_render_data.vbo_id == 0 || m_render_data.ibo_id == 0) {
        if (!send_to_gpu())
            return;
    }
#else
    assert(shader == nullptr || boost::algorithm::iends_with(shader->get_name(), "_instanced"));

    // vertex attributes
    GLint position_id = (shader != nullptr) ? shader->get_attrib_location("v_position") : -1;
    GLint normal_id = (shader != nullptr) ? shader->get_attrib_location("v_normal") : -1;
    assert(position_id != -1 && normal_id != -1);

    // instance attributes
    GLint offset_id = (shader != nullptr) ? shader->get_attrib_location("i_offset") : -1;
    GLint scales_id = (shader != nullptr) ? shader->get_attrib_location("i_scales") : -1;
    assert(offset_id != -1 && scales_id != -1);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, instances_vbo));
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    glsafe(::glVertexAttribPointer(offset_id, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (GLvoid*)0));
    glsafe(::glEnableVertexAttribArray(offset_id));
    glsafe(::glVertexAttribDivisor(offset_id, 1));

    glsafe(::glVertexAttribPointer(scales_id, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (GLvoid*)(3 * sizeof(float))));
    glsafe(::glEnableVertexAttribArray(scales_id));
    glsafe(::glVertexAttribDivisor(scales_id, 1));
#else
    if (offset_id != -1) {
        glsafe(::glVertexAttribPointer(offset_id, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (GLvoid*)0));
        glsafe(::glEnableVertexAttribArray(offset_id));
        glsafe(::glVertexAttribDivisor(offset_id, 1));
    }
    if (scales_id != -1) {
        glsafe(::glVertexAttribPointer(scales_id, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (GLvoid*)(3 * sizeof(float))));
        glsafe(::glEnableVertexAttribArray(scales_id));
        glsafe(::glVertexAttribDivisor(scales_id, 1));
    }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    const Geometry& data = m_render_data.geometry;

    GLenum mode = get_primitive_mode(data.format);
    GLenum index_type = get_index_type(data.format);

    shader->set_uniform("uniform_color", data.color);

    const size_t vertex_stride_bytes = Geometry::vertex_stride_bytes(data.format);
    const bool position = Geometry::has_position(data.format);
    const bool normal   = Geometry::has_normal(data.format);

    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, m_render_data.vbo_id));

    if (position) {
        glsafe(::glVertexAttribPointer(position_id, Geometry::position_stride_floats(data.format), GL_FLOAT, GL_FALSE, vertex_stride_bytes, (GLvoid*)Geometry::position_offset_bytes(data.format)));
        glsafe(::glEnableVertexAttribArray(position_id));
    }

    if (normal) {
        glsafe(::glVertexAttribPointer(normal_id, Geometry::normal_stride_floats(data.format), GL_FLOAT, GL_FALSE, vertex_stride_bytes, (GLvoid*)Geometry::normal_offset_bytes(data.format)));
        glsafe(::glEnableVertexAttribArray(normal_id));
    }

    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_render_data.ibo_id));
    glsafe(::glDrawElementsInstanced(mode, indices_count(), index_type, (const void*)0, instances_count));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

    if (normal)
        glsafe(::glDisableVertexAttribArray(normal_id));
    if (position)
        glsafe(::glDisableVertexAttribArray(position_id));

    glsafe(::glDisableVertexAttribArray(scales_id));
    glsafe(::glDisableVertexAttribArray(offset_id));
#else
    for (const RenderData& data : m_render_data) {
        if (data.vbo_id == 0 || data.ibo_id == 0)
            continue;

        GLenum mode;
        switch (data.type)
        {
        default:
        case EPrimitiveType::Triangles: { mode = GL_TRIANGLES; break; }
        case EPrimitiveType::Lines:     { mode = GL_LINES; break; }
        case EPrimitiveType::LineStrip: { mode = GL_LINE_STRIP; break; }
        case EPrimitiveType::LineLoop:  { mode = GL_LINE_LOOP; break; }
        }

        if (shader != nullptr)
            shader->set_uniform("uniform_color", data.color);
        else
            glsafe(::glColor4fv(data.color.data()));

        glsafe(::glBindBuffer(GL_ARRAY_BUFFER, data.vbo_id));
        if (position_id != -1) {
            glsafe(::glVertexAttribPointer(position_id, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (GLvoid*)0));
            glsafe(::glEnableVertexAttribArray(position_id));
        }
        if (normal_id != -1) {
            glsafe(::glVertexAttribPointer(normal_id, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (GLvoid*)(3 * sizeof(float))));
            glsafe(::glEnableVertexAttribArray(normal_id));
        }

        glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data.ibo_id));
        glsafe(::glDrawElementsInstanced(mode, static_cast<GLsizei>(data.indices_count), GL_UNSIGNED_INT, (const void*)0, instances_count));
        glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));

        if (normal_id != -1)
            glsafe(::glDisableVertexAttribArray(normal_id));
        if (position_id != -1)
            glsafe(::glDisableVertexAttribArray(position_id));
    }

    if (scales_id != -1)
        glsafe(::glDisableVertexAttribArray(scales_id));
    if (offset_id != -1)
        glsafe(::glDisableVertexAttribArray(offset_id));
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
}

#if ENABLE_GLBEGIN_GLEND_REMOVAL
bool GLModel::send_to_gpu()
{
    if (m_render_data.vbo_id > 0 || m_render_data.ibo_id > 0) {
        assert(false);
        return false;
    }

    Geometry& data = m_render_data.geometry;
    if (data.vertices.empty() || data.indices.empty()) {
        assert(false);
        return false;
    }

    // vertices
    glsafe(::glGenBuffers(1, &m_render_data.vbo_id));
    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, m_render_data.vbo_id));
    glsafe(::glBufferData(GL_ARRAY_BUFFER, data.vertices_size_bytes(), data.vertices.data(), GL_STATIC_DRAW));
    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));
    m_render_data.vertices_count = vertices_count();
    data.vertices = std::vector<float>();

    // indices
    glsafe(::glGenBuffers(1, &m_render_data.ibo_id));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_render_data.ibo_id));
    glsafe(::glBufferData(GL_ELEMENT_ARRAY_BUFFER, data.indices_size_bytes(), data.indices.data(), GL_STATIC_DRAW));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
    m_render_data.indices_count = indices_count();
    data.indices = std::vector<unsigned char>();

    return true;
}
#else
void GLModel::send_to_gpu(RenderData& data, const std::vector<float>& vertices, const std::vector<unsigned int>& indices)
{
    assert(data.vbo_id == 0);
    assert(data.ibo_id == 0);

    // vertex data -> send to gpu
    glsafe(::glGenBuffers(1, &data.vbo_id));
    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, data.vbo_id));
    glsafe(::glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW));
    glsafe(::glBindBuffer(GL_ARRAY_BUFFER, 0));

    // indices data -> send to gpu
    glsafe(::glGenBuffers(1, &data.ibo_id));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data.ibo_id));
    glsafe(::glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW));
    glsafe(::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
}
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

#if ENABLE_GLBEGIN_GLEND_REMOVAL
static void append_vertex(GLModel::Geometry& data, const Vec3f& position, const Vec3f& normal)
{
    data.add_vertex(position, normal);
}

static void append_triangle(GLModel::Geometry& data, unsigned short v1, unsigned short v2, unsigned short v3)
{
    data.add_ushort_index(v1);
    data.add_ushort_index(v2);
    data.add_ushort_index(v3);
}
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

GLModel::Geometry stilized_arrow(unsigned short resolution, float tip_radius, float tip_height, float stem_radius, float stem_height)
{
#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    auto append_vertex = [](GLModel::Geometry::Entity& entity, const Vec3f& position, const Vec3f& normal) {
        entity.positions.emplace_back(position);
        entity.normals.emplace_back(normal);
    };
    auto append_indices = [](GLModel::Geometry::Entity& entity, unsigned int v1, unsigned int v2, unsigned int v3) {
        entity.indices.emplace_back(v1);
        entity.indices.emplace_back(v2);
        entity.indices.emplace_back(v3);
    };
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

    resolution = std::max<unsigned short>(4, resolution);
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    resolution = std::min<unsigned short>(10922, resolution); // ensure no unsigned short overflow of indices
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    GLModel::Geometry data;
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3, GLModel::Geometry::EIndexType::USHORT };
    data.vertices.reserve((6 * resolution + 2) * GLModel::Geometry::vertex_stride_floats(data.format));
    data.indices.reserve((6 * resolution * 3) * GLModel::Geometry::index_stride_bytes(data.format));
#else
    GLModel::Geometry::Entity entity;
    entity.type = GLModel::EPrimitiveType::Triangles;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    const float angle_step = 2.0f * float(PI) / float(resolution);
    std::vector<float> cosines(resolution);
    std::vector<float> sines(resolution);

    for (unsigned short i = 0; i < resolution; ++i) {
        const float angle = angle_step * float(i);
        cosines[i] = ::cos(angle);
        sines[i] = -::sin(angle);
    }

    const float total_height = tip_height + stem_height;

    // tip vertices/normals
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    append_vertex(data, { 0.0f, 0.0f, total_height }, Vec3f::UnitZ());
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { tip_radius * sines[i], tip_radius * cosines[i], stem_height }, { sines[i], cosines[i], 0.0f });
    }

    // tip triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short v3 = (i < resolution - 1) ? i + 2 : 1;
        append_triangle(data, 0, i + 1, v3);
    }

    // tip cap outer perimeter vertices
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { tip_radius * sines[i], tip_radius * cosines[i], stem_height }, -Vec3f::UnitZ());
    }

    // tip cap inner perimeter vertices
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { stem_radius * sines[i], stem_radius * cosines[i], stem_height }, -Vec3f::UnitZ());
    }

    // tip cap triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short v2 = (i < resolution - 1) ? i + resolution + 2 : resolution + 1;
        const unsigned short v3 = (i < resolution - 1) ? i + 2 * resolution + 2 : 2 * resolution + 1;
        append_triangle(data, i + resolution + 1, v3, v2);
        append_triangle(data, i + resolution + 1, i + 2 * resolution + 1, v3);
    }

    // stem bottom vertices
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { stem_radius * sines[i], stem_radius * cosines[i], stem_height }, { sines[i], cosines[i], 0.0f });
    }

    // stem top vertices
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { stem_radius * sines[i], stem_radius * cosines[i], 0.0f }, { sines[i], cosines[i], 0.0f });
    }

    // stem triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short v2 = (i < resolution - 1) ? i + 3 * resolution + 2 : 3 * resolution + 1;
        const unsigned short v3 = (i < resolution - 1) ? i + 4 * resolution + 2 : 4 * resolution + 1;
        append_triangle(data, i + 3 * resolution + 1, v3, v2);
        append_triangle(data, i + 3 * resolution + 1, i + 4 * resolution + 1, v3);
    }

    // stem cap vertices
    append_vertex(data, Vec3f::Zero(), -Vec3f::UnitZ());
    for (unsigned short i = 0; i < resolution; ++i) {
        append_vertex(data, { stem_radius * sines[i], stem_radius * cosines[i], 0.0f }, -Vec3f::UnitZ());
    }

    // stem cap triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short v3 = (i < resolution - 1) ? i + 5 * resolution + 3 : 5 * resolution + 2;
        append_triangle(data, 5 * resolution + 1, v3, i + 5 * resolution + 2);
    }
#else
    append_vertex(entity, { 0.0f, 0.0f, total_height }, Vec3f::UnitZ());
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { tip_radius * sines[i], tip_radius * cosines[i], stem_height }, { sines[i], cosines[i], 0.0f });
    }

    // tip triangles
    for (int i = 0; i < resolution; ++i) {
        const int v3 = (i < resolution - 1) ? i + 2 : 1;
        append_indices(entity, 0, i + 1, v3);
    }

    // tip cap outer perimeter vertices
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { tip_radius * sines[i], tip_radius * cosines[i], stem_height }, -Vec3f::UnitZ());
    }

    // tip cap inner perimeter vertices
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { stem_radius * sines[i], stem_radius * cosines[i], stem_height }, -Vec3f::UnitZ());
    }

    // tip cap triangles
    for (int i = 0; i < resolution; ++i) {
        const int v2 = (i < resolution - 1) ? i + resolution + 2 : resolution + 1;
        const int v3 = (i < resolution - 1) ? i + 2 * resolution + 2 : 2 * resolution + 1;
        append_indices(entity, i + resolution + 1, v3, v2);
        append_indices(entity, i + resolution + 1, i + 2 * resolution + 1, v3);
    }

    // stem bottom vertices
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { stem_radius * sines[i], stem_radius * cosines[i], stem_height }, { sines[i], cosines[i], 0.0f });
    }

    // stem top vertices
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { stem_radius * sines[i], stem_radius * cosines[i], 0.0f }, { sines[i], cosines[i], 0.0f });
    }

    // stem triangles
    for (int i = 0; i < resolution; ++i) {
        const int v2 = (i < resolution - 1) ? i + 3 * resolution + 2 : 3 * resolution + 1;
        const int v3 = (i < resolution - 1) ? i + 4 * resolution + 2 : 4 * resolution + 1;
        append_indices(entity, i + 3 * resolution + 1, v3, v2);
        append_indices(entity, i + 3 * resolution + 1, i + 4 * resolution + 1, v3);
    }

    // stem cap vertices
    append_vertex(entity, Vec3f::Zero(), -Vec3f::UnitZ());
    for (int i = 0; i < resolution; ++i) {
        append_vertex(entity, { stem_radius * sines[i], stem_radius * cosines[i], 0.0f }, -Vec3f::UnitZ());
    }

    // stem cap triangles
    for (int i = 0; i < resolution; ++i) {
        const int v3 = (i < resolution - 1) ? i + 5 * resolution + 3 : 5 * resolution + 2;
        append_indices(entity, 5 * resolution + 1, v3, i + 5 * resolution + 2);
    }

    data.entities.emplace_back(entity);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    return data;
}

GLModel::Geometry circular_arrow(unsigned short resolution, float radius, float tip_height, float tip_width, float stem_width, float thickness)
{
#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    auto append_vertex = [](GLModel::Geometry::Entity& entity, const Vec3f& position, const Vec3f& normal) {
        entity.positions.emplace_back(position);
        entity.normals.emplace_back(normal);
    };
    auto append_indices = [](GLModel::Geometry::Entity& entity, unsigned int v1, unsigned int v2, unsigned int v3) {
        entity.indices.emplace_back(v1);
        entity.indices.emplace_back(v2);
        entity.indices.emplace_back(v3);
    };
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

    resolution = std::max<unsigned short>(2, resolution);
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    resolution = std::min<unsigned short>(8188, resolution); // ensure no unsigned short overflow of indices
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    GLModel::Geometry data;
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3, GLModel::Geometry::EIndexType::USHORT };
    data.vertices.reserve((8 * (resolution + 1) + 30) * GLModel::Geometry::vertex_stride_floats(data.format));
    data.indices.reserve(((8 * resolution + 16) * 3) * GLModel::Geometry::index_stride_bytes(data.format));
#else
    GLModel::Geometry::Entity entity;
    entity.type = GLModel::EPrimitiveType::Triangles;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    const float half_thickness = 0.5f * thickness;
    const float half_stem_width = 0.5f * stem_width;
    const float half_tip_width = 0.5f * tip_width;

    const float outer_radius = radius + half_stem_width;
    const float inner_radius = radius - half_stem_width;
    const float step_angle = 0.5f * float(PI) / float(resolution);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    // tip
    // top face vertices
    append_vertex(data, { 0.0f, outer_radius, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { 0.0f, radius + half_tip_width, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { -tip_height, radius, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { 0.0f, radius - half_tip_width, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { 0.0f, inner_radius, half_thickness }, Vec3f::UnitZ());

    // top face triangles
    append_triangle(data, 0, 1, 2);
    append_triangle(data, 0, 2, 4);
    append_triangle(data, 4, 2, 3);

    // bottom face vertices
    append_vertex(data, { 0.0f, outer_radius, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { 0.0f, radius + half_tip_width, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { -tip_height, radius, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { 0.0f, radius - half_tip_width, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { 0.0f, inner_radius, -half_thickness }, -Vec3f::UnitZ());

    // bottom face triangles
    append_triangle(data, 5, 7, 6);
    append_triangle(data, 5, 9, 7);
    append_triangle(data, 9, 8, 7);

    // side faces vertices
    append_vertex(data, { 0.0f, outer_radius, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, radius + half_tip_width, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, outer_radius, half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, radius + half_tip_width, half_thickness }, Vec3f::UnitX());

    Vec3f normal(-half_tip_width, tip_height, 0.0f);
    normal.normalize();
    append_vertex(data, { 0.0f, radius + half_tip_width, -half_thickness }, normal);
    append_vertex(data, { -tip_height, radius, -half_thickness }, normal);
    append_vertex(data, { 0.0f, radius + half_tip_width, half_thickness }, normal);
    append_vertex(data, { -tip_height, radius, half_thickness }, normal);

    normal = { -half_tip_width, -tip_height, 0.0f };
    normal.normalize();
    append_vertex(data, { -tip_height, radius, -half_thickness }, normal);
    append_vertex(data, { 0.0f, radius - half_tip_width, -half_thickness }, normal);
    append_vertex(data, { -tip_height, radius, half_thickness }, normal);
    append_vertex(data, { 0.0f, radius - half_tip_width, half_thickness }, normal);

    append_vertex(data, { 0.0f, radius - half_tip_width, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, inner_radius, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, radius - half_tip_width, half_thickness }, Vec3f::UnitX());
    append_vertex(data, { 0.0f, inner_radius, half_thickness }, Vec3f::UnitX());

    // side face triangles
    for (unsigned short i = 0; i < 4; ++i) {
        const unsigned short ii = i * 4;
        append_triangle(data, 10 + ii, 11 + ii, 13 + ii);
        append_triangle(data, 10 + ii, 13 + ii, 12 + ii);
    }

    // stem
    // top face vertices
    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        append_vertex(data, { inner_radius * ::sin(angle), inner_radius * ::cos(angle), half_thickness }, Vec3f::UnitZ());
    }

    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        append_vertex(data, { outer_radius * ::sin(angle), outer_radius * ::cos(angle), half_thickness }, Vec3f::UnitZ());
    }

    // top face triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        append_triangle(data, 26 + i, 27 + i, 27 + resolution + i);
        append_triangle(data, 27 + i, 28 + resolution + i, 27 + resolution + i);
    }

    // bottom face vertices
    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        append_vertex(data, { inner_radius * ::sin(angle), inner_radius * ::cos(angle), -half_thickness }, -Vec3f::UnitZ());
    }

    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        append_vertex(data, { outer_radius * ::sin(angle), outer_radius * ::cos(angle), -half_thickness }, -Vec3f::UnitZ());
    }

    // bottom face triangles
    for (unsigned short i = 0; i < resolution; ++i) {
        append_triangle(data, 28 + 2 * resolution + i, 29 + 3 * resolution + i, 29 + 2 * resolution + i);
        append_triangle(data, 29 + 2 * resolution + i, 29 + 3 * resolution + i, 30 + 3 * resolution + i);
    }

    // side faces vertices and triangles
    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(data, { inner_radius * s, inner_radius * c, -half_thickness }, { -s, -c, 0.0f });
    }

    for (unsigned short i = 0; i <= resolution; ++i) {
        const float angle = float(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(data, { inner_radius * s, inner_radius * c, half_thickness }, { -s, -c, 0.0f });
    }

    unsigned short first_id = 26 + 4 * (resolution + 1);
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short ii = first_id + i;
        append_triangle(data, ii, ii + 1, ii + resolution + 2);
        append_triangle(data, ii, ii + resolution + 2, ii + resolution + 1);
    }

    append_vertex(data, { inner_radius, 0.0f, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { outer_radius, 0.0f, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { inner_radius, 0.0f, half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { outer_radius, 0.0f, half_thickness }, -Vec3f::UnitY());

    first_id = 26 + 6 * (resolution + 1);
    append_triangle(data, first_id, first_id + 1, first_id + 3);
    append_triangle(data, first_id, first_id + 3, first_id + 2);

    for (short i = resolution; i >= 0; --i) {
        const float angle = float(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(data, { outer_radius * s, outer_radius * c, -half_thickness }, { s, c, 0.0f });
    }

    for (short i = resolution; i >= 0; --i) {
        const float angle = float(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(data, { outer_radius * s, outer_radius * c, +half_thickness }, { s, c, 0.0f });
    }

    first_id = 30 + 6 * (resolution + 1);
    for (unsigned short i = 0; i < resolution; ++i) {
        const unsigned short ii = first_id + i;
        append_triangle(data, ii, ii + 1, ii + resolution + 2);
        append_triangle(data, ii, ii + resolution + 2, ii + resolution + 1);
    }
#else
    // tip
    // top face vertices
    append_vertex(entity, { 0.0f, outer_radius, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, radius + half_tip_width, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { -tip_height, radius, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, radius - half_tip_width, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, inner_radius, half_thickness }, Vec3f::UnitZ());

    // top face triangles
    append_indices(entity, 0, 1, 2);
    append_indices(entity, 0, 2, 4);
    append_indices(entity, 4, 2, 3);

    // bottom face vertices
    append_vertex(entity, { 0.0f, outer_radius, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, radius + half_tip_width, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { -tip_height, radius, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, radius - half_tip_width, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { 0.0f, inner_radius, -half_thickness }, -Vec3f::UnitZ());

    // bottom face triangles
    append_indices(entity, 5, 7, 6);
    append_indices(entity, 5, 9, 7);
    append_indices(entity, 9, 8, 7);

    // side faces vertices
    append_vertex(entity, { 0.0f, outer_radius, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, radius + half_tip_width, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, outer_radius, half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, radius + half_tip_width, half_thickness }, Vec3f::UnitX());

    Vec3f normal(-half_tip_width, tip_height, 0.0f);
    normal.normalize();
    append_vertex(entity, { 0.0f, radius + half_tip_width, -half_thickness }, normal);
    append_vertex(entity, { -tip_height, radius, -half_thickness }, normal);
    append_vertex(entity, { 0.0f, radius + half_tip_width, half_thickness }, normal);
    append_vertex(entity, { -tip_height, radius, half_thickness }, normal);

    normal = Vec3f(-half_tip_width, -tip_height, 0.0f);
    normal.normalize();
    append_vertex(entity, { -tip_height, radius, -half_thickness }, normal);
    append_vertex(entity, { 0.0f, radius - half_tip_width, -half_thickness }, normal);
    append_vertex(entity, { -tip_height, radius, half_thickness }, normal);
    append_vertex(entity, { 0.0f, radius - half_tip_width, half_thickness }, normal);

    append_vertex(entity, { 0.0f, radius - half_tip_width, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, inner_radius, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, radius - half_tip_width, half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { 0.0f, inner_radius, half_thickness }, Vec3f::UnitX());

    // side face triangles
    for (int i = 0; i < 4; ++i) {
        const int ii = i * 4;
        append_indices(entity, 10 + ii, 11 + ii, 13 + ii);
        append_indices(entity, 10 + ii, 13 + ii, 12 + ii);
    }

    // stem
    // top face vertices
    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        append_vertex(entity, { inner_radius * ::sin(angle), inner_radius * ::cos(angle), half_thickness }, Vec3f::UnitZ());
    }

    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        append_vertex(entity, { outer_radius * ::sin(angle), outer_radius * ::cos(angle), half_thickness }, Vec3f::UnitZ());
    }

    // top face triangles
    for (int i = 0; i < resolution; ++i) {
        append_indices(entity, 26 + i, 27 + i, 27 + resolution + i);
        append_indices(entity, 27 + i, 28 + resolution + i, 27 + resolution + i);
    }

    // bottom face vertices
    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        append_vertex(entity, { inner_radius * ::sin(angle), inner_radius * ::cos(angle), -half_thickness }, -Vec3f::UnitZ());
    }

    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        append_vertex(entity, { outer_radius * ::sin(angle), outer_radius * ::cos(angle), -half_thickness }, -Vec3f::UnitZ());
    }

    // bottom face triangles
    for (int i = 0; i < resolution; ++i) {
        append_indices(entity, 28 + 2 * resolution + i, 29 + 3 * resolution + i, 29 + 2 * resolution + i);
        append_indices(entity, 29 + 2 * resolution + i, 29 + 3 * resolution + i, 30 + 3 * resolution + i);
    }

    // side faces vertices and triangles
    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(entity, { inner_radius * s, inner_radius * c, -half_thickness }, { -s, -c, 0.0f });
    }

    for (int i = 0; i <= resolution; ++i) {
        const float angle = static_cast<float>(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(entity, { inner_radius * s, inner_radius * c, half_thickness }, { -s, -c, 0.0f });
    }

    int first_id = 26 + 4 * (resolution + 1);
    for (int i = 0; i < resolution; ++i) {
        const int ii = first_id + i;
        append_indices(entity, ii, ii + 1, ii + resolution + 2);
        append_indices(entity, ii, ii + resolution + 2, ii + resolution + 1);
    }

    append_vertex(entity, { inner_radius, 0.0f, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { outer_radius, 0.0f, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { inner_radius, 0.0f, half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { outer_radius, 0.0f, half_thickness }, -Vec3f::UnitY());

    first_id = 26 + 6 * (resolution + 1);
    append_indices(entity, first_id, first_id + 1, first_id + 3);
    append_indices(entity, first_id, first_id + 3, first_id + 2);

    for (int i = resolution; i >= 0; --i) {
        const float angle = static_cast<float>(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(entity, { outer_radius * s, outer_radius * c, -half_thickness }, { s, c, 0.0f });
    }

    for (int i = resolution; i >= 0; --i) {
        const float angle = static_cast<float>(i) * step_angle;
        const float c = ::cos(angle);
        const float s = ::sin(angle);
        append_vertex(entity, { outer_radius * s, outer_radius * c, +half_thickness }, { s, c, 0.0f });
    }

    first_id = 30 + 6 * (resolution + 1);
    for (int i = 0; i < resolution; ++i) {
        const int ii = first_id + i;
        append_indices(entity, ii, ii + 1, ii + resolution + 2);
        append_indices(entity, ii, ii + resolution + 2, ii + resolution + 1);
    }

    data.entities.emplace_back(entity);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    return data;
}

GLModel::Geometry straight_arrow(float tip_width, float tip_height, float stem_width, float stem_height, float thickness)
{
#if !ENABLE_GLBEGIN_GLEND_REMOVAL
    auto append_vertex = [](GLModel::Geometry::Entity& entity, const Vec3f& position, const Vec3f& normal) {
        entity.positions.emplace_back(position);
        entity.normals.emplace_back(normal);
    };
    auto append_indices = [](GLModel::Geometry::Entity& entity, unsigned int v1, unsigned int v2, unsigned int v3) {
        entity.indices.emplace_back(v1);
        entity.indices.emplace_back(v2);
        entity.indices.emplace_back(v3);
    };
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

    GLModel::Geometry data;
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3, GLModel::Geometry::EIndexType::USHORT };
    data.vertices.reserve(42 * GLModel::Geometry::vertex_stride_floats(data.format));
    data.indices.reserve((24 * 3) * GLModel::Geometry::index_stride_bytes(data.format));
#else
    GLModel::Geometry::Entity entity;
    entity.type = GLModel::EPrimitiveType::Triangles;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    const float half_thickness = 0.5f * thickness;
    const float half_stem_width = 0.5f * stem_width;
    const float half_tip_width = 0.5f * tip_width;
    const float total_height = tip_height + stem_height;

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    // top face vertices
    append_vertex(data, { half_stem_width, 0.0, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { half_stem_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { half_tip_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { 0.0, total_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { -half_tip_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { -half_stem_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(data, { -half_stem_width, 0.0, half_thickness }, Vec3f::UnitZ());

    // top face triangles
    append_triangle(data, 0, 1, 6);
    append_triangle(data, 6, 1, 5);
    append_triangle(data, 4, 5, 3);
    append_triangle(data, 5, 1, 3);
    append_triangle(data, 1, 2, 3);

    // bottom face vertices
    append_vertex(data, { half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { 0.0, total_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { -half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(data, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitZ());

    // bottom face triangles
    append_triangle(data, 7, 13, 8);
    append_triangle(data, 13, 12, 8);
    append_triangle(data, 12, 11, 10);
    append_triangle(data, 8, 12, 10);
    append_triangle(data, 9, 8, 10);

    // side faces vertices
    append_vertex(data, { half_stem_width, 0.0, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { half_stem_width, stem_height, -half_thickness }, Vec3f::UnitX());
    append_vertex(data, { half_stem_width, 0.0, half_thickness }, Vec3f::UnitX());
    append_vertex(data, { half_stem_width, stem_height, half_thickness }, Vec3f::UnitX());

    append_vertex(data, { half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { half_stem_width, stem_height, half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { half_tip_width, stem_height, half_thickness }, -Vec3f::UnitY());

    Vec3f normal(tip_height, half_tip_width, 0.0f);
    normal.normalize();
    append_vertex(data, { half_tip_width, stem_height, -half_thickness }, normal);
    append_vertex(data, { 0.0, total_height, -half_thickness }, normal);
    append_vertex(data, { half_tip_width, stem_height, half_thickness }, normal);
    append_vertex(data, { 0.0, total_height, half_thickness }, normal);

    normal = { -tip_height, half_tip_width, 0.0f };
    normal.normalize();
    append_vertex(data, { 0.0, total_height, -half_thickness }, normal);
    append_vertex(data, { -half_tip_width, stem_height, -half_thickness }, normal);
    append_vertex(data, { 0.0, total_height, half_thickness }, normal);
    append_vertex(data, { -half_tip_width, stem_height, half_thickness }, normal);

    append_vertex(data, { -half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { -half_tip_width, stem_height, half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { -half_stem_width, stem_height, half_thickness }, -Vec3f::UnitY());

    append_vertex(data, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitX());
    append_vertex(data, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitX());
    append_vertex(data, { -half_stem_width, stem_height, half_thickness }, -Vec3f::UnitX());
    append_vertex(data, { -half_stem_width, 0.0, half_thickness }, -Vec3f::UnitX());

    append_vertex(data, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { -half_stem_width, 0.0, half_thickness }, -Vec3f::UnitY());
    append_vertex(data, { half_stem_width, 0.0, half_thickness }, -Vec3f::UnitY());

    // side face triangles
    for (unsigned short i = 0; i < 7; ++i) {
        const unsigned short ii = i * 4;
        append_triangle(data, 14 + ii, 15 + ii, 17 + ii);
        append_triangle(data, 14 + ii, 17 + ii, 16 + ii);
    }
#else
    // top face vertices
    append_vertex(entity, { half_stem_width, 0.0, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { half_stem_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { half_tip_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { 0.0, total_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { -half_tip_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { -half_stem_width, stem_height, half_thickness }, Vec3f::UnitZ());
    append_vertex(entity, { -half_stem_width, 0.0, half_thickness }, Vec3f::UnitZ());

    // top face triangles
    append_indices(entity, 0, 1, 6);
    append_indices(entity, 6, 1, 5);
    append_indices(entity, 4, 5, 3);
    append_indices(entity, 5, 1, 3);
    append_indices(entity, 1, 2, 3);

    // bottom face vertices
    append_vertex(entity, { half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { 0.0, total_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { -half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitZ());
    append_vertex(entity, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitZ());

    // bottom face triangles
    append_indices(entity, 7, 13, 8);
    append_indices(entity, 13, 12, 8);
    append_indices(entity, 12, 11, 10);
    append_indices(entity, 8, 12, 10);
    append_indices(entity, 9, 8, 10);

    // side faces vertices
    append_vertex(entity, { half_stem_width, 0.0, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { half_stem_width, stem_height, -half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { half_stem_width, 0.0, half_thickness }, Vec3f::UnitX());
    append_vertex(entity, { half_stem_width, stem_height, half_thickness }, Vec3f::UnitX());

    append_vertex(entity, { half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { half_stem_width, stem_height, half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { half_tip_width, stem_height, half_thickness }, -Vec3f::UnitY());

    Vec3f normal(tip_height, half_tip_width, 0.0f);
    normal.normalize();
    append_vertex(entity, { half_tip_width, stem_height, -half_thickness }, normal);
    append_vertex(entity, { 0.0, total_height, -half_thickness }, normal);
    append_vertex(entity, { half_tip_width, stem_height, half_thickness }, normal);
    append_vertex(entity, { 0.0, total_height, half_thickness }, normal);

    normal = Vec3f(-tip_height, half_tip_width, 0.0f);
    normal.normalize();
    append_vertex(entity, { 0.0, total_height, -half_thickness }, normal);
    append_vertex(entity, { -half_tip_width, stem_height, -half_thickness }, normal);
    append_vertex(entity, { 0.0, total_height, half_thickness }, normal);
    append_vertex(entity, { -half_tip_width, stem_height, half_thickness }, normal);

    append_vertex(entity, { -half_tip_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { -half_tip_width, stem_height, half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { -half_stem_width, stem_height, half_thickness }, -Vec3f::UnitY());

    append_vertex(entity, { -half_stem_width, stem_height, -half_thickness }, -Vec3f::UnitX());
    append_vertex(entity, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitX());
    append_vertex(entity, { -half_stem_width, stem_height, half_thickness }, -Vec3f::UnitX());
    append_vertex(entity, { -half_stem_width, 0.0, half_thickness }, -Vec3f::UnitX());

    append_vertex(entity, { -half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { half_stem_width, 0.0, -half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { -half_stem_width, 0.0, half_thickness }, -Vec3f::UnitY());
    append_vertex(entity, { half_stem_width, 0.0, half_thickness }, -Vec3f::UnitY());

    // side face triangles
    for (int i = 0; i < 7; ++i) {
        const int ii = i * 4;
        append_indices(entity, 14 + ii, 15 + ii, 17 + ii);
        append_indices(entity, 14 + ii, 17 + ii, 16 + ii);
    }

    data.entities.emplace_back(entity);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    return data;
}

GLModel::Geometry diamond(unsigned short resolution)
{
    resolution = std::max<unsigned short>(4, resolution);
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    resolution = std::min<unsigned short>(65534, resolution); // ensure no unsigned short overflow of indices
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    GLModel::Geometry data;
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3, GLModel::Geometry::EIndexType::USHORT };
    data.vertices.reserve((resolution + 2) * GLModel::Geometry::vertex_stride_floats(data.format));
    data.indices.reserve(((2 * (resolution + 1)) * 3) * GLModel::Geometry::index_stride_bytes(data.format));
#else
    GLModel::Geometry::Entity entity;
    entity.type = GLModel::EPrimitiveType::Triangles;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    const float step = 2.0f * float(PI) / float(resolution);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    // vertices
    for (unsigned short i = 0; i < resolution; ++i) {
        float ii = float(i) * step;
        const Vec3f p = { 0.5f * ::cos(ii), 0.5f * ::sin(ii), 0.0f };
        append_vertex(data, p, p.normalized());
    }
    Vec3f p = { 0.0f, 0.0f, 0.5f };
    append_vertex(data, p, p.normalized());
    p = { 0.0f, 0.0f, -0.5f };
    append_vertex(data, p, p.normalized());

    // triangles
    // top
    for (unsigned short i = 0; i < resolution; ++i) {
        append_triangle(data, i + 0, i + 1, resolution);
    }
    append_triangle(data, resolution - 1, 0, resolution);

    // bottom
    for (unsigned short i = 0; i < resolution; ++i) {
        append_triangle(data, i + 0, resolution + 1, i + 1);
    }
    append_triangle(data, resolution - 1, resolution + 1, 0);
#else
    // positions
    for (int i = 0; i < resolution; ++i) {
        float ii = float(i) * step;
        entity.positions.emplace_back(0.5f * ::cos(ii), 0.5f * ::sin(ii), 0.0f);
    }
    entity.positions.emplace_back(0.0f, 0.0f, 0.5f);
    entity.positions.emplace_back(0.0f, 0.0f, -0.5f);

    // normals
    for (const Vec3f& v : entity.positions) {
        entity.normals.emplace_back(v.normalized());
    }

    // triangles
    // top
    for (int i = 0; i < resolution; ++i) {
        entity.indices.push_back(i + 0);
        entity.indices.push_back(i + 1);
        entity.indices.push_back(resolution);
    }
    entity.indices.push_back(resolution - 1);
    entity.indices.push_back(0);
    entity.indices.push_back(resolution);

    // bottom
    for (int i = 0; i < resolution; ++i) {
        entity.indices.push_back(i + 0);
        entity.indices.push_back(resolution + 1);
        entity.indices.push_back(i + 1);
    }
    entity.indices.push_back(resolution - 1);
    entity.indices.push_back(resolution + 1);
    entity.indices.push_back(0);

    data.entities.emplace_back(entity);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    return data;
}

} // namespace GUI
} // namespace Slic3r
