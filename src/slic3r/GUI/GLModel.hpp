#ifndef slic3r_GLModel_hpp_
#define slic3r_GLModel_hpp_

#include "libslic3r/Point.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Color.hpp"
#include <vector>
#include <string>

struct indexed_triangle_set;

namespace Slic3r {

class TriangleMesh;
class Polygon;
using Polygons = std::vector<Polygon>;

namespace GUI {

    class GLModel
    {
    public:
#if !ENABLE_GLBEGIN_GLEND_REMOVAL
        enum class EPrimitiveType : unsigned char
        {
            Triangles,
            Lines,
            LineStrip,
            LineLoop
        };

        struct RenderData
        {
            EPrimitiveType type;
            unsigned int vbo_id{ 0 };
            unsigned int ibo_id{ 0 };
            size_t indices_count{ 0 };
            ColorRGBA color;
        };
#endif // !ENABLE_GLBEGIN_GLEND_REMOVAL

        struct Geometry
        {
#if ENABLE_GLBEGIN_GLEND_REMOVAL
            enum class EPrimitiveType : unsigned char
            {
                Points,
                Triangles,
                TriangleStrip,
                TriangleFan,
                Lines,
                LineStrip,
                LineLoop
            };

            enum class EVertexLayout : unsigned char
            {
                P2,   // position 2 floats
                P2T2, // position 2 floats + texture coords 2 floats
                P3,   // position 3 floats
                P3T2, // position 3 floats + texture coords 2 floats
                P3N3, // position 3 floats + normal 3 floats
            };

            enum class EIndexType : unsigned char
            {
                UINT,   // unsigned int
                USHORT  // unsigned short
            };

            struct Format
            {
                EPrimitiveType type{ EPrimitiveType::Triangles };
                EVertexLayout vertex_layout{ EVertexLayout::P3N3 };
                EIndexType index_type{ EIndexType::UINT };
            };

            Format format;
            std::vector<float> vertices;
            std::vector<unsigned char> indices;
            ColorRGBA color{ ColorRGBA::BLACK() };

            void reserve_vertices(size_t vertices_count);
            void reserve_indices(size_t indices_count);

            void add_vertex(const Vec2f& position);                          // EVertexLayout::P2
            void add_vertex(const Vec2f& position, const Vec2f& tex_coord);  // EVertexLayout::P2T2
            void add_vertex(const Vec3f& position);                          // EVertexLayout::P3
            void add_vertex(const Vec3f& position, const Vec2f& tex_coord);  // EVertexLayout::P3T2
            void add_vertex(const Vec3f& position, const Vec3f& normal);     // EVertexLayout::P3N3

            void add_ushort_index(unsigned short id);
            void add_uint_index(unsigned int id);

            void add_ushort_line(unsigned short id1, unsigned short id2);
            void add_uint_line(unsigned int id1, unsigned int id2);

            void add_ushort_triangle(unsigned short id1, unsigned short id2, unsigned short id3);
            void add_uint_triangle(unsigned int id1, unsigned int id2, unsigned int id3);

            Vec2f extract_position_2(size_t id) const;
            Vec3f extract_position_3(size_t id) const;
            Vec3f extract_normal_3(size_t id) const;
            Vec2f extract_tex_coord_2(size_t id) const;

            unsigned int extract_uint_index(size_t id) const;
            unsigned short extract_ushort_index(size_t id) const;

            bool is_empty() const { return vertices.empty() || indices.empty(); }

            size_t vertices_count() const { return vertices.size() / vertex_stride_floats(format); }
            size_t indices_count() const  { return indices.size() / index_stride_bytes(format); }

            size_t vertices_size_floats() const { return vertices.size(); }
            size_t vertices_size_bytes() const  { return vertices_size_floats() * sizeof(float); }
            size_t indices_size_bytes() const   { return indices.size(); }

            static size_t vertex_stride_floats(const Format& format);
            static size_t vertex_stride_bytes(const Format& format) { return vertex_stride_floats(format) * sizeof(float); }

            static size_t position_stride_floats(const Format& format);
            static size_t position_stride_bytes(const Format& format) { return position_stride_floats(format) * sizeof(float); }
            static size_t position_offset_floats(const Format& format);
            static size_t position_offset_bytes(const Format& format) { return position_offset_floats(format) * sizeof(float); }

            static size_t normal_stride_floats(const Format& format);
            static size_t normal_stride_bytes(const Format& format) { return normal_stride_floats(format) * sizeof(float); }
            static size_t normal_offset_floats(const Format& format);
            static size_t normal_offset_bytes(const Format& format) { return normal_offset_floats(format) * sizeof(float); }

            static size_t tex_coord_stride_floats(const Format& format);
            static size_t tex_coord_stride_bytes(const Format& format) { return tex_coord_stride_floats(format) * sizeof(float); }
            static size_t tex_coord_offset_floats(const Format& format);
            static size_t tex_coord_offset_bytes(const Format& format) { return tex_coord_offset_floats(format) * sizeof(float); }

            static size_t index_stride_bytes(const Format& format);

            static EIndexType index_type(size_t vertices_count);

            static bool has_position(const Format& format);
            static bool has_normal(const Format& format);
            static bool has_tex_coord(const Format& format);
#else
            struct Entity
            {
                EPrimitiveType type;
                std::vector<Vec3f> positions;
                std::vector<Vec3f> normals;
                std::vector<unsigned int> indices;
                ColorRGBA color;
            };

            std::vector<Entity> entities;

            size_t vertices_count() const;
            size_t vertices_size_floats() const { return vertices_count() * 6; }
            size_t vertices_size_bytes() const { return vertices_size_floats() * sizeof(float); }

            size_t indices_count() const;
            size_t indices_size_bytes() const { return indices_count() * sizeof(unsigned int); }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
        };

#if ENABLE_GLBEGIN_GLEND_REMOVAL
        struct RenderData
        {
            Geometry geometry;
            unsigned int vbo_id{ 0 };
            unsigned int ibo_id{ 0 };
            size_t vertices_count{ 0 };
            size_t indices_count{ 0 };
        };
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    private:
#if ENABLE_GLBEGIN_GLEND_REMOVAL
        RenderData m_render_data;
#else
        std::vector<RenderData> m_render_data;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

        BoundingBoxf3 m_bounding_box;
        std::string m_filename;

    public:
        GLModel() = default;
        virtual ~GLModel() { reset(); }

#if ENABLE_GLBEGIN_GLEND_REMOVAL
        size_t vertices_count() const { return m_render_data.vertices_count > 0 ?
            m_render_data.vertices_count : m_render_data.geometry.vertices_count(); }
        size_t indices_count() const { return m_render_data.indices_count > 0 ?
            m_render_data.indices_count : m_render_data.geometry.indices_count(); }

        size_t vertices_size_floats() const { return vertices_count() * Geometry::vertex_stride_floats(m_render_data.geometry.format); }
        size_t vertices_size_bytes() const  { return vertices_size_floats() * sizeof(float); }

        size_t indices_size_bytes() const { return indices_count() * Geometry::index_stride_bytes(m_render_data.geometry.format); }

        void init_from(Geometry&& data);
        void init_from(const TriangleMesh& mesh);
#else
        void init_from(const Geometry& data);
        void init_from(const indexed_triangle_set& its, const BoundingBoxf3& bbox);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
        void init_from(const indexed_triangle_set& its);
        void init_from(const Polygons& polygons, float z);
        bool init_from_file(const std::string& filename);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
        void set_color(const ColorRGBA& color) { m_render_data.geometry.color = color; }
        const ColorRGBA& get_color() const { return m_render_data.geometry.color; }
#else
        // if entity_id == -1 set the color of all entities
        void set_color(int entity_id, const ColorRGBA& color);
        ColorRGBA get_color(size_t entity_id = 0U) const;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

        void reset();
#if ENABLE_GLBEGIN_GLEND_REMOVAL
        void render();
        void render_instanced(unsigned int instances_vbo, unsigned int instances_count);

        bool is_initialized() const { return vertices_count() > 0 && indices_count() > 0; }
#else
        void render() const;
        void render_instanced(unsigned int instances_vbo, unsigned int instances_count) const;

        bool is_initialized() const { return !m_render_data.empty(); }
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

        const BoundingBoxf3& get_bounding_box() const { return m_bounding_box; }
        const std::string& get_filename() const { return m_filename; }

    private:
#if ENABLE_GLBEGIN_GLEND_REMOVAL
        bool send_to_gpu();
#else
        void send_to_gpu(RenderData& data, const std::vector<float>& vertices, const std::vector<unsigned int>& indices);
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
    };

    // create an arrow with cylindrical stem and conical tip, with the given dimensions and resolution
    // the origin of the arrow is in the center of the stem cap
    // the arrow has its axis of symmetry along the Z axis and is pointing upward
    // used to render bed axes and sequential marker
    GLModel::Geometry stilized_arrow(unsigned short resolution, float tip_radius, float tip_height, float stem_radius, float stem_height);

    // create an arrow whose stem is a quarter of circle, with the given dimensions and resolution
    // the origin of the arrow is in the center of the circle
    // the arrow is contained in the 1st quadrant of the XY plane and is pointing counterclockwise
    // used to render sidebar hints for rotations
    GLModel::Geometry circular_arrow(unsigned short resolution, float radius, float tip_height, float tip_width, float stem_width, float thickness);

    // create an arrow with the given dimensions
    // the origin of the arrow is in the center of the stem cap
    // the arrow is contained in XY plane and has its main axis along the Y axis
    // used to render sidebar hints for position and scale
    GLModel::Geometry straight_arrow(float tip_width, float tip_height, float stem_width, float stem_height, float thickness);

    // create a diamond with the given resolution
    // the origin of the diamond is in its center
    // the diamond is contained into a box with size [1, 1, 1]
    GLModel::Geometry diamond(unsigned short resolution);

#if ENABLE_GLBEGIN_GLEND_REMOVAL
#if ENABLE_SHOW_TOOLPATHS_COG
    // create a sphere with the given resolution and smooth normals
    // the origin of the sphere is in its center
    GLModel::Geometry smooth_sphere(unsigned short resolution, float radius);
#endif // ENABLE_SHOW_TOOLPATHS_COG
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLModel_hpp_

