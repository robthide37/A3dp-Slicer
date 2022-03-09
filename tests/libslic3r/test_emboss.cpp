#include <catch2/catch.hpp>

#include <libslic3r/Emboss.hpp>
#include <libslic3r/SVG.hpp> // only debug visualization

#include <optional>
#include <libslic3r/AABBTreeIndirect.hpp>
#include <libslic3r/Utils.hpp> // for next_highest_power_of_2()

using namespace Slic3r;

namespace Private{
        
// calculate multiplication of ray dir to intersect - inspired by
// segment_segment_intersection when ray dir is normalized retur distance from
// ray point to intersection No value mean no intersection
std::optional<double> ray_segment_intersection(const Vec2d &r_point,
                                               const Vec2d &r_dir,
                                               const Vec2d &s0,
                                               const Vec2d &s1)
{
    auto denominate = [](const Vec2d &v0, const Vec2d &v1) -> double {
        return v0.x() * v1.y() - v1.x() * v0.y();
    };

    Vec2d  segment_dir = s1 - s0;
    double d           = denominate(segment_dir, r_dir);
    if (std::abs(d) < std::numeric_limits<double>::epsilon())
        // Line and ray are collinear.
        return {};

    Vec2d  s12         = s0 - r_point;
    double s_number    = denominate(r_dir, s12);
    bool   change_sign = false;
    if (d < 0.) {
        change_sign = true;
        d           = -d;
        s_number    = -s_number;
    }

    if (s_number < 0. || s_number > d)
        // Intersection outside of segment.
        return {};

    double r_number = denominate(segment_dir, s12);
    if (change_sign) r_number = -r_number;

    if (r_number < 0.)
        // Intersection before ray start.
        return {};

    return r_number / d;
}

Vec2d get_intersection(const Vec2d &               point,
                       const Vec2d &               dir,
                       const std::array<Vec2d, 3> &triangle)
{
    std::optional<double> t;
    for (size_t i = 0; i < 3; ++i) {
        size_t i2 = i + 1;
        if (i2 == 3) i2 = 0;
        if (!t.has_value()) {
            t = ray_segment_intersection(point, dir, triangle[i],
                                         triangle[i2]);
            continue;
        }

        // small distance could be preccission inconsistance
        std::optional<double> t2 = ray_segment_intersection(point, dir,
                                                            triangle[i],
                                                            triangle[i2]);
        if (t2.has_value() && *t2 > *t) t = t2;
    }
    assert(t.has_value()); // Not found intersection.
    return point + dir * (*t);
}

Vec3d calc_hit_point(const igl::Hit &          h,
                     const Vec3i &             triangle,
                     const std::vector<Vec3f> &vertices)
{
    double c1 = h.u;
    double c2 = h.v;
    double c0 = 1.0 - c1 - c2;
    Vec3d  v0 = vertices[triangle[0]].cast<double>();
    Vec3d  v1 = vertices[triangle[1]].cast<double>();
    Vec3d  v2 = vertices[triangle[2]].cast<double>();
    return v0 * c0 + v1 * c1 + v2 * c2;
}

Vec3d calc_hit_point(const igl::Hit &h, indexed_triangle_set &its)
{
    return calc_hit_point(h, its.indices[h.id], its.vertices);
}
} // namespace Private

std::string get_font_filepath() {
    std::string resource_dir = 
        std::string(TEST_DATA_DIR) + "/../../resources/";
    return resource_dir + "fonts/NotoSans-Regular.ttf";
}

#include "imgui/imstb_truetype.h"
TEST_CASE("Read glyph C shape from font, stb library calls ONLY", "[Emboss]") {
    std::string font_path = get_font_filepath();
    char  letter   = 'C';
    
    // Read  font file
    FILE *file = fopen(font_path.c_str(), "rb");
    REQUIRE(file != nullptr);
    // find size of file
    REQUIRE(fseek(file, 0L, SEEK_END) == 0);
    size_t size = ftell(file);
    REQUIRE(size != 0);
    rewind(file);
    std::vector<unsigned char> buffer(size);
    size_t count_loaded_bytes = fread((void *) &buffer.front(), 1, size, file);
    REQUIRE(count_loaded_bytes == size);

    // Use stb true type library
    int font_offset = stbtt_GetFontOffsetForIndex(buffer.data(), 0);
    REQUIRE(font_offset >= 0);
    stbtt_fontinfo font_info;
    REQUIRE(stbtt_InitFont(&font_info, buffer.data(), font_offset) != 0);    
    int unicode_letter = (int) letter;
    int glyph_index = stbtt_FindGlyphIndex(&font_info, unicode_letter);
    REQUIRE(glyph_index != 0);
    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&font_info, glyph_index, &vertices);
    CHECK(num_verts > 0);
}

#include <libslic3r/Utils.hpp>
TEST_CASE("Convert glyph % to model", "[Emboss]") 
{
    std::string font_path = get_font_filepath();
    char  letter   = '%';
    float flatness = 2.;

    auto font = Emboss::create_font_file(font_path.c_str());
    REQUIRE(font != nullptr);

    std::optional<Emboss::Glyph> glyph = Emboss::letter2glyph(*font, letter, flatness);
    REQUIRE(glyph.has_value());

    ExPolygons shape = glyph->shape;    
    REQUIRE(!shape.empty());

    float z_depth = 1.f;
    Emboss::ProjectZ projection(z_depth);
    indexed_triangle_set its = Emboss::polygons2model(shape, projection);

    CHECK(!its.indices.empty());    
}

TEST_CASE("Test hit point", "[AABBTreeIndirect]")
{
    indexed_triangle_set its;
    its.vertices = {
        Vec3f(1, 1, 1),
        Vec3f(2, 10, 2),
        Vec3f(10, 0, 2),
    };
    its.indices = {Vec3i(0, 2, 1)};
    auto tree   = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
        its.vertices, its.indices);

    Vec3d    ray_point(8, 1, 0);
    Vec3d    ray_dir(0, 0, 1);
    igl::Hit hit;
    AABBTreeIndirect::intersect_ray_first_hit(its.vertices, its.indices, tree,
                                              ray_point, ray_dir, hit);
    Vec3d hp = Private::calc_hit_point(hit, its);
    CHECK(abs(hp.x() - ray_point.x()) < .1);
    CHECK(abs(hp.y() - ray_point.y()) < .1);
}

TEST_CASE("ray segment intersection", "[MeshBoolean]")
{
    Vec2d r_point(1, 1);
    Vec2d r_dir(1, 0);

    // colinear
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 0), Vec2d(2, 0)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(0, 0)).has_value());

    // before ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 0), Vec2d(0, 2)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 2), Vec2d(0, 0)).has_value());

    // above ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 2), Vec2d(2, 3)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 3), Vec2d(2, 2)).has_value());

    // belove ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(2, -1)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, -1), Vec2d(2, 0)).has_value());

    // intersection at [2,1] distance 1
    auto t1 = Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(2, 2));
    REQUIRE(t1.has_value());
    auto t2 = Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 2), Vec2d(2, 0));
    REQUIRE(t2.has_value());

    CHECK(abs(*t1 - *t2) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("triangle intersection", "[]")
{
    Vec2d                point(1, 1);
    Vec2d                dir(-1, 0);
    std::array<Vec2d, 3> triangle = {Vec2d(0, 0), Vec2d(5, 0), Vec2d(0, 5)};
    Vec2d                i = Private::get_intersection(point, dir, triangle);
    CHECK(abs(i.x()) < std::numeric_limits<double>::epsilon());
    CHECK(abs(i.y() - 1.) < std::numeric_limits<double>::epsilon());
}

#ifndef __APPLE__
#include <string>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
// Check function Emboss::is_italic that exist some italic and some non-italic font.
TEST_CASE("Italic check", "[Emboss]") 
{  
    std::queue<std::string> dir_paths;
#ifdef _WIN32
    dir_paths.push("C:/Windows/Fonts");
#elif defined(__linux__)
    dir_paths.push("/usr/share/fonts");
//#elif defined(__APPLE__)
//    dir_paths.push("//System/Library/Fonts");
#endif
    bool exist_italic = false;
    bool exist_non_italic = false;
    while (!dir_paths.empty()) {
        std::string dir_path = dir_paths.front();
        dir_paths.pop();
        for (const auto &entry : fs::directory_iterator(dir_path)) {
            const fs::path &act_path = entry.path();
            if (entry.is_directory()) {
                dir_paths.push(act_path.u8string());
                continue;
            }
            std::string ext = act_path.extension().u8string();
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (ext != ".ttf") continue;
            std::string path_str = act_path.u8string();
            auto        font_opt = Emboss::create_font_file(path_str.c_str());
            if (font_opt == nullptr) continue;

            unsigned int collection_number = 0;
            if (Emboss::is_italic(*font_opt, collection_number))
                exist_italic = true;
            else
                exist_non_italic = true;

            if (exist_italic && exist_non_italic) break;
            //std::cout << ((Emboss::is_italic(*font_opt)) ? "[yes] " : "[no ] ") << entry.path() << std::endl;
        }
    }
    CHECK(exist_italic);
    CHECK(exist_non_italic);
}
#endif // not __APPLE__


#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian_converter.h>

// Referencing a glyph contour (an ExPolygon) plus a vertex base of the contour.
struct GlyphContour
{
    // Index of a glyph in a vector of glyphs.
    int32_t glyph{-1};
    // Index of an ExPolygon in ExPolygons of a glyph.
    int32_t expoly{-1};
    // Index of a contour in ExPolygon.
    // 0 - outer contour, >0 - hole
    int32_t contour{-1};
    // Base of the zero'th point of a contour in text mesh.
    // There are two vertices (front and rear) created for each contour,
    // thus there are 2x more vertices in text mesh than the number of contour points.
    int32_t vertex_base{-1};
};

struct GlyphID
{
    int32_t glyph_contour{-1};
    // vertex or edge ID, where edge ID is the index of the source point.
    // There are 4 consecutive indices generated for a single glyph edge:
    // 0th - 1st text edge (straight)
    // 1th - 1st text face
    // 2nd - 2nd text edge (diagonal)
    // 3th - 2nd text face
    int32_t idx{-1};

    GlyphID &operator++()
    {
        ++idx;
        return *this;
    }
};

namespace Slic3r::MeshBoolean::cgal2 {

namespace CGALProc   = CGAL::Polygon_mesh_processing;
namespace CGALParams = CGAL::Polygon_mesh_processing::parameters;

//    using EpecKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using EpicKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using _EpicMesh  = CGAL::Surface_mesh<EpicKernel::Point_3>;
//    using _EpecMesh = CGAL::Surface_mesh<EpecKernel::Point_3>;

using CGALMesh = _EpicMesh;

// Add an indexed triangle mesh to CGAL Surface_mesh.
// Store map of CGAL face to source face index into object_face_source_id.
void triangle_mesh_to_cgal(
    const std::vector<stl_vertex>                       &V,
    const std::vector<stl_triangle_vertex_indices>      &F,
    CGALMesh                                            &out,
    CGALMesh::Property_map<CGAL::SM_Face_index, int32_t> object_face_source_id)
{
    if (F.empty()) return;

    size_t vertices_count = V.size();
    size_t edges_count    = (F.size() * 3) / 2;
    size_t faces_count    = F.size();
    out.reserve(vertices_count, edges_count, faces_count);

    for (auto &v : V)
        out.add_vertex(typename CGALMesh::Point{v.x(), v.y(), v.z()});

    using VI = typename CGALMesh::Vertex_index;
    for (auto &f : F) {
        auto fid = out.add_face(VI(f(0)), VI(f(1)), VI(f(2)));
        object_face_source_id[fid] = int32_t(&f - &F.front());
    }
}

void glyph2model(
    const ExPolygons                                     &glyph,
    int32_t                                               glyph_id,
    const Slic3r::Emboss::IProject                       &projection,
    CGALMesh                                             &out,
    std::vector<GlyphContour>                            &glyph_contours,
    CGALMesh::Property_map<CGAL::SM_Edge_index, GlyphID> &glyph_id_edge,
    CGALMesh::Property_map<CGAL::SM_Face_index, GlyphID> &glyph_id_face)
{
    std::vector<CGALMesh::Vertex_index> indices;
    auto insert_contour = [&projection, &indices, &out, glyph_id,
                           &glyph_contours, &glyph_id_edge,
                           &glyph_id_face](const Polygon &polygon,
                                           int32_t iexpoly, int32_t id) {
        indices.clear();
        indices.reserve(polygon.points.size() * 2);
        size_t  num_vertices_old = out.number_of_vertices();
        GlyphID glid{int32_t(glyph_contours.size()), 0};
        glyph_contours.push_back(
            {glyph_id, iexpoly, id, int32_t(num_vertices_old)});
        for (const Point &p2 : polygon.points) {
            auto p  = projection.project(p2);
            auto vi = out.add_vertex(typename CGALMesh::Point{p.first.x(),
                                                              p.first.y(),
                                                              p.first.z()});
            assert((size_t) vi == indices.size() + num_vertices_old);
            indices.emplace_back(vi);
            vi = out.add_vertex(typename CGALMesh::Point{p.second.x(),
                                                         p.second.y(),
                                                         p.second.z()});
            assert((size_t) vi == indices.size() + num_vertices_old);
            indices.emplace_back(vi);
        }
        for (int32_t i = 0; i < int32_t(indices.size()); i += 2) {
            int32_t j         = (i + 2) % int32_t(indices.size());
            auto    find_edge = [&out](CGALMesh::Face_index   fi,
                                    CGALMesh::Vertex_index from,
                                    CGALMesh::Vertex_index to) {
                CGALMesh::Halfedge_index hi = out.halfedge(fi);
                for (; out.target(hi) != to; hi = out.next(hi))
                    ;
                assert(out.source(hi) == from);
                assert(out.target(hi) == to);
                return hi;
            };
            auto fi = out.add_face(indices[i], indices[i + 1], indices[j]);
            glyph_id_edge[out.edge(
                find_edge(fi, indices[i], indices[i + 1]))] = glid;
            glyph_id_face[fi]                               = ++glid;
            glyph_id_edge[out.edge(
                find_edge(fi, indices[i + 1], indices[j]))] = ++glid;
            glyph_id_face[out.add_face(indices[j], indices[i + 1],
                                       indices[j + 1])]     = ++glid;
            ++glid;
        }
    };

    size_t count_point = count_points(glyph);
    out.reserve(out.number_of_vertices() + 2 * count_point,
                out.number_of_edges() + 4 * count_point,
                out.number_of_faces() + 2 * count_point);

    for (const ExPolygon &expolygon : glyph) {
        int32_t idx_contour = &expolygon - &glyph.front();
        insert_contour(expolygon.contour, idx_contour, 0);
        for (const Polygon &hole : expolygon.holes)
            insert_contour(hole, idx_contour,
                           1 + (&hole - &expolygon.holes.front()));
    }
}
} // namespace Slic3r::MeshBoolean::cgal2

bool its_write_obj(const indexed_triangle_set &its,
                   const std::vector<Vec3f>   &color,
                   const char                 *file)
{
    Slic3r::CNumericLocalesSetter locales_setter;
    FILE                         *fp = fopen(file, "w");
    if (fp == nullptr) { return false; }

    for (size_t i = 0; i < its.vertices.size(); ++i)
        fprintf(fp, "v %f %f %f %f %f %f\n", its.vertices[i](0),
                its.vertices[i](1), its.vertices[i](2), color[i](0),
                color[i](1), color[i](2));
    for (size_t i = 0; i < its.indices.size(); ++i)
        fprintf(fp, "f %d %d %d\n", its.indices[i][0] + 1,
                its.indices[i][1] + 1, its.indices[i][2] + 1);
    fclose(fp);
    return true;
}

TEST_CASE("Emboss extrude cut", "[Emboss-Cut]")
{
    std::string font_path = get_font_filepath();
    char  letter   = '%';
    float flatness = 2.;

    auto font = Emboss::create_font(font_path.c_str());
    REQUIRE(font != nullptr);

    std::optional<Emboss::Glyph> glyph = Emboss::letter2glyph(*font, letter,
                                                              flatness);
    REQUIRE(glyph.has_value());

    ExPolygons shape = glyph->shape;
    REQUIRE(!shape.empty());

    float            z_depth = 50.f;
    Emboss::ProjectZ projection(z_depth);

#if 0
    indexed_triangle_set text = Emboss::polygons2model(shape, projection);
    BoundingBoxf3 bbox = bounding_box(text);

    CHECK(!text.indices.empty());
#endif

    auto cube = its_make_cube(782 - 49 + 50, 724 + 10 + 50, 5);
    its_translate(cube, Vec3f(49 - 25, -10 - 25, 2.5));
    auto cube2 = cube;
    //    its_translate(cube2, Vec3f(0, 0, 40));
    its_translate(cube2, Vec3f(0, -40, 40));
    for (auto &face : cube2.indices)
        for (int i = 0; i < 3; ++i) face(i) += int(cube.vertices.size());
    append(cube.vertices, cube2.vertices);
    append(cube.indices, cube2.indices);

    MeshBoolean::cgal2::CGALMesh cgalcube, cgaltext;
    auto                         object_face_source_id =
        cgalcube
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Face_index,
                              int32_t>("f:object_face_source_id")
            .first;
    MeshBoolean::cgal2::triangle_mesh_to_cgal(cube.vertices, cube.indices,
                                              cgalcube,
                                              object_face_source_id);

    auto edge_glyph_id =
        cgaltext
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Edge_index,
                              GlyphID>("e:glyph_id")
            .first;
    auto face_glyph_id =
        cgaltext
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Face_index,
                              GlyphID>("f:glyph_id")
            .first;
    auto vertex_glyph_id =
        cgalcube
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Vertex_index,
                              GlyphID>("v:glyph_id")
            .first;
    std::vector<GlyphContour> glyph_contours;

    MeshBoolean::cgal2::glyph2model(shape, 0, projection, cgaltext,
                                    glyph_contours, edge_glyph_id,
                                    face_glyph_id);

    struct Visitor
    {
        using TriangleMesh = Slic3r::MeshBoolean::cgal2::CGALMesh;

        const TriangleMesh &object;
        const TriangleMesh &glyphs;
        //        const std::vector<GlyphContour> &glyph_contours;
        // Properties of the glyphs mesh:
        TriangleMesh::Property_map<CGAL::SM_Edge_index, GlyphID> glyph_id_edge;
        TriangleMesh::Property_map<CGAL::SM_Face_index, GlyphID> glyph_id_face;
        // Properties of the object mesh.
        TriangleMesh::Property_map<CGAL::SM_Face_index, int32_t>
            object_face_source_id;
        TriangleMesh::Property_map<CGAL::SM_Vertex_index, GlyphID>
            object_vertex_glyph_id;

        typedef boost::graph_traits<TriangleMesh> GT;
        typedef typename GT::face_descriptor      face_descriptor;
        typedef typename GT::halfedge_descriptor  halfedge_descriptor;
        typedef typename GT::vertex_descriptor    vertex_descriptor;

        int32_t source_face_id;

        void before_subface_creations(face_descriptor f_old,
                                      TriangleMesh   &mesh)
        {
            assert(&mesh == &object);
            source_face_id = object_face_source_id[f_old];
        }
        void after_subface_created(face_descriptor f_new, TriangleMesh &mesh)
        {
            assert(&mesh == &object);
            object_face_source_id[f_new] = source_face_id;
        }

        std::vector<const GlyphID *> intersection_point_glyph;

        // Intersecting an edge hh_edge from tm_edge with a face hh_face of tm_face.
        void intersection_point_detected(
            // ID of the intersection point, starting at 0. Ids are consecutive.
            std::size_t i_id,
            // Dimension of a simplex part of face(hh_face) that is
            // intersected by hh_edge: 0 for vertex: target(hh_face) 1 for
            // edge: hh_face 2 for the interior of face: face(hh_face)
            int simplex_dimension,
            // Edge of tm_edge, see edge_source_coplanar_with_face &
            // edge_target_coplanar_with_face whether any vertex of hh_edge is
            // coplanar with face(hh_face).
            halfedge_descriptor hh_edge,
            // Vertex, halfedge or face of tm_face intersected by hh_edge, see
            // comment at simplex_dimension.
            halfedge_descriptor hh_face,
            // Mesh containing hh_edge
            const TriangleMesh &tm_edge,
            // Mesh containing hh_face
            const TriangleMesh &tm_face,
            // source(hh_edge) is coplanar with face(hh_face).
            bool edge_source_coplanar_with_face,
            // target(hh_edge) is coplanar with face(hh_face).
            bool edge_target_coplanar_with_face)
        {
            if (i_id <= intersection_point_glyph.size()) {
                intersection_point_glyph.reserve(
                    Slic3r::next_highest_power_of_2(i_id + 1));
                intersection_point_glyph.resize(i_id + 1);
            }

            const GlyphID *glyph = nullptr;
            if (&tm_face == &glyphs) {
                assert(&tm_edge == &object);
                switch (simplex_dimension) {
                case 1:
                    // edge x edge intersection
                    glyph = &glyph_id_edge[glyphs.edge(hh_face)];
                    break;
                case 2:
                    // edge x face intersection
                    glyph = &glyph_id_face[glyphs.face(hh_face)];
                    break;
                default: assert(false);
                }
                if (edge_source_coplanar_with_face)
                    object_vertex_glyph_id[object.source(hh_edge)] = *glyph;
                if (edge_target_coplanar_with_face)
                    object_vertex_glyph_id[object.target(hh_edge)] = *glyph;
            } else {
                assert(&tm_edge == &glyphs && &tm_face == &object);
                assert(!edge_source_coplanar_with_face);
                assert(!edge_target_coplanar_with_face);
                glyph = &glyph_id_edge[glyphs.edge(hh_edge)];
                if (simplex_dimension == 0)
                    object_vertex_glyph_id[object.target(hh_face)] = *glyph;
            }
            intersection_point_glyph[i_id] = glyph;
        }

        void new_vertex_added(std::size_t         node_id,
                              vertex_descriptor   vh,
                              const TriangleMesh &tm)
        {
            assert(&tm == &object);
            assert(node_id < intersection_point_glyph.size());
            const GlyphID *glyph = intersection_point_glyph[node_id];
            assert(glyph != nullptr);
            assert(glyph->glyph_contour != -1);
            assert(glyph->idx != -1);
            object_vertex_glyph_id[vh] = glyph ? *glyph : GlyphID{};
        }

        void after_subface_creations(TriangleMesh &) {}
        void before_subface_created(TriangleMesh &) {}
        void before_edge_split(halfedge_descriptor /* h */,
                               TriangleMesh & /* tm */)
        {}
        void edge_split(halfedge_descriptor /* hnew */,
                        TriangleMesh & /* tm */)
        {}
        void after_edge_split() {}
        void add_retriangulation_edge(halfedge_descriptor /* h */,
                                      TriangleMesh & /* tm */)
        {}
    } visitor{cgalcube,
              cgaltext,
              /* glyph_contours, */ edge_glyph_id,
              face_glyph_id,
              object_face_source_id,
              vertex_glyph_id};

    auto        ecm = get(CGAL::dynamic_edge_property_t<bool>(), cgalcube);
    const auto &p =
        CGAL::Polygon_mesh_processing::parameters::throw_on_self_intersection(
            false)
            .visitor(visitor)
            .edge_is_constrained_map(ecm);
    const auto &q = CGAL::Polygon_mesh_processing::parameters::visitor(visitor)
                        .do_not_modify(true);
    //    CGAL::Polygon_mesh_processing::corefine(cgalcube, cgalcube2, p, p);

    CGAL::Polygon_mesh_processing::corefine(cgalcube, cgaltext, p, q);

    auto vertex_colors =
        cgalcube
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Vertex_index,
                              CGAL::Color>("v:color")
            .first;
    auto face_colors =
        cgalcube
            .add_property_map<MeshBoolean::cgal2::CGALMesh::Face_index,
                              CGAL::Color>("f:color")
            .first;

    const CGAL::Color marked{255, 0, 0};
    for (auto fi : cgalcube.faces()) {
        CGAL::Color color(0, 255, 0);
        auto        hi_end = cgalcube.halfedge(fi);
        auto        hi     = hi_end;
        do {
            if (get(ecm, cgalcube.edge(hi))) {
                // This face has a constrained edge.
                GlyphID g1 = vertex_glyph_id[cgalcube.source(hi)];
                GlyphID g2 = vertex_glyph_id[cgalcube.target(hi)];
                assert(g1.glyph_contour != -1 &&
                       g1.glyph_contour == g2.glyph_contour);
                assert(g1.idx != -1);
                assert(g2.idx != -1);
                const GlyphContour &glyph_contour =
                    glyph_contours[g1.glyph_contour];
                const auto &expoly = glyph->shape[glyph_contour.expoly];
                const auto &contour =
                    glyph_contour.contour == 0 ?
                        expoly.contour :
                        expoly.holes[glyph_contour.contour - 1];
                bool    inside = false;
                int32_t i1     = g1.idx / 4;
                int32_t i2     = g2.idx / 4;
                if (g1.idx == g2.idx) {
                    // Crossing both object vertices with the same glyph face.
                    int type = g1.idx % 4;
                    assert(type == 1 || type == 3);
                    const auto &p = cgalcube.point(
                        cgalcube.target(cgalcube.next(hi)));
                    int i = i1 * 2;
                    int j = (i1 + 1 == int(contour.size())) ? 0 : i + 2;
                    i += glyph_contour.vertex_base;
                    j += glyph_contour.vertex_base;
                    auto abcp =
                        type == 1 ?
                            CGAL::orientation(
                                cgaltext.point(CGAL::SM_Vertex_index(i)),
                                cgaltext.point(CGAL::SM_Vertex_index(i + 1)),
                                cgaltext.point(CGAL::SM_Vertex_index(j)), p) :
                            CGAL::orientation(
                                cgaltext.point(CGAL::SM_Vertex_index(j)),
                                cgaltext.point(CGAL::SM_Vertex_index(i + 1)),
                                cgaltext.point(CGAL::SM_Vertex_index(j + 1)),
                                p);
                    inside = abcp == CGAL::POSITIVE;
                } else if (g1.idx < g2.idx) {
                    if (i1 == 0 && i2 + 1 == contour.size()) {
                        // cw
                    } else {
                        inside = true;
                    }
                } else {
                    if (i2 == 0 && i1 + 1 == contour.size()) {
                        inside = true;
                        std::swap(g1, g2);
                        std::swap(i1, i2);
                    }
                }
                if (inside) {
                    // Is this face oriented towards p or away from p?
                    const auto &a = cgalcube.point(cgalcube.source(hi));
                    const auto &b = cgalcube.point(cgalcube.target(hi));
                    const auto &c = cgalcube.point(
                        cgalcube.target(cgalcube.next(hi)));
                    // FIXME prosim nahrad skutecnou projekci.
                    // projection.project()
                    const auto p =
                        a +
                        MeshBoolean::cgal2::EpicKernel::Vector_3(0, 0, 10);
                    auto abcp = CGAL::orientation(a, b, c, p);
                    if (abcp == CGAL::POSITIVE) color = marked;
                }
                break;
            }
            hi = cgalcube.next(hi);
        } while (hi != hi_end);
        face_colors[fi] = color;
    }

    CGAL::IO::write_OFF("c:\\data\\temp\\corefined-0.off", cgalcube);

    // Seed fill the other faces inside the region.
    std::vector<MeshBoolean::cgal2::CGALMesh::Face_index> queue;
    for (auto fi_seed : cgalcube.faces())
        if (face_colors[fi_seed] != marked) {
            // Is this face completely unconstrained?
            auto hi      = cgalcube.halfedge(fi_seed);
            auto hi_prev = cgalcube.prev(hi);
            auto hi_next = cgalcube.next(hi);
            if (!get(ecm, cgalcube.edge(hi)) &&
                !get(ecm, cgalcube.edge(hi_prev)) &&
                !get(ecm, cgalcube.edge(hi_next))) {
                queue.emplace_back(fi_seed);
                do {
                    auto fi = queue.back();
                    queue.pop_back();
                    auto hi      = cgalcube.halfedge(fi);
                    auto hi_prev = cgalcube.prev(hi);
                    auto hi_next = cgalcube.next(hi);
                    // The following condition may not apply if crossing a
                    // silhouette wrt. the glyph projection direction.
                    //                    assert(! get(ecm, cgalcube.edge(hi))
                    //                    && ! get(ecm, cgalcube.edge(hi_prev))
                    //                    && ! get(ecm, cgalcube.edge(hi_next)));
                    auto this_opposite = cgalcube.face(cgalcube.opposite(hi));
                    bool this_marked   = face_colors[this_opposite] == marked;
                    auto prev_opposite = cgalcube.face(
                        cgalcube.opposite(hi_prev));
                    bool prev_marked   = face_colors[prev_opposite] == marked;
                    auto next_opposite = cgalcube.face(
                        cgalcube.opposite(hi_next));
                    bool next_marked = face_colors[next_opposite] == marked;
                    int  num_marked = this_marked + prev_marked + next_marked;
                    if (num_marked >= 2) {
                        face_colors[fi] = marked;
                        if (num_marked == 2)
                            queue.emplace_back(!this_marked ? this_opposite :
                                               !prev_marked ? prev_opposite :
                                                              next_opposite);
                    }
                } while (!queue.empty());
            }
        }

    CGAL::IO::write_OFF("c:\\data\\temp\\corefined.off", cgalcube);

    // Mapping of its_extruded faces to source faces.
    enum class FaceState : int8_t {
        Unknown         = -1,
        Unmarked        = -2,
        UnmarkedSplit   = -3,
        Marked          = -4,
        MarkedSplit     = -5,
        UnmarkedEmitted = -6,
    };
    std::vector<FaceState> face_states(cube.indices.size(),
                                       FaceState::Unknown);
    for (auto fi_seed : cgalcube.faces()) {
        FaceState &state = face_states[object_face_source_id[fi_seed]];
        bool       m     = face_colors[fi_seed] == marked;
        switch (state) {
        case FaceState::Unknown:
            state = m ? FaceState::Marked : FaceState::Unmarked;
            break;
        case FaceState::Unmarked:
        case FaceState::UnmarkedSplit:
            state = m ? FaceState::MarkedSplit : FaceState::UnmarkedSplit;
            break;
        case FaceState::Marked:
        case FaceState::MarkedSplit: state = FaceState::MarkedSplit; break;
        default: assert(false);
        }
    }

    indexed_triangle_set its_extruded;
    its_extruded.indices.reserve(cgalcube.number_of_faces());
    its_extruded.vertices.reserve(cgalcube.number_of_vertices());
    // Mapping of its_extruded vertices (original and offsetted) to
    // cgalcuble's vertices.
    std::vector<std::pair<int32_t, int32_t>>
        map_vertices(cgalcube.number_of_vertices(),
                     std::pair<int32_t, int32_t>{-1, -1});

    Vec3f extrude_dir{0, 0, 5.f};
    for (auto fi : cgalcube.faces()) {
        const int32_t   source_face_id = object_face_source_id[fi];
        const FaceState state          = face_states[source_face_id];
        assert(state == FaceState::Unmarked ||
               state == FaceState::UnmarkedSplit ||
               state == FaceState::UnmarkedEmitted ||
               state == FaceState::Marked || state == FaceState::MarkedSplit);
        if (state == FaceState::UnmarkedEmitted) {
            // Already emitted.
        } else if (state == FaceState::Unmarked ||
                   state == FaceState::UnmarkedSplit) {
            // Just copy the unsplit source face.
            const Vec3i source_vertices = cube.indices[source_face_id];
            Vec3i       target_vertices;
            for (int i = 0; i < 3; ++i) {
                target_vertices(i) = map_vertices[source_vertices(i)].first;
                if (target_vertices(i) == -1) {
                    map_vertices[source_vertices(i)].first = target_vertices(
                        i) = int(its_extruded.vertices.size());
                    its_extruded.vertices.emplace_back(
                        cube.vertices[source_vertices(i)]);
                }
            }
            its_extruded.indices.emplace_back(target_vertices);
            face_states[source_face_id] = FaceState::UnmarkedEmitted;
        } else {
            auto hi      = cgalcube.halfedge(fi);
            auto hi_prev = cgalcube.prev(hi);
            auto hi_next = cgalcube.next(hi);
            const Vec3i
                  source_vertices{int((std::size_t) cgalcube.target(hi)),
                                int((std::size_t) cgalcube.target(hi_next)),
                                int((std::size_t) cgalcube.target(hi_prev))};
            Vec3i target_vertices;
            if (face_colors[fi] == marked) {
                // Extrude the face. Neighbor edges separating extruded face
                // from non-extruded face will be extruded.
                bool  boundary_vertex[3] = {false, false, false};
                Vec3i target_vertices_extruded{-1, -1, -1};
                for (int i = 0; i < 3; ++i) {
                    if (face_colors[cgalcube.face(cgalcube.opposite(hi))] !=
                        marked)
                        // Edge separating extruded / non-extruded region.
                        boundary_vertex[i] = boundary_vertex[(i + 2) % 3] =
                            true;
                    hi = cgalcube.next(hi);
                }
                for (int i = 0; i < 3; ++i) {
                    target_vertices_extruded(
                        i) = map_vertices[source_vertices(i)].second;
                    if (target_vertices_extruded(i) == -1) {
                        map_vertices[source_vertices(i)].second =
                            target_vertices_extruded(i) = int(
                                its_extruded.vertices.size());
                        const auto &p = cgalcube.point(cgalcube.target(hi));
                        its_extruded.vertices.emplace_back(
                            Vec3f{float(p.x()), float(p.y()), float(p.z())} +
                            extrude_dir);
                    }
                    if (boundary_vertex[i]) {
                        target_vertices(
                            i) = map_vertices[source_vertices(i)].first;
                        if (target_vertices(i) == -1) {
                            map_vertices[source_vertices(i)].first =
                                target_vertices(i) = int(
                                    its_extruded.vertices.size());
                            const auto &p = cgalcube.point(
                                cgalcube.target(hi));
                            its_extruded.vertices.emplace_back(p.x(), p.y(),
                                                               p.z());
                        }
                    }
                    hi = cgalcube.next(hi);
                }
                its_extruded.indices.emplace_back(target_vertices_extruded);
                // Add the sides.
                for (int i = 0; i < 3; ++i) {
                    int j = (i + 1) % 3;
                    assert(target_vertices_extruded[i] != -1 &&
                           target_vertices_extruded[j] != -1);
                    if (boundary_vertex[i] && boundary_vertex[j]) {
                        assert(target_vertices[i] != -1 &&
                               target_vertices[j] != -1);
                        its_extruded.indices.emplace_back(
                            Vec3i{target_vertices[i], target_vertices[j],
                                  target_vertices_extruded[i]});
                        its_extruded.indices.emplace_back(
                            Vec3i{target_vertices_extruded[i],
                                  target_vertices[j],
                                  target_vertices_extruded[j]});
                    }
                }
            } else {
                // Copy the face.
                Vec3i target_vertices;
                for (int i = 0; i < 3; ++i) {
                    target_vertices(
                        i) = map_vertices[source_vertices(i)].first;
                    if (target_vertices(i) == -1) {
                        map_vertices[source_vertices(i)].first =
                            target_vertices(i) = int(
                                its_extruded.vertices.size());
                        const auto &p = cgalcube.point(cgalcube.target(hi));
                        its_extruded.vertices.emplace_back(p.x(), p.y(),
                                                           p.z());
                    }
                    hi = cgalcube.next(hi);
                }
                its_extruded.indices.emplace_back(target_vertices);
            }
        }
    }

    its_write_obj(its_extruded, "c:\\data\\temp\\text-extruded.obj");

    indexed_triangle_set edges_its;
    std::vector<Vec3f>   edges_its_colors;
    for (auto ei : cgalcube.edges())
        if (cgalcube.is_valid(ei)) {
            const auto &p1          = cgalcube.point(cgalcube.vertex(ei, 0));
            const auto &p2          = cgalcube.point(cgalcube.vertex(ei, 1));
            bool        constrained = get(ecm, ei);
            Vec3f color = constrained ? Vec3f{1.f, 0, 0} : Vec3f{0, 1., 0};
            edges_its.indices.emplace_back(
                Vec3i(edges_its.vertices.size(), edges_its.vertices.size() + 1,
                      edges_its.vertices.size() + 2));
            edges_its.vertices.emplace_back(Vec3f(p1.x(), p1.y(), p1.z()));
            edges_its.vertices.emplace_back(Vec3f(p2.x(), p2.y(), p2.z()));
            edges_its.vertices.emplace_back(
                Vec3f(p2.x(), p2.y(), p2.z() + 0.001));
            edges_its_colors.emplace_back(color);
            edges_its_colors.emplace_back(color);
            edges_its_colors.emplace_back(color);
        }
    its_write_obj(edges_its, edges_its_colors,
                  "c:\\data\\temp\\corefined-edges.obj");

    //    MeshBoolean::cgal::minus(cube, cube2);

    //    REQUIRE(!MeshBoolean::cgal::does_self_intersect(cube));
}
