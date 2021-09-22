#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/MeshBoolean.hpp>
#include <libslic3r/SimplifyMesh.hpp>

using namespace Slic3r;

TEST_CASE("CGAL and TriangleMesh conversions", "[MeshBoolean]") {
    TriangleMesh sphere = make_sphere(1.);
    
    auto cgalmesh_ptr = MeshBoolean::cgal::triangle_mesh_to_cgal(sphere);
    
    REQUIRE(cgalmesh_ptr);
    REQUIRE(! MeshBoolean::cgal::does_self_intersect(*cgalmesh_ptr));
    
    TriangleMesh M = MeshBoolean::cgal::cgal_to_triangle_mesh(*cgalmesh_ptr);
    
    REQUIRE(M.its.vertices.size() == sphere.its.vertices.size());
    REQUIRE(M.its.indices.size() == sphere.its.indices.size());
    
    REQUIRE(M.volume() == Approx(sphere.volume()));
    
    REQUIRE(! MeshBoolean::cgal::does_self_intersect(M));
}

Vec3d calc_normal(const Vec3i &triangle, const std::vector<Vec3f> &vertices)
{
    Vec3d v0 = vertices[triangle[0]].cast<double>();
    Vec3d v1 = vertices[triangle[1]].cast<double>();
    Vec3d v2 = vertices[triangle[2]].cast<double>();
    // n = triangle normal
    Vec3d n = (v1 - v0).cross(v2 - v0);
    n.normalize();
    return n;
}

TEST_CASE("Add TriangleMeshes", "[MeshBoolean]")
{
    TriangleMesh tm1 = make_sphere(1.6, 1.6);
    Vec3f        move(5, -3, 7);
    move.normalize();
    tm1.translate(0.3 * move);
    its_write_obj(tm1.its, "tm1.obj");
    TriangleMesh tm2 = make_cube(1., 1., 1.);
    its_write_obj(tm2.its, "tm2.obj");
    MeshBoolean::cgal::plus(tm1, tm2);
    its_write_obj(tm1.its, "test_add.obj");
}

#include "libslic3r/Emboss.hpp"

ExPolygons ttf2polygons(const char * font_name, char letter, float flatness = 1.f) {
    auto font = Emboss::load_font(font_name);
    if (!font.has_value()) return ExPolygons();
    return Emboss::letter2glyph(*font, letter, flatness)->shape;
}

#include "libslic3r/SVG.hpp"
void store_to_svg(Polygons polygons,std::string file_name = "letter.svg")
{
    double scale = 1e6;
    BoundingBox bb;
    for (auto& p : polygons) {
        p.scale(scale);
        bb.merge(p.points);
    }
    SVG svg(file_name, bb);
    svg.draw(polygons);    
}

struct Plane
{
    Vec3d point = Vec3d(0., 0., 0.); // lay on plane - define zero position
    Vec3d normal = Vec3d(0., 0., 1.);// [unit vector] - define orientation
};

struct OrientedPlane: public Plane
{
    // Must be perpendiculat to normal
    Vec3d up = Vec3d(0., 1., 0.);    // [unit vector]
};

struct EmbossConfig
{
    // emboss plane must be above surface
    // point define zero for 2d polygon and must be out of model
    // normal is direction to model
    // up define orientation of polygon
    OrientedPlane projection;
    

    // Move surface distance
    // Positive value out of model (Raised)
    // Negative value into model (Engraved) 
    float height = 1.; // [in milimeters]
};

#include <optional>
#include <libslic3r/AABBTreeIndirect.hpp> 
Vec3d calc_hit_point(const igl::Hit &          h,
                     const Vec3i &             triangle,
                     const std::vector<Vec3f> &vertices)
{
    double c1  = h.u;
    double c2  = h.v;
    double c0  = 1.0 - c1 - c2;
    Vec3d v0 = vertices[triangle[0]].cast<double>();
    Vec3d v1 = vertices[triangle[1]].cast<double>();
    Vec3d v2 = vertices[triangle[2]].cast<double>();
    return v0 * c0 + v1 * c1 + v2 * c2;
}

Vec3d calc_hit_point(const igl::Hit &h, indexed_triangle_set &its)
{
    return calc_hit_point(h, its.indices[h.id], its.vertices);
}

TEST_CASE("Test hit point", "[AABBTreeIndirect]") 
{ 
    indexed_triangle_set its;
    its.vertices = {
        Vec3f(1,1,1),
        Vec3f(2, 10, 2),
        Vec3f(10, 0, 2),
    };
    its.indices = {Vec3i(0, 2, 1)};
    auto tree   = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
        its.vertices, its.indices);

    Vec3d ray_point(8, 1, 0);
    Vec3d ray_dir(0,0,1);
    igl::Hit hit;
    AABBTreeIndirect::intersect_ray_first_hit(its.vertices, its.indices, tree,
                                              ray_point, ray_dir, hit);
    Vec3d hp = calc_hit_point(hit, its);
    CHECK(abs(hp.x() - ray_point.x()) < .1);
    CHECK(abs(hp.y() - ray_point.y()) < .1);
}

// represents triangle extend by seam
struct TrianglePath
{
    // input edge, output edge, when cross triangle border
    std::optional<std::pair<char, char>> edges;

    // when edges has value than first and last point lay on triangle edge
    std::vector<Vec3f> points;

    // first point has index offset_id in result vertices
    uint32_t offset_id;
};
using TrianglePaths = std::vector<TrianglePath>;

// create transformation matrix to convert direction vectors
// do not care about up vector
// Directions are normalized
Eigen::Matrix3d create_transformation(const Vec3d &from_dir, const Vec3d &to_dir)
{
    Vec3d  axis     = from_dir.cross(to_dir);
    axis.normalize();
    double angle    = acos(from_dir.dot(to_dir));
    auto   rotation = Eigen::AngleAxisd(angle, axis);
    return rotation.matrix();
}

TEST_CASE("Transformation matrix", "[]") { 
    Vec3d d1(3, -7, 13);
    Vec3d d2(-9, 5, 1);
    d1.normalize();
    d2.normalize();
    auto tr_mat = create_transformation(d1, d2);

    Vec3d  d1_tr = tr_mat * d1;
    Vec3d  diff  = d1_tr - d2;
    double eps = 1e-12;
    for (double d : diff) 
        CHECK(abs(d) < std::numeric_limits<double>::epsilon());
}

// calculate multiplication of ray dir to intersect - inspired by segment_segment_intersection
// when ray dir is normalized retur distance from ray point to intersection
// No value mean no intersection
std::optional<double> ray_segment_intersection(
    const Vec2d& r_point, const Vec2d& r_dir, const Vec2d& s0, const Vec2d& s1){

    auto denominate = [](const Vec2d& v0,const Vec2d& v1)->double{
        return v0.x() * v1.y() - v1.x() * v0.y();
    };

    Vec2d segment_dir = s1 - s0;
    double d = denominate(segment_dir, r_dir);
    if (std::abs(d) < std::numeric_limits<double>::epsilon())
        // Line and ray are collinear.
        return {};

    Vec2d s12 = s0 - r_point;
    double s_number = denominate(r_dir, s12);
    bool change_sign = false;
    if (d < 0.) {
        change_sign = true;
        d           = -d;
        s_number    = -s_number;
    }

    if (s_number < 0.|| s_number > d)
        // Intersection outside of segment.
        return {};

    double r_number = denominate(segment_dir, s12);
    if (change_sign)
        r_number = - r_number;

    if(r_number < 0.)
        // Intersection before ray start.
        return {};

    return r_number / d;
}

TEST_CASE("ray segment intersection", "[MeshBoolean]")
{   
    Vec2d r_point(1,1);
    Vec2d r_dir(1,0);
    
    // colinear
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(0,0), Vec2d(2,0)).has_value());
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(2,0), Vec2d(0,0)).has_value());

    // before ray
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(0,0), Vec2d(0,2)).has_value());
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(0,2), Vec2d(0,0)).has_value());

    // above ray
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(2,2), Vec2d(2,3)).has_value());
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(2,3), Vec2d(2,2)).has_value());

    // belove ray
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(2,0), Vec2d(2,-1)).has_value());
    CHECK(!ray_segment_intersection(r_point, r_dir, Vec2d(2,-1), Vec2d(2,0)).has_value());

    // intersection at [2,1] distance 1
    auto t1 = ray_segment_intersection(r_point, r_dir, Vec2d(2,0), Vec2d(2,2));
    REQUIRE(t1.has_value());
    auto t2 = ray_segment_intersection(r_point, r_dir, Vec2d(2,2), Vec2d(2,0));
    REQUIRE(t2.has_value());

    CHECK(abs(*t1 - *t2) < std::numeric_limits<double>::epsilon());
}

Vec2d get_intersection(const Vec2d& point, const Vec2d& dir, const std::array<Vec2d, 3>& triangle)
{    
    std::optional<double> t;
    for (size_t i = 0; i < 3; ++i) {
        size_t i2 = i+1;
        if(i2 == 3) i2 = 0;
        if (!t.has_value()) {
            t = ray_segment_intersection(point, dir, triangle[i], triangle[i2]);
            continue;
        }

        // small distance could be preccission inconsistance
        std::optional<double> t2 = ray_segment_intersection(
            point, dir, triangle[i], triangle[i2]);
        if (t2.has_value() && *t2 > *t) t = t2;
        
    }
    assert(t.has_value()); //Not found intersection.
    return point + dir * (*t);
}

TEST_CASE("triangle intersection", "[]")
{
    Vec2d point(1,1);
    Vec2d dir(-1, 0);
    std::array<Vec2d, 3> triangle = {
        Vec2d(0, 0),
        Vec2d(5, 0), 
        Vec2d(0, 5)
    };    
    Vec2d i = get_intersection(point, dir, triangle);
    CHECK(abs(i.x()) < std::numeric_limits<double>::epsilon());
    CHECK(abs(i.y() - 1.) < std::numeric_limits<double>::epsilon());
}

#include <libslic3r/Geometry.hpp>
#include <libslic3r/Triangulation.hpp>

indexed_triangle_set emboss3d(const Polygon &shape, float height) {
    // CW order of triangle indices
    std::vector<Vec3i> shape_triangles = Triangulation::triangulate(shape);

    indexed_triangle_set result;
    const Points &pts = shape.points;
    size_t count_point = pts.size();
    result.vertices.reserve(2 * count_point);
    // top points
    std::transform(pts.begin(), pts.end(), std::back_inserter(result.vertices),
                   [](const Point &p) { return Vec3f(p.x(), p.y(), 0); });
    // bottom points
    std::transform(pts.begin(), pts.end(), std::back_inserter(result.vertices),
                   [height](const Point &p) { return Vec3f(p.x(), p.y(), height); });
       
    result.indices.reserve(shape_triangles.size() * 2 + count_point*2);
    // top triangles - change to CCW
    for (const Vec3i &t : shape_triangles) 
        result.indices.emplace_back(t.x(), t.z(), t.y());
    // bottom triangles - use CW
    for (const Vec3i &t : shape_triangles)
        result.indices.emplace_back(
            t.x() + count_point,            
            t.y() + count_point,
            t.z() + count_point
        );

    // quads around - zig zag by triangles
    for (uint32_t i = 0; i < count_point; ++i) {
        // previous index
        uint32_t ip = (i == 0) ? (count_point - 1) : (i - 1);
        // bottom indices
        uint32_t i2 = i + count_point;
        uint32_t ip2 = ip + count_point;

        result.indices.emplace_back(i, i2, ip);
        result.indices.emplace_back(ip2, ip, i2);
    }
    return result;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel   EK;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
namespace PMP    = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

// is it realy neccessary to copy all?
Mesh its_to_mesh(const indexed_triangle_set &its)
{
    Mesh mesh;
    size_t vertices_count = its.vertices.size();
    size_t edges_count    = (its.indices.size() * 3) / 2;
    size_t faces_count    = its.indices.size();
    mesh.reserve(vertices_count, edges_count, faces_count);

    for (auto &v : its.vertices)
        mesh.add_vertex(typename Mesh::Point{v.x(), v.y(), v.z()});

    using VI = typename Mesh::Vertex_index;
    for (auto &f : its.indices)
        mesh.add_face(VI(f(0)), VI(f(1)), VI(f(2)));

    return mesh;
}

indexed_triangle_set mesh_to_its(const Mesh &mesh)
{
    indexed_triangle_set res;
    res.vertices.reserve(mesh.num_vertices());
    res.indices.reserve(mesh.num_faces());

    for (auto &vi : mesh.vertices()) {
        auto &v = mesh.point(vi); // index addresing, to compare same float point
        res.vertices.emplace_back(v.x(),v.y(), v.z());
    }

    for (auto &face : mesh.faces()) {
        auto vtc = mesh.vertices_around_face(mesh.halfedge(face));
        int   i   = 0;
        Vec3i facet;
        for (auto v : vtc) {
            if (i > 2) {
                i = 0;
                break;
            }
            facet(i++) = v;
        }
        if (i == 3) res.indices.emplace_back(facet);
    }
    return res;
}

void emboss3d_(const Polygon& shape, const EmbossConfig &cfg, indexed_triangle_set &its)
{
    indexed_triangle_set shape3d = emboss3d(shape, cfg.height);

    //its_write_obj(shape3d,"shape3d.obj");

    Mesh m1 = its_to_mesh(its);
    Mesh m2 = its_to_mesh(shape3d);

    size_t in_vertices_count = its.vertices.size();
    size_t in_indices_count = its.indices.size();

    auto tt = m1.property_map<Mesh::Edge_index, bool>("e:removed");
    Mesh::Property_map<Mesh::Edge_index, bool> ecm1 = tt.first; // hope in copy

    PMP::corefine(m1, m2
        , PMP::parameters::edge_is_constrained_map(ecm1)
        //, PMP::parameters::do_not_modify(true)
    );
    //PMP::Corefinement::Intersection_of_triangle_meshes()

    size_t count_true = 0;
    for (const Mesh::Edge_index &e : edges(m1)) 
        if (ecm1[e]) {
            ++count_true;
            const Mesh::Halfedge_index &h = e.halfedge();
            const Mesh::Halfedge_index &h2 = m1.opposite(h);
            const Mesh::Vertex_index &v = m1.target(h);
            std::cout << 
                "edge " << e << 
                " halfedge " << h <<
                " target " << v <<
                " index0 " << v.idx() <<
                " index1 " << m1.target(h2).idx() <<
                std::endl;
        }
    
    
    its = mesh_to_its(m1);

    // twice and move additional vertices 
    Vec3f move = (cfg.projection.normal * cfg.height).cast<float>();
    size_t new_add_vertice = its.vertices.size() - in_vertices_count;
    its.vertices.reserve(its.vertices.size() + new_add_vertice);
    for (size_t i = in_vertices_count; i < its.vertices.size(); ++i)
        its.vertices.emplace_back(its.vertices[i] + move);            
    
    // zig zag edge collected edge
    for (size_t i = in_indices_count; i < its.indices.size(); ++i) {
        // new added triangle
        Vec3i &t = its.indices[i];

        std::array<bool, 3> is_new;
        bool is_inner = true;
        for (size_t i = 0; i < 3; ++i) {
            bool is_n = t[0] >= in_vertices_count;
            is_inner |= is_n;
            is_new[i] = is_n;
        }

    }

    its_write_obj(mesh_to_its(m1), "corefine1.obj");
    its_write_obj(mesh_to_its(m2), "corefine2.obj");

    /*MeshBoolean::cgal::CGALMesh cm1 = m1->m;

    My_visitor<MeshBoolean::cgal::CGALMesh> sm_v;

    CGAL::Polygon_mesh_processing::corefine(m1->m, m2->m, 
        CGAL::Polygon_mesh_processing::parameters::visitor(sm_v)
        , CGAL::parameters::do_not_modify(true)
        );*/

    //TriangleMesh tm1 = cgal_to_triangle_mesh(*m1);
    //store_obj("tm1.obj", &tm1);
    //TriangleMesh tm2 = cgal_to_triangle_mesh(*m2);
    //store_obj("tm2.obj", &tm1);

    its = mesh_to_its(m1);
    return;
}

void emboss3d(const Polygon& shape, const EmbossConfig &cfg, indexed_triangle_set &its)
{
    indexed_triangle_set shape3d = emboss3d(shape, abs(cfg.height));
    //its_write_obj(shape3d,"shape3d.obj");
    auto m1 = MeshBoolean::cgal::triangle_mesh_to_cgal(its);
    auto m2 = MeshBoolean::cgal::triangle_mesh_to_cgal(shape3d);

    bool is_plus = shape.is_counter_clockwise() == (cfg.height > 0.f);
    if (is_plus) {
        MeshBoolean::cgal::plus(*m1, *m2); 
    } else {
        MeshBoolean::cgal::minus(*m1, *m2); 
    }

    TriangleMesh tm1 = cgal_to_triangle_mesh(*m1);
    store_obj("tm1.obj", &tm1);
    //TriangleMesh tm2 = cgal_to_triangle_mesh(*m2);
    //store_obj("tm2.obj", &tm1);

    its = tm1.its; // copy
}

TEST_CASE("Emboss polygon example", "[MeshBoolean]")
{
    const char *font_name = "C:/windows/fonts/arialbd.ttf";
    char        letter    = '%';
    float       flatness  = 2.;
    ExPolygons  espolygons  = ttf2polygons(font_name, letter, flatness);

    //TriangleMesh tm = make_sphere(1., 1.); tm.scale(10.f);
    
    TriangleMesh tm = make_cube(10., 5., 2.);
    tm.translate(Vec3f(0, 0, 1.7));
    Polygons polygons;
    polygons        = {Polygon({{1, 1}, {1, 2}, {2, 2}, {2, 1}})}; // rectangle CW
    polygons        = {Polygon({{1, 1}, {2, 1}, {2, 2}, {1, 2}})}; // rectangle CCW

    EmbossConfig         ec;
    ec.height                = 3;
    indexed_triangle_set its = tm.its; // copy
    for (const auto& polygon : polygons) {
        // TODO: differ CW and CCW
        emboss3d_(polygon, ec, its);
    }
}

TEST_CASE("Triangulate by cgal", "[Triangulation]")
{
    Points points = {Point(1, 1), Point(2, 1), Point(2, 2), Point(1, 2)};
    Triangulation::HalfEdges edges1 = {{1, 3}};
    std::vector<Vec3i> indices1 = Triangulation::triangulate(points, edges1);

    auto check = [](int i1, int i2, Vec3i t)->bool { return true;
        return (t[0] == i1 || t[1] == i1 || t[2] == i1) &&
               (t[0] == i2 || t[1] == i2 || t[2] == i2);
    };
    REQUIRE(indices1.size() == 2);
    int i1 = edges1.begin()->first, 
        i2 = edges1.begin()->second;
    CHECK(check(i1, i2, indices1[0]));
    CHECK(check(i1, i2, indices1[1]));

    Triangulation::HalfEdges edges2 = {{0, 2}};
    std::vector<Vec3i> indices2 = Triangulation::triangulate(points, edges2);
    REQUIRE(indices2.size() == 2);
    i1 = edges2.begin()->first;
    i2 = edges2.begin()->second;
    CHECK(check(i1, i2, indices2[0]));
    CHECK(check(i1, i2, indices2[1]));
}
