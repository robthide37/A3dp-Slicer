#include "CutSurface.hpp"

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian_converter.h>

// libslic3r
#include "TriangleMesh.hpp" // its_merge
#include "Utils.hpp" // next_highest_power_of_2

namespace priv {

namespace CGALProc   = CGAL::Polygon_mesh_processing;
namespace CGALParams = CGAL::Polygon_mesh_processing::parameters;

using EpicKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CutMesh = CGAL::Surface_mesh<EpicKernel::Point_3>;
//    using EpecKernel = CGAL::Exact_predicates_exact_constructions_kernel;
//    using CutMesh = CGAL::Surface_mesh<EpecKernel::Point_3>;

using DynamicEdgeProperty = CGAL::dynamic_edge_property_t<bool>;
using SMPM = boost::property_map<priv::CutMesh, DynamicEdgeProperty>::SMPM; 
using EcmType = CGAL::internal::Dynamic<priv::CutMesh, SMPM>;

/// <summary>
/// IntersectingElement
///
/// Adress polygon inside of ExPolygon
/// Keep information about source of vertex:
///     - from face (one of 2 possible)
///     - from edge (one of 2 possible)
///
/// V1~~~~V2
/// : f1 /|
/// :   / |
/// :  /e1|
/// : /   |e2
/// :/ f2 |
/// V1'~~~V2'
///
/// | .. edge
/// / .. edge
/// : .. foreign edge - neighbor
/// ~ .. no care edge - idealy should not cross model
/// V1,V1' .. projected 2d point to 3d
/// V2,V2' .. projected 2d point to 3d
///
/// f1 .. text_face_1 (triangle face made by side of shape contour)
/// f2 .. text_face_2
/// e1 .. text_edge_1 (edge on side of face made by side of shape contour)
/// e2 .. text_edge_2
///
/// </summary>
struct IntersectingElement
{
    // Base of the zero'th point of a contour in text mesh.
    // There are two vertices (front and rear) created for each contour,
    // thus there are 2x more vertices in text mesh than the number of contour points.
    // a.k.a offset of vertex inside vertices
    int32_t vertex_base{-1};

    // index of point in Polygon contour
    int32_t point_index{-1};

    // vertex or edge ID, where edge ID is the index of the source point.
    // There are 4 consecutive indices generated for a single contour edge:
    // 0th - 1st text edge (straight)
    // 1th - 1st text face
    // 2nd - 2nd text edge (diagonal)
    // 3th - 2nd text face
    // Type of intersecting element from extruded shape( 3d )
    enum class Type {
        edge_1 = 0,
        face_1 = 1,
        edge_2 = 2,
        face_2 = 3,

        undefined = 4
    } type = Type::undefined;

    // order of point in polygon for detect place between first and last point
    bool is_first{false};
    bool is_last{false};
    
    IntersectingElement &set_type(Type t)
    {
        type = t;
        return *this;
    }
};

/// <summary>
/// Convert triangle mesh model to CGAL Surface_mesh
/// Add property map for source face index
/// </summary>
/// <param name="its">Model</param>
/// <returns>CGAL mesh - half edge mesh</returns>
CutMesh to_cgal(const indexed_triangle_set &its);

/// <summary>
/// Covert 2d shape (e.g. Glyph) to CGAL model
/// </summary>
/// <param name="shapes">2d shapes to project</param>
/// <param name="projection">Define transformation 2d point into 3d</param>
/// <param name="edge_shape_map_name">Name of property map to store conversion from edge to contour</param> 
/// <param name="face_shape_map_name">Name of property map to store conversion from face to contour</param>
/// <returns>CGAL model of extruded shape</returns>
CutMesh to_cgal(const Slic3r::ExPolygons       &shapes,
                const Slic3r::Emboss::IProject &projection,
                const std::string              &edge_shape_map_name,
                const std::string              &face_shape_map_name);

enum class FaceType {
    // face inside of the cutted shape
    inside,
    // face outside of the cutted shape
    outside,
    // face without constrained edge (In or Out)
    not_constrained
};

using FaceTypeMap = CutMesh::Property_map<CutMesh::Face_index, FaceType>;
using VertexShapeMap = CutMesh::Property_map<CutMesh::Vertex_index, IntersectingElement>;

/// <summary>
/// 
/// </summary>
/// <param name="mesh">Mesh to process</param>
/// <param name="face_type_map">Output map with type of faces</param>
/// <param name="vertex_shape_map">Keep information about source element of Face type</param>
/// <param name="ecm"></param>
/// <param name="project">projection of opoint</param>
/// <param name="shape_mesh">Vertices of mesh made by shapes</param>
void set_face_type(const CutMesh                     &mesh,
                   FaceTypeMap                       &face_type_map,
                   const VertexShapeMap              &vertex_shape_map,
                   const EcmType                     &ecm,
                   const Slic3r::Emboss::IProject    &project,
                   const CutMesh                     &shape_mesh);

void flood_fill_inner(const CutMesh &mesh, FaceTypeMap &face_type_map);

/// <summary>
/// Debug purpose store of mesh with colored face by face type
/// </summary>
/// <param name="mesh">Input mesh, could add property color</param>
/// <param name="face_type_map">Keep face type</param>
/// <param name="file">File to store</param>
void store(CutMesh &mesh, const FaceTypeMap &face_type_map, const std::string& file)
{
    auto color_prop = mesh.property_map<priv::CutMesh::Face_index, CGAL::Color>("f:color");
    if (!color_prop.second)
        color_prop = mesh.add_property_map<priv::CutMesh::Face_index, CGAL::Color>("f:color");
    auto face_colors = color_prop.first;    
    for (auto fi : mesh.faces()) { 
        auto &color = face_colors[fi];
        switch (face_type_map[fi]) {
        case FaceType::inside: color = CGAL::Color{255, 0, 0}; break;
        case FaceType::outside: color = CGAL::Color{255, 0, 255}; break;
        case FaceType::not_constrained: color = CGAL::Color{0, 255, 0}; break;
        }
    }
    CGAL::IO::write_OFF(file, mesh);
}

/// <summary>
/// Track source of intersection 
/// Anotate inner and outer face, not anotated face should be not not constrained
/// </summary>
struct Visitor {
    const CutMesh &object;
    const CutMesh &shape;

    // Properties of the shape mesh:
    CutMesh::Property_map<CutMesh::Edge_index, IntersectingElement> edge_shape_map;
    CutMesh::Property_map<CutMesh::Face_index, IntersectingElement> face_shape_map;

    // Properties of the object mesh.
    CutMesh::Property_map<CutMesh::Vertex_index, IntersectingElement> vert_shape_map;

    using GT = boost::graph_traits<CutMesh>;
    using halfedge_descriptor = GT::halfedge_descriptor;
    
    // keep source of intersection for each intersection
    std::vector<const IntersectingElement*> intersections;

    /// <summary>
    /// Called when a new intersection point is detected.
    /// The intersection is detected using a face of tm_f and an edge of tm_e.
    /// Intersecting an edge hh_edge from tm_f with a face h_e of tm_e.
    /// https://doc.cgal.org/latest/Polygon_mesh_processing/classPMPCorefinementVisitor.html#a00ee0ca85db535a48726a92414acda7f
    /// </summary>
    /// <param name="i_id">The id of the intersection point, starting at 0. Ids are consecutive.</param>
    /// <param name="sdim">Dimension of a simplex part of face(h_e) that is intersected by hh_edge:
    /// 0 for vertex: target(h_e)
    /// 1 for edge: h_e
    /// 2 for the interior of face: face(h_e) </param>
    /// <param name="h_f">
    /// A halfedge from tm_f indicating the simplex intersected: 
    /// if sdim==0 the target of h_f is the intersection point, 
    /// if sdim==1 the edge of h_f contains the intersection point in its interior,
    /// if sdim==2 the face of h_f contains the intersection point in its interior.
    /// @Vojta: Edge of tm_f, see is_target_coplanar & is_source_coplanar whether any vertex of h_f is coplanar with face(h_e).
    /// </param>
    /// <param name="h_e">A halfedge from tm_e
    /// @Vojta: Vertex, halfedge or face of tm_e intersected by h_f, see comment at sdim.
    /// </param>
    /// <param name="tm_f">Mesh containing h_f</param>
    /// <param name="tm_e">Mesh containing h_e</param>
    /// <param name="is_target_coplanar">True if the target of h_e is the intersection point
    /// @Vojta: source(h_f) is coplanar with face(made by h_e).</param>
    /// <param name="is_source_coplanar">True if the source of h_e is the intersection point
    /// @Vojta: target(h_f) is coplanar with face(h_e).</param>
    void intersection_point_detected(std::size_t         i_id,
                                     int                 sdim,
                                     halfedge_descriptor h_f,
                                     halfedge_descriptor h_e,
                                     const CutMesh      &tm_f,
                                     const CutMesh      &tm_e,
                                     bool                is_target_coplanar,
                                     bool                is_source_coplanar)
    {
        if (i_id <= intersections.size()) {
            intersections.reserve(Slic3r::next_highest_power_of_2(i_id + 1));
            intersections.resize(i_id + 1);
        }

        const IntersectingElement* intersection_ptr = nullptr;
        if (&tm_e == &shape) {
            assert(&tm_f == &object);
            switch (sdim) {
            case 1:
                // edge x edge intersection
                intersection_ptr = &edge_shape_map[shape.edge(h_e)];
                break;
            case 2:
                // edge x face intersection
                intersection_ptr = &face_shape_map[shape.face(h_e)];
                break;
            default:
                assert(false);
            }
            if (is_target_coplanar)
                vert_shape_map[object.source(h_f)] = *intersection_ptr;
            if (is_source_coplanar)
                vert_shape_map[object.target(h_f)] = *intersection_ptr;
        } else {
            assert(&tm_f == &shape && &tm_e == &object);
            assert(!is_target_coplanar);
            assert(!is_source_coplanar);
            intersection_ptr = &edge_shape_map[shape.edge(h_f)];
            if (sdim == 0)
                vert_shape_map[object.target(h_e)] = *intersection_ptr;
        }
        intersections[i_id] = intersection_ptr;
    }

    using vertex_descriptor = GT::vertex_descriptor;
    void new_vertex_added(std::size_t node_id, vertex_descriptor vh, const CutMesh &tm)
    {
        assert(&tm == &object);
        assert(node_id < intersections.size());
        const IntersectingElement * intersection_ptr = intersections[node_id];
        assert(intersection_ptr != nullptr);
        assert(intersection_ptr->point_index != -1);
        vert_shape_map[vh] = *intersection_ptr; // copy ?!?
    }

    using face_descriptor = GT::face_descriptor;
    void before_subface_creations(face_descriptor /* f_old */, CutMesh &/* mesh */){}
    void after_subface_created(face_descriptor /* f_new */, CutMesh &/* mesh */) {}
    void after_subface_creations(CutMesh&) {}
    void before_subface_created(CutMesh&) {}
    void before_edge_split(halfedge_descriptor /* h */, CutMesh& /* tm */) {}
    void edge_split(halfedge_descriptor /* hnew */, CutMesh& /* tm */) {}
    void after_edge_split() {}
    void add_retriangulation_edge(halfedge_descriptor /* h */, CutMesh& /* tm */) {}
};
} // namespace privat

using namespace Slic3r;

void Slic3r::append(SurfaceCut &sc, SurfaceCut &&sc_add)
{
    if (sc.empty()) {
        sc = std::move(sc_add);
        return;
    }

    if (!sc_add.cut.empty()) {
        SurfaceCut::Index offset = static_cast<SurfaceCut::Index>(
            sc.vertices.size());
        size_t require = sc.cut.size() + sc_add.cut.size();
        if (sc.cut.capacity() < require) sc.cut.reserve(require);
        for (std::vector<SurfaceCut::Index> &cut : sc_add.cut)
            for (SurfaceCut::Index &i : cut) i += offset;
        append(sc.cut, std::move(sc_add.cut));
    }
    its_merge(sc, std::move(sc_add));
}

SurfaceCut Slic3r::cut_surface(const indexed_triangle_set &model,
                               const ExPolygons           &shapes,
                               const Emboss::IProject     &projection)
{
    priv::CutMesh cgal_model = priv::to_cgal(model);

    std::string edge_shape_map_name = "e:IntersectingElement";
    std::string face_shape_map_name = "f:IntersectingElement";
    priv::CutMesh cgal_shape = priv::to_cgal(shapes, projection, edge_shape_map_name, face_shape_map_name);    

    auto& edge_shape_map = cgal_shape.property_map<priv::CutMesh::Edge_index, priv::IntersectingElement>(edge_shape_map_name).first;    
    auto& face_shape_map = cgal_shape.property_map<priv::CutMesh::Face_index, priv::IntersectingElement>(face_shape_map_name).first;
    
    std::string vert_shape_map_name = "v:IntersectingElement";
    auto& vert_shape_map = cgal_model.add_property_map<priv::CutMesh::Vertex_index, priv::IntersectingElement>(vert_shape_map_name).first;

    priv::Visitor visitor{cgal_model, cgal_shape, edge_shape_map, face_shape_map, vert_shape_map};

    // bool map for affected edge
    priv::EcmType ecm = get(priv::DynamicEdgeProperty(), cgal_model);

    const auto& p = CGAL::Polygon_mesh_processing::parameters::throw_on_self_intersection(false).visitor(visitor).edge_is_constrained_map(ecm);
    const auto& q = CGAL::Polygon_mesh_processing::parameters::do_not_modify(true);

    CGAL::Polygon_mesh_processing::corefine(cgal_model, cgal_shape, p, q);

    std::string face_type_map_name = "f:side";    
    priv::FaceTypeMap face_type_map = cgal_model.add_property_map<priv::CutMesh::Face_index, priv::FaceType>(face_type_map_name).first;

    // Select inside and outside face in model
    priv::set_face_type(cgal_model, face_type_map, vert_shape_map, ecm, projection, cgal_shape);
    priv::store(cgal_model, face_type_map, "C:/data/temp/constrained.off"); // only debug
    
    // Seed fill the other faces inside the region.
    priv::flood_fill_inner(cgal_model, face_type_map);
    priv::store(cgal_model, face_type_map, "C:/data/temp/filled.off"); // only debug

    SurfaceCut result;
    // for (const ExPolygon& shape : shapes)
    //     append(result, cut_surface(model, shape, projection));
    return result;
}

priv::CutMesh priv::to_cgal(const indexed_triangle_set &its)
{
    CutMesh result;
    if (its.empty()) return result;

    const std::vector<stl_vertex>                  &vertices = its.vertices;
    const std::vector<stl_triangle_vertex_indices> &indices = its.indices;

    size_t vertices_count = vertices.size();
    size_t faces_count    = indices.size();
    size_t edges_count    = (faces_count * 3) / 2;
    result.reserve(vertices_count, edges_count, faces_count);

    for (const stl_vertex &v : vertices)
        result.add_vertex(CutMesh::Point{v.x(), v.y(), v.z()});

    using VI = CutMesh::Vertex_index;
    for (const stl_triangle_vertex_indices &f : indices)
        result.add_face(static_cast<VI>(f[0]), 
                        static_cast<VI>(f[1]),
                        static_cast<VI>(f[2]));
    
    return result;
}

priv::CutMesh priv::to_cgal(const Slic3r::ExPolygons        &shapes,
                            const Slic3r::Emboss::IProject &projection,
                            const std::string           &edge_shape_map_name,
                            const std::string           &face_shape_map_name)
{
    CutMesh result;
    if (shapes.empty()) return result;

    auto edge_shape_map = result.add_property_map<CutMesh::Edge_index, IntersectingElement>(edge_shape_map_name).first;
    auto face_shape_map = result.add_property_map<CutMesh::Face_index, IntersectingElement>(face_shape_map_name).first;

    std::vector<CutMesh::Vertex_index> indices;
    auto insert_contour = [&projection, &indices, &result, 
        &edge_shape_map, &face_shape_map]
        (const Polygon &polygon) {
        indices.clear();
        indices.reserve(polygon.points.size() * 2);
        size_t  num_vertices_old = result.number_of_vertices();
        for (const Point &p2 : polygon.points) {
            auto p  = projection.project(p2);
            CutMesh::Point v_first{p.first.x(), p.first.y(), p.first.z()};
            CutMesh::Point v_second{p.second.x(), p.second.y(), p.second.z()};

            CutMesh::Vertex_index vi = result.add_vertex(v_first);
            assert(size_t(vi) == (indices.size() + num_vertices_old));
            indices.emplace_back(vi);

            vi = result.add_vertex(v_second);
            assert(size_t(vi) == (indices.size() + num_vertices_old));
            indices.emplace_back(vi);
        }

        auto find_edge = [&result](CutMesh::Face_index   fi,
                                   CutMesh::Vertex_index from,
                                   CutMesh::Vertex_index to) {
            CutMesh::Halfedge_index hi = result.halfedge(fi);
            for (; result.target(hi) != to; hi = result.next(hi));
            assert(result.source(hi) == from);
            assert(result.target(hi) == to);
            return result.edge(hi);
        };

        int32_t contour_index = 0;
        for (int32_t i = 0; i < int32_t(indices.size()); i += 2) {
            bool    is_first  = i == 0;
            bool    is_last   = (i + 2) >= indices.size();
            int32_t j = is_last ? 0 : (i + 2);
            
            auto fi1 = result.add_face(indices[i], indices[i + 1], indices[j]);
            auto ei1 = find_edge(fi1, indices[i], indices[i + 1]);
            auto ei2 = find_edge(fi1, indices[i + 1], indices[j]);
            auto fi2 = result.add_face(indices[j], indices[i + 1], indices[j + 1]);
            IntersectingElement element {num_vertices_old, contour_index, IntersectingElement::Type::undefined, is_first, is_last};
            edge_shape_map[ei1] = element.set_type(IntersectingElement::Type::edge_1);
            face_shape_map[fi1] = element.set_type(IntersectingElement::Type::face_1);
            edge_shape_map[ei2] = element.set_type(IntersectingElement::Type::edge_2);
            face_shape_map[fi2] = element.set_type(IntersectingElement::Type::face_2);
            ++contour_index;
        }
    };

    size_t count_point = count_points(shapes);
    result.reserve(result.number_of_vertices() + 2 * count_point,
                   result.number_of_edges() + 4 * count_point,
                   result.number_of_faces() + 2 * count_point);

    // Identify polygon
    for (const ExPolygon &shape : shapes) {
        insert_contour(shape.contour);
        for (const Polygon &hole : shape.holes)
            insert_contour(hole);
    }
    return result;
}

void priv::set_face_type(const CutMesh                     &mesh,
                         FaceTypeMap                       &face_type_map,
                         const VertexShapeMap              &vertex_shape_map,
                         const EcmType                     &ecm,
                         const Emboss::IProject            &project,
                         const CutMesh                     &shape_mesh)
{
    size_t count = 0;
    for (auto& fi : mesh.faces()) {
        FaceType face_type = FaceType::not_constrained;
        auto     hi_end    = mesh.halfedge(fi);
        auto     hi        = hi_end;

        do {
            CGAL::SM_Edge_index edge_index = mesh.edge(hi);
            // is edge new created - constrained?
            if (get(ecm, edge_index)) {
                // This face has a constrained edge.
                IntersectingElement shape_from = vertex_shape_map[mesh.source(hi)];
                IntersectingElement shape_to = vertex_shape_map[mesh.target(hi)];
                
                assert(shape_from.point_index != -1);
                assert(shape_from.type != IntersectingElement::Type::undefined);
                assert(shape_to.point_index != -1);
                assert(shape_to.type != IntersectingElement::Type::undefined);
                
                // assert mean: There is constrained between two shapes
                // Filip think it can't happens. 
                // consider what to do?
                assert(shape_from.vertex_base == shape_to.vertex_base);

                bool is_inside = false;

                // index into contour
                int32_t i_from = shape_from.point_index;
                int32_t i_to   = shape_to.point_index;
                if (i_from == i_to && shape_from.type == shape_to.type) {
                    // intersecting element must be face
                    assert(shape_from.type == IntersectingElement::Type::face_1 ||
                           shape_from.type == IntersectingElement::Type::face_2);

                    // count of vertices is twice as count of point in the contour
                    int i = i_from * 2;
                    // j is next contour point in vertices
                    int j = shape_from.is_last ? 0 : i + 2;
                    i += shape_from.vertex_base;
                    j += shape_from.vertex_base;

                    // opposit point(in triangle face) to edge
                    const auto &p = mesh.point(mesh.target(mesh.next(hi)));

                    // abc is source triangle face
                    auto abcp =
                        shape_from.type == IntersectingElement::Type::face_1 ?
                            CGAL::orientation(
                                shape_mesh.point(CGAL::SM_Vertex_index(i)),
                                shape_mesh.point(CGAL::SM_Vertex_index(i + 1)),
                                shape_mesh.point(CGAL::SM_Vertex_index(j)), p) :
                        //shape_from.type == IntersectingElement::Type::face_2
                            CGAL::orientation(
                                shape_mesh.point(CGAL::SM_Vertex_index(j)),
                                shape_mesh.point(CGAL::SM_Vertex_index(i + 1)),
                                shape_mesh.point(CGAL::SM_Vertex_index(j + 1)), p);
                    is_inside = abcp == CGAL::POSITIVE;
                } else if (i_from < i_to || (i_from == i_to && shape_from.type < shape_to.type)) {
                    // TODO: check that it is continous indices of contour
                    bool is_last = shape_from.is_first && shape_to.is_last;
                    if (!is_last) is_inside = true;
                } else { // i_from > i_to || (i_from == i_to && shape_from.type > shape_to.type)
                    // TODO: check that it is continous indices of contour
                    bool is_last = shape_to.is_first && shape_from.is_last;
                    if (is_last) is_inside = true;
                }

                if (is_inside) {
                    // Is this face oriented towards p or away from p?
                    const auto &a = mesh.point(mesh.source(hi));
                    const auto &b = mesh.point(mesh.target(hi));
                    const auto &c = mesh.point(mesh.target(mesh.next(hi)));
                    
                    Vec3f a_(a.x(), a.y(), a.z());
                    Vec3f p_ = project.project(a_);
                    CGAL::Epick::Point_3 p{p_.x(), p_.y(), p_.z()};
                    auto abcp = CGAL::orientation(a, b, c, p);
                    if (abcp == CGAL::POSITIVE)
                        face_type = FaceType::inside;
                    else
                        is_inside = false;
                }
                if (!is_inside) face_type = FaceType::outside;
                break;
            }
            // next half edge index inside of face
            hi = mesh.next(hi);
        } while (hi != hi_end);
        face_type_map[fi] = face_type;
    }
}

void priv::flood_fill_inner(const CutMesh &mesh, FaceTypeMap &face_type_map)
{
    for (Visitor::face_descriptor fi : mesh.faces()) {
        if (face_type_map[fi] != FaceType::not_constrained) continue;

        // check if neighbor face is inside
        Visitor::halfedge_descriptor hi     = mesh.halfedge(fi);
        Visitor::halfedge_descriptor hi_end = hi;

        bool has_inside_neighbor = false;
        std::vector<CutMesh::Face_index> queue;
        do {
            Visitor::face_descriptor fi_opposite = mesh.face(mesh.opposite(hi));
            FaceType side = face_type_map[fi_opposite];
            if (side == FaceType::inside) {
                has_inside_neighbor = true;
            } else if (side == FaceType::not_constrained) {
                queue.emplace_back(fi_opposite);
            }
            hi = mesh.next(hi);
        } while (hi != hi_end);
        if (!has_inside_neighbor) continue;
        face_type_map[fi] = FaceType::inside;
        while (!queue.empty()) {
            Visitor::face_descriptor fi = queue.back();
            queue.pop_back();
            // Do not fill twice
            if (face_type_map[fi] == FaceType::inside) continue;
            face_type_map[fi] = FaceType::inside;

            // check neighbor triangle
            Visitor::halfedge_descriptor hi = mesh.halfedge(fi);
            Visitor::halfedge_descriptor hi_end = hi;
            do {
                Visitor::face_descriptor fi_opposite = mesh.face(mesh.opposite(hi));
                FaceType side = face_type_map[fi_opposite];
                if (side == FaceType::not_constrained)
                    queue.emplace_back(fi_opposite);
                hi = mesh.next(hi);
            } while (hi != hi_end);
        }
    }
}