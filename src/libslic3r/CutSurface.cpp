#include "CutSurface.hpp"

using namespace Slic3r;

void Slic3r::append(SurfaceCut &sc, SurfaceCut &&sc_add)
{
    if (sc.empty()) {
        sc = std::move(sc_add);
        return;
    }

    if (!sc_add.contours.empty()) {
        SurfaceCut::Index offset = static_cast<SurfaceCut::Index>(
            sc.vertices.size());
        size_t require = sc.contours.size() + sc_add.contours.size();
        if (sc.contours.capacity() < require) sc.contours.reserve(require);
        for (std::vector<SurfaceCut::Index> &cut : sc_add.contours)
            for (SurfaceCut::Index &i : cut) i += offset;
        append(sc.contours, std::move(sc_add.contours));
    }
    its_merge(sc, std::move(sc_add));
}

SurfaceCut Slic3r::merge(SurfaceCuts &&cuts) { 
    SurfaceCut result; 
    for (SurfaceCut &cut : cuts) 
        append(result, std::move(cut));
    return result;
}

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian_converter.h>

// libslic3r
#include "TriangleMesh.hpp" // its_merge
#include "Utils.hpp" // next_highest_power_of_2

namespace priv {
    
using EpicKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CutMesh = CGAL::Surface_mesh<EpicKernel::Point_3>;
//    using EpecKernel = CGAL::Exact_predicates_exact_constructions_kernel;
//    using CutMesh = CGAL::Surface_mesh<EpecKernel::Point_3>;

using DynamicEdgeProperty = CGAL::dynamic_edge_property_t<bool>;
using SMPM = boost::property_map<priv::CutMesh, DynamicEdgeProperty>::SMPM; 
using EcmType = CGAL::internal::Dynamic<priv::CutMesh, SMPM>;

using VI = CGAL::SM_Vertex_index;
using HI = CGAL::SM_Halfedge_index;
using EI = CGAL::SM_Edge_index;
using FI = CGAL::SM_Face_index;

/// <summary>
/// IntersectingElement
///
/// Adress polygon inside of ExPolygon
/// Keep information about source of vertex:
///     - from face (one of 2 possible)
///     - from edge (one of 2 possible)
///
///  V1~~~~~V2
///   | f1 /:
///   |   / :
/// e1|  /e2:
///   | /   :
///   |/ f2 :
///  V1'~~~~V2'
///
/// | .. edge
/// / .. edge
/// : .. foreign edge - neighbor
/// ~ .. no care edge - idealy should not cross model
/// V1,V1' .. projected 2d point to 3d
/// V2,V2' .. projected 2d point to 3d
/// 
/// Vertex indexing
/// V1  .. i (vertex_base + 2x index of point in polygon)
/// V1' .. i + 1
/// V2  .. j = i + 2 || 0 (for last i in polygon)
/// V2' .. j + 1
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
    uint32_t vertex_base{std::numeric_limits<uint32_t>::max()};

    // index of point in Polygon contour
    uint32_t point_index{std::numeric_limits<uint32_t>::max()};

    // store together type, is_first, is_last
    unsigned char attr;

    // vertex or edge ID, where edge ID is the index of the source point.
    // There are 4 consecutive indices generated for a single contour edge:
    // 0th - 1st text edge (straight)
    // 1th - 1st text face
    // 2nd - 2nd text edge (diagonal)
    // 3th - 2nd text face
    // Type of intersecting element from extruded shape( 3d )
    // NOTE: type must be storable to 3bit -> max value is 7
    enum class Type: unsigned char {
        edge_1 = 0,
        face_1 = 1,
        edge_2 = 2,
        face_2 = 3,
        undefined = 4
    };
        
    IntersectingElement &set_type(Type t)
    {
        attr = static_cast<unsigned char>(
            attr + (int) t - (int) get_type());
        return *this;
    }
    void set_is_first(){ attr += 8; }
    void set_is_last(){ attr += 16; }
    Type get_type() const { return static_cast<Type>(attr % 8);}
    bool is_first() const { return 8 <= attr && attr < 16; }
    bool is_last() const { return attr >= 16; }
};

/// <summary>
/// Convert triangle mesh model to CGAL Surface_mesh
/// Add property map for source face index
/// </summary>
/// <param name="its">Model</param>
/// <returns>CGAL mesh - half edge mesh</returns>
CutMesh to_cgal(const indexed_triangle_set &its);

using Project = Emboss::IProject;
/// <summary>
/// Covert 2d shape (e.g. Glyph) to CGAL model
/// </summary>
/// <param name="shapes">2d shapes to project</param>
/// <param name="projection">Define transformation 2d point into 3d</param>
/// <param name="edge_shape_map_name">Name of property map to store conversion from edge to contour</param> 
/// <param name="face_shape_map_name">Name of property map to store conversion from face to contour</param>
/// <returns>CGAL model of extruded shape</returns>
CutMesh to_cgal(const ExPolygons  &shapes,
                const Project     &projection,
                const std::string &edge_shape_map_name,
                const std::string &face_shape_map_name);

using VertexShapeMap = CutMesh::Property_map<VI, const IntersectingElement *>;
/// <summary>
/// Track source of intersection 
/// Help for anotate inner and outer faces
/// </summary>
struct Visitor {
    const CutMesh &object;
    const CutMesh &shape;

    // Properties of the shape mesh:
    CutMesh::Property_map<EI, IntersectingElement> edge_shape_map;
    CutMesh::Property_map<FI, IntersectingElement> face_shape_map;

    // Properties of the object mesh.
    VertexShapeMap vert_shape_map;

    // check for anomalities
    bool* is_valid;

    // keep source of intersection for each intersection
    // used to copy data into vert_shape_map
    std::vector<const IntersectingElement*> intersections;

    /// <summary>
    /// Called when a new intersection point is detected.
    /// The intersection is detected using a face of tm_f and an edge of tm_e.
    /// Intersecting an edge hh_edge from tm_f with a face h_e of tm_e.
    /// https://doc.cgal.org/latest/Polygon_mesh_processing/classPMPCorefinementVisitor.html#a00ee0ca85db535a48726a92414acda7f
    /// </summary>
    /// <param name="i_id">The id of the intersection point, starting at 0. Ids are consecutive.</param>
    /// <param name="sdim">Dimension of a simplex part of face(h_e) that is intersected by edge(h_f):
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
    void intersection_point_detected(std::size_t    i_id,
                                     int            sdim,
                                     HI             h_f,
                                     HI             h_e,
                                     const CutMesh &tm_f,
                                     const CutMesh &tm_e,
                                     bool           is_target_coplanar,
                                     bool           is_source_coplanar);

    /// <summary>
    /// Called when a new vertex is added in tm (either an edge split or a vertex inserted in the interior of a face).
    /// Fill vertex_shape_map by intersections
    /// </summary>
    /// <param name="i_id">Order number of intersection point</param>
    /// <param name="v">New added vertex</param>
    /// <param name="tm">Affected mesh</param>
    void new_vertex_added(std::size_t i_id, VI v, const CutMesh &tm);

    // Not used visitor functions
    void before_subface_creations(FI /* f_old */, CutMesh &/* mesh */){}
    void after_subface_created(FI /* f_new */, CutMesh &/* mesh */) {}
    void after_subface_creations(CutMesh&) {}
    void before_subface_created(CutMesh&) {}
    void before_edge_split(HI /* h */, CutMesh& /* tm */) {}
    void edge_split(HI /* hnew */, CutMesh& /* tm */) {}
    void after_edge_split() {}
    void add_retriangulation_edge(HI /* h */, CutMesh& /* tm */) {}
};

/// <summary>
/// Flag for faces in CGAL mesh
/// </summary>
enum class FaceType {
    // face inside of the cutted shape
    inside,
    // face outside of the cutted shape
    outside,
    // face without constrained edge (In or Out)
    not_constrained,

    // Helper flag that inside was processed
    inside_
};
using FaceTypeMap = CutMesh::Property_map<FI, FaceType>;
/// <summary>
/// Face with constrained edge are inside/outside by type of intersection
/// Other set to not_constrained(still it could be inside/outside)
/// </summary>
/// <param name="face_type_map">[Output] property map with type of faces</param>
/// <param name="mesh">Mesh to process</param>
/// <param name="vertex_shape_map">Keep information about source element of Face
/// type</param> <param name="ecm">Dynamic Edge Constrained Map of bool</param>
/// <param name="project">projection of opoint</param>
/// <param name="shape_mesh">Vertices of mesh made by shapes</param>
void set_face_type(FaceTypeMap          &face_type_map,
                   const CutMesh        &mesh,
                   const VertexShapeMap &vertex_shape_map,
                   const EcmType        &ecm,
                   const Project        &project,
                   const CutMesh        &shape_mesh);


/// <summary>
/// Check orientation(normal direction) of face on mesh
/// </summary>
/// <param name="fi">Face index to inspect</param>
/// <param name="mesh">Mesh contained fi</param>
/// <param name="projection">Define direction of projection</param>
/// <returns>TRUE for cutted face otherwise FALSE</returns>
bool is_toward_projection(FI             fi,
                          const CutMesh &mesh,
                          const Project &projection);

/// <summary>
/// Change FaceType from not_constrained to inside
/// For neighbor(or neighbor of neighbor of ...) of inside triangles.
/// Process only not_constrained triangles
/// </summary>
/// <param name="mesh">Corefined mesh</param>
/// <param name="projection">Projection from 2d to 3d</param>
/// <param name="face_type_map">In/Out map with faces type</param>
void flood_fill_inner(const CutMesh &mesh,
                      const Project &projection,
                      FaceTypeMap   &face_type_map);

using ReductionMap = CutMesh::Property_map<VI, VI>;
/// <summary>
/// Create map to reduce unnecesary triangles,
/// Triangles are made by divided quad to two triangles
/// on side of cutting shape mesh
/// </summary>
/// <param name="reduction_map">Reduction map from vertex to vertex, 
/// when key == value than no reduction</param>
/// <param name="faces">Faces of one </param> 
/// <param name="mesh">Input object</param> 
/// <param name="face_type_map">Type of shape inside / outside</param> 
/// <param name="vert_shape_map">Source of outline vertex</param>
void create_reduce_map(ReductionMap         &reduction_map,
                       const CutMesh        &mesh,
                       const FaceTypeMap    &face_type_map,
                       const VertexShapeMap &vert_shape_map);

// connected faces(triangles) and outlines(halfEdges) for one surface cut
using CutAOI = std::pair<std::vector<FI>, std::vector<HI>>;
using CutAOIs = std::vector<CutAOI>;

/// <summary>
/// Create areas from mesh surface
/// </summary>
/// <param name="mesh">Model</param>
/// <param name="shapes">Cutted shapes</param>
/// <param name="face_type_map">Define Triangles of interest.
/// Edge between inside / outside.
/// NOTE: Not const because it need to flag proccessed faces</param>
/// <returns>Areas of interest from mesh</returns>
CutAOIs create_cut_area_of_interests(const CutMesh    &mesh,
                                     const ExPolygons &shapes,
                                     FaceTypeMap      &face_type_map);

/// <summary>
/// Filter out cuts which are behind another.
/// Prevent overlapping embossed shape in space.
/// </summary>
/// <param name="cuts">AOIs</param>
/// <param name="mesh">triangle model</param>
/// <param name="shapes">2d cutted shapes</param>
/// <param name="projection">Projection from 2d to 3d</param>
/// <param name="vert_shape_map">Identify source of intersection</param>
void filter_cuts(CutAOIs              &cuts,
                 const CutMesh        &mesh,
                 const ExPolygons     &shapes,
                 const Project        &projection,
                 const VertexShapeMap &vert_shape_map);

using ConvertMap = CutMesh::Property_map<VI, SurfaceCut::Index>;
/// <summary>
/// Create surface cuts from mesh model
/// </summary>
/// <param name="mesh">Model</param>
/// <param name="shapes">Cutted shapes</param>
/// <param name="reduction_map">Reduction of vertices</param>
/// <param name="face_type_map">Define Triangles of interest.
/// Edge between inside / outside.
/// NOTE: Not const because it need to flag proccessed faces</param>
/// <param name="convert_map">Used only inside function. 
/// Store conversion from mesh to result.</param>
/// <returns>Created surface cuts</returns>
SurfaceCuts create_surface_cuts(const CutAOIs      &cutAOIs,
                                const CutMesh      &mesh,
                                const ReductionMap &reduction_map,
                                ConvertMap         &convert_map);


/// <summary>
/// Collect connected inside faces
/// Collect outline half edges
/// </summary>
/// <param name="process">Queue of face to process - find connected</param>
/// <param name="faces">[Output] collected Face indices from mesh</param>
/// <param name="outlines">[Output] collected Halfedge indices from mesh</param>
/// <param name="face_type_map">Use flag inside / outside
/// NOTE: Modify in function: inside -> inside_</param>
/// <param name="mesh">mesh to process</param>
void collect_surface_data(std::queue<FI>  &process,
                          std::vector<FI> &faces,
                          std::vector<HI> &outlines,
                          FaceTypeMap     &face_type_map,
                          const CutMesh   &mesh);

/// <summary>
/// Copy triangles from CGAL mesh into index triangle set
/// NOTE: Skip vertices created by edge in center of Quad.
/// </summary>
/// <param name="faces">Faces to copy</param>
/// <param name="count_outlines">Count of outlines</param>
/// <param name="mesh">Source CGAL mesh</param>
/// <param name="reduction_map">Reduction of vertices</param>
/// <param name="v2v">[Output] map to convert CGAL vertex to its::vertex index</param>
/// <returns>Surface cut (Partialy filled - only index triangle set)</returns>
SurfaceCut create_index_triangle_set(const std::vector<FI> &faces,
                                     size_t                 count_outlines,
                                     const CutMesh         &mesh,
                                     const ReductionMap    &reduction_map,
                                     ConvertMap            &v2v);

/// <summary>
/// Connect outlines into closed loops
/// </summary>
/// <param name="outlines">Half edges from border of cut - Oriented</param>
/// <param name="mesh">Source CGAL mesh</param>
/// <param name="reduction_map">Reduction of vertices</param>
/// <param name="v2v">Map to convert CGAL vertex to its::vertex</param>
/// <returns>Cuts - outlines of surface</returns>
SurfaceCut::CutContour create_cut(const std::vector<HI> &outlines,
                               const CutMesh         &mesh,
                               const ReductionMap &reduction_map,
                               const ConvertMap      &v2v);

/// <summary>
/// Debug purpose store of mesh with colored face by face type
/// </summary>
/// <param name="mesh">Input mesh, could add property color
/// NOTE: Not const because need to [optionaly] append color property map</param>
/// <param name="face_type_map">Color source</param>
/// <param name="file">File to store</param>
void store(CutMesh &mesh, const FaceTypeMap &face_type_map, const std::string &file);
void store(CutMesh &mesh, const ReductionMap &reduction_map, const std::string &file);
void store(const SurfaceCuts &cut, const std::string &file_prefix);
} // namespace privat

SurfaceCuts Slic3r::cut_surface(const indexed_triangle_set &model,
                                const ExPolygons           &shapes,
                                const Emboss::IProject     &projection)
{
    priv::CutMesh cgal_model = priv::to_cgal(model);
    CGAL::IO::write_OFF("C:/data/temp/model.off", cgal_model); // only debug

    std::string edge_shape_map_name = "e:IntersectingElement";
    std::string face_shape_map_name = "f:IntersectingElement";
    priv::CutMesh cgal_shape = priv::to_cgal(shapes, projection, edge_shape_map_name, face_shape_map_name);   
    CGAL::IO::write_OFF("C:/data/temp/shape.off", cgal_shape); // only debug

    auto edge_shape_map = cgal_shape.property_map<priv::EI, priv::IntersectingElement>(edge_shape_map_name).first;    
    auto face_shape_map = cgal_shape.property_map<priv::FI, priv::IntersectingElement>(face_shape_map_name).first;
    
    std::string vert_shape_map_name = "v:IntersectingElement";
    // pointer to edge or face shape_map
    priv::VertexShapeMap vert_shape_map = cgal_model.add_property_map<priv::VI, const priv::IntersectingElement*>(vert_shape_map_name).first;
    
    // detect anomalities in visitor.
    bool is_valid = true;
    // create anotation visitor - Must be copyable
    priv::Visitor visitor{cgal_model, cgal_shape, edge_shape_map, face_shape_map, vert_shape_map, &is_valid};

    // bool map for affected edge
    priv::EcmType ecm = get(priv::DynamicEdgeProperty(), cgal_model);

    const auto &p = CGAL::parameters::visitor(visitor)
                        .edge_is_constrained_map(ecm)
                        .throw_on_self_intersection(false);
    const auto& q = CGAL::parameters::do_not_modify(true);
    CGAL::Polygon_mesh_processing::corefine(cgal_model, cgal_shape, p, q);

    if (!is_valid) return {};

    std::string face_type_map_name = "f:side";
    priv::FaceTypeMap face_type_map = cgal_model.add_property_map<priv::FI, priv::FaceType>(face_type_map_name).first;

    // Select inside and outside face in model
    priv::set_face_type(face_type_map, cgal_model, vert_shape_map, ecm, projection, cgal_shape);
    priv::store(cgal_model, face_type_map, "C:/data/temp/constrained.off"); // only debug
    
    // Seed fill the other faces inside the region.
    priv::flood_fill_inner(cgal_model, projection, face_type_map);
    priv::store(cgal_model, face_type_map, "C:/data/temp/filled.off"); // only debug

    std::string vertex_reduction_map_name = "v:reduction";
    priv::ReductionMap vertex_reduction_map = cgal_model.add_property_map<priv::VI, priv::VI>(vertex_reduction_map_name).first;
    priv::create_reduce_map(vertex_reduction_map, cgal_model, face_type_map, vert_shape_map);
    priv::store(cgal_model, vertex_reduction_map, "C:/data/temp/reduction.off"); // only debug

    priv::CutAOIs cutAOIs = create_cut_area_of_interests(cgal_model, shapes, face_type_map);

    // Filter out NO top one cuts
    priv::filter_cuts(cutAOIs, cgal_model, shapes, projection, vert_shape_map);

    // conversion map between vertex index in cgal_model and indices in result
    // used instead of std::map
    std::string vertec_convert_map_name = "v:convert";
    priv::ConvertMap vertex_convert_map = cgal_model.add_property_map<priv::VI, SurfaceCut::Index>(vertec_convert_map_name).first;    
    SurfaceCuts result = priv::create_surface_cuts(cutAOIs, cgal_model, vertex_reduction_map, vertex_convert_map);
        
    priv::store(result, "C:/data/temp/cut"); // only debug
    
    // TODO: Filter surfaceCuts to only the top most.
    return result;
}

indexed_triangle_set Slic3r::cuts2model(const SurfaceCuts      &cuts,
                                        const Emboss::IProject &projection)
{
    indexed_triangle_set result;
    size_t count_vertices = 0;
    size_t count_indices = 0;
    for (const SurfaceCut &cut : cuts) { 
        assert(!cut.empty());
        count_indices += cut.indices.size()*2;
        // indices from from zig zag
        for (const auto &c : cut.contours) {
            assert(!c.empty());
            count_indices += c.size() * 2;
        }
        count_vertices += cut.vertices.size()*2;
    }
    result.vertices.reserve(count_vertices);
    result.indices.reserve(count_indices);

    size_t indices_offset = 0;
    for (const SurfaceCut &cut : cuts) {
        // front
        for (const auto &v : cut.vertices)
            result.vertices.push_back(v);
        for (const auto &i : cut.indices) 
            result.indices.emplace_back(i.x() + indices_offset,
                                        i.y() + indices_offset,
                                        i.z() + indices_offset);

        // back
        for (const auto &v : cut.vertices) {
            Vec3f v2 = projection.project(v);
            result.vertices.push_back(v2);
        }
        size_t back_offset = indices_offset + cut.vertices.size();
        for (const auto &i : cut.indices) {
            assert(i.x() + back_offset < result.vertices.size());
            assert(i.y() + back_offset < result.vertices.size());
            assert(i.z() + back_offset < result.vertices.size());
            // Y and Z is swapped CCW triangles for back side
            result.indices.emplace_back(i.x() + back_offset,
                                        i.z() + back_offset,
                                        i.y() + back_offset);
        }

        // zig zag indices
        for (const auto &contour : cut.contours) {
            size_t prev_ci = contour.back();
            size_t prev_front_index = indices_offset + prev_ci;
            size_t prev_back_index  = back_offset + prev_ci;
            for (size_t ci : contour) {
                size_t front_index = indices_offset + ci;                
                size_t back_index  = back_offset + ci;
                assert(front_index < result.vertices.size());
                assert(prev_front_index < result.vertices.size());
                assert(back_index < result.vertices.size());
                assert(prev_back_index < result.vertices.size());
            
                result.indices.emplace_back(
                    front_index,
                    prev_front_index,
                    back_index
                );
                result.indices.emplace_back(
                    prev_front_index,
                    prev_back_index,
                    back_index
                );
                prev_front_index = front_index;
                prev_back_index  = back_index;
            }
        }

        indices_offset = result.vertices.size();
    }

    assert(count_vertices == result.vertices.size());
    assert(count_indices == result.indices.size());
    return result;
}

indexed_triangle_set Slic3r::cut2model(const SurfaceCut       &cut,
                                       const Emboss::IProject &projection)
{
    assert(!cut.empty());
    size_t count_vertices = cut.vertices.size() * 2;
    size_t count_indices  = cut.indices.size() * 2;

    // indices from from zig zag
    for (const auto &c : cut.contours) {
        assert(!c.empty());
        count_indices += c.size() * 2;
    }
    
    indexed_triangle_set result;
    result.vertices.reserve(count_vertices);
    result.indices.reserve(count_indices);

    // front
    result.vertices.insert(result.vertices.end(), 
        cut.vertices.begin(), cut.vertices.end());
    result.indices.insert(result.indices.end(), 
        cut.indices.begin(), cut.indices.end());

    // back
    for (const auto &v : cut.vertices) {
        Vec3f v2 = projection.project(v);
        result.vertices.push_back(v2);
    }

    size_t back_offset = cut.vertices.size();
    for (const auto &i : cut.indices) {
        // check range of indices in cut
        assert(i.x() >= 0 && i.x() < cut.vertices.size());
        assert(i.y() >= 0 && i.y() < cut.vertices.size());
        assert(i.z() >= 0 && i.z() < cut.vertices.size());
        // Y and Z is swapped CCW triangles for back side
        result.indices.emplace_back(i.x() + back_offset,
                                    i.z() + back_offset,
                                    i.y() + back_offset);
    }

    // zig zag indices
    for (const auto &contour : cut.contours) {
        size_t prev_front_index = contour.back();
        size_t prev_back_index  = back_offset + prev_front_index;
        for (size_t front_index : contour) {
            assert(front_index < cut.vertices.size());
            size_t back_index  = back_offset + front_index;
            result.indices.emplace_back(front_index, prev_front_index, back_index);
            result.indices.emplace_back(prev_front_index, prev_back_index, back_index);
            prev_front_index = front_index;
            prev_back_index  = back_index;
        }
    }

    assert(count_vertices == result.vertices.size());
    assert(count_indices == result.indices.size());
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

    for (const stl_triangle_vertex_indices &f : indices)
        result.add_face(static_cast<VI>(f[0]), 
                        static_cast<VI>(f[1]),
                        static_cast<VI>(f[2]));
    
    return result;
}

priv::CutMesh priv::to_cgal(const ExPolygons  &shapes,
                            const Project     &projection,
                            const std::string &edge_shape_map_name,
                            const std::string &face_shape_map_name)
{
    CutMesh result;
    if (shapes.empty()) return result;

    auto edge_shape_map = result.add_property_map<EI, IntersectingElement>(edge_shape_map_name).first;
    auto face_shape_map = result.add_property_map<FI, IntersectingElement>(face_shape_map_name).first;

    std::vector<VI> indices;
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

            VI reduction_from = result.add_vertex(v_first);
            assert(size_t(reduction_from) == (indices.size() + num_vertices_old));
            indices.emplace_back(reduction_from);

            reduction_from = result.add_vertex(v_second);
            assert(size_t(reduction_from) == (indices.size() + num_vertices_old));
            indices.emplace_back(reduction_from);
        }

        auto find_edge = [&result](FI fi, VI from, VI to) {
            HI hi = result.halfedge(fi);
            for (; result.target(hi) != to; hi = result.next(hi));
            assert(result.source(hi) == from);
            assert(result.target(hi) == to);
            return result.edge(hi);
        };

        uint32_t contour_index = 0;
        for (int32_t i = 0; i < int32_t(indices.size()); i += 2) {
            bool    is_first  = i == 0;
            bool    is_last   = size_t(i + 2) >= indices.size();
            int32_t j = is_last ? 0 : (i + 2);
            
            auto fi1 = result.add_face(indices[i], indices[i + 1], indices[j]);
            auto ei1 = find_edge(fi1, indices[i], indices[i + 1]);
            auto ei2 = find_edge(fi1, indices[i + 1], indices[j]);
            auto fi2 = result.add_face(indices[j], indices[i + 1], indices[j + 1]);
            uint32_t vertex_base = static_cast<uint32_t>(num_vertices_old);
            IntersectingElement element {vertex_base, contour_index, (unsigned char)IntersectingElement::Type::undefined};
            if (is_first) element.set_is_first();
            if (is_last) element.set_is_last();
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

void priv::set_face_type(FaceTypeMap          &face_type_map,
                         const CutMesh        &mesh,
                         const VertexShapeMap &vertex_shape_map,
                         const EcmType        &ecm,
                         const Project        &project,
                         const CutMesh        &shape_mesh)
{
    for (const FI& fi : mesh.faces()) {
        FaceType face_type = FaceType::not_constrained;

        HI hi_end = mesh.halfedge(fi);
        HI hi     = hi_end;
        do {
            EI edge_index = mesh.edge(hi);
            // is edge new created - constrained?
            if (get(ecm, edge_index)) {
                // This face has a constrained edge.
                const IntersectingElement& shape_from = *vertex_shape_map[mesh.source(hi)];
                const IntersectingElement& shape_to = *vertex_shape_map[mesh.target(hi)];
                
                assert(shape_from.point_index != std::numeric_limits<uint32_t>::max());
                assert(shape_from.attr != (unsigned char)IntersectingElement::Type::undefined);
                assert(shape_to.point_index != std::numeric_limits<uint32_t>::max());
                assert(shape_to.attr != (unsigned char)IntersectingElement::Type::undefined);
                
                // assert mean: There is constrained between two shapes
                // Filip think it can't happens. 
                // consider what to do?
                assert(shape_from.vertex_base == shape_to.vertex_base);

                bool is_inside = false;

                // index into contour
                uint32_t i_from = shape_from.point_index;
                uint32_t i_to   = shape_to.point_index;
                IntersectingElement::Type type_from = shape_from.get_type();
                IntersectingElement::Type type_to = shape_to.get_type();
                if (i_from == i_to && type_from == type_to) {
                    // intersecting element must be face
                    assert(type_from == IntersectingElement::Type::face_1 ||
                           type_from == IntersectingElement::Type::face_2);

                    // count of vertices is twice as count of point in the contour
                    uint32_t i = i_from * 2;
                    // j is next contour point in vertices
                    uint32_t j = shape_from.is_last() ? 0 : i + 2;
                    i += shape_from.vertex_base;
                    j += shape_from.vertex_base;

                    // opposit point(in triangle face) to edge
                    const auto &p = mesh.point(mesh.target(mesh.next(hi)));

                    // abc is source triangle face
                    auto abcp =
                        type_from == IntersectingElement::Type::face_1 ?
                            CGAL::orientation(
                                shape_mesh.point(VI(i)),
                                shape_mesh.point(VI(i + 1)),
                                shape_mesh.point(VI(j)), p) :
                        // type_from == IntersectingElement::Type::face_2
                            CGAL::orientation(
                                shape_mesh.point(VI(j)),
                                shape_mesh.point(VI(i + 1)),
                                shape_mesh.point(VI(j + 1)), p);
                    is_inside = abcp == CGAL::NEGATIVE;
                } else if (i_from < i_to || (i_from == i_to && type_from < type_to)) {
                    // TODO: check that it is continous indices of contour
                    bool is_last = shape_from.is_first() && shape_to.is_last() &&
                                   shape_to.vertex_base == shape_from.vertex_base;
                    if (!is_last) is_inside = true;
                } else { // i_from > i_to || (i_from == i_to && shape_from.type > shape_to.type)
                    // TODO: check that it is continous indices of contour
                    bool is_last = shape_to.is_first() && shape_from.is_last() &&
                                   shape_to.vertex_base == shape_from.vertex_base;
                    if (is_last) is_inside = true;
                }

                if (is_inside) {
                    if (is_toward_projection(fi, mesh, project))
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

bool priv::is_toward_projection(FI             fi,
                                const CutMesh &mesh,
                                const Project &projection)
{
    HI          hi = mesh.halfedge(fi);
    const auto &a  = mesh.point(mesh.source(hi));
    const auto &b  = mesh.point(mesh.target(hi));
    const auto &c  = mesh.point(mesh.target(mesh.next(hi)));

    Vec3f a_(a.x(), a.y(), a.z());
    Vec3f p_ = projection.project(a_);

    CGAL::Epick::Point_3 p{p_.x(), p_.y(), p_.z()};

    return CGAL::orientation(a, b, c, p) == CGAL::NEGATIVE;
}
 

void priv::flood_fill_inner(const CutMesh &mesh,
                            const Project &projection,
                            FaceTypeMap   &face_type_map)
{
    for (FI fi : mesh.faces()) {
        if (face_type_map[fi] != FaceType::not_constrained) continue;

        // check if neighbor(one of three in triangle) has type inside
        bool has_inside_neighbor = false;
        HI hi     = mesh.halfedge(fi);
        HI hi_end = hi;
        // list of other not constrained neighbors
        std::queue<FI> queue;
        do {
            HI hi_opposite = mesh.opposite(hi);
            // open edge doesn't have opposit half edge
            if (!hi_opposite.is_valid()) {
                hi = mesh.next(hi);
                continue;
            }
            FI fi_opposite = mesh.face(hi_opposite);
            FaceType side = face_type_map[fi_opposite];
            if (side == FaceType::inside) {
                has_inside_neighbor = true;
            } else if (side == FaceType::not_constrained) {
                queue.emplace(fi_opposite);
            }
            hi = mesh.next(hi);
        } while (hi != hi_end);


        if (!has_inside_neighbor) continue;

        face_type_map[fi] = FaceType::inside;
        while (!queue.empty()) {
            FI fi = queue.front();
            queue.pop(); 
            // Do not fill twice
            if (face_type_map[fi] == FaceType::inside) continue;
            face_type_map[fi] = FaceType::inside;

            // check neighbor triangle
            HI hi = mesh.halfedge(fi);
            HI hi_end = hi;
            do {
                HI hi_opposite = mesh.opposite(hi);
                // open edge doesn't have opposit half edge
                if (!hi_opposite.is_valid()) {
                    hi = mesh.next(hi);                    
                    continue;
                }
                FI fi_opposite = mesh.face(hi_opposite);
                FaceType &side = face_type_map[fi_opposite];
                if (side == FaceType::not_constrained) {
                    if (is_toward_projection(fi_opposite, mesh, projection)) {
                        queue.emplace(fi_opposite);
                    } else {
                        // Is in opposit direction
                        side = FaceType::outside;
                    }
                }
                hi = mesh.next(hi);
            } while (hi != hi_end);
        }
    }
}

void priv::Visitor::intersection_point_detected(std::size_t    i_id,
                                                int            sdim,
                                                HI             h_f,
                                                HI             h_e,
                                                const CutMesh &tm_f,
                                                const CutMesh &tm_e,
                                                bool is_target_coplanar,
                                                bool is_source_coplanar)
{
    if (i_id >= intersections.size()) {
        size_t capacity = Slic3r::next_highest_power_of_2(i_id + 1);
        intersections.reserve(capacity);
        intersections.resize(capacity);
    }

    const IntersectingElement *intersection_ptr = nullptr;
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
        default: assert(false);
        }
        if (is_target_coplanar)
            vert_shape_map[object.source(h_f)] = intersection_ptr;
        if (is_source_coplanar)
            vert_shape_map[object.target(h_f)] = intersection_ptr;
    } else {
        assert(&tm_f == &shape && &tm_e == &object);
        assert(!is_target_coplanar);
        assert(!is_source_coplanar);
        intersection_ptr = &edge_shape_map[shape.edge(h_f)];
        if (sdim == 0) vert_shape_map[object.target(h_e)] = intersection_ptr;
    }

    if (intersection_ptr->point_index == std::numeric_limits<uint32_t>::max()) {
        // there is unexpected intersection
        // Top (or Bottom) shape contour edge (or vertex) intersection
        // Suggest to change projection min/max limits
        *is_valid = false;
    }
    intersections[i_id] = intersection_ptr;
}

void priv::Visitor::new_vertex_added(std::size_t i_id, VI v, const CutMesh &tm)
{
    assert(&tm == &object);
    assert(i_id < intersections.size());
    const IntersectingElement *intersection_ptr = intersections[i_id];
    assert(intersection_ptr != nullptr);
    // intersection was not filled in function intersection_point_detected
    //assert(intersection_ptr->point_index != std::numeric_limits<uint32_t>::max());
    vert_shape_map[v] = intersection_ptr;
}

void priv::collect_surface_data(std::queue<FI>  &process,
                                std::vector<FI> &faces,
                                std::vector<HI> &outlines,
                                FaceTypeMap     &face_type_map,
                                const CutMesh   &mesh)
{
    while (!process.empty()) {
        FI fi_ = process.front();
        process.pop();

        // Do not process twice
        if (face_type_map[fi_] == FaceType::inside_) continue;
        assert(face_type_map[fi_] == FaceType::inside);
        // flag face as processed
        face_type_map[fi_] = FaceType::inside_;
        faces.push_back(fi_);

        // check neighbor triangle
        HI hi     = mesh.halfedge(fi_);
        HI hi_end = hi;
        do {
            HI hi_opposite = mesh.opposite(hi);
            // open edge doesn't have opposit half edge
            if (!hi_opposite.is_valid()) {
                hi = mesh.next(hi);                    
                continue;
            }
            FI fi_opposite = mesh.face(hi_opposite);
            FaceType side        = face_type_map[fi_opposite];
            if (side == FaceType::inside) {
                process.emplace(fi_opposite);
            } else if (side == FaceType::outside) {
                // store outlines
                outlines.push_back(hi);
            }
            hi = mesh.next(hi);
        } while (hi != hi_end);
    }
}

void priv::create_reduce_map(ReductionMap         &reduction_map,
                             const CutMesh        &mesh,
                             const FaceTypeMap    &face_type_map,
                             const VertexShapeMap &vert_shape_map)
{
    // IMPROVE: find better way to initialize or try use std::map
    // initialize reduction map
    for (VI reduction_from : mesh.vertices()) 
        reduction_map[reduction_from] = reduction_from;
    
    // check if vertex was made by edge_2 which is diagonal of quad
    auto is_reducible_vertex = [&vert_shape_map, &mesh](VI reduction_from) -> bool {
        const IntersectingElement *ie = vert_shape_map[reduction_from];
        if (ie == nullptr) return false;
        IntersectingElement::Type type = ie->get_type();
        return type == IntersectingElement::Type::edge_2;
    };

    /// <summary>
    /// Append reduction or change existing one.
    /// </summary>
    /// <param name="hi">HalEdge between outside and inside face.
    /// Target vertex will be reduced
    /// Source vertex left</param>
    auto add_reduction = [&reduction_map, &mesh, &is_reducible_vertex, &face_type_map]
    (HI hi) {
        VI erase = mesh.target(hi);
        VI left = mesh.source(hi);
        assert(is_reducible_vertex(erase));
        assert(!is_reducible_vertex(left));
        assert((
            FaceType::outside == face_type_map[mesh.face(hi)] && 
            FaceType::inside  == face_type_map[mesh.face(mesh.opposite(hi))] 
            ) || (
            FaceType::outside == face_type_map[mesh.face(mesh.opposite(hi))] && 
            FaceType::inside  == face_type_map[mesh.face(hi)]
        ));
        bool is_first = reduction_map[erase] == erase;
        if (is_first)
            reduction_map[erase] = left;
        // I have no better rule than take the first
        // for decide which reduction will be better
        // But it could be use only one of them
    };

    for (FI fi : mesh.faces()) {
        if (face_type_map[fi] != FaceType::inside) continue;

        // find all reducible edges
        HI hi     = mesh.halfedge(fi);
        HI hi_end = hi;
        do {
            VI reduction_from = mesh.target(hi);
            if (is_reducible_vertex(reduction_from)) {
                // reducible vertex
                VI vi_from = mesh.target(hi);

                // halfedges connected with reduction_from
                HI hi1 = hi;
                HI hi2 = mesh.next(hi);
                // faces connected with reduction_from
                FI fi1 = mesh.face(mesh.opposite(hi1));
                FI fi2 = mesh.face(mesh.opposite(hi2));

                if (face_type_map[fi1] == FaceType::outside)
                    add_reduction(hi1);
                if (face_type_map[fi2] == FaceType::outside)
                    add_reduction(mesh.opposite(hi2));
            }
            hi = mesh.next(hi);
        } while (hi != hi_end);        
    }
}

SurfaceCut priv::create_index_triangle_set(const std::vector<FI> &faces,
                                           size_t         count_outlines,
                                           const CutMesh      &mesh,
                                           const ReductionMap &reduction_map,
                                           ConvertMap &v2v)
{
    // IMPROVE: use reduced count of faces and outlines
    size_t indices_size = faces.size();
    size_t vertices_size = (indices_size * 3 - count_outlines / 2) / 2;

    SurfaceCut sc;
    sc.indices.reserve(indices_size);
    sc.vertices.reserve(vertices_size);

    for (FI fi : faces) {
        //auto reduce = get_reduce_vertex(fi);
        HI hi = mesh.halfedge(fi);
        HI hi_end = hi;

        Vec3i its_face;
        // index into its_face
        int its_face_id = 0;
        bool exist_reduction = false;
        do {
            VI vi = mesh.source(hi);

            VI vi_r = reduction_map[vi];
            if (vi_r != vi) { 
                exist_reduction = true;
                vi = vi_r;
            }

            size_t index = v2v[vi];
            if (index == std::numeric_limits<SurfaceCut::Index>::max()) {
                index = sc.vertices.size();
                const auto &p = mesh.point(vi);
                // create vertex in result
                sc.vertices.emplace_back(p.x(), p.y(), p.z());
                v2v[vi] = index;
            }
            assert(index != std::numeric_limits<SurfaceCut::Index>::max());
            its_face[its_face_id++] = index;            
            hi = mesh.next(hi);
        } while (hi != hi_end);

        // prevent add reduced triangle
        if (exist_reduction && (
            its_face[0] == its_face[1] ||
            its_face[1] == its_face[2] || 
            its_face[2] == its_face[0]
            )) continue;

        sc.indices.emplace_back(std::move(its_face));
    }

    // reduce size with respect to reduction triangles
    sc.indices.shrink_to_fit();
    sc.vertices.shrink_to_fit();
    return sc;
}


SurfaceCut::CutContour priv::create_cut(const std::vector<HI> &outlines,
                                     const CutMesh         &mesh,
                                     const ReductionMap    &reduction_map,
                                     const ConvertMap      &v2v)
{
    using Index = SurfaceCut::Index;
    SurfaceCut::CutContour cut;
    SurfaceCut::CutContour unclosed_cut;
    for (HI hi : outlines) {
        VI vi_s = mesh.source(hi);
        VI vi_t = mesh.target(hi);
        // reduced vertex
        VI vi_s_r = reduction_map[vi_s];
        VI vi_t_r = reduction_map[vi_t];
        // is reduced edge?
        if (vi_s_r == vi_t || vi_t_r == vi_s) continue;

        // source vertex (from)
        Index vi_from = v2v[vi_s_r];
        assert(vi_from != std::numeric_limits<Index>::max());

        // target vertex (to)
        Index vi_to = v2v[vi_t_r];
        assert(vi_to != std::numeric_limits<Index>::max());

        std::vector<Index> *cut_move = nullptr;
        std::vector<Index> *cut_connect = nullptr;
        for (std::vector<Index> &cut : unclosed_cut) { 
            if (cut.back() != vi_from) continue;
            if (cut.front() == vi_to) {
                // cut closing
                cut_move = &cut;
            } else {
                cut_connect = &cut;
            }
            break;
        }
        if (cut_move != nullptr) {
            // index of closed cut
            size_t index = cut_move - &unclosed_cut.front();
            // move cut to result
            cut.emplace_back(std::move(*cut_move));
            // remove it from unclosed cut
            unclosed_cut.erase(unclosed_cut.begin() + index);                
        } else if (cut_connect != nullptr) {
            // try find tail to connect cut
            std::vector<Index> *cut_tail = nullptr;
            for (std::vector<Index> &cut : unclosed_cut) {
                if (cut.front() != vi_to) continue;
                cut_tail = &cut;
                break;
            }
            if (cut_tail != nullptr) {
                // index of tail
                size_t index = cut_tail - &unclosed_cut.front();
                // move to connect vector
                cut_connect->insert(cut_connect->end(),
                    make_move_iterator(cut_tail->begin()), 
                    make_move_iterator(cut_tail->end()));
                // remove tail from unclosed cut
                unclosed_cut.erase(unclosed_cut.begin() + index);
            } else {
                cut_connect->push_back(vi_to);
            }
        } else { // not found
            bool create_cut = true;
            // try to insert to front of cut
            for (std::vector<Index> &cut : unclosed_cut) {
                if (cut.front() != vi_to) continue;
                cut.insert(cut.begin(), vi_from);
                create_cut = false;
                break;
            }
            if (create_cut)
                unclosed_cut.emplace_back(std::vector{vi_from, vi_to});
        }
    }
    assert(unclosed_cut.empty());
    return cut;
}

priv::CutAOIs priv::create_cut_area_of_interests(const CutMesh    &mesh,
                                                 const ExPolygons &shapes,
                                                 FaceTypeMap &face_type_map)
{
    // IMPROVE: Create better heuristic for count.
    size_t faces_per_cut    = mesh.faces().size() / shapes.size();
    size_t outlines_per_cut = faces_per_cut / 2;
    size_t cuts_per_model   = shapes.size() * 2;

    CutAOIs result;
    result.reserve(cuts_per_model);

    // It is faster to use one queue for all cuts
    std::queue<FI> process;
    for (FI fi : mesh.faces()) {
        if (face_type_map[fi] != FaceType::inside) continue;

        CutAOI cut;
        std::vector<FI> &faces    = cut.first;
        std::vector<HI> &outlines = cut.second;

        // faces for one surface cut
        faces.reserve(faces_per_cut);
        // outline for one surface cut
        outlines.reserve(outlines_per_cut);

        assert(process.empty());
        // Process queue of faces to separate to surface_cut
        process.push(fi);
        collect_surface_data(process, faces, outlines, face_type_map, mesh);

        assert(!faces.empty());
        assert(!outlines.empty());
        result.emplace_back(std::move(cut));
    }
    return result;
}

void priv::filter_cuts(CutAOIs              &cuts,
                       const CutMesh        &mesh,
                       const ExPolygons     &shapes,
                       const Project        &projection,
                       const VertexShapeMap &vert_shape_map)
{
    auto get_point = [&shapes](const IntersectingElement &intersection) -> Point {
        assert(intersection.vertex_base != std::numeric_limits<uint32_t>::max());
        assert(intersection.point_index != std::numeric_limits<uint32_t>::max());
        size_t offset = 0;
        for (const ExPolygon &s : shapes) {
            if (offset == intersection.vertex_base) { 
                assert(s.contour.size() > intersection.point_index);
                return s.contour[intersection.point_index];
            }
            // *2 .. see description of IntersectingElement::vertex_base
            offset += 2*s.contour.size();
            assert(offset <= intersection.vertex_base);

            for (const Polygon &h : s.holes) {
                if (offset == intersection.vertex_base) {
                    assert(h.points.size() > intersection.point_index);
                    return h.points[intersection.point_index];
                }
                // *2 .. see description of IntersectingElement::vertex_base
                offset += 2*h.points.size();
                assert(offset <= intersection.vertex_base);
            }
        }

        // index is out of shape
        assert(false);
        return Point{};
    };

    struct CutIndex
    {
        // index in vector into cuts
        size_t cut_index = std::numeric_limits<size_t>::max();
        // vertex index inside of mesh
        VI vi;
    };
    size_t count = count_points(shapes);
    // each source point from shapes could has only one nearest projection
    std::vector<CutIndex> indices(count);

    // flags which cut is not first
    std::vector<bool> del_cuts(cuts.size(), false);

    // check whether vertex is behind another cut
    auto is_behind = [&vert_shape_map, &indices, &del_cuts, &get_point,
                      &projection, &mesh]
                      (VI vi, size_t cut_index) -> bool {
        const IntersectingElement *i = vert_shape_map[vi];

        // Is vertex made by corefine?
        if (i == nullptr) return false;
        
        assert(i->vertex_base != std::numeric_limits<uint32_t>::max());
        assert(i->vertex_base%2 == 0);
        assert(i->point_index != std::numeric_limits<uint32_t>::max());
        assert(i->attr != (unsigned char)IntersectingElement::Type::undefined);

        // Use only straigh edge
        if (i->get_type() != IntersectingElement::Type::edge_1)
            return false;
        

        size_t index = i->vertex_base/2 + i->point_index;
        CutIndex &ci = indices[index];

        // is first cut for vertex OR
        // is remembred cut is deleted?
        if (ci.cut_index == std::numeric_limits<size_t>::max() || 
            del_cuts[ci.cut_index] ) {
            ci.cut_index = cut_index;
            ci.vi         = vi;
            return false;
        }

        if (ci.cut_index == cut_index) { 
            // In one connected triangles area are more points 
            // with same source point from text contour
            //assert(ci.vi == vi);
            return false;
        }

        // compare distances of vertices
        Point p = get_point(*i);
        Vec3f source_point = projection.project(p).first;
        const auto &prev = mesh.point(ci.vi);
        Vec3f prev_point(prev.x(), prev.y(), prev.z());
        float prev_sq_norm = (source_point - prev_point).squaredNorm();

        const auto &act = mesh.point(vi);
        Vec3f act_point(act.x(), act.y(), act.z());
        float act_sq_norm = (source_point - act_point).squaredNorm();
        
        if (act_sq_norm > prev_sq_norm) {
            del_cuts[cut_index] = true;
            return true;
        }

        // previous cut is behind actual one
        del_cuts[ci.cut_index] = true;
        ci.cut_index = cut_index;
        ci.vi = vi;
        return false;
    };

    // filter top one cuts
    for (const CutAOI &cut : cuts) {
        size_t cut_index = &cut - &cuts.front();
        const std::vector<HI> &outlines = cut.second;
        for (HI hi : outlines) {
            if (is_behind(mesh.source(hi), cut_index) ||
                is_behind(mesh.target(hi), cut_index))
                break;
        }
    }

    // remove flagged cuts
    for (size_t i = del_cuts.size(); i > 0; --i) { 
        size_t index = i - 1;
        if (del_cuts[index])
            cuts.erase(cuts.begin() + index);
    }
}


SurfaceCuts priv::create_surface_cuts(const CutAOIs      &cuts,
                                      const CutMesh      &mesh,
                                      const ReductionMap &reduction_map,
                                      ConvertMap         &convert_map)
{
    // initialize convert_map to MAX values
    for (VI vi : mesh.vertices())
        convert_map[vi] = std::numeric_limits<SurfaceCut::Index>::max();

    SurfaceCuts result;
    for (const CutAOI &cut : cuts) {
        const std::vector<FI>& faces = cut.first;
        const std::vector<HI> &outlines = cut.second;
        
        // convert_map could be used separately for each surface cut. 
        // But it is moore faster to use one memory allocation for them all.
        SurfaceCut sc = create_index_triangle_set(faces, outlines.size(), mesh, reduction_map, convert_map);
        
        // connect outlines
        sc.contours = create_cut(outlines, mesh, reduction_map, convert_map);
        result.emplace_back(std::move(sc));
    }
    return result;
}

// only for debug
void priv::store(CutMesh &mesh, const FaceTypeMap &face_type_map, const std::string& file)
{
    auto face_colors = mesh.add_property_map<priv::FI, CGAL::Color>("f:color").first;    
    for (FI fi : mesh.faces()) { 
        auto &color = face_colors[fi];
        switch (face_type_map[fi]) {
        case FaceType::inside: color = CGAL::Color{255, 0, 0}; break;
        case FaceType::inside_: color = CGAL::Color{150, 0, 0}; break;
        case FaceType::outside: color = CGAL::Color{255, 0, 255}; break;
        case FaceType::not_constrained: color = CGAL::Color{0, 255, 0}; break;
        default: color = CGAL::Color{127, 127, 127};
        }
    }
    CGAL::IO::write_OFF(file, mesh);
    mesh.remove_property_map(face_colors);
}

void priv::store(CutMesh &mesh, const ReductionMap &reduction_map, const std::string& file)
{
    auto vertex_colors = mesh.add_property_map<priv::VI, CGAL::Color>("v:color").first;    
    // initialize to gray color
    for (VI vi: mesh.vertices())
        vertex_colors[vi] = CGAL::Color{127, 127, 127};

    for (VI reduction_from : mesh.vertices()) {
        VI reduction_to = reduction_map[reduction_from];
        if (reduction_to != reduction_from) {
            vertex_colors[reduction_from] = CGAL::Color{255, 0, 0};
            vertex_colors[reduction_to] = CGAL::Color{0, 0, 255};
        }
    }
    CGAL::IO::write_OFF(file, mesh);
    mesh.remove_property_map(vertex_colors);
}

void priv::store(const SurfaceCuts &cut, const std::string &file_prefix) {
    for (auto &c : cut) {
        size_t index = &c - &cut.front();
        std::string file  = file_prefix + std::to_string(index) + ".obj";
        its_write_obj(c, file.c_str());  
    }
}
