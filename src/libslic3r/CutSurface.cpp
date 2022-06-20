#include "CutSurface.hpp"

/// model.off - CGAL model created from index_triangle_set
/// shape.off - CGAL model created from shapes
/// constrained.off - Visualization of inside and outside triangles
///    Green - not along constrained edge
///    Red - sure that are inside
///    Purple - sure that are outside
/// filled.off - flood fill green triangles inside of red area
///            - Same meaning of color as constrained
/// reduction.off - Visualization of reduced and non-reduced Vertices
/// aois/cutAOI{N}.obj - Cuted Area of interest from corefined model
/// cuts/cut{N}.obj - Filtered surface cuts + Reduced vertices made by e2 (text_edge_2)
#define DEBUG_OUTPUT_DIR std::string("C:/data/temp/cutSurface/")

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

using P3 = CGAL::Epick::Point_3;

using Project = Emboss::IProjection;
using Project3f = Emboss::IProject3f;

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
    // identify source point in shapes
    uint32_t shape_point_index{std::numeric_limits<uint32_t>::max()};

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
/// set true to skip map for indicies
/// </summary>
/// <param name="skip_indicies">Flag to convert triangle to cgal</param>
/// <param name="its">model</param>
/// <param name="projection">direction</param>
void set_skip_for_outward_projection(std::vector<bool> &skip_indicies,
                                     const indexed_triangle_set &its,
                                     const Project3f   &projection);

/// <summary>
/// Set true for indicies outward and almost parallel together.
/// Note: internally calculate normals
/// </summary>
/// <param name="skip_indicies">Flag to convert triangle to cgal</param>
/// <param name="its">model</param>
/// <param name="projection">Direction to measure angle</param>
/// <param name="max_angle">Maximal allowed angle between opposit normal and projection direction [in DEG]</param>
void set_skip_by_angle(std::vector<bool>          &skip_indicies,
                       const indexed_triangle_set &its,
                       const Project3f            &projection,
                       double                      max_angle = 89.);

/// <summary>
/// Set true for indices out of area of interest
/// </summary>
/// <param name="skip_indicies">Flag to convert triangle to cgal</param>
/// <param name="its">model</param>
/// <param name="projection">Convert 2d point to pair of 3d points</param>
/// <param name="shapes_bb">2d bounding box define AOI</param>
void set_skip_for_out_of_aoi(std::vector<bool>          &skip_indicies,
                             const indexed_triangle_set &its,
                             const Project              &projection,
                             const BoundingBox          &shapes_bb);

/// <summary>
/// Convert triangle mesh model to CGAL Surface_mesh
/// Filtrate out opposite triangles
/// Add property map for source face index
/// </summary>
/// <param name="its">Model</param>
/// <param name="skip_indicies">Flags that triangle should be skiped</param>
/// <returns>CGAL mesh - half edge mesh</returns>
CutMesh to_cgal(const indexed_triangle_set &its,
                const std::vector<bool> &skip_indicies);

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

/// <summary>
/// Identify contour (or hole) point from ExPolygons
/// </summary>
struct ShapePointId
{
    // index of ExPolygons
    uint32_t expolygons_index;
    // index of Polygon
    uint32_t polygon_index;
    // index of point in polygon
    uint32_t point_index;
};
/// <summary>
/// Keep conversion from ShapePointId to Index and vice versa
/// ShapePoint .. contour(or hole) poin from ExPolygons
/// Index      .. continous number
/// </summary>
class ShapePoint2index
{
    std::vector<std::vector<uint32_t>> m_offsets;
    // for check range of index
    uint32_t                           m_count;
public:
    ShapePoint2index(const ExPolygons &shapes);
    uint32_t     calc_index(const ShapePointId &id) const;
    ShapePointId calc_id(uint32_t index) const;
    uint32_t get_count() const;
};

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
    // face, inside but almost in direction of projection
    inside_parallel,
    // face outside of the cutted shape
    outside,
    // face without constrained edge (In or Out)
    not_constrained,

    // Helper flag that inside was processed
    inside_processed
};
using FaceTypeMap = CutMesh::Property_map<FI, FaceType>;
/// <summary>
/// Face with constrained edge are inside/outside by type of intersection
/// Other set to not_constrained(still it could be inside/outside)
/// </summary>
/// <param name="face_type_map">[Output] property map with type of faces</param>
/// <param name="mesh">Mesh to process</param>
/// <param name="vertex_shape_map">Keep information about source of created vertex</param>
/// <param name="ecm">Dynamic Edge Constrained Map of bool</param>
/// <param name="shape_mesh">Vertices of mesh made by shapes</param>
/// <param name="shape2index">Convert index to shape point from ExPolygons</param>
void set_face_type(FaceTypeMap          &face_type_map,
                   const CutMesh        &mesh,
                   const VertexShapeMap &vertex_shape_map,
                   const EcmType        &ecm,
                   const CutMesh          &shape_mesh,
                   const ShapePoint2index &shape2index);

void set_almost_parallel_type(FaceTypeMap              &face_type_map,
                              const CutMesh            &mesh,
                              const Project3f &projection);

/// <summary>
/// Check if face is almost parallel
/// </summary>
/// <param name="fi">Index of triangle (Face index)</param>
/// <param name="mesh">Must contain fi</param>
/// <param name="projection">Direction of projection</param>
/// <param name="threshold">Value for cos(alpha), must be greater than zero, 
/// where alpha is minimal angle between projection direction and face normal</param>
/// <returns>True when Triangle is almost parallel with direction of projection</returns>
bool is_almost_parallel(FI                        fi,
                        const CutMesh            &mesh,
                        const Project3f &projection,
                        float threshold = static_cast<float>(std::cos(80 * M_PI / 180)));

/// <summary>
/// Change FaceType from not_constrained to inside
/// For neighbor(or neighbor of neighbor of ...) of inside triangles.
/// Process only not_constrained triangles
/// </summary>
/// <param name="mesh">Corefined mesh</param>
/// <param name="face_type_map">In/Out map with faces type</param>
void flood_fill_inner(const CutMesh &mesh, FaceTypeMap &face_type_map);

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
/// To select correct area
/// </summary>
struct ProjectionDistance
{
    // index of CutAOI
    uint32_t aoi_index = std::numeric_limits<uint32_t>::max();

    // index of half edge in AOI
    uint32_t hi_index = std::numeric_limits<uint32_t>::max();

    // signed distance to projection 
    float distance = std::numeric_limits<float>::max();    
};
// addresed by ShapePoint2index
using ProjectionDistances =  std::vector<ProjectionDistance>;

/// <summary>
/// Calculate distances from CutAOI contour points to ProjectionOrigin
/// </summary>
/// <param name="cuts">AOIs</param>
/// <param name="mesh">Vertices position</param>
/// <param name="shapes_points">Count of points in shapes</param>
/// <param name="shapes_mesh">Mesh created by shapes</param>
/// <param name="source_point">Origin of projection</param>
/// <param name="vert_shape_map">Know source of new vertices</param>
/// <param name="shape_point_2_index">Convert shapepoint to index</param>
/// <returns>distances</returns>
std::vector<ProjectionDistances> create_distances(
    const CutAOIs        &cuts,
    const CutMesh        &mesh,
    uint32_t              shapes_points,
    const CutMesh        &shapes_mesh,
    float                 projection_ratio,
    const VertexShapeMap &vert_shape_map);

/// <summary>
/// Select distances in similar depth between expolygons
/// </summary>
/// <param name="distances">All distances</param>
/// <param name="shapes">Vector of letters</param>
/// <param name="shapes_bb">2d Bound of shapes</param>
/// <param name="shape_point_2_index">Convert index to addresss inside of shape</param>
/// <returns>Best projection distances</returns>
ProjectionDistances choose_best_distance(
    const std::vector<ProjectionDistances> &distances,
    const ExPolygons                      &shapes,
    const BoundingBox                     &shapes_bb,
    const ShapePoint2index                &shape_point_2_index);

/// <summary>
/// Merge 2 cuts together, cut off inner part
/// </summary>
/// <param name="cut1">[In/Out] surface cut</param>
/// <param name="cut2">Cut to intersect with</param>
/// <param name="mesh">Source of surface</param>
/// <returns>True when merge otherwise False</returns>
bool merge_cut(CutAOI &cut1, const CutAOI &cut2, const CutMesh &mesh);

/// <summary>
/// Merge cuts together
/// </summary>
/// <param name="cuts">[In/Out] cutted surface of model</param>
/// <param name="mesh">Source of surface</param>
/// <param name="use_cut_indices">Filtered indicies of Cut from best projection distances</param>
void merge_cuts(CutAOIs       &cuts,
                const CutMesh &mesh,
                const std::vector<size_t> &use_cut_indices);

/// <summary>
/// Filter out cuts which are behind another.
/// Prevent overlapping embossed shape in space.
/// </summary>
/// <param name="cuts">AOIs</param>
/// <param name="mesh">triangle model</param>
/// <param name="shapes">2d cutted shapes</param>
/// <param name="shape_point_2_index">2d cutted shapes</param>
/// <param name="projection">Projection from 2d to 3d</param>
/// <param name="vert_shape_map">Identify source of intersection</param>
void filter_cuts(CutAOIs              &cuts,
                 const CutMesh        &mesh,
                 const ExPolygons       &shapes,
                 const ShapePoint2index &shape_point_2_index,
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
/// NOTE: Modify in function: inside -> inside_processed</param>
/// <param name="mesh">mesh to process</param>
void collect_surface_data(std::queue<FI>  &process,
                          std::vector<FI> &faces,
                          std::vector<HI> &outlines,
                          FaceTypeMap     &face_type_map,
                          const CutMesh   &mesh);

/// <summary>
/// Check if contour contain at least 2 unique source contour points
/// </summary>
/// <param name="outlines">Vector of half edge define outline</param>
/// <param name="vert_shape_map">For corefine edge vertices know source of intersection</param>
/// <param name="mesh">Triangle half edge mesh</param>
/// <param name="min_count">Minimal count of unique source points creating outline</param>
/// <returns>True when outline is made by minimal or greater count of unique source point otherwise FALSE</returns>
bool has_minimal_contour_points(const std::vector<HI> &outlines,
                                const VertexShapeMap  &vert_shape_map,
                                const CutMesh         &mesh,
                                size_t                 min_count = 2);

/// <summary>
/// Check orientation of triangle
/// </summary>
/// <param name="fi">triangle</param>
/// <param name="mesh">Triangle half edge mesh</param>
/// <param name="projection">Direction to check</param>
/// <returns>True when triangle normal is toward projection otherwise FALSE</returns>
bool is_toward_projection(FI                        fi,
                          const CutMesh            &mesh,
                          const Project3f &projection);

/// <summary>
/// Check orientation of triangle defined by vertices a, b, c in CCW order
/// </summary>
/// <param name="a">Trienagle vertex</param>
/// <param name="b">Trienagle vertex</param>
/// <param name="c">Trienagle vertex</param>
/// <param name="projection">Direction to check</param>
/// <returns>True when triangle normal is toward projection otherwise FALSE</returns>
bool is_toward_projection(const Vec3f        &a,
                          const Vec3f        &b,
                          const Vec3f        &c,
                          const Project3f &projection);
/// <summary>
/// Check orientation of triangle
/// </summary>
/// <param name="t">triangle indicies into vertices vector</param>
/// <param name="vertices">model vertices</param>
/// <param name="projection">Direction to check</param>
/// <returns>True when triangle normal is toward projection otherwise FALSE</returns>
bool is_toward_projection(const stl_triangle_vertex_indices &t,
                          const std::vector<stl_vertex>     &vertices,
                          const Project3f          &projection);

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

#ifdef DEBUG_OUTPUT_DIR
indexed_triangle_set create_indexed_triangle_set(const std::vector<FI> &faces,
                                                 const CutMesh         &mesh);

/// <summary>
/// Debug purpose store of mesh with colored face by face type
/// </summary>
/// <param name="mesh">Input mesh, could add property color
/// NOTE: Not const because need to [optionaly] append color property map</param>
/// <param name="face_type_map">Color source</param>
/// <param name="file">File to store</param>
void store(CutMesh &mesh, const FaceTypeMap &face_type_map, const std::string &file);
void store(CutMesh &mesh, const ReductionMap &reduction_map, const std::string &file);
void store(const CutAOIs &aois, const CutMesh &mesh, const std::string &dir);
void store(const Vec3f &vertex, const Vec3f &normal, const std::string &file, float size = 2.f);
void store(const ProjectionDistances &pds, const CutAOIs &aois, const CutMesh &mesh, const std::string &file, float width = 0.2f/* [in mm] */);
void store(const SurfaceCuts &cut, const std::string &dir);
#endif // DEBUG_OUTPUT_DIR

} // namespace privat


SurfaceCut Slic3r::cut_surface(const indexed_triangle_set &model,
                               const ExPolygons           &shapes,
                               const Emboss::IProjection  &projection,
                               float                      projection_ratio)
{
    if (model.empty() || shapes.empty() ) return {};
#ifdef DEBUG_OUTPUT_DIR
    its_write_obj(model, (DEBUG_OUTPUT_DIR + "model_input.obj").c_str()); // only debug
#endif // DEBUG_OUTPUT_DIR

    std::vector<bool> skip_indicies(model.indices.size(), {false});
    // cut out of bounding box triangles
    BoundingBox shapes_bb = get_extents(shapes);
    priv::set_skip_for_out_of_aoi(skip_indicies, model, projection, shapes_bb);
    // cut out opposit triangles
    //priv::set_skip_for_outward_projection(skip_indicies, model, projection);
    priv::set_skip_by_angle(skip_indicies, model, projection);

    priv::CutMesh cgal_model = priv::to_cgal(model, skip_indicies);

#ifdef DEBUG_OUTPUT_DIR
    CGAL::IO::write_OFF(DEBUG_OUTPUT_DIR + "model.off", cgal_model); // only debug
#endif // DEBUG_OUTPUT_DIR

    std::string edge_shape_map_name = "e:IntersectingElement";
    std::string face_shape_map_name = "f:IntersectingElement";
    priv::CutMesh cgal_shape = priv::to_cgal(shapes, projection, edge_shape_map_name, face_shape_map_name);

#ifdef DEBUG_OUTPUT_DIR
    CGAL::IO::write_OFF(DEBUG_OUTPUT_DIR + "shape.off", cgal_shape); // only debug
#endif // DEBUG_OUTPUT_DIR

    auto edge_shape_map = cgal_shape.property_map<priv::EI, priv::IntersectingElement>(edge_shape_map_name).first;    
    auto face_shape_map = cgal_shape.property_map<priv::FI, priv::IntersectingElement>(face_shape_map_name).first;
    
    std::string vert_shape_map_name = "v:IntersectingElement";
    // pointer to edge or face shape_map
    priv::VertexShapeMap vert_shape_map = cgal_model.add_property_map<priv::VI, const priv::IntersectingElement*>(vert_shape_map_name).first;
    
    // detect anomalities in visitor.
    bool is_valid = true;
    // create anotation visitor - Must be copyable
    priv::ShapePoint2index shape_point_2_index(shapes); 

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
    priv::set_face_type(face_type_map, cgal_model, vert_shape_map, ecm,
                        cgal_shape, shape_point_2_index);
#ifdef DEBUG_OUTPUT_DIR
    priv::store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "constrained.off"); // only debug
#endif // DEBUG_OUTPUT_DIR

    // It is neccesary when almost parallel face are contained in projection
    // priv::set_almost_parallel_type(face_type_map, cgal_model, projection);
//#ifdef DEBUG_OUTPUT_DIR
//    priv::store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "constrainedWithAlmostParallel.off"); // only debug
//#endif // DEBUG_OUTPUT_DIR
    
    // flood fill the other faces inside the region.
    priv::flood_fill_inner(cgal_model, face_type_map);

#ifdef DEBUG_OUTPUT_DIR
    priv::store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "filled.off"); // only debug
#endif // DEBUG_OUTPUT_DIR

    std::string vertex_reduction_map_name = "v:reduction";
    priv::ReductionMap vertex_reduction_map = cgal_model.add_property_map<priv::VI, priv::VI>(vertex_reduction_map_name).first;
    priv::create_reduce_map(vertex_reduction_map, cgal_model, face_type_map, vert_shape_map);

#ifdef DEBUG_OUTPUT_DIR
    priv::store(cgal_model, vertex_reduction_map, DEBUG_OUTPUT_DIR + "reduction.off"); // only debug
#endif // DEBUG_OUTPUT_DIR

    // IMPROVE: AOIs area could be created during flood fill
    priv::CutAOIs cutAOIs = create_cut_area_of_interests(cgal_model, shapes, face_type_map);

#ifdef DEBUG_OUTPUT_DIR
    priv::store(cutAOIs, cgal_model, DEBUG_OUTPUT_DIR + "aois/"); // only debug
#endif // DEBUG_OUTPUT_DIR

    // calc distance to projection for all outline points of cutAOI(shape)
    // it is used for distiguish the top one 
    uint32_t shapes_points = shape_point_2_index.get_count();

    // for each point collect all projection distances
    std::vector<priv::ProjectionDistances> distances =
        priv::create_distances(cutAOIs, cgal_model, shapes_points, cgal_shape, projection_ratio, vert_shape_map);

    // NOTE: it will be fine to calc AOIs range,
    // not only outline but all vertices in direction of emboss - faster check on intersection
    
#ifdef DEBUG_OUTPUT_DIR
    auto [front,back] = projection.create_front_back(shapes_bb.center());
    Vec3f diff = back - front;
    Vec3f pos = front + diff*projection_ratio;
    priv::store(pos, diff.normalized(), DEBUG_OUTPUT_DIR + "projection_center.obj"); // only debug
#endif // DEBUG_OUTPUT_DIR

    // for each point select best projection
    priv::ProjectionDistances best_projection =
        priv::choose_best_distance(distances, shapes, shapes_bb, shape_point_2_index);    
#ifdef DEBUG_OUTPUT_DIR
    priv::store(best_projection, cutAOIs, cgal_model, DEBUG_OUTPUT_DIR + "best_projection.obj"); // only debug
#endif // DEBUG_OUTPUT_DIR

    std::vector<bool> is_best_cut(cutAOIs.size(), {false});
    for (const priv::ProjectionDistance &d : best_projection)
        is_best_cut[d.aoi_index] = true;
    std::vector<size_t> best_cut_indices;
    for (size_t i = 0; i < cutAOIs.size(); ++i)
        if (is_best_cut[i]) best_cut_indices.push_back(i);

    // cut off part + filtrate cutAOIs
    priv::merge_cuts(cutAOIs, cgal_model, best_cut_indices);
#ifdef DEBUG_OUTPUT_DIR
    priv::store(cutAOIs, cgal_model, DEBUG_OUTPUT_DIR + "merged-aois/");
    // only debug
#endif // DEBUG_OUTPUT_DIR
    
    // Filter out NO top one cuts
    priv::filter_cuts(cutAOIs, cgal_model, shapes,
                      shape_point_2_index, projection, vert_shape_map);

    // conversion map between vertex index in cgal_model and indices in result
    // used instead of std::map
    std::string vertec_convert_map_name = "v:convert";
    priv::ConvertMap vertex_convert_map = cgal_model.add_property_map<priv::VI, SurfaceCut::Index>(vertec_convert_map_name).first;    
    SurfaceCuts result = priv::create_surface_cuts(cutAOIs, cgal_model, vertex_reduction_map, vertex_convert_map);
    
#ifdef DEBUG_OUTPUT_DIR
    priv::store(result, DEBUG_OUTPUT_DIR + "cuts/"); // only debug
#endif // DEBUG_OUTPUT_DIR
    
    // TODO: fill skipped source triangles to surface of cut
    return merge(std::move(result));
}

indexed_triangle_set Slic3r::cuts2model(const SurfaceCuts        &cuts,
                                        const priv::Project3f &projection)
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

indexed_triangle_set Slic3r::cut2model(const SurfaceCut         &cut,
                                       const priv::Project3f &projection)
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
        assert(i.x() + back_offset < result.vertices.size());
        assert(i.y() + back_offset < result.vertices.size());
        assert(i.z() + back_offset < result.vertices.size());
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

bool priv::is_toward_projection(FI                        fi,
                                const CutMesh            &mesh,
                                const Project3f &projection)
{
    HI hi = mesh.halfedge(fi);

    const P3 &a = mesh.point(mesh.source(hi));
    const P3 &b = mesh.point(mesh.target(hi));
    const P3 &c = mesh.point(mesh.target(mesh.next(hi)));
        
    Vec3f a_(a.x(), a.y(), a.z());
    Vec3f p_ = projection.project(a_);
    P3 p{p_.x(), p_.y(), p_.z()};

    return CGAL::orientation(a, b, c, p) == CGAL::POSITIVE;
}

bool priv::is_toward_projection(const stl_triangle_vertex_indices &t,
                                const std::vector<stl_vertex> &vertices,
                                const Project3f &projection)
{
    return is_toward_projection(vertices[t[0]], vertices[t[1]],
                                vertices[t[2]], projection);
}

bool priv::is_toward_projection(const Vec3f              &a,
                                const Vec3f              &b,
                                const Vec3f              &c,
                                const Project3f &projection)
{
    P3 cgal_a(a.x(), a.y(), a.z());
    P3 cgal_b(b.x(), b.y(), b.z());
    P3 cgal_c(c.x(), c.y(), c.z());

    Vec3f p = projection.project(a);
    P3    cgal_p{p.x(), p.y(), p.z()};

    // is_toward_projection
    return CGAL::orientation(cgal_a, cgal_b, cgal_c, cgal_p) ==
           CGAL::POSITIVE;
}

void priv::set_skip_for_outward_projection(std::vector<bool> &skip_indicies,
                                           const indexed_triangle_set &its,
                                           const Project3f &projection)
{
    assert(skip_indicies.size() == its.indices.size());
    for (const auto &t : its.indices) {
        size_t index = &t - &its.indices.front();
        if (skip_indicies[index]) continue;
        if (is_toward_projection(t, its.vertices, projection)) continue;
        skip_indicies[index] = true;
    }
}

void priv::set_skip_by_angle(std::vector<bool>          &skip_indicies,
                             const indexed_triangle_set &its,
                             const Project3f            &projection,
                             double                      max_angle)
{
    float threshold = static_cast<float>(cos(max_angle / 180. * M_PI));
    assert(skip_indicies.size() == its.indices.size());
    for (const stl_triangle_vertex_indices& face : its.indices) {
        size_t index = &face - &its.indices.front();
        if (skip_indicies[index]) continue;
        Vec3f n = its_face_normal(its, face);
        const Vec3f v = its.vertices[face[0]];
        // Improve: For Orthogonal Projection it is same for each vertex
        Vec3f projected = projection.project(v);
        Vec3f project_dir = projected - v;
        project_dir.normalize();
        float cos_alpha = project_dir.dot(n);
        if (cos_alpha > threshold) continue;
        skip_indicies[index] = true;
    }
}

void priv::set_skip_for_out_of_aoi(std::vector<bool>          &skip_indicies,
                                   const indexed_triangle_set &its,
                                   const Project              &projection,
                                   const BoundingBox          &shapes_bb)
{
    assert(skip_indicies.size() == its.indices.size());

    //   1`*----* 2`
    //    /  2 /|
    // 1 *----* |
    //   |    | * 3`
    //   |    |/
    // 0 *----* 3
    //////////////////
    std::array<std::pair<Vec3f, Vec3f>, 4> bb;
    int index = 0;
    for (Point v :
         {shapes_bb.min, Point{shapes_bb.min.x(), shapes_bb.max.y()},
          shapes_bb.max, Point{shapes_bb.max.x(), shapes_bb.min.y()}})
        bb[index++] = projection.create_front_back(v);

    // define planes to test
    // 0 .. under
    // 1 .. left
    // 2 .. above
    // 3 .. right
    size_t prev_i = 3;
    std::array<std::pair<Vec3d, Vec3d>, 4> point_normals;
    for (size_t i = 0; i < 4; i++) {
        const Vec3f &p1 = bb[i].first;
        const Vec3f &p2 = bb[i].second;
        const Vec3f &p3 = bb[prev_i].first;
        prev_i = i;

        Vec3d v1 = (p2 - p1).cast<double>();
        v1.normalize();
        Vec3d v2 = (p3 - p1).cast<double>();
        v2.normalize();

        Vec3d normal = v2.cross(v1);
        normal.normalize(); 

        point_normals[i] = {p1.cast<double>(), normal};
    }
    // same meaning as point normal
    std::array<std::vector<bool>, 4> is_on_sides;
    for (size_t side = 0; side < 4; side++)
        is_on_sides[side] = std::vector<bool>(its.vertices.size(), {false});

    auto is_out_of = [&point_normals](int side, const Vec3d &v) -> bool {
        const auto &[p, n] = point_normals[side];
        double signed_distance = (v - p).dot(n);
        return signed_distance > 1e-5;
    };

    // inspect all vertices when it is out of bounding box
    for (size_t i = 0; i < its.vertices.size(); i++) {
        Vec3d v = its.vertices[i].cast<double>();
        // under + above
        for (int side : {0, 2}) {
            if (is_out_of(side, v)) {
                is_on_sides[side][i] = true;
                break;
            }
        }
        // left + right
        for (int side : {1, 3}) {
            if (is_out_of(side, v)) {
                is_on_sides[side][i] = true;
                break;
            }        
        }
    }

    auto is_all_on_one_side = [is_on_sides](const Vec3i &t) -> bool {
        for (size_t side = 0; side < 4; side++) {
            bool result = true;
            for (auto vi : t){ 
                if (!is_on_sides[side][vi]) { 
                    result = false; 
                    break;
                }
            }
            if (result) return true;
        }
        return false;
    };

    // inspect all triangles, when it is out of bounding box
    for (size_t i = 0; i < its.indices.size(); i++) { 
        if (skip_indicies[i]) continue;
        const auto& t = its.indices[i];
        if (is_all_on_one_side(t)) 
            skip_indicies[i] = true;
    }
}

priv::CutMesh priv::to_cgal(const indexed_triangle_set &its,
                            const std::vector<bool>    &skip_indicies)
{
    const std::vector<stl_vertex>                  &vertices = its.vertices;
    const std::vector<stl_triangle_vertex_indices> &indices = its.indices;

    std::vector<bool> use_indices(indices.size(), {false});
    std::vector<bool> use_vetices(vertices.size(), {false});

    size_t vertices_count = 0;
    size_t faces_count    = 0;
    size_t edges_count    = 0;

    for (const auto &t : indices) {     
        size_t index = &t - &indices.front();
        if (skip_indicies[index]) continue;        
        ++faces_count;
        use_indices[index] = true;
        size_t count_used_vertices = 0;
        for (const auto vi : t) {
            if (!use_vetices[vi]) {
                ++vertices_count;
                use_vetices[vi] = true;
            } else {
                ++count_used_vertices;
            }
        }
        switch (count_used_vertices) {
        case 3: break; // all edges are already counted
        case 2: edges_count += 2; break;
        case 1:
        case 0: edges_count += 3; break;
        default: assert(false);
        }        
    }
    assert(vertices_count <= vertices.size());
    assert(edges_count <= (indices.size() * 3) / 2);
    assert(faces_count <= indices.size());

    CutMesh result;
    result.reserve(vertices_count, edges_count, faces_count);

    std::vector<size_t> to_filtrated_vertices_index(vertices.size());
    size_t filtrated_vertices_index = 0;
    for (size_t i = 0; i < vertices.size(); ++i) 
        if (use_vetices[i]) { 
            to_filtrated_vertices_index[i] = filtrated_vertices_index;
            ++filtrated_vertices_index;
        }

    for (const stl_vertex& v : vertices) {
        if (!use_vetices[&v - &vertices.front()]) continue;
        result.add_vertex(CutMesh::Point{v.x(), v.y(), v.z()});
    }

    for (const stl_triangle_vertex_indices &f : indices) {
        if (!use_indices[&f - &indices.front()]) continue;
        result.add_face(static_cast<VI>(to_filtrated_vertices_index[f[0]]), 
                        static_cast<VI>(to_filtrated_vertices_index[f[1]]),
                        static_cast<VI>(to_filtrated_vertices_index[f[2]]));
    }
    
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
            auto p  = projection.create_front_back(p2);
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

        uint32_t contour_index = static_cast<uint32_t>(num_vertices_old / 2);
        for (int32_t i = 0; i < int32_t(indices.size()); i += 2) {
            bool    is_first  = i == 0;
            bool    is_last   = size_t(i + 2) >= indices.size();
            int32_t j = is_last ? 0 : (i + 2);
            
            auto fi1 = result.add_face(indices[i], indices[j], indices[i + 1]);
            auto ei1 = find_edge(fi1, indices[i + 1], indices[i]);
            auto ei2 = find_edge(fi1, indices[j], indices[i + 1]);
            auto fi2 = result.add_face(indices[j], indices[j + 1], indices[i + 1]);
            IntersectingElement element {contour_index, (unsigned char)IntersectingElement::Type::undefined};
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

void priv::set_face_type(FaceTypeMap            &face_type_map,
                         const CutMesh          &mesh,
                         const VertexShapeMap   &vertex_shape_map,
                         const EcmType          &ecm,
                         const CutMesh          &shape_mesh,
                         const ShapePoint2index &shape2index)
{
    auto get_face_type = [&mesh, &shape_mesh, &vertex_shape_map, &shape2index](HI hi) -> FaceType {
        VI vi_from = mesh.source(hi);
        VI vi_to   = mesh.target(hi);
        // This face has a constrained edge.
        const IntersectingElement &shape_from = *vertex_shape_map[vi_from];
        const IntersectingElement &shape_to   = *vertex_shape_map[vi_to];
        assert(shape_from.shape_point_index != std::numeric_limits<uint32_t>::max());
        assert(shape_from.attr != (unsigned char) IntersectingElement::Type::undefined);
        assert(shape_to.shape_point_index != std::numeric_limits<uint32_t>::max());
        assert(shape_to.attr != (unsigned char) IntersectingElement::Type::undefined);

        bool is_inside = false;
        // index into contour
        uint32_t                  i_from    = shape_from.shape_point_index;
        uint32_t                  i_to      = shape_to.shape_point_index;
        IntersectingElement::Type type_from = shape_from.get_type();
        IntersectingElement::Type type_to   = shape_to.get_type();
        if (i_from == i_to && type_from == type_to) {
            // intersecting element must be face
            assert(type_from == IntersectingElement::Type::face_1 ||
                   type_from == IntersectingElement::Type::face_2);

            // count of vertices is twice as count of point in the contour
            uint32_t i = i_from * 2;
            // j is next contour point in vertices
            uint32_t j = i + 2;
            if (shape_from.is_last()) {
                ShapePointId point_id = shape2index.calc_id(i_from);
                point_id.point_index  = 0;
                j = shape2index.calc_index(point_id)*2;
            }

            // opposit point(in triangle face) to edge
            const auto &p = mesh.point(mesh.target(mesh.next(hi)));

            // abc is source triangle face
            auto abcp = type_from == IntersectingElement::Type::face_1 ?
                            CGAL::orientation(shape_mesh.point(VI(i)),
                                              shape_mesh.point(VI(i + 1)),
                                              shape_mesh.point(VI(j)), p) :
                            // type_from == IntersectingElement::Type::face_2
                            CGAL::orientation(shape_mesh.point(VI(j)),
                                              shape_mesh.point(VI(i + 1)),
                                              shape_mesh.point(VI(j + 1)), p);
            is_inside = abcp == CGAL::POSITIVE;
        } else if (i_from < i_to || (i_from == i_to && type_from < type_to)) {
            bool is_last = shape_to.is_last() && shape_from.is_first();
            // check continuity of indicies
            assert(i_from == i_to || is_last || (i_from + 1) == i_to);
            if (!is_last) is_inside = true;
        } else {
            assert(i_from > i_to || (i_from == i_to && type_from > type_to));
            bool is_last = shape_to.is_first() && shape_from.is_last();
            // check continuity of indicies
            assert(i_from == i_to || is_last || (i_to + 1) == i_from);
            if (is_last) is_inside = true;
        }
        return (is_inside) ? FaceType::inside : FaceType::outside;
    };

    for (const FI& fi : mesh.faces()) {
        FaceType face_type = FaceType::not_constrained;
        HI hi_end = mesh.halfedge(fi);
        HI hi     = hi_end;
        do {
            // is edge new created - constrained?
            if (get(ecm, mesh.edge(hi))) {
                face_type = get_face_type(hi);
                break;
            }
            // next half edge index inside of face
            hi = mesh.next(hi);
        } while (hi != hi_end);
        face_type_map[fi] = face_type;
    }
} 

void priv::set_almost_parallel_type(FaceTypeMap              &face_type_map,
                                    const CutMesh            &mesh,
                                    const Project3f &projection)
{
    for (const FI &fi : mesh.faces()) {
        auto &type = face_type_map[fi];
        if (type != FaceType::inside) continue;
        //*
        assert(is_toward_projection(fi, mesh, projection));
        /*/
        if (!is_toward_projection(fi, mesh, projection)) {
            type = FaceType::outside;
        }else
        // */
        if (is_almost_parallel(fi, mesh, projection))
            // change type
            type = FaceType::inside_parallel;
    }
}

bool priv::is_almost_parallel(FI fi, const CutMesh &mesh, const Project3f &projection, float threshold)
{
    HI hi = mesh.halfedge(fi);
    std::array<VI, 3> vis = {
        mesh.source(hi),
        mesh.target(hi),
        mesh.target(mesh.next(hi))
    };
    std::array<Vec3f, 3> vertices;
    for (size_t i = 0; i < 3; i++) { 
        const P3 &p3 = mesh.point(vis[i]);
        vertices[i]  = Vec3f(p3.x(), p3.y(), p3.z());
    }

    Vec3f projected = projection.project(vertices[0]);
    Vec3f project_dir = projected - vertices[0];
    project_dir.normalize();
    Vec3f v1 = vertices[1] - vertices[0];
    v1.normalize();
    Vec3f v2 = vertices[2] - vertices[0];
    v2.normalize();
    // face normal
    Vec3f v_perp = v1.cross(v2);
    v_perp.normalize();
    
    float cos_alpha = project_dir.dot(v_perp);
    return cos_alpha <= threshold;
}

priv::ShapePoint2index::ShapePoint2index(const ExPolygons &shapes) {
    // prepare offsets
    m_offsets.reserve(shapes.size());
    uint32_t offset = 0;
    for (const auto &shape : shapes) {
        assert(!shape.contour.points.empty());
        std::vector<uint32_t> shape_offsets(shape.holes.size() + 1);

        shape_offsets[0] = offset;
        offset += shape.contour.points.size();

        for (uint32_t i = 0; i < shape.holes.size(); i++) { 
            shape_offsets[i + 1] = offset;
            offset += shape.holes[i].points.size();
        }
        m_offsets.push_back(std::move(shape_offsets));
    }
    m_count = offset;
}


uint32_t priv::ShapePoint2index::calc_index(const ShapePointId &id) const {
    assert(id.expolygons_index < m_offsets.size());
    const std::vector<uint32_t> &shape_offset =
        m_offsets[id.expolygons_index];
    assert(id.polygon_index < shape_offset.size());
    uint32_t res = shape_offset[id.polygon_index] + id.point_index;
    assert(res < m_count);
    return res;
}

priv::ShapePointId priv::ShapePoint2index::calc_id(uint32_t index) const {
    assert(index < m_count);
    ShapePointId result;
    // find shape index
    result.expolygons_index = 0;
    for (size_t i = 1; i < m_offsets.size(); i++) { 
        if (m_offsets[i][0] > index) break;
        result.expolygons_index = i;
    }

    // find contour index
    const std::vector<uint32_t> &shape_offset =
        m_offsets[result.expolygons_index];
    result.polygon_index = 0;
    for (size_t i = 1; i < shape_offset.size(); i++) {
        if (shape_offset[i] > index) break;
        result.polygon_index = i;     
    }

    // calculate point index
    uint32_t polygon_offset = shape_offset[result.polygon_index];
    assert(index >= polygon_offset);
    result.point_index = index - polygon_offset;
    return result;
}

uint32_t priv::ShapePoint2index::get_count() const { return m_count; }


void priv::flood_fill_inner(const CutMesh &mesh,
                            FaceTypeMap   &face_type_map)
{
    std::vector<FI> process;
    // guess count of connected not constrained triangles
    size_t guess_size = 128;
    process.reserve(guess_size);

    // check if neighbor(one of three in triangle) has type inside
    auto has_inside_neighbor = [&mesh, &face_type_map](FI fi) {
        HI hi     = mesh.halfedge(fi);
        HI hi_end = hi;
        auto exist_next = [&hi, &hi_end, &mesh]() -> bool {
            hi = mesh.next(hi);
            return hi != hi_end;
        };
        // loop over 3 half edges of face
        do {
            HI hi_opposite = mesh.opposite(hi);
            // open edge doesn't have opposit half edge
            if (!hi_opposite.is_valid()) continue;
            FI fi_opposite = mesh.face(hi_opposite);
            if (!fi_opposite.is_valid()) continue;
            if (face_type_map[fi_opposite] == FaceType::inside) return true;
        } while (exist_next());
        return false;
    };

    for (FI fi : mesh.faces()) {
        FaceType type = face_type_map[fi];
        if (type != FaceType::not_constrained &&
            type != FaceType::inside_parallel) continue;
        if (!has_inside_neighbor(fi)) continue;
        assert(process.empty());
        process.push_back(fi); 
        //store(mesh, face_type_map, DEBUG_OUTPUT_DIR + "progress.off");

        while (!process.empty()) {
            FI process_fi = process.back();
            process.pop_back();
            // Do not fill twice
            FaceType& process_type = face_type_map[process_fi];
            if (process_type == FaceType::inside) continue;
            process_type = FaceType::inside;

            // check neighbor triangle
            HI hi = mesh.halfedge(process_fi);
            HI hi_end = hi;
            auto exist_next = [&hi, &hi_end, &mesh]() -> bool {
                hi = mesh.next(hi);
                return hi != hi_end;
            };
            do {
                HI hi_opposite = mesh.opposite(hi);
                // open edge doesn't have opposit half edge
                if (!hi_opposite.is_valid()) continue;                
                FI fi_opposite = mesh.face(hi_opposite);
                if (!fi_opposite.is_valid()) continue;                
                FaceType type_opposite = face_type_map[fi_opposite];
                if (type_opposite == FaceType::not_constrained || 
                    type_opposite == FaceType::inside_parallel)
                    process.push_back(fi_opposite);
            } while (exist_next());
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

    if (intersection_ptr->shape_point_index == std::numeric_limits<uint32_t>::max()) {
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

bool priv::has_minimal_contour_points(const std::vector<HI> &outlines,
                                      const VertexShapeMap  &vert_shape_map,
                                      const CutMesh         &mesh,
                                      size_t                 min_count)
{
    // unique vector of indicies point from source contour(ExPolygon)
    std::vector<uint32_t> point_indicies;
    point_indicies.reserve(min_count);
    for (HI hi : outlines) { 
        VI vi = mesh.source(hi);
        const auto& shape = vert_shape_map[vi];
        if (shape == nullptr) continue;        
        uint32_t pi = shape->shape_point_index;
        if (pi == std::numeric_limits<uint32_t>::max()) continue;
        // is already stored in vector? 
        if (std::find(point_indicies.begin(), point_indicies.end(), pi) 
            != point_indicies.end())
            continue;
        if (point_indicies.size() == min_count) return true;
        point_indicies.push_back(pi);
    }
    // prevent weird small pieces
    return false;
}

void priv::collect_surface_data(std::queue<FI>  &process,
                                std::vector<FI> &faces,
                                std::vector<HI> &outlines,
                                FaceTypeMap     &face_type_map,
                                const CutMesh   &mesh)
{
    assert(!process.empty());
    assert(faces.empty());
    assert(outlines.empty());
    while (!process.empty()) {
        FI fi = process.front();
        process.pop();

        FaceType &fi_type = face_type_map[fi];
        // Do not process twice
        if (fi_type == FaceType::inside_processed) continue;
        assert(fi_type == FaceType::inside);
        // flag face as processed
        fi_type = FaceType::inside_processed;
        faces.push_back(fi);

        // check neighbor triangle
        HI hi     = mesh.halfedge(fi);
        HI hi_end = hi;
        do {
            HI hi_opposite = mesh.opposite(hi);
            // open edge doesn't have opposit half edge
            if (!hi_opposite.is_valid()) { 
                outlines.push_back(hi);
                hi = mesh.next(hi);
                continue; 
            }            
            FI fi_opposite = mesh.face(hi_opposite);
            if (!fi_opposite.is_valid()) {
                outlines.push_back(hi);
                hi = mesh.next(hi);
                continue;
            }
            FaceType side = face_type_map[fi_opposite];
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
    // IMPROVE: find better way to initialize
    // initialize reduction map
    for (VI reduction_from : mesh.vertices()) 
        reduction_map[reduction_from] = reduction_from;

    // check if vertex was made by edge_2 which is diagonal of quad
    auto is_reducible_vertex = [&vert_shape_map](VI reduction_from) -> bool {
        const IntersectingElement *ie = vert_shape_map[reduction_from];
        if (ie == nullptr) return false;
        IntersectingElement::Type type = ie->get_type();
        return type == IntersectingElement::Type::edge_2;
    };

    /// <summary>
    /// Append reduction or change existing one.
    /// </summary>
    /// <param name="hi">HalEdge between outside and inside face.
    /// Target vertex will be reduced, source vertex left</param>
    /// [[maybe_unused]] &face_type_map, &is_reducible_vertex are need only in debug
    auto add_reduction = [&] //&reduction_map, &mesh, &face_type_map, &is_reducible_vertex
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
    // clear v2v
    // more than one cut can share vertex and each cut need its own conversion
    for (FI fi : faces) {
        HI hi = mesh.halfedge(fi);
        for (VI vi : {mesh.source(hi), mesh.target(hi), mesh.target(mesh.next(hi))})
            v2v[vi] = std::numeric_limits<SurfaceCut::Index>::max();                
    }

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

std::vector<priv::ProjectionDistances> priv::create_distances(
    const CutAOIs        &cuts,
    const CutMesh        &mesh,
    uint32_t              shapes_points,
    const CutMesh        &shapes_mesh,
    float                 projection_ratio,
    const VertexShapeMap &vert_shape_map)
{
    // calculate distance from projection ration [in mm]
    auto calc_distance = [&mesh, &shapes_mesh, projection_ratio](uint32_t pi, VI vi) -> float {
        const P3& p = mesh.point(vi);
        
        // It is known because shapes_mesh is created inside of private space
        VI vi_start(2 * pi);
        VI vi_end(2 * pi + 1);

        // Get range for intersection
        const P3 &start = shapes_mesh.point(vi_start);
        const P3 &end   = shapes_mesh.point(vi_end);

        size_t max_i = 0;
        float  max_val = 0.f;
        for (size_t i = 0; i < 3; i++) { 
            float val = start[i] - end[i];
            // abs value
            if (val < 0.f) val *= -1.f;
            if (max_val < val) { 
                max_val = val;
                max_i   = i;
            }
        }
        return (p[max_i] - start[max_i]) - projection_ratio * (end[max_i] - start[max_i]);
    };

    std::vector<ProjectionDistances> distances(shapes_points);    
    for (const CutAOI &cut : cuts) {
        // for each half edge of outline
        for (const HI& hi : cut.second) {
            VI vi = mesh.source(hi);
            const IntersectingElement * ie = vert_shape_map[vi];
            if (ie == nullptr) continue;
            assert(ie->shape_point_index != std::numeric_limits<uint32_t>::max());
            assert(ie->attr != (unsigned char) IntersectingElement::Type::undefined);
            uint32_t pi = ie->shape_point_index;
            std::vector<ProjectionDistance> &pds = distances[pi];

            ProjectionDistance pd;
            pd.aoi_index = &cut - &cuts.front();
            pd.hi_index  = &hi - &cut.second.front();
            pd.distance = calc_distance(pi, vi);
            pds.push_back(std::move(pd));
        }        
    }
    return distances;
}

priv::ProjectionDistances priv::choose_best_distance(
    const std::vector<ProjectionDistances> &distances,
    const ExPolygons                       &shapes,
    const BoundingBox                      &shapes_bb,
    const ShapePoint2index                 &s2i)
{
    // euler square size of vector stored in just created Point
    auto calc_size_sq = [](const Point &p) -> float {
        return (float)p.x() * p.x() + (float)p.y() * p.y();
    };
    struct ClosePoint{
        // index of closest point from another shape
        uint32_t index = std::numeric_limits<uint32_t>::max();
        // squere distance to index
        float dist_sq = std::numeric_limits<float>::max();
    };
    // search in all shapes points to found closest point to given point
    auto get_closest_point_index = [&shapes, &distances, &calc_size_sq]
    (const Point &p)->uint32_t{
        ClosePoint cp;
        uint32_t id{0};
        auto get_closest = [&distances, &p, &id, &cp, &calc_size_sq]
        (const Points &pts) {
            for (const Point &p_ : pts) {
                if (distances[id].empty()) { 
                    ++id;
                    continue;
                }
                float d = calc_size_sq(p - p_);
                if (cp.dist_sq > d) {
                    cp.dist_sq = d;
                    cp.index  = id;
                }
                ++id;
            }
        };
        for (const ExPolygon &shape : shapes) {
            get_closest(shape.contour.points);
            for (const Polygon &hole : shape.holes) 
                get_closest(hole.points);
        }
        return cp.index;
    };
    // Search for closest projection to wanted distance
    auto get_closest_projection = []
    (const ProjectionDistances& distance, float wanted_distance) -> const ProjectionDistance *{
        // minimal distance
        float min_d = std::numeric_limits<float>::max();
        const ProjectionDistance *min_pd = nullptr;
        for (const ProjectionDistance &pd : distance) { 
            float d = std::fabs(pd.distance - wanted_distance);
            // There should be limit for maximal distance
            if (min_d > d) { 
                min_d = d;
                min_pd = &pd;
            }
        }
        return min_pd;
    };

    // return neighbor projection distance when exists
    auto get_next = [&get_closest_projection]
    (const ProjectionDistance &from_pd, const ProjectionDistances &from,
    const ProjectionDistances &to) -> const ProjectionDistance* {
        // exist some projection?
        if (to.empty()) return {};

        // find next same aoi (closest one)
        const ProjectionDistance* to_pd = nullptr;
        for (const ProjectionDistance &t : to) { 
            if (t.aoi_index != from_pd.aoi_index) continue;            
            if (to_pd != nullptr) {
                // when exist more than one use closest to previous
                float distance_prev = std::fabs(to_pd->distance - from_pd.distance);
                float distance = std::fabs(t.distance - from_pd.distance);
                if (distance < distance_prev) 
                    to_pd = &t;
            } else {
                to_pd = &t;
            }            
        }

        if (to_pd != nullptr) {
            // detect crossing aois
            const ProjectionDistance* cross_pd = nullptr;
            for (const ProjectionDistance &t : to) { 
                if (t.distance > to_pd->distance) continue;
                for (const ProjectionDistance &f : from) { 
                    if (f.aoi_index != t.aoi_index) continue;
                    if (f.distance < from_pd.distance) continue;
                    if (cross_pd!=nullptr) {
                        // multiple crossing
                        if (cross_pd->distance > f.distance) 
                            cross_pd = &f;
                    } else {
                        cross_pd = &f;
                    }
                }
            }
            // TODO: Detect opposit crossing - should be fixed
            if (cross_pd!=nullptr) return cross_pd;
        } else {
            // Try find another closest AOI 
            return get_closest_projection(to, from_pd.distance);
        }
        return to_pd;
    };

    ProjectionDistances result(distances.size());
    // fill result around known index inside one polygon
    auto fill_polygon_distances = [&distances, &shapes, &result, &get_next]
    (const ProjectionDistance &pd, uint32_t index, const ShapePointId& id){        
        const ExPolygon &shape = shapes[id.expolygons_index];
        const Points  & points = (id.polygon_index == 0) ?
                                                shape.contour.points :
                                                shape.holes[id.polygon_index - 1].points;
        // border of indexes for Polygon
        uint32_t first_index = index - id.point_index;
        uint32_t last_index  = first_index + points.size();        

        uint32_t act_index = index;
        const ProjectionDistance* act_pd = &pd;
        const ProjectionDistances* act_distances = &distances[act_index];

        // Copy starting pd to result
        result[act_index] = pd;

        auto exist_next = [&distances, &act_index, &act_pd, &act_distances, get_next, &result]
        (uint32_t nxt_index) {
            const ProjectionDistances* nxt_distances = &distances[nxt_index];
            const ProjectionDistance *nxt_pd = get_next(*act_pd, *act_distances, *nxt_distances);
            // exist next projection distance ?
            if (nxt_pd == nullptr) return false;

            // check no rewrite result
            assert(result[nxt_index].aoi_index == std::numeric_limits<uint32_t>::max());
            // copy founded projection to result
            result[nxt_index] = *nxt_pd; // copy

            // next
            act_index = nxt_index;
            act_pd        = &result[nxt_index];
            act_distances = nxt_distances;
            return true;
        };

        // last index in circle
        uint32_t finish_index = (index == first_index) ? (last_index - 1) :
                                                         (index - 1);
        // Positive iteration inside polygon
        do {
            uint32_t nxt_index = act_index + 1;
            // close loop of indexes inside of contour
            if (nxt_index == last_index) nxt_index = first_index;
            // check that exist next 
            if (!exist_next(nxt_index)) break;            
        } while (act_index != finish_index);

        // when all results for polygon are set no neccessary to iterate negative
        if (act_index == finish_index) return;

        act_index     = index;
        act_pd        = &pd;
        act_distances = &distances[act_index];
        // Negative iteration inside polygon
        do {
            uint32_t nxt_index = (act_index == first_index) ? 
                (last_index-1) : (act_index - 1);
            // When iterate negative it must be split to parts
            // and can't iterate in circle
            assert(nxt_index != index);
            // check that exist next 
            if (!exist_next(nxt_index)) break;
        } while (true);
    };

    std::vector<bool> finished_shapes(shapes.size(), {false});
    // choose correct cut by source point
    auto fill_shape_distances = [&distances, &s2i, &shapes, &result, &fill_polygon_distances, &calc_size_sq, &finished_shapes]
    (uint32_t known_point, const ProjectionDistance &pd) {
        const ProjectionDistance *start_pd = &pd;
        uint32_t start_index = known_point;
        uint32_t expolygons_index = s2i.calc_id(known_point).expolygons_index;
        uint32_t first_shape_index = s2i.calc_index({expolygons_index, 0, 0});
        const ExPolygon &shape = shapes[expolygons_index];
        do {
            fill_polygon_distances(*start_pd, start_index, s2i.calc_id(start_index));
            // seaching only inside shape, return index of closed finished point
            auto find_close_finished_point = [&first_shape_index, &shape, &result, &calc_size_sq]
            (const Point &p) -> ClosePoint {
                uint32_t index = first_shape_index;                
                ClosePoint cp;
                auto check_finished_points = [&cp, &result, &index, &p, &calc_size_sq]
                (const Points& pts) { 
                    for (const Point &p_ : pts) {
                        // finished point with some distances
                        if (result[index].aoi_index == std::numeric_limits<uint32_t>::max()) {
                            ++index;
                            continue;
                        }
                        float distance = calc_size_sq(p_ - p);
                        if (cp.dist_sq > distance) { 
                            cp.dist_sq = distance;
                            cp.index   = index;
                        }
                        ++index;
                    }
                };
                check_finished_points(shape.contour.points);
                for (const Polygon &h : shape.holes)
                    check_finished_points(h.points);
                return cp;
            };
                        
            // find next closest pair of points
            // (finished + unfinished) in ExPolygon
            start_index = std::numeric_limits<uint32_t>::max(); // unfinished_index
            uint32_t finished_index = std::numeric_limits<uint32_t>::max();
            float dist_sq = std::numeric_limits<float>::max();

            // first index in shape
            uint32_t index = first_shape_index;
            auto check_unfinished_points = [&index, &result, &distances, &find_close_finished_point, &dist_sq, &start_index, &finished_index]
            (const Points& pts) { 
                for (const Point &p : pts) {
                    // try find unfinished
                    if (result[index].aoi_index !=
                        std::numeric_limits<uint32_t>::max() ||
                        distances[index].empty()) {
                        ++index;
                        continue;
                    }
                    ClosePoint cp = find_close_finished_point(p);
                    if (dist_sq > cp.dist_sq) { 
                        dist_sq = cp.dist_sq;
                        start_index = index;
                        finished_index = cp.index;
                    }
                    ++index;
                }
            };
            // for each unfinished points
            check_unfinished_points(shape.contour.points);
            for (const Polygon &h : shape.holes)
                check_unfinished_points(h.points);
        } while (start_index != std::numeric_limits<uint32_t>::max());
        finished_shapes[expolygons_index] = true;
    };    

    // find close points between finished and unfinished ExPolygons
    auto find_close_point = [&shapes, &finished_shapes, &s2i, &calc_size_sq, &result]
    (const Point &p) -> ClosePoint {
        // result
        ClosePoint cp;
        // for all finished points
        for (uint32_t shape_index = 0; shape_index < shapes.size(); ++shape_index) {
            if (!finished_shapes[shape_index]) continue;
            const ExPolygon &shape = shapes[shape_index];
            uint32_t index = s2i.calc_index({shape_index, 0, 0});
            auto find_close_point_in_points = [&p, &cp, &index, &calc_size_sq, &result]
            (const Points &pts){
                for (const Point &p_ : pts) {
                    // Exist result (is finished) ?
                    if (result[index].aoi_index ==
                        std::numeric_limits<uint32_t>::max()) {
                        ++index;
                        continue;
                    }
                    float distance_sq = calc_size_sq(p - p_);
                    if (cp.dist_sq > distance_sq) { 
                        cp.dist_sq = distance_sq;
                        cp.index = index;
                    }
                    ++index;
                }
            };
            find_close_point_in_points(shape.contour.points);
            // shape could be inside of another shape's hole
            for (const Polygon& h:shape.holes)
                find_close_point_in_points(h.points);
        }
        return cp;
    };

    
    // wanted distance from ideal projection
    float wanted_distance = 0.f;
    // NOTE: it should be dependent on allign of text
    Point center = shapes_bb.center();

    // Select first point of shapes
    uint32_t unfinished_index = get_closest_point_index(center);
    // selection of closest_id should proove that pd has value
    do {
        const ProjectionDistance* pd = get_closest_projection(distances[unfinished_index], wanted_distance);
        assert(pd != nullptr);
        fill_shape_distances(unfinished_index, *pd);

        // The most close points between finished and unfinished shapes
        unfinished_index = std::numeric_limits<uint32_t>::max();
        ClosePoint best_cp; // must be finished
        
        // for each unfinished points 
        for (uint32_t shape_index = 0; shape_index < shapes.size(); ++shape_index) {
            if (finished_shapes[shape_index]) continue;
            const ExPolygon &shape  = shapes[shape_index];
            uint32_t index = s2i.calc_index({shape_index, 0, 0});
            auto find_close_point_in_points =
                [&unfinished_index, &best_cp,
                &index, &find_close_point, &distances]
            (const Points &pts) {
                for (const Point &p : pts) {
                    if (distances[index].empty()){ 
                        ++index;
                        continue;
                    }
                    ClosePoint cp = find_close_point(p);
                    if (cp.index != std::numeric_limits<uint32_t>::max() &&
                        best_cp.dist_sq > cp.dist_sq) {
                        best_cp = cp; // copy
                        unfinished_index = index;
                    }
                    ++index;
                }
            };
            find_close_point_in_points(shape.contour.points);
            // shape could be inside of another shape's hole
            for (const Polygon &h : shape.holes)
                find_close_point_in_points(h.points);
        }
        // detect finish (best doesn't have value)
        if (best_cp.index == std::numeric_limits<uint32_t>::max()) break;

        const ProjectionDistance &closest_pd = result[best_cp.index];
        // check that best_cp is finished and has result
        assert(closest_pd.aoi_index != std::numeric_limits<uint32_t>::max());
        wanted_distance = closest_pd.distance;
    } while (unfinished_index != std::numeric_limits<uint32_t>::max());
    return result;
}

// store projection center as circle
void priv::store(const Vec3f       &vertex,
                 const Vec3f       &normal,
                 const std::string &file,
                 float              size)
{
    int flatten = 20;
    size_t min_i = 0;
    for (size_t i = 1; i < 3; i++)
        if (normal[min_i] > normal[i]) 
            min_i = i;
    Vec3f up_ = Vec3f::Zero();
    up_[min_i] = 1.f;
    Vec3f side = normal.cross(up_).normalized() * size;
    Vec3f up = side.cross(normal).normalized() * size;

    indexed_triangle_set its;
    its.vertices.reserve(flatten + 1);
    its.indices.reserve(flatten);

    its.vertices.push_back(vertex);
    its.vertices.push_back(vertex + up);
    for (size_t i = 1; i < flatten; i++) {
        float angle = i * 2 * M_PI / flatten;
        Vec3f v     = vertex + sin(angle) * side + cos(angle) * up;
        its.vertices.push_back(v);
        its.indices.emplace_back(0, i, i + 1);
    }
    its.indices.emplace_back(0, flatten, 1);
    its_write_obj(its, file.c_str());
}

bool priv::merge_cut(CutAOI &cut1, const CutAOI &cut2, const CutMesh &mesh)
{
    // create cgal model and merge it together

    return false;
}

void priv::merge_cuts(CutAOIs                   &cuts,
                      const CutMesh             &mesh,
                      const std::vector<size_t> &use_cut_indices)
{
    auto create_bb = [&mesh](const CutAOI &cut) -> BoundingBoxf3 {
        Vec3f min(std::numeric_limits<float>::min(),
                  std::numeric_limits<float>::min(),
                  std::numeric_limits<float>::min());
        Vec3f max(std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max());
        for (const FI &fi : cut.first) {
            HI hi = mesh.halfedge(fi);
            for (VI vi : {mesh.source(hi), mesh.target(hi),
                          mesh.target(mesh.next(hi))}) {
                const P3 &p = mesh.point(vi);
                for (size_t i = 0; i < 3; i++) {
                    if (min[i] > p[i]) min[i] = p[i];
                    if (max[i] < p[i]) max[i] = p[i];
                }
            }
        }
        return BoundingBoxf3(min.cast<double>(), max.cast<double>());
    };

    // create bounding boxes for cuts
    std::vector<BoundingBoxf3> bbs;
    bbs.reserve(cuts.size());    
    for (const CutAOI &cut : cuts)
        bbs.push_back(create_bb(cut));
    // extend used bb by intersecting bb
    // NOTE: after merge 2 cuts could appear new intersection on surface of merged in

    std::vector<std::pair<size_t, size_t>> merge_order;
    std::vector<bool> del_cuts(cuts.size(), {true});
    for (size_t cut_index : use_cut_indices) del_cuts[cut_index] = false;

    // find intersection of cuts by Bounding boxes intersection
    for (size_t cut_index : use_cut_indices) {
        // check if cut is merged into another one
        if (del_cuts[cut_index]) continue;
        BoundingBoxf3& result_bb = bbs[cut_index];
        CutAOI &cut = cuts[cut_index];
        // all merged cuts into cut_index
        std::vector<bool> merged(cuts.size(), {false});

        // merged in last iteration
        std::vector<bool> new_merged;
        bool exist_new_extension;
        bool is_first = true;
        do {
            exist_new_extension = false;
            new_merged = std::vector<bool>(cuts.size(), {false});
            for (const BoundingBoxf3 &bb : bbs) {
                size_t bb_index = &bb - &bbs.front();
                // do not merge itself
                if (cut_index == bb_index) continue;
                if (!is_first && merged[cut_index]) continue;
                if (!bb.intersects(result_bb)) {
                    if (is_first) continue;
                    bool has_new_intersection = false;
                    for (size_t i = 0; i < cuts.size(); i++) {
                        if (!new_merged[i]) continue;
                        if (!bbs[i].intersects(bb)) continue;
                        has_new_intersection = true;
                    }
                    if (!has_new_intersection) continue;
                }
                if(!merge_cut(cut, cuts[bb_index], mesh)) continue;
                result_bb = create_bb(cut);
                merged[bb_index]     = true;
                del_cuts[bb_index]   = true;
                new_merged[bb_index] = true;
                // extend result_bb
                exist_new_extension = true;
            }
            is_first = false;
        } while (exist_new_extension);        
    }

    // remove flagged cuts
    for (size_t i = del_cuts.size(); i > 0; --i) {
        size_t index = i - 1;
        if (del_cuts[index]) cuts.erase(cuts.begin() + index);
    }
}

void priv::filter_cuts(CutAOIs              &cuts,
                       const CutMesh        &mesh,
                       const ExPolygons     &shapes,
                       const ShapePoint2index &shape_point_2_index,
                       const Project        &projection,
                       const VertexShapeMap &vert_shape_map)
{
    auto get_point = [&shapes, &shape_point_2_index]
    (const IntersectingElement &intersection) -> Point {
        assert(intersection.shape_point_index != std::numeric_limits<uint32_t>::max());
        ShapePointId point_id = shape_point_2_index.calc_id(intersection.shape_point_index);
        const ExPolygon& shape = shapes[point_id.expolygons_index];
        const Polygon   &p     = (point_id.polygon_index == 0) ?
                                     shape.contour :
                                     shape.holes[point_id.polygon_index - 1];
        return p[point_id.point_index];
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
        
        assert(i->shape_point_index != std::numeric_limits<uint32_t>::max());
        assert(i->attr != (unsigned char)IntersectingElement::Type::undefined);

        // Use only straigh edge
        if (i->get_type() != IntersectingElement::Type::edge_1)
            return false;
    
        CutIndex &ci = indices[i->shape_point_index];

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
        Vec3f source_point = projection.create_front_back(p).first;
        const auto &prev = mesh.point(ci.vi);
        Vec3f prev_point(prev.x(), prev.y(), prev.z());
        float prev_sq_norm = (source_point - prev_point).squaredNorm();

        const auto &act = mesh.point(vi);
        Vec3f act_point(act.x(), act.y(), act.z());
        float act_sq_norm = (source_point - act_point).squaredNorm();
        
        if (act_sq_norm < prev_sq_norm) {
            del_cuts[cut_index] = true;
            return true;
        }

        // previous cut is behind actual one
        del_cuts[ci.cut_index] = true;
        ci.cut_index = cut_index;
        ci.vi = vi;
        return false;
    };

    // filter small pieces
    for (const CutAOI &cut : cuts) {
        if (!has_minimal_contour_points(cut.second, vert_shape_map, mesh)) { 
            size_t index = &cut - &cuts.front();
            del_cuts[index] = true;
        }
    }

    // filter top one cuts
    for (const CutAOI &cut : cuts) {
        size_t cut_index = &cut - &cuts.front();
        if (del_cuts[cut_index]) continue;
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

#ifdef DEBUG_OUTPUT_DIR
void priv::store(CutMesh &mesh, const FaceTypeMap &face_type_map, const std::string& file)
{
    auto face_colors = mesh.add_property_map<priv::FI, CGAL::Color>("f:color").first;    
    for (FI fi : mesh.faces()) { 
        auto &color = face_colors[fi];
        switch (face_type_map[fi]) {
        case FaceType::inside: color = CGAL::Color{100, 250, 100}; break; // light green
        case FaceType::inside_parallel: color = CGAL::Color{255, 0, 0}; break; // red
        case FaceType::inside_processed: color = CGAL::Color{170, 0, 0}; break; // dark red
        case FaceType::outside: color = CGAL::Color{100, 0, 100}; break; // purple
        case FaceType::not_constrained: color = CGAL::Color{127, 127, 127}; break; // gray
        default: color = CGAL::Color{0, 0, 255}; // blue
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

indexed_triangle_set priv::create_indexed_triangle_set(
    const std::vector<FI> &faces, const CutMesh &mesh)
{
    std::vector<VI> vertices;
    vertices.reserve(faces.size() * 2);

    indexed_triangle_set its;
    its.indices.reserve(faces.size());
    for (FI fi : faces) {
        HI hi     = mesh.halfedge(fi);
        HI hi_end = hi;

        int   ti = 0;
        Vec3i t;

        do {
            VI   vi  = mesh.source(hi);
            auto res = std::find(vertices.begin(), vertices.end(), vi);
            t[ti++]  = res - vertices.begin();
            if (res == vertices.end()) vertices.push_back(vi);
            hi = mesh.next(hi);
        } while (hi != hi_end);

        its.indices.push_back(t);
    }

    its.vertices.reserve(vertices.size());
    for (VI vi : vertices) {
        const auto &p = mesh.point(vi);
        its.vertices.emplace_back(p.x(), p.y(), p.z());
    }
    return its;
}

#include <filesystem>
namespace {
static void prepare_dir(const std::string &dir)
{
    namespace fs = std::filesystem;
    if (fs::exists(dir)) {
        for (auto &path : fs::directory_iterator(dir)) fs::remove_all(path);
    } else {
        fs::create_directories(dir);
    }
}
} // namespace
void priv::store(const CutAOIs &aois, const CutMesh &mesh, const std::string &dir) {
    auto create_outline_its =
        [&mesh](const std::vector<HI> &outlines) -> indexed_triangle_set {
        static const float line_width = 0.1f;
        indexed_triangle_set its;
        its.indices.reserve(2*outlines.size());
        its.vertices.reserve(outlines.size()*4);
        for (HI hi : outlines) { 
            //FI fi = mesh.face(hi);
            VI vi_a = mesh.source(hi);
            VI vi_b = mesh.target(hi);
            VI vi_c = mesh.target(mesh.next(hi));
            P3 p3_a = mesh.point(vi_a);
            P3 p3_b = mesh.point(vi_b);
            P3 p3_c = mesh.point(vi_c);

            Vec3f a(p3_a.x(), p3_a.y(), p3_a.z());
            Vec3f b(p3_b.x(), p3_b.y(), p3_b.z());
            Vec3f c(p3_c.x(), p3_c.y(), p3_c.z());

            Vec3f v1 = b - a; // from a to b
            v1.normalize();
            Vec3f v2 = c - a; // from a to c
            v2.normalize();
            Vec3f norm = v1.cross(v2);
            norm.normalize();
            Vec3f perp_to_edge = norm.cross(v1);
            perp_to_edge.normalize();
            Vec3f dir = -perp_to_edge * line_width;

            size_t ai = its.vertices.size();
            its.vertices.push_back(a);
            size_t bi = its.vertices.size();
            its.vertices.push_back(b);
            size_t ai2 = its.vertices.size();
            its.vertices.push_back(a + dir);
            size_t bi2 = its.vertices.size();
            its.vertices.push_back(b + dir);

            its.indices.push_back(Vec3i(ai, ai2, bi));
            its.indices.push_back(Vec3i(ai2, bi2, bi));
        }
        return its;    
    };

    prepare_dir(dir);
    for (const auto &aoi : aois) {
        size_t      index = &aoi - &aois.front();
        std::string file  = dir + "aoi" + std::to_string(index) + ".obj";
        indexed_triangle_set its = create_indexed_triangle_set(aoi.first, mesh);
        its_write_obj(its, file.c_str());

        // exist some outline?
        if (aoi.second.empty()) continue;
        std::string file_outline = dir + "outline" + std::to_string(index) + ".obj";
        indexed_triangle_set outline = create_outline_its(aoi.second);
        its_write_obj(outline, file_outline.c_str());
    }
}

void priv::store(const ProjectionDistances &pds,
                 const CutAOIs             &aois,
                 const CutMesh             &mesh,
                 const std::string         &file,
                 float                      width)
{
    // create rectangle for each half edge from projection distances
    indexed_triangle_set its;
    its.vertices.reserve(4 * pds.size());
    its.indices.reserve(2 * pds.size());
    for (const ProjectionDistance &pd : pds) {
        if (pd.aoi_index == std::numeric_limits<uint32_t>::max()) continue;
        HI hi  = aois[pd.aoi_index].second[pd.hi_index];
        VI vi1 = mesh.source(hi);
        VI vi2 = mesh.target(hi);
        VI vi3 = mesh.target(mesh.next(hi));
        const P3 &p1  = mesh.point(vi1);
        const P3 &p2  = mesh.point(vi2);
        const P3 &p3  = mesh.point(vi3);
        Vec3f v1(p1.x(), p1.y(), p1.z());
        Vec3f v2(p2.x(), p2.y(), p2.z());
        Vec3f v3(p3.x(), p3.y(), p3.z());

        Vec3f v12 = v2 - v1;
        v12.normalize();
        Vec3f v13 = v3 - v1;
        v13.normalize();
        Vec3f n = v12.cross(v13);
        n.normalize();
        Vec3f side = n.cross(v12);
        side.normalize();
        side *= -width;
        
        uint32_t i = its.vertices.size();
        its.vertices.push_back(v1);
        its.vertices.push_back(v1+side);
        its.vertices.push_back(v2);
        its.vertices.push_back(v2+side);

        its.indices.emplace_back(i, i + 1, i + 2);
        its.indices.emplace_back(i + 2, i + 1, i + 3);
    }
    its_write_obj(its, file.c_str());
}

void priv::store(const SurfaceCuts &cut, const std::string &dir) {
    auto create_contour_its =
        [](const indexed_triangle_set& its, const std::vector<unsigned int> &contour)
        -> indexed_triangle_set {
        static const float line_width = 0.1f;

        auto get_triangle_tip = [&its](unsigned int vi1,
                                       unsigned int vi2) -> const Vec3f& { 
            for (const auto &t : its.indices) { 
                unsigned int tvi = std::numeric_limits<unsigned int>::max();
                for (const auto &vi : t) { 
                    if (vi == vi1) continue;
                    if (vi == vi2) continue;
                    if (tvi == std::numeric_limits<unsigned int>::max()) {
                        tvi = vi;
                    } else {
                        tvi = std::numeric_limits<unsigned int>::max();
                        break;
                    }
                }
                if (tvi != std::numeric_limits<unsigned int>::max())
                    return its.vertices[tvi];
            }
            // triangle with indices vi1 and vi2 doesnt exist
            assert(false);
        };

        indexed_triangle_set result;
        result.vertices.reserve((contour.size() + 1) * 4);
        result.indices.reserve((contour.size() + 1) * 2);
        unsigned int prev_vi = contour.back();
        for (unsigned int vi : contour) {
            const Vec3f &a = its.vertices[vi];
            const Vec3f &b = its.vertices[prev_vi];
            const Vec3f &c = get_triangle_tip(vi, prev_vi);

            Vec3f v1 = b - a; // from a to b
            v1.normalize();
            Vec3f v2 = c - a; // from a to c
            v2.normalize();
            // triangle normal
            Vec3f norm = v1.cross(v2);
            norm.normalize();
            // perpendiculat to edge lay on triangle 
            Vec3f perp_to_edge = norm.cross(v1);
            perp_to_edge.normalize();

            Vec3f dir = -perp_to_edge * line_width;

            size_t ai = result.vertices.size();
            result.vertices.push_back(a);
            size_t bi = result.vertices.size();
            result.vertices.push_back(b);
            size_t ai2 = result.vertices.size();
            result.vertices.push_back(a + dir);
            size_t bi2 = result.vertices.size();
            result.vertices.push_back(b + dir);

            result.indices.push_back(Vec3i(ai, bi, ai2));
            result.indices.push_back(Vec3i(ai2, bi, bi2));
            prev_vi = vi;
        }
        return result;
    };
    prepare_dir(dir);
    for (const auto &c : cut) {
        size_t index = &c - &cut.front();
        std::string file  = dir + "cut" + std::to_string(index) + ".obj";
        its_write_obj(c, file.c_str());
        for (const auto& contour : c.contours) {
            size_t c_index = &contour - &c.contours.front();
            std::string c_file = dir + "cut" + std::to_string(index) + 
                "contour" + std::to_string(c_index) + ".obj";
            indexed_triangle_set c_its = create_contour_its(c, contour);
            its_write_obj(c_its, c_file.c_str());
        }
    }
}
#endif // DEBUG_OUTPUT_DIR
