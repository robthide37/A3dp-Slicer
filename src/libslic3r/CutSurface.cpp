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

SurfaceCut Slic3r::merge(SurfaceCuts &&cuts)
{ 
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
using CutMeshes = std::vector<CutMesh>;
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

// stored in model made by shape
using EdgeShapeMap = CutMesh::Property_map<EI, IntersectingElement>;
using FaceShapeMap = CutMesh::Property_map<FI, IntersectingElement>;

// stored in surface source - pointer to EdgeShapeMap | FaceShapeMap
using VertexShapeMap = CutMesh::Property_map<VI, const IntersectingElement *>;

// stored in model made by shape
const std::string edge_shape_map_name = "e:IntersectingElement";
const std::string face_shape_map_name = "f:IntersectingElement";

// stored in surface source
const std::string vert_shape_map_name = "v:IntersectingElement";

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
/// Convert triangle mesh model to CGAL Surface_mesh
/// </summary>
/// <param name="its">Input mesh model</param>
/// <param name="edges_count">It depends on count of opened edge</param>
/// <returns>CGAL mesh - half edge mesh</returns>
CutMesh to_cgal(const indexed_triangle_set &its, size_t edges_count = 0);

/// <summary>
/// Covert 2d shape (e.g. Glyph) to CGAL model
/// NOTE: internaly create 
/// edge_shape_map .. Property map to store conversion from edge to contour
/// face_shape_map .. Property map to store conversion from face to contour
/// </summary>
/// <param name="shapes">2d shapes to project</param>
/// <param name="projection">Define transformation 2d point into 3d</param>
/// <returns>CGAL model of extruded shape</returns>
CutMesh to_cgal(const ExPolygons &shapes,
                const Project    &projection);

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

/// <summary>
/// Track source of intersection 
/// Help for anotate inner and outer faces
/// </summary>
struct Visitor {
    const CutMesh &object;
    const CutMesh &shape;

    // Properties of the shape mesh:
    EdgeShapeMap edge_shape_map;
    FaceShapeMap face_shape_map;

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
const std::string face_type_map_name = "f:side";

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
/// Note: also use from mesh (have to be created)
/// face_type_map .. Type of shape inside / outside
/// vert_shape_map .. Source of outline vertex
/// </summary>
/// <param name="reduction_map">Reduction map from vertex to vertex, 
/// when key == value than no reduction</param>
/// <param name="faces">Faces of one </param> 
/// <param name="mesh">Input object</param>
void create_reduce_map(ReductionMap &reduction_map, const CutMesh &meshes);

// Patch made by Cut area of interest from model
// connected faces(triangles) and outlines(halfEdges) for one surface cut
using CutAOI = std::pair<std::vector<FI>, std::vector<HI>>;
// vector of Cutted Area of interest cutted from one CGAL model
using CutAOIs = std::vector<CutAOI>;
// vector of CutAOIs for each model
using VCutAOIs = std::vector<CutAOIs>;

struct ModelCutId
{
    // index of model
    uint32_t model_index;
    // index of cut inside model
    uint32_t cut_index;
};

/// <summary>
/// Keep conversion from VCutAOIs to Index and vice versa
/// Model_index .. contour(or hole) poin from ExPolygons
/// Index      .. continous number
/// </summary>
class ModelCut2index
{
    std::vector<uint32_t> m_offsets;
    // for check range of index
    uint32_t m_count;
public:
    ModelCut2index(const VCutAOIs &cuts);
    uint32_t   calc_index(const ModelCutId &id) const;
    ModelCutId calc_id(uint32_t index) const;
    uint32_t   get_count() const;
};

/// <summary>
/// Cut surface from model
/// </summary>
/// <param name="cgal_model">Input model converted to CGAL</param>
/// <param name="shapes">2d contours</param>
/// <param name="cgal_shape">[const]Model made by shapes
/// NOTE: Can't be definde as const because of corefine function input definition,
/// but it is.</param> <param name="s2i">Convert index to shape point from
/// ExPolygons</param> <returns>Patches from model surface</returns>
CutAOIs cut_from_model(CutMesh                &cgal_model,
                       const ExPolygons       &shapes,
                       /*const*/ CutMesh      &cgal_shape,
                       float                   projection_ratio,
                       const ShapePoint2index &s2i);

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
    // index of source model
    uint32_t model_index = std::numeric_limits<uint32_t>::max();

    // index of CutAOI
    uint32_t aoi_index = std::numeric_limits<uint32_t>::max();

    // index of half edge in AOI
    uint32_t hi_index = std::numeric_limits<uint32_t>::max();

    // signed distance to projection 
    float distance = std::numeric_limits<float>::max();    
};
// addresed by ShapePoint2index
using ProjectionDistances =  std::vector<ProjectionDistance>;

// each point in shapes has its ProjectionDistances
using VDistances = std::vector<ProjectionDistances>;

/// <summary>
/// Calculate distances from CutAOI contour points to ProjectionOrigin
/// NOTE:
/// each model has to have "vert_shape_map" .. Know source of new vertices
/// </summary>
/// <param name="cuts">models AOIs</param>
/// <param name="count_shapes_points">Count of contour points in shapes</param>
/// <param name="models">Vertices position</param>
/// <param name="shapes_mesh">Mesh created by shapes</param>
/// <param name="projection_ratio">Define best distnace</param>
/// <param name="vert_shape_map">Know source of new vertices</param>
/// <returns>Projection distances of cutted shape points</returns>
VDistances calc_distances(const VCutAOIs       &cuts,
                          size_t                count_shapes_points,
                          const CutMeshes      &models,
                          const CutMesh        &shapes_mesh,
                          float                 projection_ratio);
    /// <summary>
/// Select distances in similar depth between expolygons
/// </summary>
/// <param name="distances">All distances</param>
/// <param name="shapes">Vector of letters</param>
/// <param name="shape_point_2_index">Convert index to addresss inside of shape</param>
/// <returns>Best projection distances</returns>
ProjectionDistances choose_best_distance(
    const std::vector<ProjectionDistances> &distances,
    const ExPolygons                      &shapes,
    const ShapePoint2index                &shape_point_2_index);

using ConvertMap = CutMesh::Property_map<VI, SurfaceCut::Index>;
/// <summary>
/// Create surface cuts from mesh model
/// </summary>
/// <param name="cutAOIs">Areas of interests from model surface</param>
/// <param name="mesh">Model - can't be const because of create temporary property map</param>
/// <param name="reduction_map">Reduction of vertices</param>
/// <returns>Created surface cuts</returns>
SurfaceCuts create_surface_cuts(const CutAOIs      &cutAOIs,
                                CutMesh            &mesh,
                                const ReductionMap &reduction_map);

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

/// <summary>
/// Self Intersection of surface cuts are made by 
/// damaged models OR multi volumes emboss
/// </summary>
/// <param name="cuts">Surface cuts to merge
/// NOTE: Merge process move data from cuts to result</param>
/// <param name="use_cut">Mask for wanted cuts [Same size as cuts]</param>
/// <returns>Merged all surface cuts into one</returns>
SurfaceCut merge_intersections(SurfaceCuts &cuts, const CutAOIs& cutAOIs, const std::vector<bool>& use_cut);

/// <summary>
/// Differenciate other models
/// </summary>
/// <param name="cuts">Patches from meshes</param>
/// <param name="use_cut">Define wanted cuts (addressing by m2i)</param>
/// <param name="m2i">Convert model_index and cut_index into one index</param>
/// <param name="cut_models">Source points for Cutted AOIs</param>
/// <param name="models">Original models without cut modifications used for differenciate</param>
void diff_models(VCutAOIs                &cuts,
                 const std::vector<bool> &use_cut,
                 const ModelCut2index    &m2i,
                 const CutMeshes         &cut_models,
                 /*const*/ CutMeshes     &models);

// keep CGAL Mesh for next processing
struct SurfaceCutWithMesh : public SurfaceCut{
    CutMesh cgalMesh = CutMesh();
};


/// <summary>
/// Merge 2 Cuts when has intersection
/// </summary>
/// <param name="cut1">In/Out cut to merge into</param>
/// <param name="cut2">Cut to merge from</param>
/// <returns>Has intersection</returns>
bool merge_intersection(SurfaceCut &cut1, const SurfaceCut &cut2);

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
void store(const ProjectionDistances &pds, const VCutAOIs &aois, const CutMeshes &meshes, const std::string &file, float width = 0.2f/* [in mm] */);
void store(const SurfaceCuts &cut, const std::string &dir);

void store(const std::vector<indexed_triangle_set> &models, const std::string &obj_filename);
void store(const std::vector<CutMesh>&models, const std::string &off_filename);
void store(const Emboss::IProjection &projection, const Point &point_to_project, float projection_ratio, const std::string &obj_filename);
#endif // DEBUG_OUTPUT_DIR
} // namespace privat

SurfaceCut Slic3r::cut_surface(const ExPolygons &shapes,
                               const std::vector<indexed_triangle_set> &models,
                               const Emboss::IProjection &projection,
                               float projection_ratio)
{
    assert(!models.empty());
    assert(!shapes.empty());
    if (models.empty() || shapes.empty() ) return {};

#ifdef DEBUG_OUTPUT_DIR
    priv::store(models, DEBUG_OUTPUT_DIR + "models_input.obj");
#endif // DEBUG_OUTPUT_DIR

    // for filter out triangles out of bounding box
    BoundingBox shapes_bb = get_extents(shapes);
    Point projection_center = shapes_bb.center();
#ifdef DEBUG_OUTPUT_DIR
    priv::store(projection, projection_center, projection_ratio,
                DEBUG_OUTPUT_DIR + "projection_center.obj");
#endif // DEBUG_OUTPUT_DIR

    // for filttrate opposite triangles and a little more
    const float max_angle = 89.f;

    priv::CutMeshes cgal_models;
    cgal_models.reserve(models.size());
    for (const indexed_triangle_set &its : models) {
        std::vector<bool> skip_indicies(its.indices.size(), {false});
        priv::set_skip_for_out_of_aoi(skip_indicies, its, projection, shapes_bb);
        // cut out opposit triangles
        // priv::set_skip_for_outward_projection(skip_indicies, model, projection);
        // cut out more than only opposit triangles 
        priv::set_skip_by_angle(skip_indicies, its, projection, max_angle);
        cgal_models.push_back(priv::to_cgal(its, skip_indicies));
    }

#ifdef DEBUG_OUTPUT_DIR
    priv::store(cgal_models, DEBUG_OUTPUT_DIR + "model");// model[0-N].off
#endif // DEBUG_OUTPUT_DIR

    priv::CutMesh cgal_shape = priv::to_cgal(shapes, projection);

#ifdef DEBUG_OUTPUT_DIR
    CGAL::IO::write_OFF(DEBUG_OUTPUT_DIR + "shape.off", cgal_shape); // only debug
#endif // DEBUG_OUTPUT_DIR
    
    // create tool for convert index to shape Point adress and vice versa
    priv::ShapePoint2index s2i(shapes);

    // create copy of source models for cut intersections
    priv::CutMeshes cgal_models_copy = cgal_models;
    priv::VCutAOIs model_cuts;
    // cut shape from each cgal model
    for (priv::CutMesh &cgal_model : cgal_models) { 
        priv::CutAOIs cutAOIs = priv::cut_from_model(
            cgal_model, shapes, cgal_shape, projection_ratio, s2i);
#ifdef DEBUG_OUTPUT_DIR
        size_t index = &cgal_model - &cgal_models.front();
        priv::store(cutAOIs, cgal_model, DEBUG_OUTPUT_DIR + "aois" + std::to_string(index) + "/"); // only debug
#endif // DEBUG_OUTPUT_DIR
        model_cuts.push_back(std::move(cutAOIs));
    }

    // calc distance to projection for all outline points of cutAOI(shape)
    // it is used for distiguish the top one
    uint32_t shapes_points = s2i.get_count();
    // for each point collect all projection distances
    priv::VDistances distances = priv::calc_distances(
        model_cuts, shapes_points, cgal_models, cgal_shape, projection_ratio);

    // for each point select best projection
    priv::ProjectionDistances best_projection =
        priv::choose_best_distance(distances, shapes, s2i);

#ifdef DEBUG_OUTPUT_DIR
    store(best_projection, model_cuts, cgal_models, DEBUG_OUTPUT_DIR + "best_projection.obj"); // only debug
#endif // DEBUG_OUTPUT_DIR

    priv::ModelCut2index m2i(model_cuts);
    // Create mask for wanted AOIs
    std::vector<bool> is_best_cut(m2i.get_count(), {false});  
    for (const priv::ProjectionDistance& d : best_projection)
        if (d.model_index != std::numeric_limits<uint32_t>::max())
            is_best_cut[m2i.calc_index({d.model_index, d.aoi_index})] = true;    
    priv::diff_models(model_cuts, is_best_cut, m2i, cgal_models, cgal_models_copy);
    
    // IMPROVE: create reduce map on demand - may be model do not need it (when it is not used for result)
    // Reduction prepare
    std::string vertex_reduction_map_name = "v:reduction";
    for (size_t model_index = 0; model_index < models.size(); model_index++) {
        priv::CutMesh &cgal_model = cgal_models[model_index];
        priv::ReductionMap vertex_reduction_map = cgal_model.add_property_map<priv::VI, priv::VI>(vertex_reduction_map_name).first;
        priv::create_reduce_map(vertex_reduction_map, cgal_model);
#ifdef DEBUG_OUTPUT_DIR
        priv::store(cgal_model, vertex_reduction_map, DEBUG_OUTPUT_DIR + "reduction" + std::to_string(model_index) + ".off");
#endif // DEBUG_OUTPUT_DIR
    }

    // NOTE: it will be fine to calc AOIs range,
    // not only outline but all vertices in direction of emboss - faster check
    // on intersection


    // NOTE: It is not neccessary to convert all AOIs to SurfaceCut
    // but exist case where AOI intersect best AOI without edge intersection
    // So filtering is made from SurfaceCuts in merge function

    //SurfaceCuts surface_cuts = create_surface_cuts(cutAOIs, cgal_model);
//#ifdef DEBUG_OUTPUT_DIR
//    store(surface_cuts, DEBUG_OUTPUT_DIR + "cuts/"); // only debug
//#endif                                               // DEBUG_OUTPUT_DIR
//
//    return surface_cuts;


    // TODO: cut of inside part of cuts
    
    //priv::merge_aois(cutAOIs, is_best_cut, cgal_model);

    // Self Intersection of surface cuts are
    // made by damaged models AND multi volumes
    //SurfaceCut result = priv::merge_intersections(surface_cuts, cutAOIs, is_best_cut);

    SurfaceCut result;
    //for (SurfaceCuts &cuts : model_cuts)
    //    for (SurfaceCut &cut : cuts) append(result, std::move(cut));
//#ifdef DEBUG_OUTPUT_DIR
//    its_write_obj(result, (DEBUG_OUTPUT_DIR + "resultCut.obj").c_str()); // only debug
//#endif // DEBUG_OUTPUT_DIR
    return result;
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

priv::CutMesh priv::to_cgal(const indexed_triangle_set &its, size_t edges_count)
{
    CutMesh result;
    if (its.empty()) return result;

    const std::vector<stl_vertex>                  &vertices = its.vertices;
    const std::vector<stl_triangle_vertex_indices> &indices  = its.indices;

    size_t vertices_count = vertices.size();
    size_t faces_count    = indices.size();
    if (edges_count == 0) edges_count = (faces_count * 3) / 2;
    result.reserve(vertices_count, edges_count, faces_count);

    for (const stl_vertex &v : vertices)
        result.add_vertex(CutMesh::Point{v.x(), v.y(), v.z()});

    for (const stl_triangle_vertex_indices &f : indices)
        result.add_face(static_cast<VI>(f[0]), static_cast<VI>(f[1]),
                        static_cast<VI>(f[2]));

    return result;
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
                            const Project     &projection)
{
    CutMesh result;
    if (shapes.empty()) return result;
        
    EdgeShapeMap& edge_shape_map = result.add_property_map<EI, IntersectingElement>(edge_shape_map_name).first;
    FaceShapeMap& face_shape_map = result.add_property_map<FI, IntersectingElement>(face_shape_map_name).first;

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


priv::CutAOIs priv::cut_from_model(CutMesh                &cgal_model,
                                   const ExPolygons       &shapes,
                                   CutMesh                &cgal_shape,
                                   float                   projection_ratio,
                                   const ShapePoint2index &s2i)
{
    // pointer to edge or face shape_map
    VertexShapeMap vert_shape_map = cgal_model.add_property_map<VI, const IntersectingElement*>(vert_shape_map_name).first;
    
    // detect anomalities in visitor.
    bool is_valid = true;
    // NOTE: map are created when convert shapes to cgal model
    const EdgeShapeMap& edge_shape_map = cgal_shape.property_map<EI, IntersectingElement>(edge_shape_map_name).first;
    const FaceShapeMap& face_shape_map = cgal_shape.property_map<FI, IntersectingElement>(face_shape_map_name).first;
    Visitor visitor{cgal_model, cgal_shape, edge_shape_map, face_shape_map, vert_shape_map, &is_valid};

    // bool map for affected edge
    EcmType ecm = get(DynamicEdgeProperty(), cgal_model);

    const auto &p = CGAL::parameters::visitor(visitor)
                        .edge_is_constrained_map(ecm)
                        .throw_on_self_intersection(false);
    const auto& q = CGAL::parameters::do_not_modify(true);
    CGAL::Polygon_mesh_processing::corefine(cgal_model, cgal_shape, p, q);

    if (!is_valid) return {};

    FaceTypeMap face_type_map = cgal_model.add_property_map<FI, FaceType>(face_type_map_name).first;

    // Select inside and outside face in model
    set_face_type(face_type_map, cgal_model, vert_shape_map, ecm, cgal_shape, s2i);
#ifdef DEBUG_OUTPUT_DIR
    store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "constrained.off"); // only debug
#endif // DEBUG_OUTPUT_DIR

    // It is neccesary when almost parallel face are contained in projection
    // set_almost_parallel_type(face_type_map, cgal_model, projection);
//#ifdef DEBUG_OUTPUT_DIR
//    store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "constrainedWithAlmostParallel.off"); // only debug
//#endif // DEBUG_OUTPUT_DIR
    
    // flood fill the other faces inside the region.
    flood_fill_inner(cgal_model, face_type_map);

#ifdef DEBUG_OUTPUT_DIR
    store(cgal_model, face_type_map, DEBUG_OUTPUT_DIR + "filled.off"); // only debug
#endif // DEBUG_OUTPUT_DIR
        
    // IMPROVE: AOIs area could be created during flood fill
    return create_cut_area_of_interests(cgal_model, shapes, face_type_map);
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

priv::ModelCut2index::ModelCut2index(const VCutAOIs &cuts)
{
    // prepare offsets
    m_offsets.reserve(cuts.size());
    uint32_t offset = 0;
    for (const CutAOIs &model_cuts: cuts) {
        m_offsets.push_back(offset);
        offset += model_cuts.size();
    }
    m_count = offset;
}

uint32_t priv::ModelCut2index::calc_index(const ModelCutId &id) const
{
    assert(id.model_index < m_offsets.size());
    uint32_t offset = m_offsets[id.model_index];
    uint32_t res = offset + id.cut_index;
    assert(((id.model_index+1) < m_offsets.size() && res < m_offsets[id.model_index+1]) ||
           ((id.model_index+1) == m_offsets.size() && res < m_count));
    return res;
}

priv::ModelCutId priv::ModelCut2index::calc_id(uint32_t index) const
{
    assert(index < m_count);
    ModelCutId result;
    // find shape index
    result.model_index = 0;
    for (size_t model_index = 1; model_index < m_offsets.size(); ++model_index) {
        if (m_offsets[model_index] > index) break;
        result.model_index = model_index;
    }
    result.cut_index = index - m_offsets[result.model_index];
    return result;
}

uint32_t priv::ModelCut2index::get_count() const { return m_count; }


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

void priv::create_reduce_map(ReductionMap &reduction_map, const CutMesh &mesh)
{
    const FaceTypeMap &face_type_map = mesh.property_map<FI, FaceType>(face_type_map_name).first;
    const VertexShapeMap &vert_shape_map = mesh.property_map<VI, const IntersectingElement*>(vert_shape_map_name).first;

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

priv::VDistances priv::calc_distances(const std::vector<CutAOIs> &cuts,
                                      size_t count_shapes_points,
                                      const std::vector<CutMesh> &models,
                                      const CutMesh              &shapes_mesh,
                                      float                 projection_ratio)
{
    // calculate distance from projection ration [in mm]
    auto calc_distance = [&models, &shapes_mesh, projection_ratio](uint32_t pi, VI vi, uint32_t model_index) -> float {
        const P3& p = models[model_index].point(vi);
        
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

    priv::VDistances result(count_shapes_points);
    for (uint32_t model_index = 0; model_index < cuts.size(); model_index++) {
        const CutAOIs& model_cuts = cuts[model_index];
        // map is created during intersection by corefine visitor
        const VertexShapeMap &vert_shape_map = 
            models[model_index].property_map<VI, const IntersectingElement*>(vert_shape_map_name).first;
        for (const CutAOI &cut : model_cuts) {
            // for each half edge of outline
            for (const HI& hi : cut.second) {
                VI vi = models[model_index].source(hi);
                const IntersectingElement * ie = vert_shape_map[vi];
                if (ie == nullptr) continue;
                assert(ie->shape_point_index != std::numeric_limits<uint32_t>::max());
                assert(ie->attr != (unsigned char) IntersectingElement::Type::undefined);
                uint32_t pi = ie->shape_point_index;
                assert(pi <= count_shapes_points);
                std::vector<ProjectionDistance> &pds = result[pi];
                uint32_t aoi_index = &cut - &model_cuts.front();
                uint32_t hi_index = &hi - &cut.second.front();
                float distance = calc_distance(pi, vi, model_index);
                pds.push_back({model_index, aoi_index, hi_index, distance});
            }        
        }
    }
    return result;
}

// functions for choose_best_distance
namespace priv {

// euler square size of vector stored in Point
float calc_size_sq(const Point &p);

// structure to store index and distance together
struct ClosePoint
{
    // index of closest point from another shape
    uint32_t index = std::numeric_limits<uint32_t>::max();
    // squere distance to index
    float dist_sq = std::numeric_limits<float>::max();
};

// search in all shapes points to found closest point to given point
uint32_t get_closest_point_index(const Point &p, const ExPolygons &shapes, const VDistances &distances);

// Search for closest projection to wanted distance
const ProjectionDistance *get_closest_projection(const ProjectionDistances &distance, float wanted_distance);

// return neighbor projection distance when exists
const ProjectionDistance *get_next(const ProjectionDistance  &from_pd,
                                   const ProjectionDistances &from,
                                   const ProjectionDistances &to);

// fill result around known index inside one polygon
void fill_polygon_distances(const ProjectionDistance &pd, uint32_t index, const ShapePointId &id, ProjectionDistances & result, const ExPolygons &shapes, const VDistances &distances);

// choose correct cut by source point
void fill_shape_distances(uint32_t known_point, const ProjectionDistance &pd, ProjectionDistances &result, std::vector<bool>& finished_shapes, const ShapePoint2index &s2i, const ExPolygons &shapes, const VDistances &distances);

// find close points between finished and unfinished ExPolygons
ClosePoint find_close_point(const Point &p, ProjectionDistances &result, std::vector<bool>& finished_shapes,const ShapePoint2index &s2i, const ExPolygons &shapes);

}

float priv::calc_size_sq(const Point &p){
    return (float) p.x() * p.x() + (float) p.y() * p.y();
}

uint32_t priv::get_closest_point_index(const Point      &p,
                                       const ExPolygons &shapes,
                                       const VDistances &distances)
{
    ClosePoint cp;
    uint32_t id{0};
    auto get_closest = [&distances, &p, &id, &cp]
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
}

const priv::ProjectionDistance *priv::get_closest_projection(
    const ProjectionDistances &distance, float wanted_distance)
{
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
}

const priv::ProjectionDistance *priv::get_next(
    const ProjectionDistance  &from_pd,
    const ProjectionDistances &from,
    const ProjectionDistances &to)
{
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
}

void priv::fill_polygon_distances(const ProjectionDistance &pd,
                                  uint32_t                  index,
                                  const ShapePointId       &id,
                                  ProjectionDistances      &result,
                                  const ExPolygons         &shapes,
                                  const VDistances         &distances)
{
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

    auto exist_next = [&distances, &act_index, &act_pd, &act_distances, &result]
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
}

void priv::fill_shape_distances(uint32_t                  known_point,
                                const ProjectionDistance &pd,
                                ProjectionDistances      &result,
                                std::vector<bool>        &finished_shapes,
                                const ShapePoint2index   &s2i,
                                const ExPolygons         &shapes,
                                const VDistances         &distances)
{
    const ProjectionDistance *start_pd = &pd;
    uint32_t start_index = known_point;
    uint32_t expolygons_index = s2i.calc_id(known_point).expolygons_index;
    uint32_t first_shape_index = s2i.calc_index({expolygons_index, 0, 0});
    const ExPolygon &shape = shapes[expolygons_index];
    do {
        fill_polygon_distances(*start_pd, start_index, s2i.calc_id(start_index),result, shapes, distances);
        // seaching only inside shape, return index of closed finished point
        auto find_close_finished_point = [&first_shape_index, &shape, &result]
        (const Point &p) -> ClosePoint {
            uint32_t index = first_shape_index;                
            ClosePoint cp;
            auto check_finished_points = [&cp, &result, &index, &p]
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
}

priv::ClosePoint priv::find_close_point(const Point         &p,
                                        ProjectionDistances &result,
                                        std::vector<bool>   &finished_shapes,
                                        const ShapePoint2index &s2i,
                                        const ExPolygons       &shapes)
{
    // result
    ClosePoint cp;
    // for all finished points
    for (uint32_t shape_index = 0; shape_index < shapes.size(); ++shape_index) {
        if (!finished_shapes[shape_index]) continue;
        const ExPolygon &shape = shapes[shape_index];
        uint32_t index = s2i.calc_index({shape_index, 0, 0});
        auto find_close_point_in_points = [&p, &cp, &index, &result]
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
}

// IMPROVE: create better structure to find closest points e.g. Tree
priv::ProjectionDistances priv::choose_best_distance(
    const std::vector<ProjectionDistances> &distances,
    const ExPolygons                       &shapes,
    const ShapePoint2index                 &s2i)
{
    // collect one closest projection for each outline point
    ProjectionDistances result(distances.size());

    // store info about finished shapes
    std::vector<bool> finished_shapes(shapes.size(), {false});
    
    // wanted distance from ideal projection
    // Distances are relative to projection distance
    // so first wanted distance is the closest one (ZERO)
    float wanted_distance = 0.f;

    // NOTE: Shapes are centered to respect allign of text
    Point center(0, 0);
    // Select point from shapes(text contour) which is closest to center (all in 2d)
    uint32_t unfinished_index = get_closest_point_index(center, shapes, distances);

    do {
        const ProjectionDistance* pd = get_closest_projection(distances[unfinished_index], wanted_distance);
        // selection of closest_id should proove that pd has value 
        // (functions: get_closest_point_index and find_close_point_in_points)
        assert(pd != nullptr);
        fill_shape_distances(unfinished_index, *pd, result, finished_shapes, s2i, shapes, distances);

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
                &index, &result, &finished_shapes, &distances, &s2i, &shapes]
            (const Points &pts) {
                for (const Point &p : pts) {
                    if (distances[index].empty()){ 
                        ++index;
                        continue;
                    }
                    ClosePoint cp = find_close_point(p, result, finished_shapes, s2i, shapes);
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

// functions to help 'merge_intersection'
namespace priv {
const VI default_vi(std::numeric_limits<uint32_t>::max());

struct Source
{
    HI  hi;
    int sdim;
};
using Sources = std::vector<Source>;
const std::string vertex_source_map_name = "v:SourceIntersecting";
using VertexSourceMap = CutMesh::Property_map<VI, Source>;

/// <summary>
/// Corefine visitor
/// Store intersection source for vertices of constrained edge of tm1
/// Must be used with corefine flag no modification of tm2
/// </summary>
struct IntersectionSources
{
    const CutMesh *patch; // patch
    const CutMesh *model; // const model

    VertexSourceMap vmap;

    // keep sources from call intersection_point_detected
    // until call new_vertex_added
    Sources m_sources;

    // count intersections
    void intersection_point_detected(std::size_t i_id, int sdim, HI h_f, HI h_e, 
        const CutMesh &tm_f, const CutMesh &tm_e, bool is_target_coplanar, bool is_source_coplanar){
        Source source;
        if (&tm_e == model) {
            source = {h_e, sdim};            
            // check other CGAL model that is patch
            assert(&tm_f == patch);
            if (is_target_coplanar) vmap[tm_f.source(h_f)] = source;
            if (is_source_coplanar) vmap[tm_f.target(h_f)] = source;

            // clear source to be able check that this intersection source is not used any more
            if (is_source_coplanar || is_target_coplanar) source = {};
        } else {
            source = {h_f, sdim};
            assert(&tm_f == model && &tm_e == patch);
            assert(!is_target_coplanar);
            assert(!is_source_coplanar);
            //if (is_target_coplanar) vmap[tm_e.source(h_e)] = source;
            //if (is_source_coplanar) vmap[tm_e.target(h_e)] = source;
            //if (sdim == 0)
            //    vmap[tm_e.target(h_e)] = source;
        }

        // add source of intersection
        m_sources.push_back(source);

        // By documentation i_id is consecutive.
        // check id goes in a row, without skips
        assert(m_sources.size() == i_id);
    }

    /// <summary>
    /// Store VI to intersections by i_id
    /// </summary>
    /// <param name="i_id">Order number of intersection point</param>
    /// <param name="v">New added vertex</param>
    /// <param name="tm">Affected mesh</param>
    void new_vertex_added(std::size_t i_id, VI v, const CutMesh &tm) {
        // check that it is first insertation into item of vmap
        assert(!vmap[v].hi.is_valid());
        // check valid addresing into sources
        assert(i_id < m_sources.size());
        // check that source has value
        assert(m_sources[i_id].hi.is_valid());
        vmap[v] = m_sources[i_id];
    }

    // Not used visitor functions
    void before_subface_creations(FI /* f_old */, CutMesh & /* mesh */) {}
    void after_subface_created(FI /* f_new */, CutMesh & /* mesh */) {}
    void after_subface_creations(CutMesh &) {}
    void before_subface_created(CutMesh &) {}
    void before_edge_split(HI /* h */, CutMesh & /* tm */) {}
    void edge_split(HI /* hnew */, CutMesh & /* tm */) {}
    void after_edge_split() {}
    void add_retriangulation_edge(HI /* h */, CutMesh & /* tm */) {}
};

/// <summary>
/// Create map1 and map2
/// </summary>
/// <param name="map">Convert tm1.face to type</param>
/// <param name="tm1">Corefined mesh</param>
/// <param name="tm2">Source of intersection</param>
/// <param name="ecm1">Identify constrainde edge</param>
/// <param name="sources">Convert tm1.face to type</param>
void create_face_types(FaceTypeMap& map, const CutMesh &tm1,const CutMesh &tm2, const EcmType &ecm, const VertexSourceMap& sources);

/// <summary>
/// Create map1 and map2
/// </summary>
/// <param name="tm1">First mesh</param>
/// <param name="map1">First map</param>
/// <param name="tm2">Second mesh</param>
/// <param name="map2">Second map</param>
/// <param name="vi1_to_vi2">Cvt vertex 1 to vertex 2</param>
/// <param name="ecm1">Identify constrainde edge</param>
//void create_face_types(const CutMesh &tm1, FaceTypeMap& map1, const CutMesh &tm2, const EcmType &ecm1);

//// Convert vertex index to indexed_triangle_set::index (index of vertext)
//using V2I_map = CutMesh::Property_map<VI, SurfaceCut::Index>;
//
//// Create its from merged CGAL models(without outline)
//indexed_triangle_set create_merged_its(
//    const CutMesh &tm1, const FaceTypeMap &map1, V2I_map &vi1_to_res,
//    const CutMesh &tm2, const FaceTypeMap &map2, V2I_map &vi2_to_res, const V2V_map &vi1_to_vi2);
}

//priv::V2V_map priv::create_map_from_vi1_to_vi2(
//    CutMesh                     &tm1,
//    const Store_VI_pairs::Pairs &intersections,
//    const std::string           &map_name)
//{
//    // Create conversion map from tm1.vertex to tm2.vertex on constrained edge
//    V2V_map vi2vi_map = tm1.add_property_map<VI, VI>(map_name).first;
//    // initialize to default value
//    for (VI vi : tm1.vertices()) vi2vi_map[vi] = default_vi;
//    // fill data by intersections
//    for (const std::pair<VI, VI> &intersection : intersections) {
//        // weird intersection?? what does it mean ??
//        if (intersection.first == default_vi ||
//            intersection.second == default_vi)
//            continue;
//
//        // is intersection point used only once?
//        assert(vi2vi_map[intersection.first] == default_vi);
//        vi2vi_map[intersection.first] = intersection.second;
//    }
//    return vi2vi_map;
//}

//void priv::create_face_types(const CutMesh &tm1,
//                             FaceTypeMap   &map1,
//                             const CutMesh &tm2,
//                             const EcmType &ecm1)
//{
//    // initialize maps
//    for (FI fi : tm1.faces()) map1[fi] = FaceType::not_constrained;
//    for (EI ei1 : tm1.edges()) { 
//        if (!get(ecm1, ei1)) continue;
//
//        // get faces from tm1 (f1a + f1b)
//        HI hi1 = tm1.halfedge(ei1);
//        assert(hi1.is_valid());
//        FI f1a = tm1.face(hi1);
//        assert(f1a.is_valid());
//        HI hi_op = tm1.opposite(hi1);
//        assert(hi_op.is_valid());
//        FI f1b = tm1.face(hi_op);
//        assert(f1b.is_valid());
//
//        // get faces from tm2 (f2a + f2b)
//        VI vi1_source = tm1.source(hi1);
//        assert(vi1_source.is_valid());
//        VI vi1_target = tm1.target(hi1);
//        assert(vi1_target.is_valid());
//        
//        // Need to know source triangle of intersection-----------------------------<<<<<<<<<<<<<<<
//
//        VI vi2_source = vi1_to_vi2[vi1_source];
//        assert(vi2_source != default_vi);
//        assert(vi2_source.is_valid());
//        VI vi2_target = vi1_to_vi2[vi1_target];
//        assert(vi2_target != default_vi);
//        assert(vi2_target.is_valid());
//
//        // target(halfedge(v)) == v
//        HI hi2_target = tm2.halfedge(vi2_target);
//        HI hi2; // half edge on tm2 in same position as hi
//        for (HI hi : tm2.halfedges_around_target(hi2_target)) { 
//            if (tm2.source(hi) != vi2_source) continue;            
//            hi2 = hi;
//            break;            
//        }
//
//        assert(hi2.is_valid());
//        FI f2a = tm2.face(hi2);
//        assert(f2a.is_valid());
//        HI hi2_op = tm2.opposite(hi2);
//        assert(hi2_op.is_valid());
//        FI f2b = tm2.face(hi2_op);
//        assert(f2a.is_valid());
//
//        // check orientation
//        const P3 &a = tm1.point(vi1_source);
//        const P3 &b = tm1.point(vi1_target);
//        // triangle tip from face f1a
//        VI vi1a_tip = tm1.target(tm1.next(hi1));
//        assert(vi1a_tip.is_valid());
//        const P3 &c = tm1.point(vi1a_tip);
//        // triangle tip from face f2a
//        VI vi2a_tip = tm2.target(tm2.next(hi2));
//        const P3 &p = tm2.point(vi2a_tip);
//
//        // check if f1a is behinde f2a
//        // inside mean it will be used
//        // outside will be discarded
//        if (CGAL::orientation(a, b, c, p) == CGAL::POSITIVE) {
//            map1[f1a] = FaceType::outside;
//            map1[f1b] = FaceType::inside;
//            map2[f2a] = FaceType::inside;
//            map2[f2b] = FaceType::outside;
//        } else {
//            map1[f1a] = FaceType::inside;
//            map1[f1b] = FaceType::outside;
//            map2[f2a] = FaceType::outside;
//            map2[f2b] = FaceType::inside;
//        }
//    }
//}

//indexed_triangle_set priv::create_merged_its(
//    const CutMesh &tm1, const FaceTypeMap &map1, V2I_map &vi1_to_res,
//    const CutMesh &tm2, const FaceTypeMap &map2, V2I_map &vi2_to_res, const V2V_map &vi1_to_vi2)
//{
//    // clear for sure
//    const SurfaceCut::Index invalid = std::numeric_limits<SurfaceCut::Index>::max();
//    for (VI vi:tm1.vertices()) vi1_to_res[vi] = invalid;    
//    for (VI vi:tm2.vertices()) vi2_to_res[vi] = invalid;
//
//    // count result triangles
//    size_t indices_size = 0;
//    for (FI fi : tm1.faces()) if (map1[fi] == FaceType::inside) ++indices_size;
//    for (FI fi : tm2.faces()) if (map2[fi] == FaceType::inside) ++indices_size;
//    
//    // approx maximal count
//    size_t vertices_count = 3 * indices_size;
//
//    indexed_triangle_set result;
//    result.indices.reserve(indices_size);
//    result.vertices.reserve(vertices_count);
//
//    // Start with tm2 because of convert map vi1_to_vi2 not vice versa
//    for (FI fi : tm2.faces()) {
//        if (map2[fi] != FaceType::inside) continue;
//        Vec3i t;
//        int   i = 0;
//        for (VI vi : tm2.vertices_around_face(tm2.halfedge(fi))) {
//            SurfaceCut::Index &si = vi2_to_res[vi];
//            if (si == invalid) {
//                si = result.vertices.size();
//                const P3 &p = tm2.point(vi);
//                Vec3f p_(p.x(), p.y(), p.z());
//                result.vertices.push_back(p_);
//            }
//            t[i++] = si;
//        }
//        result.indices.push_back(t);
//    }
//
//    // convert tm1 without same points from tm2 defined in map vi1_to_vi2
//    for (FI fi : tm1.faces()) {
//        if (map1[fi] != FaceType::inside) continue;
//        Vec3i t;
//        int i = 0;
//        for (VI vi : tm1.vertices_around_face(tm1.halfedge(fi))) {
//            SurfaceCut::Index &si = vi1_to_res[vi];
//            if (si == invalid) {
//                // check if it is not already created from tm2
//                VI vi2 = vi1_to_vi2[vi];
//                if (vi2 == default_vi) {
//                    si = result.vertices.size();
//                    const P3 &p = tm1.point(vi);
//                    Vec3f p_(p.x(), p.y(), p.z());
//                    result.vertices.push_back(p_);
//                } else {
//                    // constrained vertex copied from tm2
//                    si = vi2_to_res[vi2];
//                    assert(si != invalid);
//                }
//            }
//            t[i++] = si;
//        }
//        result.indices.push_back(t);
//    }
//    // fix approx of vertices count
//    result.vertices.shrink_to_fit();
//    return result;
//}
//
bool priv::merge_intersection(SurfaceCut &cut1, const SurfaceCut &cut2) {
    return false;
    //auto cut_to_cgal = [](const SurfaceCut &cut) {
    //    size_t count_edges = (cut.indices.size() * 3 + cut.contours.size()) / 2;
    //    return to_cgal(cut, count_edges);
    //};

    //CutMesh tm1 = cut_to_cgal(cut1);
    //CutMesh tm2 = cut_to_cgal(cut2);

    ////Store_VI_pairs::Pairs intersections;
    //Store_Source visitor = {&tm1, &tm2};

    //// bool map for affected edge
    //EcmType ecm1 = get(DynamicEdgeProperty(), tm1);
    //const auto &p = CGAL::parameters::visitor(visitor)
    //    .edge_is_constrained_map(ecm1)
    //    .throw_on_self_intersection(false);
    //const auto &q = CGAL::parameters::throw_on_self_intersection(false);

    //CGAL::Polygon_mesh_processing::corefine(tm1, tm2, p, q);
    //// when no intersection detected than no result surface cut
    //if (intersections.empty()) return false;

    //// Create conversion map from tm1.vertex to tm2.vertex on constrained edge
    //std::string vi1_to_vi2_name = "v:vi1_to_vi2";
    //V2V_map vi1_to_vi2 = create_map_from_vi1_to_vi2(tm1, intersections, vi1_to_vi2_name);

    //std::string face_type_map_name = "f:side";
    //FaceTypeMap face_type_map1 = tm1.add_property_map<FI, FaceType>(face_type_map_name).first;
    //FaceTypeMap face_type_map2 = tm2.add_property_map<FI, FaceType>(face_type_map_name).first;
    //create_face_types(tm1, face_type_map1, tm2, face_type_map2, vi1_to_vi2, ecm1);
    //
    //std::string dir = "C:/data/temp/out/";
    //store(tm1, face_type_map1, dir + "tm1_constrained.off");
    //store(tm2, face_type_map2, dir + "tm2_constrained.off");
    //flood_fill_inner(tm1, face_type_map1);
    //flood_fill_inner(tm2, face_type_map2);
    //store(tm1, face_type_map1, dir + "tm1_filled.off");
    //store(tm2, face_type_map2, dir + "tm2_filled.off");

    //std::string vi_to_res_name = "v:vertex_to_result";
    //V2I_map vi1_to_res = tm1.add_property_map<VI, SurfaceCut::Index>(vi_to_res_name).first;
    //V2I_map vi2_to_res = tm2.add_property_map<VI, SurfaceCut::Index>(vi_to_res_name).first;
    //indexed_triangle_set its = create_merged_its(
    //    tm1, face_type_map1, vi1_to_res,
    //    tm2, face_type_map2, vi2_to_res, vi1_to_vi2);
    //its_write_obj(its, (dir + "merged_result.obj").c_str());

    //// calculate contours:
    //


    //// set result into cut1
    //cut1.indices = std::move(its.indices);
    //cut1.vertices = std::move(its.vertices);
    //return true;
}

namespace priv {
struct ExtendAOI{
    // source for extend
    const CutAOI *source;
    // converted cut to CGAL mesh
    CutMesh mesh;
};

void diff_model(ExtendAOI &cut, /*const*/ CutMesh &another_model);



BoundingBoxf3 bounding_box(const CutAOI &cut, const CutMesh &mesh);
BoundingBoxf3 bounding_box(const ExtendAOI &ecut);

ExtendAOI create_extend_aoi(CutAOI &cut, const CutMesh &mesh);
} // namespace priv

void priv::create_face_types(FaceTypeMap           &map,
                             const CutMesh         &tm1,
                             const CutMesh         &tm2,
                             const EcmType         &ecm,
                             const VertexSourceMap &sources)
{
    auto get_intersection_source = [&tm2](const Source& s1, const Source& s2)->FI{        
        // when one of sources is face than return it
        FI fi1 = tm2.face(s1.hi);
        if (s1.sdim == 2) return fi1;
        FI fi2 = tm2.face(s2.hi);
        if (s2.sdim == 2) return fi2;
        // both vertices are made by same source triangle
        if (fi1 == fi2) return fi1;

        // when one from sources is edge second one decide side of triangle triangle
        HI hi1_opposit = tm2.opposite(s1.hi);
        FI fi1_opposit;
        if (hi1_opposit.is_valid())
            fi1_opposit = tm2.face(hi1_opposit);
        if (fi2 == fi1_opposit) return fi2;

        HI hi2_opposit = tm2.opposite(s2.hi);
        FI fi2_opposit;
        if (hi2_opposit.is_valid())
            fi2_opposit = tm2.face(hi2_opposit);
        if (fi1 == fi2_opposit) return fi1;
        if (fi1_opposit.is_valid() && fi1_opposit == fi2_opposit)
            return fi1_opposit;

        // when intersection is vertex need loop over neighbor
        for (FI fi_around_hi1 : tm2.faces_around_target(s1.hi)) {
            for (FI fi_around_hi2 : tm2.faces_around_target(s2.hi)) { 
                if (fi_around_hi1 == fi_around_hi2) 
                    return fi_around_hi1;
            }
        }
        return FI();
    };

    for (FI fi : tm1.faces()) map[fi] = FaceType::not_constrained;
    for (EI ei1 : tm1.edges()) {
        if (!get(ecm, ei1)) continue;

        // get faces from tm1 (f1a + f1b)
        HI hi1 = tm1.halfedge(ei1);
        assert(hi1.is_valid());
        FI f1a = tm1.face(hi1);
        assert(f1a.is_valid());
        HI hi_op = tm1.opposite(hi1);
        assert(hi_op.is_valid());
        FI f1b = tm1.face(hi_op);
        assert(f1b.is_valid());

        // get faces from tm2 (f2a + f2b)
        VI vi1_source = tm1.source(hi1);
        assert(vi1_source.is_valid());
        VI vi1_target = tm1.target(hi1);
        assert(vi1_target.is_valid());

        // Need to know source triangle of intersection-----------------------------<<<<<<<<<<<<<<<
        const Source &s_s = sources[vi1_source];
        const Source &s_t = sources[vi1_target];
        FI fi2 = get_intersection_source(s_s, s_t);
        HI hi2 = tm2.halfedge(fi2);
        std::array<const P3 *, 3> t;
        size_t ti =0;
        for (VI vi2 : tm2.vertices_around_face(hi2))
            t[ti++] = &tm2.point(vi2);

        // triangle tip from face f1a
        VI vi1a_tip = tm1.target(tm1.next(hi1));
        assert(vi1a_tip.is_valid());
        const P3 &p = tm1.point(vi1a_tip);

        // check if f1a is behinde f2a
        // inside mean it will be used
        // outside will be discarded
        if (CGAL::orientation(*t[0], *t[1], *t[2], p) == CGAL::POSITIVE) {
            map[f1a] = FaceType::outside;
            map[f1b] = FaceType::inside;
        } else {
            map[f1a] = FaceType::inside;
            map[f1b] = FaceType::outside;
        }
    }
}


void priv::diff_model(ExtendAOI &cut, CutMesh &another_model)
{
    CutMesh& tm1 = cut.mesh;
    CutMesh &tm2 = another_model;

    //Store_VI_pairs::Pairs intersections;
    std::string vertex_source_map_name = "v:source_intersections";
    VertexSourceMap vmap = tm1.add_property_map<VI, Source>(vertex_source_map_name).first;
    IntersectionSources visitor = {&tm1, &tm2, vmap};

    //// bool map for affected edge
    EcmType  ecm = get(DynamicEdgeProperty(), tm1);
    const auto &p = CGAL::parameters::visitor(visitor)
                        .edge_is_constrained_map(ecm)
                        .throw_on_self_intersection(false);
    const auto &q = CGAL::parameters::do_not_modify(true)
                        .throw_on_self_intersection(false);
    CGAL::Polygon_mesh_processing::corefine(tm1, tm2, p, q);


    //// when no intersection detected than no result surface cut
    //if (intersections.empty()) return;

    //// Create conversion map from tm1.vertex to tm2.vertex on constrained edge
    //std::string vi1_to_vi2_name = "v:vi1_to_vi2";
    //V2V_map vi1_to_vi2 = create_map_from_vi1_to_vi2(tm1, intersections, vi1_to_vi2_name);

    ////std::string face_type_map_name = "f:side";
    //FaceTypeMap face_type_map = tm1.add_property_map<FI, FaceType>(face_type_map_name).first;
    //create_face_types(tm1, face_type_map, tm2, ecm);

    //std::string dir = "C:/data/temp/out/";
    //store(tm1, face_type_map, dir + "constrained.off");
    //flood_fill_inner(tm1, face_type_map);
    //store(tm1, face_type_map, dir + "filled.off");

    //std::string vi_to_res_name = "v:vertex_to_result";
    //V2I_map vi1_to_res = tm1.add_property_map<VI, SurfaceCut::Index>(vi_to_res_name).first;

    // create result Surface cut
    
    //indexed_triangle_set its = create_merged_its(tm1, face_type_map1,
    //                                             vi1_to_res, tm2,
    //                                             face_type_map2, vi2_to_res,
    //                                             vi1_to_vi2);
    //its_write_obj(its, (dir + "merged_result.obj").c_str());

    // calculate contours:

    // set result into cut1
    //cut1.indices  = std::move(its.indices);
    //cut1.vertices = std::move(its.vertices);
}

BoundingBoxf3 priv::bounding_box(const CutAOI &cut, const CutMesh &mesh) {
    const P3& p_from_cut = mesh.point(mesh.target(mesh.halfedge(cut.first.front())));
    Vec3d min(p_from_cut.x(), p_from_cut.y(), p_from_cut.z());
    Vec3d max = min;
    for (FI fi : cut.first) { 
        for(VI vi: mesh.vertices_around_face(mesh.halfedge(fi))){
            const P3& p = mesh.point(vi);
            for (size_t i = 0; i < 3; ++i) { 
                if (min[i] > p[i]) min[i] = p[i];
                if (max[i] < p[i]) max[i] = p[i];
            }
        } 
    }
    return BoundingBoxf3(min, max);
}

BoundingBoxf3 priv::bounding_box(const ExtendAOI &ecut) {
    const CutMesh& mesh = ecut.mesh;
    const P3      &p_from_cut = *mesh.points().begin();
    Vec3d min(p_from_cut.x(), p_from_cut.y(), p_from_cut.z());
    Vec3d max = min;
    for (VI vi: mesh.vertices()) { 
        const P3& p = mesh.point(vi);
        for (size_t i = 0; i < 3; ++i) { 
            if (min[i] > p[i]) min[i] = p[i];
            if (max[i] < p[i]) max[i] = p[i];
        }
    }
    return BoundingBoxf3(min, max);
}

priv::ExtendAOI priv::create_extend_aoi(CutAOI &cut, const CutMesh &mesh)
{
    std::vector<bool> is_counted(mesh.vertices().size(), {false});
    uint32_t count_vertices = 0;
    for (FI fi : cut.first) { 
        for (VI vi : mesh.vertices_around_face(mesh.halfedge(fi))) { 
            if (!is_counted[vi.idx()]) { 
                is_counted[vi.idx()] = true;
                ++count_vertices;
            }
        }
    }

    uint32_t count_faces = cut.first.size();    
    // NOTE: It is more than neccessary, guess it better
    uint32_t count_edges = count_faces*3; 

    CutMesh cm;
    cm.reserve(count_vertices, count_edges, count_faces);

    // vertex conversion function
    constexpr uint32_t def_val = std::numeric_limits<uint32_t>::max();       
    std::vector<uint32_t> v_cvt(mesh.vertices().size(), {def_val});
    for (FI fi : cut.first) {
        std::array<VI, 3> t;
        int index = 0;
        for (VI vi : mesh.vertices_around_face(mesh.halfedge(fi))) {
            uint32_t &cvt = v_cvt[vi.idx()];
            if (cvt == def_val) {
                cvt = cm.vertices().size(); 
                cm.add_vertex(mesh.point(vi));
            }
            t[index++] = VI(cvt);
        }
        cm.add_face(t[0], t[1], t[2]);
    }
    return {&cut, cm};
}

void priv::diff_models(VCutAOIs                &cuts,
                       const std::vector<bool> &use_cut,
                       const ModelCut2index    &m2i,
                       const CutMeshes         &cut_models,
                       /*const*/ CutMeshes     &models)
{
    // create bounding boxes for cuts
    std::vector<BoundingBoxf3> bbs;    
    bbs.reserve(m2i.get_count());
    for (size_t model_index = 0; model_index < models.size(); ++model_index) {
        const CutMesh &cut_model = cut_models[model_index];
        const CutAOIs &cutAOIs = cuts[model_index];
        for (size_t cut_index = 0; cut_index < cutAOIs.size(); ++cut_index) {
            const CutAOI &cut = cutAOIs[cut_index];
            bbs.push_back(bounding_box(cut, cut_model));
        }
    }

    // keep converted AOI to extend from
    std::vector<std::optional<ExtendAOI>> ecuts(m2i.get_count());

    // NOTE: bad model could have self intersection but this can't solve it
    // find intersection of cuts by Bounding boxes intersection
    size_t index = 0;
    for (size_t model_index = 0; model_index < models.size(); ++model_index) {
        CutAOIs &model_cuts = cuts[model_index];
        for (size_t cut_index = 0; cut_index < model_cuts.size(); ++cut_index, ++index)
        {
            if (!use_cut[index]) continue;
            BoundingBoxf3 &result_bb = bbs[index];
            std::optional<ExtendAOI> &ecut = ecuts[index];
            if (!ecut.has_value()) {
                CutAOI &cut = model_cuts[cut_index];
                const CutMesh &cut_model = cut_models[model_index];
                ecut = create_extend_aoi(cut, cut_model);
            }

            // all differenced models from this model
            std::vector<bool> differenced(models.size(), {false});
            // do not merge itself 
            differenced[model_index] = true;
                        
            // check when exist intersection with result_bb
            size_t index2 = 0;
            for (size_t model_index2 = 0; model_index2 < models.size(); ++model_index2) {
                if (differenced[model_index2]) {
                    if ((model_index2+1) < models.size())
                        index2 = m2i.calc_index({uint32_t(model_index2 + 1), 0});
                    continue; 
                }
                size_t count_cuts = cuts[model_index2].size();
                for (size_t cut_index2 = 0; cut_index2 < count_cuts; ++cut_index2, ++index2){
                    const BoundingBoxf3 &bb = bbs[index2];
                    if (!bb.intersects(result_bb)) continue;
                    priv::diff_model(*ecut, models[model_index2]);
                    differenced[model_index2] = true;
                    if ((model_index2+1) < models.size())
                        index2 = m2i.calc_index({uint32_t(model_index2 + 1), 0});
                    break;
                }
            }
        }    
    }

    // TODO: convert ecuts to CutSurface

    // TODO: merge AOIs without intersection - only append
    
    // Cuts merged in are signed in finished vector as TRUE
    // All rest cuts must be merged simple way
    //SurfaceCut result;
    //for (size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
    //    if (finished[cut_index]) continue;
    //    append(result, std::move(cuts[cut_index]));
    //}
    //return result;
}

bool Slic3r::merge_intersection(SurfaceCut& sc1, const SurfaceCut& sc2) {
    return priv::merge_intersection(sc1, sc2);
}

SurfaceCut priv::merge_intersections(
    SurfaceCuts &cuts, const CutAOIs& cutAOIs, const std::vector<bool> &use_cut)
{
    // create bounding boxes for cuts
    std::vector<BoundingBoxf3> bbs;
    bbs.reserve(cuts.size());
    for (const SurfaceCut &cut : cuts) bbs.push_back(bounding_box(cut));

    // extend used bb by intersecting bb
    // NOTE: after merge 2 cuts could appears
    // new intersection on surface of merged in

    std::vector<bool> finished(cuts.size(), {false});

    // find intersection of cuts by Bounding boxes intersection
    for (size_t cut_index = 0; cut_index < cuts.size(); ++cut_index)
    {
        if (finished[cut_index]) continue;
        if (!use_cut[cut_index]) continue;

        BoundingBoxf3 &result_bb = bbs[cut_index];
        SurfaceCut    &cut       = cuts[cut_index];

        // all merged cuts into cut_index
        std::vector<bool> merged(cuts.size(), {false});

        // merged in last iteration
        std::vector<bool> new_merged;

        bool exist_new_extension;
        bool is_first = true;
        // while exist bb intersection
        do {
            exist_new_extension = false;
            new_merged = std::vector<bool>(cuts.size(), {false});
            // check when exist intersection with result_bb
            for (const BoundingBoxf3 &bb : bbs) {
                size_t bb_index = &bb - &bbs.front();
                // do not merge itself
                if (cut_index == bb_index) continue;
                if (!is_first && merged[bb_index]) continue;
                if (!bb.intersects(result_bb)) {
                    if (is_first) continue;
                    bool has_new_intersection = false;
                    for (size_t i = 0; i < cuts.size(); i++) {
                        if (!new_merged[i]) continue;
                        if (!bbs[i].intersects(bb)) continue;
                        // TODO: check that really intersect by merging
                        has_new_intersection = true;
                    }
                    if (!has_new_intersection) continue;
                }
                if (!priv::merge_intersection(cut, cuts[bb_index])) continue;

                result_bb            = bounding_box(cut);
                merged[bb_index]     = true;
                finished[bb_index]   = true;
                new_merged[bb_index] = true;
                exist_new_extension  = true;
            }
            is_first = false;
        } while (exist_new_extension);
    }

    // Cuts merged in are signed in finished vector as TRUE
    // All rest cuts must be merged simple way
    SurfaceCut result;
    for (size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
        if (finished[cut_index]) continue;
        append(result, std::move(cuts[cut_index]));
    }
    return result;
}

SurfaceCuts priv::create_surface_cuts(const CutAOIs      &cuts,
                                      CutMesh            &mesh,
                                      const ReductionMap &reduction_map)
{    
    // conversion map between vertex index in cgal_model and indices in result
    // used instead of std::map
    // NOTE: can't be used outside because it is rewrited during create_index_triangle_set
    std::string convert_map_name = "v:convert";
    priv::ConvertMap convert_map = mesh.add_property_map<priv::VI, SurfaceCut::Index>(convert_map_name).first;

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

    mesh.remove_property_map(convert_map);
    return result;
}

#ifdef DEBUG_OUTPUT_DIR
// store projection center as circle
void priv::store(const Vec3f       &vertex,
                 const Vec3f       &normal,
                 const std::string &file,
                 float              size)
{
    int    flatten = 20;
    size_t min_i   = 0;
    for (size_t i = 1; i < 3; i++)
        if (normal[min_i] > normal[i]) min_i = i;
    Vec3f up_  = Vec3f::Zero();
    up_[min_i] = 1.f;
    Vec3f side = normal.cross(up_).normalized() * size;
    Vec3f up   = side.cross(normal).normalized() * size;

    indexed_triangle_set its;
    its.vertices.reserve(flatten + 1);
    its.indices.reserve(flatten);

    its.vertices.push_back(vertex);
    its.vertices.push_back(vertex + up);
    size_t max_i = static_cast<size_t>(flatten);
    for (size_t i = 1; i < max_i; i++) {
        float angle = i * 2 * M_PI / flatten;
        Vec3f v     = vertex + sin(angle) * side + cos(angle) * up;
        its.vertices.push_back(v);
        its.indices.emplace_back(0, i, i + 1);
    }
    its.indices.emplace_back(0, flatten, 1);
    its_write_obj(its, file.c_str());
}

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
                 const VCutAOIs            &aois,
                 const CutMeshes           &meshes,
                 const std::string         &file,
                 float                      width)
{
    // create rectangle for each half edge from projection distances
    indexed_triangle_set its;
    its.vertices.reserve(4 * pds.size());
    its.indices.reserve(2 * pds.size());
    for (const ProjectionDistance &pd : pds) {
        if (pd.aoi_index == std::numeric_limits<uint32_t>::max()) continue;
        HI hi  = aois[pd.model_index][pd.aoi_index].second[pd.hi_index];
        const CutMesh &mesh = meshes[pd.model_index];
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
                    unsigned int vi_ = static_cast<unsigned int>(vi);
                    if (vi_ == vi1) continue;
                    if (vi_ == vi2) continue;
                    if (tvi == std::numeric_limits<unsigned int>::max()) {
                        tvi = vi_;
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
            return its.vertices[vi1];
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

void priv::store(const std::vector<indexed_triangle_set> &models,
                 const std::string                       &obj_filename)
{
    indexed_triangle_set merged_model;
    for (const indexed_triangle_set &model : models)
        its_merge(merged_model, model);
    its_write_obj(merged_model, obj_filename.c_str());
}

void priv::store(const std::vector<priv::CutMesh> &models,
                 const std::string                &off_filename)
{
    size_t model_index = 0;
    for (const priv::CutMesh& model : models) {
        std::string filename = off_filename + std::to_string(model_index++) + ".off";
        CGAL::IO::write_OFF(filename, model);
    }
}

// store projection center
void priv::store(const Emboss::IProjection &projection,
                 const Point               &point_to_project,
                 float                      projection_ratio,
                 const std::string         &obj_filename)
{
    auto [front, back] = projection.create_front_back(point_to_project);
    Vec3f diff         = back - front;
    Vec3f pos          = front + diff * projection_ratio;
    priv::store(pos, diff.normalized(),
                DEBUG_OUTPUT_DIR + "projection_center.obj"); // only debug
}

#endif // DEBUG_OUTPUT_DIR
