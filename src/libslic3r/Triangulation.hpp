#ifndef libslic3r_Triangulation_hpp_
#define libslic3r_Triangulation_hpp_

#include <vector>
#include <set>
#include <libslic3r/Point.hpp>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/ExPolygon.hpp>

namespace Slic3r {

class Triangulation
{
public:
    Triangulation() = delete;

    // define oriented connection of 2 vertices(defined by its index)
    using HalfEdge  = std::pair<uint32_t, uint32_t>;
    using HalfEdges = std::vector<HalfEdge>;
    using Indices   = std::vector<Vec3i>;

    /// <summary>
    /// Connect points by triangulation to create filled surface by triangles
    /// Input points have to be unique
    /// Inspiration for make unique points is Emboss::dilate_to_unique_points
    /// </summary>
    /// <param name="points">Points to connect</param>
    /// <param name="edges">Constraint for edges, pair is from point(first) to
    /// point(second), sorted lexicographically</param> 
    /// <returns>Triangles</returns>
    static Indices triangulate(const Points &points,
                               const HalfEdges &half_edges);
    static Indices triangulate(const Polygon &polygon);
    static Indices triangulate(const Polygons &polygons);
    static Indices triangulate(const ExPolygons &expolygons);

    /// <summary>
    /// Filter out triagles without both side edge or inside half edges
    /// Main purpose: Filter out triangles which lay outside of ExPolygon
    /// given to triangulation
    /// </summary>
    /// <param name="indices">Triangles</param>
    /// <param name="half_edges">Only outer edges</param>
    static void remove_outer(Indices &indices, const HalfEdges &half_edges);

};

} // namespace Slic3r
#endif // libslic3r_Triangulation_hpp_