#ifndef slic3r_CutSurface_hpp_
#define slic3r_CutSurface_hpp_

#include <vector>
#include <optional>
#include <memory>
#include <admesh/stl.h> // indexed_triangle_set
#include "Polygon.hpp"
#include "ExPolygon.hpp"
#include "Emboss.hpp"

namespace Slic3r{

/// <summary>
/// Address of contour point in ExPolygon
/// </summary>
struct ExPolygonPoint
{
    // Index of Polygon in ExPolygon
    // 0 .. ExPolygon::contour
    // N .. ExPolygon::hole[N-1]
    size_t poly_id;

    // Index of point in Polygon
    size_t index;
};

/// <summary>
/// Address of contour point in ExPolygons
/// </summary>
struct ExPolygonsPoint : public ExPolygonPoint
{
    // Index of ExPolygon in ExPolygons
    size_t expoly_id;
};

/// <summary>
/// Represents cutted surface from object
/// Extend index triangle set by outlines
/// </summary>
struct SurfaceCut : public indexed_triangle_set
{
    // connected cutted surface
    indexed_triangle_set mesh;

    // verticex index(index to mesh vertices)
    using Index = unsigned int;
    using CutType = std::vector<std::vector<Index>>;
    // list of circulated open surface
    CutType cut;

    // conversion map from vertex index to contour point
    // std::map<Index, ExPolygonsPoint> vertex2contour;
};
using SurfaceCuts = std::vector<SurfaceCut>;

/// <summary>
/// Merge two surface cuts together
/// Added surface cut will be consumed
/// </summary>
/// <param name="sc">Surface cut to extend</param>
/// <param name="sc_add">Surface cut to consume</param>
void append(SurfaceCut &sc, SurfaceCut &&sc_add);

/// <summary>
/// Cut surface shape from model.
/// IMPROVE1: It is possible to prefiltrate model triangles befor cut.(AABB)
/// IMPROVE2: Make a cut by quad. Two triangles are possible slower but it is question for profiler.
/// </summary>
/// <param name="model">Mesh to cut</param>
/// <param name="shapes">Multiple shape to cut from model</param>
/// <param name="projection">Define transformation from 2d coordinate of shape to 3d</param>
/// <returns>Cutted surface from model</returns>
SurfaceCuts cut_surface(const indexed_triangle_set &model,
                        const ExPolygons           &shapes,
                        const Emboss::IProject     &projection);

} // namespace Slic3r
#endif // slic3r_CutSurface_hpp_
