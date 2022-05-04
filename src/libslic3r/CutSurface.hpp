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
/// Represents cutted surface from object
/// Extend index triangle set by outlines
/// </summary>
struct SurfaceCut : public indexed_triangle_set
{
    // connected cutted surface
    indexed_triangle_set mesh;

    // vertex indices(index to mesh vertices)
    using Index = unsigned int;
    using CutContour = std::vector<std::vector<Index>>;
    // list of circulated open surface
    CutContour contours;

    // Conversion map from vertex index to contour point
    // Could be used for filtration of surface cuts
    // Still I don't have an idea how to filtrate it.
    // What is wanted result on wave?
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
/// Merge surface cuts int one
/// </summary>
/// <param name="cuts">input</param>
SurfaceCut merge(SurfaceCuts&& cuts);

/// <summary>
/// Cut surface shape from model.
/// IMPROVE1: It is possible to prefiltrate model triangles befor cut.(AABB)
/// IMPROVE2: Make a cut by quad. Two triangles are possible slower but it is question for profiler.
/// </summary>
/// <param name="model">Mesh to cut</param>
/// <param name="shapes">Multiple shape to cut from model</param>
/// <param name="projection">Define transformation from 2d coordinate of shape to 3d</param>
/// <returns>Cutted surface from model</returns>
SurfaceCut cut_surface(const indexed_triangle_set &model,
                       const ExPolygons           &shapes,
                       const Emboss::IProjection  &projection);

/// <summary>
/// Create model from surface cuts by projection
/// </summary>
/// <param name="cuts">Surfaces from model</param>
/// <param name="projection">Way of emboss</param>
/// <returns>Mesh</returns>
indexed_triangle_set cuts2model(const SurfaceCuts        &cuts,
                                const Emboss::IProject3f &projection);
/// <summary>
/// Create model from surface cuts by projection
/// </summary>
/// <param name="cut">Surface from model with outlines</param>
/// <param name="projection">Way of emboss</param>
/// <returns>Mesh</returns>
indexed_triangle_set cut2model(const SurfaceCut         &cut,
                               const Emboss::IProject3f &projection);


} // namespace Slic3r
#endif // slic3r_CutSurface_hpp_
