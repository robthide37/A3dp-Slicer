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
    // connected cutted surface --> inheritance is used
    //indexed_triangle_set mesh;

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

// call private function with same name to test it
bool merge_intersection(SurfaceCut &sc1, const SurfaceCut &sc2);


/// <summary>
/// Merge surface cuts int one
/// </summary>
/// <param name="cuts">input</param>
SurfaceCut merge(SurfaceCuts&& cuts);

/// <summary>
/// Cut surface shape from models.
/// </summary>
/// <param name="shapes">Multiple shape to cut from model</param>
/// <param name="models">Multi mesh to cut, need to be in same coordinate system</param>
/// <param name="projection">Define transformation 2d shape into 3d</param>
/// <param name="projection_ratio">Define ideal ratio between front and back projection to cut
/// 0 .. means use closest to front projection
/// 1 .. means use closest to back projection
/// value from <0, 1>
/// </param>
/// <returns>Cutted surface from model</returns>
SurfaceCut cut_surface(const ExPolygons                        &shapes,
                       const std::vector<indexed_triangle_set> &models,
                       const Emboss::IProjection               &projection,
                       float projection_ratio);

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
