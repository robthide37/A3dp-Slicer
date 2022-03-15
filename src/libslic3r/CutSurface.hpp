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
    using Index = unsigned int;
    // cutted surface
    indexed_triangle_set mesh;

    // list of circulated open surface
    std::vector<std::vector<Index>> cut;
};

/// <summary>
/// Merge two surface cuts together
/// Added surface cut will be consumed
/// </summary>
/// <param name="sc">Surface cut to extend</param>
/// <param name="sc_add">Surface cut to consume</param>
void append(SurfaceCut &sc, SurfaceCut &&sc_add);

/// <summary>
/// Cut surface shape from model
/// </summary>
/// <param name="model">Mesh to cut</param>
/// <param name="shapes">Multi shapes to cut from model</param>
/// <param name="projection">Define transformation from 2d shape to 3d</param>
/// <returns>Cutted surface from model</returns>
SurfaceCut cut_surface(const indexed_triangle_set &model,
                       const ExPolygons           &shapes,
                       const Emboss::IProject     &projection);

} // namespace Slic3r
#endif // slic3r_CutSurface_hpp_
