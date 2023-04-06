#ifndef slic3r_NSVGUtils_hpp_
#define slic3r_NSVGUtils_hpp_

#include "Polygon.hpp"
#include "ExPolygon.hpp"
#include "nanosvg/nanosvg.h"    // load SVG file

// Helper function to work with nano svg
namespace Slic3r { 

/// <summary>
/// Convert cubic curve to lines
/// Inspired by nanosvgrast.h function nsvgRasterize->nsvg__flattenShape 
/// </summary>
/// <param name="polygon">Result points</param>
/// <param name="tessTol">Tesselation tolerance</param>
/// <param name="p1">Curve point</param>
/// <param name="p2">Curve point</param>
/// <param name="p3">Curve point</param>
/// <param name="p4">Curve point</param>
/// <param name="level"></param>
void flatten_cubic_bez(Polygon &polygon, float tessTol, Vec2f p1, Vec2f p2, Vec2f p3, Vec2f p4, int level);

/// <summary>
/// Convert .svg opened by nanoSvg to Polygons
/// </summary>
/// <param name="image">Opened file</param>
/// <param name="tessTol">Tesselation tolerance</param>
/// <param name="max_level">Maximal depth</param>
/// <param name="is_y_negative">Flag is y negative, when true than y coor is multiplied by -1</param>
/// <returns>Polygons extracted from svg</returns>
Polygons   to_polygons(NSVGimage *image, float tessTol = 10., int max_level = 10, bool is_y_negative = true);
ExPolygons to_expolygons(NSVGimage *image, float tessTol = 10., int max_level = 10, bool is_y_negative = true);

} // namespace Slic3r
#endif // slic3r_NSVGUtils_hpp_
