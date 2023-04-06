#ifndef slic3r_NSVGUtils_hpp_
#define slic3r_NSVGUtils_hpp_

#include "Polygon.hpp"
#include "ExPolygon.hpp"
#include "nanosvg/nanosvg.h"    // load SVG file

// Helper function to work with nano svg
namespace Slic3r {

/// <summary>
/// Convert .svg opened by nanoSvg to Polygons
/// </summary>
/// <param name="image">Opened file</param>
/// <param name="tessTol">Tesselation tolerance 
/// NOTE: Value is in image scale</param>
/// <param name="max_level">Maximal depth</param>
/// <param name="scale">Multiplicator of point coors
/// NOTE: Every point coor from image(float) is multiplied by scale and rounded to integer</param>
/// <param name="is_y_negative">Flag is y negative, when true than y coor is multiplied by -1</param>
/// <returns>Polygons extracted from svg</returns>
Polygons   to_polygons(NSVGimage *image, float tessTol = 10., int max_level = 10, float scale = 1.f, bool is_y_negative = true);
ExPolygons to_expolygons(NSVGimage *image, float tessTol = 10., int max_level = 10, float scale = 1.f, bool is_y_negative = true);

} // namespace Slic3r
#endif // slic3r_NSVGUtils_hpp_
