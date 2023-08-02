#ifndef slic3r_NSVGUtils_hpp_
#define slic3r_NSVGUtils_hpp_

#include <memory>
#include <string>
#include <sstream>
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
Polygons   to_polygons(const NSVGimage &image, float tessTol = 10., int max_level = 10, float scale = 1.f, bool is_y_negative = true);
ExPolygons to_expolygons(const NSVGimage &image, float tessTol = 10., int max_level = 10, float scale = 1.f, bool is_y_negative = true);

void bounds(const NSVGimage &image, Vec2f &min, Vec2f &max);

// read text data from file
std::unique_ptr<char[]> read_from_disk(const std::string &path);

using NSVGimage_ptr = std::unique_ptr<NSVGimage, void (*)(NSVGimage*)>;
NSVGimage_ptr nsvgParseFromFile(const std::string &svg_file_path, const char *units = "mm", float dpi = 96.0f);
NSVGimage_ptr nsvgParse(const std::shared_ptr<char[]> file_data, const char *units = "mm", float dpi = 96.0f);

void save(const NSVGimage &image, std::ostream &data);
bool save(const NSVGimage &image, const std::string &svg_file_path);
} // namespace Slic3r
#endif // slic3r_NSVGUtils_hpp_
