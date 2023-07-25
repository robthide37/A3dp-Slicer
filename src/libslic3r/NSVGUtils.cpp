#include "NSVGUtils.hpp"
#include "ClipperUtils.hpp"

namespace {
Slic3r::Polygons to_polygons(const NSVGshape &shape, float tessTol, int max_level, float scale, bool is_y_negative);
} // namespace

namespace Slic3r {

Polygons to_polygons(const NSVGimage &image, float tessTol, int max_level, float scale, bool is_y_negative)
{
    Polygons polygons;
    for (NSVGshape *shape = image.shapes; shape != NULL; shape = shape->next)
        polygons_append(polygons, ::to_polygons(*shape, tessTol, max_level, scale, is_y_negative));
    return polygons;
}

ExPolygons to_expolygons(const NSVGimage &image, float tessTol, int max_level, float scale, bool is_y_negative){
    ExPolygons expolygons;
    for (NSVGshape *shape = image.shapes; shape != NULL; shape = shape->next) {
        Polygons polygons = ::to_polygons(*shape, tessTol, max_level, scale, is_y_negative);
        if (polygons.empty())
            continue;
        expolygons_append(expolygons, union_ex(polygons));
    }
    return union_ex(expolygons);
}

NSVGimage_ptr nsvgParseFromFile(const std::string &filename, const char *units, float dpi)
{
    NSVGimage *image = ::nsvgParseFromFile(filename.c_str(), units, dpi);
    return {image, ::nsvgDelete};
}

bool save(const NSVGimage &image, const std::string &svg_file_path) 
{
    //BoundingBox bb = get_extents(to_polygons(image));
    FILE * file = boost::nowide::fopen(svg_file_path.c_str(), "w");
    if (file == NULL)
        return false;

    fprintf(file, "<svg xmlns=\"http://www.w3.org/2000/svg\">\n");    
    for (const NSVGshape *shape = image.shapes; shape != NULL; shape = shape->next) {
        std::string d = "M "; // move on start point
        for (const NSVGpath *path = shape->paths; path != NULL; path = path->next) {
            size_t path_size = (path->npts > 1) ? static_cast<size_t>(path->npts - 1) : 0;
            for (size_t i = 0; i < path_size; i += 3) {
                const float *p = &path->pts[i * 2];
                d += std::to_string(p[0]) + "," + std::to_string(p[1]) + " ";
            }
        }
        d += "Z"; // closed path
        fprintf(file, "   <path fill=\"#D2D2D2\" d=\"%s\" />\n", d.c_str());
    }
    fprintf(file, "</svg>\n");
    fclose(file);
}
} // namespace Slic3r

namespace {
using namespace Slic3r; // Polygon + Vec2f

bool is_useable(const NSVGshape &shape)
{
    if (!(shape.flags & NSVG_FLAGS_VISIBLE))
        return false;
    if (shape.fill.type == NSVG_PAINT_NONE)
        return false;
    return true;
}

Point::coord_type to_coor(float val, float scale) { return static_cast<Point::coord_type>(std::round(val * scale)); }



bool need_flattening(float tessTol, const Vec2f &p1, const Vec2f &p2, const Vec2f &p3, const Vec2f &p4) {
    // f .. first
    // s .. second
    auto det = [](const Vec2f &f, const Vec2f &s) {
        return std::fabs(f.x() * s.y() - f.y() * s.x()); 
    };

    Vec2f pd  = (p4 - p1);
    Vec2f pd2 = (p2 - p4);
    float d2  = det(pd2, pd);
    Vec2f pd3 = (p3 - p4);
    float d3  = det(pd3, pd);
    float d23 = d2 + d3;

    return (d23 * d23) >= tessTol * pd.squaredNorm();
}

// see function nsvg__lineTo(NSVGparser* p, float x, float y)
bool is_line(const float *p, float precision = 1e-4){
    //Vec2f p1(p[0], p[1]);
    //Vec2f p2(p[2], p[3]);
    //Vec2f p3(p[4], p[5]);
    //Vec2f p4(p[6], p[7]);
    float dx_3 = (p[6] - p[0]) / 3.f;
    float dy_3 = (p[7] - p[1]) / 3.f;

    return 
        is_approx(p[2], p[0] + dx_3, precision) && 
        is_approx(p[4], p[6] - dx_3, precision) && 
        is_approx(p[3], p[1] + dy_3, precision) &&
        is_approx(p[5], p[7] - dy_3, precision);
}

/// <summary>
/// Convert cubic curve to lines
/// Inspired by nanosvgrast.h function nsvgRasterize -> nsvg__flattenShape -> nsvg__flattenCubicBez
/// https://github.com/memononen/nanosvg/blob/f0a3e1034dd22e2e87e5db22401e44998383124e/src/nanosvgrast.h#L335
/// </summary>
/// <param name="polygon">Result points</param>
/// <param name="tessTol">Tesselation tolerance</param>
/// <param name="p1">Curve point</param>
/// <param name="p2">Curve point</param>
/// <param name="p3">Curve point</param>
/// <param name="p4">Curve point</param>
/// <param name="level">Actual depth of recursion</param>
void flatten_cubic_bez(Polygon &polygon, float tessTol, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3, const Vec2f& p4, int level)
{
    if (!need_flattening(tessTol, p1, p2, p3, p4)) {
        Point::coord_type x = static_cast<Point::coord_type>(std::round(p4.x()));
        Point::coord_type y = static_cast<Point::coord_type>(std::round(p4.y()));
        polygon.points.emplace_back(x, y);
        return;
    }

    --level;
    if (level == 0)
        return;
    
    Vec2f p12  = (p1 + p2) * 0.5f;
    Vec2f p23  = (p2 + p3) * 0.5f;
    Vec2f p34  = (p3 + p4) * 0.5f;
    Vec2f p123 = (p12 + p23) * 0.5f;
    Vec2f p234  = (p23 + p34) * 0.5f;
    Vec2f p1234 = (p123 + p234) * 0.5f;
    flatten_cubic_bez(polygon, tessTol, p1, p12, p123, p1234, level);
    flatten_cubic_bez(polygon, tessTol, p1234, p234, p34, p4, level);
}

Polygons to_polygons(NSVGpath *first_path, float tessTol, int max_level, float scale)
{
    Polygons polygons;
    Polygon  polygon;
    for (NSVGpath *path = first_path; path != NULL; path = path->next) {
        // Flatten path
        Point::coord_type x = to_coor(path->pts[0], scale);
        Point::coord_type y = to_coor(path->pts[1], scale);
        polygon.points.emplace_back(x, y);
        size_t path_size = (path->npts > 1) ? static_cast<size_t>(path->npts - 1) : 0;
        for (size_t i = 0; i < path_size; i += 3) {
            const float *p = &path->pts[i * 2];
            if (is_line(p)) {
                // point p4
                Point::coord_type x = to_coor(p[6], scale);
                Point::coord_type y = to_coor(p[7], scale);
                polygon.points.emplace_back(x, y);
                continue;
            }
            Vec2f p1(p[0], p[1]);
            Vec2f p2(p[2], p[3]);
            Vec2f p3(p[4], p[5]);
            Vec2f p4(p[6], p[7]);
            flatten_cubic_bez(polygon, tessTol, p1 * scale, p2 * scale, p3 * scale, p4 * scale, max_level);
        }
        if (path->closed && !polygon.empty()) {
            polygons.push_back(polygon);
            polygon = Polygon();
        }
    }
    if (!polygon.empty())
        polygons.push_back(polygon);

    remove_same_neighbor(polygons);
    return polygons;
}

Polygons to_polygons(const NSVGshape &shape, float tessTol, int max_level, float scale, bool is_y_negative)
{
    if (!is_useable(shape))
        return {};

    Polygons polygons = to_polygons(shape.paths, tessTol, max_level, scale);
    
    if (is_y_negative)
        for (Polygon &polygon : polygons)
            for (Point &p : polygon.points)
                p.y() = -p.y();

    return polygons;
}

} // namespace