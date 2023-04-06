#include "NSVGUtils.hpp"
#include "ClipperUtils.hpp"

namespace {
using namespace Slic3r; // Polygon + Vec2f
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
/// <param name="level">Actual depth of recursion</param>
/// <param name="scale">Scale of point - multiplicator
/// NOTE: increase preccission by number greater than 1.</param>
void flatten_cubic_bez(Polygon &polygon, float tessTol, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3, const Vec2f& p4, int level, float scale);

Point::coord_type to_coor(float val, float scale) { return static_cast<Point::coord_type>(std::round(val * scale)); }

} // namespace

namespace Slic3r {
Polygons to_polygons(NSVGimage *image, float tessTol, int max_level, float scale, bool is_y_negative)
{
    Polygons polygons;
    for (NSVGshape *shape = image->shapes; shape != NULL; shape = shape->next) {
        if (!(shape->flags & NSVG_FLAGS_VISIBLE))
            continue;
        if (shape->fill.type == NSVG_PAINT_NONE)
            continue;

        Polygon polygon;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
            // Flatten path
            Point::coord_type x = to_coor(path->pts[0], scale);
            Point::coord_type y = to_coor(path->pts[1], scale);
            polygon.points.emplace_back(x, y);
            size_t path_size = (path->npts > 1) ? static_cast<size_t>(path->npts - 1) : 0;
            for (size_t i = 0; i < path_size; i += 3) {
                const float *p = &path->pts[i * 2];
                Vec2f p1(p[0], p[1]);
                Vec2f p2(p[2], p[3]);
                Vec2f p3(p[4], p[5]);
                Vec2f p4(p[6], p[7]);
                flatten_cubic_bez(polygon, tessTol, p1, p2, p3, p4, max_level, scale);
            }
            if (path->closed && !polygon.empty()) {
                polygons.push_back(polygon);
                polygon = Polygon();
            }
        }
        if (!polygon.empty()) 
            polygons.push_back(polygon);
    }

    if (is_y_negative)
        for (Polygon &polygon : polygons)
            for (Point &p : polygon.points)
                p.y() = -p.y();    

    return polygons;
}

ExPolygons to_expolygons(NSVGimage *image, float tessTol, int max_level, float scale, bool is_y_negative){
    return union_ex(to_polygons(image, tessTol, max_level, scale, is_y_negative));
}

} // namespace Slic3r

namespace {
// inspired by nanosvgrast.h function nsvgRasterize -> nsvg__flattenShape -> nsvg__flattenCubicBez
// https://github.com/memononen/nanosvg/blob/f0a3e1034dd22e2e87e5db22401e44998383124e/src/nanosvgrast.h#L335
void flatten_cubic_bez(Polygon &polygon, float tessTol, const Vec2f& p1, const Vec2f& p2, const Vec2f& p3, const Vec2f& p4, int level, float scale)
{
    // f .. first
    // s .. second
    auto det = [](const Vec2f &f, const Vec2f &s) {
        return std::fabs(f.x() * s.y() - f.y() * s.x());
    };

    Vec2f pd  = p4 - p1;
    Vec2f pd2 = p2 - p4;
    float d2  = det(pd2, pd);
    Vec2f pd3 = p3 - p4;
    float d3  = det(pd3, pd);
    float d23 = d2 + d3;

    if ((d23 * d23) < tessTol * pd.squaredNorm()) {
        Point::coord_type x = to_coor(p4.x(), scale);
        Point::coord_type y = to_coor(p4.y(), scale);
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
    flatten_cubic_bez(polygon, tessTol, p1, p12, p123, p1234, level, scale);
    flatten_cubic_bez(polygon, tessTol, p1234, p234, p34, p4, level, scale);
}
} // namespace