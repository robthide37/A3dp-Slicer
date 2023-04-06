#include "NSVGUtils.hpp"
#include "ClipperUtils.hpp"
namespace Slic3r {
// inspired by nanosvgrast.h function nsvgRasterize -> nsvg__flattenShape -> nsvg__flattenCubicBez
// https://github.com/memononen/nanosvg/blob/f0a3e1034dd22e2e87e5db22401e44998383124e/src/nanosvgrast.h#L335
void flatten_cubic_bez(Polygon &polygon, float tessTol, Vec2f p1, Vec2f p2, Vec2f p3, Vec2f p4, int level)
{
    Vec2f p12  = (p1 + p2) * 0.5f;
    Vec2f p23  = (p2 + p3) * 0.5f;
    Vec2f p34  = (p3 + p4) * 0.5f;
    Vec2f p123 = (p12 + p23) * 0.5f;

    Vec2f pd  = p4 - p1;
    Vec2f pd2 = p2 - p4;
    float d2  = std::abs(pd2.x() * pd.y() - pd2.y() * pd.x());
    Vec2f pd3 = p3 - p4;
    float d3  = std::abs(pd3.x() * pd.y() - pd3.y() * pd.x());
    float d23 = d2 + d3;

    if ((d23 * d23) < tessTol * (pd.x() * pd.x() + pd.y() * pd.y())) {
        polygon.points.emplace_back(p4.x(), p4.y());
        return;
    }

    --level;
    if (level == 0)
        return;
    Vec2f p234  = (p23 + p34) * 0.5f;
    Vec2f p1234 = (p123 + p234) * 0.5f;
    flatten_cubic_bez(polygon, tessTol, p1, p12, p123, p1234, level);
    flatten_cubic_bez(polygon, tessTol, p1234, p234, p34, p4, level);
}

Polygons to_polygons(NSVGimage *image, float tessTol, int max_level, bool is_y_negative)
{
    const float y_sign = is_y_negative ? -1.f : 1.f;
    Polygons    polygons;
    for (NSVGshape *shape = image->shapes; shape != NULL; shape = shape->next) {
        if (!(shape->flags & NSVG_FLAGS_VISIBLE))
            continue;
        if (shape->fill.type == NSVG_PAINT_NONE)
            continue;

        Polygon polygon;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
            // Flatten path
            polygon.points.emplace_back(path->pts[0], y_sign * path->pts[1]);
            size_t path_size = (path->npts > 1) ? static_cast<size_t>(path->npts - 1) : 0;
            for (size_t i = 0; i < path_size; i += 3) {
                const float *p = &path->pts[i * 2];
                Vec2f        p1(p[0], p[1]), p2(p[2], p[3]), p3(p[4], p[5]), p4(p[6], p[7]);
                flatten_cubic_bez(polygon, tessTol, p1, p2, p3, p4, max_level);
            }
            if (path->closed && !polygon.empty()) {
                polygons.push_back(polygon);
                polygon = Polygon();
            }
        }
        if (!polygon.empty())
            polygons.push_back(polygon);
    }
    return polygons;
}

ExPolygons to_expolygons(NSVGimage *image, float tessTol, int max_level, bool is_y_negative)
{
    return union_ex(to_polygons(image, tessTol, max_level, is_y_negative));
}

} // namespace Slic3r