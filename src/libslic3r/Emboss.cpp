#include "Emboss.hpp"
#include <stdio.h>

#define STB_TRUETYPE_IMPLEMENTATION // force following include to generate implementation
#include "imgui/imstb_truetype.h" // stbtt_fontinfo

using namespace Slic3r;

// do not expose out of this file stbtt_ data types
class Privat
{
public:
    Privat() = delete;

    static std::optional<stbtt_fontinfo> load_font_info(const Emboss::Font &font);

    struct Glyph
    {
        Polygons polygons;
        int      advance_width, left_side_bearing;
    };
    static std::optional<Glyph> get_glyph(stbtt_fontinfo &font_info, int unicode_letter, float flatness = 2.f);
};

std::optional<stbtt_fontinfo> Privat::load_font_info(const Emboss::Font &font)
{
    int font_offset = stbtt_GetFontOffsetForIndex(font.buffer.data(), font.index);
    if (font_offset < 0) {
        std::cerr << "Font index("<<font.index<<") doesn't exist.";
        return {};        
    }
    stbtt_fontinfo font_info;
    if (stbtt_InitFont(&font_info, font.buffer.data(), font_offset) == 0) {
        std::cerr << "Can't initialize font.";
        return {};
    }
    return font_info;
}

std::optional<Privat::Glyph> Privat::get_glyph(stbtt_fontinfo &font_info, int unicode_letter, float flatness)
{
    int glyph_index = stbtt_FindGlyphIndex(&font_info, unicode_letter);
    if (glyph_index == 0) { 
        std::cerr << "Character codepoint(" << unicode_letter 
            << " = '" << (char) unicode_letter << "') is not defined in the font.";
        return {};
    }

    Privat::Glyph glyph;
    stbtt_GetGlyphHMetrics(&font_info, glyph_index, &glyph.advance_width, &glyph.left_side_bearing);

    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&font_info, glyph_index, &vertices);
    if (num_verts <= 0) return glyph; // no shape

    int *         contour_lengths = NULL;
    int           num_countour    = 0;
    stbtt__point *points          = stbtt_FlattenCurves(vertices, num_verts,
                                               flatness,
                                               &contour_lengths,
                                               &num_countour,
                                               font_info.userdata);

    glyph.polygons.reserve(num_countour);
    size_t pi = 0; // point index
    for (size_t ci = 0; ci < num_countour; ++ci) {
        int    length = contour_lengths[ci];
        if (length <= 0) continue;
        Points pts;
        pts.reserve(length);
        for (size_t i = 0; i < length; i++) {
            const stbtt__point &point = points[pi];
            pi++;
            pts.emplace_back(point.x, point.y);
        }
        if (pts.front() == pts.back()) { 
            pts.pop_back(); 
        } else {
            int j = 42;
        }

        glyph.polygons.emplace_back(pts);
    }

    // inner ccw
    // outer cw
    return glyph;
}

std::optional<Emboss::Font> Emboss::load_font(const char *file_path)
{
    FILE *file = fopen(file_path, "rb");
    if (file == nullptr) {
        std::cerr << "Couldn't open " << file_path << " for reading.";
        return {};
    }

    // find size of file
    if (fseek(file, 0L, SEEK_END) != 0) {
        std::cerr << "Couldn't fseek file " << file_path << " for size measure.";
        return {};
    }
    size_t size = ftell(file);
    if (size == 0) {
        std::cerr << "Size of font file is zero. Can't read.";
        return {};    
    }
    rewind(file);

    Font res;
    res.buffer = std::vector<unsigned char>(size);
    size_t count_loaded_bytes = fread((void *) &res.buffer.front(), 1, size, file);

    unsigned int index = 0;
    int font_offset = 0;
    while (font_offset >= 0) {
        font_offset = stbtt_GetFontOffsetForIndex(res.buffer.data(), index++);
    }
    // at least one font must be inside collection
    if (index < 1) {
        std::cerr << "There is no font collection inside file.";
        return {};        
    }
    // select default font on index 0
    res.index = 0;
    res.count = index;

    auto font_info = Privat::load_font_info(res);
    if (!font_info.has_value()) return {};

    // load information about line gap
    stbtt_GetFontVMetrics(&(*font_info), &res.ascent, &res.descent, &res.linegap);
    return res;
}

Polygons Emboss::letter2polygons(const Font &font, char letter)
{
    auto font_info_opt = Privat::load_font_info(font);
    if (!font_info_opt.has_value()) return Polygons();
    stbtt_fontinfo *font_info = &(*font_info_opt);

    auto glyph_opt = Privat::get_glyph(*font_info_opt, (int) letter, font.flatness);
    if (!glyph_opt.has_value()) return Polygons();

    return glyph_opt->polygons;
}

Polygons Emboss::text2polygons(const Font &font, const std::string &text)
{
    auto font_info_opt = Privat::load_font_info(font);
    if (!font_info_opt.has_value()) return Polygons();
    stbtt_fontinfo *font_info = &(*font_info_opt);

    Point    cursor(0, 0);
    Polygons result;
    for (const char &letter : text) {
        if (letter == '\0') break;

        auto glyph_opt = Privat::get_glyph(*font_info_opt, (int) letter, font.flatness);
        if (!glyph_opt.has_value()) continue;

        // move glyph to cursor position
        Polygons polygons = glyph_opt->polygons; // copy
        for (Polygon &polygon : polygons) 
            for (Point &p : polygon.points) p += cursor;
        
        cursor.x() += glyph_opt->advance_width;

        polygons_append(result, polygons);
    }
    return result;
}

std::vector<Vec3i> its_create_neighbors_index_2(const indexed_triangle_set &its)
{
    std::vector<Vec3i> out(its.indices.size(), Vec3i(-1, -1, -1));

    // Create a mapping from triangle edge into face.
    struct EdgeToFace
    {
        // Index of the 1st vertex of the triangle edge. vertex_low <= vertex_high.
        int vertex_low;
        // Index of the 2nd vertex of the triangle edge.
        int vertex_high;
        // Index of a triangular face.
        int face;
        // Index of edge in the face, starting with 1. Negative indices if the
        // edge was stored reverse in (vertex_low, vertex_high).
        int  face_edge;
        bool operator==(const EdgeToFace &other) const
        {
            return vertex_low == other.vertex_low &&
                   vertex_high == other.vertex_high;
        }
        bool operator<(const EdgeToFace &other) const
        {
            return vertex_low < other.vertex_low ||
                   (vertex_low == other.vertex_low &&
                    vertex_high < other.vertex_high);
        }
    };
    std::vector<EdgeToFace> edges_map;
    edges_map.assign(its.indices.size() * 3, EdgeToFace());
    for (uint32_t facet_idx = 0; facet_idx < its.indices.size(); ++facet_idx)
        for (int i = 0; i < 3; ++i) {
            EdgeToFace &e2f = edges_map[facet_idx * 3 + i];
            e2f.vertex_low  = its.indices[facet_idx][i];
            e2f.vertex_high = its.indices[facet_idx][(i + 1) % 3];
            e2f.face        = facet_idx;
            // 1 based indexing, to be always strictly positive.
            e2f.face_edge = i + 1;
            if (e2f.vertex_low > e2f.vertex_high) {
                // Sort the vertices
                std::swap(e2f.vertex_low, e2f.vertex_high);
                // and make the face_edge negative to indicate a flipped edge.
                e2f.face_edge = -e2f.face_edge;
            }
        }
    std::sort(edges_map.begin(), edges_map.end());

    // Assign a unique common edge id to touching triangle edges.
    int num_edges = 0;
    for (size_t i = 0; i < edges_map.size(); ++i) {
        EdgeToFace &edge_i = edges_map[i];
        if (edge_i.face == -1)
            // This edge has been connected to some neighbor already.
            continue;
        // Unconnected edge. Find its neighbor with the correct orientation.
        size_t j;
        bool   found = false;
        for (j = i + 1; j < edges_map.size() && edge_i == edges_map[j]; ++j)
            if (edge_i.face_edge * edges_map[j].face_edge < 0 &&
                edges_map[j].face != -1) {
                // Faces touching with opposite oriented edges and none of the
                // edges is connected yet.
                found = true;
                break;
            }
        if (!found) {
            // FIXME Vojtech: Trying to find an edge with equal orientation.
            // This smells.
            // admesh can assign the same edge ID to more than two facets (which is
            // still topologically correct), so we have to search for a
            // duplicate of this edge too in case it was already seen in this
            // orientation
            for (j = i + 1; j < edges_map.size() && edge_i == edges_map[j];
                 ++j)
                if (edges_map[j].face != -1) {
                    // Faces touching with equally oriented edges and none of
                    // the edges is connected yet.
                    found = true;
                    break;
                }
        }
        // Assign an edge index to the 1st face.
        //        out[edge_i.face](std::abs(edge_i.face_edge) - 1) = num_edges;
        if (found) {
            EdgeToFace &edge_j                               = edges_map[j];
            out[edge_i.face](std::abs(edge_i.face_edge) - 1) = edge_j.face;
            out[edge_j.face](std::abs(edge_j.face_edge) - 1) = edge_i.face;
            // Mark the edge as connected.
            edge_j.face = -1;
        }
        ++num_edges;
    }
    return out;
}

void its_remove_edge_triangles(indexed_triangle_set &its)
{
    ////                       start, count
    //std::vector<std::pair<uint32_t, uint32_t>> neighbors_vertices(its.vertices.size(), {0,0});
    //// calc counts
    //for (const auto &i : its.indices) {
    //    for (size_t j = 0; j < 3; j++) 
    //        ++neighbors_vertices[i[j]].second;
    //}
    //uint32_t triangle_start = 0;
    //for (auto &neighbor : neighbors_vertices) {
    //    neighbor.first = triangle_start;
    //    triangle_start += neighbor.second; 
    //    neighbor.second = 0;
    //}
    //std::vector<uint32_t> neighbors_data(its.indices.size()*3);
    //for (const auto &i : its.indices) {
    //    uint32_t index = &i - &its.indices.front();
    //    for (size_t j = 0; j < 3; j++) {
    //        auto & neighbor = neighbors_vertices[i[j]];
    //        size_t index_data = neighbor.second + neighbor.first;
    //        neighbors_data[index_data] = index;
    //        ++neighbor.second;
    //    }
    //}
    
    auto neighbors = its_create_neighbors_index_2(its);
    std::set<uint32_t> remove;
    std::queue<uint32_t> insert;
    int no_value = -1;
    for (const auto &neighbor : neighbors) {
        uint32_t index = &neighbor - &neighbors.front();
        auto     it    = remove.find(index);
        if (it != remove.end()) continue; // already removed

        if (neighbor[0] != no_value && 
            neighbor[1] != no_value &&
            neighbor[2] != no_value)
            continue;

        insert.push(index);
        while (!insert.empty()) {
            uint32_t i = insert.front();
            insert.pop();
            if (remove.find(i) != remove.end()) continue;
            remove.insert(i);
            for (size_t j = 0; j < 3; j++) {
                if (neighbor[j] == no_value) continue;
                uint32_t i2 = static_cast<uint32_t>(neighbor[j]);
                insert.push(i2);
            }            
        }        
    }
    std::vector<uint32_t> rem(remove.begin(), remove.end());
    std::sort(rem.begin(), rem.end());
    uint32_t offset = 0;
    for (uint32_t i : rem) { 
        its.indices.erase(its.indices.begin() + i - offset);
        ++offset;
    }
}

indexed_triangle_set Emboss::polygons2model(const Polygons &shape2d,
                                            const IProject &projection)
{
    indexed_triangle_set result;
    size_t count_point = count_points(shape2d);
    result.vertices.reserve(2 * count_point);

    std::vector<Vec3f> &front_points = result.vertices;
    std::vector<Vec3f>  back_points;
    back_points.reserve(count_point);

    for (const Polygon &polygon : shape2d) {
        for (const Point &p : polygon.points) {
            auto p2 = projection.project(p);
            front_points.emplace_back(p2.first);
            back_points.emplace_back(p2.second);
        }
    }
    // insert back points, front are already in
    result.vertices.insert(result.vertices.end(),
                           std::make_move_iterator(back_points.begin()),
                           std::make_move_iterator(back_points.end()));

    // CW order of triangle indices
    std::vector<Vec3i> shape_triangles = triangulate(shape2d);
    result.indices.reserve(shape_triangles.size() * 2 + count_point * 2);
    // top triangles - change to CCW
    for (const Vec3i &t : shape_triangles)
        result.indices.emplace_back(t.x(), t.z(), t.y());
    // bottom triangles - use CW
    for (const Vec3i &t : shape_triangles)
        result.indices.emplace_back(t.x() + count_point, t.y() + count_point,
                                    t.z() + count_point);

    // quads around - zig zag by triangles
    size_t polygon_offset = 0;
    for (const Polygon &polygon : shape2d) {
        uint32_t polygon_points = polygon.points.size();
        for (uint32_t p = 0; p < polygon_points; p++) { 
            uint32_t i = polygon_offset + p;
            // previous index
            uint32_t ip = (p == 0) ? (polygon_offset + polygon_points - 1) : (i - 1);
            // bottom indices
            uint32_t i2  = i + count_point;
            uint32_t ip2 = ip + count_point;

            result.indices.emplace_back(i, i2, ip);
            result.indices.emplace_back(ip2, ip, i2);
        }
        polygon_offset += polygon_points;
    }

    // remove bad triangulated faces
    its_remove_edge_triangles(result);
    return result;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
std::vector<Vec3i> Emboss::triangulate(
    const Points& points,
    const std::set<std::pair<uint32_t, uint32_t>> &edges)
{
    // IMPROVE use int point insted of float !!!

    // use cgal triangulation
    using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Itag = CGAL::Exact_predicates_tag;
    using CDT  = CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>;
    using Point = CDT::Point;

    // construct a constrained triangulation
    CDT                               cdt;
    std::map<CDT::Vertex_handle, int> map;             // for indices
    std::vector<CDT::Vertex_handle>   vertices_handle; // for constriants
    vertices_handle.reserve(points.size());
    for (const auto& p: points) {
        Point  cdt_p(p.x(), p.y());
        auto   handl = cdt.insert(cdt_p);
        vertices_handle.push_back(handl);
        size_t i = &p - &points.front();
        map[handl] = i;
    }

    // triangle can not contain forbiden edge
    for (const std::pair<uint32_t, uint32_t> &edge : edges) {
        const CDT::Vertex_handle& vh1 = vertices_handle[edge.first];
        const CDT::Vertex_handle& vh2 = vertices_handle[edge.second];
        cdt.insert_constraint(vh1, vh2);
    }

    auto               faces = cdt.finite_face_handles();
    std::vector<Vec3i> indices;
    indices.reserve(faces.size());
    for (CDT::Face_handle face : faces) {
        auto v0 = face->vertex(0);
        auto v1 = face->vertex(1);
        auto v2 = face->vertex(2);
        uint32_t i0 = map[v0];
        uint32_t i1 = map[v1];
        uint32_t i2 = map[v2];

        // check forbiden triangle edge - opposit order
        if (edges.find(std::make_pair(i0, i1)) != edges.end()) continue;
        if (edges.find(std::make_pair(i1, i2)) != edges.end()) continue;
        if (edges.find(std::make_pair(i2, i0)) != edges.end()) continue;

        indices.emplace_back(map[v0], map[v1], map[v2]);
    }
    return indices;
}

std::vector<Vec3i> Emboss::triangulate(const Polygon &polygon)
{
    const Points &                          pts = polygon.points;
    std::set<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t i = 1; i < pts.size(); ++i) edges.insert({i - 1, i});
    edges.insert({(uint32_t)pts.size() - 1, uint32_t(0)});
    return triangulate(pts, edges);
}

std::vector<Vec3i> Emboss::triangulate(const Polygons &polygons)
{
    size_t count = count_points(polygons);
    Points points;
    points.reserve(count);
    for (const Polygon &polygon : polygons)
        points.insert(points.end(), polygon.points.begin(),
                      polygon.points.end());

    std::set<std::pair<uint32_t, uint32_t>> edges;
    uint32_t offset = 0;
    for (const Polygon& polygon : polygons) {
        const Points &pts = polygon.points;
        for (uint32_t i = 1; i < pts.size(); ++i) {
            uint32_t i2 = i + offset;
            edges.insert({i2 - 1, i2}); 
        }
        uint32_t size = static_cast<uint32_t>(pts.size());
        // add connection from first to last point
        edges.insert({offset + size - 1, offset});
        offset += size;
    }
    return triangulate(points, edges);
}

std::pair<Vec3f, Vec3f> Emboss::ProjectZ::project(const Point &p) const
{
    Vec3f front(p.x(),p.y(),0.f);
    Vec3f back = front; // copy
    back.z() = m_depth;
    return std::make_pair(front, back);
}
