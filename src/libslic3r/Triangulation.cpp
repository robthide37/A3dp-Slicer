#include "Triangulation.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

using namespace Slic3r;

namespace Slic3r::Private {

inline void insert_points(Points& points, const Polygon &polygon) {
    points.insert(points.end(), polygon.points.begin(), polygon.points.end());
}

inline void insert_edge(Triangulation::HalfEdges &edges, uint32_t &offset, const Polygon &polygon) {
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

} // namespace Private

Triangulation::Indices Triangulation::triangulate(const Points &   points,
                                                  const HalfEdges &half_edges,
                                                  bool allow_opposite_edge)
{
    // IMPROVE use int point insted of float !!!

    // use cgal triangulation
    using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Itag = CGAL::Exact_predicates_tag;
    using CDT =
        CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>;
    using Point = CDT::Point;

    // construct a constrained triangulation
    CDT                                    cdt;
    std::map<CDT::Vertex_handle, uint32_t> map;             // for indices
    std::vector<CDT::Vertex_handle>        vertices_handle; // for constriants
    vertices_handle.reserve(points.size());
    for (const auto &p : points) {
        Point cdt_p(p.x(), p.y());
        auto handle = cdt.insert(cdt_p);
        vertices_handle.push_back(handle);
        // point index
        uint32_t pi = &p - &points.front();
        map[handle]  = pi;
    }

    // triangle can not contain forbiden edge
    for (const HalfEdge &edge : half_edges) {
        const CDT::Vertex_handle &vh1 = vertices_handle[edge.first];
        const CDT::Vertex_handle &vh2 = vertices_handle[edge.second];
        cdt.insert_constraint(vh1, vh2);
    }

    auto faces = cdt.finite_face_handles();
    std::vector<Vec3i> indices;
    indices.reserve(faces.size());
    for (CDT::Face_handle face : faces) {
        // point indices
        std::array<uint32_t, 3> pi;
        for (size_t i = 0; i < 3; ++i) pi[i] = map[face->vertex(i)];

        // Do not use triangles with opposit edges
        if (!allow_opposite_edge) {
            if (half_edges.find(std::make_pair(pi[1], pi[0])) != half_edges.end()) continue;
            if (half_edges.find(std::make_pair(pi[2], pi[1])) != half_edges.end()) continue;
            if (half_edges.find(std::make_pair(pi[0], pi[2])) != half_edges.end()) continue;
        }

        indices.emplace_back(pi[0], pi[1], pi[2]);
    }
    return indices;
}

Triangulation::Indices Triangulation::triangulate(const Polygon &polygon)
{
    const Points &pts = polygon.points;
    HalfEdges edges;
    uint32_t offset = 0;
    Private::insert_edge(edges, offset, polygon);
    Triangulation::Indices indices = triangulate(pts, edges);
    remove_outer(indices, edges);
    return indices;
}

Triangulation::Indices Triangulation::triangulate(const Polygons &polygons)
{
    size_t count = count_points(polygons);
    Points points;
    points.reserve(count);

    HalfEdges edges;
    uint32_t  offset = 0;

    for (const Polygon &polygon : polygons) {
        Private::insert_points(points, polygon);
        Private::insert_edge(edges, offset, polygon);
    }

    Triangulation::Indices indices = triangulate(points, edges);
    remove_outer(indices, edges);
    return indices;
}

Triangulation::Indices Triangulation::triangulate(const ExPolygons &expolygons)
{
    size_t count = count_points(expolygons);
    Points points;
    points.reserve(count);

    HalfEdges edges;
    uint32_t offset = 0;

    for (const ExPolygon &expolygon : expolygons) { 
        Private::insert_points(points, expolygon.contour);
        Private::insert_edge(edges, offset, expolygon.contour);
        for (const Polygon &hole : expolygon.holes) {
            Private::insert_points(points, hole);
            Private::insert_edge(edges, offset, hole);        
        }
    }        

    Triangulation::Indices indices = triangulate(points, edges);
    remove_outer(indices, edges);
    return indices;
}

void Triangulation::remove_outer(Indices &indices, const HalfEdges &half_edges)
{
    // convert half edge to triangle index
    std::map<HalfEdge, uint32_t> edge2triangle;
    // triangles with all edges out of half_edge, candidate to remove
    std::vector<uint32_t> triangles_to_check;
    triangles_to_check.reserve(indices.size() / 3);

    // triangle half edge
    auto create_half_edge = [](const Vec3i& triangle, size_t index) {
        size_t next_index = (index == 2) ? 0 : (index + 1);
        return HalfEdge(triangle[index], triangle[next_index]);
    };

    // neighbor triangle half edge
    auto create_opposit_half_edge = [](const Vec3i &triangle, size_t index) {
        size_t prev_index = (index == 0) ? 2 : (index - 1);
        return HalfEdge(triangle[index], triangle[prev_index]);
    };

    for (const auto &t : indices) {
        uint32_t index     = &t - &indices.front();
        bool     is_border = false;
        for (size_t j = 0; j < 3; ++j) {
            HalfEdge he = create_half_edge(t, j);
            if (half_edges.find(he) != half_edges.end())
                is_border = true;
            else {
                // half edge already exists
                assert(edge2triangle.find(he) == edge2triangle.end());
                edge2triangle[he] = index;
            }
        }
        if (!is_border) triangles_to_check.push_back(index);
    }

    std::set<uint32_t>   remove;
    std::queue<uint32_t> insert;
    for (uint32_t index : triangles_to_check) {
        auto it = remove.find(index);
        if (it != remove.end()) continue; // already removed

        bool is_edge = false;
        const Vec3i &t = indices[index];
        // check if exist opposit triangle
        for (size_t j = 0; j < 3; ++j) {
            HalfEdge he = create_opposit_half_edge(t, j);
            auto it = edge2triangle.find(he);
            if (it == edge2triangle.end()) {
                is_edge = true;
                break;
            }
            assert(it->second != index);
        }
        if (!is_edge) continue; // correct
        // when there is no opposit edge, remove triangle and all neighbors recursive

        insert.push(index);
        do {
            uint32_t triangle_index = insert.front();
            insert.pop();
            if (remove.find(triangle_index) != remove.end()) continue; // already removed
            remove.insert(triangle_index);

            for (size_t j = 0; j < 3; ++j) {
                const Vec3i &t  = indices[triangle_index];
                HalfEdge he = create_opposit_half_edge(t, j);
                auto     it = edge2triangle.find(he);
                if (it == edge2triangle.end()) continue; // edge triangle without neighbor
                // process neighbor
                uint32_t neighbor_index = it->second;
                assert(neighbor_index != triangle_index);
                insert.push(neighbor_index);
            }
        } while (!insert.empty());
    }

    // remove indices
    std::vector<uint32_t> rem(remove.begin(), remove.end());
    std::sort(rem.begin(), rem.end());
    uint32_t offset = 0;
    for (uint32_t i : rem) {
        indices.erase(indices.begin() + (i - offset));
        ++offset;
    }
}