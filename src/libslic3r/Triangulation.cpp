#include "Triangulation.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

using namespace Slic3r;

namespace Slic3r::Private {

inline void insert_points(Points& points, const Polygon &polygon) {
    points.insert(points.end(), polygon.points.begin(), polygon.points.end());
}

inline void insert_edges(Triangulation::HalfEdges &edges, uint32_t &offset, const Polygon &polygon) {
    const Points &pts = polygon.points;
    for (uint32_t i = 1; i < pts.size(); ++i) {
        uint32_t i2 = i + offset;
        edges.push_back({i2 - 1, i2});
    }
    uint32_t size = static_cast<uint32_t>(pts.size());
    // add connection from first to last point
    edges.push_back({offset + size - 1, offset});
    offset += size;
}

} // namespace Private

Triangulation::Indices Triangulation::triangulate(const Points    &points,
                                                  const HalfEdges &constrained_half_edges)
{
    // IMPROVE use int point insted of float !!!

    // use cgal triangulation
    using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb   = CGAL::Triangulation_vertex_base_with_info_2<uint32_t, K>;
    using Fb   = CGAL::Constrained_triangulation_face_base_2<K>;
    using Tds  = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using CDT =
        CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;

    // construct a constrained triangulation
    CDT cdt;
    {
        std::vector<CDT::Vertex_handle> vertices_handle; // for constriants
        vertices_handle.reserve(points.size());
        for (const auto &p : points) {
            uint32_t pi = &p - &points.front();
            auto  handle = cdt.insert({ p.x(), p.y() });
            handle->info() = pi;
            vertices_handle.push_back(handle);
        }
        // Constrain the triangulation.
        for (const HalfEdge &edge : constrained_half_edges)
            cdt.insert_constraint(vertices_handle[edge.first], vertices_handle[edge.second]);
    }

    auto faces = cdt.finite_face_handles();

    // Unmark constrained edges of outside faces.
    size_t num_faces = 0;
    for (CDT::Face_handle fh : faces) {
        for (int i = 0; i < 3; ++ i) {
            if (fh->is_constrained(i)) {
                auto key = std::make_pair(fh->vertex((i + 2) % 3)->info(), fh->vertex((i + 1) % 3)->info());
                if (auto it = std::lower_bound(constrained_half_edges.begin(), constrained_half_edges.end(), key); it != constrained_half_edges.end() && *it == key) {
                    // This face contains a constrained edge and it is outside.
                    for (int j = 0; j < 3; ++ j)
                        fh->set_constraint(j, false);
                    goto end;
                }
            }
        }
        ++ num_faces;
    end:;
    }

    // Propagate inside the constrained regions.
    std::vector<CDT::Face_handle> queue;
    queue.reserve(num_faces);
    auto inside = [](CDT::Face_handle &fh) { return fh->neighbor(0) != fh && (fh->is_constrained(0) || fh->is_constrained(1) || fh->is_constrained(2)); };
    for (CDT::Face_handle seed : faces)
        if (inside(seed)) {
            // Seed fill to neighbor faces.
            queue.emplace_back(seed);
            while (! queue.empty()) {
                CDT::Face_handle fh = queue.back();
                queue.pop_back();
                for (int i = 0; i < 3; ++ i)
                    if (! fh->is_constrained(i)) {
                        // Propagate along this edge.
                        fh->set_constraint(i, true);
                        CDT::Face_handle nh = fh->neighbor(i);
                        bool was_inside = inside(nh);
                        // Mark the other side of this edge.
                        nh->set_constraint(nh->index(fh), true);
                        if (! was_inside)
                            queue.push_back(nh);
                    }
            }
        }

    std::vector<Vec3i> indices;
    indices.reserve(num_faces);
    for (CDT::Face_handle fh : faces)
        if (inside(fh))
            indices.emplace_back(fh->vertex(0)->info(), fh->vertex(1)->info(), fh->vertex(2)->info());
    return indices;
}

Triangulation::Indices Triangulation::triangulate(const Polygon &polygon)
{
    const Points &pts = polygon.points;
    HalfEdges edges;
    edges.reserve(pts.size());
    uint32_t offset = 0;
    Private::insert_edges(edges, offset, polygon);
    std::sort(edges.begin(), edges.end());
    return triangulate(pts, edges);
}

Triangulation::Indices Triangulation::triangulate(const Polygons &polygons)
{
    size_t count = count_points(polygons);
    Points points;
    points.reserve(count);

    HalfEdges edges;
    edges.reserve(count);
    uint32_t  offset = 0;

    for (const Polygon &polygon : polygons) {
        Private::insert_points(points, polygon);
        Private::insert_edges(edges, offset, polygon);
    }

    std::sort(edges.begin(), edges.end());
    return triangulate(points, edges);
}

Triangulation::Indices Triangulation::triangulate(const ExPolygons &expolygons)
{
    size_t count = count_points(expolygons);
    Points points;
    points.reserve(count);

    HalfEdges edges;
    edges.reserve(count);
    uint32_t offset = 0;

    for (const ExPolygon &expolygon : expolygons) { 
        Private::insert_points(points, expolygon.contour);
        Private::insert_edges(edges, offset, expolygon.contour);
        for (const Polygon &hole : expolygon.holes) {
            Private::insert_points(points, hole);
            Private::insert_edges(edges, offset, hole);        
        }
    }        

    std::sort(edges.begin(), edges.end());
    return triangulate(points, edges);
}
