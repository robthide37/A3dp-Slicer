///|/ Copyright (c) Prusa Research 2021 Vojtěch Bubník @bubnikv
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_Geometry_Voronoi_hpp_
#define slic3r_Geometry_Voronoi_hpp_

#include "../Line.hpp"
#include "../Polyline.hpp"

#ifdef _MSC_VER
// Suppress warning C4146 in OpenVDB: unary minus operator applied to unsigned type, result still unsigned
#pragma warning(push)
#pragma warning(disable : 4146)
#endif // _MSC_VER
#include "boost/polygon/voronoi.hpp"
#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

namespace Slic3r::Geometry {

class VoronoiDiagram
{
public:
    using coord_type           = double;
    using voronoi_diagram_type = boost::polygon::voronoi_diagram<coord_type>;
    using point_type           = boost::polygon::point_data<voronoi_diagram_type::coordinate_type>;
    using segment_type         = boost::polygon::segment_data<voronoi_diagram_type::coordinate_type>;
    using rect_type            = boost::polygon::rectangle_data<voronoi_diagram_type::coordinate_type>;

    using coordinate_type = voronoi_diagram_type::coordinate_type;
    using vertex_type     = voronoi_diagram_type::vertex_type;
    using edge_type       = voronoi_diagram_type::edge_type;
    using cell_type       = voronoi_diagram_type::cell_type;

    using const_vertex_iterator = voronoi_diagram_type::const_vertex_iterator;
    using const_edge_iterator   = voronoi_diagram_type::const_edge_iterator;
    using const_cell_iterator   = voronoi_diagram_type::const_cell_iterator;

    using vertex_container_type = voronoi_diagram_type::vertex_container_type;
    using edge_container_type   = voronoi_diagram_type::edge_container_type;
    using cell_container_type   = voronoi_diagram_type::cell_container_type;

    VoronoiDiagram() = default;

    virtual ~VoronoiDiagram() = default;

    void clear() { m_voronoi_diagram.clear(); }

    const cell_container_type &cells() const { return m_voronoi_diagram.cells(); }

    const vertex_container_type &vertices() const { return m_voronoi_diagram.vertices(); }

    const edge_container_type &edges() const { return m_voronoi_diagram.edges(); }

    std::size_t num_cells() const { return m_voronoi_diagram.num_cells(); }

    std::size_t num_edges() const { return m_voronoi_diagram.num_edges(); }

    std::size_t num_vertices() const { return m_voronoi_diagram.num_vertices(); }

    template<typename SegmentIterator>
    typename boost::polygon::enable_if<
        typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
            typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
        void>::type
    construct_voronoi(const SegmentIterator first, const SegmentIterator last)
    {
        boost::polygon::construct_voronoi(first, last, &m_voronoi_diagram);
    }

    template<typename PointIterator>
    typename boost::polygon::enable_if<
        typename boost::polygon::gtl_if<typename boost::polygon::is_point_concept<
            typename boost::polygon::geometry_concept<typename std::iterator_traits<PointIterator>::value_type>::type>::type>::type,
        void>::type
    construct_voronoi(const PointIterator first, const PointIterator last)
    {
        boost::polygon::construct_voronoi(first, last, &m_voronoi_diagram);
    }

    template<typename PointIterator, typename SegmentIterator>
    typename boost::polygon::enable_if<
        typename boost::polygon::gtl_and<
            typename boost::polygon::gtl_if<typename boost::polygon::is_point_concept<
                typename boost::polygon::geometry_concept<typename std::iterator_traits<PointIterator>::value_type>::type>::type>::type,
            typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<typename boost::polygon::geometry_concept<
                typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type>::type,
        void>::type
    construct_voronoi(const PointIterator p_first, const PointIterator p_last, const SegmentIterator s_first, const SegmentIterator s_last)
    {
        boost::polygon::construct_voronoi(p_first, p_last, s_first, s_last, &m_voronoi_diagram);
    }

private:
    voronoi_diagram_type m_voronoi_diagram;
};

} // namespace Slic3r::Geometry

#endif // slic3r_Geometry_Voronoi_hpp_
