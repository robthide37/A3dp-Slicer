#include "Voronoi.hpp"

#include "libslic3r/Arachne/utils/PolygonsSegmentIndex.hpp"
#include "libslic3r/Geometry/VoronoiUtils.hpp"
#include "libslic3r/Geometry/VoronoiUtilsCgal.hpp"
#include "libslic3r/MultiMaterialSegmentation.hpp"

#include <boost/log/trivial.hpp>

namespace Slic3r::Geometry {

using PolygonsSegmentIndexConstIt = std::vector<Arachne::PolygonsSegmentIndex>::const_iterator;
using LinesIt                     = Lines::iterator;
using ColoredLinesIt              = ColoredLines::iterator;

// Explicit template instantiation.
template void VoronoiDiagram::construct_voronoi(LinesIt, LinesIt, bool);
template void VoronoiDiagram::construct_voronoi(ColoredLinesIt, ColoredLinesIt, bool);
template void VoronoiDiagram::construct_voronoi(PolygonsSegmentIndexConstIt, PolygonsSegmentIndexConstIt, bool);

template<typename SegmentIterator>
typename boost::polygon::enable_if<
    typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
        typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
    void>::type
VoronoiDiagram::construct_voronoi(const SegmentIterator segment_begin, const SegmentIterator segment_end, const bool try_to_repair_if_needed) {
    boost::polygon::construct_voronoi(segment_begin, segment_end, &m_voronoi_diagram);
    if (try_to_repair_if_needed) {
        if (m_issue_type = detect_known_issues(*this, segment_begin, segment_end); m_issue_type != IssueType::NO_ISSUE_DETECTED) {
            if (m_issue_type == IssueType::MISSING_VORONOI_VERTEX) {
                BOOST_LOG_TRIVIAL(warning) << "Detected missing Voronoi vertex, input polygons will be rotated back and forth.";
            } else if (m_issue_type == IssueType::NON_PLANAR_VORONOI_DIAGRAM) {
                BOOST_LOG_TRIVIAL(warning) << "Detected non-planar Voronoi diagram, input polygons will be rotated back and forth.";
            } else if (m_issue_type == IssueType::VORONOI_EDGE_INTERSECTING_INPUT_SEGMENT) {
                BOOST_LOG_TRIVIAL(warning) << "Detected Voronoi edge intersecting input segment, input polygons will be rotated back and forth.";
            } else if (m_issue_type == IssueType::FINITE_EDGE_WITH_NON_FINITE_VERTEX) {
                BOOST_LOG_TRIVIAL(warning) << "Detected finite Voronoi vertex with non finite vertex, input polygons will be rotated back and forth.";
            } else {
                BOOST_LOG_TRIVIAL(error) << "Detected unknown Voronoi diagram issue, input polygons will be rotated back and forth.";
            }

            if (m_issue_type = try_to_repair_degenerated_voronoi_diagram(segment_begin, segment_end); m_issue_type != IssueType::NO_ISSUE_DETECTED) {
                if (m_issue_type == IssueType::MISSING_VORONOI_VERTEX) {
                    BOOST_LOG_TRIVIAL(error) << "Detected missing Voronoi vertex even after the rotation of input.";
                } else if (m_issue_type == IssueType::NON_PLANAR_VORONOI_DIAGRAM) {
                    BOOST_LOG_TRIVIAL(error) << "Detected non-planar Voronoi diagram even after the rotation of input.";
                } else if (m_issue_type == IssueType::VORONOI_EDGE_INTERSECTING_INPUT_SEGMENT) {
                    BOOST_LOG_TRIVIAL(error) << "Detected Voronoi edge intersecting input segment even after the rotation of input.";
                } else if (m_issue_type == IssueType::FINITE_EDGE_WITH_NON_FINITE_VERTEX) {
                    BOOST_LOG_TRIVIAL(error) << "Detected finite Voronoi vertex with non finite vertex even after the rotation of input.";
                } else {
                    BOOST_LOG_TRIVIAL(error) << "Detected unknown Voronoi diagram issue even after the rotation of input.";
                }

                m_state = State::REPAIR_UNSUCCESSFUL;
            } else {
                m_state = State::REPAIR_SUCCESSFUL;
            }
        } else {
            m_state      = State::REPAIR_NOT_NEEDED;
            m_issue_type = IssueType::NO_ISSUE_DETECTED;
        }
    } else {
        m_state      = State::UNKNOWN;
        m_issue_type = IssueType::UNKNOWN;
    }
}

void VoronoiDiagram::clear()
{
    if (m_is_modified) {
        m_vertices.clear();
        m_edges.clear();
        m_cells.clear();
        m_is_modified = false;
    } else {
        m_voronoi_diagram.clear();
    }

    m_state      = State::UNKNOWN;
    m_issue_type = IssueType::UNKNOWN;
}

template<typename SegmentIterator>
typename boost::polygon::enable_if<
    typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
        typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
    VoronoiDiagram::IssueType>::type
VoronoiDiagram::detect_known_issues(const VoronoiDiagram &voronoi_diagram, SegmentIterator segment_begin, SegmentIterator segment_end)
{
    if (has_finite_edge_with_non_finite_vertex(voronoi_diagram)) {
        return IssueType::FINITE_EDGE_WITH_NON_FINITE_VERTEX;
    } else if (const IssueType cell_issue_type = detect_known_voronoi_cell_issues(voronoi_diagram, segment_begin, segment_end); cell_issue_type != IssueType::NO_ISSUE_DETECTED) {
        return cell_issue_type;
    } else if (!VoronoiUtilsCgal::is_voronoi_diagram_planar_angle(voronoi_diagram, segment_begin, segment_end)) {
        // Detection of non-planar Voronoi diagram detects at least GH issues #8474, #8514 and #8446.
        return IssueType::NON_PLANAR_VORONOI_DIAGRAM;
    }

    return IssueType::NO_ISSUE_DETECTED;
}

template<typename SegmentIterator>
typename boost::polygon::enable_if<
    typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
        typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
    VoronoiDiagram::IssueType>::type
VoronoiDiagram::detect_known_voronoi_cell_issues(const VoronoiDiagram &voronoi_diagram,
                                                 const SegmentIterator segment_begin,
                                                 const SegmentIterator segment_end)
{
    using Segment          = typename std::iterator_traits<SegmentIterator>::value_type;
    using Point            = typename boost::polygon::segment_point_type<Segment>::type;
    using SegmentCellRange = SegmentCellRange<Point>;

    for (VD::cell_type cell : voronoi_diagram.cells()) {
        if (cell.is_degenerate() || !cell.contains_segment())
            continue; // Skip degenerated cell that has no spoon. Also, skip a cell that doesn't contain a segment.

        if (const SegmentCellRange cell_range = VoronoiUtils::compute_segment_cell_range(cell, segment_begin, segment_end); cell_range.is_valid()) {
            // Detection if Voronoi edge is intersecting input segment.
            // It detects this type of issue at least in GH issues #8446, #8474 and #8514.

            const Segment &source_segment      = Geometry::VoronoiUtils::get_source_segment(cell, segment_begin, segment_end);
            const Vec2d    source_segment_from = boost::polygon::segment_traits<Segment>::get(source_segment, boost::polygon::LOW).template cast<double>();
            const Vec2d    source_segment_to   = boost::polygon::segment_traits<Segment>::get(source_segment, boost::polygon::HIGH).template cast<double>();
            const Vec2d    source_segment_vec  = source_segment_to - source_segment_from;

            // All Voronoi vertices must be on the left side of the source segment, otherwise the Voronoi diagram is invalid.
            for (const VD::edge_type *edge = cell_range.edge_begin; edge != cell_range.edge_end; edge = edge->next()) {
                if (edge->is_infinite()) {
                    // When there is a missing Voronoi vertex, we may encounter an infinite Voronoi edge.
                    // This happens, for example, in GH issue #8846.
                    return IssueType::MISSING_VORONOI_VERTEX;
                } else if (const Vec2d edge_v1(edge->vertex1()->x(), edge->vertex1()->y()); Slic3r::cross2(source_segment_vec, edge_v1 - source_segment_from) < 0) {
                    return IssueType::VORONOI_EDGE_INTERSECTING_INPUT_SEGMENT;
                }
            }
        } else {
            // When there is a missing Voronoi vertex (especially at one of the endpoints of the input segment),
            // the returned cell_range is marked as invalid.
            // It detects this type of issue at least in GH issue #8846.
            return IssueType::MISSING_VORONOI_VERTEX;
        }
    }

    return IssueType::NO_ISSUE_DETECTED;
}

bool VoronoiDiagram::has_finite_edge_with_non_finite_vertex(const VoronoiDiagram &voronoi_diagram)
{
    for (const voronoi_diagram_type::edge_type &edge : voronoi_diagram.edges()) {
        if (edge.is_finite()) {
            assert(edge.vertex0() != nullptr && edge.vertex1() != nullptr);
            if (edge.vertex0() == nullptr || edge.vertex1() == nullptr || !VoronoiUtils::is_finite(*edge.vertex0()) || !VoronoiUtils::is_finite(*edge.vertex1()))
                return true;
        }
    }
    return false;
}

template<typename SegmentIterator>
typename boost::polygon::enable_if<
    typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
        typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
    VoronoiDiagram::IssueType>::type
VoronoiDiagram::try_to_repair_degenerated_voronoi_diagram(const SegmentIterator segment_begin, const SegmentIterator segment_end)
{
    return this->m_issue_type;
}

template<typename SegmentIterator>
typename boost::polygon::enable_if<
    typename boost::polygon::gtl_if<typename boost::polygon::is_segment_concept<
        typename boost::polygon::geometry_concept<typename std::iterator_traits<SegmentIterator>::value_type>::type>::type>::type,
    VoronoiDiagram::IssueType>::type
VoronoiDiagram::try_to_repair_degenerated_voronoi_diagram_by_rotation(const SegmentIterator segment_begin,
                                                                      const SegmentIterator segment_end,
                                                                      const double          fix_angle)
{
    return this->m_issue_type;
}

} // namespace Slic3r::Geometry
