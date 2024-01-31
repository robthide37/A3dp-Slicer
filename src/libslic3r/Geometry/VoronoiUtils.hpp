#ifndef slic3r_VoronoiUtils_hpp_
#define slic3r_VoronoiUtils_hpp_

#include "libslic3r/Geometry/Voronoi.hpp"

using VD = Slic3r::Geometry::VoronoiDiagram;

namespace Slic3r::Geometry {

class VoronoiUtils
{
public:
    static bool is_finite(const VD::vertex_type &vertex);

    template<typename T> static bool is_in_range(double value)
    {
        return double(std::numeric_limits<T>::lowest()) <= value && value <= double(std::numeric_limits<T>::max());
    }

    template<typename T> static bool is_in_range(const VD::vertex_type &vertex)
    {
        return VoronoiUtils::is_finite(vertex) && is_in_range<T>(vertex.x()) && is_in_range<T>(vertex.y());
    }

    template<typename T> static bool is_in_range(const VD::edge_type &edge)
    {
        if (edge.vertex0() == nullptr || edge.vertex1() == nullptr)
            return false;

        return is_in_range<T>(*edge.vertex0()) && is_in_range<T>(*edge.vertex1());
    }
};

} // namespace Slic3r::Geometry

#endif // slic3r_VoronoiUtils_hpp_
