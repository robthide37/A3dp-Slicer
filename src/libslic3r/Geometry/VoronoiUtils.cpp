#include "VoronoiUtils.hpp"

namespace Slic3r::Geometry {

bool VoronoiUtils::is_finite(const VD::vertex_type &vertex)
{
    return std::isfinite(vertex.x()) && std::isfinite(vertex.y());
}

} // namespace Slic3r::Geometry