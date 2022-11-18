#ifndef slic3r_PrincipalComponents2D_hpp_
#define slic3r_PrincipalComponents2D_hpp_

#include "AABBTreeLines.hpp"
#include "BoundingBox.hpp"
#include "libslic3r.h"
#include <vector>
#include "Polygon.hpp"

namespace Slic3r {

// returns two eigenvectors of the area covered by given polygons. The vectors are sorted by their corresponding eigenvalue, largest first
std::tuple<Vec2d, Vec2d> compute_principal_components(const Polygons &polys, const double unscaled_resolution);

}

#endif