#ifndef SRC_LIBSLIC3R_ALGORITHM_REGION_EXPANSION_HPP_
#define SRC_LIBSLIC3R_ALGORITHM_REGION_EXPANSION_HPP_

#include <cstdint>

namespace Slic3r {

class Polygon;
using Polygons = std::vector<Polygon>;
class ExPolygon;
using ExPolygons = std::vector<ExPolygon>;

namespace Algorithm {

std::vector<Polygons> expand_expolygons(const ExPolygons &src, const ExPolygons &boundary,
    // Scaled expansion value
    float expansion, 
    // Expand by waves of expansion_step size (expansion_step is scaled).
    float expansion_step,
    // Don't take more than max_nr_steps for small expansion_step.
    size_t max_nr_steps);

} // Algorithm
} // Slic3r

#endif /* SRC_LIBSLIC3R_ALGORITHM_REGION_EXPANSION_HPP_ */
