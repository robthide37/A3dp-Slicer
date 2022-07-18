#include "ClosestPoint.hpp"
#include "Point.hpp"

// Compile template for specific data type
template<>
size_t Slic3r::find_closest_in_sorted<Slic3r::Point>(const Slic3r::Point &p, const Slic3r::Points &pts);