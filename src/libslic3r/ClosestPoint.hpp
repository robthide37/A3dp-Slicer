#ifndef slic3r_ClosestPoint_hpp_
#define slic3r_ClosestPoint_hpp_

#include "Point.hpp"
#include "ExPolygon.hpp"

namespace Slic3r {

/// <summary>
/// Sort points and use sweep line algo to find closest point in points
/// </summary>
/// <param name="p">Seach for closest index to this point</param>
/// <param name="pts">Search inside of thoose points</param>
/// <returns>Index of closest point from sorted_pts</returns>
//size_t find_closest(const Point &p, const Points& pts);

/// <summary>
/// Use a plane sweep algorithm to find closest point in sorted points
/// https://www.cs.mcgill.ca/~cs251/ClosestPair/ClosestPairPS.html
/// </summary>
/// <param name="p">Seach for closest index to this point</param>
/// <param name="sorted_pts">Sorted points by X coordinate</param>
/// <returns>Index of closest point from sorted_pts</returns>
size_t find_closest_in_sorted(const Point &p, const Points &sorted_pts);

/// <summary>
/// Use a plane sweep algorithm to find closest point from pts in sorted_pts
/// https://www.cs.mcgill.ca/~cs251/ClosestPair/ClosestPairPS.html
/// </summary>
/// <param name="pts">Seach for closest index to thoose points</param>
/// <param name="sorted_pts">Sorted points by X coordinate</param>
/// <returns>Index of closest point from sorted_pts</returns>
//size_t find_closest_in_sorted(const Point &pts, const Points &sorted_pts);

namespace closestPoint {
/// <summary>
/// Function used to sort points before searching for closest
/// </summary>
/// <param name="p1">First point</param>
/// <param name="p2">Second Point</param>
/// <returns>True when, p1.x < p2.x </returns>
bool sort_fnc(const Point &p1, const Point &p2);

/// <summary> Function used to find upper bound in sorted points. </summary>
bool upper_fnc(coord_t value, const Point &p);
/// <summary> Function used to find lower bound in sorted points. </summary>
bool lower_fnc(const Point &p, coord_t value);
/// <summary> Calc manhatn size of point. Mainly to explain meaning</summary>
uint32_t manhattan_size(const Point &p);

} // namespace closestPoint
} // namespace Slic3r
#endif // slic3r_ClosestPoint_hpp_
