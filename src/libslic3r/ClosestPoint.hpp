#ifndef slic3r_ClosestPoint_hpp_
#define slic3r_ClosestPoint_hpp_

#include <vector>

namespace Slic3r {

/// <summary>
/// Sort points and use sweep line algo to find closest point in points
/// </summary>
/// <param name="p">Seach for closest index to this point</param>
/// <param name="pts">Search inside of thoose points</param>
/// <returns>Index of closest point from sorted_pts</returns>
// size_t find_closest(const Point &p, const Points& pts);

/// <summary>
/// Use a plane sweep algorithm to find closest point in sorted points
/// https://www.cs.mcgill.ca/~cs251/ClosestPair/ClosestPairPS.html
/// </summary>
/// <param name="p">Seach for closest index to this point</param>
/// <param name="sorted_pts">Sorted points by X coordinate</param>
/// <returns>Index of closest point from sorted_pts</returns>
template<class P>
size_t find_closest_in_sorted(const P &p, const std::vector<P> &sorted_pts);

/// <summary>
/// Use a plane sweep algorithm to find closest point from pts in sorted_pts
/// https://www.cs.mcgill.ca/~cs251/ClosestPair/ClosestPairPS.html
/// </summary>
/// <param name="pts">Seach for closest index to thoose points</param>
/// <param name="sorted_pts">Sorted points by X coordinate</param>
/// <returns>Index of closest point from sorted_pts</returns>
// size_t find_closest_in_sorted(const Point &pts, const Points &sorted_pts);

namespace closestPoint {
/// <summary>
/// Function used to sort points before searching for closest
/// </summary>
/// <param name="p1">First point</param>
/// <param name="p2">Second Point</param>
/// <returns>True when, p1.x < p2.x </returns>
template<class P> bool sort_fnc(const P &p1, const P &p2){ return p1.x() < p2.x(); }

/// <summary> Function used to find upper bound in sorted points. </summary>
template<class P, typename V> bool upper_fnc(V value, const P &p){ return value < p.x(); }
/// <summary> Function used to find lower bound in sorted points. </summary>
template<class P, typename V> bool lower_fnc(const P &p, V value){ return value > p.x(); }
/// <summary> Calc manhatn size of point. Mainly to explain meaning</summary>
template<class P> uint32_t manhattan_size(const P &p){ return std::abs(p.x()) + abs(p.y()); }
} // namespace closestPoint
} // namespace Slic3r

template<class P>
size_t Slic3r::find_closest_in_sorted(const P &p, const std::vector<P> &pts)
{
    using namespace closestPoint;
    // check that input is really sorted
    assert(std::is_sorted(pts.begin(), pts.end(), sort_fnc<P>));

    // check input
    if (pts.empty()) return std::numeric_limits<size_t>::max();
    if (pts.size() == 1) return 0;

    using V = decltype(p.x());

    // closest point node in X
    Points::const_iterator it_x = std::upper_bound(pts.begin(), pts.end(), p.x(), upper_fnc<P,V>);
    bool is_it_x_end = it_x == pts.end();
    // it_x can't pointing to end so change to last point
    if (is_it_x_end) --it_x;
    // manhatn distance to closest point
    uint32_t manhattan_dist = manhattan_size(*it_x - p);

    // node for lower bound
    Points::const_iterator it_l;
    if (it_x == pts.begin()) { 
        it_l = it_x;
    } else {
        it_l = std::lower_bound(pts.begin(), it_x, p.x() - manhattan_dist, lower_fnc<P,V>);
        for (auto it = it_x - 1; it > it_l; --it) {
            uint32_t diff_y = std::abs(it->y() - p.y());
            if (diff_y > manhattan_dist) continue;
            uint32_t diff_x   = std::abs(it->x() - p.x());
            uint32_t act_dist = diff_y + diff_x;
            if (manhattan_dist > act_dist) {
                manhattan_dist = act_dist;
                it_l = std::lower_bound(it_l, it_x, p.x() - manhattan_dist, lower_fnc<P,V>);
            }
        }
    }

    // node for upper bound
    Points::const_iterator it_u;
    if (is_it_x_end) { 
        it_u = pts.end();
    } else {
        it_u = std::upper_bound(it_x, pts.end(), p.x() + manhattan_dist, upper_fnc<P,V>);
        for (auto it = it_x + 1; it < it_u; ++it) { 
            uint32_t diff_y = std::abs(it->y() - p.y());
            if (diff_y > manhattan_dist) continue;
            uint32_t diff_x = std::abs(it->x() - p.x());
            uint32_t act_dist = diff_y + diff_x;
            if (manhattan_dist > act_dist) {
                // IMPROVE: calc euclid distance when e.g. (diff_Biggery < 2*diff_smaller)
                manhattan_dist = act_dist;
                it_u = std::upper_bound(it_x, it_u, p.x() + manhattan_dist, upper_fnc<P,V>);
            }
        }
    }

    // find closest by squer distance
    float dist_sq = std::numeric_limits<float>::max();
    size_t result = it_x - pts.begin();
    for (Points::const_iterator it = it_l; it < it_u; ++it) {
        uint32_t diff_y = std::abs(it->y() - p.y());
        if (diff_y > manhattan_dist) continue;
        float diff_x = it->x() - p.x();
        // calculate square distance
        float d = (float) diff_y * diff_y + diff_x * diff_x;
        if (dist_sq > d) {
            dist_sq = d;
            result  = it - pts.begin();
        }
    }
    return result;
}

#endif // slic3r_ClosestPoint_hpp_
