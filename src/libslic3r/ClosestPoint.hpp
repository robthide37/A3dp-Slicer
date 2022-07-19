#ifndef slic3r_ClosestPoint_hpp_
#define slic3r_ClosestPoint_hpp_

#include <vector>
#include <type_traits>

namespace Slic3r {

/// <summary>
/// Sort points and use find_closest_in_sorted
/// </summary>
/// <param name="p">Seach for closest index to this point</param>
/// <param name="pts">Search inside of thoose points</param>
/// <returns>Index of closest point from sorted_pts</returns>
template<class P>
// /*SFINAE*/ typename std::enable_if<std::is_member_function_pointer_v<decltype(&P::x)>>::type * = nullptr >
size_t find_closest(const P &p, const std::vector<P> &pts){
    // check input
    if (pts.empty()) return std::numeric_limits<size_t>::max();
    if (pts.size() == 1) return 0;

    using V = decltype(p.x());
    // extend P with order
    struct PP
    {
        size_t ord;
        const P* p;
        V x() const { return p->x(); }
        V y() const { return p->y(); }
    };
    std::vector<PP> pts_ord;
    pts_ord.reserve(pts.size());
    size_t ord = 0;
    for (const P &pp : pts) pts_ord.push_back({ord++, &pp});
    std::stable_sort(pts_ord.begin(), pts_ord.end(), closestPoint::sort_fnc<PP>);
    PP pp{0, &p};
    size_t closest_index = find_closest_in_sorted(pp, pts_ord);
    return pts_ord[closest_index].ord;
}

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
/// <param name="pts">Seach for closest index from thoose points</param>
/// <param name="sorted_pts">Seach for closest index to thoose points
/// Sorted points by X coordinate(by function closestPoint::sort_fnc)</param>
/// <returns>Indices to pts(first) and sorted_pts(second)</returns>
template<class P>
std::pair<size_t, size_t> find_closest_in_sorted(
    const std::vector<P> &pts, const std::vector<P> &sorted_pts);

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
template<class P> uint32_t manhattan_size(const P &p1, const P &p2)
{ return std::abs(p1.x()-p2.x()) + abs(p1.y()-p2.y()); }

template<class P>
typename std::vector<P>::const_iterator get_lower_bound(uint32_t& dist, const P &p, 
    typename std::vector<P>::const_iterator it_x, const std::vector<P>& pts)
{
    if(it_x == pts.begin()) return pts.begin();
    assert(it_x > pts.begin() && it_x < pts.end());

    using V = decltype(p.x());
    using It = std::vector<P>::const_iterator;
    It res = std::lower_bound(pts.begin(), it_x, p.x() - dist, lower_fnc<P,V>);
    for (auto it = it_x - 1; it > res; --it) {
        uint32_t diff_y = std::abs(it->y() - p.y());
        if (diff_y > dist) continue;
        uint32_t diff_x   = std::abs(it->x() - p.x());
        uint32_t act_dist = diff_y + diff_x;
        if (dist > act_dist) {
            dist = act_dist;
            res = std::lower_bound(res, it_x, p.x() - dist, lower_fnc<P,V>);
        }
    }
    return res;
}

template<class P>
typename std::vector<P>::const_iterator get_uppper_bound(uint32_t& dist, const P &p, 
    typename std::vector<P>::const_iterator it_x, const std::vector<P>& pts)
{
    assert(it_x >= pts.begin() && it_x < pts.end());
    
    using V = decltype(p.x());
    using It = std::vector<P>::const_iterator;
    It res = std::upper_bound(it_x, pts.end(), p.x() + dist, upper_fnc<P,V>);
    for (auto it = it_x + 1; it < res; ++it) { 
        uint32_t diff_y = std::abs(it->y() - p.y());
        if (diff_y > dist) continue;
        uint32_t diff_x = std::abs(it->x() - p.x());
        uint32_t act_dist = diff_y + diff_x;
        if (dist > act_dist) {
            // IMPROVE: calc euclid distance when e.g. (diff_Biggery < 2*diff_smaller)
            dist = act_dist;
            res = std::upper_bound(it_x, res, p.x() + dist, upper_fnc<P,V>);
        }
    }
    return res;
}

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
    using It = std::vector<P>::const_iterator;
    // closest point node in X
    It it_x = std::upper_bound(pts.begin(), pts.end(), p.x(), upper_fnc<P,V>);
    bool is_it_x_end = it_x == pts.end();
    // it_x can't pointing to end so change to last point
    if (is_it_x_end) --it_x;
    // manhatn distance to closest point
    uint32_t manhattan_dist = manhattan_size(*it_x, p);
    It it_l = get_lower_bound(manhattan_dist,p, it_x, pts);    

    // node for upper bound
    It it_u = (is_it_x_end) ? pts.end() : 
        get_uppper_bound(manhattan_dist, p, it_x, pts);
   
    // find closest by squer distance
    float dist_sq = std::numeric_limits<float>::max();
    size_t result = it_x - pts.begin();
    for (It it = it_l; it < it_u; ++it) {
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

template<class P>
std::pair<size_t, size_t> Slic3r::find_closest_in_sorted(
    const std::vector<P> &pts, const std::vector<P> &sorted_pts)
{
    using namespace closestPoint;
    // check that input is really sorted
    assert(std::is_sorted(sorted_pts.begin(), sorted_pts.end(), sort_fnc<P>));

    std::pair<size_t, size_t> res = {std::numeric_limits<size_t>::max(),
                                     std::numeric_limits<size_t>::max()};
    // check inputs
    if (pts.empty() || sorted_pts.empty()) return res;
    // Speed up for one element in vector
    if (pts.size() == 1)
        return {0, find_closest_in_sorted(pts.front(), sorted_pts)};
    if (sorted_pts.size() == 1) 
        return {find_closest(sorted_pts.front(), pts), 0};

    using V  = decltype(p.x());
    using It = std::vector<P>::const_iterator;
        
    uint32_t res_dist = std::numeric_limits<uint32_t>::max();
    float res_dist_sq = std::numeric_limits<float>::max();
    auto  set_res_dist_sq = [&res_dist_sq, &res_dist](float dist_sq) {
        res_dist_sq = dist_sq;
        float dist  = sqrt(dist_sq);
        // crop dist to unsigned int
        uint32_t dist_ui = static_cast<uint32_t>(std::ceil(dist));
        if (res_dist < dist_ui) res_dist = dist_ui;
    };
    // find limit for manhatn distance
    bool is_first = true;
    for (const P &p : pts) {
        if (is_first) {
            is_first = false;
            size_t index = find_closest_in_sorted(p, sorted_pts);
            res          = {0, index};
            const P &p2  = sorted_pts[index];
            V d_x = p.x()-p2.x();
            V d_y = p.y()-p2.y();
            set_res_dist_sq(d_x * d_x + d_y * d_y);
            continue;
        }
        // closest point node in X
        It it_x = std::upper_bound(pts.begin(), pts.end(), p.x(), upper_fnc<P, V>);        
        bool is_it_x_end = it_x == pts.end();
        // it_x can't pointing to end so change to last point
        if (is_it_x_end) --it_x;
        // manhatn distance to closest point
        uint32_t p_dist = manhattan_size(*it_x, p);
        
        It it_l = get_lower_bound(res_dist, p, it_x, pts);
                
        // node for upper bound
        It it_u = pts.end();
        if (!is_it_x_end) it_u = get_uppper_bound(res_dist, p, it_x, pts);

        // find closest by squer distance
        float dist_sq = res_dist_sq;
        for (It it = it_l; it < it_u; ++it) {
            uint32_t diff_y = std::abs(it->y() - p.y());
            if (diff_y > res_dist) continue;
            float diff_x = it->x() - p.x();
            // calculate square distance
            float d = (float) diff_y * diff_y + diff_x * diff_x;
            if (dist_sq > d) {
                dist_sq = d;
                size_t i1 = &p - &pts.begin();
                size_t i2 = it - pts.begin();
                res    = {i1, i2};
            }
        }

        // found closer?
        if (res_dist_sq > dist_sq)
            set_res_dist_sq(dist_sq);
    }
    return res;
}

#endif // slic3r_ClosestPoint_hpp_
