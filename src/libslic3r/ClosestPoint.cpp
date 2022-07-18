#include "ClosestPoint.hpp"

size_t Slic3r::find_closest_in_sorted(const Point &p, const Points &pts)
{
    using namespace closestPoint;
    // check that input is really sorted
    assert(std::is_sorted(pts.begin(), pts.end(), sort_fnc));

    // check input
    if (pts.empty()) return std::numeric_limits<size_t>::max();
    if (pts.size() == 1) return 0;

    // closest point node in X
    Points::const_iterator it_x = std::upper_bound(pts.begin(), pts.end(), p.x(), upper_fnc);
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
        it_l = std::lower_bound(pts.begin(), it_x, p.x() - manhattan_dist, lower_fnc);
        for (auto it = it_x - 1; it > it_l; --it) {
            uint32_t diff_y = std::abs(it->y() - p.y());
            if (diff_y > manhattan_dist) continue;
            uint32_t diff_x   = std::abs(it->x() - p.x());
            uint32_t act_dist = diff_y + diff_x;
            if (manhattan_dist > act_dist) {
                manhattan_dist = act_dist;
                it_l = std::lower_bound(it_l, it_x, p.x() - manhattan_dist, lower_fnc);
            }
        }
    }

    // node for upper bound
    Points::const_iterator it_u;
    if (is_it_x_end) { 
        it_u = pts.end();
    } else {
        it_u = std::upper_bound(it_x, pts.end(), p.x() + manhattan_dist, upper_fnc);
        for (auto it = it_x + 1; it < it_u; ++it) { 
            uint32_t diff_y = std::abs(it->y() - p.y());
            if (diff_y > manhattan_dist) continue;
            uint32_t diff_x = std::abs(it->x() - p.x());
            uint32_t act_dist = diff_y + diff_x;
            if (manhattan_dist > act_dist) {
                // IMPROVE: calc euclid distance when e.g. (diff_Biggery < 2*diff_smaller)
                manhattan_dist = act_dist;
                it_u = std::upper_bound(it_x, it_u, p.x() + manhattan_dist, upper_fnc);
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

using namespace Slic3r;

bool closestPoint::sort_fnc(const Point &p1, const Point &p2)
{
    return p1.x() < p2.x();
}

bool closestPoint::upper_fnc(coord_t value, const Point &p)
{
    return value < p.x();
}

bool closestPoint::lower_fnc(const Point &p, coord_t value)
{
    return value > p.x();
}

uint32_t closestPoint::manhattan_size(const Point &p)
{
    return std::abs(p.x()) + abs(p.y());
}