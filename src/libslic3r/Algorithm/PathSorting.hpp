#ifndef SRC_LIBSLIC3R_PATH_SORTING_HPP_
#define SRC_LIBSLIC3R_PATH_SORTING_HPP_

#include "AABBTreeLines.hpp"
#include "ankerl/unordered_dense.h"
#include <algorithm>
#include <iterator>
#include <libslic3r/Point.hpp>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/ExPolygon.hpp>
#include <type_traits>
#include <unordered_set>

namespace Slic3r {
namespace Algorithm {

//Sorts the paths such that all paths between begin and last_seed are printed first, in some order. The rest of the paths is sorted
// such that the paths that are touching some of the already printed are printed first, sorted secondary by the distance to the last point of the last 
// printed path.
// begin, end, and last_seed are random access iterators. touch_limit_distance is used to check if the paths are touching - if any part of the path gets this close
// to the second, then they touch.
// convert_to_lines is a lambda that should accept the path as argument and return it as Lines vector, in correct order.
template<typename RandomAccessIterator, typename ToLines>
void sort_paths(RandomAccessIterator begin, RandomAccessIterator end, RandomAccessIterator last_seed, double touch_limit_distance, ToLines convert_to_lines)
{
    size_t paths_count = std::distance(begin, end);
    if (paths_count <= 1)
        return;

    std::vector<AABBTreeLines::LinesDistancer<Line>> distancers(paths_count);
    for (size_t path_idx = 0; path_idx < paths_count; path_idx++) {
        distancers[path_idx] = AABBTreeLines::LinesDistancer<Line>{convert_to_lines(*std::next(begin, path_idx))};
    }

    auto paths_touch = [touch_limit_distance](const AABBTreeLines::LinesDistancer<Line> &left,
                                              const AABBTreeLines::LinesDistancer<Line> &right) {
        for (const Line &l : left.get_lines()) {
            if (right.distance_from_lines<false>(l.a) < touch_limit_distance) {
                return true;
            }
        }
        if (right.distance_from_lines<false>(left.get_lines().back().b) < touch_limit_distance) {
            return true;
        }

        for (const Line &l : right.get_lines()) {
            if (left.distance_from_lines<false>(l.a) < touch_limit_distance) {
                return true;
            }
        }
        if (left.distance_from_lines<false>(right.get_lines().back().b) < touch_limit_distance) {
            return true;
        }
        return false;
    };

    std::vector<std::unordered_set<size_t>> dependencies(paths_count);
    for (size_t path_idx = 0; path_idx < paths_count; path_idx++) {
        for (size_t prev_path_idx = 0; prev_path_idx < path_idx; prev_path_idx++) {
            if (paths_touch(distancers[path_idx], distancers[prev_path_idx])) {
                dependencies[path_idx].insert(prev_path_idx);
                dependencies[prev_path_idx].insert(path_idx);
            }
        }
    }

    size_t index_of_last_fixed = std::distance(begin, last_seed);

    std::vector<bool> processed(paths_count, false);
    for (size_t path_idx = 0; path_idx <= index_of_last_fixed; path_idx++) {
        processed[path_idx] = true;
    }

    for (size_t i = index_of_last_fixed + 1; i < paths_count; i++) {
        bool change = false;
        for (size_t path_idx = index_of_last_fixed + 1; path_idx < paths_count; path_idx++) {
            if (processed[path_idx])
                continue;
            auto processed_dep = std::find_if(dependencies[path_idx].begin(), dependencies[path_idx].end(),
                                              [&](size_t dep) { return processed[dep]; });
            if (processed_dep != dependencies[path_idx].end()) {
                for (auto it = dependencies[path_idx].begin(); it != dependencies[path_idx].end();) {
                    if (!processed[*it]) {
                        dependencies[*it].insert(path_idx);
                        dependencies[path_idx].erase(it++);
                    } else {
                        ++it;
                    }
                }
                processed[path_idx] = true;
                change              = true;
            }
        }
        if (!change) {
            break;
        }
    }

    Point current_point = distancers.begin()->get_lines().begin()->a;

    size_t null_idx     = size_t(-1);
    size_t unsorted_idx = 0;
    size_t next_idx     = null_idx;
    bool   reverse      = false;
    while (true) {
        if (next_idx == null_idx) { // find next pidx to print
            double dist = std::numeric_limits<double>::max();
            for (size_t path_idx = 0; path_idx < paths_count; path_idx++) {
                if (!dependencies[path_idx].empty())
                    continue;
                const auto& lines = distancers[path_idx].get_lines();
                double      dist_a = (lines.front().a - current_point).cast<double>().squaredNorm();
                if (dist_a < dist) {
                    dist     = dist_a;
                    next_idx = path_idx;
                    reverse  = false;
                }
                double dist_b = (lines.back().b - current_point).cast<double>().squaredNorm();
                if (dist_b < dist) {
                    dist     = dist_b;
                    next_idx = path_idx;
                    reverse  = true;
                }
            }
            if (next_idx == null_idx) {
                       break;
            }
        } else {
            // we have valid next_idx, sort it, update dependencies, update current point and potentialy set new next_idx
            std::iter_swap(std::next(begin, unsorted_idx), std::next(begin, next_idx)); // next_path is now at sorted spot
            if (reverse) {
                std::next(begin, unsorted_idx)->reverse();
            }
            unsorted_idx++;

            assert(dependencies[next_idx].empty());
            dependencies[next_idx].insert(null_idx);
            current_point = reverse ? distancers[next_idx].get_lines().front().a : distancers[next_idx].get_lines().back().b;
            for (size_t path_idx = 0; path_idx < paths_count; path_idx++) {
                dependencies[path_idx].erase(next_idx);
            }
            double dist = std::numeric_limits<double>::max();
            next_idx    = null_idx;

            for (size_t path_idx = next_idx + 1; path_idx < paths_count; path_idx++) {
                if (!dependencies[path_idx].empty()) {
                    continue;
                }
                const auto &lines  = distancers[path_idx].get_lines();
                double      dist_a = (lines.front().a - current_point).cast<double>().squaredNorm();
                if (dist_a < dist) {
                    dist     = dist_a;
                    next_idx = path_idx;
                    reverse  = false;
                }
                double dist_b = (lines.back().b - current_point).cast<double>().squaredNorm();
                if (dist_b < dist) {
                    dist     = dist_b;
                    next_idx = path_idx;
                    reverse  = true;
                }
            }
            if (dist > scaled(5.0)) {
                next_idx = null_idx;
            }
        }
    }
}


}
}

#endif /*SRC_LIBSLIC3R_PATH_SORTING_HPP_*/