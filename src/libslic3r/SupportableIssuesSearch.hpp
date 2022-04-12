#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"

namespace Slic3r {

namespace SupportableIssues {

struct Params {
    float bridge_distance = 10.0f;
    float limit_curvature = 0.3f;

    float max_unsupported_distance_factor = 0.0f;
    float max_ex_perim_unsupported_distance_factor = 1.0f;
    float bridge_distance_decrease_by_curvature_factor = 5.0f;

    float perimeter_length_diff_tolerance = 8.0f;
    float centroid_offset_tolerance = 1.0f;
};


struct Issues {
    std::vector<Vec3f> supports_nedded;
    std::vector<Vec3f> curling_up;

    void add(const Issues &layer_issues);
    bool empty() const;
};

std::vector<size_t> quick_search(const PrintObject *po, const Params &params = Params { });
Issues full_search(const PrintObject *po, const Params &params = Params { });

}

}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
