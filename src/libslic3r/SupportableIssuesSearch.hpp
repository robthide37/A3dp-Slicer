#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"

namespace Slic3r {

namespace SupportableIssues {

struct Params {
    float bridge_distance = 10.0f;
    float limit_curvature = 0.15f; // used to detect curling issues

    float max_first_ex_perim_unsupported_distance_factor = 0.0f; // if external perim first, return tighter max allowed distance from previous layer extrusion
    float max_unsupported_distance_factor = 1.0f; // For internal perimeters, infill, bridges etc, allow gap of [extrusion width] size, these extrusions have usually something to stick to.
    float bridge_distance_decrease_by_curvature_factor = 5.0f; // allowed bridge distance = bridge_distance / ( 1 + this factor * (curvature / PI) )

    float base_adhesion = 100.0f; // adhesion per mm^2 of first layer; Force needed to remove the object from the bed, divided by the adhesion area, so in Pascal units N/mm^2
    float support_adhesion = 300.0f; // adhesion per mm^2 of support interface layer
    float support_points_inteface_area = 5.0f; // mm^2
    float max_acceleration = 1000.0f; // mm/s^2 ; max acceleration in XY
    float filament_density = 1.25f * 0.001f; // g/mm^3

};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float weight);
    explicit SupportPoint(const Vec3f &position);
    Vec3f position;
    float weight;
};

struct CurledFilament {
    CurledFilament(const Vec3f &position, float estimated_height);
    explicit CurledFilament(const Vec3f &position);
    Vec3f position;
    float estimated_height;
};

struct Issues {
    std::vector<SupportPoint> supports_nedded;
    std::vector<CurledFilament> curling_up;

    void add(const Issues &layer_issues);
    bool empty() const;
};

std::vector<size_t> quick_search(const PrintObject *po, const Params &params = Params { });
Issues full_search(const PrintObject *po, const Params &params = Params { });

}

}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
