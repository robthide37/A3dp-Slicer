#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"

namespace Slic3r {

namespace SupportSpotsGenerator {

struct Params {
    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.

    float bridge_distance = 10.0f; //mm
    float bridge_distance_decrease_by_curvature_factor = 5.0f; // >0 REQUIRED; allowed bridge distance = bridge_distance / (this factor * (curvature / PI) )

    float min_distance_between_support_points = 0.5f;

    // Adhesion computation : from experiment, PLA holds about 3g per mm^2 of base area (with reserve); So it can withstand about 3*gravity_constant force per mm^2
    float base_adhesion = 3.0f * gravity_constant; // adhesion per mm^2 of first layer
    float support_adhesion = 1.0f * gravity_constant; // adhesion per mm^2 of support interface layer

    float support_points_interface_radius = 0.5f; // mm

    float max_acceleration = 1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY
    float filament_density = 1.25f * 0.001f; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important

    float tolerable_extruder_conflict_force = 50.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... ); current value corresponds to weight of X grams
    float max_curled_conflict_extruder_force = 300.0f * gravity_constant; // for areas with possible high layered curled filaments, max force to account for;
};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float weight);
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
