#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"

namespace Slic3r {

namespace SupportSpotsGenerator {

struct Params {
    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.

    const float bridge_distance = 12.0f; //mm
    const float bridge_distance_decrease_by_curvature_factor = 5.0f; // allowed bridge distance = bridge_distance / (this factor * (curvature / PI) )

    const float min_distance_between_support_points = 3.0f;

    // Adhesion computation : from experiment, PLA holds about 3g per mm^2 of base area (with reserve); So it can withstand about 3*gravity_constant force per mm^2
    const float base_adhesion = 3.0f * gravity_constant; // adhesion per mm^2 of first layer
    const float support_adhesion = 1.0f * gravity_constant; // adhesion per mm^2 of support interface layer

    const float support_points_interface_radius = 0.6f; // mm

    const float max_acceleration = 9*1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY (NOTE: The max hit is received by the object in the jerk phase, so the usual machine limits are too low)
    const float filament_density = 1.25f * 0.001f; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important
    const float yield_strength = 33000.0f; // mN/mm^2;    33 MPa is yield strength of ABS, which has the lowest yield strength from common materials.
    const float standard_extruder_conflict_force = 1.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... ); current value corresponds to weight of X grams
    const float malformations_additive_conflict_extruder_force = 100.0f * gravity_constant; // for areas with possible high layered curled filaments

};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float force,const Vec3f& direction);
    Vec3f position;
    float force;
    Vec3f direction;
};

struct Issues {
    std::vector<SupportPoint> support_points;
};

std::vector<size_t> quick_search(const PrintObject *po, const Params &params = Params { });
Issues full_search(const PrintObject *po, const Params &params = Params { });

}

}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
