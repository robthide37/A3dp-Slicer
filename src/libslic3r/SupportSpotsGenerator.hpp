#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"

namespace Slic3r {

namespace SupportSpotsGenerator {

struct Params {
    // the algorithm should use the following units for all computations: distance [mm], mass [g], time [s], force [g*mm/s^2]
    const float bridge_distance = 12.0f; //mm
    const float bridge_distance_decrease_by_curvature_factor = 3.0f; // allowed bridge distance = bridge_distance / (this factor * (curvature / PI) )

    const float min_distance_between_support_points = 3.0f; //mm
    const float support_points_interface_radius = 0.6f; // mm

    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.
    const float max_acceleration = 9*1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY (NOTE: The max hit is received by the object in the jerk phase, so the usual machine limits are too low)
    const float filament_density = 1.25e-3f ; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important
    const float bed_adhesion_yield_strength = 0.128f * 1e6f; //MPa * 1e^6 = (g*mm/s^2)/mm^2 = g/(mm*s^2); yield strength of the bed surface
    const float material_yield_strength = 15.0f * 1e6f; // (g*mm/s^2)/mm^2; 33 MPa is yield strength of ABS, which has the lowest yield strength from common materials.
    const float standard_extruder_conflict_force = 20.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... );
    const float malformations_additive_conflict_extruder_force = 300.0f * gravity_constant; // for areas with possible high layered curled filaments
};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float force, const Vec3f& direction);
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
