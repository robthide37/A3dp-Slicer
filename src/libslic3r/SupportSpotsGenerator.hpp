#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "libslic3r/Print.hpp"
#include <boost/log/trivial.hpp>

namespace Slic3r {

namespace SupportSpotsGenerator {

struct Params {
    Params(const std::vector<std::string> &filament_types) {
        if (filament_types.size() > 1) {
            BOOST_LOG_TRIVIAL(warning)
            << "SupportSpotsGenerator does not currently handle different materials properly, only first will be used";
        }
        if (filament_types.empty() || filament_types[0].empty()) {
            BOOST_LOG_TRIVIAL(error)
            << "SupportSpotsGenerator error: empty filament_type";
            filament_type = std::string("PLA");
        } else {
            filament_type = filament_types[0];
        }
    }

    // the algorithm should use the following units for all computations: distance [mm], mass [g], time [s], force [g*mm/s^2]
    const float bridge_distance = 12.0f; //mm
    const float bridge_distance_decrease_by_curvature_factor = 5.0f; // allowed bridge distance = bridge_distance / (1 + this factor * (curvature / PI) )
    const float overhang_angle_deg = 80.0f;
    const std::pair<float,float> malformation_angle_span_deg = std::pair<float, float> { 45.0f, 80.0f };

    const float min_distance_between_support_points = 3.0f; //mm
    const float support_points_interface_radius = 2.0f; // mm

    // NOTE: Currently disabled, does not work correctly due to inability of the algorithm to correctly detect islands at each layer
    const float supportable_volume_threshold = 0.0f; // mm^3

    std::string filament_type;
    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.
    const float max_acceleration = 9 * 1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY (NOTE: The max hit is received by the object in the jerk phase, so the usual machine limits are too low)
    const float filament_density = 1.25e-3f; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important
    const float material_yield_strength = 33.0f * 1e6f; // (g*mm/s^2)/mm^2; 33 MPa is yield strength of ABS, which has the lowest yield strength from common materials.
    const float standard_extruder_conflict_force = 20.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... );
    const float malformations_additive_conflict_extruder_force = 300.0f * gravity_constant; // for areas with possible high layered curled filaments

    // MPa * 1e^6 = (g*mm/s^2)/mm^2 = g/(mm*s^2); yield strength of the bed surface
    float get_bed_adhesion_yield_strength() const {
        if (filament_type == "PLA") {
            return 0.018f * 1e6f;
        } else if (filament_type == "PET" || filament_type == "PETG") {
            return 0.3f * 1e6f;
        } else { //PLA default value - defensive approach, PLA has quite low adhesion
            return 0.018f * 1e6f;
        }
    }

    //just return PLA adhesion value as value for supports
    float get_support_spots_adhesion_strength() const {
         return 0.018f * 1e6f; 
    }
};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float force, const Vec3f &direction);
    Vec3f position;
    float force;
    Vec3f direction;
};

struct Issues {
    std::vector<SupportPoint> support_points;
};

// std::vector<size_t> quick_search(const PrintObject *po, const Params &params);
Issues full_search(const PrintObject *po, const Params &params);

}

}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
