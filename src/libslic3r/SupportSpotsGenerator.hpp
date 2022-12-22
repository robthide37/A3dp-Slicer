#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "Layer.hpp"
#include "Line.hpp"
#include <boost/log/trivial.hpp>
#include <vector>

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
    const float bridge_distance_decrease_by_curvature_factor = 10.0f; // allowed bridge distance = bridge_distance / (1 + this factor * curvature )
    const std::pair<float,float> malformation_distance_factors = std::pair<float, float> { 0.4, 1.2 };
    const float max_malformation_factor = 10.0f;

    const float min_distance_between_support_points = 3.0f; //mm
    const float support_points_interface_radius = 1.5f; // mm
    const float connections_min_considerable_area = 1.5f; //mm^2
    const float min_distance_to_allow_local_supports = 2.0f; //mm

    std::string filament_type;
    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.
    const float max_acceleration = 9 * 1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY (NOTE: The max hit is received by the object in the jerk phase, so the usual machine limits are too low)
    const double filament_density = 1.25e-3f; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important
    const double material_yield_strength = 33.0f * 1e6f; // (g*mm/s^2)/mm^2; 33 MPa is yield strength of ABS, which has the lowest yield strength from common materials.
    const float standard_extruder_conflict_force = 20.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... );
    const float malformations_additive_conflict_extruder_force = 300.0f * gravity_constant; // for areas with possible high layered curled filaments

    // MPa * 1e^6 = (g*mm/s^2)/mm^2 = g/(mm*s^2); yield strength of the bed surface
    double get_bed_adhesion_yield_strength() const {
        if (filament_type == "PLA") {
            return 0.018 * 1e6;
        } else if (filament_type == "PET" || filament_type == "PETG") {
            return 0.3 * 1e6;
        } else { //PLA default value - defensive approach, PLA has quite low adhesion
            return 0.018 * 1e6;
        }
    }

    //just return PLA adhesion value as value for supports
    double get_support_spots_adhesion_strength() const {
         return 0.018f * 1e6; 
    }
};

struct SupportPoint {
    SupportPoint(const Vec3f &position, float force, float spot_radius, const Vec2f &direction);
    Vec3f position;
    float force;
    float spot_radius;
    Vec2f direction;
};

using SupportPoints = std::vector<SupportPoint>;

struct Malformations {
    std::vector<Lines> layers; //for each layer
};

// std::vector<size_t> quick_search(const PrintObject *po, const Params &params);
SupportPoints full_search(const PrintObject *po, const Params &params);

void estimate_supports_malformations(std::vector<SupportLayer*> &layers, float supports_flow_width, const Params &params);
void estimate_malformations(std::vector<Layer*> &layers, const Params &params);

} // namespace SupportSpotsGenerator
}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
