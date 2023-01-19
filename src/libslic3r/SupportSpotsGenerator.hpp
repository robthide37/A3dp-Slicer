#ifndef SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_
#define SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_

#include "Layer.hpp"
#include "Line.hpp"
#include "PrintBase.hpp"
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
             BOOST_LOG_TRIVIAL(debug)
            << "SupportSpotsGenerator: applying filament type: " << filament_type;
        }
    }

    // the algorithm should use the following units for all computations: distance [mm], mass [g], time [s], force [g*mm/s^2]
    const float bridge_distance = 12.0f; //mm
    const std::pair<float,float> malformation_distance_factors = std::pair<float, float> { 0.4, 1.2 };
    const float max_curled_height_factor = 10.0f;

    const float min_distance_between_support_points = 3.0f; //mm
    const float support_points_interface_radius = 1.5f; // mm
    const float connections_min_considerable_area = 1.5f; //mm^2
    const float min_distance_to_allow_local_supports = 1.0f; //mm

    std::string filament_type;
    const float gravity_constant = 9806.65f; // mm/s^2; gravity acceleration on Earth's surface, algorithm assumes that printer is in upwards position.
    const float max_acceleration = 9 * 1000.0f; // mm/s^2 ; max acceleration of object (bed) in XY (NOTE: The max hit is received by the object in the jerk phase, so the usual machine limits are too low)
    const double filament_density = 1.25e-3f; // g/mm^3  ; Common filaments are very lightweight, so precise number is not that important
    const double material_yield_strength = 33.0f * 1e6f; // (g*mm/s^2)/mm^2; 33 MPa is yield strength of ABS, which has the lowest yield strength from common materials.
    const float standard_extruder_conflict_force = 20.0f * gravity_constant; // force that can occasionally push the model due to various factors (filament leaks, small curling, ... );
    const float malformations_additive_conflict_extruder_force = 100.0f * gravity_constant; // for areas with possible high layered curled filaments

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

// The support points are generated for two reasons:
// 1. Local extrusion support for extrusions that are printed in the air and would not 
//    withstand on their own (too long bridges, sharp turns in large overhang, concave bridge holes, etc.)
//    These points have negative force (-EPSILON) and Vec2f::Zero() direction
//    The algorithm still expects that these points will be supported and accounts for them in the global stability check
// 2. Global stability support points are generated at each spot, where the algorithm detects that extruding the current line 
//    may cause separation of the object part from the bed and/or its support spots or crack in the weak connection of the object parts
//    The generated point's direction is the estimated falling direction of the object part, and the force is equal to te difference 
//    between forces that destabilize the object (extruder conflicts with curled filament, weight if instable center of mass, bed movements etc)
//    and forces that stabilize the object (bed adhesion, other support spots adhesion, weight if stable center of mass)
//    Note that the force is only the difference - the amount needed to stabilize the object again.
struct SupportPoint {
    SupportPoint(const Vec3f &position, float force, float spot_radius, const Vec2f &direction);
    bool is_local_extrusion_support() const { return force < 0; }
    bool is_global_object_support() const { return !is_local_extrusion_support(); }

    //position is in unscaled coords. The z coordinate is aligned with the layers bottom_z coordiantes
    Vec3f position;
    // force that destabilizes the object to the point of falling/breaking. It is in g*mm/s^2 units
    // values gathered from large XL print: Min : 0 | Max : 18713800 | Average : 1361186 | Median : 329103
    // For reference 18713800 is weight of 1.8 Kg object, 329103 is weight of 0.03 Kg
    // The final printed object weight was approx 0.5 Kg
    float force;
    // Expected spot size. The support point strength is calculated from the area defined by this value. 
    // Currently equal to the support_points_interface_radius parameter above
    float spot_radius;
    // direction of the fall of the object (z part is neglected)
    Vec2f direction;
};

using SupportPoints = std::vector<SupportPoint>;

struct Malformations {
    std::vector<Lines> layers; //for each layer
};

// std::vector<size_t> quick_search(const PrintObject *po, const Params &params);
SupportPoints full_search(const PrintObject *po, const PrintTryCancel& cancel_func, const Params &params);

void estimate_supports_malformations(std::vector<SupportLayer*> &layers, float supports_flow_width, const Params &params);
void estimate_malformations(std::vector<Layer*> &layers, const Params &params);

} // namespace SupportSpotsGenerator
}

#endif /* SRC_LIBSLIC3R_SUPPORTABLEISSUESSEARCH_HPP_ */
