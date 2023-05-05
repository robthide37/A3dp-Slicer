#ifndef slic3r_SupportParameters_hpp_
#define slic3r_SupportParameters_hpp_

#include "../libslic3r.h"
#include "../Flow.hpp"

namespace Slic3r {

class PrintObject;
enum InfillPattern : int;

struct SupportParameters {
	SupportParameters(const PrintObject &object);

	// Flow at the 1st print layer.
	Flow 					first_layer_flow;
	// Flow at the support base (neither top, nor bottom interface).
	// Also flow at the raft base with the exception of raft interface and contact layers.
	Flow 					support_material_flow;
	// Flow at the top interface and contact layers.
	Flow 					support_material_interface_flow;
	// Flow at the bottom interfaces and contacts.
	Flow 					support_material_bottom_interface_flow;
	// Flow at raft inteface & contact layers.
	Flow    				raft_interface_flow;
	// Is merging of regions allowed? Could the interface & base support regions be printed with the same extruder?
	bool 					can_merge_support_regions;

    coordf_t 				support_layer_height_min;
//	coordf_t				support_layer_height_max;

	coordf_t				gap_xy;

    float    				base_angle;
    float    				interface_angle;

    // Density of the top / bottom interface and contact layers.
    coordf_t 				interface_density;
    // Density of the raft interface and contact layers.
    coordf_t 				raft_interface_density;
    // Density of the base support layers.
    coordf_t 				support_density;

    // Pattern of the sparse infill including sparse raft layers.
    InfillPattern           base_fill_pattern;
    // Pattern of the top / bottom interface and contact layers.
    InfillPattern           interface_fill_pattern;
    // Pattern of the raft interface and contact layers.
    InfillPattern           raft_interface_fill_pattern;
    // Pattern of the contact layers.
    InfillPattern 			contact_fill_pattern;
    // Shall the sparse (base) layers be printed with a single perimeter line (sheath) for robustness?
    bool                    with_sheath;

    float 					raft_angle_1st_layer;
    float 					raft_angle_base;
    float 					raft_angle_interface;

    // Produce a raft interface angle for a given SupportLayer::interface_id()
    float 					raft_interface_angle(size_t interface_id) const 
    	{ return this->raft_angle_interface + ((interface_id & 1) ? float(- M_PI / 4.) : float(+ M_PI / 4.)); }
};

} // namespace Slic3r

#endif /* slic3r_SupportParameters_hpp_ */
