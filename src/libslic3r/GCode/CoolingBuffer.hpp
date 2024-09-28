///|/ Copyright (c) Prusa Research 2017 - 2022 Vojtěch Bubník @bubnikv
///|/ Copyright (c) Slic3r 2016 Alessandro Ranellucci @alranel
///|/
///|/ ported from lib/Slic3r/GCode/CoolingBuffer.pm:
///|/ Copyright (c) Prusa Research 2016 - 2017 Vojtěch Bubník @bubnikv
///|/ Copyright (c) Slic3r 2013 - 2016 Alessandro Ranellucci @alranel
///|/ Copyright (c) 2016 Chow Loong Jin @hyperair
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_CoolingBuffer_hpp_
#define slic3r_CoolingBuffer_hpp_

#include "../libslic3r.h"
#include <map>
#include <string>

namespace Slic3r {

class GCodeGenerator;
class Layer;
struct PerExtruderAdjustments;

// A standalone G-code filter, to control cooling of the print.
// The G-code is processed per layer. Once a layer is collected, fan start / stop commands are edited
// and the print is modified to stretch over a minimum layer time.
//
// The simple it sounds, the actual implementation is significantly more complex.
// Namely, for a multi-extruder print, each material may require a different cooling logic.
// For example, some materials may not like to print too slowly, while with some materials 
// we may slow down significantly.
//
class CoolingBuffer {
public:
    CoolingBuffer(GCodeGenerator &gcodegen);
    void        reset(const Vec3d &position);
    void        set_current_extruder(unsigned int extruder_id) { m_current_extruder = extruder_id; }
    /// process the layer: check the time and apply fan / speed change
    /// append_time_only: if the layer is only support, then you can put this at true to not process the layer but just append its time to the next one.
    std::string process_layer(std::string &&gcode, size_t layer_id, bool flush, bool append_time_only = false);
    std::string process_layer(const std::string &gcode, size_t layer_id, bool flush)
        { return this->process_layer(std::string(gcode), layer_id, flush); }

private:
	CoolingBuffer& operator=(const CoolingBuffer&) = delete;
    std::vector<PerExtruderAdjustments> parse_layer_gcode(const std::string &gcode, std::array<float, 7> &current_pos) const;
    float       calculate_layer_slowdown(std::vector<PerExtruderAdjustments> &per_extruder_adjustments);
    // Apply slow down over G-code lines stored in per_extruder_adjustments, enable fan if needed.
    // Returns the adjusted G-code.
    std::string apply_layer_cooldown(const std::string &gcode, size_t layer_id, float layer_time, std::vector<PerExtruderAdjustments> &per_extruder_adjustments);

    // G-code snippet cached for the support layers preceding an object layer.
    std::string                 m_gcode;
    // Internal data.
    std::vector<char>           m_axis;
    enum AxisIdx : int {
        X = 0, Y, Z, E, F, I, J, K, R, Count
    };
    std::array<float, 7>        m_current_pos;
    // Current known fan speed or -1 if not known yet.
    int                         m_fan_speed;
    // Cached from GCodeWriter.
    // Printing extruder IDs, zero based.
    std::vector<unsigned int>   m_extruder_ids;
    // Highest of m_extruder_ids plus 1.
    uint16_t                    m_num_extruders { 0 };
    // Referencs GCodeGenerator::m_config, which is FullPrintConfig. While the PrintObjectConfig slice of FullPrintConfig is being modified,
    // the PrintConfig slice of FullPrintConfig is constant, thus no thread synchronization is required.
    const FullPrintConfig      &m_config;
    uint16_t                    m_current_extruder;
    uint16_t                    m_previous_extruder;

    //saved previous unslowed layer 
    std::map<size_t, float> saved_layer_time_support;
    std::map<size_t, float> saved_layer_time_object;


    // Old logic: proportional.
    bool                        m_cooling_logic_proportional = false;
};

}

#endif
