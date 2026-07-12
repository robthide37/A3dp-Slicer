#ifndef SLIC3R_FORMAT_HFP_HPP_
#define SLIC3R_FORMAT_HFP_HPP_

#include "../PrintConfig.hpp"
#include "../GCode/ThumbnailData.hpp"
#include <functional>
#include <string>
#include <nlohmann/json.hpp>
#include <vector>
#include <any>
#include <utility>

#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include "../libslic3r.h"
#include "../Exception.hpp"
#include "../Model.hpp"
#include "../Utils.hpp"
#include "../LocalesUtils.hpp"
#include "../GCode.hpp"
#include "../Geometry.hpp"

namespace Slic3r {

class Model;
class DynamicPrintConfig;

class HFP
{
public:
    HFP();
    struct Filament
    {
        std::string Brand;
        std::string Color;
        std::string Name;
        bool Owned;
        double Transmissivity;
        std::string Type;
        std::string uuid;
    };

    bool valid_hfp() const;
    bool load_hfp(const std::string &input_file);

    // Getter functions
    const std::vector<Filament> &get_filament_set() const { return m_filament_set; }
    const std::vector<int> &get_slider_values() const { return m_slider_values; }
    double get_base_layer_height() const { return m_base_layer_height; }
    double get_layer_height() const { return m_layer_height; }
    std::string get_stl_path() { return m_stl_path; }

    // apply
    void update_config(DynamicPrintConfig &fff_print_config) const;
    void set_custom_gcode_z(Model &model) const;

private:
    nlohmann::json m_json_data;
    double m_base_layer_height;
    double m_layer_height;
    std::vector<Filament> m_filament_set;
    // always increase slider_values by +1
    std::vector<int> m_slider_values;
    std::string m_stl_path;
};

} // namespace Slic3r

#endif /* SLIC3R_FORMAT_HFP_HPP_ */
