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

namespace Slic3r {

class Model;
class DynamicPrintConfig;

class HFP {

public:
    HFP();
    ~HFP();

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
    bool load_hfp(const std::string& input_file, const DynamicPrintConfig* config);
    bool apply_to_config();  // NEW FUNCTION TO APPLY VALUES TO CONFIG

    // Getter functions
    const std::vector<Filament> &get_filament_set() const;
    const std::vector<int> &get_slider_values() const;
    const float *get_base_layer_height() const;
    const float *get_layer_height() const;

private:
    std::string file_path;
    nlohmann::json json_data;
    DynamicPrintConfig* cfg;
    Model* m_model;
    const float* m_base_layer_height;
    const float* m_layer_height;
    std::vector<Filament> m_filament_set;
    // always increase slider_values by +1
    std::vector<int> m_slider_values;

};

} // namespace Slic3r

#endif /* SLIC3R_FORMAT_HFP_HPP_ */
