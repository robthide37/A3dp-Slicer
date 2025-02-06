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
    HFP() = default;    
    
    bool valid_hfp() const;
    bool load_hfp(const std::string& input_file, const DynamicPrintConfig* config);
    std::vector<std::pair<std::string, std::any>> get_hfp_values() const;
    bool apply_to_config();  // NEW FUNCTION TO APPLY VALUES TO CONFIG

private:
    std::string file_path;
    nlohmann::json json_data;
    DynamicPrintConfig* cfg;
};

} // namespace Slic3r

#endif /* SLIC3R_FORMAT_HFP_HPP_ */
