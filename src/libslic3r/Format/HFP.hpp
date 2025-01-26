#ifndef slic3r_Format_HFP_hpp_
#define slic3r_Format_HFP_hpp_

#include "../PrintConfig.hpp"
#include "../GCode/ThumbnailData.hpp"
#include <functional>
#include <string>
#include <nlohmann/json.hpp>
#include <vector.h>

namespace Slic3r {

class Model;
class DynamicPrintConfig;

class HFP {

    bool valid_hfp();
    bool load_hfp(std::string file_path);
    std::vector<std::pair<std::string, std::any>> get_hfp_values() const;



private:
    std::string file_path;
    nlohmann::json json_data;
    DynamicPrintConfig* cfg;


}


} // namespace Slic3r

#endif /* slic3r_Format_HFP_hpp_ */
