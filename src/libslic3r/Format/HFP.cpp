#include "../libslic3r.h"
#include "../Model.hpp"
#include "../TriangleMesh.hpp"
#include "HFP.hpp"
#include "STL.hpp"

#include <string>


std::string KEYS = {
    "base_layer_height",
    "layer_height",
    "filament_set",
    "slider_values",
    "reverse_litho",
    "stl"
}

namespace Slic3r {


bool HFP::load_hfp(std::string& file_path) {
    std::ifstream input_file(file_path);
    if (!inputFile.is_open()) {
        return false;
    }

    try {
        input_file >> json_data
    } catch (nhlohmann::json::parse_error& err) {
        std::cerr << "parse error: " << e.what() << std::endl;
        return false
    }

    input_file.close();
    return true

}


}