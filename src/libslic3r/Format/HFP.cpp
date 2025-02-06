#include "HFP.hpp"
#include <fstream>
#include <iostream>

namespace Slic3r {

bool HFP::valid_hfp() const {
    return json_data.is_object();
}

bool HFP::load_hfp(const std::string& input_file, const DynamicPrintConfig* config) {
    std::ifstream file(input_file);
    
    if (!file.is_open()) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open HFP file: " << input_file;
        return false;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();  // Read the entire file into a buffer

    std::string file_content = buffer.str();  // Convert buffer to a string

    file.close();  // Close file after reading

    try {
        // Process the file_content based on its format
        std::cout << "HFP file content:\n" << file_content;
        
        // Example: If the file contains key-value pairs, process them here
        // Example: Parsing line by line
        std::istringstream stream(file_content);
        std::string line;
        while (std::getline(stream, line)) {
            BOOST_LOG_TRIVIAL(info) << "Processing line: " << line;
            // Process each line as needed
            
            
            
        }

    } catch (const std::exception& e) {
        BOOST_LOG_TRIVIAL(error) << "Error processing HFP file: " << e.what();
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Successfully loaded HFP file: " << input_file;
    return true;
}

bool HFP::apply_to_config() {
    if (!cfg) {
        BOOST_LOG_TRIVIAL(error) << "DynamicPrintConfig is null!";
        return false;
    }

    if (!valid_hfp()) {
        BOOST_LOG_TRIVIAL(error) << "Invalid HFP file structure!";
        return false;
    }

    // Define a mapping from HFP keys to Slic3r config keys
    std::unordered_map<std::string, std::string> key_mapping = {
        {"layer_height", "layer_height"},
        {"nozzle_diameter", "nozzle_diameter"},
        {"print_speed", "speed"},
        {"infill_density", "fill_density"},
        {"extruder_temperature", "temperature"},
        {"bed_temperature", "bed_temperature"}
    };

    for (const auto& [hfp_key, slic3r_key] : key_mapping) {
        if (json_data.contains(hfp_key)) {
            try {
                std::any value = json_data[hfp_key];

                if (value.type() == typeid(int)) {
                    cfg->set(slic3r_key, std::any_cast<int>(value));
                } else if (value.type() == typeid(double)) {
                    cfg->set(slic3r_key, std::any_cast<double>(value));
                } else if (value.type() == typeid(std::string)) {
                    cfg->set(slic3r_key, std::any_cast<std::string>(value));
                } else {
                    BOOST_LOG_TRIVIAL(error) << "Unsupported data type for key: " << hfp_key;
                }

                BOOST_LOG_TRIVIAL(info) << "Applied " << hfp_key << " -> " << slic3r_key;
            } catch (const std::exception& e) {
                BOOST_LOG_TRIVIAL(error) << "Failed to apply " << hfp_key << ": " << e.what();
            }
        }
    }

    return true;
}

} // namespace Slic3r
