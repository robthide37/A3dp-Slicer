#include "HFP.hpp"
#include <fstream>
#include <iostream>

namespace Slic3r {


HFP::HFP() {
    m_base_layer_height = new float(0.0f);
    m_layer_height = new float(0.0f);
}


bool HFP::valid_hfp() const {
    return json_data.is_object();
}


HFP::~HFP() {
    if (m_base_layer_height) {
        delete m_base_layer_height;
        m_base_layer_height = nullptr;
    }
    if (m_layer_height) {
        delete m_layer_height;
        m_layer_height = nullptr;
    }
}


bool HFP::load_hfp(const std::string &input_file, const DynamicPrintConfig *config) {
    std::ifstream file(input_file);

    if (!file.is_open()) {
        BOOST_LOG_TRIVIAL(error) << "Failed to open HFP file: " << input_file;
        return false;
    }

    std::stringstream buffer;
    buffer << file.rdbuf();                  // Read the entire file into a buffer
    std::string file_content = buffer.str(); // Convert buffer to a string
    file.close();                            // Close file after reading

    BOOST_LOG_TRIVIAL(info) << "Processing HFP file: " << input_file;

    try {
        // Check if content is JSON format
        if (file_content.find("{") != std::string::npos) {
            BOOST_LOG_TRIVIAL(info) << "Detected JSON format.";
            json_data = nlohmann::json::parse(file_content);

            // Load Base Layer Height
            if (json_data.contains("base_layer_height")) {
                float base_layer_height = json_data.value("base_layer_height", 0.0f);
                m_base_layer_height = new float(base_layer_height);
                BOOST_LOG_TRIVIAL(info) << "Base Layer Height: " << *m_base_layer_height;
            }

            // Load Layer Height
            if (json_data.contains("layer_height")) {
                float layer_height = json_data.value("layer_height", 0.0f);
                m_layer_height = new float(layer_height);
                BOOST_LOG_TRIVIAL(info) << "Layer Height: " << *m_layer_height;
            }

            // Load Filament Set
            if (json_data.contains("filament_set") && json_data["filament_set"].is_array()) {
                BOOST_LOG_TRIVIAL(info) << "Loading filament set from JSON...";
                m_filament_set.clear();

                for (const auto &filament_json : json_data["filament_set"]) {
                    Filament filament;
                    filament.Brand = filament_json.value("Brand", "");
                    filament.Color = filament_json.value("Color", "");
                    filament.Name = filament_json.value("Name", "");
                    filament.Owned = filament_json.value("Owned", false);
                    filament.Transmissivity = filament_json.value("Transmissivity", 0.0);
                    filament.Type = filament_json.value("Type", "");
                    filament.uuid = filament_json.value("uuid", "");

                    m_filament_set.push_back(filament);
                    BOOST_LOG_TRIVIAL(info) << "Loaded filament: " << filament.Brand << " (" << filament.Name << ")";
                }
            } else {
                BOOST_LOG_TRIVIAL(warning) << "No 'filament_set' found in JSON.";
            }

            // Load Slider Values
            if (json_data.contains("slider_values") && json_data["slider_values"].is_array()) {
                BOOST_LOG_TRIVIAL(info) << "Loading slider values...";
                m_slider_values.clear();

                for (const auto &value : json_data["slider_values"]) {
                    int slider_value = value.get<int>() + 1; // Always increase by 1
                    m_slider_values.push_back(slider_value);
                }

                BOOST_LOG_TRIVIAL(info) << "Slider values loaded with +1 increment.";
            } else {
                BOOST_LOG_TRIVIAL(warning) << "No 'slider_values' found in JSON.";
            }
        }

        // Otherwise, process as key-value pairs (plain text format)
        else {
            BOOST_LOG_TRIVIAL(info) << "Detected Key-Value format.";
            std::istringstream stream(file_content);
            std::string line;

            while (std::getline(stream, line)) {
                BOOST_LOG_TRIVIAL(info) << "Processing line: " << line;

                size_t delimiter_pos = line.find(":");
                if (delimiter_pos != std::string::npos) {
                    std::string key = line.substr(0, delimiter_pos);
                    std::string value = line.substr(delimiter_pos + 1);

                    // Trim spaces
                    key.erase(0, key.find_first_not_of(" \t"));
                    key.erase(key.find_last_not_of(" \t") + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    value.erase(value.find_last_not_of(" \t") + 1);

                    // Process each key dynamically
                    if (key == "base_layer_height") {
                        m_base_layer_height = new float(std::stof(value));
                        BOOST_LOG_TRIVIAL(info) << "Base Layer Height: " << *m_base_layer_height;
                    } else if (key == "layer_height") {
                        m_layer_height = new float(std::stof(value));
                        BOOST_LOG_TRIVIAL(info) << "Layer Height: " << *m_layer_height;
                    } else if (key == "slider_values") {
                        std::istringstream value_stream(value);
                        int slider_value;
                        while (value_stream >> slider_value) {
                            m_slider_values.push_back(slider_value + 1); // Always increment by 1
                        }
                        BOOST_LOG_TRIVIAL(info) << "Loaded slider values with +1 increment.";
                    }
                    // Handle filaments dynamically (if stored as key-value in this format)
                    else if (key.find("filament_") == 0) {
                        Filament filament;
                        filament.Name = key;
                        filament.Brand = value;
                        m_filament_set.push_back(filament);
                        BOOST_LOG_TRIVIAL(info) << "Loaded filament: " << filament.Brand;
                    }
                }
            }
        }
    } catch (const std::exception &e) {
        BOOST_LOG_TRIVIAL(error) << "Error processing HFP file: " << e.what();
        return false;
    }

    return true;
}

const std::vector<HFP::Filament> &HFP::get_filament_set() const { return m_filament_set; }

const std::vector<int> &HFP::get_slider_values() const { return m_slider_values; }

const float *HFP::get_base_layer_height() const { return m_base_layer_height; }

const float *HFP::get_layer_height() const { return m_layer_height; }

bool HFP::apply_to_config() {
    if (!cfg) {
        BOOST_LOG_TRIVIAL(error) << "DynamicPrintConfig is null!";
        return false;
    }

    if (!valid_hfp()) {
        BOOST_LOG_TRIVIAL(error) << "Invalid HFP file structure!";
        return false;
    }



    return true;
}

} // namespace Slic3r
