#ifndef slic3r_Format_AMF_hpp_
#define slic3r_Format_AMF_hpp_

namespace Slic3r {

class Model;
class DynamicPrintConfig;

// Load the content of an amf file into the given model and configuration.
extern bool load_amf(const char* path, DynamicPrintConfig* config, ConfigSubstitutionContext* config_substitutions, Model* model, bool check_version);

struct OptionStoreAmf {
    bool fullpath_sources = true;
    bool export_config = true;
    bool export_modifiers = true;
    OptionStoreAmf& set_fullpath_sources(bool use_fullpath_sources) { fullpath_sources = use_fullpath_sources; return *this; }
    OptionStoreAmf& set_export_config(bool use_export_config) { export_config = use_export_config; return *this; }
    OptionStoreAmf& set_export_modifiers(bool use_export_modifiers) { export_modifiers = use_export_modifiers; return *this; }
};

// Save the given model and the config data into an amf file.
// The model could be modified during the export process if meshes are not repaired or have no shared vertices
extern bool store_amf(std::string &path, Model *model, const DynamicPrintConfig *config, OptionStoreAmf options);

} // namespace Slic3r

#endif /* slic3r_Format_AMF_hpp_ */
