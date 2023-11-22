#include "ScriptExecutor.hpp"
#include "GUI_App.hpp"
#include "Plater.hpp"
#include "Tab.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Print.hpp"

#include <string>

#include <angelscript/source/as_config.h>
#include <angelscript/add_on/autowrapper/aswrappedcall.h>
#include <angelscript/add_on/scriptarray/scriptarray.h>
#include <angelscript/add_on/scriptbuilder/scriptbuilder.h>
#include <angelscript/add_on/scriptstdstring/scriptstdstring.h>
#include <angelscript/add_on/scriptmath/scriptmath.h>

using namespace gw;

namespace Slic3r {  namespace GUI {

void as_message_callback(const AngelScript::asSMessageInfo* msg, void* param)
{
    const char* type = "ERR ";
    if (msg->type == AngelScript::asMSGTYPE_WARNING)
        type = "WARN";
    else if (msg->type == AngelScript::asMSGTYPE_INFORMATION)
        type = "INFO";

    printf("%s (%d, %d) : %s : %s\n", msg->section, msg->row, msg->col, type, msg->message);
}

//FIXME put script methods in the ScriptContainer class so there isn't this dangerous global var here.
ScriptContainer* current_script;
void as_print(std::string& str)
{
    std::cout << str;
}
void as_print_float(float f)
{
    std::cout << f;
}
//void as_register_key(std::string& key) {
//    if (watching_keys != nullptr)
//        watching_keys->push_back(key);
//}
std::pair<const PresetCollection*, const ConfigOption*> get_coll(const std::string& str) {
    const PresetCollection* coll = nullptr;
    const ConfigOption* opt = nullptr;
    if(opt == nullptr && (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0) {
        coll = &current_script->tab()->m_preset_bundle->fff_prints;
        opt = coll->get_edited_preset().config.option(str);
        if (opt == nullptr) {
            coll = &current_script->tab()->m_preset_bundle->filaments;
            opt = coll->get_edited_preset().config.option(str);
        }
    }
    if (opt == nullptr && (current_script->tab()->get_printer_technology() & PrinterTechnology::ptSLA) != 0) {
        coll = &current_script->tab()->m_preset_bundle->sla_prints;
        opt = coll->get_edited_preset().config.option(str);
        if (opt == nullptr) {
            coll = &current_script->tab()->m_preset_bundle->sla_materials;
            opt = coll->get_edited_preset().config.option(str);
        }
    }
    if (opt == nullptr) {
        coll = &current_script->tab()->m_preset_bundle->printers;
        opt = coll->get_edited_preset().config.option(str);
    }
    return std::pair<const PresetCollection*, const ConfigOption*>{ coll,  opt };
}
//PresetCollection* get_coll(int preset_type) {
//    if (preset_type <= 0)
//        return current_script->tech() == (PrinterTechnology::ptFFF)
//        ? &current_script->tab()->m_preset_bundle->fff_prints
//        : &current_script->tab()->m_preset_bundle->sla_prints;
//    else if (preset_type == 1)
//        return current_script->tech() == (PrinterTechnology::ptFFF)
//        ? &current_script->tab()->m_preset_bundle->filaments
//        : &current_script->tab()->m_preset_bundle->sla_materials;
//    else return &current_script->tab()->m_preset_bundle->printers;
//}
bool as_get_bool(std::string& key)
{
    const ConfigOption* opt = get_coll(key).second;
    if (opt == nullptr || opt->type() != ConfigOptionType::coBool)
        throw NoDefinitionException("error, can't find bool option " + key);
    if (opt->is_vector()) {
        const ConfigOptionVectorBase* vector = static_cast<const ConfigOptionVectorBase*>(opt);
        return vector->getFloat(0) != 0;
    } else {
        return opt->getBool();
    }
}
void as_set_bool(std::string& key, bool b)
{
    if (!current_script->can_set()) return;
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw NoDefinitionException("error, can't find bool option " + key);
    DynamicPrintConfig& conf = current_script->to_update()[result.first->type()];
    if (result.second->type() == ConfigOptionType::coBool) {
        conf.set_key_value(key, new ConfigOptionBool(b));
    } else if (result.second->type() == ConfigOptionType::coBools) {
        ConfigOptionBools* new_val = static_cast<ConfigOptionBools*>(result.second->clone());
        new_val->set_at(b, 0);
        conf.set_key_value(key, new_val);
    }
}
int32_t as_get_int(std::string& key)
{
    const ConfigOption* opt = get_coll(key).second;
    if (opt == nullptr || (opt->type() != ConfigOptionType::coInt && opt->type() != ConfigOptionType::coEnum))
        throw NoDefinitionException("error, can't find int option " + key);
    if (opt->is_vector()) {
        const ConfigOptionVectorBase* vector = static_cast<const ConfigOptionVectorBase*>(opt);
        return (int32_t)vector->getFloat(0);
    } else {
        return (int32_t)(opt->getInt());
    }
}
void as_set_int(std::string& key, int val)
{
    if (!current_script->can_set()) return;
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw NoDefinitionException("error, can't find int option " + key);
    DynamicPrintConfig& conf = current_script->to_update()[result.first->type()];
    if (result.second->type() == ConfigOptionType::coInt) {
        conf.set_key_value(key, new ConfigOptionInt(val));
    } else if (result.second->type() == ConfigOptionType::coInts) {
        ConfigOptionInts* new_val = static_cast<ConfigOptionInts*>(result.second->clone());
        new_val->set_at(val, 0);
        conf.set_key_value(key, new_val);
    } else if (result.second->type() == ConfigOptionType::coEnum) {
        //const ConfigOptionDef* def = result.first->get_edited_preset().config.get_option_def(key);
        //if (!def->enum_values.empty()) {
        //    std::string key = "";
        //    for (const auto& entry : *def->enum_keys_map) {
        //        if (entry.second == val) {
        //            key = entry.first;
        //        }
        //    }
        //    int32_t value = -1;
        //    for (int i = 0; i < def->enum_values.size(); i++) {
        //        if (def->enum_values[i] == key) {
        //            value = i;
        //            break;
        //        }
        //    }
        //    if (value >= 0 && value < def->enum_values.size()) {
                ConfigOption* copy = result.second->clone();
                copy->setInt(val);
                conf.set_key_value(key, copy);
        //        return;
        //    }
        //}
        //BOOST_LOG_TRIVIAL(error) << "Error, can't access enum '" << key << "'";
    }
}
float as_get_float(std::string& key)
{
    const ConfigOption* opt = get_coll(key).second;
    if (opt == nullptr) //TODO check if  float, etc..
        throw NoDefinitionException("error, can't find float option " + key);
    if (opt->is_vector()) {
        const ConfigOptionVectorBase* vector = static_cast<const ConfigOptionVectorBase*>(opt);
        return (float)vector->getFloat(0);
    } else {
        return (float)(opt->getFloat());
    }
}

double round(float f) {
    std::stringstream ss;
    double dbl_val;
    ss << f;
    ss >> dbl_val;
    return dbl_val;
}

void as_set_float(std::string& key, float f_val)
{
    if (!current_script->can_set()) return;
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw NoDefinitionException("error, can't find float option " + key);

    DynamicPrintConfig& conf = current_script->to_update()[result.first->type()];
    if (result.second->type() == ConfigOptionType::coFloat) {
        double old_value = result.second->getFloat();
        double new_val = round(f_val);
        // only update if difference is significant
        if (std::abs(old_value - new_val) / std::abs(old_value) < 0.0000001)
            new_val = old_value; // don't return int these check, as it can escpae a refresh of the scripted widget
        conf.set_key_value(key, new ConfigOptionFloat(new_val));
    } else if (result.second->type() == ConfigOptionType::coFloats) {
        ConfigOptionFloats* new_opt = static_cast<ConfigOptionFloats*>(result.second->clone());
        double new_val = round(f_val);
        if (!new_opt->values.empty()) {
            // only update if difference is significant
            double old_value = new_opt->values.front();
            if (std::abs(old_value - new_val) / std::abs(old_value) < 0.0000001)
                new_val = old_value;
        }
        new_opt->set_at(new_val, 0);
        conf.set_key_value(key, new_opt);
    } else if (result.second->type() == ConfigOptionType::coPercent) {
        double percent_f = floor(f_val * 100000. + 0.5) / 1000.;
        // only update if difference is significant
        double old_value = result.second->getFloat();
        if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
            percent_f = old_value;
        conf.set_key_value(key, new ConfigOptionPercent(percent_f));
    } else if (result.second->type() == ConfigOptionType::coPercents) {
        ConfigOptionPercents* new_opt = static_cast<ConfigOptionPercents*>(result.second->clone());
        double percent_f = floor(f_val * 100000. + 0.5) / 1000.;
        if (!new_opt->values.empty()) {
            // only update if difference is significant
            double old_value = new_opt->values.front();
            if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
                percent_f = old_value;
        }
        new_opt->set_at(percent_f, 0);
        conf.set_key_value(key, new_opt);
    } else if (result.second->type() == ConfigOptionType::coFloatOrPercent) {
        double new_val = round(f_val);
        if (!static_cast<const ConfigOptionFloatOrPercent*>(result.second)->percent) {
            // only update if difference is significant
            double old_value = result.second->getFloat();
            if (std::abs(old_value - new_val) / std::abs(old_value) < 0.0000001)
                new_val = old_value;
        }
        conf.set_key_value(key, new ConfigOptionFloatOrPercent(new_val, false));
    } else if (result.second->type() == ConfigOptionType::coFloatsOrPercents) {
        ConfigOptionFloatsOrPercents* new_opt = static_cast<ConfigOptionFloatsOrPercents*>(result.second->clone());
        double new_val = round(f_val);
        if (!new_opt->values.empty() && !new_opt->values.front().percent) {
            // only update if difference is significant
            double old_value = new_opt->values.front().value;
            if (std::abs(old_value - new_val) / std::abs(old_value) < 0.0000001)
                new_val = old_value;
        }
        new_opt->set_at(FloatOrPercent{ new_val, false}, 0);
        conf.set_key_value(key, new_opt);
    }
}
bool as_is_percent(std::string& key)
{
    const ConfigOption* opt = get_coll(key).second;
    if (opt == nullptr)
        throw NoDefinitionException("error, can't find percent option " + key);
    return (opt->type() == ConfigOptionType::coPercent) || (opt->type() == ConfigOptionType::coPercents) 
        || (opt->type() == ConfigOptionType::coFloatOrPercent && ((ConfigOptionFloatOrPercent*)opt)->percent) 
        || (opt->type() == ConfigOptionType::coFloatsOrPercents && ((ConfigOptionFloatsOrPercents*)opt)->get_at(0).percent);
}
void as_set_percent(std::string& key, float f_val)
{
    if (!current_script->can_set()) return;
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw NoDefinitionException("error, can't find percent option " + key);
    double percent_f = floor(f_val * 1000. + 0.5) / 1000.;
    DynamicPrintConfig& conf = current_script->to_update()[result.first->type()];
    if (result.second->type() == ConfigOptionType::coFloat) {
        // only update if difference is significant
        double old_value = result.second->getFloat() * 100;
        if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
            percent_f = old_value; // don't return int these check, as it can escpae a refresh of the scripted widget
        conf.set_key_value(key, new ConfigOptionFloat(percent_f / 100.));
    } else if (result.second->type() == ConfigOptionType::coFloats) {
        ConfigOptionFloats* new_opt = static_cast<ConfigOptionFloats*>(result.second->clone());
        if (!new_opt->values.empty()) {
            // only update if difference is significant
            double old_value = new_opt->values.front() * 100;
            if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
                percent_f = old_value;
        }
        new_opt->set_at(percent_f / 100., 0);
        conf.set_key_value(key, new_opt);
    } else if (result.second->type() == ConfigOptionType::coPercent) {
        // only update if difference is significant
        double old_value = get_coll(key).second->getFloat();
        if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
            percent_f = old_value;
        conf.set_key_value(key, new ConfigOptionPercent(percent_f));
    } else if (result.second->type() == ConfigOptionType::coPercents) {
        ConfigOptionPercents* new_opt = static_cast<ConfigOptionPercents*>(result.second->clone());
        if (!new_opt->values.empty()) {
            // only update if difference is significant
            double old_value = new_opt->values.front();
            if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
                percent_f = old_value;
        }
        new_opt->set_at(percent_f, 0);
        conf.set_key_value(key, new_opt);
    } else if (result.second->type() == ConfigOptionType::coFloatOrPercent) {
        if (static_cast<const ConfigOptionFloatOrPercent*>(result.second)->percent) {
            // only update if difference is significant
            double old_value = result.second->getFloat();
            if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
                percent_f = old_value;
        }
        conf.set_key_value(key, new ConfigOptionFloatOrPercent(percent_f, true));
    } else if (result.second->type() == ConfigOptionType::coFloatsOrPercents) {
        ConfigOptionFloatsOrPercents* new_opt = static_cast<ConfigOptionFloatsOrPercents*>(result.second->clone());
        if (!new_opt->values.empty() && new_opt->values.front().percent) {
            // only update if difference is significant
            double old_value = new_opt->values.front().value;
            if (std::abs(old_value - percent_f) / std::abs(old_value) < 0.0000001)
                percent_f = old_value;
        }
        new_opt->set_at(FloatOrPercent{ percent_f, true }, 0);
        conf.set_key_value(key, new_opt);
    }
}
void as_get_string(std::string& key, std::string& val)
{
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    const ConfigOption* opt = result.second;
    if (opt == nullptr) //TODO check if  float, etc..
        throw NoDefinitionException("error, can't find string option " + key);
    if (opt->type() == ConfigOptionType::coString) {
        val = ((ConfigOptionString*)opt)->value;
    } else if (opt->type() == ConfigOptionType::coStrings) {
        val = ((ConfigOptionStrings*)opt)->get_at(0);
    } else if (opt->type() == ConfigOptionType::coEnum) {
        val = opt->serialize();
    }
}
void as_set_string(std::string& key, std::string& val)
{
    if (!current_script->can_set()) return;
    std::pair<const PresetCollection*, const ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw NoDefinitionException("error, can't find string option " + key);
    DynamicPrintConfig& conf = current_script->to_update()[result.first->type()];
    if (result.second->type() == ConfigOptionType::coString) {
        conf.set_key_value(key, new ConfigOptionString(val));
    } else if (result.second->type() == ConfigOptionType::coStrings) {
        ConfigOptionStrings* new_val = (ConfigOptionStrings*)result.second->clone();
        new_val->set_at(val, 0);
        conf.set_key_value(key, new_val);
    } else if (result.second->type() == ConfigOptionType::coEnum) {
        const ConfigOptionDef* def = result.first->get_edited_preset().config.get_option_def(key);
        int idx = 0;
        for (; idx < def->enum_values.size() && def->enum_values[idx] != val; idx++) {}
        if (idx >= 0 && idx < def->enum_values.size()) {
            ConfigOption* copy = result.second->clone();
            copy->setInt(idx);
            conf.set_key_value(key, copy);
        }
    }
}

/////// custom vars ////////

std::string get_custom_var_option(int preset_type) {
    if (preset_type <= 0)
        return (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
        ? current_script->tab()->m_preset_bundle->fff_prints.get_edited_preset().config.opt_string("print_custom_variables")
        : current_script->tab()->m_preset_bundle->sla_prints.get_edited_preset().config.opt_string("print_custom_variables");
    else if (preset_type == 1) {
        return (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
            ? current_script->tab()->m_preset_bundle->filaments.get_edited_preset().config.opt_string("filament_custom_variables", (unsigned int)(0))
            : current_script->tab()->m_preset_bundle->sla_materials.get_edited_preset().config.opt_string("filament_custom_variables", (unsigned int)(0));
    } else return current_script->tab()->m_preset_bundle->printers.get_edited_preset().config.opt_string("printer_custom_variables");
}
std::string get_custom_value(std::string custom_var_field, const std::string& opt_key) {
    if (custom_var_field.find(opt_key) != std::string::npos) {
        boost::erase_all(custom_var_field, "\r");
        std::vector<std::string> lines;
        boost::algorithm::split(lines, custom_var_field, boost::is_any_of("\n"));
        for (const std::string& line : lines) {
            size_t equal_pos = line.find_first_of('=');
            if (equal_pos != std::string::npos) {
                std::string name = line.substr(0, equal_pos);
                std::string value = line.substr(equal_pos + 1);
                boost::algorithm::trim(name);
                if (name == opt_key) {
                    boost::algorithm::trim(value);
                    return value;
                }
            }
        }
    }
    return "";
}
void set_custom_option(int preset_type, std::string new_value) {
    if (!current_script->can_set()) return;
    if (preset_type <= 0) {
        PresetCollection& coll = (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
            ? current_script->tab()->m_preset_bundle->fff_prints
            : current_script->tab()->m_preset_bundle->sla_prints;
        DynamicPrintConfig& conf = current_script->to_update()[coll.type()];
        conf.set_key_value("print_custom_variables", new ConfigOptionString(new_value));
    } else if (preset_type == 1) {
        const PresetCollection& coll = (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
            ? current_script->tab()->m_preset_bundle->filaments
            : current_script->tab()->m_preset_bundle->sla_materials;
        DynamicPrintConfig& conf = current_script->to_update()[coll.type()];
        conf.set_key_value("filament_custom_variables", new ConfigOptionString(new_value));
    } else {
        DynamicPrintConfig& conf = current_script->to_update()[Preset::Type::TYPE_PRINTER];
        conf.set_key_value("printer_custom_variables", new ConfigOptionString(new_value));
    }
}
void set_custom_value(std::string& custom_var_field, const std::string& opt_key, const std::string& new_value) {
    boost::erase_all(custom_var_field, "\r");
    bool found = false;
    //iterate onlines,until you find the good one
    std::vector<std::string> lines;
    boost::algorithm::split(lines, custom_var_field, boost::is_any_of("\n"));
    for (auto line = lines.begin(); line != lines.end(); ++line) {
        size_t equal_pos = line->find_first_of('=');
        if (equal_pos != std::string::npos) {
            std::string name = line->substr(0, equal_pos);
            std::string value = line->substr(equal_pos + 1);
            boost::algorithm::trim(name);
            if (name == opt_key) {
                if (new_value == "") {
                    //delete
                    lines.erase(line);
                } else {
                    *line = name + " = " + new_value;
                }
                found = true;
                break;
            }
        }
    }
    if (!found) {
        //put it at the end
        lines.push_back(opt_key + " = " + new_value);
    }
    //reconstruct from lines
    custom_var_field = "";
    for (std::string& line : lines)
        if(!line.empty())
            custom_var_field += line + "\n";
}
bool as_get_custom_bool(int preset, std::string& key, bool& result)
{
    std::string serialized_value = get_custom_value(get_custom_var_option(preset), key);
    if (serialized_value.empty())
        return false;

    if (serialized_value == "true") {
        result = true;
        return true;
    } else if (serialized_value == "false") {
        result = false;
        return true;
    }
    return false;
}
void as_set_custom_bool(int preset, std::string& key, bool b)
{
    if (!current_script->can_set()) return;
    std::string serialized_vars = get_custom_var_option(preset);
    set_custom_value(serialized_vars, key, b ? "true" : "false");
    set_custom_option(preset, serialized_vars);
}
bool as_get_custom_int(int preset, std::string& key, int32_t& result)
{
    std::string serialized_value = get_custom_value(get_custom_var_option(preset), key);
    if (serialized_value.empty())
        return false;

    try {
        result = boost::lexical_cast<int32_t>(serialized_value);
        return true;
    }
    catch (boost::bad_lexical_cast&) {
    }
    return false;
}
void as_set_custom_int(int preset, std::string& key, int32_t val)
{
    if (!current_script->can_set()) return;
    std::string serialized_vars = get_custom_var_option(preset);
    set_custom_value(serialized_vars, key, std::to_string(val));
    set_custom_option(preset, serialized_vars);
}
bool as_get_custom_float(int preset, std::string& key, float& result)
{
    std::string serialized_value = get_custom_value(get_custom_var_option(preset), key);
    if (serialized_value.empty())
        return false;

    try {
        result = boost::lexical_cast<float>(serialized_value);
        return true;
    }
    catch (boost::bad_lexical_cast&) {
    }
    return false;
}
void as_set_custom_float(int preset, std::string& key, float f)
{
    if (!current_script->can_set()) return;
    std::string serialized_vars = get_custom_var_option(preset);
    set_custom_value(serialized_vars, key, std::to_string(f)); //FIXME: 0.0 format, with only 9 significant numbers
    set_custom_option(preset, serialized_vars);
}
bool as_get_custom_string(int preset, std::string& key, std::string& result)
{
    std::string serialized_value = get_custom_value(get_custom_var_option(preset), key);
    if (serialized_value.empty())
        return false;

    result = serialized_value;
    return true;
}
void as_set_custom_string(int preset, std::string& key, std::string& val)
{
    if (!current_script->can_set()) return;
    std::string serialized_vars = get_custom_var_option(preset);
    set_custom_value(serialized_vars, key, val);
    set_custom_option(preset, serialized_vars);
}

////// others //////

class ConfigAdapter : public ConfigBase
{
public:
    const ConfigBase* real_storage;
    ConfigAdapter(const ConfigBase* conf) : real_storage(conf) { this->parent = nullptr; }
    ConfigAdapter(const ConfigBase* conf, const ConfigBase* conf_parent) : real_storage(conf) { this->parent = conf_parent;  }
    virtual ~ConfigAdapter() { if (this->parent != nullptr) delete parent; }

    virtual const ConfigDef* def() const override { return real_storage->def(); }
    virtual ConfigOption* optptr(const t_config_option_key& opt_key, bool create = false) { return nullptr; }
    virtual t_config_option_keys keys() const override { return real_storage->keys(); };

    const ConfigOption* optptr(const t_config_option_key& opt_key) const override
    {
        const ConfigOption* opt = real_storage->optptr(opt_key);
            //if not find, try with the parent config.
        if (opt == nullptr && parent != nullptr)
            opt = parent->optptr(opt_key);
        return opt;
    }


};

float as_get_computed_float(std::string& key)
{
    if ((current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0) {
        ConfigAdapter fullconfig(
            &current_script->tab()->m_preset_bundle->fff_prints.get_edited_preset().config, 
            new ConfigAdapter(
                &current_script->tab()->m_preset_bundle->filaments.get_edited_preset().config, 
                new ConfigAdapter(&current_script->tab()->m_preset_bundle->printers.get_edited_preset().config)));
        try {
            return (float)fullconfig.get_computed_value(key, 0);
            //return (float)wxGetApp().plater()->fff_print().full_print_config().get_computed_value(key, 0);
        }
        catch (Exception e) {
            if(wxGetApp().initialized())
                BOOST_LOG_TRIVIAL(error) << "Error, can't compute fff option '" << key << "'";
        }

    } else {
        ConfigAdapter fullconfig(
            &current_script->tab()->m_preset_bundle->sla_prints.get_edited_preset().config,
            new ConfigAdapter(
                &current_script->tab()->m_preset_bundle->sla_materials.get_edited_preset().config,
                new ConfigAdapter(&current_script->tab()->m_preset_bundle->printers.get_edited_preset().config)));
        try {
            return (float)fullconfig.get_computed_value(key, 0);
            //return (float)wxGetApp().plater()->fff_print().full_print_config().get_computed_value(key, 0);
        }
        catch (Exception e) {
            if (wxGetApp().initialized())
                BOOST_LOG_TRIVIAL(error) << "Error, can't compute sla option '" << key << "'";
        }

    }
    return 0;
}

void as_ask_for_refresh()
{
    current_script->request_refresh();
}

//function to reset a field
void as_back_initial_value(std::string& key) {
    current_script->add_to_reset(key);
}
void as_back_custom_initial_value(int preset_type, std::string& key) {
    if (!current_script->can_set()) return;
    std::string initial_serialized_vars;
    if (preset_type <= 0)
        initial_serialized_vars = (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
        ? current_script->tab()->m_preset_bundle->fff_prints.get_selected_preset().config.opt_string("print_custom_variables")
        : current_script->tab()->m_preset_bundle->sla_prints.get_selected_preset().config.opt_string("print_custom_variables");
    else if (preset_type == 1) {
        initial_serialized_vars = (current_script->tab()->get_printer_technology() & PrinterTechnology::ptFFF) != 0
            ? current_script->tab()->m_preset_bundle->filaments.get_selected_preset().config.opt_string("filament_custom_variables", (unsigned int)(0))
            : current_script->tab()->m_preset_bundle->sla_materials.get_selected_preset().config.opt_string("filament_custom_variables", (unsigned int)(0));
    } else initial_serialized_vars = current_script->tab()->m_preset_bundle->printers.get_selected_preset().config.opt_string("printer_custom_variables");
    std::string serialized_value = get_custom_value(initial_serialized_vars, key);
    std::string serialized_vars = get_custom_var_option(preset_type);
    set_custom_value(serialized_vars, key, serialized_value);
    set_custom_option(preset_type, serialized_vars);
}

/////// main script fucntions //////

//TODO: add "unset" function, that revert to last value (befoer a scripted set) if a set has been made since last not-scripted change.
void ScriptContainer::init(const std::string& tab_key, Tab* tab)
{
    m_tab = tab;
    const boost::filesystem::path ui_script_file = Slic3r::GUI::get_app_config()->layout_config_path() / (tab_key + ".as");
    if (boost::filesystem::exists(ui_script_file)) {
        //launch the engine if not yet
        if (m_script_engine.get() == nullptr) {

            m_script_engine.reset(AngelScript::asCreateScriptEngine());
            // Create the script engine
            if (m_script_engine.get() == nullptr)
            {
                throw ScriptError("Failed to create script engine.");
            }
            // The script compiler will send any compiler messages to the callback function
#ifdef AS_MAX_PORTABILITY
            m_script_engine->SetMessageCallback(WRAP_FN(as_message_callback), 0, AngelScript::asCALL_GENERIC);
#else
            m_script_engine->SetMessageCallback(AngelScript::asFUNCTION(as_message_callback), 0, AngelScript::asCALL_CDECL);
#endif
            // Configure the script engine with the callback function
            AngelScript::RegisterScriptArray(m_script_engine.get(), false);
            AngelScript::RegisterStdString(m_script_engine.get());
            AngelScript::RegisterStdStringUtils(m_script_engine.get());
            AngelScript::RegisterScriptMath(m_script_engine.get());

            //if (use_generic) {
            //    r = engine->RegisterGlobalFunction("void print(string & in)", WRAP_FN(print), asCALL_GENERIC); assert(r >= 0);
            //} else {
            //    r = engine->RegisterGlobalFunction("void print(string & in)", asFUNCTION(print), asCALL_CDECL); assert(r >= 0);
            //}
#ifdef AS_MAX_PORTABILITY
            m_script_engine.get()->RegisterGlobalFunction("void print(string &in)", WRAP_FN(as_print), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void print_float(float)", WRAP_FN(as_print_float), AngelScript::asCALL_GENERIC);
            //m_script_engine.get()->RegisterGlobalFunction("void register_key(string &in)", WRAP_FN(as_register_key), AngelScript::asCALL_GENERIC);

            m_script_engine.get()->RegisterGlobalFunction("bool get_bool(string &in)", WRAP_FN(as_get_bool), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_bool(string &in, bool new_val)", WRAP_FN(as_set_bool), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("int get_int(string &in)", WRAP_FN(as_get_int), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_int(string &in, int new_val)", WRAP_FN(as_set_int), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("float get_float(string &in)", WRAP_FN(as_get_float), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_float(string &in, float new_val)", WRAP_FN(as_set_float), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("bool is_percent(string &in)", WRAP_FN(as_is_percent), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_percent(string &in, float new_val)", WRAP_FN(as_set_percent), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void get_string(string &in, string &out get_val)", WRAP_FN(as_get_string), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_string(string &in, string &in new_val)", WRAP_FN(as_set_string), AngelScript::asCALL_GENERIC);


            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_bool(int, string &in, bool &out)", WRAP_FN(as_get_custom_bool), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_bool(int, string &in, bool)", WRAP_FN(as_set_custom_bool), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_int(int, string &in, int &out)", WRAP_FN(as_get_custom_int), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_int(int, string &in, int)", WRAP_FN(as_set_custom_int), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_float(int, string &in, float &out)", WRAP_FN(as_get_custom_float), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_float(int, string &in, float)", WRAP_FN(as_set_custom_float), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_string(int, string &in, string &out)", WRAP_FN(as_get_custom_string), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_string(int, string &in, string &in)", WRAP_FN(as_set_custom_string), AngelScript::asCALL_GENERIC);

            m_script_engine.get()->RegisterGlobalFunction("float get_computed_float(string &in)", WRAP_FN(as_get_computed_float), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void back_initial_value(string &in)", WRAP_FN(as_back_initial_value), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void back_custom_initial_value(int, string &in)", WRAP_FN(as_back_custom_initial_value), AngelScript::asCALL_GENERIC);
            m_script_engine.get()->RegisterGlobalFunction("void ask_for_refresh()", WRAP_FN(as_ask_for_refresh), AngelScript::asCALL_GENERIC);

#else
            m_script_engine.get()->RegisterGlobalFunction("void print(string &in)",     AngelScript::asFUNCTION(as_print),          AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void print_float(float)",    AngelScript::asFUNCTION(as_print_float),    AngelScript::asCALL_CDECL);
            //m_script_engine.get()->RegisterGlobalFunction("void register_key(string &in)", AngelScript::asFUNCTION(as_register_key), AngelScript::asCALL_CDECL);

            m_script_engine.get()->RegisterGlobalFunction("bool get_bool(string &in)",                          AngelScript::asFUNCTION(as_get_bool),   AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_bool(string &in, bool new_val)",            AngelScript::asFUNCTION(as_set_bool),   AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("int get_int(string &in)",                            AngelScript::asFUNCTION(as_get_int),    AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_int(string &in, int new_val)",              AngelScript::asFUNCTION(as_set_int),    AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("float get_float(string &in)",                        AngelScript::asFUNCTION(as_get_float),  AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_float(string &in, float new_val)",          AngelScript::asFUNCTION(as_set_float),  AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool is_percent(string &in)",                        AngelScript::asFUNCTION(as_is_percent), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_percent(string &in, float new_val)",        AngelScript::asFUNCTION(as_set_percent),AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void get_string(string &in, string &out get_val)",   AngelScript::asFUNCTION(as_get_string), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_string(string &in, string &in new_val)",    AngelScript::asFUNCTION(as_set_string), AngelScript::asCALL_CDECL);

            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_bool(int, string &in, bool &out)",       AngelScript::asFUNCTION(as_get_custom_bool),    AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_bool(int, string &in, bool)",            AngelScript::asFUNCTION(as_set_custom_bool),    AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_int(int, string &in, int &out)",         AngelScript::asFUNCTION(as_get_custom_int),     AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_int(int, string &in, int)",              AngelScript::asFUNCTION(as_set_custom_int),     AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_float(int, string &in, float &out)",     AngelScript::asFUNCTION(as_get_custom_float),   AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_float(int, string &in, float)",          AngelScript::asFUNCTION(as_set_custom_float),   AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool get_custom_string(int, string &in, string &out)",   AngelScript::asFUNCTION(as_get_custom_string),  AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_custom_string(int, string &in, string &in)",    AngelScript::asFUNCTION(as_set_custom_string),  AngelScript::asCALL_CDECL);

            m_script_engine.get()->RegisterGlobalFunction("float get_computed_float(string &in)",   AngelScript::asFUNCTION(as_get_computed_float), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void back_initial_value(string &in)",    AngelScript::asFUNCTION(as_back_initial_value), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void back_custom_initial_value(int, string &in)",    AngelScript::asFUNCTION(as_back_custom_initial_value), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void ask_for_refresh()",                 AngelScript::asFUNCTION(as_ask_for_refresh),    AngelScript::asCALL_CDECL);
#endif
        }

        //m_script_module = m_script_engine->GetModule(tab_key.c_str(), AngelScript::asGM_CREATE_IF_NOT_EXISTS);
        AngelScript::CScriptBuilder builder;
        int res = builder.StartNewModule(m_script_engine.get(), tab_key.c_str());
        if (res < 0) throw CompileErrorException("Error, can't build the script for tab " + tab_key);
        // Let the builder load the script, and do the necessary pre-processing (include files, etc)
        //res = builder.AddSectionFromFile(ui_script_file.string().c_str()); //seems to be problematic on cyrillic locale
        {
            std::string all_file;
            boost::filesystem::load_string_file(ui_script_file, all_file);
            res = builder.AddSectionFromMemory(ui_script_file.string().c_str(), all_file.c_str(), (unsigned int)(all_file.length()), 0);
        }
        if (res < 0) throw CompileErrorException("Error, can't build the script for tab " + tab_key);
        res = builder.BuildModule();
        if (res < 0) throw CompileErrorException("Error, can't build the script for tab " + tab_key);
        m_script_module = m_script_engine->GetModule(tab_key.c_str(), AngelScript::asGM_ONLY_IF_EXISTS);
        //AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl("void main()");
        //AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
        //ctx->Prepare(func);

        //script_current_tab = this;
        //res = ctx->Execute();
        //script_current_tab = nullptr;
        //std::cout << "\nres is " << res << "\n";
        m_initialized = true;
    } else {
        BOOST_LOG_TRIVIAL(warning) << "Warning, can't find file script '" << ui_script_file.string() << "', is it needed? Even if something need it, it won't execute.";
        m_script_module = nullptr;
        disable();
    }
}

std::string get_type_name(ConfigOptionType type)
{
    switch (type) {
    case coFloat: return "float";
    case coFloats: return "float";
    case coInt: return "int";
    case coInts: return "int";
    case coString: return "string &in";
    case coStrings: return "string &in";
    case coPercent: return "float";
    case coPercents: return "float";
    case coFloatOrPercent: return "float, bool";
    case coFloatsOrPercents: return "float, bool";
    case coPoint: return "float, float";
    case coPoints: return "float, float";
    case coPoint3: return "float, float, float";
    case coBool: return "bool";
    case coBools: return "bool";
    case coEnum: return "string &in, int idx";
    default: return "";
    }
}

void ScriptContainer::call_script_function_set(const ConfigOptionDef& def, const boost::any& value)
{
    if (value.empty() || !is_intialized())
        return;
    std::string func_name = ("void " + def.opt_key + "_set(" + get_type_name(def.type) + ")");
    AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl(func_name.c_str());
    if (func == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't find function '" << func_name << "' in the script file";
        return;
    }
    AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
    if (ctx == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't create script context for function '" << func_name << "'";
        return;
    }
    ctx->Prepare(func);
    std::string str_arg;
    switch (def.type) {
    case coBool:
    case coBools: ctx->SetArgByte(0, boost::any_cast<bool>(value)); break;
    case coInt:
    case coInts: ctx->SetArgDWord(0, boost::any_cast<int>(value)); break;
    case coPercent:
    case coPercents:
    case coFloat:
    case coFloats: ctx->SetArgFloat(0, (float)boost::any_cast<double>(value)); break;
    case coFloatOrPercent:
    case coFloatsOrPercents: {
        std::string flOrPercent = boost::any_cast<std::string>(value);
        float val = 0;
        bool is_percent = false;
        if (flOrPercent[flOrPercent.size() - 1] == '%') {
            flOrPercent = flOrPercent.substr(0, flOrPercent.size() - 1);
            val = std::stof(flOrPercent);
            is_percent = true;
        } else {
            val = std::stof(flOrPercent);
        }
        ctx->SetArgDWord(0, boost::any_cast<float>((float)val));
        ctx->SetArgByte(1, boost::any_cast<bool>(is_percent));
        break;
    }
    case coPoint:
    case coPoints: { ctx->SetArgFloat(0, (float)boost::any_cast<double>(value)); ctx->SetArgFloat(1, (float)boost::any_cast<double>(value)); break; } //FIXME
    case coPoint3: { ctx->SetArgFloat(0, (float)boost::any_cast<double>(value)); ctx->SetArgFloat(1, (float)boost::any_cast<double>(value)); ctx->SetArgFloat(2, (float)boost::any_cast<double>(value)); break; }
    case coString:
    case coStrings: {
        str_arg = boost::any_cast<std::string>(value);
        ctx->SetArgAddress(0, &str_arg);
        break;
    }
    case coEnum: {
        int32_t enum_idx = boost::any_cast<std::int32_t>(value);
        if (enum_idx >= 0 && enum_idx < def.enum_values.size()) {
            str_arg = def.enum_values[enum_idx];
            ctx->SetArgAddress(0, &str_arg);
            ctx->SetArgDWord(1, enum_idx);
        }
        break;
    }
    }
    // init globals for script exec (TODO find a way to change that)
    current_script = this;
    m_need_refresh = false;
    m_to_update.clear();
    for (Tab* tab : wxGetApp().tabs_list)
        if (tab->completed())
            m_to_update[tab->type()] = {};
    m_can_set = true;
    // exec
    /*int res = */ctx->Execute();
    m_can_set = false;
    std::map<Preset::Type, DynamicPrintConfig> to_update = m_to_update;
    m_to_update.clear();
    auto to_reset = m_to_reset_initial;
    m_to_reset_initial.clear();

    //update the tabs from the results
    for (auto& data : to_update) {
        Tab* tab = wxGetApp().get_tab(data.first);
        //also reset
        if (!to_reset.empty()) {
            const DynamicPrintConfig& initial_conf = tab->m_presets->get_selected_preset().config;
            for (size_t key_idx = 0; key_idx != to_reset.size(); ++key_idx) {
                const std::string& key = to_reset[key_idx];
                if (initial_conf.has(key)) {
                    data.second.set_key_value(key, initial_conf.option(key)->clone());
                    to_reset.erase(to_reset.begin() + key_idx);
                    key_idx--;
                }
            }
        }
        tab->load_config(data.second);
    }
    //also call for value_changed, as it's not really a load but a change
    for (const auto& data : to_update) {
        Tab* tab = wxGetApp().get_tab(data.first);
        for (auto opt_key : data.second.keys()) {
            tab->on_value_change(opt_key, data.second.option(opt_key)->getAny());
        }
    }
    // refresh the field if needed
    if (m_need_refresh && m_tab) {
        Field* f = m_tab->get_field(def.opt_key);
        if (f != nullptr) {
            f->set_value(call_script_function_get_value(def), false);
        }
    }
}

bool ScriptContainer::call_script_function_reset(const ConfigOptionDef& def)
{
    std::string func_name = ("void " + def.opt_key + "_reset()");
    AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl(func_name.c_str());
    if (func == nullptr) {
        return false;
    }
    AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
    if (ctx == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't create script context for function '" << func_name << "'";
        return false;
    }
    ctx->Prepare(func);
    m_can_set = true;
    // exec
    /*int res = */ctx->Execute();
    m_can_set = false;
    std::map<Preset::Type, DynamicPrintConfig> to_update = m_to_update;
    m_to_update.clear();
    auto to_reset = m_to_reset_initial;
    m_to_reset_initial.clear();

    //update the tabs from the results
    for (auto& data : to_update) {
        Tab* tab = wxGetApp().get_tab(data.first);
        //also reset
        if (!to_reset.empty()) {
            const DynamicPrintConfig& initial_conf = tab->m_presets->get_selected_preset().config;
            for (size_t key_idx = 0; key_idx != to_reset.size(); ++key_idx) {
                const std::string& key = to_reset[key_idx];
                if (initial_conf.has(key)) {
                    data.second.set_key_value(key, initial_conf.option(key)->clone());
                    to_reset.erase(to_reset.begin() + key_idx);
                    key_idx--;
                }
            }
        }
        tab->load_config(data.second);
    }
    //also call for value_changed, as it's not really a load but a change
    for (const auto& data : to_update) {
        Tab* tab = wxGetApp().get_tab(data.first);
        for (auto opt_key : data.second.keys()) {
            tab->on_value_change(opt_key, data.second.option(opt_key)->getAny());
        }
    }
    // refresh the field if needed
    if (m_need_refresh && m_tab) {
        Field* f = m_tab->get_field(def.opt_key);
        if (f != nullptr) {
            f->set_value(call_script_function_get_value(def), false);
        }
    }
    return true;
}
 
//void ScriptContainer::call_script_function_refresh(const std::string& def_id)
//{
//    std::string func_name = ("int " + def_id + "_refresh()");
//    AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl(func_name.c_str());
//    if (func == nullptr) {
//        BOOST_LOG_TRIVIAL(error) << "Error, can't find function '" << func_name << "' in the script file";
//        return;
//    }
//    AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
//    if (ctx == nullptr) {
//        BOOST_LOG_TRIVIAL(error) << "Error, can't create script context for function '" << func_name << "'";
//        return;
//    }
//    ctx->Prepare(func);
//    // init globals for script exec (TODO find a way to change that)
//    script_current_tab = m_tab;
//    current_script->tech() = m_tech;
//    // exec
//    int res = ctx->Execute();
//    int ret = ctx->GetReturnDWord();
//    if (ret >= 0) {
//        m_tab->set_value(def_id, unsigned char(ret));
//    } else {
//        m_tab->set_value(def_id, unsigned char(2));
//    }
//    //TODO: add the keyt into a collection of dirty script-widget in our tab. Then, ask for update_dirty() and add the code to use that collection in update_changed_ui
//    m_tab->add_dirty_setting(def_id);
//    //m_tab->update_dirty();
//    //m_tab->
//}

boost::any ScriptContainer::call_script_function_get_value(const ConfigOptionDef& def)
{
    if (!is_intialized())
        return boost::any();

    std::string func_name;

    switch (def.type) {
    case coBool:
    case coBools:
    case coInt:
    case coInts: func_name = "int"; break;
    case coPercent:
    case coPercents:
    case coFloat:
    case coFloats: 
    case coFloatOrPercent:
    case coFloatsOrPercents: func_name = "float"; break;
    case coPoint:
    case coPoints: func_name = "float"; break; //FIXME
    case coPoint3: func_name = "float"; break;
    case coString:
    case coStrings:func_name = "void"; break;
    case coEnum: func_name = "int"; break;
    }
    func_name += (" " + def.opt_key + "_get(");
    switch (def.type) {
    case coFloatOrPercent:
    case coFloatsOrPercents: func_name += "bool &out"; break;
    case coString:
    case coStrings:
    case coEnum: func_name += "string &out"; break;
    }
    func_name += ")";
    AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl(func_name.c_str());
    if (func == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't find function '" << func_name << "' in the script file";
        return boost::any{};
    }
    AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
    if (ctx == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't create script context for function '" << func_name << "'";
        return boost::any{};
    }
    ctx->Prepare(func);
    std::string ret_str;
    bool ret_percent = false;
    switch (def.type) {
    case coFloatOrPercent:
    case coFloatsOrPercents: ctx->SetArgAddress(0, &ret_percent); break;
    case coString:
    case coStrings:
    case coEnum: ctx->SetArgObject(0, &ret_str); break;
    }
    // init globals for script exec (TODO find a way to change that)
    current_script = this;
    m_need_refresh = false;
    // exec
    int res = ctx->Execute();
    int32_t ret_int;
    float ret_float;
    boost::any field_val;
    boost::any opt_val;
    switch (def.type) {
    case coBool:
    case coBools: { ret_int = ctx->GetReturnDWord(); field_val = uint8_t(ret_int < 0 ? 2 : ret_int); opt_val = uint8_t((ret_int > 0)?1:0); break; } //CheckBox
    case coInt:
    case coInts: { ret_int = ctx->GetReturnDWord(); field_val = int32_t(ret_int); opt_val = int(ret_int); break; } //SpinCtrl
    case coString:
    case coStrings: { field_val = from_u8(ret_str); opt_val = ret_str; break; } //TextCtrl
    case coPercent:
    case coPercents: ret_percent = true;
    case coFloat:
    case coFloats: opt_val = double(ctx->GetReturnFloat());
    case coFloatOrPercent:
    case coFloatsOrPercents:
    {
        ret_float = ctx->GetReturnFloat();
        wxString ret_wstring = double_to_string(ret_float);
        if (ret_percent)
            ret_wstring += '%';
        field_val = ret_wstring; //TextCtrl
        if (opt_val.empty()) { opt_val = ret_wstring.ToStdString(); }
        break;
    }
    case coPoint:
    case coPoints: { ret_float = ctx->GetReturnFloat(); field_val = Vec2d{ ret_float, ret_float }; opt_val = double(ctx->GetReturnFloat()); break; } //FIXME PointCtrl
    case coPoint3: { ret_float = ctx->GetReturnFloat(); field_val = Vec3d{ ret_float, ret_float, ret_float }; opt_val = double(ctx->GetReturnFloat());  break; }
    case coEnum: { 
        ret_int = ctx->GetReturnDWord();
        if (ret_int >= 0 && ret_int < def.enum_values.size()) {
            field_val = int32_t(ret_int);
        } else {
            field_val = int32_t(0);
            for (size_t i = 0; i < def.enum_values.size(); i++) {
                if (ret_str == def.enum_values[i])
                    field_val = int32_t(i);
            }
        }
        opt_val = field_val;
        break; //Choice
    }
    }
    if (m_need_refresh) {
        refresh(def, opt_val);
    }
    if (field_val.empty()) {
        std::cout << "Error nullptr for script\n";
    }
    return field_val;
}

void ScriptContainer::refresh(const ConfigOptionDef& def, boost::any value)
{
    auto it = std::find(m_currently_reset.begin(), m_currently_reset.end(), def.opt_key);
    // Don't refresh while in refresh
    if (it != m_currently_reset.end())
        return;
    //if bool, need to change the int in bool
    if (def.type == coBool || def.type == coBools) {
        value = bool(boost::any_cast<uint8_t>(value) != 0);
    }
    m_currently_reset.push_back(def.opt_key);
    call_script_function_set(def, value);
    //remove opt_key
    it = std::find(m_currently_reset.begin(), m_currently_reset.end(), def.opt_key);
    if (it != m_currently_reset.end())
        m_currently_reset.erase(it);
}

//TODO find a way to use the depends_on to add the same lock & points as real configoption in the gui

} }//namespace Slic3r Gui

