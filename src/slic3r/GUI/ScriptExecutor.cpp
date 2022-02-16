#include "ScriptExecutor.hpp"
#include "GUI_App.hpp"
#include "Tab.hpp"
#include "libslic3r/PresetBundle.hpp"

#include <string>

#include <angelscript/add_on/scriptarray/scriptarray.h>
#include <angelscript/add_on/scriptbuilder/scriptbuilder.h>
#include <angelscript/add_on/scriptstdstring/scriptstdstring.h>

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
Tab* script_current_tab;
PrinterTechnology current_tech;
bool can_set = false;
std::map<Preset::Type, DynamicPrintConfig> to_update;
std::vector<std::string>* watching_keys;
void as_print(std::string& str)
{
    std::cout << str;
}
void as_print_float(float f)
{
    std::cout << f;
}
void as_register_key(std::string& key) {
    if (watching_keys != nullptr)
        watching_keys->push_back(key);
}
std::pair<PresetCollection*, ConfigOption*> get_coll(const std::string& str) {
    PresetCollection* coll = current_tech == (PrinterTechnology::ptFFF)
        ? &script_current_tab->m_preset_bundle->fff_prints
        : &script_current_tab->m_preset_bundle->sla_prints;
    ConfigOption* opt = coll->get_edited_preset().config.option(str);
    if (opt == nullptr) {
        coll = current_tech == (PrinterTechnology::ptFFF)
            ? &script_current_tab->m_preset_bundle->filaments
            : &script_current_tab->m_preset_bundle->sla_materials;
        opt = coll->get_edited_preset().config.option(str);
    }
    if (opt == nullptr) {
        coll = &script_current_tab->m_preset_bundle->printers;
        opt = coll->get_edited_preset().config.option(str);
    }
    return { coll,  opt };
}
bool as_get_bool(std::string& str)
{
    const ConfigOption* opt = get_coll(str).second;
    if (opt == nullptr || opt->type() != ConfigOptionType::coBool)
        throw Exception("error, can't find bool option " + str);
    return opt->getBool();
}
void as_set_bool(std::string& str, bool b)
{
    if (!can_set) return;
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(str);
    if (result.second == nullptr)
        throw Exception("error, can't find bool option " + str);
    DynamicPrintConfig& conf = to_update[result.first->type()];
    if (result.second->type() == ConfigOptionType::coBool)
        conf.set_key_value(str, new ConfigOptionBool(b));
}
int32_t as_get_int(std::string& str)
{
    const ConfigOption* opt = get_coll(str).second;
    if (opt == nullptr || (opt->type() != ConfigOptionType::coInt && opt->type() != ConfigOptionType::coEnum))
        throw Exception("error, can't find int option " + str);
    return (int32_t)(opt->getInt());
}
void as_set_int(std::string& key, int val)
{
    if (!can_set) return;
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw Exception("error, can't find int option " + key);
    DynamicPrintConfig& conf = to_update[result.first->type()];
    if (result.second->type() == ConfigOptionType::coInt)
        conf.set_key_value(key, new ConfigOptionInt(val));
    else if (result.second->type() == ConfigOptionType::coEnum) {
        const ConfigOptionDef* def = result.first->get_edited_preset().config.get_option_def(key);
        if (val >= 0 && val < def->enum_values.size()) {
            ConfigOption* copy = result.second->clone();
            copy->setInt(val);
            conf.set_key_value(key, copy);
        }
    }
}
float as_get_float(std::string& str)
{
    const ConfigOption* opt = get_coll(str).second;
    if (opt == nullptr) //TODO check if  float, etc..
        throw Exception("error, can't find float option " + str);
    return (float)(opt->getFloat());
}
void as_set_float(std::string& str, float f)
{
    if (!can_set) return;
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(str);
    if (result.second == nullptr)
        throw Exception("error, can't find float option " + str);
    DynamicPrintConfig& conf = to_update[result.first->type()];
    if (result.second->type() == ConfigOptionType::coFloat)
        conf.set_key_value(str, new ConfigOptionFloat(f));
    else if (result.second->type() == ConfigOptionType::coPercent)
        conf.set_key_value(str, new ConfigOptionPercent(f/100.));
    else if (result.second->type() == ConfigOptionType::coFloatOrPercent)
        conf.set_key_value(str, new ConfigOptionFloatOrPercent(f, false));
}
bool as_is_percent(std::string& str)
{
    const ConfigOption* opt = get_coll(str).second;
    if (opt == nullptr)
        throw Exception("error, can't find percent option " + str);
    return (opt->type() == ConfigOptionType::coPercent) || (opt->type() == ConfigOptionType::coFloatOrPercent &&
        ((ConfigOptionFloatOrPercent*)opt)->percent);
}
void as_set_percent(std::string& str, float f)
{
    if (!can_set) return;
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(str);
    if (result.second == nullptr)
        throw Exception("error, can't find percent option " + str);
    DynamicPrintConfig& conf = to_update[result.first->type()];
    if (result.second->type() == ConfigOptionType::coFloat)
        conf.set_key_value(str, new ConfigOptionFloat(f/100.));
    else if (result.second->type() == ConfigOptionType::coPercent)
        conf.set_key_value(str, new ConfigOptionPercent(f));
    else if (result.second->type() == ConfigOptionType::coFloatOrPercent)
        conf.set_key_value(str, new ConfigOptionFloatOrPercent(f, true));
}
void as_get_string(std::string& key, std::string& val)
{
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(key);
    const ConfigOption* opt = result.second;
    if (opt == nullptr) //TODO check if  float, etc..
        throw Exception("error, can't find string option " + key);
    if (opt->type() == ConfigOptionType::coString)
        val = ((ConfigOptionString*)opt)->value;
    else if (opt->type() == ConfigOptionType::coEnum) {
        int idx = opt->getInt();
        const ConfigOptionDef* def = result.first->get_edited_preset().config.get_option_def(key);
        if (idx >= 0 && idx < def->enum_values.size())
            val = def->enum_values[idx];
    }
}
void as_set_string(std::string& key, std::string& val)
{
    if (!can_set) return;
    std::pair<PresetCollection*, ConfigOption*> result = get_coll(key);
    if (result.second == nullptr)
        throw Exception("error, can't find string option " + key);
    DynamicPrintConfig& conf = to_update[result.first->type()];
    if (result.second->type() == ConfigOptionType::coString)
        conf.set_key_value(key, new ConfigOptionString(val));
    else if (result.second->type() == ConfigOptionType::coEnum) {
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
//function to reset a field
void as_back_initial_value(std::string& key) {
    for (Tab* tab : wxGetApp().tabs_list) {
        if (tab != nullptr && tab->supports_printer_technology(current_tech)) {
            // more optimal to create a tab::back_to_initial_value() that call the optgroup... here we are going down to the optgroup to get the field to get back to the opgroup via on_back_to_initial_value
            Field* f = tab->get_field(key);
            if (f != nullptr) {
                f->on_back_to_initial_value();
            }
        }
    }
}
//TODO: add "unset" function, that revert to last value (befoer a scripted set) if a set has been made since last not-scripted change.
void ScriptContainer::init(const std::string& resource_dir, const std::string& tab_key, Tab* tab, PrinterTechnology current_tech)
{
    m_tab = tab;
    m_tech = current_tech;
    const boost::filesystem::path ui_script_file = (boost::filesystem::path(resource_dir) / "ui_layout" / (tab_key + ".as")).make_preferred();
    if (boost::filesystem::exists(ui_script_file)) {
        //launch the engine if not yet
        if (m_script_engine.get() == nullptr) {

            m_script_engine.reset(AngelScript::asCreateScriptEngine());
            // Create the script engine
            if (m_script_engine.get() == nullptr)
            {
                std::cout << "Failed to create script engine." << std::endl;
                throw Exception("Failed to create script engine.");
            }
            // The script compiler will send any compiler messages to the callback function
            m_script_engine->SetMessageCallback(AngelScript::asFUNCTION(as_message_callback), 0, AngelScript::asCALL_CDECL);
            // Configure the script engine with the callback function
            AngelScript::RegisterScriptArray(m_script_engine.get(), false);
            AngelScript::RegisterStdString(m_script_engine.get());
            AngelScript::RegisterStdStringUtils(m_script_engine.get());
            m_script_engine.get()->RegisterGlobalFunction("void print(string &in)", AngelScript::asFUNCTION(as_print), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void print_float(float)", AngelScript::asFUNCTION(as_print_float), AngelScript::asCALL_CDECL);
            //m_script_engine.get()->RegisterGlobalFunction("void register_key(string &in)", AngelScript::asFUNCTION(as_register_key), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool get_bool(string &in)", AngelScript::asFUNCTION(as_get_bool), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_bool(string &in, bool new_val)", AngelScript::asFUNCTION(as_set_bool), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("int get_int(string &in)", AngelScript::asFUNCTION(as_get_int), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_int(string &in, int new_val)", AngelScript::asFUNCTION(as_set_int), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("float get_float(string &in)", AngelScript::asFUNCTION(as_get_float), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_float(string &in, float new_val)", AngelScript::asFUNCTION(as_set_float), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("bool is_percent(string &in)", AngelScript::asFUNCTION(as_is_percent), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_percent(string &in, float new_val)", AngelScript::asFUNCTION(as_set_percent), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void get_string(string &in, string &out get_val)", AngelScript::asFUNCTION(as_get_string), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void set_string(string &in, string &in new_val)", AngelScript::asFUNCTION(as_set_string), AngelScript::asCALL_CDECL);
            m_script_engine.get()->RegisterGlobalFunction("void back_initial_value(string &in)", AngelScript::asFUNCTION(as_back_initial_value), AngelScript::asCALL_CDECL);
        }

        //m_script_module = m_script_engine->GetModule(tab_key.c_str(), AngelScript::asGM_CREATE_IF_NOT_EXISTS);
        AngelScript::CScriptBuilder builder;
        int res = builder.StartNewModule(m_script_engine.get(), tab_key.c_str());
        if (res < 0) throw Exception("Error, can't build the script for tab " + tab_key);
        // Let the builder load the script, and do the necessary pre-processing (include files, etc)
        res = builder.AddSectionFromFile(ui_script_file.string().c_str());
        if (res < 0) throw Exception("Error, can't build the script for tab " + tab_key);
        res = builder.BuildModule();
        if (res < 0) throw Exception("Error, can't build the script for tab " + tab_key);
        m_script_module = m_script_engine->GetModule(tab_key.c_str(), AngelScript::asGM_ONLY_IF_EXISTS);
        //AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl("void main()");
        //AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
        //ctx->Prepare(func);

        //script_current_tab = this;
        //res = ctx->Execute();
        //script_current_tab = nullptr;
        //std::cout << "\nres is " << res << "\n";
    } else {
        m_script_module = nullptr;
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
void ScriptContainer::call_script_function_set(const Option& def, const boost::any& value)
{
    if (value.empty())
        return;
    std::string func_name = ("void " + def.opt.opt_key + "_set(" + get_type_name(def.opt.type) + ")");
    AngelScript::asIScriptFunction* func = m_script_module->GetFunctionByDecl(func_name.c_str());
    if (func == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't find function '" << func_name<<"' in the script file";
        return;
    }
    AngelScript::asIScriptContext* ctx = m_script_engine->CreateContext();
    if (ctx == nullptr) {
        BOOST_LOG_TRIVIAL(error) << "Error, can't create script context for function '" << func_name << "'";
        return;
    }
    ctx->Prepare(func);
    std::string str_arg;
    switch (def.opt.type) {
    case coBool:
    case coBools: ctx->SetArgByte(0, boost::any_cast<bool>(value)); break;
    case coInt:
    case coInts: ctx->SetArgDWord(0, boost::any_cast<int>(value)); break;
    case coPercent:
    case coPercents:
    case coFloat:
    case coFloats: ctx->SetArgFloat(0, boost::any_cast<float>(value)); break;
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
    case coPoints: { ctx->SetArgFloat(0, boost::any_cast<float>(value)); ctx->SetArgFloat(1, boost::any_cast<float>(value)); break; } //FIXME
    case coPoint3: { ctx->SetArgFloat(0, boost::any_cast<float>(value)); ctx->SetArgFloat(1, boost::any_cast<float>(value)); ctx->SetArgFloat(2, boost::any_cast<float>(value)); break; }
    case coString:
    case coStrings: {
        str_arg = boost::any_cast<std::string>(value);
        ctx->SetArgAddress(0, &str_arg);
        break;
    }
    case coEnum: { 
        int32_t enum_idx = boost::any_cast<std::int32_t>(value);
        if (enum_idx >= 0 && enum_idx < def.opt.enum_values.size()) {
            str_arg = def.opt.enum_values[enum_idx];
            ctx->SetArgAddress(0, &str_arg);
            ctx->SetArgDWord(1, enum_idx);
        }
        break; 
    }
    }
    // init globals for script exec (TODO find a way to change that)
    script_current_tab = m_tab;
    current_tech = m_tech;
    to_update.clear();
    can_set = true;
    // exec
    int res = ctx->Execute();
    can_set = false;
    //update the tabs from the results
    for (const auto& data : to_update) {
        Tab* tab = wxGetApp().get_tab(data.first);
        tab->load_config(data.second);
    }
    //int ret = ctx->GetReturnDWord();
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
//    current_tech = m_tech;
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
    script_current_tab = m_tab;
    current_tech = m_tech;
    // exec
    int res = ctx->Execute();
    int32_t ret_int;
    float ret_float;
    switch (def.type) {
    case coBool:
    case coBools: { ret_int = ctx->GetReturnDWord(); return uint8_t(ret_int < 0 ? 2 : ret_int); } //CheckBox
    case coInt:
    case coInts: { ret_int = ctx->GetReturnDWord(); return int32_t(ret_int); } //SpinCtrl
    case coString:
    case coStrings: { return from_u8(ret_str); break; } //TextCtrl
    case coPercent:
    case coPercents: ret_percent = true;
    case coFloatOrPercent:
    case coFloatsOrPercents:
    case coFloat:
    case coFloats: {
        ret_float = ctx->GetReturnFloat();
        wxString ret_wstring = double_to_string(ret_float);
        if (ret_percent)
            ret_wstring += '%';
        return ret_wstring; //TextCtrl
    }
    case coPoint:
    case coPoints: { ret_float = ctx->GetReturnFloat(); return Vec2d{ ret_float, ret_float };  } //FIXME PointCtrl
    case coPoint3: { ret_float = ctx->GetReturnFloat(); return Vec3d{ ret_float, ret_float, ret_float }; }
    case coEnum: { 
        ret_int = ctx->GetReturnDWord();
        if (ret_int >= 0 && ret_int < def.enum_values.size()) {
            return int32_t(ret_int);
        } else {
            for (size_t i = 0; i < def.enum_values.size(); i++) {
                if (ret_str == def.enum_values[i])
                    return int32_t(i);
            }
        }
        return int32_t(0); } //Choice
    }
    return boost::any{};
}

//TODO find a way to use the depends_on to add the same lock & points as real configoption in the gui

} }//namespace Slic3r Gui

