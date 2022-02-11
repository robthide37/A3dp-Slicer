#ifndef slic3r_ScriptExecutor_hpp_
#define slic3r_ScriptExecutor_hpp_

#include "libslic3r/Config.hpp"
#include "OptionsGroup.hpp"

#include <angelscript/include/angelscript.h>
#include <angelscript/add_on/scriptbuilder/scriptbuilder.h>

namespace Slic3r { namespace GUI {

class Tab;

class ScriptContainer
{
    Tab* m_tab;
    PrinterTechnology   m_tech = ptFFF;
public:
    ScriptContainer() {}
    inline static AngelScript::PtrRelease<AngelScript::asIScriptEngine> m_script_engine;
    AngelScript::asIScriptModule* m_script_module{ nullptr };

    void init(const std::string& resource_dir, const std::string& tab_key, Tab* tab, PrinterTechnology current_tech);
    void call_script_function_set(const Option& def, const boost::any& value);
    //void call_script_function_refresh(const std::string& def_id);
    boost::any call_script_function_get_value(const ConfigOptionDef& def);
};

} } //namespace Slic3r Gui

#endif //slic3r_ScriptExecutor_hpp_
