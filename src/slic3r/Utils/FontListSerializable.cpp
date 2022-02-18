#include "FontListSerializable.hpp"

#include <libslic3r/AppConfig.hpp>
#include "WxFontUtils.hpp"

using namespace Slic3r;
using namespace Slic3r::GUI;

const std::string FontListSerializable::APP_CONFIG_FONT_NAME        = "name";
const std::string FontListSerializable::APP_CONFIG_FONT_DESCRIPTOR  = "descriptor";
const std::string FontListSerializable::APP_CONFIG_FONT_LINE_HEIGHT = "line_height";
const std::string FontListSerializable::APP_CONFIG_FONT_DEPTH       = "depth";
const std::string FontListSerializable::APP_CONFIG_FONT_BOLDNESS    = "boldness";
const std::string FontListSerializable::APP_CONFIG_FONT_SKEW        = "skew";
const std::string FontListSerializable::APP_CONFIG_FONT_DISTANCE    = "distance";
const std::string FontListSerializable::APP_CONFIG_FONT_ANGLE       = "angle";
const std::string FontListSerializable::APP_CONFIG_FONT_CHAR_GAP    = "char_gap";
const std::string FontListSerializable::APP_CONFIG_FONT_LINE_GAP    = "line_gap";

std::string FontListSerializable::create_section_name(unsigned index)
{
    return AppConfig::SECTION_FONT + ':' + std::to_string(index);
}

#include "fast_float/fast_float.h"
bool FontListSerializable::read(const std::map<std::string, std::string>& section, const std::string& key, float& value){
    auto item = section.find(key);
    if (item == section.end()) return false;
    const std::string &data = item->second;
    if (data.empty()) return false;
    float value_;
    fast_float::from_chars(data.c_str(), data.c_str() + data.length(), value_);
    // read only non zero value
    if (fabs(value_) <= std::numeric_limits<float>::epsilon()) return false;

    value = value_;
    return true;
}

bool FontListSerializable::read(const std::map<std::string, std::string>& section, const std::string& key, std::optional<int>& value){
    auto item = section.find(key);
    if (item == section.end()) return false;
    const std::string &data = item->second;
    if (data.empty()) return false;
    int value_ = std::atoi(data.c_str());
    if (value_ == 0) return false;

    value = value_;
    return true;
}

bool FontListSerializable::read(const std::map<std::string, std::string>& section, const std::string& key, std::optional<float>& value){
    auto item = section.find(key);
    if (item == section.end()) return false;
    const std::string &data = item->second;
    if (data.empty()) return false;
    float value_;
    fast_float::from_chars(data.c_str(), data.c_str() + data.length(), value_);
    // read only non zero value
    if (fabs(value_) <= std::numeric_limits<float>::epsilon()) return false;

    value = value_;
    return true;
}

std::optional<FontItem> FontListSerializable::load_font_item(
    const std::map<std::string, std::string> &app_cfg_section)
{
    auto path_it = app_cfg_section.find(APP_CONFIG_FONT_DESCRIPTOR);
    if (path_it == app_cfg_section.end()) return {};
    const std::string &path = path_it->second;

    auto name_it = app_cfg_section.find(APP_CONFIG_FONT_NAME);
    static const std::string default_name = "font_name";
    const std::string &name = 
        (name_it == app_cfg_section.end()) ?
        default_name : name_it->second;
        
    FontProp fp;
    read(app_cfg_section, APP_CONFIG_FONT_LINE_HEIGHT, fp.size_in_mm);
    read(app_cfg_section, APP_CONFIG_FONT_DEPTH, fp.emboss);
    read(app_cfg_section, APP_CONFIG_FONT_BOLDNESS, fp.boldness);
    read(app_cfg_section, APP_CONFIG_FONT_SKEW, fp.skew);
    read(app_cfg_section, APP_CONFIG_FONT_DISTANCE, fp.distance);
    read(app_cfg_section, APP_CONFIG_FONT_ANGLE, fp.angle);
    read(app_cfg_section, APP_CONFIG_FONT_CHAR_GAP, fp.char_gap);
    read(app_cfg_section, APP_CONFIG_FONT_LINE_GAP, fp.line_gap);

    FontItem::Type type = WxFontUtils::get_actual_type();
    return FontItem(name, path, type, fp);
}

void FontListSerializable::store_font_item(AppConfig &     cfg,
                                           const FontItem &fi,
                                           unsigned        index)
{
    std::string section_name = create_section_name(index);
    cfg.clear_section(section_name);
    cfg.set(section_name, APP_CONFIG_FONT_NAME, fi.name);
    cfg.set(section_name, APP_CONFIG_FONT_DESCRIPTOR, fi.path);
    const FontProp &fp = fi.prop;
    cfg.set(section_name, APP_CONFIG_FONT_LINE_HEIGHT, std::to_string(fp.size_in_mm));
    cfg.set(section_name, APP_CONFIG_FONT_DEPTH, std::to_string(fp.emboss));
    if (fp.boldness.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_BOLDNESS, std::to_string(*fp.boldness));
    if (fp.skew.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_SKEW, std::to_string(*fp.skew));
    if (fp.distance.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_DISTANCE, std::to_string(*fp.distance));
    if (fp.angle.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_ANGLE, std::to_string(*fp.angle));
    if (fp.char_gap.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_CHAR_GAP, std::to_string(*fp.char_gap));
    if (fp.line_gap.has_value())
        cfg.set(section_name, APP_CONFIG_FONT_LINE_GAP, std::to_string(*fp.line_gap));
}