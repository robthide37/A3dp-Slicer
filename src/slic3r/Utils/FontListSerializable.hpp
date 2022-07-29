#ifndef slic3r_FontListSerializable_hpp_
#define slic3r_FontListSerializable_hpp_

#include <string>
#include <map>
#include <optional>
#include <libslic3r/TextConfiguration.hpp> // FontList+FontItem

namespace Slic3r {
class AppConfig;
}

namespace Slic3r::GUI {

/// <summary>
/// For store/load font list to/from AppConfig
/// </summary>
class FontListSerializable
{
    static const std::string APP_CONFIG_FONT_NAME;
    static const std::string APP_CONFIG_FONT_DESCRIPTOR;
    static const std::string APP_CONFIG_FONT_LINE_HEIGHT;
    static const std::string APP_CONFIG_FONT_DEPTH;
    static const std::string APP_CONFIG_FONT_USE_SURFACE;
    static const std::string APP_CONFIG_FONT_BOLDNESS;
    static const std::string APP_CONFIG_FONT_SKEW;
    static const std::string APP_CONFIG_FONT_DISTANCE;
    static const std::string APP_CONFIG_FONT_ANGLE;
    static const std::string APP_CONFIG_FONT_COLLECTION;
    static const std::string APP_CONFIG_FONT_CHAR_GAP;
    static const std::string APP_CONFIG_FONT_LINE_GAP;

    static const std::string APP_CONFIG_ACTIVE_FONT;
public:
    FontListSerializable() = delete;

    static void store_font_index(AppConfig &cfg, unsigned index);
    static std::optional<size_t> load_font_index(const AppConfig &cfg);

    static FontList load_font_list(const AppConfig &cfg);
    static void     store_font_list(AppConfig &cfg, const FontList font_list);

private:
    static std::string create_section_name(unsigned index);
    static std::optional<FontItem> load_font_item(const std::map<std::string, std::string> &app_cfg_section);
    static void store_font_item(AppConfig &cfg, const FontItem &fi, unsigned index);

    // TODO: move to app config like read from section
    static bool read(const std::map<std::string, std::string>& section, const std::string& key, bool& value);
    static bool read(const std::map<std::string, std::string>& section, const std::string& key, float& value);
    static bool read(const std::map<std::string, std::string>& section, const std::string& key, std::optional<int>& value);
    static bool read(const std::map<std::string, std::string>& section, const std::string& key, std::optional<unsigned int>& value);
    static bool read(const std::map<std::string, std::string>& section, const std::string& key, std::optional<float>& value);    
};
} // namespace Slic3r

#endif // #define slic3r_FontListSerializable_hpp_

