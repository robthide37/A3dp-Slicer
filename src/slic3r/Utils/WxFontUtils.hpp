#ifndef slic3r_WxFontUtils_hpp_
#define slic3r_WxFontUtils_hpp_

#include <memory>
#include <optional>
#include <wx/font.h>
#include "libslic3r/Emboss.hpp"

namespace Slic3r::GUI {

// Help class to  work with wx widget font object( wxFont )
class WxFontUtils
{
public:
    // only static functions
    WxFontUtils() = delete;

    // check if exist file for wxFont
    // return pointer on data or nullptr when can't load
    static void *can_load(const wxFont &font);

    // os specific load of wxFont
    static std::unique_ptr<Slic3r::Emboss::FontFile> create_font_file(const wxFont &font);

    static FontItem::Type get_actual_type();
    static FontItem       get_font_item(const wxFont &font);

    // load font used by Operating system as default GUI
    static FontItem    get_os_font();
    static std::string get_human_readable_name(const wxFont &font);

    // serialize / deserialize font
    static std::string           store_wxFont(const wxFont &font);
    static std::optional<wxFont> load_wxFont(const std::string &font_descriptor);

    // Try to create similar font, loaded from 3mf from different Computer
    static std::optional<wxFont> create_wxFont(const FontItem &fi,
                                               const FontProp &fp);
    // update font property by wxFont
    static void update_property(FontProp &font_prop, const wxFont &font);

    static bool is_italic(const wxFont &font);
    static bool is_bold(const wxFont &font);

    // Font could not support italic than return FALSE.
    // For check of support is neccessary font file pointer.
    // Font file is optional it could be created inside of function, but it slow down.
    // To not load font file twice on success font_file contain new created font file.
    static bool set_italic(wxFont &font, std::shared_ptr<Emboss::FontFile>& font_file = std::shared_ptr<Emboss::FontFile>());

    // Font could not support bold than return FALSE.
    // For check of support is neccessary font file pointer.
    // Font file is optional it could be created inside of function, but it slow down.
    // To not load font file twice on success font_file contain new created font file.
    static bool set_bold(wxFont &font, std::shared_ptr<Emboss::FontFile>& font_file = std::shared_ptr<Emboss::FontFile>());

    // map to convert wxFont type to string and vice versa
    static const std::map<wxFontFamily, std::string> from_family;
    static const std::map<std::string, wxFontFamily> to_family;

    static const std::map<wxFontStyle, std::string> from_style;
    static const std::map<std::string, wxFontStyle> to_style;

    static const std::map<wxFontWeight, std::string> from_weight;
    static const std::map<std::string, wxFontWeight> to_weight;
};

} // namespace Slic3r::GUI
#endif // slic3r_WxFontUtils_hpp_
