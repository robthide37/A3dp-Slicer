#include "WxFontUtils.hpp"
#include "libslic3r/MapUtils.hpp"

#if defined(__APPLE__)
#include <wx/uri.h>
#include <CoreText/CTFont.h>
#include <wx/osx/core/cfdictionary.h>
#elif defined(__linux__)
#include "slic3r/Utils/FontConfigHelp.hpp"
#endif

using namespace Slic3r;
using namespace Slic3r::GUI;

void *WxFontUtils::can_load(const wxFont &font)
{
    if (!font.IsOk()) return nullptr;
#ifdef _WIN32
    return Emboss::can_load(font.GetHFONT());
#elif defined(__APPLE__)
    // use file path
    return font.GetNativeFontInfo();
#elif defined(__linux__)
    // TODO: find better way
    static FontConfigHelp help;
    std::string font_path = help.get_font_path(font);
    if (font_path.empty()) return nullptr;
    return Emboss::load_font(font_path.c_str());
#endif
    return nullptr;
}

std::unique_ptr<Emboss::FontFile> WxFontUtils::create_font_file(const wxFont &font)
{
#ifdef _WIN32
    return Emboss::load_font(font.GetHFONT());
#elif defined(__APPLE__)
    // use file path
    const wxNativeFontInfo *info = font.GetNativeFontInfo();
    if (info == nullptr) return nullptr;
    CTFontDescriptorRef descriptor = info->GetCTFontDescriptor();
    CFURLRef            typeref    = (CFURLRef)
        CTFontDescriptorCopyAttribute(descriptor, kCTFontURLAttribute);
    CFStringRef url = CFURLGetString(typeref);
    if (url == NULL) return nullptr;
    wxString file_uri;
    wxCFTypeRef(url).GetValue(file_uri);
    std::string file_path(wxURI::Unescape(file_uri).c_str());
    size_t      start = std::string("file://").size();
    if (file_path.empty() || file_path.size() <= start) return nullptr;
    file_path = file_path.substr(start, file_path.size() - start);
    return Emboss::load_font(file_path.c_str());
#elif defined(__linux__)
    static FontConfigHelp help;
    std::string font_path = help.get_font_path(font);
    if (font_path.empty()) return nullptr;
    return Emboss::load_font(font_path.c_str());
#else
    // HERE is place to add implementation for another platform
    // to convert wxFont to font data as windows or font file path as linux
    return nullptr;
#endif
}

FontItem::Type WxFontUtils::get_actual_type()
{
#ifdef _WIN32
    return FontItem::Type::wx_win_font_descr;
#elif defined(__APPLE__)
    return FontItem::Type::wx_mac_font_descr;
#elif defined(__linux__)
    return FontItem::Type::wx_lin_font_descr;
#else
    return FontItem::Type::undefined;
#endif
}

FontItem WxFontUtils::get_font_item(const wxFont &font)
{
    std::string    name     = get_human_readable_name(font);
    std::string    fontDesc = store_wxFont(font);
    FontItem::Type type     = get_actual_type();

    // synchronize font property with actual font
    FontProp font_prop;    
    WxFontUtils::update_property(font_prop, font);
    return FontItem(name, fontDesc, type, font_prop);
}

FontItem WxFontUtils::get_os_font()
{
    wxSystemFont system_font = wxSYS_DEFAULT_GUI_FONT;
    wxFont       font        = wxSystemSettings::GetFont(system_font);
    FontItem     fi          = get_font_item(font);
    fi.name += std::string(" (OS default)");
    return get_font_item(font);
}

std::string WxFontUtils::get_human_readable_name(const wxFont &font)
{
    if (!font.IsOk()) return "Font is NOT ok.";
    // Face name is optional in wxFont
    if (!font.GetFaceName().empty()) {
        return std::string(font.GetFaceName().c_str());
    } else {
        return std::string((font.GetFamilyString() + " " +
                            font.GetStyleString() + " " +
                            font.GetWeightString())
                               .c_str());
    }
}

std::string WxFontUtils::store_wxFont(const wxFont &font)
{
    // wxString os = wxPlatformInfo::Get().GetOperatingSystemIdName();
    wxString font_descriptor = font.GetNativeFontInfoDesc();
    return std::string(font_descriptor.c_str());
}

std::optional<wxFont> WxFontUtils::load_wxFont(
    const std::string &font_descriptor)
{
    wxString font_descriptor_wx(font_descriptor);
    wxFont   wx_font(font_descriptor_wx);
    if (!wx_font.IsOk()) return {};
    return wx_font;
}

const std::map<wxFontFamily, std::string> WxFontUtils::from_family(
    {{wxFONTFAMILY_DEFAULT, "default"},
     {wxFONTFAMILY_DECORATIVE, "decorative"},
     {wxFONTFAMILY_ROMAN, "roman"},
     {wxFONTFAMILY_SCRIPT, "script"},
     {wxFONTFAMILY_SWISS, "swiss"},
     {wxFONTFAMILY_MODERN, "modern"},
     {wxFONTFAMILY_TELETYPE, "teletype"},
     {wxFONTFAMILY_MAX, "max"},
     {wxFONTFAMILY_UNKNOWN, "unknown"}});
const std::map<std::string, wxFontFamily> WxFontUtils::to_family =
    MapUtils::create_oposit(WxFontUtils::from_family);

const std::map<wxFontStyle, std::string> WxFontUtils::from_style(
    {{wxFONTSTYLE_ITALIC, "italic"},
     {wxFONTSTYLE_SLANT, "slant"},
     {wxFONTSTYLE_NORMAL, "normal"}});
const std::map<std::string, wxFontStyle> WxFontUtils::to_style =
    MapUtils::create_oposit(WxFontUtils::from_style);

const std::map<wxFontWeight, std::string> WxFontUtils::from_weight(
    {{wxFONTWEIGHT_THIN, "thin"},
     {wxFONTWEIGHT_EXTRALIGHT, "extraLight"},
     {wxFONTWEIGHT_LIGHT, "light"},
     {wxFONTWEIGHT_NORMAL, "normal"},
     {wxFONTWEIGHT_MEDIUM, "medium"},
     {wxFONTWEIGHT_SEMIBOLD, "semibold"},
     {wxFONTWEIGHT_BOLD, "bold"},
     {wxFONTWEIGHT_EXTRABOLD, "extraBold"},
     {wxFONTWEIGHT_HEAVY, "heavy"},
     {wxFONTWEIGHT_EXTRAHEAVY, "extraHeavy"}});
const std::map<std::string, wxFontWeight> WxFontUtils::to_weight =
    MapUtils::create_oposit(WxFontUtils::from_weight);

std::optional<wxFont> WxFontUtils::create_wxFont(const FontItem &fi,
                                                 const FontProp &fp)
{
    double     point_size = static_cast<double>(fp.size_in_mm);
    wxFontInfo info(point_size);
    if (fp.family.has_value()) {
        auto it = to_family.find(*fp.style);
        if (it != to_family.end()) info.Family(it->second);
    }
    if (fp.face_name.has_value()) {
        wxString face_name(*fp.face_name);
        info.FaceName(face_name);
    }
    if (fp.style.has_value()) {
        auto it = to_style.find(*fp.style);
        if (it != to_style.end()) info.Style(it->second);
    }
    if (fp.weight.has_value()) {
        auto it = to_weight.find(*fp.weight);
        if (it != to_weight.end()) info.Weight(it->second);
    }

    // Improve: load descriptor instead of store to font property to 3mf
    // switch (fi.type) {
    // case FontItem::Type::wx_lin_font_descr:
    // case FontItem::Type::wx_win_font_descr:
    // case FontItem::Type::wx_mac_font_descr:
    // case FontItem::Type::file_path:
    // case FontItem::Type::undefined:
    // default:
    //}

    wxFont font(info);
    if (!font.IsOk()) return {};
    return font;
}

void WxFontUtils::update_property(FontProp &font_prop, const wxFont &font)
{
    // The point size is defined as 1/72 of the Anglo-Saxon inch (25.4 mm): it
    // is approximately 0.0139 inch or 352.8 um. But it is too small, so I
    // decide use point size as mm for emboss
    font_prop.size_in_mm = font.GetPointSize(); // *0.3528f;
    font_prop.emboss     = font_prop.size_in_mm / 2.f;

    wxString    wx_face_name = font.GetFaceName();
    std::string face_name((const char *) wx_face_name.ToUTF8());
    if (!face_name.empty()) font_prop.face_name = face_name;

    wxFontFamily wx_family = font.GetFamily();
    if (wx_family != wxFONTFAMILY_DEFAULT) {
        auto it = from_family.find(wx_family);
        if (it != from_family.end()) font_prop.family = it->second;
    }

    wxFontStyle wx_style = font.GetStyle();
    if (wx_style != wxFONTSTYLE_NORMAL) {
        auto it = from_style.find(wx_style);
        if (it != from_style.end()) font_prop.style = it->second;
    }

    wxFontWeight wx_weight = font.GetWeight();
    if (wx_weight != wxFONTWEIGHT_NORMAL) {
        auto it = from_weight.find(wx_weight);
        if (it != from_weight.end()) font_prop.weight = it->second;
    }
}

bool WxFontUtils::is_italic(const wxFont &font) {
    wxFontStyle wx_style = font.GetStyle();
    return wx_style == wxFONTSTYLE_ITALIC || 
        wx_style == wxFONTSTYLE_SLANT;
}

bool WxFontUtils::is_bold(const wxFont &font) {
    wxFontWeight wx_weight = font.GetWeight();
    return wx_weight != wxFONTWEIGHT_NORMAL;
}

std::unique_ptr<Emboss::FontFile> WxFontUtils::set_italic(wxFont &font, const Emboss::FontFile &font_file)
{
    static std::vector<wxFontStyle> italic_styles = {
        wxFontStyle::wxFONTSTYLE_ITALIC,
        wxFontStyle::wxFONTSTYLE_SLANT
    };

    for (wxFontStyle style : italic_styles) { 
        font.SetStyle(style);
        std::unique_ptr<Emboss::FontFile> new_font_file =
            WxFontUtils::create_font_file(font);
        
        // can create italic font?
        if (new_font_file == nullptr) continue;

        // is still same font file pointer?
        if (font_file == *new_font_file) continue;

        return new_font_file;
    }
    return nullptr;
}

std::unique_ptr<Emboss::FontFile> WxFontUtils::set_bold(wxFont &font, const Emboss::FontFile& font_file)
{
    static std::vector<wxFontWeight> bold_weight = {
        wxFontWeight::wxFONTWEIGHT_BOLD,
        wxFontWeight::wxFONTWEIGHT_HEAVY,
        wxFontWeight::wxFONTWEIGHT_EXTRABOLD,
        wxFontWeight::wxFONTWEIGHT_EXTRAHEAVY
    };

    for (wxFontWeight weight : bold_weight) { 
        font.SetWeight(weight);
        std::unique_ptr<Emboss::FontFile> new_font_file =
            WxFontUtils::create_font_file(font);

        // can create bold font file?
        if (new_font_file == nullptr) continue;

        // is still same font file pointer?
        if (font_file == *new_font_file) continue;

        return new_font_file;
    }
    return nullptr;
}