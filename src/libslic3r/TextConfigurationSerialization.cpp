#include "TextConfigurationSerialization.hpp"
#include "utils.hpp"

using namespace Slic3r;
// Convert map
const std::map<FontItem::Type, std::string> TextConfigurationSerialization::to_string = {
    {FontItem::Type::file_path, "file_path"},
    {FontItem::Type::wx_font_descr, "wx_font_descriptor"}
};

const char TextConfigurationSerialization::separator = '|';

const std::map<std::string, FontItem::Type> TextConfigurationSerialization::to_type =
    TextConfigurationSerialization::create_oposit_map(TextConfigurationSerialization::to_string);

std::string TextConfigurationSerialization::serialize(const TextConfiguration &text_configuration)
{
    auto twice_separator = [](const std::string& data) {
        if (data.empty()) return std::string(" ");
        std::string::size_type pos = data.find(separator);
        if (pos == data.npos) return data;
        // twice all separator inside data
        std::string copy = data;
        do { 
            copy.insert(pos, 1, separator);
            pos += 2;
        } while (copy.npos != (pos = copy.find(separator, pos)));
        return copy;
    };

    const FontItem &font_item = text_configuration.font_item;
    const FontProp &font_prop = text_configuration.font_prop;
    return twice_separator(text_configuration.text) + separator +
           twice_separator(font_item.name) + separator +
           twice_separator(font_item.path) + separator +
           to_string.at(font_item.type) + separator +
           std::to_string(font_prop.emboss) + separator +
           std::to_string(font_prop.flatness) + separator +
           std::to_string(font_prop.size_in_mm) + separator +
           std::to_string(font_prop.char_gap) + separator +
           std::to_string(font_prop.line_gap);
}

void TextConfigurationSerialization::to_xml(
    const TextConfiguration &text_configuration,
    unsigned                 count_indent,
    std::stringstream &      stream)
{
    // Because of back compatibility must be main tag metadata (source: 3mf.cpp)
    const std::string_view main_tag     = "metadata";
    const std::string_view main_name_attr = "name";
    const std::string_view main_name_value= "TextConfiguration";

    const std::string_view font_item_tag = "fontItem";
    const std::string_view name_attr     = "name";
    const std::string_view path_attr     = "path";
    const std::string_view type_attr     = "type";

    const std::string_view font_prop_tag = "fontProp";
    const std::string_view char_gap_attr = "charGap";
    const std::string_view line_gap_attr = "lineGap";
    const std::string_view flatness_attr = "flatness";
    const std::string_view height_attr   = "height";
    const std::string_view depth_attr    = "depth";
        
    const std::string_view text_tag = "text";

    auto get_path = [&text_configuration]() {
        const std::string &path = text_configuration.font_item.path;   
        return xml_escape(
            (text_configuration.font_item.type == FontItem::Type::file_path) ?
                path.substr(path.find_last_of("/\\") + 1) : path);
    };

    std::string indent = std::string(count_indent, ' ');
    std::string indent2 = std::string(count_indent+1, ' '); 

    stream << indent << "<" << main_tag << " " << main_name_attr << "=\"" << main_name_value << "\">\n";
    stream << indent2 << "<" << font_item_tag 
        << ' ' << name_attr << "=\"" << xml_escape(text_configuration.font_item.name) << '"' 
        << ' ' << path_attr << "=\"" << get_path() << '"' 
        << ' ' << type_attr << "=\"" << to_string.at(text_configuration.font_item.type) << '"' 
        << "/>\n";
    stream << indent2 << "<" << font_prop_tag 
        << ' ' << char_gap_attr << "=\"" << text_configuration.font_prop.char_gap << '"' 
        << ' ' << line_gap_attr << "=\"" << text_configuration.font_prop.line_gap << '"' 
        << ' ' << flatness_attr << "=\"" << text_configuration.font_prop.flatness << '"' 
        << ' ' << height_attr << "=\"" << text_configuration.font_prop.size_in_mm << '"' 
        << ' ' << depth_attr << "=\"" << text_configuration.font_prop.emboss << '"' 
        << "/>\n";
    stream << indent2 << "<" << text_tag << ">";
    stream << xml_escape(text_configuration.text);
    stream << indent2 << "</" << text_tag << ">\n";
    stream << indent << "</" << main_tag << ">\n";
}

std::optional<TextConfiguration> TextConfigurationSerialization::from_string(const std::string &data)
{
    TextConfiguration tc;
    return tc;
}