#include "TextConfigurationSerialization.hpp"
#include "Utils.hpp"

using namespace Slic3r;
// Convert map
const std::map<FontItem::Type, std::string> TextConfigurationSerialization::to_string = {
    {FontItem::Type::file_path, "file_path"},
    {FontItem::Type::wx_font_descr, "wx_font_descriptor"},
    {FontItem::Type::undefined, "unknown"}
};

const std::map<std::string, FontItem::Type> TextConfigurationSerialization::to_type =
    TextConfigurationSerialization::create_oposit_map(TextConfigurationSerialization::to_string);

const char TextConfigurationSerialization::separator = '|';

std::string TextConfigurationSerialization::serialize(const TextConfiguration &text_configuration)
{
    size_t size = 1 + 3 + 5;
    std::vector<std::string> columns;
    columns.reserve(size);
    columns.emplace_back(text_configuration.text);
    to_columns(text_configuration.font_item, columns);
    to_columns(text_configuration.font_prop, columns);
    assert(columns.size() == size);
    return serialize(columns, separator);
}

std::optional<TextConfiguration> TextConfigurationSerialization::deserialize(const std::string &data)
{
    size_t size = 1 + 3 + 5;
    std::vector<std::string> columns;
    columns.reserve(size);
    deserialize(data, separator, columns);
    assert(columns.size() == size);

    const std::string& text = columns[0];
    std::optional<FontItem> font_item = get_font_item(columns, 1);
    if (!font_item.has_value()) return {};
    std::optional<FontProp> font_prop = get_font_prop(columns, 4);
    if (!font_prop.has_value()) return {};

    return TextConfiguration(*font_item, *font_prop, text);
}

std::string TextConfigurationSerialization::serialize(const FontList &font_list)
{
    std::vector<std::string> columns;
    columns.reserve(3 * font_list.size());
    for (const FontItem &fi : font_list) to_columns(fi, columns);
    return serialize(columns, separator);
}

FontList TextConfigurationSerialization::deserialize_font_list(const std::string &data)
{
    std::vector<std::string> columns;
    deserialize(data, separator, columns);
    if ((columns.size() % 3) != 0) return {};
    size_t   count = columns.size() / 3;
    FontList fl;
    fl.reserve(count);
    for (size_t i = 0; i < count; i++) 
    {
        std::optional<FontItem> fi = get_font_item(columns, i * 3);
        if (!fi.has_value()) return {};
        fl.emplace_back(*fi);
    }
    return fl;
}

std::string TextConfigurationSerialization::serialize(const std::vector<std::string> &columns, char separator)
{
    std::string result;
    const std::string separator_str = std::string(" ") + separator + ' ';
    bool is_first = true;
    for (const std::string& column : columns) {
        if (is_first) is_first = false;
        else result += separator_str;
        result += twice(column, separator);
    }
    return result;
}

void TextConfigurationSerialization::deserialize(const std::string &data, char separator, std::vector<std::string> &columns)
{
    if (data.empty()) return;

    size_t position = 0;
    while (position < data.size()) {
        size_t start = position;
        position     = find_odd(data, position+1, separator);
        size_t size  = position - start;

        // is not last column
        if (position != data.size()) {
            assert(size != 0);
            --size; 
        }
        // is not first column
        if (start != 0) {
            // previous separator + space = 2
            start+=2;
            assert(size >=2);
            size-=2;
        }

        std::string column = data.substr(start, size);

        // heal column 
        columns.emplace_back(reduce(column, separator));
    }
}

std::string TextConfigurationSerialization::serialize(const FontItem::Type &type)
{
    auto it = to_string.find(type);
    assert(it != to_string.end());
    if (it == to_string.end()) 
        return serialize(FontItem::Type::undefined);
    return it->second;
}

FontItem::Type TextConfigurationSerialization::deserialize_type(const std::string &type)
{
    auto it = to_type.find(type);
    assert(it != to_type.end());
    if (it == to_type.end()) return FontItem::Type::undefined;
    return it->second;
}

void TextConfigurationSerialization::to_columns(
    const FontItem &font_item, std::vector<std::string> &columns)
{
    columns.emplace_back(font_item.name);
    columns.emplace_back(font_item.path);
    columns.emplace_back(serialize(font_item.type));
}

std::optional<FontItem> TextConfigurationSerialization::get_font_item(
    const std::vector<std::string> &columns, size_t offset)
{
    if (columns.size() <= (offset + 2)) return {}; // no enough columns
    FontItem::Type type = deserialize_type(columns[offset + 2]);
    return FontItem(columns[offset], columns[offset + 1], type);
}

void TextConfigurationSerialization::to_columns(
    const FontProp& font_prop, std::vector<std::string> &columns)
{
    columns.emplace_back(std::to_string(font_prop.emboss));
    columns.emplace_back(std::to_string(font_prop.flatness));
    columns.emplace_back(std::to_string(font_prop.size_in_mm));
    columns.emplace_back(std::to_string(font_prop.char_gap));
    columns.emplace_back(std::to_string(font_prop.line_gap));
}

std::optional<FontProp> TextConfigurationSerialization::get_font_prop(
    const std::vector<std::string> &columns, size_t offset)
{
    if (columns.size() <= (offset + 4)) return {}; // no enough columns
    FontProp font_prop;
    font_prop.emboss     = static_cast<float>(std::atof(columns[offset].c_str()));
    font_prop.flatness   = static_cast<float>(std::atof(columns[offset+1].c_str()));
    font_prop.size_in_mm = static_cast<float>(std::atof(columns[offset+2].c_str()));
    font_prop.char_gap   = static_cast<float>(std::atof(columns[offset+3].c_str()));
    font_prop.line_gap   = std::atoi(columns[offset+4].c_str());
    return font_prop;
}