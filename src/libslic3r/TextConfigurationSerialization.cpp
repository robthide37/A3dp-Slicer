#include "TextConfigurationSerialization.hpp"
#include "utils.hpp"

using namespace Slic3r;
// Convert map
const std::map<FontItem::Type, std::string> TextConfigurationSerialization::to_string = {
    {FontItem::Type::file_path, "file_path"},
    {FontItem::Type::wx_font_descr, "wx_font_descriptor"}
};

const std::map<std::string, FontItem::Type> TextConfigurationSerialization::to_type =
    TextConfigurationSerialization::create_oposit_map(TextConfigurationSerialization::to_string);

const char TextConfigurationSerialization::separator = '|';

std::string TextConfigurationSerialization::serialize(const TextConfiguration &text_configuration)
{
    // IMPROVE: make general and move to string utils
    auto twice_separator = [](const std::string& data) {
        // no value is one space
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

std::optional<TextConfiguration> TextConfigurationSerialization::deserialize(const std::string &data)
{
    // IMPROVE: make general and move to string utils
    auto find_separator = [&data](std::string::size_type pos) {
        pos = data.find(separator, pos);
        while (pos != data.npos && data[pos + 1] == separator)
            pos = data.find(separator, pos + 2);
        return pos;
    };
    // IMPROVE: make general and move to string utils
    auto reduce_separator = [](const std::string& item) {
        std::string::size_type pos = item.find(separator);
        if (pos == item.npos) return item;
        std::string copy = item;
        do {
            assert(copy[pos] == separator);
            assert(copy[pos+1] == separator);
            copy.erase(pos, size_t(1));
            pos = copy.find(separator, pos + 1);
        } while (pos != copy.npos);
        return copy;
    };

    std::string::size_type start = 0;
    std::string::size_type size  = find_separator(start);
    auto reduce_and_move = [&data, &start, &size, &reduce_separator, &find_separator]() {
        if (start == data.npos) return std::string();
        std::string res = reduce_separator(data.substr(start, size));
        start = size + 1;
        size  = find_separator(start) - start;
        return res;
    };
    auto get_float_and_move = [&data, &start, &size, &find_separator]() {
        if (start == data.npos) return 0.f;
        float res = std::atof(data.substr(start, size).c_str());
        start     = size + 1;
        size      = find_separator(start) - start;
        return res;
    };
    auto get_int_and_move = [&data, &start, &size, &find_separator]() {
        if (start == data.npos) return 0;
        int res = std::atoi(data.substr(start, size).c_str());
        start   = size + 1;
        size    = find_separator(start) - start;
        return res;
    };

    std::string text = reduce_and_move();
    std::string name = reduce_and_move();
    std::string path = reduce_and_move();
    std::string type = reduce_and_move();
    auto it = to_type.find(type);
    if (it == to_type.end()) return {}; // no valid type
    FontItem font_item(name,path,it->second);

    FontProp font_prop;
    font_prop.emboss     = get_float_and_move();
    font_prop.flatness   = get_float_and_move();
    font_prop.size_in_mm = get_float_and_move();
    font_prop.char_gap   = get_int_and_move();
    if (start == data.npos) return {}; // no valid data
    font_prop.line_gap   = get_int_and_move();
    return TextConfiguration(font_item, font_prop, text);
}