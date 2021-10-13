#ifndef slic3r_TextConfigurationSerialization_hpp_
#define slic3r_TextConfigurationSerialization_hpp_

#include "TextConfiguration.hpp"
#include <optional>
#include <string_view>
#include <sstream>
#include <map>

namespace Slic3r {

// utility for convert TextConfiguration into string and vice versa
class TextConfigurationSerialization
{
public:
    TextConfigurationSerialization() = delete; // only static functions

    // store / load TextConfiguration - .3mf files
    static std::string serialize(const TextConfiguration &text_configuration);
    static std::optional<TextConfiguration> deserialize(const std::string &data);

    // store / load FontList - AppConfig
    static std::string serialize(const FontList &font_list);
    static FontList deserialize_font_list(const std::string &data);

    static std::string serialize(const FontItem::Type &type);
    static FontItem::Type deserialize_type(const std::string& type);

private:
    // convert type to string and vice versa
    static const std::map<std::string, FontItem::Type> to_type;
    static const std::map<FontItem::Type, std::string> to_string;
    static const char separator;
    
    // store / load general connection of string into one string
    static std::string serialize(const std::vector<std::string> &columns, char separator);
    static void deserialize(const std::string& data, char separator, std::vector<std::string> &columns);// columns vector should be reserved on valid count

    static void to_columns(const FontItem& font_item, std::vector<std::string> &columns);
    static std::optional<FontItem> get_font_item(const std::vector<std::string> &columns, size_t offset);

    static void to_columns(const FontProp& font_prop, std::vector<std::string> &columns);
    static std::optional<FontProp> get_font_prop(const std::vector<std::string> &columns, size_t offset);

    /// <summary>
    /// Twice all appearance of character in data.
    /// Twiced character could be used as separator for this data.
    /// IMPROVE: move to string utils
    /// </summary>
    /// <param name="data">input data to twice character</param>
    /// <param name="letter">Specify character to be twiced in data</param>
    /// <returns>String conatin only pair continous count of specified character</returns>
    static std::string twice(const std::string &data, char letter) 
    {
        // no value is one space
        if (data.empty()) return std::string(" ");
        std::string::size_type pos = data.find(letter);
        if (pos == data.npos) return data;
        // twice all separator inside data
        std::string copy = data; // copy
        do {
            copy.insert(pos, 1, letter);
            pos += 2;
        } while (copy.npos != (pos = copy.find(letter, pos)));
        return copy;
    };

    /// <summary>
    /// Reduce all twice appearance of character in data.
    /// Main purpose heal after twice function.
    /// IMPROVE: move to string utils
    /// </summary>
    /// <param name="data">input data to reduce character</param>
    /// <param name="letter">Specify character to be reduced</param>
    /// <returns>String conatining only half letter in row</returns>
    static std::string reduce(const std::string &data, char letter)
    {
        std::string::size_type pos = data.find(letter);
        if (pos == data.npos) return data;
        std::string copy = data; // copy
        do {
            std::string::size_type pos_plus_one = pos + 1;
            // is data endig with odd number of letters
            if (pos_plus_one == copy.npos) return copy;
            // is pair count of letter - should not appear after twice
            if (copy[pos_plus_one] != letter) continue;
            
            // reduce by removing first
            copy.erase(pos, size_t(1));
            pos = copy.find(letter, pos_plus_one);
        } while (pos != copy.npos);
        return copy;
    };

    /// <summary>
    /// Find odd position of letter in text data
    /// Used with combination of twice and reduce
    /// IMPROVE: move to string utils
    /// </summary>
    /// <param name="data">Input text</param>
    /// <param name="pos">Start index into data for searching odd letter</param>
    /// <param name="letter">Character to find</param>
    /// <returns>Index to data with next odd appearance of letter 
    /// OR size of data when NO next one exists</returns>
    static size_t find_odd(const std::string &data, size_t pos, char letter)
    {
        pos = data.find(letter, pos);
        while ((pos+1) <= data.size() && data[pos + 1] == letter)
            pos = data.find(letter, pos + 2);
        return pos;
    }

    /// <summary>
    /// Create map with swaped key-value
    /// IMPROVE: Move to map utils
    /// </summary>
    /// <param name="map">Input map</param>
    /// <returns>Map with swapped key-value</returns>
    template<typename Key, typename Value>
    static std::map<Value, Key> create_oposit_map(
        const std::map<Key, Value> &map)
    {
        std::map<Value, Key> result;
        for (const auto &it : map) result[it.second] = it.first;
        return result;
    }
};
} // namespace Slic3r

#endif // slic3r_TextConfigurationSerialization_hpp_
