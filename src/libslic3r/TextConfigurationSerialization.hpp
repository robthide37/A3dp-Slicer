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
    static std::string serialize(const TextConfiguration &text_configuration);
    static std::optional<TextConfiguration> deserialize(const std::string &data);

    // convert type to string and vice versa
    static const std::map<std::string, FontItem::Type> to_type;
    static const std::map<FontItem::Type, std::string> to_string;
    
    static const char separator;

    // Move to map utils
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
