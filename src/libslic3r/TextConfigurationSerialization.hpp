#ifndef slic3r_TextConfigurationSerialization_hpp_
#define slic3r_TextConfigurationSerialization_hpp_

#include "TextConfiguration.hpp"
#include <optional>

namespace Slic3r {

// utility for convert TextConfiguration into string and vice versa
class TextConfigurationSerialization
{
public:
    TextConfigurationSerialization() = delete; // only static functions
    static std::string to_string(const TextConfiguration &text_configuration);
    static std::optional<TextConfiguration> from_string(const std::string &data);

private:
    static const std::string font_item;
    static const std::string font_prop;
    static const std::string text;
};
} // namespace Slic3r

#endif // slic3r_TextConfigurationSerialization_hpp_
