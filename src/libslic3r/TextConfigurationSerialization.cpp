#include "TextConfigurationSerialization.hpp"

#include <cereal/archives/xml.hpp>
#include <strstream>

using namespace Slic3r;

const std::string TextConfigurationSerialization::font_item = "FontItem";
const std::string TextConfigurationSerialization::font_prop = "FontProp";
const std::string TextConfigurationSerialization::text = "Text";

namespace cereal {
template<class Archive>
void serialize(Archive &archive, FontItem &font_item)
{
    archive(
        CEREAL_NVP(font_item.name), 
        CEREAL_NVP(font_item.path)
        //,CEREAL_NVP(font_item.type)
    );
}

template<class Archive>
void serialize(Archive &archive, FontProp &font_prop)
{
    archive(
        CEREAL_NVP(font_prop.char_gap),
        CEREAL_NVP(font_prop.line_gap),
        CEREAL_NVP(font_prop.flatness),
        CEREAL_NVP(font_prop.size_in_mm),
        CEREAL_NVP(font_prop.emboss)
    );
}

template<class Archive>
void serialize(Archive &archive, TextConfiguration &text_configuration)
{
    archive(CEREAL_NVP(text_configuration.font_item),
            CEREAL_NVP(text_configuration.font_prop),
            CEREAL_NVP(text_configuration.text));
}
} // namespace cereal

std::string TextConfigurationSerialization::to_string(
    const TextConfiguration &text_configuration)
{
    std::strstream ss;
    {
        cereal::XMLOutputArchive archive(ss);
        // CEREAL_NVP - Names the output the same as the variable name
        archive(CEREAL_NVP(text_configuration));
    }
    std::string result = ss.str();    
    static size_t start  = std::string("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<cereal>\n\t").size();
    size_t end = result.find("\n</cereal>", start);
    result = result.substr(start, end - start);
    return result;
}

std::optional<TextConfiguration> TextConfigurationSerialization::from_string(const std::string &data)
{
    std::strstream ss;
    ss << data;
    cereal::XMLInputArchive archive(ss);
    TextConfiguration       tc;
    archive(tc);
    return tc;
}