#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include "StopableJob.hpp"
#include "libslic3r/Emboss.hpp"
//#include "libslic3r/ObjectID.hpp"

namespace Slic3r {
class ModelVolume;
}

namespace Slic3r::GUI {

struct EmbossData
{
    // Pointer on Data of font (glyph shapes)
    std::shared_ptr<Emboss::Font> font;
    // font item is not used for create object
    TextConfiguration text_configuration;
    // new volume name created from text
    std::string volume_name;

    // unique identifier of volume to change
    // I can't proove of alive pointer
    ModelVolume *volume;

    // unique identifier of volume to change
    // Change of volume change id, last change could disapear
    //ObjectID     volume_id;

    EmbossData(std::shared_ptr<Emboss::Font> font,
               TextConfiguration             text_configuration,
               std::string                   volume_name,
               ModelVolume *                 volume)
        : font(std::move(font))
        , text_configuration(text_configuration)
        , volume_name(volume_name)
        , volume(volume)
    {}
};

class EmbossJob : public StopableJob<EmbossData>
{
protected:
    void process(std::unique_ptr<EmbossData> input) override;
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
