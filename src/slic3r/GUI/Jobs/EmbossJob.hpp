#ifndef slic3r_EmbossJob_hpp_
#define slic3r_EmbossJob_hpp_

#include "StopableJob.hpp"
#include "libslic3r/Emboss.hpp"

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
    // when volume_ptr == nullptr than new volume will be created
    ModelVolume *volume_ptr;
    // when volume_ptr == nullptr && object_idx < 0 than new object will be created
    int object_idx;
};

class EmbossJob : public StopableJob<EmbossData>
{
protected:
    void process(std::unique_ptr<EmbossData> input, StopCondition is_stop) override;
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
