#ifndef slic3r_CreateFontStyleImagesJob_hpp_
#define slic3r_CreateFontStyleImagesJob_hpp_

#include <vector>
#include <string>
#include <libslic3r/Emboss.hpp>
#include "slic3r/Utils/FontManager.hpp"
#include "Job.hpp"

namespace Slic3r::GUI {

/// <summary>
/// Create texture with name of styles written by its style
/// NOTE: Access to glyph cache is possible only from job
/// </summary>
class CreateFontStyleImagesJob : public Job
{
    EmbossStyleManager::StyleImagesData m_input;

    // Output data
    // texture size
    int width, height;
    // texture data
    std::vector<unsigned char> pixels; 
    // descriptors of sub textures
    std::vector<EmbossStyleManager::StyleImage> images;

public:
    CreateFontStyleImagesJob(EmbossStyleManager::StyleImagesData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};

} // namespace Slic3r::GUI

#endif // slic3r_CreateFontStyleImagesJob_hpp_
