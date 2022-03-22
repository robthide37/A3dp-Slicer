#ifndef slic3r_CreateFontStyleImagesJob_hpp_
#define slic3r_CreateFontStyleImagesJob_hpp_

#include <vector>
#include <string>
#include <libslic3r/Emboss.hpp>
#include "slic3r/Utils/FontManager.hpp"
#include "Job.hpp"

namespace Slic3r::GUI {

/// <summary>
/// Data needed to create Font Style Images
/// </summary>
struct StyleImagesData
{
    struct Item
    {
        Emboss::FontFileWithCache font;
        std::string               text;
        FontProp                  prop;
    };
    using Items = std::vector<Item>;

    // Keep styles to render
    Items styles;

    // maximal width in pixels of image
    int max_width;

    // is used in finalize to set result
    // and I Can't proof of alive
    FontManager *mng;
};

/// <summary>
/// Create texture with name of styles written by its style
/// NOTE: Access to glyph cache is possible only from job
/// </summary>
class CreateFontStyleImagesJob : public Job
{
    StyleImagesData m_input;

    // Output data
    // texture size
    int width, height;
    // texture data
    std::vector<unsigned char> pixels; 
    // descriptors of sub textures
    std::vector<FontManager::StyleImage> images;

public:
    CreateFontStyleImagesJob(StyleImagesData &&input);
    void process(Ctl &ctl) override;
    void finalize(bool canceled, std::exception_ptr &) override;
};

} // namespace Slic3r::GUI

#endif // slic3r_EmbossJob_hpp_
