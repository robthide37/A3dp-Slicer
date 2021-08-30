#ifndef slic3r_Emboss_hpp_
#define slic3r_Emboss_hpp_

#include <vector>
#include <optional>
#include <admesh/stl.h> // indexed_triangle_set
#include "Polygon.hpp"

namespace Slic3r {

/// <summary>
/// class with only static function add ability to engraved OR raised
/// text OR polygons onto model surface
/// </summary>
class Emboss
{
public:
    Emboss() = delete;

    /// <summary>
    /// keep information from file about font
    /// </summary>
    struct Font
    {
        // loaded data from font file
        std::vector<unsigned char> buffer;

        unsigned int index=0; // index of actual file info in collection
        unsigned int count=0; // count of fonts in file collection

        // vertical position is "scale*(ascent - descent + lineGap)"
        int ascent=0, descent=0, linegap=0;
        // user defined unscaled char space
        int extra_char_space = 0;

        // unscaled precision of lettter outline curve in conversion to lines
        float flatness = 2.;

        // change size of font
        float scale = 1.;
    };

    /// <summary>
    /// Load font file into buffer
    /// </summary>
    /// <param name="file_path">Location of .ttf or .ttc font file</param>
    /// <returns>Font object when loaded.</returns>
    static std::optional<Font> load_font(const char *file_path);

    /// <summary>
    /// convert letter into polygons
    /// </summary>
    /// <param name="font">Define fonts</param>
    /// <param name="letter">Character to convert</param>
    /// <returns>inner polygon ccw(outer cw)</returns>
    static Polygons letter2polygons(const Font &font, char letter);    

    static Polygons text2polygons(const Font &font, const std::string &text);

    /// <summary>
    /// Create triangle model for text
    /// </summary>
    /// <param name="text"></param>
    /// <param name="font"></param>
    /// <param name="z_size"></param>
    /// <returns></returns>
    static indexed_triangle_set create_model(const std::string &text,
                                             const Font &       font,
                                             float              z_size);
};

} // namespace Slic3r
#endif // slic3r_Emboss_hpp_
