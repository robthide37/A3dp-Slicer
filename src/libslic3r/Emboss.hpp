#ifndef slic3r_Emboss_hpp_
#define slic3r_Emboss_hpp_

#include <vector>
#include <set>
#include <optional>
#include <memory>
#include <admesh/stl.h> // indexed_triangle_set
#include "Polygon.hpp"
#include "TextConfiguration.hpp"

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
    /// Collect fonts registred inside OS
    /// </summary>
    /// <returns>OS registred TTF font files(full path) with names</returns>
    static FontList get_font_list();
#ifdef _WIN32
    static FontList get_font_list_by_register();
    static FontList get_font_list_by_enumeration();
    static FontList get_font_list_by_folder();
#endif

    /// <summary>
    /// OS dependent function to get location of font by its name descriptor
    /// </summary>
    /// <param name="font_face_name">Unique identificator for font</param>
    /// <returns>File path to font when found</returns>
    static std::optional<std::wstring> get_font_path(const std::wstring &font_face_name);

    // description of one letter
    struct Glyph
    {
        ExPolygons shape;
        int advance_width, left_side_bearing;
    };
    // cache for glyph by unicode
    using Glyphs = std::map<int, Glyph>;
        
    /// <summary>
    /// keep information from file about font
    /// </summary>
    struct Font
    {
        // loaded data from font file
        std::vector<unsigned char> buffer;

        unsigned int index = 0; // index of actual file info in collection
        unsigned int count = 0; // count of fonts in file collection

        // vertical position is "scale*(ascent - descent + lineGap)"
        int            ascent = 0, descent = 0, linegap = 0;

        Emboss::Glyphs cache; // cache of glyphs
    };

    /// <summary>
    /// Load font file into buffer
    /// </summary>
    /// <param name="file_path">Location of .ttf or .ttc font file</param>
    /// <returns>Font object when loaded.</returns>
    static std::optional<Font> load_font(const char *file_path);
    // data = raw file data
    static std::optional<Font> load_font(std::vector<unsigned char> data);
#ifdef _WIN32
    // fix for unknown pointer HFONT
    using HFONT = void*;
    static std::optional<Font> load_font(HFONT hfont);
#endif // _WIN32

    /// <summary>
    /// convert letter into polygons
    /// </summary>
    /// <param name="font">Define fonts</param>
    /// <param name="letter">One character defined by unicode codepoint</param>
    /// <param name="flatness">Precision of lettter outline curve in conversion to lines</param>
    /// <returns>inner polygon cw(outer ccw)</returns>
    static std::optional<Glyph> letter2glyph(const Font &font, int letter, float flatness);

    /// <summary>
    /// Convert text into polygons
    /// </summary>
    /// <param name="font">Define fonts + cache, which could extend</param>
    /// <param name="text">Characters to convert</param>
    /// <param name="font_prop">User defined property of the font</param>
    /// <returns>Inner polygon cw(outer ccw)</returns>
    static ExPolygons text2shapes(Font &          font,
                                  const char *    text,
                                  const FontProp &font_prop);

    /// <summary>
    /// Project 2d point into space
    /// Could be plane, sphere, cylindric, ...
    /// </summary>
    class IProject
    {
    public:
        virtual ~IProject() = default;
        /// <summary>
        /// convert 2d point to 3d point
        /// </summary>
        /// <param name="p">2d coordinate</param>
        /// <returns>
        /// first - front spatial point
        /// second - back spatial point
        /// </returns>
        virtual std::pair<Vec3f, Vec3f> project(const Point &p) const = 0;
    };

    /// <summary>
    /// Create triangle model for text
    /// </summary>
    /// <param name="shape2d">text or image</param>
    /// <param name="projection">Define transformation from 2d to 3d(orientation, position, scale, ...)</param>
    /// <returns>Projected shape into space</returns>
    static indexed_triangle_set polygons2model(const ExPolygons &shape2d, const IProject& projection);

    class ProjectZ : public IProject
    {
    public:
        ProjectZ(float depth) : m_depth(depth) {}
        // Inherited via IProject
        virtual std::pair<Vec3f, Vec3f> project(const Point &p) const override;
        float m_depth;
    };

    class ProjectScale : public IProject
    {
        std::unique_ptr<IProject> core;
    public:
        ProjectScale(std::unique_ptr<IProject> core, float scale)
            : m_scale(scale)
            , core(std::move(core))
        {}

        // Inherited via IProject
        virtual std::pair<Vec3f, Vec3f> project(const Point &p) const override
        {
            auto res = core->project(p);
            return std::make_pair(res.first * m_scale, res.second * m_scale);
        }

        float m_scale;
    };
};

} // namespace Slic3r
#endif // slic3r_Emboss_hpp_
