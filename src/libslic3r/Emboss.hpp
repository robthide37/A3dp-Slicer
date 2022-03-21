#ifndef slic3r_Emboss_hpp_
#define slic3r_Emboss_hpp_

#include <vector>
#include <set>
#include <optional>
#include <memory>
#include <admesh/stl.h> // indexed_triangle_set
#include "Polygon.hpp"
#include "ExPolygon.hpp"
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
    
    // every glyph's shape point is divided by SHAPE_SCALE - increase precission of fixed point value
    static double SHAPE_SCALE;    

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
        // NOTE: shape is scaled by SHAPE_SCALE 
        // to be able store points without floating points
        ExPolygons shape;

        // values are in font points
        int advance_width=0, left_side_bearing=0;
    };
    // cache for glyph by unicode
    using Glyphs = std::map<int, Glyph>;
        
    /// <summary>
    /// keep information from file about font 
    /// (store file data itself)
    /// + cache data readed from buffer
    /// </summary>
    struct FontFile
    {
        // loaded data from font file
        // must store data size for imgui rasterization
        // To not store data on heap and To prevent unneccesary copy
        // data are stored inside unique_ptr
        std::unique_ptr<std::vector<unsigned char>> data;

        // count of fonts when data are collection of fonts
        unsigned int count; 

        // vertical position is "scale*(ascent - descent + lineGap)"
        int ascent, descent, linegap;

        // for convert font units to pixel
        int unit_per_em;

        FontFile(std::unique_ptr<std::vector<unsigned char>> data,
                 unsigned int                     count,
                 int                              ascent,
                 int                              descent,
                 int                              linegap,
                 int                              unit_per_em)
            : data(std::move(data))
            , count(count)
            , ascent(ascent)
            , descent(descent)
            , linegap(linegap)
            , unit_per_em(unit_per_em)
        {
            assert(this->data != nullptr);
        }
        bool operator==(const FontFile &other) const {
            return count == other.count && ascent == other.ascent &&
                   descent == other.descent && linegap == other.linegap &&
                   data->size() == other.data->size();
                //&& *data == *other.data;            
        }
    };

    /// <summary>
    /// Add caching for shape of glyphs
    /// </summary>
    struct FontFileWithCache
    {
        std::shared_ptr<const FontFile> font_file;
        // cache for glyph shape
        std::shared_ptr<Emboss::Glyphs> cache;

        FontFileWithCache() : font_file(nullptr), cache(nullptr) {}
        FontFileWithCache(std::unique_ptr<FontFile> font_file)
            : font_file(std::move(font_file))
            , cache(std::make_shared<Emboss::Glyphs>())
        {}
        bool has_value() const { return font_file != nullptr && cache != nullptr; }
    };

    /// <summary>
    /// Load font file into buffer
    /// </summary>
    /// <param name="file_path">Location of .ttf or .ttc font file</param>
    /// <returns>Font object when loaded.</returns>
    static std::unique_ptr<FontFile> create_font_file(const char *file_path);
    // data = raw file data
    static std::unique_ptr<FontFile> create_font_file(std::unique_ptr<std::vector<unsigned char>> data);
#ifdef _WIN32
    // fix for unknown pointer HFONT
    using HFONT = void*;
    static void * can_load(HFONT hfont);
    static std::unique_ptr<FontFile> create_font_file(HFONT hfont);
#endif // _WIN32

    /// <summary>
    /// convert letter into polygons
    /// </summary>
    /// <param name="font">Define fonts</param>
    /// <param name="letter">One character defined by unicode codepoint</param>
    /// <param name="flatness">Precision of lettter outline curve in conversion to lines</param>
    /// <returns>inner polygon cw(outer ccw)</returns>
    static std::optional<Glyph> letter2glyph(const FontFile &font, int letter, float flatness);

    /// <summary>
    /// Convert text into polygons
    /// </summary>
    /// <param name="font">Define fonts + cache, which could extend</param>
    /// <param name="text">Characters to convert</param>
    /// <param name="font_prop">User defined property of the font</param>
    /// <returns>Inner polygon cw(outer ccw)</returns>
    static ExPolygons text2shapes(FontFileWithCache &font,
                                  const char *    text,
                                  const FontProp &font_prop);

    /// <summary>
    /// Use data from font property to modify transformation
    /// </summary>
    /// <param name="font_prop">Z-move as surface distance(FontProp::distance)
    /// Z-rotation as angle to Y axis(FontProp::angle)</param>
    /// <param name="transformation">In / Out transformation to modify by property</param>
    static void apply_transformation(const FontProp &font_prop,
                                     Transform3d    &transformation);

    /// <summary>
    /// Read information from naming table of font file
    /// search for italic (or oblique), bold italic (or bold oblique)
    /// </summary>
    /// <param name="font">Selector of font</param>
    /// <param name="font_index">Index of font in collection</param>
    /// <returns>True when the font description contains italic/obligue otherwise False</returns>
    static bool is_italic(const FontFile &font, unsigned int font_index);

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

        /// <summary>
        /// Move point with respect to projection direction
        /// e.g. Orthogonal projection will move with point by direction
        /// e.g. Spherical projection need to use center of projection
        /// </summary>
        /// <param name="point">Spatial point coordinate</param>
        /// <returns>Projected spatial point</returns>
        virtual Vec3f project(const Vec3f &point) const = 0;
    };

    /// <summary>
    /// Create triangle model for text
    /// </summary>
    /// <param name="shape2d">text or image</param>
    /// <param name="projection">Define transformation from 2d to 3d(orientation, position, scale, ...)</param>
    /// <returns>Projected shape into space</returns>
    static indexed_triangle_set polygons2model(const ExPolygons &shape2d, const IProject& projection);
        
    /// <summary>
    /// Create transformation for emboss text object to lay on surface point
    /// </summary>
    /// <param name="position">Position of surface point</param>
    /// <param name="normal">Normal of surface point</param>
    /// <param name="up_limit">Is compared with normal.z to suggest up direction</param>
    /// <returns>Transformation onto surface point</returns>
    static Transform3d create_transformation_onto_surface(
        const Vec3f &position, const Vec3f &normal, float up_limit = 0.9f);

    class ProjectZ : public IProject
    {
    public:
        ProjectZ(float depth) : m_depth(depth) {}
        // Inherited via IProject
        std::pair<Vec3f, Vec3f> project(const Point &p) const override;
        Vec3f project(const Vec3f &point) const override;
        float m_depth;
    };

    class ProjectScale : public IProject
    {
        std::unique_ptr<IProject> core;
        float m_scale;
    public:
        ProjectScale(std::unique_ptr<IProject> core, float scale)
            : core(std::move(core)), m_scale(scale)
        {}

        // Inherited via IProject
        std::pair<Vec3f, Vec3f> project(const Point &p) const override
        {
            auto res = core->project(p);
            return std::make_pair(res.first * m_scale, res.second * m_scale);
        }
        Vec3f project(const Vec3f &point) const override{
            return core->project(point);
        }

    };
};

} // namespace Slic3r
#endif // slic3r_Emboss_hpp_
