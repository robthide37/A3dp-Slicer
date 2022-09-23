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
    static EmbossStyles get_font_list();
#ifdef _WIN32
    static EmbossStyles get_font_list_by_register();
    static EmbossStyles get_font_list_by_enumeration();
    static EmbossStyles get_font_list_by_folder();
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

        struct Info
        {
            // vertical position is "scale*(ascent - descent + lineGap)"
            int ascent, descent, linegap;

            // for convert font units to pixel
            int unit_per_em;
        };
        // info for each font in data
        std::vector<Info> infos;

        FontFile(std::unique_ptr<std::vector<unsigned char>> data,
                 std::vector<Info>                         &&infos)
            : data(std::move(data)), infos(std::move(infos))
        {
            assert(this->data != nullptr);
            assert(!this->data->empty());
        }

        bool operator==(const FontFile &other) const {
            if (data->size() != other.data->size())
                return false;
            //if(*data != *other.data) return false;
            for (size_t i = 0; i < infos.size(); i++) 
                if (infos[i].ascent != other.infos[i].ascent ||
                    infos[i].descent == other.infos[i].descent ||
                    infos[i].linegap == other.infos[i].linegap)
                    return false;
            return true;
        }
    };

    /// <summary>
    /// Add caching for shape of glyphs
    /// </summary>
    struct FontFileWithCache
    {
        // Pointer on data of the font file
        std::shared_ptr<const FontFile> font_file;

        // Cache for glyph shape
        // IMPORTANT: accessible only in plater job thread !!!
        // main thread only clear cache by set to another shared_ptr
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
    /// <param name="font_index">Index of font in collection</param>
    /// <param name="letter">One character defined by unicode codepoint</param>
    /// <param name="flatness">Precision of lettter outline curve in conversion to lines</param>
    /// <returns>inner polygon cw(outer ccw)</returns>
    static std::optional<Glyph> letter2glyph(const FontFile &font, unsigned int font_index, int letter, float flatness);

    /// <summary>
    /// Convert text into polygons
    /// </summary>
    /// <param name="font">Define fonts + cache, which could extend</param>
    /// <param name="text">Characters to convert</param>
    /// <param name="font_prop">User defined property of the font</param>
    /// <param name="was_canceled">Way to interupt processing</param>
    /// <returns>Inner polygon cw(outer ccw)</returns>
    static ExPolygons text2shapes(FontFileWithCache &font, const char *text, const FontProp &font_prop, std::function<bool()> was_canceled = nullptr);

    /// <summary>
    /// Fix intersections and self intersections in polygons glyph shape 
    /// </summary>
    /// <param name="shape">Input shape to heal</param>
    /// <returns>Healed shapes</returns>
    static ExPolygons heal_shape(const Polygons &shape);

    /// <summary>
    /// Use data from font property to modify transformation
    /// </summary>
    /// <param name="font_prop">Z-move as surface distance(FontProp::distance)
    /// Z-rotation as angle to Y axis(FontProp::angle)</param>
    /// <param name="transformation">In / Out transformation to modify by property</param>
    static void apply_transformation(const FontProp &font_prop, Transform3d &transformation);

    /// <summary>
    /// Read information from naming table of font file
    /// search for italic (or oblique), bold italic (or bold oblique)
    /// </summary>
    /// <param name="font">Selector of font</param>
    /// <param name="font_index">Index of font in collection</param>
    /// <returns>True when the font description contains italic/obligue otherwise False</returns>
    static bool is_italic(const FontFile &font, unsigned int font_index);

    /// <summary>
    /// Create unique character set from string with filtered from text with only character from font
    /// </summary>
    /// <param name="text">Source vector of glyphs</param>
    /// <param name="font">Font descriptor</param>
    /// <param name="font_index">Define font in collection</param>
    /// <param name="exist_unknown">True when text contain glyph unknown in font</param>
    /// <returns>Unique set of character from text contained in font</returns>
    static std::string create_range_text(const std::string &text, const FontFile &font, unsigned int font_index, bool* exist_unknown = nullptr);    

    /// <summary>
    /// calculate scale for glyph shape convert from shape points to mm
    /// </summary>
    /// <param name="fp"></param>
    /// <param name="ff"></param>
    /// <returns>Conversion to mm</returns>
    static double get_shape_scale(const FontProp &fp, const FontFile &ff);

    /// <summary>
    /// Project spatial point
    /// </summary>
    class IProject3f
    {
    public:
        virtual ~IProject3f() = default;
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
    /// Project 2d point into space
    /// Could be plane, sphere, cylindric, ...
    /// </summary>
    class IProjection : public IProject3f
    {
    public:
        virtual ~IProjection() = default;

        /// <summary>
        /// convert 2d point to 3d points
        /// </summary>
        /// <param name="p">2d coordinate</param>
        /// <returns>
        /// first - front spatial point
        /// second - back spatial point
        /// </returns>
        virtual std::pair<Vec3f, Vec3f> create_front_back(const Point &p) const = 0;

        virtual std::optional<Point> unproject(const Vec3f &p) const = 0;
    };

    /// <summary>
    /// Create triangle model for text
    /// </summary>
    /// <param name="shape2d">text or image</param>
    /// <param name="projection">Define transformation from 2d to 3d(orientation, position, scale, ...)</param>
    /// <returns>Projected shape into space</returns>
    static indexed_triangle_set polygons2model(const ExPolygons &shape2d, const IProjection& projection);
        
    /// <summary>
    /// Create transformation for emboss text object to lay on surface point
    /// </summary>
    /// <param name="position">Position of surface point</param>
    /// <param name="normal">Normal of surface point</param>
    /// <param name="up_limit">Is compared with normal.z to suggest up direction</param>
    /// <returns>Transformation onto surface point</returns>
    static Transform3d create_transformation_onto_surface(
        const Vec3f &position, const Vec3f &normal, float up_limit = 0.9f);

    class ProjectZ : public IProjection
    {
    public:
        ProjectZ(float depth) : m_depth(depth) {}
        // Inherited via IProject
        std::pair<Vec3f, Vec3f> create_front_back(const Point &p) const override;
        Vec3f project(const Vec3f &point) const override;
        std::optional<Point> unproject(const Vec3f &p) const override;
        float m_depth;
    };

    class ProjectScale : public IProjection
    {
        std::unique_ptr<IProjection> core;
        float m_scale;
    public:
        ProjectScale(std::unique_ptr<IProjection> core, float scale)
            : core(std::move(core)), m_scale(scale)
        {}

        // Inherited via IProject
        std::pair<Vec3f, Vec3f> create_front_back(const Point &p) const override
        {
            auto res = core->create_front_back(p);
            return std::make_pair(res.first * m_scale, res.second * m_scale);
        }
        Vec3f project(const Vec3f &point) const override{
            return core->project(point);
        }
        std::optional<Point> unproject(const Vec3f &p) const override {
            return core->unproject(p / m_scale);
        }
    };

    class OrthoProject3f : public Emboss::IProject3f
    {
        // size and direction of emboss for ortho projection
        Vec3f m_direction;
    public:
        OrthoProject3f(Vec3f direction) : m_direction(direction) {}
        Vec3f project(const Vec3f &point) const override{ return point + m_direction;}
    };

    class OrthoProject: public Emboss::IProjection {
        Transform3d m_matrix;
        // size and direction of emboss for ortho projection
        Vec3f       m_direction;
        Transform3d m_matrix_inv;
    public:
        OrthoProject(Transform3d matrix, Vec3f direction)
            : m_matrix(matrix), m_direction(direction), m_matrix_inv(matrix.inverse())
        {}
        // Inherited via IProject
        std::pair<Vec3f, Vec3f> create_front_back(const Point &p) const override;
        Vec3f project(const Vec3f &point) const override;
        std::optional<Point> unproject(const Vec3f &p) const override;     
    };
};

} // namespace Slic3r
#endif // slic3r_Emboss_hpp_
