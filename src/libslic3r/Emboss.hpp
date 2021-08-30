#ifndef slic3r_Emboss_hpp_
#define slic3r_Emboss_hpp_

#include <vector>
#include <set>
#include <optional>
#include <memory>
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

        // enum class Align: center/left/right
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

    /// <summary>
    /// Convert text into polygons
    /// </summary>
    /// <param name="font">Define fonts</param>
    /// <param name="letter">Character to convert</param>
    /// <returns>inner polygon ccw(outer cw)</returns>
    static Polygons text2polygons(const Font &font, const std::string &text);

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
    static indexed_triangle_set polygons2model(const Polygons &shape2d, const IProject& projection);

    /// <summary>
    /// Connect points by triangulation to create filled surface by triangle indices
    /// </summary>
    /// <param name="points">Points to connect</param>
    /// <param name="edges">Constraint for edges, pair is from point(first) to point(second)</param>
    /// <returns>Triangles</returns>
    static std::vector<Vec3i> triangulate(const Points &points, const std::set<std::pair<uint32_t, uint32_t>> &edges);
    static std::vector<Vec3i> triangulate(const Polygon &polygon);
    static std::vector<Vec3i> triangulate(const Polygons &polygons);

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

        float                           m_scale;
    };
};

} // namespace Slic3r
#endif // slic3r_Emboss_hpp_
