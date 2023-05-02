#ifndef slic3r_TextLines_hpp_
#define slic3r_TextLines_hpp_

#include <vector>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/Point.hpp>
#include <libslic3r/Emboss.hpp>
#include "slic3r/GUI/GLModel.hpp"


namespace Slic3r::GUI {

class Selection;

/// <summary>
/// Define polygon for draw letters
/// </summary>
struct TextLine
{
    // slice of object
    Polygon polygon;

    // index to point in polygon which starts line, which is closest to zero
    size_t start_index;

    // Point on line closest to zero
    Point start_point;
};
using TextLines = std::vector<TextLine>;

class TextLinesModel
{
public:
    // line_height in mm
    void init(const Selection &selection, double line_height);
    void render(const Transform3d &text_world);

    bool is_init() const { return m_model.is_initialized(); }
    void reset() { m_model.reset(); }
    const TextLines &get_lines() const { return m_lines; }
    static double calc_line_height(const Slic3r::Emboss::FontFile& ff, const FontProp& fp);
private:
    TextLines m_lines;

    // Keep model for visualization text lines
    GLModel m_model;
};

} // namespace Slic3r::GUI
#endif // slic3r_TextLines_hpp_