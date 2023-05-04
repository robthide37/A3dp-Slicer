#ifndef slic3r_TextLines_hpp_
#define slic3r_TextLines_hpp_

#include <vector>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/Point.hpp>
#include <libslic3r/Emboss.hpp>
#include "slic3r/GUI/GLModel.hpp"

namespace Slic3r::GUI {
class Selection;
class TextLinesModel
{
public:
    /// <summary>
    /// Initialize model and lines
    /// </summary>
    /// <param name="selection">Must be selected text volume</param>
    /// <param name="line_height">Height of text line with spacing [in mm]</param>
    /// <param name="line_offset">Offset of base line from center [in mm]</param>
    void init(const Selection &selection, double line_height, double line_offset = 0.);
    void render(const Transform3d &text_world);

    bool is_init() const { return m_model.is_initialized(); }
    void reset() { m_model.reset(); }
    const Slic3r::Emboss::TextLines &get_lines() const { return m_lines; }
    static double calc_line_height(const Slic3r::Emboss::FontFile& ff, const FontProp& fp);
private:
    Slic3r::Emboss::TextLines m_lines;

    // Keep model for visualization text lines
    GLModel m_model;
};

} // namespace Slic3r::GUI
#endif // slic3r_TextLines_hpp_