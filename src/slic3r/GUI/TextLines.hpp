#ifndef slic3r_TextLines_hpp_
#define slic3r_TextLines_hpp_

#include <vector>
#include <libslic3r/Polygon.hpp>
#include <libslic3r/Point.hpp>
#include <libslic3r/Emboss.hpp>
#include "slic3r/GUI/GLModel.hpp"

namespace Slic3r {
class ModelVolume;
typedef std::vector<ModelVolume *> ModelVolumePtrs;
}

namespace Slic3r::GUI {
class TextLinesModel
{
public:
    // line offset in y direction (up/down)
    float offset = 0;

    /// <summary>
    /// Initialize model and lines
    /// </summary>
    /// <param name="text_tr">Transformation of text volume inside object (aka inside of instance)</param>
    /// <param name="volumes_to_slice">Vector of volumes to be sliced</param>
    /// <param name="align">Vertical (Y) align of the text</param>
    /// <param name="line_height">Distance between lines [in mm]</param>
    /// <param name="count_lines">Count lines(slices over volumes)</param>
    void init(const Transform3d &text_tr, const ModelVolumePtrs& volumes_to_slice, FontProp::Align align, double line_height, unsigned count_lines);

    void render(const Transform3d &text_world);

    bool is_init() const { return m_model.is_initialized(); }
    void reset() { m_model.reset(); m_lines.clear(); }
    const Slic3r::Emboss::TextLines &get_lines() const { return m_lines; }

    static double calc_line_height(const Slic3r::Emboss::FontFile& ff, const FontProp& fp); // return lineheight in mm
private:
    Slic3r::Emboss::TextLines m_lines;

    // Keep model for visualization text lines
    GLModel m_model;
};

} // namespace Slic3r::GUI
#endif // slic3r_TextLines_hpp_