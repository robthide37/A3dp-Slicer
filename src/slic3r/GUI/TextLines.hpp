#ifndef slic3r_TextLines_hpp_
#define slic3r_TextLines_hpp_

#include <vector>
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
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
    void init(const Selection &selection);
    void render(const Transform3d &text_world);

    bool is_init() const { return model.is_initialized(); }
    void reset() { model.reset(); }
    const TextLines &get_lines() const { return lines; }
private:
    TextLines lines;

    // Keep model for visualization text lines
    GLModel model;
};

} // namespace Slic3r::GUI
#endif // slic3r_TextLines_hpp_