#include "Emboss.hpp"
#include <stdio.h>

#define STB_TRUETYPE_IMPLEMENTATION // force following include to generate implementation
#include "imgui/imstb_truetype.h" // stbtt_fontinfo

using namespace Slic3r;

std::optional<Emboss::Font> Emboss::load_font(const char *file_path)
{
    FILE *file = fopen(file_path, "rb");
    if (file == nullptr) {
        std::cerr << "Couldn't open " << file_path << " for reading.";
        return {};
    }

    // find size of file
    if (fseek(file, 0L, SEEK_END) != 0) {
        std::cerr << "Couldn't fseek file " << file_path << " for size measure.";
        return {};
    }
    size_t size = ftell(file);
    if (size == 0) {
        std::cerr << "Size of font file is zero. Can't read.";
        return {};    
    }
    rewind(file);

    Font res;
    res.buffer = std::vector<unsigned char>(size);
    size_t count_loaded_bytes = fread((void *) &res.buffer.front(), 1, size, file);

    unsigned int index = 0;
    int font_offset = 0;
    while (font_offset >= 0) {
        font_offset = stbtt_GetFontOffsetForIndex(res.buffer.data(), index++);
    }
    // at least one font must be inside collection
    if (index < 1) {
        std::cerr << "There is no font collection inside file.";
        return {};        
    }
    // select default font on index 0
    res.index = 0;
    res.count = index;
    // first fonst has offset zero
    font_offset = 0;

    stbtt_fontinfo font_info;
    if (stbtt_InitFont(&font_info, res.buffer.data(), font_offset) == 0) {
        std::cerr << "Can't initialize font.";
        return {};    
    }

    // load information about line gap
    stbtt_GetFontVMetrics(&(font_info), &res.ascent, &res.descent, &res.linegap);
    return res;
}

Polygons Emboss::letter2polygons(const Font &font, char letter)
{
    int font_offset = stbtt_GetFontOffsetForIndex(font.buffer.data(), font.index);
    stbtt_fontinfo font_info;
    if (stbtt_InitFont(&font_info, font.buffer.data(), font_offset) == 0)
        return Polygons();

    int glyph_index = stbtt_FindGlyphIndex(&(font_info), letter);
    int advanceWidth, leftSideBearing;
    stbtt_GetGlyphHMetrics(&(font_info), glyph_index, &advanceWidth,
                           &leftSideBearing);

    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&(font_info), glyph_index, &vertices);
    if (num_verts < 0) return {}; // no shape

    int *         contour_lengths = NULL;
    int           num_countour    = 0;
    stbtt__point *points = stbtt_FlattenCurves(vertices, num_verts, font.flatness,
                                               &contour_lengths,
                                               &num_countour,
                                               (font_info).userdata);

    Polygons result;
    result.reserve(num_countour);
    size_t pi = 0; // point index
    for (size_t ci = 0; ci < num_countour; ++ci) {
        int    length = contour_lengths[ci];
        Points pts;
        pts.reserve(length);
        for (size_t i = 0; i < length; i++) {
            const stbtt__point &point = points[pi];
            pi++;
            pts.emplace_back(point.x, point.y);
        }
        result.emplace_back(pts);
    }

    BoundingBox bb;
    for (auto &r : result) bb.merge(r.points);

    // inner ccw
    // outer cw
    return result;
}

indexed_triangle_set Emboss::create_model(const std::string &text,
                                          const Font &       font,
                                          float              z_size)
{
    return indexed_triangle_set();
}
