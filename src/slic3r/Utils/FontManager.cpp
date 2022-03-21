#include "FontManager.hpp"
#include <wx/font.h>
#include <GL/glew.h> // Imgui texture
#include <imgui/imgui_internal.h> // ImTextCharFromUtf8
#include "WxFontUtils.hpp"
#include "libslic3r/Utils.hpp" // ScopeGuard
#include "slic3r/GUI/3DScene.hpp" // ::glsafe

using namespace Slic3r;
using namespace Slic3r::GUI;

FontManager::FontManager(const ImWchar *language_glyph_range)
    : m_imgui_init_glyph_range(language_glyph_range)
    , m_font_selected(std::numeric_limits<size_t>::max())
    , m_exist_style_images(false)
{}

FontManager::~FontManager() { 
    free_imgui_fonts();
    free_style_images();
}

void FontManager::swap(size_t i1, size_t i2) {
    if (i1 >= m_font_list.size() || 
        i2 >= m_font_list.size()) return;
    std::swap(m_font_list[i1], m_font_list[i2]);

    // fix selected index
    if (!is_activ_font()) return;
    if (m_font_selected == i1) {
        m_font_selected = i2;
    } else if (m_font_selected == i2) {
        m_font_selected = i1;
    }
}

void FontManager::duplicate() { duplicate(m_font_selected); }
void FontManager::duplicate(size_t index) {
    if (index >= m_font_list.size()) return;
    Item item = m_font_list[index]; // copy
    make_unique_name(item.font_item.name);
    item.truncated_name.clear();    

    // take original font imgui pointer
    //ImFont *imgui_font = get_imgui_font(index);
    //if (imgui_font != nullptr)
    //    m_font_list[index].imgui_font_index.reset();

    m_font_list.insert(m_font_list.begin() + index, item);
    // fix selected index
    if (!is_activ_font()) return;
    if (index < m_font_selected) ++m_font_selected;
}

void FontManager::erase(size_t index) {
    if (index >= m_font_list.size()) return;
    //ImFont *imgui_font = get_imgui_font(index);
    //if (imgui_font != nullptr)
    //    IM_DELETE(imgui_font);

    // fix selected index
    if (is_activ_font() && index < m_font_selected)
        --m_font_selected;

    m_font_list.erase(m_font_list.begin() + index);
}

bool FontManager::wx_font_changed(std::unique_ptr<Emboss::FontFile> font_file)
{
    if (!is_activ_font()) return false;
    auto &wx_font = get_wx_font();
    if (!wx_font.has_value()) return false;

    if (font_file == nullptr) {        
        font_file = WxFontUtils::create_font_file(*wx_font);
        if (font_file == nullptr) return false;
    }
    m_font_list[m_font_selected].font_file_with_cache = 
        Emboss::FontFileWithCache(std::move(font_file));
        
    auto &fi = get_font_item();
    fi.type  = WxFontUtils::get_actual_type();
    fi.path = WxFontUtils::store_wxFont(*wx_font);
    clear_imgui_font();
    free_style_images();
    return true;
}

bool FontManager::load_font(size_t font_index)
{
    if (font_index >= m_font_list.size()) return false;
    std::swap(font_index, m_font_selected);
    bool is_loaded = load_active_font();
    if (!is_loaded) std::swap(font_index, m_font_selected);
    return is_loaded;
}

bool FontManager::load_font(size_t font_index, const wxFont &font)
{
    if (font_index >= m_font_list.size()) return false;
    std::swap(font_index, m_font_selected);
    bool is_loaded = set_wx_font(font);
    if (!is_loaded) std::swap(font_index, m_font_selected);
    return is_loaded;
}

static std::string get_file_name(const std::string &file_path)
{
    size_t pos_last_delimiter = file_path.find_last_of('\\');
    size_t pos_point          = file_path.find_last_of('.');
    size_t offset             = pos_last_delimiter + 1;
    size_t count              = pos_point - pos_last_delimiter - 1;
    return file_path.substr(offset, count);
}

bool FontManager::load_active_font()
{ 
    return set_up_font_file(m_font_selected); 
}

bool FontManager::is_activ_font() {
    return m_font_selected < m_font_list.size();
}

bool FontManager::load_first_valid_font() {
    while (!m_font_list.empty()) {
        if (load_font(0)) return true;
        // can't load so erase it from list
        m_font_list.erase(m_font_list.begin());
    }
    return false;
}

void FontManager::add_font(FontItem font_item)
{ 
    make_unique_name(font_item.name);
    Item item;
    item.font_item = font_item;
    m_font_list.push_back(item); 
}

void FontManager::add_fonts(FontList font_list)
{
    for (const FontItem &fi : font_list) 
        add_font(fi);
}

std::shared_ptr<const Emboss::FontFile> &FontManager::get_font_file()
{
    // TODO: fix not selected font
    //if (!is_activ_font()) return nullptr;
    return m_font_list[m_font_selected].font_file_with_cache.font_file;
}

const FontItem &FontManager::get_font_item() const
{
    // TODO: fix not selected font
    //if (!is_activ_font()) return nullptr;
    return m_font_list[m_font_selected].font_item;
}

FontItem &FontManager::get_font_item()
{
    // TODO: fix not selected font
    //if (!is_activ_font()) return nullptr;
    return m_font_list[m_font_selected].font_item;
}

const FontProp &FontManager::get_font_prop() const
{
    // TODO: fix not selected font
    //if (!is_activ_font()) return nullptr;
    return m_font_list[m_font_selected].font_item.prop;
}

FontProp &FontManager::get_font_prop()
{
    // TODO: fix not selected font
    return m_font_list[m_font_selected].font_item.prop;
}

std::optional<wxFont> &FontManager::get_wx_font()
{
    //if (!is_activ_font()) return {};
    return m_font_list[m_font_selected].wx_font;

    //std::optional<wxFont> &wx_font = m_font_list[m_font_selected].wx_font;
    //if (wx_font.has_value()) return wx_font;

    //const FontItem &fi = get_font_item();
    //if (fi.type != WxFontUtils::get_actual_type()) return {};

    //wx_font = WxFontUtils::load_wxFont(fi.path);
    //return wx_font;
}

bool FontManager::set_wx_font(const wxFont &wx_font) { 
    return set_wx_font(m_font_selected, wx_font);
}

std::string &FontManager::get_truncated_name()
{
    return m_font_list[m_font_selected].truncated_name;
}

const std::optional<wxFont> &FontManager::get_wx_font() const
{
    return m_font_list[m_font_selected].wx_font;
}

void FontManager::clear_glyphs_cache()
{
    if (!is_activ_font()) return;
    Emboss::FontFileWithCache &ff = m_font_list[m_font_selected].font_file_with_cache;
    if (!ff.has_value()) return;
    ff.cache = std::make_shared<Emboss::Glyphs>();
}

void FontManager::clear_imgui_font() { 
    // TODO: improove to clear only actual font
    if (!is_activ_font()) return;
    free_imgui_fonts();
    return;    
    ImFont *imgui_font = get_imgui_font(m_font_selected);
    m_font_list[m_font_selected].imgui_font_index.reset();
    if (imgui_font != nullptr) IM_DELETE(imgui_font);
}

ImFont *FontManager::get_imgui_font(const std::string &text)
{
    return get_imgui_font(m_font_selected, text);
}

ImFont *FontManager::get_imgui_font(size_t item_index, const std::string &text)
{    
    Item &item = m_font_list[item_index];
    // check is already loaded
    if (!item.imgui_font_index.has_value())
        return load_imgui_font(item_index, text);

    size_t index = *item.imgui_font_index;
    auto & fonts = m_imgui_font_atlas.Fonts;

    // check correct index
    int f_size = fonts.size();
    assert(f_size > 0 && index < (size_t)f_size);
    if (f_size <= 0 || index >= (size_t) f_size) return nullptr;
    ImFont *font = fonts[index];
    if (font == nullptr) return nullptr;
    if (!font->IsLoaded()) return nullptr;
    if (font->Scale <= 0.f) return nullptr;
    // automatic extend range
    if (!text.empty() && !is_text_in_ranges(font, text)) 
        return extend_imgui_font_range(item_index, text);

    return font;
}

void FontManager::free_except_active_font() { 
    free_imgui_fonts();

    // free font_files
    const Item &act_item = m_font_list[m_font_selected];
    for (auto &item : m_font_list) {
        if (&item == &act_item) continue; // keep alive actual font file
        item.font_file_with_cache = {};
    }
}

const std::vector<FontManager::Item> &FontManager::get_fonts() const
{
    return m_font_list;
}

const FontManager::Item &FontManager::get_font() const
{
    return m_font_list[m_font_selected];
}

void FontManager::make_unique_name(std::string &name)
{
    auto is_unique = [&](const std::string &name) -> bool {
        for (const Item &it : m_font_list)
            if (it.font_item.name == name) return false;
        return true;
    };

    if (name.empty()) name = "font";
    if (is_unique(name)) return;

    auto pos = name.find(" (");
    if (pos != std::string::npos && name.find(")", pos) != std::string::npos) {
        // short name by ord number
        name = name.substr(0, pos);
    }

    int         order = 1; // start with value 2 to represents same font name
    std::string new_name;
    do {
        new_name = name + " (" + std::to_string(++order) + ")";
    } while (!is_unique(new_name));
    name = new_name;
}

bool FontManager::set_up_font_file(size_t item_index)
{ 
    Item &item = m_font_list[item_index];
    FontItem &fi = item.font_item;
    if (fi.type == FontItem::Type::file_path) {
        // fill font name after load from .3mf
        if (fi.name.empty()) fi.name = get_file_name(fi.path);
        std::unique_ptr<Emboss::FontFile> font_ptr =
            Emboss::create_font_file(fi.path.c_str());
        if (font_ptr == nullptr) return false;
        item.font_file_with_cache = 
            Emboss::FontFileWithCache(std::move(font_ptr));
        return true;
    }
    if (fi.type != WxFontUtils::get_actual_type()) return false;
    if (!item.wx_font.has_value())
        item.wx_font = WxFontUtils::load_wxFont(fi.path);
    if (!item.wx_font.has_value()) return false;
    return set_wx_font(item_index, *item.wx_font);
}

bool FontManager::is_text_in_ranges(const ImFont *font, const std::string &text)
{
    if (font == nullptr) return false;
    if (!font->IsLoaded()) return false;
    const ImFontConfig *fc = font->ConfigData;
    if (fc == nullptr) return false;
    return is_text_in_ranges(fc->GlyphRanges, text);
}

bool FontManager::is_char_in_ranges(const ImWchar *ranges, unsigned int letter)
{
    for (const ImWchar *range = ranges; range[0] && range[1]; range += 2) {
        ImWchar from = range[0];
        ImWchar to   = range[1];
        if (from <= letter && letter <= to) return true;
        if (letter < to) return false; // ranges should be sorted
    }
    return false;
};

bool FontManager::is_text_in_ranges(const ImWchar *ranges, const std::string &text)
{
    const char *text_char_ptr = text.c_str();
    while (*text_char_ptr) {
        unsigned int c = 0;
        // UTF-8 to 32-bit character need imgui_internal
        int c_len = ImTextCharFromUtf8(&c, text_char_ptr, NULL);
        text_char_ptr += c_len;
        if (c_len == 0) break;
        if (!is_char_in_ranges(ranges, c)) return false;
    }
    return true;
}

ImFont* FontManager::extend_imgui_font_range(size_t index, const std::string& text)
{
    auto &font_index_opt = m_font_list[m_font_selected].imgui_font_index;
    if (!font_index_opt.has_value()) 
        return load_imgui_font(index, text);

    // TODO: start using merge mode
    // ImFontConfig::MergeMode = true;

    free_imgui_fonts();
    return load_imgui_font(index, text);
}

#include "libslic3r/SLA/AGGRaster.hpp"
void FontManager::create_texture(size_t index, const std::string &text, GLuint& tex_id, ImVec2& tex_size)
{
    if (index >= m_font_list.size()) return;
    Item &item = m_font_list[index];
    if (!item.font_file_with_cache.has_value() && !set_up_font_file(index))
        return;

    const FontProp &font_prop = item.font_item.prop;
    ExPolygons      shapes    = Emboss::text2shapes(item.font_file_with_cache,
                                                    text.c_str(), font_prop);

    BoundingBox bb;
    for (ExPolygon &shape : shapes) bb.merge(BoundingBox(shape.contour.points));
    for (ExPolygon &shape : shapes) shape.translate(-bb.min);        

    double scale = font_prop.size_in_mm;
    BoundingBoxf bb2 = unscaled(bb);
    bb2.scale(scale);
    tex_size.x       = bb2.max.x() - bb2.min.x();
    tex_size.y       = bb2.max.y() - bb2.min.y();
    sla::Resolution resolution(tex_size.x,tex_size.y);
    sla::PixelDim   dim(1/scale, 1/scale);
    const double no_gamma = 1.;
    std::unique_ptr<sla::RasterBase> r =
        sla::create_raster_grayscale_aa(resolution, dim, no_gamma);
    for (const ExPolygon &shape : shapes) r->draw(shape);      
    // reserve texture on GPU
    glGenTextures(1, &tex_id);
    glBindTexture(GL_TEXTURE_2D, tex_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    sla::RasterEncoder encoder = [](const void *ptr, size_t w, size_t h, size_t num_components) -> sla::EncodedRaster {
        GLsizei width = w, height = h;
        GLint   border = 0;
        glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, width, height, border, GL_ALPHA, GL_UNSIGNED_BYTE, ptr);
        return sla::EncodedRaster();
    };
    r->encode(encoder);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// for get DPI
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/MainFrame.hpp"

void FontManager::init_style_images(int max_width) {
    // check already initialized
    if (m_exist_style_images) return;

    // create shapes and calc size (bounding boxes)
    std::vector<ExPolygons> name_shapes(m_font_list.size());
    std::vector<double> scales(m_font_list.size());
    for (Item &item : m_font_list) {
        FontItem &        font_item = item.font_item;
        const FontProp &  font_prop = font_item.prop;
        size_t index = &item - &m_font_list.front();
        if (!item.font_file_with_cache.has_value() && !set_up_font_file(index))
            continue;

        ExPolygons &shapes = name_shapes[index];
        shapes = Emboss::text2shapes(item.font_file_with_cache, font_item.name.c_str(), font_prop);

        // create image description
        item.image = StyleImage();
        StyleImage &image = *item.image;

        BoundingBox &bounding_box = image.bounding_box;
        for (ExPolygon &shape : shapes)
            bounding_box.merge(BoundingBox(shape.contour.points));
        for (ExPolygon &shape : shapes) shape.translate(-bounding_box.min);
        
        // calculate conversion from FontPoint to screen pixels by size of font
        auto   mf  = wxGetApp().mainframe;
        // dot per inch for monitor
        int    dpi = get_dpi_for_window(mf);
        double ppm = dpi / 25.4; // pixel per milimeter
        double unit_per_em = item.font_file_with_cache.font_file->unit_per_em;
        double scale = font_prop.size_in_mm / unit_per_em * Emboss::SHAPE_SCALE * ppm;
        scales[index] = scale;

        //double scale = font_prop.size_in_mm * SCALING_FACTOR;
        BoundingBoxf bb2(bounding_box.min.cast<double>(),
                         bounding_box.max.cast<double>());
        bb2.scale(scale);
        image.tex_size.x = std::ceil(bb2.max.x() - bb2.min.x());
        image.tex_size.y = std::ceil(bb2.max.y() - bb2.min.y());
        // crop image width
        if (image.tex_size.x > max_width) 
            image.tex_size.x = max_width;
    }

    // arrange bounding boxes
    int offset_y = 0;
    int width    = 0;
    for (Item &item : m_font_list) {
        if (!item.image.has_value()) continue;
        StyleImage &image = *item.image;
        image.offset.y() = offset_y;
        offset_y += image.tex_size.y+1;
        if (width < image.tex_size.x) 
            width = image.tex_size.x;
    }
    int height = offset_y;
    for (Item &item : m_font_list) {
        if (!item.image.has_value()) continue;
        StyleImage &image = *item.image;
        const Point &o = image.offset;
        const ImVec2 &s = image.tex_size;
        image.uv0 = ImVec2(o.x() / (double) width, 
                           o.y() / (double) height);
        image.uv1 = ImVec2((o.x() + s.x) / (double) width,
                           (o.y() + s.y) / (double) height);
    }

    // reserve texture on GPU
    GLuint tex_id;
    GLenum target = GL_TEXTURE_2D, format = GL_ALPHA, type = GL_UNSIGNED_BYTE;
    GLint level = 0, border = 0;
    glsafe(::glGenTextures(1, &tex_id));
    glsafe(::glBindTexture(target, tex_id));
    glsafe(::glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
    glsafe(::glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
    // texture size
    GLint w = width, h = height;
    glsafe(::glTexImage2D(target, level, GL_ALPHA, w, h, border, format, type, nullptr));

    // set up texture id
    void *texture_id = (void *)(intptr_t) tex_id;
    for (Item &item : m_font_list)
        if (item.image.has_value())
            item.image->texture_id = texture_id;

    // upload sub textures
    for (Item &item : m_font_list) {
        if (!item.image.has_value()) continue;
        StyleImage &image = *item.image;
        sla::Resolution resolution(image.tex_size.x, image.tex_size.y);

        size_t index = &item - &m_font_list.front();
        double pixel_dim = SCALING_FACTOR / scales[index];
        sla::PixelDim dim(pixel_dim, pixel_dim);
        double gamma = 1.;
        std::unique_ptr<sla::RasterBase> r = sla::create_raster_grayscale_aa(resolution, dim, gamma);
        for (const ExPolygon &shape : name_shapes[index]) r->draw(shape);
        const Point& offset = image.offset;
        sla::RasterEncoder encoder = 
            [offset, target, level, format, type]
            (const void *ptr, size_t w, size_t h, size_t num_components) {
            GLint sub_w = w, sub_h = h, xoffset = offset.x(), yoffset = offset.y();
            glsafe(::glTexSubImage2D(target, level, xoffset, yoffset, sub_w, sub_h, format, type, ptr));
            return sla::EncodedRaster();
        };
        // upload texture data to GPU
        r->encode(encoder);
    }

    // bind default texture
    GLuint no_texture_id = 0;
    glsafe(::glBindTexture(target, no_texture_id));

    m_exist_style_images = true;
}

void FontManager::free_style_images() {
    if (!is_activ_font()) return;
    if (!m_exist_style_images) return;
    GLuint tex_id = (GLuint) (intptr_t) m_font_list.front().image->texture_id;
    for (Item &it : m_font_list) it.image.reset();

    glsafe(::glDeleteTextures(1, &tex_id));
    m_exist_style_images = false;
}


void FontManager::free_imgui_fonts()
{ 
    for (auto &item : m_font_list) 
        item.imgui_font_index.reset();
    m_imgui_font_atlas.Clear();
}

float FontManager::min_imgui_font_size = 18.f;
float FontManager::max_imgui_font_size = 60.f;
float FontManager::get_imgui_font_size(const FontProp         &prop,
                                       const Emboss::FontFile &file)
{
    // coeficient for convert line height to font size
    float c1 = (file.ascent - file.descent + file.linegap) / (float) file.unit_per_em;

    // The point size is defined as 1/72 of the Anglo-Saxon inch (25.4 mm):
    // It is approximately 0.0139 inch or 352.8 um.
    return c1 * std::abs(prop.size_in_mm) / 0.3528f;
}

ImFont * FontManager::load_imgui_font(size_t index, const std::string &text)
{
    free_imgui_fonts(); // TODO: remove it after correct initialization

    if (index >= m_font_list.size()) return nullptr;
    Item &item = m_font_list[index];
    if (!item.font_file_with_cache.has_value()) return nullptr;
    const Emboss::FontFile &font_file = *item.font_file_with_cache.font_file;

    // TODO: Create glyph range
    ImFontGlyphRangesBuilder builder;
    builder.AddRanges(m_imgui_init_glyph_range);
    if (!text.empty())
        builder.AddText(text.c_str());
    item.font_ranges.clear();

    builder.BuildRanges(&item.font_ranges);
        
    m_imgui_font_atlas.Flags |= ImFontAtlasFlags_NoMouseCursors |
                                ImFontAtlasFlags_NoPowerOfTwoHeight;

    const FontProp &font_prop = item.font_item.prop;
    float font_size = get_imgui_font_size(font_prop, font_file);
    if (font_size < min_imgui_font_size)
        font_size = min_imgui_font_size;
    if (font_size > max_imgui_font_size)
        font_size = max_imgui_font_size;

    ImFontConfig font_config;
    // TODO: start using merge mode
    //font_config.MergeMode = true;
    if (font_prop.char_gap.has_value()) {
        float coef = font_size / (double) font_file.unit_per_em;
        font_config.GlyphExtraSpacing.x = coef * (*font_prop.char_gap);
    }
    if (font_prop.line_gap.has_value()) {
        float coef = font_size / (double) font_file.unit_per_em;
        font_config.GlyphExtraSpacing.y = coef * (*font_prop.line_gap);
    }

    font_config.FontDataOwnedByAtlas = false;

    const std::vector<unsigned char> &buffer = *font_file.data;
    m_imgui_font_atlas.AddFontFromMemoryTTF(
        (void *) buffer.data(), buffer.size(), font_size, &font_config, item.font_ranges.Data);

    unsigned char *pixels;
    int            width, height;
    m_imgui_font_atlas.GetTexDataAsAlpha8(&pixels, &width, &height);

    // Upload texture to graphics system
    GLint last_texture;
    glsafe(::glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture));
    ScopeGuard sg([last_texture]() {
        glsafe(::glBindTexture(GL_TEXTURE_2D, last_texture));
    });

    GLuint font_texture;
    glsafe(::glGenTextures(1, &font_texture));
    glsafe(::glBindTexture(GL_TEXTURE_2D, font_texture));
    glsafe(::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
    glsafe(::glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
    glsafe(::glPixelStorei(GL_UNPACK_ROW_LENGTH, 0));
    glsafe(::glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, width, height, 0,
                          GL_ALPHA, GL_UNSIGNED_BYTE, pixels));

    // Store our identifier
    m_imgui_font_atlas.TexID = (ImTextureID) (intptr_t) font_texture;
    assert(!m_imgui_font_atlas.Fonts.empty());
    if (m_imgui_font_atlas.Fonts.empty()) return nullptr;
    item.imgui_font_index = m_imgui_font_atlas.Fonts.size() - 1;
    return m_imgui_font_atlas.Fonts.back();
}

bool FontManager::set_wx_font(size_t item_index, const wxFont &wx_font) {
    std::unique_ptr<Emboss::FontFile> font_file = 
        WxFontUtils::create_font_file(wx_font);
    if (font_file == nullptr) return false;

    Item &item = m_font_list[item_index];
    item.font_file_with_cache = 
        Emboss::FontFileWithCache(std::move(font_file));
    item.wx_font = wx_font;

    FontItem &fi = item.font_item;
    fi.type      = WxFontUtils::get_actual_type();

    // fill font name after load from .3mf
    if (fi.name.empty())
        fi.name = WxFontUtils::get_human_readable_name(*item.wx_font);

    clear_imgui_font();
    return true;
}
