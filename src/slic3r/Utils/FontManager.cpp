#include "FontManager.hpp"
#include <wx/font.h>
#include <GL/glew.h> // Imgui texture
#include <imgui/imgui_internal.h> // ImTextCharFromUtf8
#include "WxFontUtils.hpp"
#include "libslic3r/Utils.hpp" // ScopeGuard

#include "slic3r/GUI/3DScene.hpp" // ::glsafe
#include "slic3r/GUI/Jobs/CreateFontStyleImagesJob.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp" // check of font ranges

using namespace Slic3r;
using namespace Slic3r::GUI;

FontManager::FontManager(const ImWchar *language_glyph_range)
    : m_imgui_init_glyph_range(language_glyph_range)
    , m_font_selected(std::numeric_limits<size_t>::max())
    , m_exist_style_images(false)
    , m_temp_style_images(nullptr)
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
    size_t pos_last_delimiter = file_path.find_last_of("/\\");
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

ImFont *FontManager::get_imgui_font()
{
    return get_imgui_font(m_font_selected);
}

ImFont *FontManager::create_imgui_font(const std::string &text)
{
    return create_imgui_font(m_font_selected, text);
}

ImFont *FontManager::get_imgui_font(size_t item_index)
{
    if (!is_activ_font()) return nullptr;
    Item &item = m_font_list[item_index];
    // check is already loaded
    if (!item.imgui_font_index.has_value()) return nullptr;

    size_t              index = *item.imgui_font_index;
    ImVector<ImFont *> &fonts = m_imgui_font_atlas.Fonts;

    // check correct index
    int f_size = fonts.size();
    assert(f_size > 0 && index < (size_t) f_size);
    if (f_size <= 0 || index >= (size_t) f_size) return nullptr;
    ImFont *font = fonts[index];
    if (font == nullptr) return nullptr;
    if (!font->IsLoaded()) return nullptr;
    if (font->Scale <= 0.f) return nullptr;
    // Symbol fonts doesn't have atlas because their glyph range is out of language range
    if (font->ContainerAtlas == nullptr) return nullptr;
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

const std::vector<FontManager::Item> &FontManager::get_styles() const
{
    return m_font_list;
}

const FontManager::Item &FontManager::get_activ_style() const
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

ImFont* FontManager::extend_imgui_font_range(size_t index, const std::string& text)
{
    auto &font_index_opt = m_font_list[m_font_selected].imgui_font_index;
    if (!font_index_opt.has_value()) 
        return create_imgui_font(index, text);

    // TODO: start using merge mode
    // ImFontConfig::MergeMode = true;

    free_imgui_fonts();
    return create_imgui_font(index, text);
}


void FontManager::init_trunc_names(float max_width) { 
    for (auto &s : m_font_list) {
        if (s.truncated_name.empty())
            s.truncated_name = ImGuiWrapper::trunc(s.font_item.name, max_width);
    }
}

#include "slic3r/GUI/Jobs/CreateFontStyleImagesJob.hpp"

// for access to worker
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp" 

void FontManager::init_style_images(const Vec2i &max_size,
                                    const std::string &text)
{
    // check already initialized
    if (m_exist_style_images) return;

    // check is initializing
    if (m_temp_style_images != nullptr) {
        // is initialization finished
        if (!m_temp_style_images->styles.empty()) { 
            assert(m_temp_style_images->images.size() ==
                   m_temp_style_images->styles.size());
            // copy images into styles
            for (FontManager::StyleImage &image : m_temp_style_images->images){
                size_t index = &image - &m_temp_style_images->images.front();
                StyleImagesData::Item &style = m_temp_style_images->styles[index];

                // find style in font list and copy to it
                for (auto &it : m_font_list) {
                    if (it.font_item.name != style.text ||
                        !(it.font_item.prop == style.prop))
                        continue;
                    it.image = image;
                    break;
                }
            }
            m_temp_style_images = nullptr;
            m_exist_style_images = true;
            return;
        }
        // in process of initialization inside of job
        return;
    }

    // create job for init images
    m_temp_style_images = std::make_shared<StyleImagesData::StyleImages>();
    StyleImagesData::Items styles;
    styles.reserve(m_font_list.size());
    for (Item &item : m_font_list) {
        size_t index = &item - &m_font_list.front();
        if (!item.font_file_with_cache.has_value() && !set_up_font_file(index))
            continue;
        styles.push_back({
            item.font_file_with_cache, 
            item.font_item.name,
            item.font_item.prop
        });
    }

    auto &worker = wxGetApp().plater()->get_ui_job_worker();
    StyleImagesData data{std::move(styles), max_size, text, m_temp_style_images};
    queue_job(worker, std::make_unique<CreateFontStyleImagesJob>(std::move(data)));
}

void FontManager::free_style_images() {
    if (!m_exist_style_images) return;
    if (!is_activ_font()) return;

    GLuint tex_id = 0;
    
    for (Item &it : m_font_list) {
        if (tex_id == 0 && it.image.has_value())
            tex_id = (GLuint)(intptr_t) it.image->texture_id;
        it.image.reset();
    }
    if (tex_id != 0)
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
    const auto  &cn = prop.collection_number;
    unsigned int font_index = (cn.has_value()) ? *cn : 0;
    const auto  &font_info  = file.infos[font_index];
    // coeficient for convert line height to font size
    float c1 = (font_info.ascent - font_info.descent + font_info.linegap) /
               (float) font_info.unit_per_em;

    // The point size is defined as 1/72 of the Anglo-Saxon inch (25.4 mm):
    // It is approximately 0.0139 inch or 352.8 um.
    return c1 * std::abs(prop.size_in_mm) / 0.3528f;
}

ImFont *FontManager::create_imgui_font(size_t index, const std::string &text)
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

    const auto  &cn = font_prop.collection_number;
    unsigned int font_index = (cn.has_value()) ? *cn : 0;
    const auto  &font_info  = font_file.infos[font_index];
    if (font_prop.char_gap.has_value()) {
        float coef = font_size / (double) font_info.unit_per_em;
        font_config.GlyphExtraSpacing.x = coef * (*font_prop.char_gap);
    }
    if (font_prop.line_gap.has_value()) {
        float coef = font_size / (double) font_info.unit_per_em;
        font_config.GlyphExtraSpacing.y = coef * (*font_prop.line_gap);
    }

    font_config.FontDataOwnedByAtlas = false;

    const std::vector<unsigned char> &buffer = *font_file.data;
    ImFont * font = m_imgui_font_atlas.AddFontFromMemoryTTF(
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
    assert(font == m_imgui_font_atlas.Fonts.back());
    item.imgui_font_index = m_imgui_font_atlas.Fonts.size() - 1;
    assert(font->IsLoaded());
    return font;
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
