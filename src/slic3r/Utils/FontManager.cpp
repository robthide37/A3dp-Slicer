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
    : m_imgui_init_glyph_range(language_glyph_range), 
    m_font_selected(0)
{}

void FontManager::select(size_t index)
{
    if (index < m_font_list.size())
        m_font_selected = index;
}

void FontManager::duplicate(size_t index) {
    if (index >= m_font_list.size()) return;
    Item item = m_font_list[index]; // copy
    make_unique_name(item.font_item.name);

    m_font_list.insert(m_font_list.begin() + index, item);
    // fix selected index
    if (index < m_font_selected) ++m_font_selected;
}

void FontManager::erase(size_t index) {
    if (index >= m_font_list.size()) return;
    m_font_list.erase(m_font_list.begin() + index);
    // fix selected index
    if (index < m_font_selected) --m_font_selected;
}

bool FontManager::load_font(size_t font_index)
{
    if (font_index >= m_font_list.size()) return false;
    std::swap(font_index, m_font_selected);
    bool is_loaded = load_font();
    if (!is_loaded) std::swap(font_index, m_font_selected);
    return is_loaded;
}

bool FontManager::load_font(size_t font_index, const wxFont &font)
{
    if (font_index >= m_font_list.size()) return false;
    std::swap(font_index, m_font_selected);
    bool is_loaded = load_font(font);
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

bool FontManager::load_font()
{
    // next condition may be safely removed
    if (m_font_selected >= m_font_list.size()) return false; 

    Item &item = m_font_list[m_font_selected];
    FontItem &fi = item.font_item;
    if (fi.type == FontItem::Type::file_path) {
        // fill font name after load from .3mf
        if (fi.name.empty())
            fi.name = get_file_name(fi.path);
        std::unique_ptr<Emboss::FontFile> font_ptr = Emboss::load_font(
            fi.path.c_str());
        if (font_ptr == nullptr) return false;
        item.font_file = std::move(font_ptr);
        load_imgui_font();
        return true;
    }
    if (fi.type != WxFontUtils::get_actual_type()) return false;
    std::optional<wxFont> wx_font = WxFontUtils::load_wxFont(fi.path);
    if (!wx_font.has_value()) return false;

    // fill font name after load from .3mf
    if (fi.name.empty())
        fi.name = WxFontUtils::get_human_readable_name(*wx_font);
    return load_font(*wx_font);
}

bool FontManager::load_first_valid_font() {
    // try to load valid font
    m_font_selected     = 0;
    bool is_font_loaded = load_font();
    while (!is_font_loaded && !m_font_list.empty()) {
        // can't load so erase it from list
        m_font_list.erase(m_font_list.begin());
        is_font_loaded = load_font();
    }
    return !m_font_list.empty();
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

std::shared_ptr<Emboss::FontFile> &FontManager::get_font_file()
{
    return m_font_list[m_font_selected].font_file;
}

const FontItem &FontManager::get_font_item() const
{
    return m_font_list[m_font_selected].font_item;
}

FontItem &FontManager::get_font_item()
{
    return m_font_list[m_font_selected].font_item;
}

ImFont *FontManager::get_imgui_font()
{
    return m_font_list[m_font_selected].imgui_font;
}

const std::vector<FontManager::Item> &Slic3r::GUI::FontManager::get_fonts() const
{
    return m_font_list;
}

std::vector<FontManager::Item> &FontManager::get_fonts()
{
    return m_font_list;
}

const FontManager::Item &FontManager::get_font() const
{
    return m_font_list[m_font_selected];
}

FontManager::Item &FontManager::get_font(size_t index)
{
    return m_font_list[index];
}

const FontManager::Item &FontManager::get_font(size_t index) const
{
    return m_font_list[index];
}

void FontManager::make_unique_name(std::string &name) {
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

void FontManager::check_imgui_font_range(const std::string& text)
{
    const ImFont *font = m_imgui_font_atlas.Fonts.front();
    if (!font->IsLoaded()) {
        // when create font no one letter in text was inside font
        // check text again
        load_imgui_font();
        return;
    }
    if (font->ConfigData == nullptr) return;
    const ImWchar *ranges       = font->ConfigData->GlyphRanges;
    auto           is_in_ranges = [ranges](unsigned int letter) -> bool {
        for (const ImWchar *range = ranges; range[0] && range[1]; range += 2) {
            ImWchar from = range[0];
            ImWchar to   = range[1];
            if (from <= letter && letter <= to) return true;
            if (letter < to) return false; // ranges should be sorted
        }
        return false;
    };

    bool exist_unknown = false;
    const char *text_char_ptr = text.c_str();
    while (*text_char_ptr) {
        unsigned int c = 0;
        // UTF-8 to 32-bit character need imgui_internal
        int c_len = ImTextCharFromUtf8(&c, text_char_ptr, NULL);
        text_char_ptr += c_len;
        if (c_len == 0) break;
        if (!is_in_ranges(c)) {
            exist_unknown = true;
            break;
        }
    }
    if (exist_unknown) load_imgui_font();
}

bool FontManager::load_font(const wxFont &font)
{
    auto font_ptr = WxFontUtils::load_font(font);
    if (font_ptr == nullptr) return false;
    m_font_list[m_font_selected].font_file = std::move(font_ptr);
    load_imgui_font();
    return true;
}

void FontManager::load_imgui_font(const std::string &text)
{
    Item &item = m_font_list[m_font_selected];
    if (item.font_file == nullptr) return;

    // TODO: Create glyph range
    ImFontGlyphRangesBuilder builder;
    builder.AddRanges(m_imgui_init_glyph_range);
    if (!text.empty())
        builder.AddText(text.c_str());
    item.font_ranges.clear();

    builder.BuildRanges(&item.font_ranges);
    const FontProp &font_prop = item.font_item.prop;
    int             font_size = static_cast<int>(
        std::round(std::abs(font_prop.size_in_mm / 0.3528)));
    if (font_size < m_cfg.min_imgui_font_size)
        font_size = m_cfg.min_imgui_font_size;
    if (font_size > m_cfg.max_imgui_font_size)
        font_size = m_cfg.max_imgui_font_size;

    ImFontConfig font_config;
    font_config.FontDataOwnedByAtlas = false;
    m_imgui_font_atlas.Flags |= ImFontAtlasFlags_NoMouseCursors |
                                ImFontAtlasFlags_NoPowerOfTwoHeight;
    m_imgui_font_atlas.Clear();

    const std::vector<unsigned char> &buffer = item.font_file->buffer;
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
    item.imgui_font = m_imgui_font_atlas.Fonts.front();
}
