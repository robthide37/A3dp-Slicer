#include "EmbossStyleManager.hpp"
#include <wx/font.h>
#include <GL/glew.h> // Imgui texture
#include <imgui/imgui_internal.h> // ImTextCharFromUtf8
#include "WxFontUtils.hpp"
#include "libslic3r/Utils.hpp" // ScopeGuard

#include "slic3r/GUI/3DScene.hpp" // ::glsafe
#include "slic3r/GUI/Jobs/CreateFontStyleImagesJob.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp" // check of font ranges

#include "slic3r/Utils/EmbossStylesSerializable.hpp"

using namespace Slic3r;
using namespace Slic3r::GUI;

EmbossStyleManager::EmbossStyleManager(const ImWchar *language_glyph_range)
    : m_imgui_init_glyph_range(language_glyph_range)
    , m_exist_style_images(false)
    , m_temp_style_images(nullptr)
    , m_app_config(nullptr)
    , m_last_style_index(std::numeric_limits<size_t>::max())
{}

EmbossStyleManager::~EmbossStyleManager() { 
    clear_imgui_font();
    free_style_images();
    store_styles_to_app_config(false, false);
    if (m_app_config != nullptr &&
        m_last_style_index != std::numeric_limits<size_t>::max())
        EmbossStylesSerializable::store_style_index(*m_app_config, m_last_style_index);
}

void EmbossStyleManager::init(AppConfig *app_config, const EmbossStyles &default_styles)
{
    m_app_config = app_config;
    EmbossStyles styles = (app_config != nullptr) ?
        EmbossStylesSerializable::load_styles(*app_config) :
        default_styles;
    if (styles.empty()) styles = default_styles;
    for (EmbossStyle &style : styles) {
        make_unique_name(style.name);
        m_style_items.push_back({style});
    }

    std::optional<size_t> activ_index_opt = (app_config != nullptr) ? 
        EmbossStylesSerializable::load_style_index(*app_config) : 
        std::optional<size_t>{};

    size_t activ_index = 0;
    if (activ_index_opt.has_value()) activ_index = *activ_index_opt;    
    if (activ_index >= m_style_items.size()) activ_index = 0;
    
    // find valid font item
    if (!load_style(activ_index)) {
        m_style_items.erase(m_style_items.begin() + activ_index);
        activ_index = 0;
        while (m_style_items.empty() || !load_style(activ_index))
            m_style_items.erase(m_style_items.begin());
        // no one style from config is loadable
        if (m_style_items.empty()) {
            // set up default font list
            for (EmbossStyle style : default_styles) {
                make_unique_name(style.name);
                m_style_items.push_back({std::move(style)});
            }
            // try to load first default font
            bool loaded = load_style(activ_index);
            assert(loaded);
        }
    }
}

bool EmbossStyleManager::store_styles_to_app_config(bool use_modification,
                                                    bool store_activ_index)
{
    assert(m_app_config != nullptr);
    if (m_app_config == nullptr) return false;
    if (use_modification) {
        if (exist_stored_style()) {
            // update stored item
            m_style_items[m_style_cache.style_index].style = m_style_cache.style;
        } else {
            // add new into stored list
            EmbossStyle &style = m_style_cache.style;
            make_unique_name(style.name);
            m_style_cache.truncated_name.clear();
            m_style_cache.style_index = m_style_items.size();
            m_style_items.push_back({style});
            m_style_cache.stored_wx_font = m_style_cache.wx_font;
        }
    }
    if (store_activ_index)
        EmbossStylesSerializable::store_style_index(*m_app_config, m_style_cache.style_index);

    EmbossStyles styles;
    styles.reserve(m_style_items.size());
    for (const Item &item : m_style_items) styles.push_back(item.style);
    EmbossStylesSerializable::store_styles(*m_app_config, styles);
    return true;
}

void EmbossStyleManager::add_style(const std::string &name) {
    EmbossStyle& style = m_style_cache.style;
    style.name = name;
    make_unique_name(style.name);
    m_style_cache.style_index = m_style_items.size();
    m_style_cache.stored_wx_font = m_style_cache.wx_font;
    m_style_cache.truncated_name.clear();
    m_style_items.push_back({style});
}

void EmbossStyleManager::swap(size_t i1, size_t i2) {
    if (i1 >= m_style_items.size() || 
        i2 >= m_style_items.size()) return;
    std::swap(m_style_items[i1], m_style_items[i2]);
    // fix selected index
    if (!exist_stored_style()) return;
    if (m_style_cache.style_index == i1) {
        m_style_cache.style_index = i2;
    } else if (m_style_cache.style_index == i2) {
        m_style_cache.style_index = i1;
    }
}

void EmbossStyleManager::discard_style_changes() {
    if (exist_stored_style()) {
        if (load_style(m_style_cache.style_index))
            return; // correct reload style
    } else {
        if(load_style(m_last_style_index))
            return; // correct load last used style
    }

    // try to save situation by load some font
    load_first_valid_font();
}

void EmbossStyleManager::erase(size_t index) {
    if (index >= m_style_items.size()) return;

    // fix selected index
    if (exist_stored_style()) {
        size_t &i = m_style_cache.style_index;
        if (index < i) --i;
        else if (index == i) i = std::numeric_limits<size_t>::max();
    }

    m_style_items.erase(m_style_items.begin() + index);
}

void EmbossStyleManager::rename(const std::string& name) {
    m_style_cache.style.name = name;
    m_style_cache.truncated_name.clear();
    if (exist_stored_style()) { 
        Item &it = m_style_items[m_style_cache.style_index];
        it.style.name = name;
        it.truncated_name.clear();
    }
}

bool EmbossStyleManager::load_style(size_t style_index)
{
    if (style_index >= m_style_items.size()) return false;
    if (!load_style(m_style_items[style_index].style)) return false;
    m_style_cache.style_index    = style_index;
    m_style_cache.stored_wx_font = m_style_cache.wx_font; // copy
    m_last_style_index           = style_index;
    return true;
}

bool EmbossStyleManager::load_style(const EmbossStyle &style) {
    if (style.type == EmbossStyle::Type::file_path) {
        std::unique_ptr<Emboss::FontFile> font_ptr =
            Emboss::create_font_file(style.path.c_str());
        if (font_ptr == nullptr) return false;
        m_style_cache.wx_font = {};
        m_style_cache.font_file = 
            Emboss::FontFileWithCache(std::move(font_ptr));
        m_style_cache.style          = style; // copy
        m_style_cache.style_index    = std::numeric_limits<size_t>::max();
        m_style_cache.stored_wx_font = {};
        return true;
    }
    if (style.type != WxFontUtils::get_actual_type()) return false;
    std::optional<wxFont> wx_font_opt = WxFontUtils::load_wxFont(style.path);
    if (!wx_font_opt.has_value()) return false;
    return load_style(style, *wx_font_opt);
}

bool EmbossStyleManager::load_style(const EmbossStyle &style, const wxFont &font)
{
    if (!set_wx_font(font)) return false;
    m_style_cache.style      = style; // copy
    m_style_cache.style_index     = std::numeric_limits<size_t>::max();
    m_style_cache.stored_wx_font = {};
    m_style_cache.truncated_name.clear();
    return true;
}

bool EmbossStyleManager::is_activ_font() { return m_style_cache.font_file.has_value(); }

bool EmbossStyleManager::load_first_valid_font() {
    while (!m_style_items.empty()) {
        if (load_style(0)) return true;
        // can't load so erase it from list
        m_style_items.erase(m_style_items.begin());
    }
    return false;
}

const EmbossStyle* EmbossStyleManager::get_stored_style() const
{
    if (m_style_cache.style_index >= m_style_items.size()) return nullptr;
    return &m_style_items[m_style_cache.style_index].style;
}

void EmbossStyleManager::clear_glyphs_cache()
{
    Emboss::FontFileWithCache &ff = m_style_cache.font_file;
    if (!ff.has_value()) return;
    ff.cache = std::make_shared<Emboss::Glyphs>();
}

void EmbossStyleManager::clear_imgui_font() { m_style_cache.atlas.Clear(); }

ImFont *EmbossStyleManager::get_imgui_font()
{
    if (!is_activ_font()) return nullptr;
    
    ImVector<ImFont *> &fonts = m_style_cache.atlas.Fonts;
    if (fonts.empty()) return nullptr;

    // check correct index
    int f_size = fonts.size();
    assert(f_size == 1);
    if (f_size != 1) return nullptr;
    ImFont *font = fonts.front();
    if (font == nullptr) return nullptr;
    return font;
}

const std::vector<EmbossStyleManager::Item> &EmbossStyleManager::get_styles() const{ return m_style_items; }

ImFont* EmbossStyleManager::extend_imgui_font_range(size_t index, const std::string& text)
{
    // TODO: start using merge mode
    // ImFontConfig::MergeMode = true;
    return create_imgui_font(text);
}

void EmbossStyleManager::make_unique_name(std::string &name)
{
    auto is_unique = [&](const std::string &name) -> bool {
        for (const Item &it : m_style_items)
            if (it.style.name == name) return false;
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

void EmbossStyleManager::init_trunc_names(float max_width) { 
    for (auto &s : m_style_items)
        if (s.truncated_name.empty())
            s.truncated_name = ImGuiWrapper::trunc(s.style.name, max_width);
}

#include "slic3r/GUI/Jobs/CreateFontStyleImagesJob.hpp"

// for access to worker
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp" 

void EmbossStyleManager::init_style_images(const Vec2i &max_size,
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
            for (EmbossStyleManager::StyleImage &image : m_temp_style_images->images){
                size_t index = &image - &m_temp_style_images->images.front();
                StyleImagesData::Item &style = m_temp_style_images->styles[index];

                // find style in font list and copy to it
                for (auto &it : m_style_items) {
                    if (it.style.name != style.text ||
                        !(it.style.prop == style.prop))
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
    styles.reserve(m_style_items.size());
    for (const Item &item : m_style_items) {
        const EmbossStyle &style = item.style;
        std::optional<wxFont> wx_font_opt = WxFontUtils::load_wxFont(style.path);
        if (!wx_font_opt.has_value()) continue;
        std::unique_ptr<Emboss::FontFile> font_file =
            WxFontUtils::create_font_file(*wx_font_opt);
        if (font_file == nullptr) continue;
        styles.push_back({
            Emboss::FontFileWithCache(std::move(font_file)), 
            style.name,
            style.prop
        });
    }
    auto &worker = wxGetApp().plater()->get_ui_job_worker();
    StyleImagesData data{std::move(styles), max_size, text, m_temp_style_images};
    queue_job(worker, std::make_unique<CreateFontStyleImagesJob>(std::move(data)));
}

void EmbossStyleManager::free_style_images() {
    if (!m_exist_style_images) return;
    if (!is_activ_font()) return;

    GLuint tex_id = 0;
    
    for (Item &it : m_style_items) {
        if (tex_id == 0 && it.image.has_value())
            tex_id = (GLuint)(intptr_t) it.image->texture_id;
        it.image.reset();
    }
    if (tex_id != 0)
        glsafe(::glDeleteTextures(1, &tex_id));
    m_exist_style_images = false;
}

float EmbossStyleManager::min_imgui_font_size = 18.f;
float EmbossStyleManager::max_imgui_font_size = 60.f;
float EmbossStyleManager::get_imgui_font_size(const FontProp         &prop,
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

ImFont *EmbossStyleManager::create_imgui_font(const std::string &text)
{
    auto& ff = m_style_cache.font_file;
    if (!ff.has_value()) return nullptr;
    const Emboss::FontFile &font_file = *ff.font_file;

    // TODO: Create glyph range
    ImFontGlyphRangesBuilder builder;
    builder.AddRanges(m_imgui_init_glyph_range);
    if (!text.empty())
        builder.AddText(text.c_str());

    ImVector<ImWchar> &ranges = m_style_cache.ranges;
    ranges.clear();
    builder.BuildRanges(&ranges);
        
    m_style_cache.atlas.Flags |= ImFontAtlasFlags_NoMouseCursors |
                                ImFontAtlasFlags_NoPowerOfTwoHeight;

    const FontProp &font_prop = m_style_cache.style.prop;
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
    ImFont * font = m_style_cache.atlas.AddFontFromMemoryTTF(
        (void *) buffer.data(), buffer.size(), font_size, &font_config, m_style_cache.ranges.Data);

    unsigned char *pixels;
    int            width, height;
    m_style_cache.atlas.GetTexDataAsAlpha8(&pixels, &width, &height);

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
    m_style_cache.atlas.TexID = (ImTextureID) (intptr_t) font_texture;
    assert(!m_style_cache.atlas.Fonts.empty());
    if (m_style_cache.atlas.Fonts.empty()) return nullptr;
    assert(font == m_style_cache.atlas.Fonts.back());
    if (!font->IsLoaded()) return nullptr;
    assert(font->IsLoaded());
    return font;
}

bool EmbossStyleManager::set_wx_font(const wxFont &wx_font) {
    std::unique_ptr<Emboss::FontFile> font_file = 
        WxFontUtils::create_font_file(wx_font);
    return set_wx_font(wx_font, std::move(font_file));
}

bool EmbossStyleManager::set_wx_font(const wxFont &wx_font, std::unique_ptr<Emboss::FontFile> font_file)
{
    if (font_file == nullptr) return false;
    m_style_cache.wx_font   = wx_font; // copy
    m_style_cache.font_file = 
        Emboss::FontFileWithCache(std::move(font_file));

    EmbossStyle &style = m_style_cache.style;
    style.type = WxFontUtils::get_actual_type();
    // update string path
    style.path = WxFontUtils::store_wxFont(wx_font);
    clear_imgui_font();
    free_style_images();
    return true;
}