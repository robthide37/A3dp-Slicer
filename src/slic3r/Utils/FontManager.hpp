#ifndef slic3r_FontManager_hpp_
#define slic3r_FontManager_hpp_

#include <optional>
#include <imgui/imgui.h>
#include "libslic3r/Emboss.hpp"

class wxFont; 

namespace Slic3r::GUI {

/// <summary>
/// GUI list of loaded fonts
/// Keep pointer to ImGui font pointers
/// Keep file data of TTF files
/// </summary>
class FontManager
{
public:
    FontManager(const ImWchar *language_glyph_range);
    ~FontManager();

    void select(size_t index);
    void duplicate(size_t index);
    void erase(size_t index);
    
    // load actual selected font
    bool load_font();
    // try to select and load font_index
    bool load_font(size_t font_index);
    // fastering load font on index by wxFont
    bool load_font(size_t font_index, const wxFont &font);
    
    void load_imgui_font(const std::string &text = "");
    void load_imgui_font(size_t index, const std::string &text);

    // erase font when not possible to load
    bool load_first_valid_font();

    // add font into manager
    void add_font(FontItem font_item);
    // add multiple font into manager
    void add_fonts(FontList font_list);

    // getter on active font file for access to glyphs
    std::shared_ptr<Emboss::FontFile> &get_font_file();

    // getter on active font item for access to font property
    const FontItem &get_font_item() const;
    FontItem &get_font_item();

    // getter on active font property
    const FontProp &get_font_prop() const;
    FontProp &get_font_prop();

    // getter on activ wx font
    const std::optional<wxFont> &get_wx_font() const;
    std::optional<wxFont> &get_wx_font();

    // getter on acitve font pointer for imgui
    // text could extend font atlas when not in glyph range
    ImFont *get_imgui_font(const std::string &text = "");

    // getter on index selected font pointer for imgui
    // text could extend font atlas when not in glyph range
    ImFont *get_imgui_font(size_t item_index, const std::string &text = "");

    // free used memory and font file data
    void free_except_active_font();

    struct Item;
    // access to all managed fonts
    const std::vector<Item> &get_fonts() const;

    std::vector<Item> &get_fonts();
    const Item &get_font() const;
    Item &get_font();
    const Item &get_font(size_t index) const;
    Item &get_font(size_t index);

    struct Item
    {
        FontItem                          font_item;

        // cache for view font name with maximal width in imgui
        std::string                       truncated_name; 

        // share font file data with emboss job thread
        std::shared_ptr<Emboss::FontFile> font_file = nullptr;

        std::optional<size_t> imgui_font_index;

        // must live same as imgui_font inside of atlas
        ImVector<ImWchar> font_ranges;

        // wx widget font
        std::optional<wxFont> wx_font;
    };   

    // TODO: make private
    ImFontAtlas m_imgui_font_atlas;

private:
    // extend actual imgui font when exist unknown char in text
    // NOTE: imgui_font has to be unused
    void extend_imgui_font_range(size_t font_index, const std::string &text);

    // Move to imgui utils
    static bool is_text_in_ranges(const ImFont *font, const std::string &text);
    static bool is_text_in_ranges(const ImWchar *ranges, const std::string &text);
    static bool is_char_in_ranges(const ImWchar *ranges, unsigned int letter);

    bool check_imgui_font_range(ImFont *font, const std::string &text);
    void free_imgui_fonts();

    struct Configuration
    {
        // limits for imgui loaded font
        // Value out of limits is crop
        int min_imgui_font_size = 18;
        int max_imgui_font_size = 60;
    } m_cfg;

    // load actual font by wx font
    bool load_font(const wxFont &font);

    void make_unique_name(std::string &name);

    // Privat member
    std::vector<Item> m_font_list;
    size_t            m_font_selected; // index to m_font_list

    // store all font GLImages
    //ImFontAtlas    m_imgui_font_atlas;
    const ImWchar *m_imgui_init_glyph_range;
};

} // namespace Slic3r

#endif // slic3r_FontManager_hpp_
