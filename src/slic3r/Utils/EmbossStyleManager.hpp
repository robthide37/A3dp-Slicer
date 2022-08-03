#ifndef slic3r_EmbossStyleManager_hpp_
#define slic3r_EmbossStyleManager_hpp_

#include <memory>
#include <optional>
#include <string>
#include <imgui/imgui.h>
#include <GL/glew.h>
#include <libslic3r/Emboss.hpp>
#include <libslic3r/AppConfig.hpp>

class wxFont; 

namespace Slic3r::GUI {
/// <summary>
/// GUI list of loaded fonts
/// Keep pointer to ImGui font pointers
/// Keep file data of TTF files
/// Cache wx font objects
/// </summary>
class EmbossStyleManager
{
    friend class CreateFontStyleImagesJob; // access to StyleImagesData

public:
    EmbossStyleManager(const ImWchar *language_glyph_range);
    ~EmbossStyleManager();

    /// <summary>
    /// Load font style list from config
    /// Also select actual activ font
    /// </summary>
    /// <param name="cfg">Application configuration loaded from file "PrusaSlicer.ini"</param>
    /// <param name="default_font_list">Used when list is not loadable from config</param>
    void init(const AppConfig *cfg, const FontList &default_font_list);
    
    /// <summary>
    /// Write font list into AppConfig
    /// </summary>
    /// <param name="cfg">Stor into this configuration</param>
    /// <param name="item_to_store">Configuration</param>
    bool store_font_list_to_app_config(AppConfig *cfg);

    /// <summary>
    /// Append actual style to style list and store
    /// </summary>
    /// <param name="name">New name for style</param>
    void store_style(const std::string& name);

    /// <summary>
    /// Change order of style item in m_style_items.
    /// Fix selected font index when (i1 || i2) == m_font_selected 
    /// </summary>
    /// <param name="i1">First index to m_style_items</param>
    /// <param name="i2">Second index to m_style_items</param>
    void swap(size_t i1, size_t i2);

    /// <summary>
    /// Track using of swap between saves
    /// </summary>
    /// <returns>True when swap was call after save otherwise false</returns>
    bool is_style_order_changed() const;

    /// <summary>
    /// Check that actual selected style is same as activ style stored in "PrusaSlicer.ini"
    /// </summary>
    /// <returns>True when actual selection is not stored otherwise False</returns>
    bool is_activ_style_changed() const;

    /// <summary>
    /// Remove style from m_style_items.
    /// Fix selected font index when index is under m_font_selected
    /// </summary>
    /// <param name="index">Index of style to be removed</param>
    void erase(size_t index);

    /// <summary>
    /// Rename actual selected font item
    /// </summary>
    /// <param name="name">New name</param>
    void rename(const std::string &name);

    /// <summary>
    /// Actual wx font was changed
    /// Clear caches
    /// </summary>
    /// <param name="font_file">font file created by WxFontUtils::create_font_file(wx_font)</param>
    bool wx_font_changed(std::unique_ptr<Emboss::FontFile> font_file = nullptr);
    
    /// <summary>
    /// Change active font
    /// When font not loaded roll back activ font
    /// </summary>
    /// <param name="font_index">New font index(from m_style_items range)</param>
    /// <returns>True on succes. False on fail load font</returns>
    bool load_font(size_t font_index);
    // load font style not stored in list
    bool load_font(const FontItem &fi);
    // fastering load font on index by wxFont, ignore type and descriptor
    bool load_font(const FontItem &fi, const wxFont &font);
    
    // clear actual selected glyphs cache
    void clear_glyphs_cache();

    // remove cached imgui font for actual selected font
    void clear_imgui_font();

    // erase font when not possible to load
    // used at initialize phaze - fonts could be modified in appConfig file by user
    bool load_first_valid_font();

    // getter on stored fontItem
    const FontItem *get_stored_font_item() const;

    // getter on stored wxFont
    const std::optional<wxFont> &get_stored_wx_font() const;

    // getter on active font item for access to font property
    const FontItem &get_font_item() const;
    FontItem &get_font_item();

    // getter on active font property
    const FontProp &get_font_prop() const;
    FontProp &get_font_prop();

    // getter on activ wx font
    const std::optional<wxFont> &get_wx_font() const;
    std::optional<wxFont> &get_wx_font();

    bool exist_stored_style() const;
    size_t get_style_index() const;

    // getter on font file with cache
    Emboss::FontFileWithCache &get_font_file_with_cache();

    // Getter for cached trucated name for style list selector
    std::string &get_truncated_name();
        
    // setter of font for actual selection
    bool set_wx_font(const wxFont &wx_font);


    // Getter on acitve font pointer for imgui
    // Initialize imgui font(generate texture) when doesn't exist yet.
    // Extend font atlas when not in glyph range
    ImFont *get_imgui_font();
    // initialize font range by unique symbols in text
    ImFont *create_imgui_font(const std::string& text);
    
    // init truncated names of styles
    void init_trunc_names(float max_width);

    /// <summary>
    /// Initialization texture with rendered font style
    /// </summary>
    /// <param name="max_size">Maximal width and height of one style texture</param>
    /// <param name="text">Text to render by style</param>
    void init_style_images(const Vec2i& max_size, const std::string &text);
    void free_style_images();
    
    struct Item;
    // access to all managed font styles
    const std::vector<Item> &get_styles() const;

    /// <summary>
    /// Describe image in GPU to show settings of style
    /// </summary>
    struct StyleImage
    {
        void* texture_id = 0; // GLuint
        BoundingBox bounding_box;
        ImVec2 tex_size, uv0, uv1;
        Point  offset     = Point(0, 0);
        StyleImage()      = default;
    };

    /// <summary>
    /// All connected with one style 
    /// keep temporary data and caches for style
    /// </summary>
    struct Item
    {
        // define font, style and other property of text
        FontItem font_item;

        // cache for view font name with maximal width in imgui
        std::string truncated_name; 

        // visualization of style
        std::optional<StyleImage> image;
    };

    // check if exist selected font style in manager
    bool is_activ_font();

    // Limits for imgui loaded font size
    // Value out of limits is crop
    static float min_imgui_font_size;
    static float max_imgui_font_size;
    static float get_imgui_font_size(const FontProp& prop, const Emboss::FontFile& file);

private:
    /// <summary>
    /// Cache data from style to reduce amount of:
    /// 1) loading font from file
    /// 2) Create atlas of symbols for imgui
    /// 3) Keep loaded(and modified by style) glyphs from font
    /// </summary>
    struct StyleCache
    {
        // share font file data with emboss job thread
        Emboss::FontFileWithCache font_file = {};

        // must live same as imgui_font inside of atlas
        ImVector<ImWchar> ranges = {};

        // Keep only actual style in atlas
        ImFontAtlas atlas = {};

        // wx widget font
        std::optional<wxFont> wx_font = {};

        // cache for view font name with maximal width in imgui
        std::string truncated_name; 

        // actual used font item
        FontItem font_item = {};

        // cache for stored wx font to not create every frame
        std::optional<wxFont> stored_wx_font;

        // index into m_style_items
        size_t font_index = std::numeric_limits<size_t>::max();

    } m_style_cache;
            
    // extend actual imgui font when exist unknown char in text
    // NOTE: imgui_font has to be unused
    // return true when extend range otherwise FALSE
    ImFont *extend_imgui_font_range(size_t font_index, const std::string &text);

    void make_unique_name(std::string &name);

    // Privat member
    std::vector<Item> m_style_items;
    bool m_change_order = false;
    size_t m_stored_activ_index;

    /// <summary>
    /// Keep data needed to create Font Style Images in Job
    /// </summary>
    struct StyleImagesData
    {
        struct Item
        {
            Emboss::FontFileWithCache font;
            std::string               text;
            FontProp                  prop;
        };
        using Items = std::vector<Item>;

        // Keep styles to render
        Items styles;
        // Maximal width and height in pixels of image
        Vec2i max_size;
        // Text to render
        std::string text;

        /// <summary>
        /// Result of job
        /// </summary>
        struct StyleImages
        {
            // vector of inputs
            StyleImagesData::Items styles;
            // job output
            std::vector<EmbossStyleManager::StyleImage> images;
        };

        // place to store result in main thread in Finalize
        std::shared_ptr<StyleImages> result;
    };
    std::shared_ptr<StyleImagesData::StyleImages> m_temp_style_images;
    bool m_exist_style_images;

    // store all font GLImages
    //ImFontAtlas    m_imgui_font_atlas;
    const ImWchar *m_imgui_init_glyph_range;
};

} // namespace Slic3r

#endif // slic3r_EmbossStyleManager_hpp_
