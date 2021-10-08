#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "slic3r/GUI/GLTexture.hpp"

#include "admesh/stl.h" // indexed_triangle_set
#include <optional>
#include <memory>

#include "libslic3r/Emboss.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Model.hpp"

namespace Slic3r::GUI {
class GLGizmoEmboss : public GLGizmoBase
{    
public:
    GLGizmoEmboss(GLCanvas3D& parent);
    virtual ~GLGizmoEmboss();

    void set_volume_type(ModelVolumeType volume_type) { m_volume_type = volume_type; }
protected:
    virtual bool on_init() override;
    virtual std::string on_get_name() const override;
    virtual void on_render() override;
    virtual void on_render_for_picking() override;    
    virtual void on_render_input_window(float x, float y, float bottom_limit) override;
    virtual bool on_is_activable() const override { return true; }
    virtual bool on_is_selectable() const override { return false; }
    virtual void on_set_state() override;    

private:
    void initialize();
    void load_font_list();
    void store_font_list();
    static FontList create_default_font_list();
    void set_default_configuration();
    void check_selection();
    // more general function --> move to select
    ModelVolume *get_selected_volume();
    static ModelVolume *get_selected_volume(const Selection &selection, const ModelObjectPtrs objects);
    // create volume from text - main functionality
    bool process();
    bool add_volume(const std::string& name, indexed_triangle_set& its);
    void close();
    void draw_window();
    void draw_font_list();
    void draw_advanced();

    bool load_font();
    // try to set font_index
    bool load_font(size_t font_index);

    bool choose_font_by_wxdialog();
    bool choose_true_type_file();
    bool choose_svg_file();

    // Create object described how to make a Volume
    TextConfiguration create_configuration();
    bool load_configuration(ModelVolume *volume);

    // Create notification when unknown font type is used
    bool notify_unknown_font_type(ModelVolume *volume);
    void notify_cant_load_font(const FontItem &font_item);

    // TODO: Move to imgui utils
    std::string imgui_trunc(const std::string &text, float width);

    std::string create_volume_name();

    // This configs holds GUI layout size given by translated texts.
    // etc. When language changes, GUI is recreated and this class constructed again,
    // so the change takes effect. (info by GLGizmoFdmSupports.hpp)
    struct GuiCfg
    {
        size_t max_count_char_in_volume_name = 20;
        // Zero means it is calculated in init function
        float combo_font_width = 0.f;
        float rename_pos_x = 0.f;
        float delete_pos_x = 0.f;
        float max_font_name_width = 0.f;
        float icon_width = 0.;
        float icon_width_with_spacing = 0.;
        GuiCfg() = default;
    };
    std::optional<GuiCfg> m_gui_cfg;

    FontList m_font_list;    
    size_t   m_font_selected;// index to m_font_list

    std::optional<Emboss::Font> m_font;

    std::string m_text;

    FontProp m_font_prop;

    // text position
    struct Orientation
    {
        Vec3f origin = Vec3f(0.f, 0.f, 0.f);
        Vec3f normal = Vec3f(0.f, 0.f, 1.f);
        Vec3f up     = Vec3f(0.f, 1.f, 0.f);
        Orientation() = default;
    };
    Orientation m_orientation;

    // actual volume
    ModelVolume    *m_volume; 

    // Only for new created volume
    ModelVolumeType m_volume_type; 

    // initialize when GL is accessible
    bool m_is_initialized;

    // drawing icons
    GLTexture m_icons_texture;
    bool init_icons();
    enum class IconType: unsigned { rename = 0, erase /*1*/};
    enum class IconState: unsigned { activable = 0, hovered /*1*/, disabled /*2*/};
    void draw_icon(IconType icon, IconState state);
    bool draw_button(IconType icon, bool disable = false);

    static const std::string M_APP_CFG_FONT_LIST;
    // only temporary solution
    static const std::string M_ICON_FILENAME;
};

} // namespace Slic3r::GUI

#endif // slic3r_GLGizmoEmboss_hpp_
