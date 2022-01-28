#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/GLTexture.hpp"
#include "slic3r/Utils/RaycastManager.hpp"
#include "slic3r/Utils/FontManager.hpp"

#include "admesh/stl.h" // indexed_triangle_set
#include <optional>
#include <memory>
#include <mutex>
#include <thread>

#include "libslic3r/Emboss.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Model.hpp"

#include <imgui/imgui.h>

class wxFont;
namespace Slic3r{
    class AppConfig;
    class GLVolume;
}

namespace Slic3r::GUI {
class MeshRaycaster;

class GLGizmoEmboss : public GLGizmoBase
{
public:
    GLGizmoEmboss(GLCanvas3D& parent);

    /// <summary>
    /// Create new embossed text volume by type on position of mouse
    /// </summary>
    /// <param name="volume_type">Object part / Negative volume / Modifier</param>
    /// <param name="mouse_pos">Define position of new volume</param>
    void create_volume(ModelVolumeType volume_type, const Vec2d &mouse_pos = Vec2d(-1,-1));

    /// <summary>
    /// Move window for edit emboss text near to embossed object
    /// NOTE: embossed object must be selected
    /// </summary>
    void set_fine_position();
        
    /// <summary>
    /// Rotate by text on dragging rotate grabers
    /// </summary>
    /// <param name="mouse_event">Information about mouse</param>
    /// <returns>Propagete normaly return false.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;
protected:
    bool on_init() override;
    std::string on_get_name() const override;
    void on_render() override;
    void on_render_for_picking() override;    
    void on_render_input_window(float x, float y, float bottom_limit) override;
    bool on_is_activable() const override { return true; }
    bool on_is_selectable() const override { return false; }
    void on_set_state() override;    

    void on_set_hover_id() override{ m_rotate_gizmo.set_hover_id(m_hover_id); }
    void on_enable_grabber(unsigned int id) override { m_rotate_gizmo.enable_grabber(0); }
    void on_disable_grabber(unsigned int id) override { m_rotate_gizmo.disable_grabber(0); }
    void on_update(const UpdateData &data) override { m_rotate_gizmo.update(data); }
    void on_start_dragging() override;
    void on_stop_dragging() override;

private:
    void initialize();
    void set_default_text();

    void check_selection();
    // more general function --> move to select
    ModelVolume *get_selected_volume();
    static ModelVolume *get_model_volume(const GLVolume *gl_volume, const ModelObjectPtrs& objects);
    static ModelVolume *get_selected_volume(const Selection &selection, const ModelObjectPtrs& objects);
    // create volume from text - main functionality
    bool process();
    void close();
    void draw_window();
    void draw_text_input();
    void draw_model_type();
    void draw_style_list();
    void draw_rename_style(bool start_rename);
    void draw_font_list();
    void draw_style_edit();
    bool italic_button();
    bool bold_button();
    void draw_advanced();

    void set_minimal_window_size(bool is_edit_style, bool is_advance_edit_style);
    const ImVec2 &get_minimal_window_size() const;
    // process mouse event
    bool on_mouse_for_rotation(const wxMouseEvent &mouse_event);
    bool on_mouse_for_translate(const wxMouseEvent &mouse_event);

    bool choose_font_by_wxdialog();
    bool choose_true_type_file();
    bool choose_svg_file();

    // Create object described how to make a Volume
    TextConfiguration create_configuration();
    bool load_configuration(ModelVolume *volume);

    // Create notification when unknown font type is used
    bool m_exist_notification;
    void create_notification_not_valid_font(const TextConfiguration& tc);
    void remove_notification_not_valid_font();

    std::string create_volume_name();

    // This configs holds GUI layout size given by translated texts.
    // etc. When language changes, GUI is recreated and this class constructed again,
    // so the change takes effect. (info by GLGizmoFdmSupports.hpp)
    struct GuiCfg
    {
        size_t max_count_char_in_volume_name = 20;
        // Zero means it is calculated in init function
        ImVec2 minimal_window_size              = ImVec2(0, 0);
        ImVec2 minimal_window_size_with_edit    = ImVec2(0, 0);
        ImVec2 minimal_window_size_with_advance = ImVec2(0, 0);
        float advanced_input_width    = 0.f;
        float combo_font_width        = 0.f;
        float delete_pos_x            = 0.f;
        float max_font_name_width     = 0.f;
        float icon_width              = 0.f;
        
        float min_style_image_height = 0.f;
        int   max_style_image_width   = 0.f;

        float style_edit_text_width   = 0.f;

        ImVec2 text_size;

        // Only translations needed for calc GUI size
        struct Translations
        {
            std::string font;
            std::string size;
            std::string depth;
        };
        Translations translations;
        GuiCfg() = default;
    };
    std::optional<const GuiCfg> m_gui_cfg;
    // setted only when wanted to use - not all the time
    std::optional<ImVec2> m_set_window_offset;
    bool m_is_edit_style = false;
    bool m_is_advanced_edit_style = false;

    FontManager m_font_manager;

    //FontList m_font_list;    
    //size_t   m_font_selected;// index to m_font_list
    
    std::string m_text;

    // actual volume
    ModelVolume *m_volume; 

    // Rotation gizmo
    GLGizmoRotate m_rotate_gizmo;

    // TODO: it should be accessible by other gizmo too.
    // May be move to plater?
    RaycastManager m_raycast_manager;

    // Only when drag text object it stores world position
    std::optional<Transform3d> m_temp_transformation;

    // initialize when GL is accessible
    bool m_is_initialized;

    // drawing icons
    GLTexture m_icons_texture;
    bool init_icons();
    enum class IconType : unsigned {
        rename = 0,
        erase,
        duplicate,
        save,
        undo,
        italic,
        unitalic,
        bold,
        unbold,
        system_selector,
        open_file,
        _count /* automatic calc of icon size */
    };
    enum class IconState: unsigned { activable = 0, hovered /*1*/, disabled /*2*/};
    void draw_icon(IconType icon, IconState state);
    bool draw_button(IconType icon, bool disable = false);

    // load / store appConfig
    static FontList load_font_list_from_app_config(const AppConfig *cfg);
    void store_font_list_to_app_config() const;
    void store_font_item_to_app_config() const;

    // only temporary solution
    static const std::string M_ICON_FILENAME;

public:
    // TODO: move to file utils
    static std::string get_file_name(const std::string &file_path);

};

} // namespace Slic3r::GUI

#endif // slic3r_GLGizmoEmboss_hpp_
