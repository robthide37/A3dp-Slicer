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
#include <atomic>

#include "libslic3r/Emboss.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TextConfiguration.hpp"

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
    void on_enable_grabber(unsigned int id) override { m_rotate_gizmo.enable_grabber(); }
    void on_disable_grabber(unsigned int id) override { m_rotate_gizmo.disable_grabber(); }
    void on_start_dragging() override;
    void on_stop_dragging() override;
    void on_dragging(const UpdateData &data) override;    

    /// <summary>
    /// Rotate by text on dragging rotate grabers
    /// </summary>
    /// <param name="mouse_event">Information about mouse</param>
    /// <returns>Propagete normaly return false.</returns>
    bool on_mouse(const wxMouseEvent &mouse_event) override;

    bool wants_enter_leave_snapshots() const override { return true; }
    std::string get_gizmo_entering_text() const override { return _u8L("Enter emboss gizmo"); }
    std::string get_gizmo_leaving_text() const override { return _u8L("Leave emboss gizmo"); }
    std::string get_action_snapshot_name() override { return _u8L("Embossing actions"); }
private:
    void initialize();
    static FontList create_default_font_list();
    // Could exist systems without installed font so last chance is used own file
    static FontItem create_default_font();
    void set_default_text();

    bool start_volume_creation(ModelVolumeType volume_type, const Vec2d &screen_coor);

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
    void draw_delete_style_button();
    void draw_undo_style_button(bool is_stored, bool is_changed);
    void draw_revert_all_styles_button();
    void draw_rename_style(bool start_rename);
    void draw_font_list();
    void draw_style_edit();
    bool italic_button();
    bool bold_button();
    void draw_advanced();

    bool select_facename(const wxString& facename);
    void init_face_names();

    void do_translate(const Vec3d& relative_move);
    void do_rotate(float relative_z_angle);

    /// <summary>
    /// Choose valid source Volume to project on(cut surface from).
    /// </summary>
    /// <returns>ModelVolume to project on</returns>
    const ModelVolume *get_volume_to_cut_surface_from();

    /// <summary>
    /// Reversible input float with option to restor default value
    /// TODO: make more general, static and move to ImGuiWrapper 
    /// </summary>
    /// <returns>True when value changed otherwise FALSE.</returns>
    bool rev_input(const std::string &name, float &value, float *default_value, 
        const std::string &undo_tooltip, float step, float step_fast, const char *format, 
        ImGuiInputTextFlags flags = 0);
    bool rev_checkbox(const std::string &name, bool &value, bool* default_value, const std::string  &undo_tooltip);
    bool rev_slider(const std::string &name, std::optional<int>& value, std::optional<int> *default_value,
        const std::string &undo_tooltip, int v_min, int v_max, const std::string &format, const wxString &tooltip);
    bool rev_slider(const std::string &name, std::optional<float>& value, std::optional<float> *default_value,
        const std::string &undo_tooltip, float v_min, float v_max, const std::string &format, const wxString &tooltip);
    bool rev_slider(const std::string &name, float &value, float *default_value, 
        const std::string &undo_tooltip, float v_min, float v_max, const std::string &format, const wxString &tooltip);
    template<typename T, typename Draw>
    bool revertible(const std::string &name, T &value, T *default_value, bool exist_change, const std::string &undo_tooltip, float undo_offset, Draw draw);

    void set_minimal_window_size(bool is_advance_edit_style);
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
        ImVec2 minimal_window_size_with_advance = ImVec2(0, 0);
        float        input_width                      = 0.f;
        float        delete_pos_x                     = 0.f;
        float        max_font_name_width              = 0.f;
        unsigned int icon_width                       = 0;

        float min_style_image_height = 0.f;
        int   max_style_image_width  = 0.f;

        float style_offset          = 0.f;
        float input_offset          = 0.f;
        float advanced_input_offset = 0.f;

        ImVec2 text_size;

        // Only translations needed for calc GUI size
        struct Translations
        {
            std::string type;
            std::string style;
            std::string font;
            std::string size;
            std::string depth;
            std::string use_surface;

            // advanced
            std::string char_gap;
            std::string line_gap;
            std::string boldness;
            std::string italic;
            std::string surface_distance;
            std::string angle;
            std::string collection;
        };
        Translations translations;
                
        GuiCfg() = default;
    };
    std::optional<const GuiCfg> m_gui_cfg;
    // setted only when wanted to use - not all the time
    std::optional<ImVec2> m_set_window_offset;
    bool m_is_advanced_edit_style = false;

    FontManager m_font_manager;

    // Keep sorted list of loadable face names
    struct Facenames
    {
        bool                  is_init = false;
        std::vector<wxString> names;
        wxFontEncoding        encoding;
    } m_face_names;

    // Track stored values in AppConfig
    std::optional<FontItem> m_stored_font_item;
    std::optional<wxFont> m_stored_wx_font; // cache for stored wx font to not create every frame
    std::map<std::string, FontItem> m_stored_font_items;
    void fill_stored_font_items();
    void select_stored_font_item();

    // Text to emboss
    std::string m_text;
    // True when m_text contain character unknown by selected font
    bool m_text_contain_unknown_glyph = false;

    // cancel for previous update of volume to cancel finalize part
    std::shared_ptr<std::atomic<bool>> m_update_job_cancel;

    // actual volume
    ModelVolume *m_volume; 

    // Rotation gizmo
    GLGizmoRotate m_rotate_gizmo;

    // when draging with text object hold screen offset of cursor from object center
    std::optional<Vec2d> m_dragging_mouse_offset;

    // TODO: it should be accessible by other gizmo too.
    // May be move to plater?
    RaycastManager m_raycast_manager;

    // Only when drag text object it stores world position
    std::optional<Transform3d> m_temp_transformation;

    // drawing icons
    GLTexture m_icons_texture;
    void init_icons();
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
        revert_all,
        // VolumeType icons
        part,
        negative,
        modifier,
        // automatic calc of icon's count
        _count
    };
    enum class IconState: unsigned { activable = 0, hovered /*1*/, disabled /*2*/};
    void draw_icon(IconType icon, IconState state, ImVec2 size = ImVec2(0,0));
    void draw_transparent_icon();
    bool draw_clickable(IconType icon, IconState state, IconType hover_icon, IconState hover_state);
    bool draw_button(IconType icon, bool disable = false);

    // load / store appConfig
    static FontList load_font_list_from_app_config(const AppConfig *cfg, size_t &activ_font_index);
    void store_font_list_to_app_config();
    //void store_font_item_to_app_config() const;

    // only temporary solution
    static const std::string M_ICON_FILENAME;

public:
    /// <summary>
    /// Check if text is last solid part of object
    /// TODO: move to emboss gui utils
    /// </summary>
    /// <param name="text">Model volume of Text</param>
    /// <returns>True when object otherwise False</returns>
    static bool is_text_object(const ModelVolume *text);

    // TODO: move to file utils
    static std::string get_file_name(const std::string &file_path);
};

} // namespace Slic3r::GUI

#endif // slic3r_GLGizmoEmboss_hpp_
