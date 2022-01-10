#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/GLTexture.hpp"
#include "slic3r/Utils/RaycastManager.hpp"

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
class EmbossJob;
class MeshRaycaster;

class GLGizmoEmboss : public GLGizmoBase
{  
    friend EmbossJob;
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
    void set_default_configuration();
    TriangleMesh create_default_mesh();
    TriangleMesh create_mesh();

    /// <summary>
    /// Create mesh from text
    /// </summary>
    /// <param name="text">Text to convert on mesh</param>
    /// <param name="font">Define shape of characters. 
    /// NOTE: Can't be const cache glyphs</param>
    /// <param name="font_prop">Property of font</param>
    /// <returns>Triangle mesh model</returns>
    static TriangleMesh create_mesh(const char *    text,
                                    Emboss::Font &  font,
                                    const FontProp &font_prop);

    void check_selection();
    // more general function --> move to select
    ModelVolume *get_selected_volume();
    static ModelVolume *get_model_volume(const GLVolume *gl_volume, const ModelObjectPtrs& objects);
    static ModelVolume *get_selected_volume(const Selection &selection, const ModelObjectPtrs& objects);
    // create volume from text - main functionality
    bool process();
    void close();
    void draw_window();
    void draw_font_list();
    void draw_text_input();
    void draw_advanced();
    // process mouse event
    bool on_mouse_for_rotation(const wxMouseEvent &mouse_event);
    bool on_mouse_for_translate(const wxMouseEvent &mouse_event);

    bool load_font();
    // try to set font_index
    bool load_font(size_t font_index);
    bool load_font(const wxFont &font);
    void load_imgui_font();
    void check_imgui_font_range();

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
        int    count_line_of_text            = 6;
        // limits for font size inside gizmo
        // When out of limits no change in size will appear in text input
        int min_imgui_font_size = 18;
        int max_imgui_font_size = 60;

        bool draw_advanced = false;

        // setted only when wanted to use - not all the time
        std::optional<ImVec2> offset;

        // Zero means it is calculated in init function
        ImVec2 minimal_window_size = ImVec2(0, 0);
        ImVec2 minimal_window_size_with_advance = ImVec2(0, 0);
        float advanced_input_width    = 0.f;
        float combo_font_width        = 0.f;
        float rename_pos_x            = 0.f;
        float delete_pos_x            = 0.f;
        float max_font_name_width     = 0.f;
        float icon_width              = 0.f;
        float icon_width_with_spacing = 0.f;
        ImVec2 text_size;
        GuiCfg() = default;
    };
    std::optional<GuiCfg> m_gui_cfg;

    FontList m_font_list;    
    size_t   m_font_selected;// index to m_font_list

    // to share data with job thread
    std::shared_ptr<Emboss::Font> m_font;
    std::string m_text;

    // actual volume
    ModelVolume    *m_volume; 

    // Rotation gizmo
    GLGizmoRotate m_rotate_gizmo;

    // TODO: it should be accessible by other gizmo too.
    // May be move to plater?
    RaycastManager m_raycast_manager;

    // Only when drag text object it stores world position
    std::optional<Transform3d> m_temp_transformation;

    // initialize when GL is accessible
    bool m_is_initialized;

    // imgui font
    ImFontAtlas m_imgui_font_atlas;
    // must live same as font in atlas
    ImVector<ImWchar> m_imgui_font_ranges;

    // drawing icons
    GLTexture m_icons_texture;
    bool init_icons();
    enum class IconType: unsigned { rename = 0, erase /*1*/};
    enum class IconState: unsigned { activable = 0, hovered /*1*/, disabled /*2*/};
    void draw_icon(IconType icon, IconState state);
    bool draw_button(IconType icon, bool disable = false);

    // load / store appConfig
    void load_font_list_from_app_config();
    void store_font_list_to_app_config() const;
    void store_font_item_to_app_config() const;

    // call after functions to work outside of drawing
    static void create_emboss_object(TriangleMesh &&   mesh,
                                     std::string       name,
                                     TextConfiguration cfg);
    static void create_emboss_volume(TriangleMesh &&   mesh,
                                     Transform3d       transformation,
                                     std::string       name,
                                     TextConfiguration cfg,
                                     ModelVolumeType   type,
                                     size_t            object_idx);
    static void update_emboss_volume(TriangleMesh &&          mesh,
                                     const std::string &      name,
                                     const TextConfiguration &cfg,
                                     ModelVolume *            volume);

    // only temporary solution
    static const std::string M_ICON_FILENAME;

public:
    // TODO: move to file utils
    static std::string get_file_name(const std::string &file_path);

};

} // namespace Slic3r::GUI

#endif // slic3r_GLGizmoEmboss_hpp_
