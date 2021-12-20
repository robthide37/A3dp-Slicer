#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "GLGizmoRotate.hpp"
#include "slic3r/GUI/GLTexture.hpp"

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

    void create_volume(ModelVolumeType volume_type, const Vec2d &mouse_pos = Vec2d(-1,-1));
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
    CommonGizmosDataID on_get_requirements() const override;
    void on_start_dragging() override;
    void on_stop_dragging() override;

private:
    void initialize();
    static FontList create_default_font_list();
    void set_default_configuration();
    void check_selection();
    // more general function --> move to select
    ModelVolume *get_selected_volume();
    static ModelVolume *get_model_volume(const GLVolume *gl_volume, const ModelObjectPtrs objects);
    static ModelVolume *get_selected_volume(const Selection &selection, const ModelObjectPtrs objects);
    // create volume from text - main functionality
    bool process();
    void close();
    void draw_window();
    void preview_positon();
    void draw_font_list();
    void draw_text_input();
    void draw_advanced();

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

    static Transform3d get_emboss_transformation(const Vec3f &position,
                                                 const Vec3f &emboss_dir);
    std::optional<Transform3d> transform_on_surface(const Vec2d &mouse_pos,
        const std::vector<const MeshRaycaster *> &raycasters,
        const std::vector<Transform3d> &          raycasters_tr
    );

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
    FontProp m_font_prop;
    TriangleMesh m_default_mesh; // when add new text this shape is used

    // actual volume
    ModelVolume    *m_volume; 

    // Rotation gizmo
    GLGizmoRotate m_rotate_gizmo;

    // preview position
    GLModel     m_preview;
    Transform3d m_preview_trmat;

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
    void load_font_list();
    void store_font_list();
    static const std::string APP_CONFIG_FONT_NAME;
    static const std::string APP_CONFIG_FONT_DESCRIPTOR;
    std::string get_app_config_font_section(unsigned index);
    std::optional<FontItem> get_font_item(const std::map<std::string, std::string> &app_cfg_section);
    void set_font_item(AppConfig &cfg, const FontItem &fi, unsigned index);

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
