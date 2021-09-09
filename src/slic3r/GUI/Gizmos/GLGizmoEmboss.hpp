#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "GLGizmosCommon.hpp"

#include "admesh/stl.h" // indexed_triangle_set
#include <optional>
#include <memory>

#include "libslic3r/Emboss.hpp"
#include "libslic3r/Point.hpp"

namespace Slic3r {
class ModelVolume;
namespace GUI {

class GLGizmoEmboss : public GLGizmoBase
{    
public:
    GLGizmoEmboss(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);
    virtual ~GLGizmoEmboss();
    // pseudo virtual function, no inheritance
    virtual bool gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down);
protected:
    virtual bool on_init() override;
    virtual std::string on_get_name() const override;
    virtual void on_render() override;
    virtual void on_render_for_picking() override;    
    virtual void on_render_input_window(float x, float y, float bottom_limit) override;
    virtual bool on_is_activable() const override;
    virtual bool on_is_selectable() const override { return true; }
    virtual void on_set_state() override;    

private:
    void process();
    void close();
    void draw_add_button();
    bool load_font();
    void sort_fonts();
    void add_fonts(const Emboss::FontList &font_list);

    // This configs holds GUI layout size given by translated texts.
    // etc. When language changes, GUI is recreated and this class constructed again,
    // so the change takes effect. (info by GLGizmoFdmSupports.hpp)
    struct GuiCfg
    {
        const size_t max_font_name = 20; // count characters
        GuiCfg() = default;
    };
    std::optional<GuiCfg> m_gui_cfg;

    Emboss::FontList m_font_list;    
    size_t           m_font_selected;// index to m_font_list

    std::optional<Emboss::Font> m_font;
    Emboss::Glyphs              m_font_glyph_cache;

    size_t                  m_text_size;
    std::unique_ptr<char[]> m_text;

    Emboss::FontProp m_font_prop;

    // text position
    struct Orientation
    {
        Vec3f origin = Vec3f(0.f, 0.f, 0.f);
        Vec3f normal = Vec3f(0.f, 0.f, 1.f);
        Vec3f up     = Vec3f(0.f, 1.f, 0.f);
        Orientation() = default;
    };
    Orientation m_orientation;

    float m_scale;
    float m_emboss;

    // actual volume
    ModelVolume *m_volume; 
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoEmboss_hpp_
