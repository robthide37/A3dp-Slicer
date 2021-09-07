#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "admesh/stl.h" // indexed_triangle_set
#include <optional>
#include <memory>

#include "libslic3r/Emboss.hpp"

namespace Slic3r {
class ModelVolume;
namespace GUI {

class GLGizmoEmboss : public GLGizmoBase
{    
public:
    GLGizmoEmboss(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id);
    virtual ~GLGizmoEmboss();
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

    };
    std::optional<GuiCfg> m_gui_cfg;

    Emboss::FontList m_font_list;    
    size_t           m_font_selected;// index to m_font_list

    std::optional<Emboss::Font> m_font;

    size_t                  m_text_size;
    std::unique_ptr<char[]> m_text;

    Emboss::FontProp m_font_prop;

    float m_scale;
    float m_emboss;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoEmboss_hpp_
