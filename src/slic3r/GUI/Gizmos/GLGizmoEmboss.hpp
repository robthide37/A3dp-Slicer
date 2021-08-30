#ifndef slic3r_GLGizmoEmboss_hpp_
#define slic3r_GLGizmoEmboss_hpp_

// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code,
// which overrides our localization "L" macro.
#include "GLGizmoBase.hpp"
#include "admesh/stl.h" // indexed_triangle_set
#include <optional>
#include <memory>

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
    void close();

    struct MyFont;
    MyFont createFont(const char * fileName);

    struct Configuration
    {
        bool use_count = false;
        // minimal triangle count
        float    decimate_ratio = 50.f; // in percent
        uint32_t wanted_count   = 0; // initialize by percents

        // maximal quadric error
        float max_error = 1.;

        void fix_count_by_ratio(size_t triangle_count)
        {
            wanted_count = static_cast<uint32_t>(
                std::round(triangle_count * (100.f-decimate_ratio) / 100.f));
        }
    } m_configuration;

    // This configs holds GUI layout size given by translated texts.
    // etc. When language changes, GUI is recreated and this class constructed again,
    // so the change takes effect. (info by GLGizmoFdmSupports.hpp)
    struct GuiCfg
    {

    };
    std::optional<GuiCfg> m_gui_cfg;

    struct MyFont
    {
        std::string name;
        std::string file_path;

        MyFont(const std::string &name, const std::string &file_path)
            : name(name), file_path(file_path)
        {}
    };

    std::vector<MyFont> m_fonts;    
    size_t              m_selected;// index to m_fonts


    size_t m_text_size; // allocated size by m_text
    std::unique_ptr<char[]> m_text;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoEmboss_hpp_
