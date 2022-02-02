#ifndef slic3r_Preferences_hpp_
#define slic3r_Preferences_hpp_

#include "GUI.hpp"
#include "GUI_Utils.hpp"
#include "wxExtensions.hpp"

#include <wx/dialog.h>
#include <wx/timer.h>
#include <vector>
#include <map>

class wxColourPickerCtrl;
class wxBookCtrlBase;

namespace Slic3r {

	enum  NotifyReleaseMode {
		NotifyReleaseAll,
		NotifyReleaseOnly,
		NotifyReleaseNone
	};

namespace GUI {

class ConfigOptionsGroup;
class OG_CustomCtrl;

class PreferencesDialog : public DPIDialog
{
	std::map<std::string, std::string>	m_values;
	std::shared_ptr<ConfigOptionsGroup>	m_optgroup_general;
	std::shared_ptr<ConfigOptionsGroup>	m_optgroup_camera;
	std::shared_ptr<ConfigOptionsGroup>	m_optgroup_gui;
#ifdef _WIN32
	std::shared_ptr<ConfigOptionsGroup>	m_optgroup_dark_mode;
#endif //_WIN32
#if ENABLE_ENVIRONMENT_MAP
	std::shared_ptr<ConfigOptionsGroup>	m_optgroup_render;
#endif // ENABLE_ENVIRONMENT_MAP
	wxSizer*                            m_icon_size_sizer;
	wxColourPickerCtrl*					m_sys_colour {nullptr};
	wxColourPickerCtrl*					m_mod_colour {nullptr};
	wxBookCtrlBase*						tabs {nullptr};

    bool                                isOSX {false};
	bool								m_settings_layout_changed {false};
	bool								m_seq_top_layer_only_changed{ false };
	bool								m_recreate_GUI{false};

public:
	explicit PreferencesDialog(wxWindow* paren);
	~PreferencesDialog() = default;

	bool settings_layout_changed() const { return m_settings_layout_changed; }
	bool seq_top_layer_only_changed() const { return m_seq_top_layer_only_changed; }
	bool recreate_GUI() const { return m_recreate_GUI; }
	void	build();
	void	update_ctrls_alignment();
	void	accept(wxEvent&);
	void	show(const std::string& highlight_option = std::string(), const std::string& tab_name = std::string());

protected:
	void msw_rescale();
	void on_dpi_changed(const wxRect& suggested_rect) override { msw_rescale(); }
	void on_sys_color_changed() override;
    void layout();
    void create_icon_size_slider();
    void create_settings_mode_widget();
    void create_settings_text_color_widget();
	void init_highlighter(const t_config_option_key& opt_key);
	std::vector<ConfigOptionsGroup*> optgroups();

	HighlighterForWx						m_highlighter;
	std::map<std::string, BlinkingBitmap*>	m_blinkers;
};

} // GUI
} // Slic3r


#endif /* slic3r_Preferences_hpp_ */
