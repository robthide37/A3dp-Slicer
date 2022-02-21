#ifndef slic3r_GUI_Tags_hpp_
#define slic3r_GUI_Tags_hpp_

#include "libslic3r/Config.hpp"

#include <wx/checklst.h>
#include <wx/combo.h>

#include "wxExtensions.hpp"

namespace Slic3r {
namespace GUI {


// ----------------------------------------------------------------------------
// ModeButton
// ----------------------------------------------------------------------------

class ModeButton : public ScalableButton
{
public:
    ModeButton(
        wxWindow*           parent,
        wxWindowID          id,
        const std::string&  icon_name = "",
        const wxString&     mode = wxEmptyString,
        const wxSize&       size = wxDefaultSize,
        const wxPoint&      pos = wxDefaultPosition);

    ModeButton(
        wxWindow*           parent,
        const wxString&     mode,
        const std::string&  icon_name,
        int                 px_cnt);

    ModeButton(wxWindow* parent,
        const wxString&     mode,
        wxBitmap*           bitmap,
        int                 px_cnt);

    ~ModeButton() {}

    void Init(const wxString& mode);

    void    OnButton(wxCommandEvent& event);
    void    OnEnterBtn(wxMouseEvent& event) { focus_button(true); event.Skip(); }
    void    OnLeaveBtn(wxMouseEvent& event) { focus_button(m_is_selected); event.Skip(); }

    void    SetState(const bool state);
    bool    is_selected() { return m_is_selected; }

protected:
    void    focus_button(const bool focus);

private:
    bool        m_is_selected = false;

    wxString    m_tt_selected;
    wxString    m_tt_focused;
};



// ----------------------------------------------------------------------------
// ModeSizer
// ----------------------------------------------------------------------------

class ModeSizer : public wxFlexGridSizer
{
public:
    ModeSizer( wxWindow *parent, int hgap, int max_col);
    ~ModeSizer() {}

    void SetMode(const ConfigOptionMode mode);

    void set_items_flag(int flag);
    void set_items_border(int border);

    void msw_rescale();
    const std::vector<ModeButton*>& get_btns() { return m_mode_btns; }

private:
    std::vector<ConfigOptionMode> m_bt_mode;
    std::vector<ModeButton*> m_mode_btns;
    wxWindow*                m_parent {nullptr};
    double                   m_hgap_unscaled;
};

} // namespace slic3r
} // namespace GUI
#endif // slic3r_GUI_Tags_hpp_
