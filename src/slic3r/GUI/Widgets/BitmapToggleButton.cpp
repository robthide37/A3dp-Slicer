#include "BitmapToggleButton.hpp"

#include <wx/settings.h>

BitmapToggleButton::BitmapToggleButton(wxWindow* parent, const wxString& label, wxWindowID id)
{
    if (label.IsEmpty())
        wxBitmapToggleButton::Create(parent, id, wxNullBitmap, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE | wxBU_EXACTFIT);
    else {
#ifdef __linux__
        wxSize def_size = wxSize(parent->GetTextExtent(label).GetX() + 20, 20);
#else
        wxSize def_size = wxDefaultSize;
#endif
        // Call Create() from wxToggleButton instead of wxBitmapToggleButton to allow add Label text under Linux
        wxToggleButton::Create(parent, id, label, wxDefaultPosition, def_size, wxBORDER_NONE | wxBU_EXACTFIT);
    }

#ifdef __WXMSW__
	if (parent) {
		SetBackgroundColour(parent->GetBackgroundColour());
		SetForegroundColour(parent->GetForegroundColour());
	}
#elif __linux__
    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW));
#endif

    Bind(wxEVT_TOGGLEBUTTON, [this](auto& e) {
	    update();

	    wxCommandEvent evt(wxEVT_CHECKBOX);
	    evt.SetInt(int(GetValue()));
	    wxPostEvent(this, evt);

	    e.Skip();
	});
}

void BitmapToggleButton::update_size()
{
#ifdef __linux__
    wxSize bmp_sz = GetBitmap().GetSize();
    wxSize sz = GetSize();
	if (GetLabel().IsEmpty())
        SetSize(bmp_sz);
    else
        SetSize(sz.x, bmp_sz.y);
#else
    wxSize best_sz = GetBestSize();
    SetSize(best_sz);
#endif
}
