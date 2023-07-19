#include "BitmapToggleButton.hpp"

#include <wx/settings.h>

BitmapToggleButton::BitmapToggleButton(wxWindow* parent, const wxString& label, wxWindowID id)
{
#ifdef __linux__
    long style = wxBORDER_NONE | wxBU_EXACTFIT;
	if (label.IsEmpty())
		style = style | wxBU_NOTEXT;
	// Call Create() from wxToggleButton instead of wxBitmapToggleButton to allow add Label text under Linux
	wxToggleButton::Create(parent, id, label, wxDefaultPosition, wxDefaultSize, style);
#else
	wxBitmapToggleButton::Create(parent, id, wxNullBitmap, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE | wxBU_EXACTFIT);
	if (!label.IsEmpty())
		SetLabel(label);
#endif

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
	if (GetLabel().IsEmpty())
		SetSize(GetBitmap().GetSize());
	else
#endif
	    SetSize(GetBestSize());
}
