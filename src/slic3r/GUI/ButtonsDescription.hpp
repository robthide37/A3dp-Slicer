#ifndef slic3r_ButtonsDescription_hpp
#define slic3r_ButtonsDescription_hpp

#include <wx/dialog.h>
#include <vector>

#include <wx/bmpbndl.h>

#include "BitmapComboBox.hpp"

class ScalableBitmap;
class wxColourPickerCtrl;

namespace Slic3r {
namespace GUI {

class BitmapCache;

// ---------------------------------
// ***  PaletteComboBox  ***
// ---------------------------------

// BitmapComboBox used to palets list in GUI Preferences
class ModePaletteComboBox : public BitmapComboBox
{
public:
    ModePaletteComboBox(wxWindow* parent);
	~ModePaletteComboBox() = default;

	void UpdateSelection(const std::vector<wxColour>& palette_in);

protected:
    // Caching bitmaps for the all bitmaps, used in preset comboboxes
    static BitmapCache&		bitmap_cache();
    wxBitmapBundle*			get_bmp( const std::vector<std::string>& palette);
};


class ButtonsDescription : public wxDialog
{
	wxColourPickerCtrl* sys_colour{ nullptr };
	wxColourPickerCtrl* mod_colour{ nullptr };

	wxColourPickerCtrl* simple    { nullptr };
	wxColourPickerCtrl* advanced  { nullptr };
	wxColourPickerCtrl* expert    { nullptr };

	std::vector<wxColour> mode_palette;
public:
	struct Entry {
		Entry(ScalableBitmap *bitmap, const std::string &symbol, const std::string &explanation) : bitmap(bitmap), symbol(symbol), explanation(explanation) {}

		ScalableBitmap *bitmap;
		std::string     symbol;
		std::string   	explanation;
	};

	ButtonsDescription(wxWindow* parent, const std::vector<Entry> &entries);
	~ButtonsDescription() {}

	static void FillSizerWithTextColorDescriptions(wxSizer* sizer, wxWindow* parent, wxColourPickerCtrl** sys_colour, wxColourPickerCtrl** mod_colour);
	static void FillSizerWithModeColorDescriptions(wxSizer* sizer, wxWindow* parent, 
		                                            std::vector<wxColourPickerCtrl**> clr_pickers, 
		                                            std::vector<wxColour>& mode_palette);

private:
	std::vector<Entry> m_entries;
};

} // GUI
} // Slic3r


#endif 

