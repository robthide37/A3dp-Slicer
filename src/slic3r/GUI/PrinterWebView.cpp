// Webview Panel created by Steve.
#include "PrinterWebView.hpp"

#include "slic3r/GUI/wxExtensions.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "libslic3r_version.h"

#include <wx/sizer.h>
#include <wx/string.h>
#include <wx/toolbar.h>
#include <wx/textdlg.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/bmpcbox.h>
#include <wx/statbox.h>
#include <wx/statbmp.h>
#include <wx/filedlg.h>
#include <wx/dnd.h>
#include <wx/progdlg.h>
#include <wx/wupdlock.h>
#include <wx/numdlg.h>
#include <wx/debug.h>
#include <wx/busyinfo.h>

#include "slic3r/Utils/WebViewHelper.hpp"
#include <wx/webview.h>
#include <wx/msw/webview_edge.h>
#include <wx/msw/webview_ie.h>
#include "slic3r/GUI/PresetComboBoxes.hpp"

// Todo add if mac

#if _MSW_
#include <WebView2.h>
#endif

namespace Slic3r {

namespace GUI {

WebViewPanel::WebViewPanel(wxWindow *parent) : wxPanel(parent, wxID_ANY)
{
    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);
    this->Bind(wxEVT_COMBOBOX, &WebViewPanel::on_select_preset, this);

    m_webView = CreateWebView(this);  // Correct parent is `this`, not `parent`
    
    auto init_combo = [this](DevicePresetComboBox **combo, wxString label, Preset::Type preset_type, bool filament) {
        *combo = new DevicePresetComboBox(this, preset_type);
    };
    
    init_combo(&m_combo_printer, "Printer", Preset::TYPE_PRINTER, false);
    m_combo_printer->update();
    m_combo_printer->SetMaxSize(wxSize(600, 100));
    
    sizer->Add(m_combo_printer, 0, wxEXPAND | wxALL, 5);  // Add padding around the combo box
    sizer->Add(m_webView, 1, wxEXPAND | wxALL, 5);  // Add padding around the web view and make it expand to fill the remaining space

    this->SetSizer(sizer);
}

WebViewPanel::~WebViewPanel(){}

wxWebView* WebViewPanel::CreateWebView(wxWindow* parent) {

    wxWebView *webview = WebView::CreateWebView(parent, "");

    webview->EnableContextMenu(true);
    webview->EnableAccessToDevTools(true);

    return webview;
}

void WebViewPanel::on_select_preset(wxCommandEvent &evt)
{
    DevicePresetComboBox *combo       = static_cast<DevicePresetComboBox *>(evt.GetEventObject());
    Preset::Type          preset_type = combo->get_type();

    // see https://github.com/prusa3d/PrusaSlicer/issues/3889
    // Under OSX: in case of use of a same names written in different case (like "ENDER" and "Ender"),
    // m_presets_choice->GetSelection() will return first item, because search in PopupListCtrl is case-insensitive.
    // So, use GetSelection() from event parameter
    int selection = evt.GetSelection();

    auto idx = combo->get_extruder_idx();

    std::string preset_name =
        wxGetApp().preset_bundle->get_preset_name_by_alias(preset_type,
                                                           Preset::remove_suffix_modified(combo->GetString(selection).ToUTF8().data()));

    if (preset_type == Preset::TYPE_FFF_FILAMENT) {
        wxGetApp().preset_bundle->set_filament_preset(idx, preset_name);
    }

    bool select_preset = !combo->selection_is_changed_according_to_physical_printers();
    // TODO: ?
    DynamicPrintConfig phys = wxGetApp().preset_bundle->physical_printers.get_selected_printer().config;
    
    this->m_webView->LoadURL("http://" + phys.opt_string("print_host"));
    combo->update();
    
#ifdef __WXMSW__
    // From the Win 2004 preset combobox lose a focus after change the preset selection
    // and that is why the up/down arrow doesn't work properly
    // (see https://github.com/prusa3d/PrusaSlicer/issues/5531 ).
    // So, set the focus to the combobox explicitly
    combo->SetFocus();
#endif
}

void WebViewPanel::load_url(wxString &url)
{
    if (m_webView == nullptr)
        return;

    m_webView->LoadURL(url);

    if (url.IsEmpty()) {
        std::cout << "URL is empty";
    }
}


} // namespace GUI

} // namespace Slic3r
