#ifndef slic3r_PrinterWebView_hpp_
#define slic3r_PrinterWebView_hpp_

#include <wx/wx.h>

#if __MSW__
#include <WebView2.h>
#endif

#include <wx/webview.h>
#include <wx/msw/webview_edge.h>
#include <wx/msw/webview_ie.h>
#include "slic3r/GUI/PresetComboBoxes.hpp"

namespace Slic3r { namespace GUI {

class WebViewPanel : public wxPanel
{
public:
    WebViewPanel(wxWindow *parent);
    ~WebViewPanel();

    wxWebView *m_webView;
    DevicePresetComboBox* m_combo_printer;
    void on_select_preset(wxCommandEvent &);

    void       load_url(wxString &url);

private:
    wxWebView *CreateWebView(wxWindow *parent);

};
}
}     // namespace Slic3r::GUI
#endif // slic3r_PrinterWebView_hpp
