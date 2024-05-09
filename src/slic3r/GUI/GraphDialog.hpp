#ifndef _GRAPH_DIALOG_H_
#define _GRAPH_DIALOG_H_

#include <wx/spinctrl.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/checkbox.h>
#include <wx/msgdlg.h>

#include "RammingChart.hpp"

namespace Slic3r { namespace GUI {

class GraphPanel : public wxPanel
{
public:
    GraphPanel(wxWindow *parent, const std::string &data);
    std::string get_parameters();

private:
    Chart *           m_chart           = nullptr;
    wxSpinCtrl *      m_widget_speed    = nullptr;
    wxSpinCtrlDouble *m_widget_min_flow = nullptr;
    wxSpinCtrlDouble *m_widget_max_flow = nullptr;
    int               m_last_speed      = 120;
};

class GraphDialog : public wxDialog
{
public:
    GraphDialog(wxWindow *parent, const std::string &parameters);
    std::string get_parameters() { return m_output_data; }

private:
    GraphPanel *m_panel_graph = nullptr;
    std::string m_output_data;
};

}}     // namespace Slic3r::GUI
#endif  // _GRAPH_DIALOG_H_