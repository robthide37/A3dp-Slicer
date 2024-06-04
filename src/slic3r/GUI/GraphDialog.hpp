#ifndef _GRAPH_DIALOG_H_
#define _GRAPH_DIALOG_H_

#include <wx/spinctrl.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/checkbox.h>
#include <wx/msgdlg.h>

#include "RammingChart.hpp"

namespace Slic3r { namespace GUI {

struct GraphSettings
{
    std::string title;
    std::string description;
    std::string y_label;
    std::string x_label;
    std::string null_label;
    double min_x, max_x, step_x;
    double min_y, max_y, step_y;
    std::string label_min_x;
    std::string label_max_x;
    std::string label_min_y;
    std::string label_max_y;
    std::vector<GraphData::GraphType> allowed_types;
    GraphData reset_vals;
};

class GraphPanel : public wxPanel
{
public:
    GraphPanel(wxWindow *parent, GraphData data,const GraphSettings &settings);
    GraphData get_data();
    bool      is_disabled();

private:
    Chart *           m_chart        = nullptr;
    wxSpinCtrlDouble *m_widget_min_x = nullptr;
    wxSpinCtrlDouble *m_widget_max_x = nullptr;
    wxSpinCtrlDouble *m_widget_min_y = nullptr;
    wxSpinCtrlDouble *m_widget_max_y = nullptr;
    double            m_last_min_x   = 0.f;
    double            m_last_max_x   = 1.f;
    double            m_last_min_y   = 0.f;
    double            m_last_max_y   = 1.f;
};

class GraphDialog : public wxDialog
{
public:
    GraphDialog(wxWindow *parent, const GraphData &parameters, const GraphSettings &settings);
    GraphData get_data() { return m_output_data; }
    bool      is_disabled() { return m_disabled; }

private:
    GraphPanel *m_panel_graph = nullptr;
    GraphData m_output_data;
    bool m_disabled = false;
};

}}     // namespace Slic3r::GUI
#endif  // _GRAPH_DIALOG_H_