#include <algorithm>
#include <sstream>
#include "GraphDialog.hpp"
#include "BitmapCache.hpp"
#include "GUI.hpp"
#include "I18N.hpp"
#include "GUI_App.hpp"
#include "MsgDialog.hpp"

#include <wx/sizer.h>
#include <wx/spinctrl.h>

namespace Slic3r { namespace GUI {

int scale(const int val) { return val * Slic3r::GUI::wxGetApp().em_unit(); }
#ifdef __WXGTK3__
int ITEM_WIDTH() { return scale(10); }
#else
int ITEM_WIDTH() { return scale(6); }
#endif

static void update_ui(wxWindow *window) { Slic3r::GUI::wxGetApp().UpdateDarkUI(window); }

GraphDialog::GraphDialog(wxWindow *parent, const std::string &parameters)
    : wxDialog(parent,
               wxID_ANY,
               _(L("Extrusion multiplier per extrusion speed")),
               wxDefaultPosition,
               wxDefaultSize,
               wxDEFAULT_DIALOG_STYLE /* | wxRESIZE_BORDER*/)
{
    update_ui(this);
    m_panel_graph = new GraphPanel(this, parameters);

    // Not found another way of getting the background colours of GraphDialog, GraphPanel and Chart correct than
    // setting them all explicitely. Reading the parent colour yielded colour that didn't really match it, no
    // wxSYS_COLOUR_... matched colour used for the dialog. Same issue (and "solution") here :
    // https://forums.wxwidgets.org/viewtopic.php?f=1&t=39608 Whoever can fix this, feel free to do so.
#ifndef _WIN32
    this->SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_FRAMEBK));
    m_panel_graph->SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_FRAMEBK));
#endif
    m_panel_graph->Show(true);
    this->Show();

    auto main_sizer = new wxBoxSizer(wxVERTICAL);
    main_sizer->Add(m_panel_graph, 1, wxEXPAND | wxTOP | wxLEFT | wxRIGHT, 5);
    main_sizer->Add(CreateButtonSizer(wxOK | wxCANCEL), 0, wxALIGN_CENTER_HORIZONTAL | wxTOP | wxBOTTOM, 10);
    SetSizer(main_sizer);
    main_sizer->SetSizeHints(this);

    update_ui(static_cast<wxButton *>(this->FindWindowById(wxID_OK, this)));
    update_ui(static_cast<wxButton *>(this->FindWindowById(wxID_CANCEL, this)));

    this->Bind(wxEVT_CLOSE_WINDOW, [this](wxCloseEvent &e) { 
        EndModal(wxCANCEL);
    });

    this->Bind(
        wxEVT_BUTTON,
        [this](wxCommandEvent &) {
            m_output_data = m_panel_graph->get_parameters();
            EndModal(wxID_OK);
        },
        wxID_OK);
    this->Show();
    // Slic3r::GUI::MessageDialog dlg(this, _(L("Graph .")), _(L("Warning")), wxOK | wxICON_EXCLAMATION);
    // dlg.ShowModal();
}

#ifdef _WIN32
#define style wxSP_ARROW_KEYS | wxBORDER_SIMPLE
#else
#define style wxSP_ARROW_KEYS
#endif



GraphPanel::GraphPanel(wxWindow *parent, const std::string &parameters)
    : wxPanel(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize /*,wxPoint(50,50), wxSize(800,350),wxBORDER_RAISED*/)
{
    update_ui(this);
    auto sizer_chart = new wxBoxSizer(wxVERTICAL);

    std::stringstream stream{parameters};
    //stream >> m_graph_line_width_multiplicator >> m_graph_step_multiplicator;
    int   Graph_speed_size = 0;
    float dummy            = 0.f;
    float min              = 2.f;
    float max              = 0.f;
    while (stream >> dummy) {
        ++Graph_speed_size;
        if (dummy > 0 && dummy <= 2) {
            min = std::min(min, dummy);
            max = std::max(max, dummy);
        }
    }
    stream.clear();
    stream.get();

    if (min >= max) {
        max = 1.2f;
        min = 0.9f;
    } else {
        min = int(min * 10 - 1 + EPSILON) / 10.f;
        max = int(1.9f + max * 10 - EPSILON) / 10.f;
    }

    std::vector<std::pair<float, float>> buttons;
    float                                x = 0.f;
    float                                y = 0.f;
    while (stream >> x >> y) buttons.push_back(std::make_pair(x, y));

    m_chart = new Chart(this, wxRect(scale(1), scale(1), scale(64), scale(36)), buttons, scale(1));
    m_chart->set_manual_points_manipulation(true);
    m_chart->set_xy_range(0, min, Graph_speed_size * 10.f, max);
    m_chart->set_x_label(_L("Print speed") + " ("+_L("mm/s")+")", 1.f);
    m_chart->set_y_label(_L("Extrusion multiplier"), 0.001f);
    m_chart->set_no_point_label(_L("No compensation"));
#ifdef _WIN32
    update_ui(m_chart);
#else
    m_chart->SetBackgroundColour(parent->GetBackgroundColour()); // see comment in GraphDialog constructor
#endif
    sizer_chart->Add(new wxStaticText(this, wxID_ANY, 
        _L("Choose the extrusion multipler value for multiple speeds.\nYou can add/remove points with a right clic.")));
    sizer_chart->Add(m_chart, 0, wxALL, 5);
    
    m_last_speed   = Graph_speed_size * 10;
    m_widget_speed = new wxSpinCtrl(this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize(ITEM_WIDTH(), -1),
                                    style | wxTE_PROCESS_ENTER, 10, 2000, m_last_speed);
    // note: wxTE_PROCESS_ENTER allow the wxSpinCtrl to receive wxEVT_TEXT_ENTER events

    
    m_widget_min_flow = new wxSpinCtrlDouble(this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize(ITEM_WIDTH(), -1),
                                    style, 0.1, 2, min, 0.1f);
    m_widget_max_flow = new wxSpinCtrlDouble(this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize(ITEM_WIDTH(), -1),
                                    style, 0.1, 2, max, 0.1f);

#ifdef _WIN32
    update_ui(m_widget_speed);
    update_ui(m_widget_min_flow);
    update_ui(m_widget_max_flow);
#endif
    // line for speed max & reset
    wxBoxSizer *size_line = new wxBoxSizer(wxHORIZONTAL);
    size_line->Add(new wxStaticText(this, wxID_ANY, wxString(_L("Graph max speed") + " :")), 0, wxALIGN_CENTER_VERTICAL);
    size_line->Add(m_widget_speed);
    size_line->AddSpacer(80);
    wxButton *bt_reset = new wxButton(this, wxID_ANY, _L("Reset"));
    bt_reset->SetToolTip(_L("Reset all values to 1. Also reset all points to defaults."));
    size_line->Add(bt_reset);
    sizer_chart->Add(size_line);

    //line for y min & max
    size_line = new wxBoxSizer(wxHORIZONTAL);
    size_line->Add(new wxStaticText(this, wxID_ANY, wxString(_L("Minimum flow") + " :")), 0, wxALIGN_CENTER_VERTICAL);
    size_line->Add(m_widget_min_flow);
    size_line->AddSpacer(20);
    size_line->Add(new wxStaticText(this, wxID_ANY, wxString(_L("Maximum flow") + " :")), 0, wxALIGN_CENTER_VERTICAL);
    size_line->Add(m_widget_max_flow);
    sizer_chart->Add(size_line);
    
    sizer_chart->SetSizeHints(this);
    SetSizer(sizer_chart);

    bt_reset->Bind(wxEVT_BUTTON, ([this](wxCommandEvent& e) {
        std::vector<std::pair<float,float>> buttons;// = m_chart->get_buttons();
        //for (std::pair<float, float> &button : buttons) {
        //    button.second = 1.f;
        //}
        //buttons.emplace_back(5,1.f);
        buttons.emplace_back(10,1.f);
        //buttons.emplace_back(15,1.f);
        buttons.emplace_back(20,1.f);
        buttons.emplace_back(30,1.f);
        buttons.emplace_back(40,1.f);
        buttons.emplace_back(60,1.f);
        buttons.emplace_back(80,1.f);
        buttons.emplace_back(120,1.f);
        buttons.emplace_back(160,1.f);
        buttons.emplace_back(240,1.f);
        buttons.emplace_back(480,1.f);
        buttons.emplace_back(640,1.f);
        buttons.emplace_back(960,1.f);
        buttons.emplace_back(1280,1.f);
        m_chart->set_buttons(buttons);
    }));

    m_widget_speed->Bind(wxEVT_TEXT_ENTER, [this](wxCommandEvent &) {
        int old_speed = m_last_speed;
        m_last_speed = 10 * ((m_widget_speed->GetValue() + 5) / 10);
        m_last_speed = std::min(std::max(m_last_speed, 20), 2000);
        m_widget_speed->SetValue(m_last_speed);
        if (old_speed < m_last_speed) {
            if (old_speed < 1000 && m_last_speed >= 1000)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 100.f);
            else if (old_speed < 100 && m_last_speed >= 100)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 10.f);
        } else {
            if (old_speed >= 100 && m_last_speed < 100)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 1.f);
            else if (old_speed >= 1000 && m_last_speed < 1000)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 10.f);
        }
        m_chart->set_xy_range(0, -1, m_last_speed, -1);
    });
    m_widget_speed->Bind(wxEVT_SPINCTRL, [this](wxSpinEvent &evt) {
        assert(evt.GetInt() != m_last_speed);
        int incr = 10;
        if (m_last_speed >= 60)
            incr = 20;
        if (m_last_speed >= 100)
            incr = 50;
        if (m_last_speed >= 300)
            incr = 100;
        if (m_last_speed >= 600)
            incr = 200;
        if (m_last_speed >= 1000)
            incr = 500;
        if (evt.GetInt() > m_last_speed) {
            if (m_last_speed < 100 && m_last_speed + incr >= 100)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 10.f);
            else if (m_last_speed < 1000 && m_last_speed + incr >= 1000)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 100.f);
            m_last_speed += incr;
        } else {
            if (m_last_speed >= 100 && m_last_speed - incr < 100)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 1.f);
            else if (m_last_speed >= 1000 && m_last_speed - incr < 1000)
                m_chart->set_x_label(_L("Print speed") + " (" + _L("mm/s") + ")", 10.f);
            m_last_speed -= incr;
        }
        m_last_speed = std::min(std::max(m_last_speed, 20), 2000);
        m_widget_speed->SetValue(m_last_speed);
        m_chart->set_xy_range(0, 0, m_last_speed, -1);
    });
    // thses don't work, i don't know why.
    //m_widget_speed->Bind(wxEVT_SPIN_UP, [this](wxSpinEvent &evt) {
    //        std::cout<<"up";
    //    });
    //m_widget_speed->Bind(wxEVT_SPIN_DOWN, [this](wxSpinEvent &evt) {
    //        std::cout<<"down";
    //    });
    
    m_widget_min_flow->Bind(wxEVT_TEXT, [this](wxCommandEvent &) {
        m_chart->set_xy_range(-1, m_widget_min_flow->GetValue(), -1, m_widget_max_flow->GetValue());
    });
    m_widget_min_flow->Bind(wxEVT_CHAR, [](wxKeyEvent &) {});   // do nothing - prevents the user to change the value
    m_widget_max_flow->Bind(wxEVT_TEXT, [this](wxCommandEvent &) {
        m_chart->set_xy_range(-1, m_widget_min_flow->GetValue(), -1, m_widget_max_flow->GetValue());
    });
    m_widget_max_flow->Bind(wxEVT_CHAR, [](wxKeyEvent &) {});   // do nothing - prevents the user to change the value
    Bind(EVT_WIPE_TOWER_CHART_CHANGED, [this](wxCommandEvent &) {
        int nb_samples = m_chart->get_speed(10.f).size();
        m_last_speed = 10 * nb_samples;
        m_widget_speed->SetValue(m_last_speed);
    });
    Refresh(true); // erase background
}

std::string GraphPanel::get_parameters()
{
    std::vector<float>                   flow_rates  = m_chart->get_speed(10.f);
    std::vector<std::pair<float, float>> buttons = m_chart->get_buttons();

    //write string
    std::stringstream                    stream;
    //stream << m_graph_line_width_multiplicator << " " << m_graph_step_multiplicator;
    //if all are at 1, then set the first to 0 so it's "disabled"
    bool disabled = true;
    for (const float &flow_rate : flow_rates)
        if (flow_rate != 1.)
            disabled = false;
    for (size_t i = 0; i < flow_rates.size() ; i++) {
        const float &flow_rate = flow_rates[i];
        if (0 == i) {
            stream << (disabled ? 0.f : flow_rate);
        } else {
            stream << " " << flow_rate;
        }
    }
    stream << "|";
    for (const auto &button : buttons) stream << " " << button.first << " " << button.second;
    return stream.str();
}

}} // namespace Slic3r::GUI
