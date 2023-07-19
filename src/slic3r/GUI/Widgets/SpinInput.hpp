#ifndef slic3r_GUI_SpinInput_hpp_
#define slic3r_GUI_SpinInput_hpp_

#include <wx/textctrl.h>
#include "StaticBox.hpp"

class Button;

class SpinInput : public wxNavigationEnabled<StaticBox>
{
    wxSize labelSize;
    StateColor   label_color;
    StateColor   text_color;
    wxTextCtrl * text_ctrl{nullptr};
    Button * button_inc {nullptr};
    Button * button_dec {nullptr};
    wxTimer timer;

    int val;
    int min;
    int max;
    int delta;

    static const int SpinInputWidth = 200;
    static const int SpinInputHeight = 50;

public:
    SpinInput();

    SpinInput(wxWindow *     parent,
              wxString       text,
              wxString       label = "",
              const wxPoint &pos   = wxDefaultPosition,
              const wxSize & size  = wxDefaultSize,
              long           style = 0,
              int min = 0, int max = 100, int initial = 0);

    void Create(wxWindow *     parent,
              wxString       text,
              wxString       label   = "",
              const wxPoint &pos     = wxDefaultPosition,
              const wxSize & size    = wxDefaultSize,
              long           style   = 0,
              int            min     = 0,
              int            max     = 100,
              int            initial = 0);

    void SetCornerRadius(double radius);

    void SetLabel(const wxString &label) wxOVERRIDE;

    void SetLabelColor(StateColor const &color);

    void SetTextColor(StateColor const &color);

    void SetSize(wxSize const &size);

    void Rescale();

    virtual bool Enable(bool enable = true) wxOVERRIDE;

    wxTextCtrl * GetText() { return text_ctrl; }

    void SetValue(const wxString &text);

    void SetValue (int value);

    int GetValue () const;
    wxString GetTextValue() const;

    void SetRange(int min, int max);

    bool SetFont(wxFont const& font) override;

    bool SetBackgroundColour(const wxColour& colour) override;
    bool SetForegroundColour(const wxColour& colour) override;
    void SetBorderColor(StateColor const& color);

    int GetMin() const { return this->min; }
    int GetMax() const { return this->max; }
    void SetSelection(long from, long to);

protected:
    void DoSetToolTipText(wxString const &tip) override;

private:
    void paintEvent(wxPaintEvent& evt);

    void render(wxDC& dc);

    void messureSize();

    Button *createButton(bool inc);

    // some useful events
    void mouseWheelMoved(wxMouseEvent& event);
    void keyPressed(wxKeyEvent& event);
    void onTimer(wxTimerEvent &evnet);
    void onTextLostFocus(wxEvent &event);
    void onText(wxCommandEvent &event);
    void onTextEnter(wxCommandEvent &event);

    void sendSpinEvent();

    DECLARE_EVENT_TABLE()
};

#endif // !slic3r_GUI_SpinInput_hpp_
