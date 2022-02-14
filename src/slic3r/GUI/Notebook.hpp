#ifndef slic3r_Notebook_hpp_
#define slic3r_Notebook_hpp_

#ifdef _WIN32

#include <wx/bookctrl.h>

class ModeSizer;
class ScalableButton;

// custom message the ButtonsListCtrl sends to its parent (Notebook) to notify a selection change:
wxDECLARE_EVENT(wxCUSTOMEVT_NOTEBOOK_SEL_CHANGED, wxCommandEvent);
wxDECLARE_EVENT(wxCUSTOMEVT_NOTEBOOK_BT_PRESSED, wxCommandEvent);

class ButtonsListCtrl : public wxControl
{
public:
    ButtonsListCtrl(wxWindow* parent, bool add_mode_buttons = false);
    ~ButtonsListCtrl() {}

    void OnPaint(wxPaintEvent&);
    void SetSelection(int sel);
    int  GetSelection() { return m_selection; }
    void UpdateMode();
    void Rescale();
    bool InsertPage(size_t n, const wxString& text, bool bSelect = false, const std::string& bmp_name = "", const int bmp_size = 16);
    void RemovePage(size_t n);
    bool InsertSpacer(size_t n, int size);
    void RemoveSpacer(size_t n);
    bool SetPageImage(size_t n, const std::string& bmp_name, const int bmp_size = -1) const;
    bool SetPageImage(size_t n, const wxBitmap& bmp) const;
    void SetPageText(size_t n, const wxString& strText);
    wxString GetPageText(size_t n) const;
    ScalableButton* GetPageButton(size_t n);

private:
    wxWindow*                       m_parent;
    wxFlexGridSizer*                m_buttons_sizer;
    wxBoxSizer*                     m_sizer;
    std::vector<ScalableButton*>    m_pageButtons;
    std::vector<bool>               m_spacers;
    int                             m_selection {-1};
    int                             m_btn_margin;
    int                             m_line_margin;
    ModeSizer*                      m_mode_sizer {nullptr};
};

// A tabpane but with custom buttons.
// The top buttons are stored in ButtonsListCtrl.
// The panes are stored here, in the wxBookCtrlBase.
// It's possible to add "fake" button that link to an existing pane.
// This notebook selection and tab count include the fake buttons.
// So it's possible to return the same pane for two different index.
// It's possible to set icons for the buttons, but they have no index as it's managed in the ButtonsListCtrl
class Notebook: public wxBookCtrlBase
{
public:
    Notebook(wxWindow * parent,
                 wxWindowID winid = wxID_ANY,
                 const wxPoint & pos = wxDefaultPosition,
                 const wxSize & size = wxDefaultSize,
                 long style = 0,
                 bool add_mode_buttons = false)
    {
        Init();
        Create(parent, winid, pos, size, style, add_mode_buttons);
    }

    bool Create(wxWindow * parent,
                wxWindowID winid = wxID_ANY,
                const wxPoint & pos = wxDefaultPosition,
                const wxSize & size = wxDefaultSize,
                long style = 0,
                bool add_mode_buttons = false)
    {
        if (!wxBookCtrlBase::Create(parent, winid, pos, size, style | wxBK_TOP))
            return false;

        m_bookctrl = new ButtonsListCtrl(this, add_mode_buttons);

        wxSizer* mainSizer = new wxBoxSizer(IsVertical() ? wxVERTICAL : wxHORIZONTAL);

        if (style & wxBK_RIGHT || style & wxBK_BOTTOM)
            mainSizer->Add(0, 0, 1, wxEXPAND, 0);

        m_controlSizer = new wxBoxSizer(IsVertical() ? wxHORIZONTAL : wxVERTICAL);
        m_controlSizer->Add(m_bookctrl, wxSizerFlags(1).Expand());
        wxSizerFlags flags;
        if (IsVertical())
            flags.Expand();
        else
            flags.CentreVertical();
        mainSizer->Add(m_controlSizer, flags.Border(wxALL, m_controlMargin));
        SetSizer(mainSizer);

        this->Bind(wxCUSTOMEVT_NOTEBOOK_SEL_CHANGED, [this](wxCommandEvent& evt)
        {                    
            if (int page_idx = evt.GetId(); page_idx >= 0)
                SetSelection(page_idx);
        });

        this->Bind(wxEVT_NAVIGATION_KEY, &Notebook::OnNavigationKey, this);

        return true;
    }


    // Methods specific to this class.

    // A method allowing to add a new page without any label (which is unused
    // by this control) and show it immediately.
    bool ShowNewPage(wxWindow * page)
    {
        return AddPage(page, wxString(), ""/*true *//* select it */);
    }


    // Set effect to use for showing/hiding pages.
    void SetEffects(wxShowEffect showEffect, wxShowEffect hideEffect)
    {
        m_showEffect = showEffect;
        m_hideEffect = hideEffect;
    }

    // Or the same effect for both of them.
    void SetEffect(wxShowEffect effect)
    {
        SetEffects(effect, effect);
    }

    // And the same for time outs.
    void SetEffectsTimeouts(unsigned showTimeout, unsigned hideTimeout)
    {
        m_showTimeout = showTimeout;
        m_hideTimeout = hideTimeout;
    }

    void SetEffectTimeout(unsigned timeout)
    {
        SetEffectsTimeouts(timeout, timeout);
    }


    // Implement base class pure virtual methods.

    // adds a new page to the control
    bool AddPage(wxWindow* page,
                 const wxString& text,
                 const std::string& bmp_name,
                 bool bSelect = false)
    {
        DoInvalidateBestSize();
        return InsertPage(GetPageCount(), page, text, bmp_name, bSelect);
    }

    // Page management
    virtual bool InsertPage(size_t n,
                            wxWindow * page,
                            const wxString & text,
                            bool bSelect = false,
                            int imageId = NO_IMAGE) override
    {
        //update btidx_to_tabpage
        while (btidx_to_tabpage.size() < n)
            btidx_to_tabpage.push_back(-1);
        int16_t last = -1;
        for (int i = 0; i < n; i++)
            if (btidx_to_tabpage[i] >= 0)
                last = i;
        for (int i = 0; i < btidx_to_tabpage.size(); i++)
            if (btidx_to_tabpage[i] > last)
                btidx_to_tabpage[i]++;
        btidx_to_tabpage.insert(btidx_to_tabpage.begin() + n, last + 1);

        if (!wxBookCtrlBase::InsertPage(last + 1, page, text, bSelect, imageId))
            return false;

        GetBtnsListCtrl()->InsertPage(n, text, bSelect);

        if (!DoSetSelectionAfterInsertion(n, bSelect))
            page->Hide();
        else
            EmitEventSelChanged(last + 1);

        return true;
    }

    void EmitEventSelChanged(int16_t new_sel);

    bool InsertPage(size_t n,
        wxWindow* page,
        const wxString& text,
        const std::string& bmp_name = "",
        const int bmp_size = 16,
        bool bSelect = false)
    {

        //update btidx_to_tabpage
        while (btidx_to_tabpage.size() < n)
            btidx_to_tabpage.push_back(-1);
        int16_t last = -1;
        for (int i = 0; i < n; i++)
            if (btidx_to_tabpage[i] >= 0)
                last = i;
        for (int i = 0; i < btidx_to_tabpage.size(); i++)
            if (btidx_to_tabpage[i] > last)
                btidx_to_tabpage[i]++;
        btidx_to_tabpage.insert(btidx_to_tabpage.begin() + n, last + 1);

        //create the tab
        if (!wxBookCtrlBase::InsertPage(last + 1, page, text, bSelect))
            return false;

        //create the button
        GetBtnsListCtrl()->InsertPage(n, text, bSelect, bmp_name, bmp_size);

        if (bSelect)
            SetSelection(n);

        return true;
    }

    void InsertFakePage(size_t n,
        const int page_idx,
        const wxString& text,
        const std::string& bmp_name = "",
        const int bmp_size = 16,
        bool bSelect = false)
    {
        //update btidx_to_tabpage
        while (btidx_to_tabpage.size() < n)
            btidx_to_tabpage.push_back(-1);
        btidx_to_tabpage.insert(btidx_to_tabpage.begin() + n, page_idx);

        //create the button
        GetBtnsListCtrl()->InsertPage(n, text, bSelect, bmp_name, bmp_size);


        if (bSelect)
            SetSelection(n);

    }

    virtual int SetSelection(size_t n) override
    {
        GetBtnsListCtrl()->SetSelection(n);
        int16_t real_page = n< btidx_to_tabpage.size() ? btidx_to_tabpage[n] : -1;
        if (real_page >= 0) {
            int ret = DoSetSelection(real_page, SetSelection_SendEvent);

            // check that only the selected page is visible and others are hidden:
            for (size_t page = 0; page < m_pages.size(); page++)
                if (page != real_page)
                    m_pages[page]->Hide();

            EmitEventSelChanged(n);

            return ret;
        }
        return -1;
    }

    virtual int ChangeSelection(size_t n) override
    {
        GetBtnsListCtrl()->SetSelection(n);
        int16_t real_page = n < btidx_to_tabpage.size() ? btidx_to_tabpage[n] : -1;
        if (real_page >= 0) {
            int ret = DoSetSelection(real_page);
            EmitEventSelChanged(n);
            return ret;
        }
        return -1;
    }

    virtual int GetButtonSelection() const { return GetBtnsListCtrl()->GetSelection(); }

    // Neither labels nor images are supported but we still store the labels
    // just in case the user code attaches some importance to them.
    virtual bool SetPageText(size_t n, const wxString & strText) override
    {
        wxCHECK_MSG(n < GetPageCount(), false, wxS("Invalid page"));

        GetBtnsListCtrl()->SetPageText(n, strText);

        return true;
    }

    virtual wxString GetPageText(size_t n) const override
    {
        wxCHECK_MSG(n < GetPageCount(), wxString(), wxS("Invalid page"));
        return GetBtnsListCtrl()->GetPageText(n);
    }

    virtual bool SetPageImage(size_t n, int imageId) override
    {
        return GetBtnsListCtrl()->SetPageImage(n, this->GetImageList()->GetBitmap(imageId));
    }

    virtual int GetPageImage(size_t WXUNUSED(n)) const override
    {
        return NO_IMAGE;
    }

    bool SetPageImage(size_t n, const std::string& bmp_name, const int bmp_size)
    {
        return GetBtnsListCtrl()->SetPageImage(n, bmp_name, bmp_size);
    }

    // Override some wxWindow methods too.
    virtual void SetFocus() override
    {
        wxWindow* const page = GetCurrentPage();
        if (page)
            page->SetFocus();
    }

    ButtonsListCtrl* GetBtnsListCtrl() const { return static_cast<ButtonsListCtrl*>(m_bookctrl); }

    void UpdateMode()
    {
        GetBtnsListCtrl()->UpdateMode();
    }

    void Rescale()
    {
        GetBtnsListCtrl()->Rescale();
    }

    void Notebook::OnNavigationKey(wxNavigationKeyEvent& event)
    {
        if (event.IsWindowChange()) {
            // change pages
            AdvanceSelection(event.GetDirection());
        }
        else {
            // we get this event in 3 cases
            //
            // a) one of our pages might have generated it because the user TABbed
            // out from it in which case we should propagate the event upwards and
            // our parent will take care of setting the focus to prev/next sibling
            //
            // or
            //
            // b) the parent panel wants to give the focus to us so that we
            // forward it to our selected page. We can't deal with this in
            // OnSetFocus() because we don't know which direction the focus came
            // from in this case and so can't choose between setting the focus to
            // first or last panel child
            //
            // or
            //
            // c) we ourselves (see MSWTranslateMessage) generated the event
            //
            wxWindow* const parent = GetParent();

            // the wxObject* casts are required to avoid MinGW GCC 2.95.3 ICE
            const bool isFromParent = event.GetEventObject() == (wxObject*)parent;
            const bool isFromSelf = event.GetEventObject() == (wxObject*)this;
            const bool isForward = event.GetDirection();

            if (isFromSelf && !isForward)
            {
                // focus is currently on notebook tab and should leave
                // it backwards (Shift-TAB)
                event.SetCurrentFocus(this);
                parent->HandleWindowEvent(event);
            }
            else if (isFromParent || isFromSelf)
            {
                // no, it doesn't come from child, case (b) or (c): forward to a
                // page but only if entering notebook page (i.e. direction is
                // backwards (Shift-TAB) comething from out-of-notebook, or
                // direction is forward (TAB) from ourselves),
                if (m_selection != wxNOT_FOUND &&
                    (!event.GetDirection() || isFromSelf))
                {
                    // so that the page knows that the event comes from it's parent
                    // and is being propagated downwards
                    event.SetEventObject(this);

                    wxWindow* page = m_pages[m_selection];
                    if (!page->HandleWindowEvent(event))
                    {
                        page->SetFocus();
                    }
                    //else: page manages focus inside it itself
                }
                else // otherwise set the focus to the notebook itself
                {
                    SetFocus();
                }
            }
            else
            {
                // it comes from our child, case (a), pass to the parent, but only
                // if the direction is forwards. Otherwise set the focus to the
                // notebook itself. The notebook is always the 'first' control of a
                // page.
                if (!isForward)
                {
                    SetFocus();
                }
                else if (parent)
                {
                    event.SetCurrentFocus(this);
                    parent->HandleWindowEvent(event);
                }
            }
        }
    }

protected:
    virtual void UpdateSelectedPage(size_t newsel) override
    {
        // even if Nothing to do here, but must be overridden to avoid the assert in
        // the base class version.
    }

    virtual wxBookCtrlEvent * CreatePageChangingEvent() const override
    {
        return new wxBookCtrlEvent(wxEVT_BOOKCTRL_PAGE_CHANGING,
                                   GetId());
    }

    virtual void MakeChangedEvent(wxBookCtrlEvent & event) override
    {
        event.SetEventType(wxEVT_BOOKCTRL_PAGE_CHANGED);
    }

    virtual wxWindow * DoRemovePage(size_t page) override
    {
        wxWindow* const win = wxBookCtrlBase::DoRemovePage(page);
        if (win)
        {
            int16_t real_page = page < btidx_to_tabpage.size() ? btidx_to_tabpage[page] : -1;
            GetBtnsListCtrl()->RemovePage(page);
            btidx_to_tabpage.erase(btidx_to_tabpage.begin() + page);
            if (real_page >= 0) {
                DoSetSelectionAfterRemoval(real_page);
                for (int i = page; i < btidx_to_tabpage.size(); i++)
                    if (btidx_to_tabpage[i] > real_page)
                        btidx_to_tabpage[i]--;
                    else if(btidx_to_tabpage[i] == real_page)
                        btidx_to_tabpage[i] = -1;
            }
        }

        return win;
    }

    virtual void DoSize() override
    {
        wxWindow* const page = GetCurrentPage();
        if (page)
            page->SetSize(GetPageRect());
    }

    virtual void DoShowPage(wxWindow * page, bool show) override
    {
        if (show)
            page->ShowWithEffect(m_showEffect, m_showTimeout);
        else
            page->HideWithEffect(m_hideEffect, m_hideTimeout);
    }

private:
    void Init()
    {
        // We don't need any border as we don't have anything to separate the
        // page contents from.
        SetInternalBorder(0);

        // No effects by default.
        m_showEffect =
        m_hideEffect = wxSHOW_EFFECT_NONE;

        m_showTimeout =
        m_hideTimeout = 0;
    }

    wxShowEffect m_showEffect,
                 m_hideEffect;

    unsigned m_showTimeout,
             m_hideTimeout;

    std::vector<int16_t> btidx_to_tabpage;

    ButtonsListCtrl* m_ctrl{ nullptr };

};
#endif // _WIN32
#endif // slic3r_Notebook_hpp_
