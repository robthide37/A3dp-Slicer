#include "WifiConfigDialog.hpp"

#include "GUI_App.hpp"
#include "GUI.hpp"
#include "I18N.hpp"
#include "format.hpp"
#include "RemovableDriveManager.hpp"
#include "MsgDialog.hpp"

#include <wx/stattext.h>
#include <wx/button.h>
#include <boost/nowide/convert.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>

namespace Slic3r {
namespace GUI {

WifiConfigDialog::WifiConfigDialog(wxWindow* parent, std::string& file_path, RemovableDriveManager* removable_manager)
     : DPIDialog(parent, wxID_ANY, _L("Physical Printer Instalation"), wxDefaultPosition, wxDefaultSize/*wxSize(25 * wxGetApp().em_unit(), 20 * wxGetApp().em_unit())*/, wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER)
     , m_wifi_scanner(new WifiScanner())
     , out_file_path(file_path)
     , m_removable_manager(removable_manager)
{
    const int em = GUI::wxGetApp().em_unit();

    wxPanel* panel = new wxPanel(this);
    wxBoxSizer* vsizer = new wxBoxSizer(wxVERTICAL);
    panel->SetSizer(vsizer);

    auto* ssid_sizer = new wxBoxSizer(wxHORIZONTAL);
    // TRN SSID of WiFi network.
    wxStaticText* ssid_label = new wxStaticText(panel, wxID_ANY, GUI::format_wxstr("%1%:", _L("SSID")));
    m_ssid_combo = new wxComboBox(panel, wxID_ANY);
#if __APPLE__
    m_ssid_combo->SetToolTip(_L("On some versions of MacOS, this only loads SSID of connencted network."));
#endif // __APPLE__
    rescan_networks(false);
    // TRN Text of button to rescan visible networks in Wifi Config dialog.
    wxButton* ssid_button = new wxButton(panel, wxID_ANY, _(L("Rescan")));
    ssid_sizer->Add(m_ssid_combo, 1, wxALIGN_CENTER_VERTICAL, 10);
    ssid_sizer->Add(ssid_button, 0);

    auto* pass_sizer = new wxBoxSizer(wxHORIZONTAL);
    // TRN Password of WiFi network.
    wxStaticText* password_label = new wxStaticText(panel, wxID_ANY, GUI::format_wxstr("%1%:", _L("Password")));
    m_pass_textctrl = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxDefaultSize);
    pass_sizer->Add(m_pass_textctrl, 1, wxALIGN_CENTER_VERTICAL, 10);
#if __APPLE__
    // TRN Text of button to retrieve password from keychain in Wifi Config dialog. Only on Mac.
    wxButton* pass_button = new wxButton(panel, wxID_ANY, _(L("Retrieve")));
    pass_sizer->Add(pass_button, 0);
    pass_button->Bind(wxEVT_BUTTON, &WifiConfigDialog::on_retrieve_password, this);
#endif // __APPLE__
    // show password if current ssid was selected already
    fill_password();

    auto* drive_sizer = new wxBoxSizer(wxHORIZONTAL);
    // TRN description of Combo Box with path to USB drive.
    wxStaticText* drive_label = new wxStaticText(panel, wxID_ANY, GUI::format_wxstr("%1%:", _L("Drive")));
    m_drive_combo = new wxComboBox(panel, wxID_ANY);
    rescan_drives(preffered_drive);
    // TRN Text of button to rescan connect usb drives in Wifi Config dialog.
    wxButton* drive_button = new wxButton(panel, wxID_ANY, _(L("Rescan")));
    drive_sizer->Add(m_drive_combo, 1, wxALIGN_CENTER_VERTICAL, 10);
    drive_sizer->Add(drive_button, 0);

    // TRN Text of button to write config file in Wifi Config dialog.
    wxButton* ok_button = new wxButton(panel, wxID_OK, _L("Write"));
    wxButton* cancel_button = new wxButton(panel, wxID_CANCEL);

    auto* grid = new wxFlexGridSizer(2, 15, 15);
    grid->AddGrowableCol(1);

    grid->Add(ssid_label, 0, wxALIGN_CENTER_VERTICAL);
    grid->Add(ssid_sizer, 0, wxEXPAND);

    grid->Add(password_label, 0, wxALIGN_CENTER_VERTICAL);
    grid->Add(pass_sizer, 0, wxEXPAND);

    grid->Add(drive_label, 0, wxALIGN_CENTER_VERTICAL);
    grid->Add(drive_sizer, 0, wxEXPAND);

    vsizer->Add(grid, 0, wxEXPAND | wxTOP | wxBOTTOM, 15);

    wxBoxSizer* buttons_sizer = new wxBoxSizer(wxHORIZONTAL);
    buttons_sizer->Add(ok_button, 1, wxLEFT);
    buttons_sizer->AddStretchSpacer();
    buttons_sizer->Add(cancel_button, 1, wxRIGHT);

    vsizer->Add(buttons_sizer, 0, wxEXPAND);

    wxBoxSizer* topsizer = new wxBoxSizer(wxVERTICAL);
    topsizer->Add(panel, 1, wxEXPAND | wxALL, 15);
    SetSizerAndFit(topsizer);

    ok_button->Bind(wxEVT_BUTTON, &WifiConfigDialog::on_ok, this);
    m_ssid_combo->Bind(wxEVT_TEXT, &WifiConfigDialog::on_combo, this);
    drive_button->Bind(wxEVT_BUTTON, &WifiConfigDialog::on_rescan_drives, this);
    ssid_button->Bind(wxEVT_BUTTON, &WifiConfigDialog::on_rescan_networks, this);
}

WifiConfigDialog::~WifiConfigDialog()
{
    if (m_wifi_scanner)
        delete m_wifi_scanner;
}

void WifiConfigDialog::on_combo(wxCommandEvent& e)
{
    fill_password();
}

void WifiConfigDialog::fill_password()
{
    assert(m_ssid_combo && m_pass_textctrl && m_wifi_scanner);
    if (auto it = m_wifi_scanner->get_map().find(m_ssid_combo->GetValue()); it != m_wifi_scanner->get_map().end())
        m_pass_textctrl->SetValue(boost::nowide::widen(it->second));
}

void WifiConfigDialog::on_retrieve_password(wxCommandEvent& e)
{
    if (m_ssid_combo->GetValue().empty()) {
        return;
    }
    
    std::string psk = m_wifi_scanner->get_psk(boost::nowide::narrow(m_ssid_combo->GetValue()));
    if (psk.empty()) {
        // TRN Alert message when retrieving password for wifi from keychain. Probably will display only on Apple so keychain is MacOS term.
        wxString msg = _L("No password in the keychain for given SSID.");
        show_info(nullptr, msg);
        return;
    }
    m_pass_textctrl->SetValue(boost::nowide::widen(psk));
}

void WifiConfigDialog::on_rescan_drives(wxCommandEvent& e)
{
   rescan_drives();
}

void WifiConfigDialog::rescan_drives()
{
    assert(m_drive_combo && m_removable_manager);
    m_drive_combo->Clear();
    std::vector<DriveData> ddata = m_removable_manager->get_drive_list();
    for (const auto& data : ddata) {
        m_drive_combo->Append(boost::nowide::widen(data.path));
    }
    if (m_drive_combo->GetCount() > 0) {
        m_drive_combo->Select(0);
    }
}

void WifiConfigDialog::on_rescan_networks(wxCommandEvent& e)
{
    rescan_networks(true);
}

void WifiConfigDialog::rescan_networks(bool select)
{
    assert(m_ssid_combo && m_wifi_scanner);
    m_wifi_scanner->scan();
    std::string current = m_wifi_scanner->get_current_ssid();
    const auto& map = m_wifi_scanner->get_map();
    m_ssid_combo->Clear();
    for (const auto pair : map) {
        m_ssid_combo->Append(pair.first);
        // select ssid of current network (if connected)
        if (current == pair.first)
            m_ssid_combo->Select(m_ssid_combo->GetCount() - 1);
    }
    if (m_ssid_combo->GetSelection() == -1 && m_ssid_combo->GetCount() > 0)
        m_ssid_combo->Select(0);
    if (select && m_ssid_combo->GetSelection() != -1)
        fill_password();
}

void WifiConfigDialog::on_ok(wxCommandEvent& e)
{
    assert(m_ssid_combo && m_pass_textctrl);
    if (m_ssid_combo->GetValue().empty()) {
        // TRN Alert message when writing WiFi configuration file to usb drive.
        wxString msg = _L("SSID field is empty.");
        show_info(nullptr, msg);
        return;
    }
   
    std::string selected_path = boost::nowide::narrow(m_drive_combo->GetValue());

    if (selected_path.empty()) {
        // TRN Alert message when writing WiFi configuration file to usb drive.
        wxString msg = _L("Drive field is empty.");
        show_info(nullptr, msg);
        return;
    }

    boost::filesystem::path file_path = boost::filesystem::path(selected_path) / "prusa_printer_settings.ini";

    bool path_on_removable_media = m_removable_manager->set_and_verify_last_save_path(file_path.string());
    if (!path_on_removable_media) {
        // TRN Alert message when writing WiFi configuration file to usb drive.
        wxString msg = _L("Selected path is not on removable media.");
        show_info(nullptr, msg);
        return;
    }

    if (boost::filesystem::exists(file_path))
    {
        wxString msg_text = GUI::format_wxstr("%1% already exists. Do you want to rewrite it?", file_path.string());
        WarningDialog dialog(m_parent, msg_text, _L("Warning"), wxYES | wxNO);
        if (dialog.ShowModal() == wxID_NO)
        {
            return;
        }
    }
   
    wxString ssid = m_ssid_combo->GetValue();
    wxString pass = m_pass_textctrl->GetValue();
    wxString data = GUI::format(
        "[wifi]\n"
        "ssid=%1%\n"
        "psk=%2%\n"
        , ssid
        , pass
    );

    FILE* file;
    file = fopen(file_path.string().c_str(), "w");
    if (file == NULL) {
        BOOST_LOG_TRIVIAL(error) << "Failed to write to file " << file_path;
        // TODO show error
        show_error(nullptr,  _L("Failed to open file for writing."));
        return;
    }
    fwrite(data.c_str(), 1, data.size(), file);
    fclose(file);
    out_file_path = file_path.string();
    this->EndModal(wxID_OK);
}

void WifiConfigDialog::on_dpi_changed(const wxRect& suggested_rect)
{
    SetFont(wxGetApp().normal_font());

    const int em = em_unit();

    //msw_buttons_rescale(this, em, { wxID_OK, wxID_CANCEL });

    Fit();
    Refresh();
}
}}// Slicer::GUI