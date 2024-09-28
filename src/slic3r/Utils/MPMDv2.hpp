///|/ Copyright (c) Superslicer 2021 - 2024 Durand Rémi @supermerill
///|/ Copyright (c) 2021 Alexander Bachler Jansson
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_MPMDv2_hpp_
#define slic3r_MPMDv2_hpp_

#include <optional>
#include <string>
#include <wx/string.h>
#include <wx/arrstr.h>

#include "PrintHost.hpp"
#include "libslic3r/PrintConfig.hpp"


namespace Slic3r {

class DynamicPrintConfig;
class Http;

class MPMDv2 : public PrintHost
{
public:
    MPMDv2(DynamicPrintConfig *config);
    ~MPMDv2() override = default;

    const char* get_name() const override;

    bool test(wxString &curl_msg) const override;
    wxString get_test_ok_msg () const override;
    wxString get_test_failed_msg (wxString &msg) const override;
    bool upload(PrintHostUpload upload_data, ProgressFn prorgess_fn, ErrorFn error_fn, InfoFn info_fn) const override;
    bool has_auto_discovery() const override { return false; }
    bool can_test() const override { return true; }
    PrintHostPostUploadActions get_post_upload_actions() const override { return PrintHostPostUploadAction::StartPrint; }
    std::string get_host() const override { return host; }

protected:
    virtual bool validate_version_text(const std::optional<std::string> &version_text) const;

private:
    std::string host;

    std::string make_url(const std::string &path) const;
};

}

#endif
