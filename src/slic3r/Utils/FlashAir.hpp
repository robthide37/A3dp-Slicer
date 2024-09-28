///|/ Copyright (c) Prusa Research 2018 - 2022 David Kocík @kocikdav, Vojtěch Bubník @bubnikv, Vojtěch Král @vojtechkral
///|/ Copyright (c) 2019 Stephan Reichhelm @stephanr
///|/ Copyright (c) 2018 Martin Loidl @LoidlM
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_FlashAir_hpp_
#define slic3r_FlashAir_hpp_

#include <string>
#include <wx/string.h>

#include "PrintHost.hpp"


namespace Slic3r {

class DynamicPrintConfig;
class Http;

class FlashAir : public PrintHost
{
public:
	FlashAir(DynamicPrintConfig *config);
	~FlashAir() override = default;

	const char* get_name() const override;

	bool test(wxString &curl_msg) const override;
	wxString get_test_failed_msg(wxString &msg) const override;
	bool upload(PrintHostUpload upload_data, ProgressFn prorgess_fn, ErrorFn error_fn, InfoFn info_fn) const override;
	bool has_auto_discovery() const override { return false; }
	bool can_test() const override { return true; }
    PrintHostPostUploadActions get_post_upload_actions() const override { return {}; }
	std::string get_host() const override { return host; }
    
private:
	std::string host;

	std::string timestamp_str() const;
	std::string make_url(const std::string &path) const;
	std::string make_url(const std::string &path, const std::string &arg, const std::string &val) const;
};

}

#endif
