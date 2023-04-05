#include "Mainsail.hpp"

#include <algorithm>
#include <sstream>
#include <exception>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/nowide/convert.hpp>

#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/I18N.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/format.hpp"
#include "Http.hpp"

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;
namespace Slic3r {

Mainsail::Mainsail(DynamicPrintConfig *config) :
    m_host(config->opt_string("print_host")),
    m_apikey(config->opt_string("printhost_apikey")),
    m_cafile(config->opt_string("printhost_cafile")),
    m_ssl_revoke_best_effort(config->opt_bool("printhost_ssl_ignore_revoke"))
{}

const char* Mainsail::get_name() const { return "Mainsail"; }

wxString Mainsail::get_test_ok_msg () const
{
    return _(L("Connection to Mainsail works correctly."));
}

wxString Mainsail::get_test_failed_msg (wxString &msg) const
{
    return GUI::format_wxstr("%s: %s"
        , _L("Could not connect to Mainsail")
        , msg);
}

bool Mainsail::test(wxString& msg) const
{
    // GET /server/info

    // Since the request is performed synchronously here,
    // it is ok to refer to `msg` from within the closure
    const char* name = get_name();

    bool res = true;
    auto url = make_url("server/info");

    BOOST_LOG_TRIVIAL(info) << boost::format("%1%: Get version at: %2%") % name % url;

    auto http = Http::get(std::move(url));
    set_auth(http);
    http.on_error([&](std::string body, std::string error, unsigned status) {
        BOOST_LOG_TRIVIAL(error) << boost::format("%1%: Error getting version: %2%, HTTP %3%, body: `%4%`") % name % error % status % body;
        res = false;
        msg = format_error(body, error, status);
    })
    .on_complete([&, this](std::string body, unsigned) {
        BOOST_LOG_TRIVIAL(debug) << boost::format("%1%: Got server/info: %2%") % name % body;
            
        try {
            // All successful HTTP requests will return a json encoded object in the form of :
            // {result: <response data>}
            std::stringstream ss(body);
            pt::ptree ptree;
            pt::read_json(ss, ptree);
            if (ptree.front().first != "result") {
                msg = "Could not parse server response";
                res = false;
                return;
            }
            if (!ptree.front().second.get_optional<std::string>("moonraker_version")) {
                msg = "Could not parse server response";
                res = false;
                return;
            }
            BOOST_LOG_TRIVIAL(info) << boost::format("%1%: Got version: %2%") % name % ptree.front().second.get_optional<std::string>("moonraker_version");
        } catch (const std::exception&) {
            res = false;
            msg = "Could not parse server response";
        }
    })
    .perform_sync();

    return res;
}

bool Mainsail::upload(PrintHostUpload upload_data, ProgressFn prorgess_fn, ErrorFn error_fn, InfoFn info_fn) const
{
    // POST /server/files/upload

    const char* name = get_name();
    const auto upload_filename = upload_data.upload_path.filename();
    const auto upload_parent_path = upload_data.upload_path.parent_path();

    // If test fails, test_msg_or_host_ip contains the error message.
    wxString test_msg_or_host_ip;
    if (!test(test_msg_or_host_ip)) {
        error_fn(std::move(test_msg_or_host_ip));
        return false;
    }

    std::string url;
    bool res = true;

    url = make_url("server/files/upload");
    
    BOOST_LOG_TRIVIAL(info) << boost::format("%1%: Uploading file %2% at %3%, filename: %4%, path: %5%, print: %6%")
        % name
        % upload_data.source_path
        % url
        % upload_filename.string()
        % upload_parent_path.string()
        % (upload_data.post_action == PrintHostPostUploadAction::StartPrint ? "true" : "false");
    /*
    The file must be uploaded in the request's body multipart/form-data (ie: <input type="file">). The following arguments may also be added to the form-data:
    root: The root location in which to upload the file.Currently this may be gcodes or config.If not specified the default is gcodes.
    path : This argument may contain a path(relative to the root) indicating a subdirectory to which the file is written.If a path is present the server will attempt to create any subdirectories that do not exist.
    checksum : A SHA256 hex digest calculated by the client for the uploaded file.If this argument is supplied the server will compare it to its own checksum calculation after the upload has completed.A checksum mismatch will result in a 422 error.
    Arguments available only for the gcodes root :
    print: If set to "true", Klippy will attempt to start the print after uploading.Note that this value should be a string type, not boolean.This provides compatibility with OctoPrint's upload API.
    */
    auto http = Http::post(std::move(url));
    set_auth(http);
    
    http.form_add("root", "gcodes");
    if (!upload_parent_path.empty())
        http.form_add("path", upload_parent_path.string());
    if (upload_data.post_action == PrintHostPostUploadAction::StartPrint)
        http.form_add("print", "true");
   
    http.form_add_file("file", upload_data.source_path.string(), upload_filename.string())
        .on_complete([&](std::string body, unsigned status) {
            BOOST_LOG_TRIVIAL(debug) << boost::format("%1%: File uploaded: HTTP %2%: %3%") % name % status % body;
        })
        .on_error([&](std::string body, std::string error, unsigned status) {
            BOOST_LOG_TRIVIAL(error) << boost::format("%1%: Error uploading file: %2%, HTTP %3%, body: `%4%`") % name % error % status % body;
            error_fn(format_error(body, error, status));
            res = false;
        })
        .on_progress([&](Http::Progress progress, bool& cancel) {
            prorgess_fn(std::move(progress), cancel);
            if (cancel) {
                // Upload was canceled
                BOOST_LOG_TRIVIAL(info) << name << ": Upload canceled";
                res = false;
            }
        })
#ifdef WIN32
        .ssl_revoke_best_effort(m_ssl_revoke_best_effort)
#endif
        .perform_sync();

    return res;
}

void Mainsail::set_auth(Http &http) const
{
    if (!m_apikey.empty())
        http.header("X-Api-Key", m_apikey);
    if (!m_cafile.empty())
        http.ca_file(m_cafile);
}

std::string Mainsail::make_url(const std::string &path) const
{
    if (m_host.find("http://") == 0 || m_host.find("https://") == 0) {
        if (m_host.back() == '/') {
            return (boost::format("%1%%2%") % m_host % path).str();
        } else {
            return (boost::format("%1%/%2%") % m_host % path).str();
        }
    } else {
        return (boost::format("http://%1%/%2%") % m_host % path).str();
    }
}


}
