#include "A3DPWebView.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/Utils/MacDarkMode.hpp"
#include <boost/log/trivial.hpp>
#include <slic3r/GUI/Widgets/StateColor.hpp>
#include <wx/webviewarchivehandler.h>
#include <wx/webviewfshandler.h>
#if wxUSE_WEBVIEW_EDGE
#include <wx/msw/webview_edge.h>
#elif defined(__WXMAC__)
#include <wx/osx/webview_webkit.h>
#endif
#include <wx/uri.h>
#if defined(__WIN32__) || defined(__WXMAC__)
#include "wx/private/jsscriptwrapper.h"
#endif
#ifdef __WIN32__
#include <WebView2.h>
#include <Shellapi.h>
#include <slic3r/Utils/Http.hpp>
#elif defined __linux__
#include <gtk/gtk.h>
#define WEBKIT_API struct WebKitWebView; struct WebKitJavascriptResult; extern "C" { WEBKIT_API void webkit_web_view_run_javascript ( WebKitWebView *web_view, const gchar *script, GCancellable *cancellable, GAsyncReadyCallback callback, gpointer user_data ) ; WEBKIT_API WebKitJavascriptResult * webkit_web_view_run_javascript_finish ( WebKitWebView *web_view, GAsyncResult *result, GError **error ) ; WEBKIT_API void webkit_javascript_result_unref ( WebKitJavascriptResult *js_result ) ; }
#endif

namespace fs = boost::filesystem;

#ifdef __WIN32__
// Run Download and Install in another thread so we don't block the UI thread
DWORD DownloadAndInstallWV2RT() {
    int returnCode = 2; // Download failed
    fs::path target_file_path = (fs::temp_directory_path() / "MicrosoftEdgeWebview2Setup.exe");
    bool downloaded = false;
    Slic3r::Http::get("https://go.microsoft.com/fwlink/p/?LinkId=2124703")
        .on_error([](std::string body, std::string error, unsigned http_status) {})
        .on_complete([&downloaded, target_file_path](std::string body, unsigned http_status) {
            fs::fstream file(target_file_path, std::ios::out | std::ios::binary | std::ios::trunc);
            file.write(body.c_str(), body.size());
            file.flush();
            file.close();
            downloaded = true;
        })
        .perform_sync();
    // Sleep for 1 second to wait for the buffer written into disk
    std::this_thread::sleep_for(std::chrono::seconds(1));
    if (downloaded) {
        // Either Package the WebView2 Bootstrapper with your app or download it using fwlink
        // Then invoke install at Runtime.
        SHELLEXECUTEINFOW shExInfo = {0};
        shExInfo.cbSize = sizeof(shExInfo);
        shExInfo.fMask = SEE_MASK_NOCLOSEPROCESS;
        shExInfo.hwnd = NULL;
        shExInfo.lpVerb = L"runas";
        shExInfo.lpFile = target_file_path.wstring().c_str();
        shExInfo.lpParameters = L"";
        shExInfo.lpDirectory = NULL;
        shExInfo.nShow = SW_SHOWNORMAL;
        shExInfo.hInstApp = NULL;
        if (ShellExecuteExW(&shExInfo)) {
            WaitForSingleObject(shExInfo.hProcess, INFINITE);
            DWORD exitCode;
            GetExitCodeProcess(shExInfo.hProcess, &exitCode);
            CloseHandle(shExInfo.hProcess);
            returnCode = exitCode;
        }
    }
    return returnCode;
}
#endif

wxWebView* WebView::CreateWebView(wxWindow *parent, wxString const &url) {
    #if wxUSE_WEBVIEW_EDGE
    wxFileName edgeFixedDir(wxStandardPaths::Get().GetExecutablePath());
    edgeFixedDir.SetFullName("");
    edgeFixedDir.AppendDir("edge_fixed");
    if (edgeFixedDir.DirExists()) {
        wxWebViewEdge::MSWSetBrowserExecutableDir(edgeFixedDir.GetFullPath());
        wxLogMessage("Using fixed edge version");
    }
    #endif
    auto url2 = url;
    #ifdef __WIN32__
    url2.Replace("\\\\", "/");
    #endif
    if (!url2.empty()) {
        url2 = wxURI(url2).BuildURI();
    }
    BOOST_LOG_TRIVIAL(trace) << __FUNCTION__ << ": " << url2.ToUTF8();
    #ifdef __WIN32__
    wxWebView* webView = new WebViewEdge;
    #elif defined(__WXOSX__)
    wxWebView *webView = new WebViewWebKit;
    #else
    auto webView = wxWebView::New();
    #endif
    if (webView) {
        webView->SetBackgroundColour(wxColour(*wxWHITE));
        #ifdef __WIN32__
        webView->SetUserAgent(wxString::Format("BBL-Slicer/v%s (%s) Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36 Edg/107.0.1418.52", SLIC3R_VERSION, Slic3r::GUI::wxGetApp().dark_mode() ? "dark" : "light"));
        webView->Create(parent, wxID_ANY, url2, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE);
        webView->RegisterHandler(wxSharedPtr<wxWebViewHandler>(new wxWebViewArchiveHandler("bbl")));
        webView->RegisterHandler(wxSharedPtr<wxWebViewHandler>(new wxWebViewFSHandler("memory")));
        #else
        webView->RegisterHandler(wxSharedPtr<wxWebViewHandler>(new wxWebViewArchiveHandler("wxfs")));
        webView->RegisterHandler(wxSharedPtr<wxWebViewHandler>(new wxWebViewFSHandler("memory")));
        webView->Create(parent, wxID_ANY, url2, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE);
        webView->SetUserAgent(wxString::Format("BBL-Slicer/v%s (%s) Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko)", SLIC3R_VERSION, Slic3r::GUI::wxGetApp().dark_mode() ? "dark" : "light"));
        #endif
        webView->EnableContextMenu(true);
    } else {
        BOOST_LOG_TRIVIAL(info) << __FUNCTION__ << ": failed. Use fake web view.";
        webView = new FakeWebView;
    }
    webView->SetRefData(new WebViewRef(webView));
    g_webviews.push_back(webView);
    return webView;
}

void WebView::LoadUrl(wxWebView *webView, wxString const &url) {
    auto url2 = url;
    #ifdef __WIN32__
    url2.Replace("\\\\", "/");
    #endif
    if (!url2.empty()) {
        url2 = wxURI(url2).BuildURI();
    }
    BOOST_LOG_TRIVIAL(trace) << __FUNCTION__ << ": " << url2.ToUTF8();
    webView->LoadURL(url2);
}

#if wxUSE_WEBVIEW_EDGE
bool WebView::CheckWebViewRuntime() {
    wxWebViewFactoryEdge factory;
    auto wxVersion = factory.GetVersionInfo();
    return wxVersion.GetMajor() != 0;
}

bool WebView::DownloadAndInstallWebViewRuntime() {
    return DownloadAndInstallWV2RT() == 0;
}
#endif
