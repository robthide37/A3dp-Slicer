#include "AppUpdater.hpp"

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/ini_parser.hpp> 
#include <curl/curl.h>

#include "slic3r/GUI/format.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/Utils/Http.hpp"

#ifdef _WIN32
#include <shellapi.h>
#include <Shlobj_core.h>
#include <windows.h>
#include <KnownFolders.h>
#include <shlobj.h>
#endif // _WIN32


namespace Slic3r {

namespace {
	
#ifdef _WIN32
	bool run_file(const boost::filesystem::path& path)
	{
		// find updater exe
		if (boost::filesystem::exists(path)) {
			// run updater. Original args: /silent -restartapp prusa-slicer.exe -startappfirst

			// Using quoted string as mentioned in CreateProcessW docs, silent execution parameter.
			std::wstring wcmd = L"\"" + path.wstring();

			// additional information
			STARTUPINFOW si;
			PROCESS_INFORMATION pi;

			// set the size of the structures
			ZeroMemory(&si, sizeof(si));
			si.cb = sizeof(si);
			ZeroMemory(&pi, sizeof(pi));

			// start the program up
			if (CreateProcessW(NULL,   // the path
				wcmd.data(),    // Command line
				NULL,           // Process handle not inheritable
				NULL,           // Thread handle not inheritable
				FALSE,          // Set handle inheritance to FALSE
				0,              // No creation flags
				NULL,           // Use parent's environment block
				NULL,           // Use parent's starting directory 
				&si,            // Pointer to STARTUPINFO structure
				&pi             // Pointer to PROCESS_INFORMATION structure (removed extra parentheses)
			)) {
				// Close process and thread handles.
				CloseHandle(pi.hProcess);
				CloseHandle(pi.hThread);
				return true;
			}
			else {
				BOOST_LOG_TRIVIAL(error) << "Failed to run " << wcmd;
			}
		}
		return false;
	}

	std::string get_downloads_path()
	{
		std::string ret;
		PWSTR path = NULL;
		HRESULT hr = SHGetKnownFolderPath(FOLDERID_Downloads, 0, NULL, &path);
		if (SUCCEEDED(hr)) {
			ret = boost::nowide::narrow(path);
		}
		CoTaskMemFree(path);
		return ret;
	}

	bool open_folder(const boost::filesystem::path& path)
	{
		// this command can run the installer exe as well, but is it better than CreateProcessW?
		ShellExecuteW(NULL, NULL, path.wstring().c_str(), NULL, NULL, SW_SHOWNORMAL);
		return true;
	}

#elif __linux__
	bool run_file(const boost::filesystem::path& path)
	{	
		return false;
	}

	std::string get_downloads_path()
	{
		wxString command = "xdg-user-dir DOWNLOAD";
		wxArrayString output;
		 
		//Check if we're running in an AppImage container, if so, we need to remove AppImage's env vars,
		// because they may mess up the environment expected by the file manager.
		// Mostly this is about LD_LIBRARY_PATH, but we remove a few more too for good measure.
		if (wxGetEnv("APPIMAGE", nullptr)) {
			// We're running from AppImage
			wxEnvVariableHashMap env_vars;
			wxGetEnvMap(&env_vars);

			env_vars.erase("APPIMAGE");
			env_vars.erase("APPDIR");
			env_vars.erase("LD_LIBRARY_PATH");
			env_vars.erase("LD_PRELOAD");
			env_vars.erase("UNION_PRELOAD");

			wxExecuteEnv exec_env;
			exec_env.env = std::move(env_vars);

			wxString owd;
			if (wxGetEnv("OWD", &owd)) {
				// This is the original work directory from which the AppImage image was run,
				// set it as CWD for the child process:
					exec_env.cwd = std::move(owd);
			}

			::wxExecute(command, output, 0, &exec_env);
				
		} else {
			// Looks like we're NOT running from AppImage, we'll make no changes to the environment.
			::wxExecute(command, output);
		}
		if (output.GetCount() > 0) {
			return boost::nowide::narrow(output[0]);
		}
		return std::string();
	}

	bool open_folder(const boost::filesystem::path& path)
	{
		if (boost::filesystem::is_directory(path)) {
			const char *argv[] = { "xdg-open", path.string().c_str(), nullptr };

			// Check if we're running in an AppImage container, if so, we need to remove AppImage's env vars,
			// because they may mess up the environment expected by the file manager.
			// Mostly this is about LD_LIBRARY_PATH, but we remove a few more too for good measure.
			if (wxGetEnv("APPIMAGE", nullptr)) {
				// We're running from AppImage
				wxEnvVariableHashMap env_vars;
				wxGetEnvMap(&env_vars);

				env_vars.erase("APPIMAGE");
				env_vars.erase("APPDIR");
				env_vars.erase("LD_LIBRARY_PATH");
				env_vars.erase("LD_PRELOAD");
				env_vars.erase("UNION_PRELOAD");

				wxExecuteEnv exec_env;
				exec_env.env = std::move(env_vars);

				wxString owd;
				if (wxGetEnv("OWD", &owd)) {
					// This is the original work directory from which the AppImage image was run,
					// set it as CWD for the child process:
					exec_env.cwd = std::move(owd);
				}

				::wxExecute(const_cast<char**>(argv), wxEXEC_ASYNC, nullptr, &exec_env);

			} else {
				// Looks like we're NOT running from AppImage, we'll make no changes to the environment.
				::wxExecute(const_cast<char**>(argv), wxEXEC_ASYNC, nullptr, nullptr);
			}
			return true;
		}
		return false;
	}

#elif  __APPLE__
	bool run_file(const boost::filesystem::path& path)
	{
		if (boost::filesystem::exists(path)) {
			// attach downloaded dmg file
            const char* argv1[] = { "hdiutil", "attach", path.string().c_str(), nullptr };
            ::wxExecute(const_cast<char**>(argv1), wxEXEC_ASYNC, nullptr);
            // open inside attached as a folder in finder
            const char* argv2[] = { "open", "/Volumes/PrusaSlicer", nullptr };
			::wxExecute(const_cast<char**>(argv2), wxEXEC_ASYNC, nullptr);
			return true;
		}
		return false;
	}

	std::string get_downloads_path()
	{
		// call objective-c implementation
		return get_downloads_path_mac();
	}

	bool open_folder(const boost::filesystem::path& path)
	{

		if (boost::filesystem::is_directory(path)) {
			const char* argv[] = { "open", path.string().c_str(), nullptr };
			::wxExecute(const_cast<char**>(argv), wxEXEC_ASYNC, nullptr);
			return true;
		}
		return false;
	}
#endif 
} // namespace

wxDEFINE_EVENT(EVT_SLIC3R_VERSION_ONLINE, wxCommandEvent);
wxDEFINE_EVENT(EVT_SLIC3R_EXPERIMENTAL_VERSION_ONLINE, wxCommandEvent);
wxDEFINE_EVENT(EVT_SLIC3R_APP_DOWNLOAD_PROGRESS, wxCommandEvent);
wxDEFINE_EVENT(EVT_SLIC3R_APP_DOWNLOAD_FAILED, wxCommandEvent);

// priv handles all operations in separate thread
// 1) download version file and parse it.
// 2) download new app file and open in folder / run it.
struct AppUpdater::priv {
	priv();
	// Download file. What happens with the data is specified in completefn.
	bool http_get_file(const std::string& url
		, size_t size_limit
		, std::function<bool(Http::Progress)> progress_fn
		, std::function<bool(std::string /*body*/, std::string& error_message)> completefn
		, std::string& error_message
	) const;

	// Download installer / app
	boost::filesystem::path download_file(const DownloadAppData& data) const;
	// Run file in m_last_dest_path
	bool run_downloaded_file(boost::filesystem::path path);
	// gets version file via http
	void version_check(const std::string& version_check_url);
#if 0
	// parsing of Prusaslicer.version2 
	void parse_version_string_old(const std::string& body) const;
#endif
	// parses ini tree of version file, saves to m_online_version_data and queue event(s) to UI 
	void parse_version_string(const std::string& body);
	// thread
	std::thread				m_thread;
	std::atomic_bool        m_cancel;
	std::mutex				m_data_mutex;
	// read only variable used to init m_online_version_data.target_path
	boost::filesystem::path m_default_dest_folder; // readonly
	// DownloadAppData read / write needs to be locked by m_data_mutex
	DownloadAppData			m_online_version_data;
	DownloadAppData get_app_data();
	void            set_app_data(DownloadAppData data);	
};

AppUpdater::priv::priv() :
	m_cancel (false)
#ifdef __linux__
    , m_default_dest_folder (boost::filesystem::path("/tmp"))
#else
	, m_default_dest_folder (boost::filesystem::path(data_dir()) / "cache")
#endif //_WIN32
{	
	boost::filesystem::path downloads_path = boost::filesystem::path(get_downloads_path());
	if (!downloads_path.empty()) {
		m_default_dest_folder = std::move(downloads_path);
	}
	BOOST_LOG_TRIVIAL(error) << "Default download path: " << m_default_dest_folder;
	
}

bool  AppUpdater::priv::http_get_file(const std::string& url, size_t size_limit, std::function<bool(Http::Progress)> progress_fn, std::function<bool(std::string /*body*/, std::string& error_message)> complete_fn, std::string& error_message) const
{
	bool res = false;
	Http::get(url)
		.size_limit(size_limit)
		.on_progress([&, progress_fn](Http::Progress progress, bool& cancel) {
			// progress function returns true as success (to continue) 
			cancel = (this->m_cancel ? true : !progress_fn(std::move(progress)));
			if (cancel) {
				error_message = GUI::format("Error getting: `%1%`: Download was canceled.",
					url);
				BOOST_LOG_TRIVIAL(debug) << "AppUpdater::priv::http_get_file message: "<< error_message;
			}
		})
		.on_error([&](std::string body, std::string error, unsigned http_status) {
			error_message = GUI::format("Error getting: `%1%`: HTTP %2%, %3%",
				url,
				http_status,
				error);
			BOOST_LOG_TRIVIAL(error) << error_message;
		})
		.on_complete([&](std::string body, unsigned /* http_status */) {
			assert(complete_fn != nullptr);
			res = complete_fn(body, error_message);
		})
		.perform_sync();
	
	return res;
}

boost::filesystem::path AppUpdater::priv::download_file(const DownloadAppData& data) const
{
	boost::filesystem::path dest_path;
	size_t last_gui_progress = 0;
	size_t expected_size = data.size;
	dest_path = data.target_path;
	assert(!dest_path.empty());
	if (dest_path.empty())
	{
		BOOST_LOG_TRIVIAL(error) << "Download from " << data.url << " could not start. Destination path is empty.";
		return boost::filesystem::path();
	}
	std::string error_message;
	bool res = http_get_file(data.url, 80 * 1024 * 1024 //TODO: what value here
		// on_progress
		, [&last_gui_progress, expected_size](Http::Progress progress) {
			// size check
			if (progress.dltotal > 0 && progress.dltotal > expected_size) {
				std::string message = GUI::format("Downloading new %1% has failed. The file has incorrect file size. Aborting download.\nExpected size: %2%\nDownload size: %3%", SLIC3R_APP_NAME, expected_size, progress.dltotal);
				BOOST_LOG_TRIVIAL(error) << message;
				wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_APP_DOWNLOAD_FAILED);
				evt->SetString(message);
				GUI::wxGetApp().QueueEvent(evt);
				return false;
			} else if (progress.dltotal > 0 && progress.dltotal < expected_size) {
				BOOST_LOG_TRIVIAL(error) << GUI::format("Downloading new %1% has incorrect size. The download will continue. \nExpected size: %2%\nDownload size: %3%", SLIC3R_APP_NAME, expected_size, progress.dltotal);;
			}
			// progress event
			size_t gui_progress = progress.dltotal > 0 ? 100 * progress.dlnow / progress.dltotal : 0;
			//BOOST_LOG_TRIVIAL(error) << "App download " << gui_progress << "% " << progress.dlnow << " of " << progress.dltotal;
			if (last_gui_progress < gui_progress && (last_gui_progress != 0 || gui_progress != 100)) {
				last_gui_progress = gui_progress;
				wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_APP_DOWNLOAD_PROGRESS);
				evt->SetString(GUI::from_u8(std::to_string(gui_progress)));
				GUI::wxGetApp().QueueEvent(evt);
			}
			return true;
		}
		// on_complete
		, [dest_path, expected_size](std::string body, std::string& error_message){
			// Size check. Does always 1 char == 1 byte?
			size_t body_size = body.size(); 
			if (body_size != expected_size) {
				BOOST_LOG_TRIVIAL(error) << "Downloaded file has wrong size. Expected size: " <<  expected_size << " Downloaded size: " << body_size;
				return false;
			}
			boost::filesystem::path tmp_path = dest_path;
			tmp_path += format(".%1%%2%", get_current_pid(), ".download");
			try
			{
				boost::filesystem::fstream file(tmp_path, std::ios::out | std::ios::binary | std::ios::trunc);
				file.write(body.c_str(), body.size());
				file.close();
				boost::filesystem::rename(tmp_path, dest_path);
			}
			catch (const std::exception&)
			{
				BOOST_LOG_TRIVIAL(error) << "Failed to write and move " << tmp_path << " to " << dest_path;
				return false;
			}
			return true;
		}
		, error_message
	);
	if (!res)
	{
		if (this->m_cancel)
		{
			BOOST_LOG_TRIVIAL(error) << error_message;
		} else {
			std::string message = GUI::format("Downloading new %1% has failed:\n%2%", SLIC3R_APP_NAME, error_message);
			BOOST_LOG_TRIVIAL(error) << message;
			wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_APP_DOWNLOAD_FAILED);
			evt->SetString(message);
			GUI::wxGetApp().QueueEvent(evt);
		}
		return boost::filesystem::path();
	}
	
	return dest_path;
}

bool AppUpdater::priv::run_downloaded_file(boost::filesystem::path path)
{
	assert(!path.empty());
	bool res = run_file(path);
	BOOST_LOG_TRIVIAL(error) << "Running "<< path.string() << " was " << res;
	return res;
}

void AppUpdater::priv::version_check(const std::string& version_check_url) 
{
	assert(!version_check_url.empty());
	std::string error_message;
	bool res = http_get_file(version_check_url, 1024
		// on_progress
		, [](Http::Progress progress) { return true; }
		// on_complete
		, [&](std::string body, std::string& error_message) {
			boost::trim(body);
			parse_version_string(body);
			return true;
		}
		, error_message
	);
	if (!res)
		BOOST_LOG_TRIVIAL(error) << "Failed to download version file: " << error_message;
}

void AppUpdater::priv::parse_version_string(const std::string& body)
{
	size_t start = body.find('[');
	if (start == std::string::npos) {
#if 0
		BOOST_LOG_TRIVIAL(error) << "Could not find property tree in version file. Starting old parsing.";
		parse_version_string_old(body);
		return;
#endif // 0
		BOOST_LOG_TRIVIAL(error) << "Could not find property tree in version file. Checking for application update has failed.";
		return;
	}
	std::string tree_string = body.substr(start);
	boost::property_tree::ptree tree;
	std::stringstream ss(tree_string);
	try {
		boost::property_tree::read_ini(ss, tree);
	} catch (const boost::property_tree::ini_parser::ini_parser_error& err) {
		//throw Slic3r::RuntimeError(format("Failed reading version file property tree Error: \"%1%\" at line %2%. \nTree:\n%3%", err.message(), err.line(), tree_string).c_str());
		BOOST_LOG_TRIVIAL(error) << format("Failed reading version file property tree Error: \"%1%\" at line %2%. \nTree:\n%3%", err.message(), err.line(), tree_string);
		return;
	}

	DownloadAppData new_data;

	for (const auto& section : tree) {
		std::string section_name = section.first;

		// online release version info
		if (section_name == 
#ifdef _WIN32
			"release:win64"
#elif __linux__
			"release:linux"
#else
			"release:osx"
#endif
			) {
			for (const auto& data : section.second) {
				if (data.first == "url") {
					new_data.url = data.second.data();
					new_data.target_path = m_default_dest_folder / AppUpdater::get_filename_from_url(new_data.url);
					BOOST_LOG_TRIVIAL(error) << format("parsing version string: url: %1%", new_data.url);
				} else if (data.first == "size"){
					new_data.size = std::stoi(data.second.data());
					BOOST_LOG_TRIVIAL(error) << format("parsing version string: expected size: %1%", new_data.size);
				}
			}
		}

		// released versions - to be send to UI layer
		if (section_name == "common") {
			std::vector<std::string> prerelease_versions;
			for (const auto& data : section.second) {
				// release version - save and send to UI layer
				if (data.first == "release") {
					std::string version = data.second.data();
					boost::optional<Semver> release_version = Semver::parse(version);
					if (!release_version) {
						BOOST_LOG_TRIVIAL(error) << format("Received invalid contents from version file: Not a correct semver: `%1%`", version);
						return;
					}
					new_data.version = release_version;
					// Send after all data is read
					/*
					BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, version);
					wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_VERSION_ONLINE);
					evt->SetString(GUI::from_u8(version));
					GUI::wxGetApp().QueueEvent(evt);
					*/
				// prerelease versions - write down to be sorted and send to UI layer
				} else if (data.first == "alpha") {
					prerelease_versions.emplace_back(data.second.data());
				} else if (data.first == "beta") {
					prerelease_versions.emplace_back(data.second.data());
				} else if (data.first == "rc") {
					prerelease_versions.emplace_back(data.second.data());
				}
			}
			// find recent version that is newer than last full release.
			boost::optional<Semver> recent_version;
			std::string				version_string;
			for (const std::string& ver_string : prerelease_versions) {
				boost::optional<Semver> ver = Semver::parse(ver_string);
				if (ver && *new_data.version < *ver && ((recent_version && *recent_version < *ver) || !recent_version)) {
					recent_version = ver;
					version_string = ver_string;
				}
			}
			// send prerelease version to UI layer
			if (recent_version) {
				BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, version_string);
				wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_EXPERIMENTAL_VERSION_ONLINE);
				evt->SetString(GUI::from_u8(version_string));
				GUI::wxGetApp().QueueEvent(evt);
			}
		}
	}
	assert(!new_data.url.empty());
	assert(new_data.version);
	// save
	set_app_data(new_data);
	// send
	std::string version = new_data.version.get().to_string();
	BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, version);
	wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_VERSION_ONLINE);
	evt->SetString(GUI::from_u8(version));
	GUI::wxGetApp().QueueEvent(evt);
}

#if 0
void AppUpdater::priv::parse_version_string_old(const std::string& body) const
{

	// release version
	std::string version;
	const auto first_nl_pos = body.find_first_of("\n\r");
	if (first_nl_pos != std::string::npos)
		version = body.substr(0, first_nl_pos);
	else
		version = body;
	boost::optional<Semver> release_version = Semver::parse(version);
	if (!release_version) {
		BOOST_LOG_TRIVIAL(error) << format("Received invalid contents from `%1%`: Not a correct semver: `%2%`", SLIC3R_APP_NAME, version);
		return;
	}
	BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, version);
	wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_VERSION_ONLINE);
	evt->SetString(GUI::from_u8(version));
	GUI::wxGetApp().QueueEvent(evt);

	// alpha / beta version
	std::vector<std::string> prerelease_versions;
	size_t nexn_nl_pos = first_nl_pos;
	while (nexn_nl_pos != std::string::npos && body.size() > nexn_nl_pos + 1) {
		const auto last_nl_pos = nexn_nl_pos;
		nexn_nl_pos = body.find_first_of("\n\r", last_nl_pos + 1);
		std::string line;
		if (nexn_nl_pos == std::string::npos)
			line = body.substr(last_nl_pos + 1);
		else
			line = body.substr(last_nl_pos + 1, nexn_nl_pos - last_nl_pos - 1);

		// alpha
		if (line.substr(0, 6) == "alpha=") {
			version = line.substr(6);
			if (!Semver::parse(version)) {
				BOOST_LOG_TRIVIAL(error) << format("Received invalid contents for alpha release from `%1%`: Not a correct semver: `%2%`", SLIC3R_APP_NAME, version);
				return;
			}
			prerelease_versions.emplace_back(version);
			// beta
		}
		else if (line.substr(0, 5) == "beta=") {
			version = line.substr(5);
			if (!Semver::parse(version)) {
				BOOST_LOG_TRIVIAL(error) << format("Received invalid contents for beta release from `%1%`: Not a correct semver: `%2%`", SLIC3R_APP_NAME, version);
				return;
			}
			prerelease_versions.emplace_back(version);
		}
	}
	// find recent version that is newer than last full release.
	boost::optional<Semver> recent_version;
	for (const std::string& ver_string : prerelease_versions) {
		boost::optional<Semver> ver = Semver::parse(ver_string);
		if (ver && *release_version < *ver && ((recent_version && *recent_version < *ver) || !recent_version)) {
			recent_version = ver;
			version = ver_string;
		}
	}
	if (recent_version) {
		BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, version);
		wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_EXPERIMENTAL_VERSION_ONLINE);
		evt->SetString(GUI::from_u8(version));
		GUI::wxGetApp().QueueEvent(evt);
	}
}
#endif // 0

DownloadAppData AppUpdater::priv::get_app_data()
{
	const std::lock_guard<std::mutex> lock(m_data_mutex);
	DownloadAppData ret_val(m_online_version_data);
	return ret_val;
}

void AppUpdater::priv::set_app_data(DownloadAppData data)
{
	const std::lock_guard<std::mutex> lock(m_data_mutex);
	m_online_version_data = data;
}

AppUpdater::AppUpdater()
	:p(new priv())
{
}
AppUpdater::~AppUpdater()
{
	if (p && p->m_thread.joinable()) {
		// This will stop transfers being done by the thread, if any.
		// Cancelling takes some time, but should complete soon enough.
		p->m_cancel = true;
		p->m_thread.join();
	}
}
void AppUpdater::sync_download()
{
	assert(p);
	// join thread first - it could have been in sync_version
	if (p->m_thread.joinable()) {
		// This will stop transfers being done by the thread, if any.
		// Cancelling takes some time, but should complete soon enough.
		p->m_cancel = true;
		p->m_thread.join();
	}
	p->m_cancel = false;

	DownloadAppData input_data = p->get_app_data();
	assert(!input_data.url.empty());

 	p->m_thread = std::thread(
		[this, input_data]() {
			if (boost::filesystem::path dest_path = p->download_file(input_data); boost::filesystem::exists(dest_path)){
				if (input_data.start_after) {
					p->run_downloaded_file(std::move(dest_path));
				} else {
					open_folder(dest_path.parent_path());
				}
			}
		});
}

void AppUpdater::sync_version(const std::string& version_check_url)
{
	assert(p);
	// join thread first - it could have been in sync_download
	if (p->m_thread.joinable()) {
		// This will stop transfers being done by the thread, if any.
		// Cancelling takes some time, but should complete soon enough.
		p->m_cancel = true;
		p->m_thread.join();
	}
	p->m_cancel = false;
	p->m_thread = std::thread(
		[this, version_check_url]() {
			p->version_check(version_check_url);
		});
}

void AppUpdater::cancel()
{
	p->m_cancel = true;
}
bool AppUpdater::cancel_callback()
{
	cancel();
	return true;
}

std::string AppUpdater::get_default_dest_folder()
{
	return p->m_default_dest_folder.string();
}

std::string AppUpdater::get_filename_from_url(const std::string& url)
{
	size_t slash = url.rfind('/');
	return (slash != std::string::npos ? url.substr(slash + 1) : url);
}

std::string AppUpdater::get_file_extension_from_url(const std::string& url)
{
	size_t dot = url.rfind('.');
	return (dot != std::string::npos ? url.substr(dot) : url);
}

void AppUpdater::set_app_data(DownloadAppData data)
{
	p->set_app_data(std::move(data));
}

DownloadAppData AppUpdater::get_app_data()
{
	return p->get_app_data();
}


} //namespace Slic3r 
