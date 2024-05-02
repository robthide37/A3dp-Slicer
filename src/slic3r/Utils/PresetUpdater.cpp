#include "PresetUpdater.hpp"

#include <algorithm>
#include <thread>
#include <unordered_map>
#include <ostream>
#include <utility>
#include <stdexcept>
#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp> 

#include <wx/app.h>
#include <wx/msgdlg.h>

#include "libslic3r/libslic3r.h"
#include "libslic3r/format.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/I18N.hpp"
#include "slic3r/GUI/UpdateDialogs.hpp"
#include "slic3r/GUI/ConfigWizard.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/format.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/Utils/Http.hpp"
#include "slic3r/Config/Version.hpp"
#include "slic3r/Config/Snapshot.hpp"

namespace fs = boost::filesystem;
using Slic3r::GUI::Config::Index;
using Slic3r::GUI::Config::Version;
using Slic3r::GUI::Config::Snapshot;
using Slic3r::GUI::Config::SnapshotDB;



// FIXME: Incompat bundle resolution doesn't deal with inherited user presets


namespace Slic3r {


static const char *INDEX_FILENAME = "index.idx";
static const char *TMP_EXTENSION = ".download";


void copy_file_fix(const fs::path &source, const fs::path &target)
{
	BOOST_LOG_TRIVIAL(debug) << format("PresetUpdater: Copying %1% -> %2%", source, target);
	std::string error_message;
	CopyFileResult cfr = copy_file(source.string(), target.string(), error_message, false);
	if (cfr != CopyFileResult::SUCCESS) {
		BOOST_LOG_TRIVIAL(error) << "Copying failed(" << cfr << "): " << error_message;
		throw Slic3r::CriticalException(GUI::format(
				_L("Copying of file %1% to %2% failed: %3%"),
				source, target, error_message));
	}
	// Permissions should be copied from the source file by copy_file(). We are not sure about the source
	// permissions, let's rewrite them with 644.
	static constexpr const auto perms = fs::owner_read | fs::owner_write | fs::group_read | fs::others_read;
	fs::permissions(target, perms);
}

struct Update
{
	fs::path source;
	fs::path target;

	Version version;
	std::string vendor;
	std::string changelog_url;

	bool forced_update;

	Update() {}
	Update(fs::path &&source, fs::path &&target, const Version &version, std::string vendor, std::string changelog_url, bool forced = false)
		: source(std::move(source))
		, target(std::move(target))
		, version(version)
		, vendor(std::move(vendor))
		, changelog_url(std::move(changelog_url))
		, forced_update(forced)
	{}

	void install() const
	{
		copy_file_fix(source, target);
	}

	friend std::ostream& operator<<(std::ostream& os, const Update &self)
	{
		os << "Update(" << self.source.string() << " -> " << self.target.string() << ')';
		return os;
	}
};

struct Incompat
{
	fs::path bundle;
	Version version;
	std::string vendor;

	Incompat(fs::path &&bundle, const Version &version, std::string vendor)
		: bundle(std::move(bundle))
		, version(version)
		, vendor(std::move(vendor))
	{}

	void remove() {
		// Remove the bundle file
		fs::remove(bundle);

		// Look for an installed index and remove it too if any
		const fs::path installed_idx = bundle.replace_extension("idx");
		if (fs::exists(installed_idx)) {
			fs::remove(installed_idx);
		}
	}

	friend std::ostream& operator<<(std::ostream& os , const Incompat &self) {
		os << "Incompat(" << self.bundle.string() << ')';
		return os;
	}
};

struct Updates
{
	std::vector<Incompat> incompats;
	std::vector<Update> updates;
};


wxDEFINE_EVENT(EVT_SLIC3R_VERSION_ONLINE, wxCommandEvent);
wxDEFINE_EVENT(EVT_SLIC3R_EXPERIMENTAL_VERSION_ONLINE, wxCommandEvent);

struct PresetUpdater::priv
{
	std::vector<Index> index_db;

	bool enabled_version_check;
	bool enabled_config_update;
	std::string version_check_url;

	fs::path cache_path;
	fs::path rsrc_path;
	fs::path vendor_path;

	bool cancel;
	std::thread thread;

	bool has_waiting_updates { false };
	Updates waiting_updates;

	priv();

	void set_download_prefs(AppConfig *app_config);
	bool get_file(const std::string &url, const fs::path &target_path) const;
	void prune_tmps() const;
	void sync_version() const;
	void parse_version_string(const std::string& body) const;
	void sync_config(const VendorMap vendors);

	void check_install_indices() const;
	Updates get_config_updates(const Semver& old_slic3r_version) const;
	bool perform_updates(Updates &&updates, bool snapshot = true) const;
	void set_waiting_updates(Updates u);
};

PresetUpdater::priv::priv()
	: cache_path(fs::path(Slic3r::data_dir()) / "cache")
	, rsrc_path(fs::path(resources_dir()) / "profiles")
	, vendor_path(fs::path(Slic3r::data_dir()) / "vendor")
	, cancel(false)
{
	set_download_prefs(GUI::wxGetApp().app_config.get());
	// Install indicies from resources. Only installs those that are either missing or older than in resources.
	check_install_indices();
	// Load indices from the cache directory.
	index_db = Index::load_db();
}

// Pull relevant preferences from AppConfig
void PresetUpdater::priv::set_download_prefs(AppConfig *app_config)
{
	enabled_version_check = app_config->get("notify_release") != "none";
	version_check_url = app_config->version_check_url();
	enabled_config_update = app_config->get("preset_update") == "1" && !app_config->legacy_datadir();
}

// Downloads a file (http get operation). Cancels if the Updater is being destroyed.
bool PresetUpdater::priv::get_file(const std::string &url, const fs::path &target_path) const
{
	bool res = false;
	fs::path tmp_path = target_path;
	tmp_path += format(".%1%%2%", get_current_pid(), TMP_EXTENSION);

	BOOST_LOG_TRIVIAL(info) << format("Get: `%1%`\n\t-> `%2%`\n\tvia tmp path `%3%`",
		url,
		target_path.string(),
		tmp_path.string());

	Http::get(url)
        .on_progress([](Http::Progress, bool &cancel) {
			if (cancel) { cancel = true; }
		})
		.on_error([&](std::string body, std::string error, unsigned http_status) {
			(void)body;
			BOOST_LOG_TRIVIAL(error) << format("Error getting: `%1%`: HTTP %2%, %3%",
				url,
				http_status,
				error);
		})
		.on_complete([&](std::string body, unsigned /* http_status */) {
			fs::fstream file(tmp_path, std::ios::out | std::ios::binary | std::ios::trunc);
			file.write(body.c_str(), body.size());
			file.close();
			fs::rename(tmp_path, target_path);
			res = true;
		})
		.perform_sync();

	return res;
}

// Remove leftover paritally downloaded files, if any.
void PresetUpdater::priv::prune_tmps() const
{
    for (auto &dir_entry : boost::filesystem::directory_iterator(cache_path))
		if (is_plain_file(dir_entry) && dir_entry.path().extension() == TMP_EXTENSION) {
			BOOST_LOG_TRIVIAL(debug) << "Cache prune: " << dir_entry.path().string();
			fs::remove(dir_entry.path());
		}
}

//parse the string, if it doesn't contain a valid version string, return invalid version.
Semver get_version(const std::string &str, const std::regex &regexp) {
	std::smatch match;
	if (std::regex_match(str, match, regexp)) {
		std::string version_cleaned = match[0];
		const boost::optional<Semver> version = Semver::parse(version_cleaned);
		if (version.has_value()) {
			return *version;
		}
	}
	return Semver::invalid();
}

// Get Slic3rPE version available online, save in AppConfig.
void PresetUpdater::priv::sync_version() const
{
	if (! enabled_version_check) { return; }

	BOOST_LOG_TRIVIAL(info) << format("Downloading %1% online version from: `%2%`", SLIC3R_APP_NAME, version_check_url);

	Http::get(version_check_url)
		.on_progress([this](Http::Progress, bool &cancel) {
			cancel = this->cancel;
		})
		.on_error([&](std::string body, std::string error, unsigned http_status) {
			(void)body;
			BOOST_LOG_TRIVIAL(error) << format("Error getting: `%1%`: HTTP %2%, %3%",
				version_check_url,
				http_status,
				error);
		})
		.on_complete([&](std::string body, unsigned /* http_status */) {

			parse_version_string(body);
		})
		.perform_sync();
}

// Parses version string obtained in sync_version() and sends events to UI thread.
// Version string must contain release version on first line. Follows non-mandatory alpha / beta releases on following lines (alpha=2.0.0-alpha1).
void PresetUpdater::priv::parse_version_string(const std::string& constbody) const
{
/* PRUSASLICER
	std::string body = constbody;
	boost::trim(body);
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
*/

// github verison

	boost::property_tree::ptree root;
	std::stringstream json_stream(constbody);
	boost::property_tree::read_json(json_stream, root);
	bool i_am_pre = false;
	//at least two number, use '.' as separator. can be followed by -Az23 for prereleased and +Az42 for metadata
	std::regex matcher("[0-9]+\\.[0-9]+(\\.[0-9]+)*(-[A-Za-z0-9]+)?(\\+[A-Za-z0-9]+)?");

	Semver current_version(SLIC3R_VERSION_FULL);
	Semver best_pre(1,0,0,0);
	Semver best_release(1, 0, 0, 0);
	std::string best_pre_url;
	std::string best_release_url;
	const std::regex reg_num("([0-9]+)");
	for (auto json_version : root) {
		std::string tag = json_version.second.get<std::string>("tag_name");
		for (std::regex_iterator it = std::sregex_iterator(tag.begin(), tag.end(), reg_num); it != std::sregex_iterator(); ++it) {

		}
		Semver tag_version = get_version(tag, matcher);
		if (current_version == tag_version)
			i_am_pre = json_version.second.get<bool>("prerelease");
		if (json_version.second.get<bool>("prerelease")) {
			if (best_pre < tag_version) {
				best_pre = tag_version;
				best_pre_url = json_version.second.get<std::string>("html_url");
			}
		} else {
			if (best_release < tag_version) {
				best_release = tag_version;
				best_release_url = json_version.second.get<std::string>("html_url");
			}
		}
	}

	//if release is more recent than beta, use release anyway
	if (best_pre < best_release) {
		best_pre = best_release;
		best_pre_url = best_release_url;
	}
	//if we're the most recent, don't do anything
	if ((i_am_pre ? best_pre : best_release) <= current_version)
		return;

	BOOST_LOG_TRIVIAL(info) << format("Got %1% online version: `%2%`. Sending to GUI thread...", SLIC3R_APP_NAME, i_am_pre? best_pre:best_release);

	wxCommandEvent* evt = new wxCommandEvent(EVT_SLIC3R_VERSION_ONLINE);
	evt->SetString((i_am_pre ? best_pre : best_release).to_string());
	GUI::wxGetApp().QueueEvent(evt);
}

// Download vendor indices. Also download new bundles if an index indicates there's a new one available.
// Both are saved in cache.
void PresetUpdater::priv::sync_config(const VendorMap vendors)
{
	BOOST_LOG_TRIVIAL(info) << "Syncing configuration cache";

	if (!enabled_config_update) { return; }

	// Donwload vendor preset bundles
	// Over all indices from the cache directory:
	for (auto &index : index_db) {
		if (cancel) { return; }

		const auto vendor_it = vendors.find(index.vendor());
		if (vendor_it == vendors.end()) {
			BOOST_LOG_TRIVIAL(warning) << "No such vendor: " << index.vendor();
			continue;
		}

		const VendorProfile &vendor = vendor_it->second;
		if (vendor.config_update_url.empty()) {
			BOOST_LOG_TRIVIAL(info) << "Vendor has no config_update_url: " << vendor.name;
			continue;
		}

		// Download a fresh index
		BOOST_LOG_TRIVIAL(info) << "Downloading index for vendor: " << vendor.name;
		const auto idx_url = vendor.config_update_url + "/" + INDEX_FILENAME;
		const std::string idx_path = (cache_path / (vendor.id + ".idx")).string();
		const std::string idx_path_temp = idx_path + "-update";
		//check if idx_url is leading to a safe site 
		//if (! boost::starts_with(idx_url, "http://files.my_company.com/wp-content/uploads/repository/")
		//	&& ! boost::starts_with(idx_url, "https://files.my_company.com/wp-content/uploads/repository/"))
		//{
		//	BOOST_LOG_TRIVIAL(warning) << "unsafe url path for vendor \"" << vendor.name << "\" rejected: " << idx_url;
		//	continue;
		//}
		if (!get_file(idx_url, idx_path_temp)) { continue; }
		if (cancel) { return; }

		// Load the fresh index up
		{
			Index new_index;
			try {
				new_index.load(idx_path_temp);
			} catch (const std::exception & /* err */) {
				BOOST_LOG_TRIVIAL(error) << format("Could not load downloaded index %1% for vendor %2%: invalid index?", idx_path_temp, vendor.name);
				continue;
			}
			if (new_index.version() < index.version()) {
				BOOST_LOG_TRIVIAL(warning) << format("The downloaded index %1% for vendor %2% is older than the active one. Ignoring the downloaded index.", idx_path_temp, vendor.name);
				continue;
			}
			Slic3r::rename_file(idx_path_temp, idx_path);
			//if we rename path we need to change it in Index object too or create the object again
			//index = std::move(new_index);
			try {
				index.load(idx_path);
			}
			catch (const std::exception& /* err */) {
				BOOST_LOG_TRIVIAL(error) << format("Could not load downloaded index %1% for vendor %2%: invalid index?", idx_path, vendor.name);
				continue;
			}
			if (cancel)
				return;
		}

		// See if a there's a new version to download
		const auto recommended_it = index.recommended();
		if (recommended_it == index.end()) {
			BOOST_LOG_TRIVIAL(error) << format("No recommended version for vendor: %1%, invalid index?", vendor.name);
			continue;
		}

		const auto recommended = recommended_it->config_version;

		BOOST_LOG_TRIVIAL(debug) << format("Got index for vendor: %1%: current version: %2%, recommended version: %3%",
			vendor.name,
			vendor.config_version.to_string(),
			recommended.to_string());

		if (vendor.config_version >= recommended) { continue; }

		// Download a fresh bundle
		BOOST_LOG_TRIVIAL(info) << "Downloading new bundle for vendor: " << vendor.name;
		const std::string bundle_url = format("%1%/%2%.ini", vendor.config_update_url, recommended.to_string());
		const fs::path bundle_path = cache_path / (vendor.id + ".ini");
		if (! get_file(bundle_url, bundle_path)) { continue; }
		if (cancel) { return; }
	}
}

// Install indicies from resources. Only installs those that are either missing or older than in resources.
void PresetUpdater::priv::check_install_indices() const
{
	BOOST_LOG_TRIVIAL(info) << "Checking if indices need to be installed from resources...";
	if (!fs::exists(rsrc_path))
		return;
	for (auto& dir_entry : boost::filesystem::directory_iterator(rsrc_path)) {
		if (is_idx_file(dir_entry)) {
			const auto& path = dir_entry.path();
			const auto path_in_cache = cache_path / path.filename();

			if (!fs::exists(path_in_cache)) {
				BOOST_LOG_TRIVIAL(info) << "Install index from resources: " << path.filename();
				copy_file_fix(path, path_in_cache);
			} else {
				Index idx_rsrc, idx_cache;
				idx_rsrc.load(path);
				idx_cache.load(path_in_cache);
				fs::path bundle_path = vendor_path / (idx_cache.vendor() + ".ini");

				//test if the cache file is bad and the resource one is good.
				if (fs::exists(bundle_path)) {
					Semver version = VendorProfile::from_ini(bundle_path, false).config_version;
					const auto ver_from_cache = idx_cache.find(version);
					const auto ver_from_resource = idx_rsrc.find(version);
					if (ver_from_resource != idx_rsrc.end()) {
						if (idx_cache.version() < idx_rsrc.version()) {
							if (fs::exists(bundle_path)) {
								BOOST_LOG_TRIVIAL(info) << "Update index from resources (new version): " << path.filename();
								copy_file_fix(path, path_in_cache);
							}
						} else if (ver_from_cache == idx_cache.end()) {
							BOOST_LOG_TRIVIAL(info) << "Update index from resources (only way to have a consistent idx): " << path.filename();
							copy_file_fix(path, path_in_cache);
						}
					}
				} else if (idx_cache.version() < idx_rsrc.version() || idx_cache.configs().back().max_slic3r_version < idx_rsrc.configs().back().max_slic3r_version) {
					//not installed, force update the .idx from resource
					BOOST_LOG_TRIVIAL(info) << "Update index from resources (uninstalled & more up-to-date): " << path.filename();
					copy_file_fix(path, path_in_cache);
				}
			}
		}
	}
}

// Generates a list of bundle updates that are to be performed.
// Version of slic3r that was running the last time and which was read out from PrusaSlicer.ini is provided
// as a parameter.
Updates PresetUpdater::priv::get_config_updates(const Semver &old_slic3r_version) const
{
	Updates updates;

	BOOST_LOG_TRIVIAL(info) << "Checking for cached configuration updates...";

	// Over all indices from the cache directory:
    for (const Index& idx : index_db) {
		auto bundle_path = vendor_path / (idx.vendor() + ".ini");
		auto bundle_path_idx = vendor_path / idx.path().filename();

		if (! fs::exists(bundle_path)) {
			BOOST_LOG_TRIVIAL(info) << format("Confing bundle not installed for vendor %1%, skipping: ", idx.vendor());
			continue;
		}

		// Perform a basic load and check the version of the installed preset bundle.
		VendorProfile vp = VendorProfile::from_ini(bundle_path, false);

		// Getting a recommended version from the latest index, wich may have been downloaded
		// from the internet, or installed / updated from the installation resources.
		auto recommended = idx.recommended();
		if (recommended == idx.end()) {
			BOOST_LOG_TRIVIAL(error) << format("No recommended version for vendor: %1%, invalid index? Giving up.", idx.vendor());
			// XXX: what should be done here?
			continue;
		}

		const auto ver_current = idx.find(vp.config_version);
		const bool ver_current_found = ver_current != idx.end();

		BOOST_LOG_TRIVIAL(debug) << format("Vendor: %1%, version installed: %2%%3%, version cached: %4%",
			vp.name,
			vp.config_version.to_string(),
			(ver_current_found ? "" : " (not found in index!)"),
			recommended->config_version.to_string());

		if (! ver_current_found) {
			// Config bundle inside the resources directory.
			fs::path path_in_rsrc = rsrc_path / (idx.vendor() + ".ini");
			fs::path path_idx_in_rsrc = rsrc_path / (idx.vendor() + ".idx");
			if (fs::exists(path_idx_in_rsrc)) {
				Index rsrc_idx;
				rsrc_idx.load(path_idx_in_rsrc);
			
				// Any published config shall be always found in the latest config index.
				std::string message = format("Preset bundle `%1%` version not found in index: %2%, do we force the update to the version %3%? ", idx.vendor(), vp.config_version.to_string(), rsrc_idx.version().to_string());
				wxMessageDialog msg_wingow(nullptr, message, wxString(SLIC3R_APP_NAME " - ") + (_L("Notice")), wxYES | wxNO | wxICON_INFORMATION);
				if (msg_wingow.ShowModal() == wxID_YES) {
					//copy idx
					copy_file_fix(path_idx_in_rsrc, idx.path());
					//copy profile
					copy_file_fix(path_in_rsrc, bundle_path);
				}
				continue;
			} else {
				// Any published config shall be always found in the latest config index.
				std::string message = format("Preset bundle `%1%` version not found in index: %2%", idx.vendor(), vp.config_version.to_string());
				BOOST_LOG_TRIVIAL(error) << message;
				GUI::show_error(nullptr, message);
			}
			continue;
		}

		bool current_not_supported = false; //if slcr is incompatible but situation is not downgrade, we do forced updated and this bool is information to do it 

		if (ver_current_found && !ver_current->is_current_slic3r_supported()){
			if(ver_current->is_current_slic3r_downgrade()) {
				// "Reconfigure" situation.
				BOOST_LOG_TRIVIAL(warning) << "Current Slic3r incompatible with installed bundle: " << bundle_path.string();
				updates.incompats.emplace_back(std::move(bundle_path), *ver_current, vp.name);
				continue;
			}
		current_not_supported = true;
		}

		if (recommended->config_version < vp.config_version) {
			BOOST_LOG_TRIVIAL(warning) << format("Recommended config version for the currently running %1% is older than the currently installed config for vendor %2%. This should not happen.", SLIC3R_APP_NAME, idx.vendor());
			continue;
		}

		if (recommended->config_version == vp.config_version) {
			// The recommended config bundle is already installed.
			continue;
		}

		// Config bundle update situation. The recommended config bundle version for this PrusaSlicer version from the index from the cache is newer
		// than the version of the currently installed config bundle.

		// The config index inside the cache directory (given by idx.path()) is one of the following:
		// 1) The last config index downloaded by any previously running PrusaSlicer instance
		// 2) The last config index installed by any previously running PrusaSlicer instance (older or newer) from its resources.
		// 3) The last config index installed by the currently running PrusaSlicer instance from its resources.
		// The config index is always the newest one (given by its newest config bundle referenced), and older config indices shall fully contain
		// the content of the older config indices.

		// Config bundle inside the cache directory.
		fs::path path_in_cache 		= cache_path / (idx.vendor() + ".ini");
		// Config bundle inside the resources directory.
		fs::path path_in_rsrc 		= rsrc_path  / (idx.vendor() + ".ini");
		// Config index inside the resources directory.
		fs::path path_idx_in_rsrc 	= rsrc_path  / (idx.vendor() + ".idx");

		// Search for a valid config bundle in the cache directory.
		bool 		found = false;
		Update    	new_update;
		fs::path 	bundle_path_idx_to_install;
		if (fs::exists(path_in_cache)) {
			try {
				VendorProfile new_vp = VendorProfile::from_ini(path_in_cache, false);
				if (new_vp.config_version == recommended->config_version) {
					// The config bundle from the cache directory matches the recommended version of the index from the cache directory.
					// This is the newest known recommended config. Use it.
					new_update = Update(std::move(path_in_cache), std::move(bundle_path), *recommended, vp.name, vp.changelog_url, current_not_supported);
					// and install the config index from the cache into vendor's directory.
					bundle_path_idx_to_install = idx.path();
					found = true;
				}
			} catch (const std::exception &ex) {
				BOOST_LOG_TRIVIAL(info) << format("Failed to load the config bundle `%1%`: %2%", path_in_cache.string(), ex.what());
			}
		}

		// Keep the rsrc_idx outside of the next block, as we will reference the "recommended" version by an iterator.
		Index rsrc_idx;
		if (! found && fs::exists(path_in_rsrc) && fs::exists(path_idx_in_rsrc)) {
			// Trying the config bundle from resources (from the installation).
			// In that case, the recommended version number has to be compared against the recommended version reported by the config index from resources as well, 
			// as the config index in the cache directory may already be newer, recommending a newer config bundle than available in cache or resources.
			VendorProfile rsrc_vp;
			try {
				rsrc_vp = VendorProfile::from_ini(path_in_rsrc, false);
			} catch (const std::exception &ex) {
				BOOST_LOG_TRIVIAL(info) << format("Cannot load the config bundle at `%1%`: %2%", path_in_rsrc.string(), ex.what());
			}
			if (rsrc_vp.valid()) {
				try {
					rsrc_idx.load(path_idx_in_rsrc);
				} catch (const std::exception &ex) {
					BOOST_LOG_TRIVIAL(info) << format("Cannot load the config index at `%1%`: %2%", path_idx_in_rsrc.string(), ex.what());
				}
				recommended = rsrc_idx.recommended();
				if (recommended != rsrc_idx.end() && recommended->config_version == rsrc_vp.config_version && recommended->config_version > vp.config_version) {
					new_update = Update(std::move(path_in_rsrc), std::move(bundle_path), *recommended, vp.name, vp.changelog_url, current_not_supported);
					bundle_path_idx_to_install = path_idx_in_rsrc;
					found = true;
				} else {
					BOOST_LOG_TRIVIAL(warning) << format("The recommended config version for vendor `%1%` in resources does not match the recommended\n"
			                                             " config version for this version of `%2%`. Corrupted installation?", idx.vendor(), SLIC3R_APP_NAME);
				}
			}
		}

		if (found) {
			// Load 'installed' idx, if any.
			// 'Installed' indices are kept alongside the bundle in the `vendor` subdir
			// for bookkeeping to remember a cancelled update and not offer it again.
			if (fs::exists(bundle_path_idx)) {
				Index existing_idx;
				try {
					existing_idx.load(bundle_path_idx);
					// Find a recommended config bundle version for the slic3r version last executed. This makes sure that a config bundle update will not be missed
					// when upgrading an application. On the other side, the user will be bugged every time he will switch between slic3r versions.
                    /*const auto existing_recommended = existing_idx.recommended(old_slic3r_version);
                    if (existing_recommended != existing_idx.end() && recommended->config_version == existing_recommended->config_version) {
						// The user has already seen (and presumably rejected) this update
						BOOST_LOG_TRIVIAL(info) << format("Downloaded index for `%1%` is the same as installed one, not offering an update.",idx.vendor());
						continue;
					}*/
				} catch (const std::exception &err) {
					BOOST_LOG_TRIVIAL(error) << format("Cannot load the installed index at `%1%`: %2%", bundle_path_idx, err.what());
				}
			}

			// Check if the update is already present in a snapshot
			if(!current_not_supported)
			{
				const auto recommended_snap = SnapshotDB::singleton().snapshot_with_vendor_preset(vp.name, recommended->config_version);
				if (recommended_snap != SnapshotDB::singleton().end()) {
					BOOST_LOG_TRIVIAL(info) << format("Bundle update %1% %2% already found in snapshot %3%, skipping...",
						vp.name,
						recommended->config_version.to_string(),
						recommended_snap->id);
					continue;
				}
			}

			updates.updates.emplace_back(std::move(new_update));
			// 'Install' the index in the vendor directory. This is used to memoize
			// offered updates and to not offer the same update again if it was cancelled by the user.
			copy_file_fix(bundle_path_idx_to_install, bundle_path_idx);
		} else {
			BOOST_LOG_TRIVIAL(warning) << format("Index for vendor %1% indicates update (%2%) but the new bundle was found neither in cache nor resources",
				idx.vendor(),
				recommended->config_version.to_string());
		}
	}

	return updates;
}

bool PresetUpdater::priv::perform_updates(Updates &&updates, bool snapshot) const
{
	if (updates.incompats.size() > 0) {
		if (snapshot) {
			BOOST_LOG_TRIVIAL(info) << "Taking a snapshot...";
			if (! GUI::Config::take_config_snapshot_cancel_on_error(*GUI::wxGetApp().app_config, Snapshot::SNAPSHOT_DOWNGRADE, "",
				_u8L("Continue and install configuration updates?")))
				return false;
		}
		
		BOOST_LOG_TRIVIAL(info) << format("Deleting %1% incompatible bundles", updates.incompats.size());

		for (auto &incompat : updates.incompats) {
			BOOST_LOG_TRIVIAL(info) << '\t' << incompat;
			incompat.remove();
		}

		
	} else if (updates.updates.size() > 0) {
		
		if (snapshot) {
			BOOST_LOG_TRIVIAL(info) << "Taking a snapshot...";
			if (! GUI::Config::take_config_snapshot_cancel_on_error(*GUI::wxGetApp().app_config, Snapshot::SNAPSHOT_UPGRADE, "",
				_u8L("Continue and install configuration updates?")))
				return false;
		}

		BOOST_LOG_TRIVIAL(info) << format("Performing %1% updates", updates.updates.size());

		for (const auto &update : updates.updates) {
			BOOST_LOG_TRIVIAL(info) << '\t' << update;

			update.install();

			PresetBundle bundle;
			// Throw when parsing invalid configuration. Only valid configuration is supposed to be provided over the air.
			bundle.load_configbundle(update.source.string(), PresetBundle::LoadConfigBundleAttribute::LoadSystem, ForwardCompatibilitySubstitutionRule::Disable);

			BOOST_LOG_TRIVIAL(info) << format("Deleting %1% conflicting presets", bundle.fff_prints.size() + bundle.filaments.size() + bundle.printers.size());

			auto preset_remover = [](const Preset &preset) {
				BOOST_LOG_TRIVIAL(info) << '\t' << preset.file;
				fs::remove(preset.file);
			};

			for (const auto &preset : bundle.fff_prints){ preset_remover(preset); }
			for (const auto &preset : bundle.filaments) { preset_remover(preset); }
			for (const auto &preset : bundle.printers)  { preset_remover(preset); }

			// Also apply the `obsolete_presets` property, removing obsolete ini files

			BOOST_LOG_TRIVIAL(info) << format("Deleting %1% obsolete presets",
				bundle.obsolete_presets.fff_prints.size() + bundle.obsolete_presets.filaments.size() + bundle.obsolete_presets.printers.size());

			auto obsolete_remover = [](const char *subdir, const std::string &preset) {
				auto path = fs::path(Slic3r::data_dir()) / subdir / preset;
				path += ".ini";
				BOOST_LOG_TRIVIAL(info) << '\t' << path.string();
				fs::remove(path);
			};

			for (const auto &name : bundle.obsolete_presets.fff_prints)    { obsolete_remover("print", name); }
			for (const auto &name : bundle.obsolete_presets.filaments) { obsolete_remover("filament", name); }
			for (const auto &name : bundle.obsolete_presets.sla_prints) { obsolete_remover("sla_print", name); } 
			for (const auto &name : bundle.obsolete_presets.sla_materials/*filaments*/) { obsolete_remover("sla_material", name); } 
			for (const auto &name : bundle.obsolete_presets.printers)  { obsolete_remover("printer", name); }
		}
	}

	return true;
}

void PresetUpdater::priv::set_waiting_updates(Updates u)
{
	waiting_updates = u;
	has_waiting_updates = true;
}

PresetUpdater::PresetUpdater() :
	p(new priv())
{}


// Public

PresetUpdater::~PresetUpdater()
{
	if (p && p->thread.joinable()) {
		// This will stop transfers being done by the thread, if any.
		// Cancelling takes some time, but should complete soon enough.
		p->cancel = true;
		p->thread.join();
	}
}

void PresetUpdater::sync(PresetBundle *preset_bundle)
{
	p->set_download_prefs(GUI::wxGetApp().app_config.get());
	if (!p->enabled_version_check && !p->enabled_config_update) { return; }

	// Copy the whole vendors data for use in the background thread
	// Unfortunatelly as of C++11, it needs to be copied again
	// into the closure (but perhaps the compiler can elide this).
	VendorMap vendors = preset_bundle->vendors;

    p->thread = std::thread([this, vendors]() {
		this->p->prune_tmps();
		this->p->sync_version();
		this->p->sync_config(std::move(vendors));
    });
}

void PresetUpdater::slic3r_update_notify()
{
	if (! p->enabled_version_check)
		return;
	auto* app_config = GUI::wxGetApp().app_config.get();
	const auto ver_online_str = app_config->get("version_online");
	const auto ver_online = Semver::parse(ver_online_str);
	const auto ver_online_seen = Semver::parse(app_config->get("version_online_seen"));

	if (ver_online) {
		// Only display the notification if the version available online is newer AND if we haven't seen it before
		if (*ver_online > Slic3r::SEMVER && (! ver_online_seen || *ver_online_seen < *ver_online)) {
			GUI::MsgUpdateSlic3r notification(Slic3r::SEMVER, *ver_online);
			notification.ShowModal();
			if (notification.disable_version_check()) {
				app_config->set("notify_release", "none");
				p->enabled_version_check = false;
			}
		}

		app_config->set("version_online_seen", ver_online_str);
	}
}

static bool reload_configs_update_gui()
{
	wxString header = _L("Configuration Updates causes a lost of preset modification.\n"
						 "So, check unsaved changes and save them if necessary.");
	if (!GUI::wxGetApp().check_and_save_current_preset_changes(_L("Updating"), header, false ))
		return false;

	// Reload global configuration
	auto* app_config = GUI::wxGetApp().app_config.get();
	// System profiles should not trigger any substitutions, user profiles may trigger substitutions, but these substitutions
	// were already presented to the user on application start up. Just do substitutions now and keep quiet about it.
	// However throw on substitutions in system profiles, those shall never happen with system profiles installed over the air.
	GUI::wxGetApp().preset_bundle->load_presets(*app_config, ForwardCompatibilitySubstitutionRule::EnableSilentDisableSystem);
	GUI::wxGetApp().load_current_presets();
	GUI::wxGetApp().plater()->set_bed_shape();

	return true;
}

PresetUpdater::UpdateResult PresetUpdater::config_update(const Semver& old_slic3r_version, UpdateParams params) const
{
	if (! p->enabled_config_update) { return R_NOOP; }

	auto updates = p->get_config_updates(old_slic3r_version);
	if (updates.incompats.size() > 0) {
		BOOST_LOG_TRIVIAL(info) << format("%1% bundles incompatible. Asking for action...", updates.incompats.size());

		std::unordered_map<std::string, wxString> incompats_map;
		for (const auto &incompat : updates.incompats) {
			const auto min_slic3r = incompat.version.min_slic3r_version;
			const auto max_slic3r = incompat.version.max_slic3r_version;
			wxString restrictions;
			if (min_slic3r != Semver::zero() && max_slic3r != Semver::inf()) {
                restrictions = GUI::format_wxstr(_L("requires min. %s and max. %s"),
                    min_slic3r.to_string(),
                    max_slic3r.to_string());
			} else if (min_slic3r != Semver::zero()) {
				restrictions = GUI::format_wxstr(_L("requires min. %s"), min_slic3r.to_string());
				BOOST_LOG_TRIVIAL(debug) << "Bundle is not downgrade, user will now have to do whole wizard. This should not happen.";
			} else {
                restrictions = GUI::format_wxstr(_L("requires max. %s"), max_slic3r.to_string());
			}

			incompats_map.emplace(std::make_pair(incompat.vendor, std::move(restrictions)));
		}

		GUI::MsgDataIncompatible dlg(std::move(incompats_map));
		const auto res = dlg.ShowModal();
		if (res == wxID_REPLACE) {
			BOOST_LOG_TRIVIAL(info) << "User wants to re-configure...";

			// This effectively removes the incompatible bundles:
			// (snapshot is taken beforehand)
			if (! p->perform_updates(std::move(updates)) ||
				! GUI::wxGetApp().run_wizard(GUI::ConfigWizard::RR_DATA_INCOMPAT))
				return R_INCOMPAT_EXIT;

			return R_INCOMPAT_CONFIGURED;
		}
		else {
			BOOST_LOG_TRIVIAL(info) << "User wants to exit Slic3r, bye...";
			return R_INCOMPAT_EXIT;
		}

	} else if (updates.updates.size() > 0) {

		bool incompatible_version = false;
		for (const auto& update : updates.updates) {
			incompatible_version = (update.forced_update ? true : incompatible_version);
			//td::cout << update.forced_update << std::endl;
			//BOOST_LOG_TRIVIAL(info) << format("Update requires higher version.");
		}

		//forced update
		if (incompatible_version)
		{
			BOOST_LOG_TRIVIAL(info) << format("Update of %1% bundles available. At least one requires higher version of Slicer.", updates.updates.size());

			std::vector<GUI::MsgUpdateForced::Update> updates_msg;
			for (const auto& update : updates.updates) {
				std::string changelog_url = update.version.config_version.prerelease() == nullptr ? update.changelog_url : std::string();
				updates_msg.emplace_back(update.vendor, update.version.config_version, update.version.comment, std::move(changelog_url));
			}

			GUI::MsgUpdateForced dlg(updates_msg);

			const auto res = dlg.ShowModal();
			if (res == wxID_OK) {
				BOOST_LOG_TRIVIAL(info) << "User wants to update...";
				if (! p->perform_updates(std::move(updates)) ||
					! reload_configs_update_gui())
					return R_INCOMPAT_EXIT;
				return R_UPDATE_INSTALLED;
			}
			else {
				BOOST_LOG_TRIVIAL(info) << "User wants to exit Slic3r, bye...";
				return R_INCOMPAT_EXIT;
			}
		}

		// regular update
		if (params == UpdateParams::SHOW_NOTIFICATION) {
			p->set_waiting_updates(updates);
			GUI::wxGetApp().plater()->get_notification_manager()->push_notification(GUI::NotificationType::PresetUpdateAvailable);
		}
		else {
			BOOST_LOG_TRIVIAL(info) << format("Update of %1% bundles available. Asking for confirmation ...", p->waiting_updates.updates.size());

			std::vector<GUI::MsgUpdateConfig::Update> updates_msg;
			for (const auto& update : updates.updates) {
				std::string changelog_url = update.version.config_version.prerelease() == nullptr ? update.changelog_url : std::string();
				updates_msg.emplace_back(update.vendor, update.version.config_version, update.version.comment, std::move(changelog_url));
			}

			GUI::MsgUpdateConfig dlg(updates_msg, params == UpdateParams::FORCED_BEFORE_WIZARD);

			const auto res = dlg.ShowModal();
			if (res == wxID_OK) {
				BOOST_LOG_TRIVIAL(debug) << "User agreed to perform the update";
				if (! p->perform_updates(std::move(updates)) ||
					! reload_configs_update_gui())
					return R_ALL_CANCELED;
				return R_UPDATE_INSTALLED;
			}
			else {
				BOOST_LOG_TRIVIAL(info) << "User refused the update";
				if (params == UpdateParams::FORCED_BEFORE_WIZARD && res == wxID_CANCEL)
					return R_ALL_CANCELED;
				return R_UPDATE_REJECT;
			}
		}
		
		// MsgUpdateConfig will show after the notificaation is clicked
	} else {
		BOOST_LOG_TRIVIAL(info) << "No configuration updates available.";
	}

	return R_NOOP;
}

bool PresetUpdater::install_bundles_rsrc(std::vector<std::string> bundles, bool snapshot) const
{
	Updates updates;

	BOOST_LOG_TRIVIAL(info) << format("Installing %1% bundles from resources ...", bundles.size());

	for (const auto &bundle : bundles) {
		auto path_in_rsrc = (p->rsrc_path / bundle).replace_extension(".ini");
		auto path_in_vendors = (p->vendor_path / bundle).replace_extension(".ini");
		updates.updates.emplace_back(std::move(path_in_rsrc), std::move(path_in_vendors), Version(), "", "");
	}

	return p->perform_updates(std::move(updates), snapshot);
}

void PresetUpdater::on_update_notification_confirm()
{
	if (!p->has_waiting_updates)
		return;
	BOOST_LOG_TRIVIAL(info) << format("Update of %1% bundles available. Asking for confirmation ...", p->waiting_updates.updates.size());

	std::vector<GUI::MsgUpdateConfig::Update> updates_msg;
	for (const auto& update : p->waiting_updates.updates) {
		std::string changelog_url = update.version.config_version.prerelease() == nullptr ? update.changelog_url : std::string();
		updates_msg.emplace_back(update.vendor, update.version.config_version, update.version.comment, std::move(changelog_url));
	}

	GUI::MsgUpdateConfig dlg(updates_msg);

	const auto res = dlg.ShowModal();
	if (res == wxID_OK) {
		BOOST_LOG_TRIVIAL(debug) << "User agreed to perform the update";
		if (p->perform_updates(std::move(p->waiting_updates)) &&
			reload_configs_update_gui()) {
			p->has_waiting_updates = false;
		}
	}
	else {
		BOOST_LOG_TRIVIAL(info) << "User refused the update";
	}	
}

bool PresetUpdater::version_check_enabled() const
{
	return p->enabled_version_check;
}

}
