///|/ Copyright (c) Prusa Research 2018 - 2023 David Kocík @kocikdav, Lukáš Matěna @lukasmatena, Vojtěch Bubník @bubnikv, Oleksandra Iushchenko @YuSanka, Vojtěch Král @vojtechkral
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_PresetUpdate_hpp_
#define slic3r_PresetUpdate_hpp_

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "libslic3r/Semver.hpp"
#include "libslic3r/Preset.hpp"

#include <wx/event.h>

namespace Slic3r {


class AppConfig;
class PresetBundle;
class Semver;
class PresetUpdater;

#define USE_GTHUB_PRESET_UPDATE 1

struct VendorAvailable
{
    Semver config_version;
    Semver slicer_version;
    std::string local_file;
    std::string url_zip;
    std::string commit_sha;
    std::string commit_url;
    std::string tag;
    std::string notes;
};
struct VendorSync
{
    VendorProfile profile;
    bool is_installed = false;
    bool is_synch = false;
    bool synch_in_progress = false;
    bool synch_failed = false;
    bool can_upgrade = false;
    // true if "profile" is inside datadir's vendor directory 
    std::vector<VendorAvailable> available_profiles = {};
    VendorAvailable *best = nullptr;


    bool parse_tags(const std::string &json);
    void sort_available();
    // reurn error message (empty if install succeed)
    std::string install_vendor_config(const VendorAvailable &to_install, PresetUpdater& api_slot);
    bool uninstall_vendor_config();
    void reset(const VendorProfile &profile, bool installed);
};
class PresetUpdater
{
protected:
public:
    bool cancel = false; //TODO
    std::thread thread;

    // max request: 60/h per ip
    // half to allow two config update from same ip at the same time.
    std::atomic_int max_request = 25;
    std::time_t next_time_slot = 0;

    //std::vector<VendorSync> installed_vendors;
    //std::vector<VendorSync> unused_vendors;
    std::recursive_mutex all_vendors_mutex;
    std::map<std::string, VendorSync> all_vendors;

    std::atomic_int profiles_synch = 0;
    //std::atomic_int profiles_to_update = 0;
    std::mutex callback_update_preset_mutex;
    std::function<void(int)> callback_update_preset;
    std::atomic_bool synch_process_ongoing = false;
    bool is_synch = false;

    std::mutex callback_update_changelog_mutex;
    std::atomic_int changelog_synch = 0;

    PresetUpdater();
    PresetUpdater(PresetUpdater &&) = delete;
    PresetUpdater(const PresetUpdater &) = delete;
    PresetUpdater &operator=(PresetUpdater &&) = delete;
    PresetUpdater &operator=(const PresetUpdater &) = delete;
    ~PresetUpdater();

    enum PresetUpdaterResult
    {
        PUR_Nothing,
        PUR_Updates
    };
    
    int get_profile_count_to_update();
    size_t count_installed();
    VendorSync *get_vendor(std::string id) {
        auto it = all_vendors.find(id);
        return it == all_vendors.end() ? nullptr : &it->second;
    }

    bool has_api_request_slot(const std::string &url);

    void set_installed_vendors(const PresetBundle *preset_bundle);
    void reload_all_vendors();

    // If either version check or config updating is enabled, get the appropriate data in the background and cache it.
    // call the handler with PresetUpdaterResult if somethign isntalle
    void sync_async(std::function<void(int)> callback_update_preset, bool force = false);
    void download_logs(const std::string &vendor_id, std::function<void(bool)> callback_result, bool force = false);

    // Show the window to select what to download
    // emit a EVT_CONFIG_UPDATER_SYNC_DONE into evt_handler when done (if not null)
    void show_synch_window(wxWindow *parent,
                           const PresetBundle *preset_bundle,
                           const wxString &message,
                           std::function<void(bool)> callback_dialog_closed);

    void download_new_repo(const std::string &github_org_repo, std::function<void(bool)> callback_result);
    void uninstall_vendor(const std::string &vendor_id, std::function<void(bool)> callback_result);
    void install_vendor(const std::string &vendor_id, const VendorAvailable &version, std::function<void(const std::string &)> callback_result);

    void uninstall_all_vendors(std::function<void(bool)> callback_result);
    void install_all_vendors(std::function<void(const std::string &)> callback_result);
    void upgrade_all_installed_vendors(std::function<void(const std::string &)> callback_result);
protected:
    void load_unused_vendors(std::set<std::string> &vendors_id, const boost::filesystem::path vendor_dir, bool is_installed);

    void update_vendor(VendorSync &vendor, bool force = false);
    // must be call by each update_vendor call at some point.
    void end_updating();
};

wxDECLARE_EVENT(EVT_CONFIG_UPDATER_SYNC_DONE, wxCommandEvent);
}
#endif
