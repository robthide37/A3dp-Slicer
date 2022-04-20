#include "ZipperArchiveImport.hpp"

#include "libslic3r/miniz_extension.hpp"
#include "libslic3r/Exception.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>

namespace Slic3r {

static const constexpr char *CONFIG_FNAME  = "config.ini";
static const constexpr char *PROFILE_FNAME = "prusaslicer.ini";

namespace {

boost::property_tree::ptree read_ini(const mz_zip_archive_file_stat &entry,
                                     MZ_Archive                     &zip)
{
    std::string buf(size_t(entry.m_uncomp_size), '\0');

    if (!mz_zip_reader_extract_file_to_mem(&zip.arch, entry.m_filename,
                                           buf.data(), buf.size(), 0))
        throw Slic3r::FileIOError(zip.get_errorstr());

    boost::property_tree::ptree tree;
    std::stringstream           ss(buf);
    boost::property_tree::read_ini(ss, tree);
    return tree;
}

EntryBuffer read_entry(const mz_zip_archive_file_stat &entry,
                       MZ_Archive                     &zip,
                       const std::string              &name)
{
    std::vector<uint8_t> buf(entry.m_uncomp_size);

    if (!mz_zip_reader_extract_file_to_mem(&zip.arch, entry.m_filename,
                                           buf.data(), buf.size(), 0))
        throw Slic3r::FileIOError(zip.get_errorstr());

    return {std::move(buf), (name.empty() ? entry.m_filename : name)};
}

} // namespace

ZipperArchive read_zipper_archive(const std::string &zipfname,
                                  const std::vector<std::string> &includes,
                                  const std::vector<std::string> &excludes)
{
    ZipperArchive arch;

    // Little RAII
    struct Arch : public MZ_Archive
    {
        Arch(const std::string &fname)
        {
            if (!open_zip_reader(&arch, fname))
                throw Slic3r::FileIOError(get_errorstr());
        }

        ~Arch() { close_zip_reader(&arch); }
    } zip(zipfname);

    mz_uint num_entries = mz_zip_reader_get_num_files(&zip.arch);

    for (mz_uint i = 0; i < num_entries; ++i) {
        mz_zip_archive_file_stat entry;

        if (mz_zip_reader_file_stat(&zip.arch, i, &entry)) {
            std::string name = entry.m_filename;
            boost::algorithm::to_lower(name);

            if (!std::any_of(includes.begin(), includes.end(),
                            [&name](const std::string &incl) {
                                return boost::algorithm::contains(name, incl);
                            }))
                continue;

            if (std::any_of(excludes.begin(), excludes.end(),
                            [&name](const std::string &excl) {
                                return boost::algorithm::contains(name, excl);
                            }))
                continue;

            if (name == CONFIG_FNAME) arch.config = read_ini(entry, zip);
            if (name == PROFILE_FNAME) arch.profile = read_ini(entry, zip);

            auto it = std::lower_bound(
                arch.entries.begin(), arch.entries.end(),
                EntryBuffer{{}, name},
                [](const EntryBuffer &r1, const EntryBuffer &r2) {
                    return std::less<std::string>()(r1.fname, r2.fname);
                });

            arch.entries.insert(it, read_entry(entry, zip, name));
        }
    }

    return arch;
}

} // namespace Slic3r
