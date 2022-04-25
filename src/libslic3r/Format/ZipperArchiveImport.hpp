#ifndef ZIPPERARCHIVEIMPORT_HPP
#define ZIPPERARCHIVEIMPORT_HPP

#include <vector>
#include <string>
#include <cstdint>

#include <boost/property_tree/ptree.hpp>

#include "libslic3r/PrintConfig.hpp"

namespace Slic3r {

struct EntryBuffer
{
    std::vector<uint8_t> buf;
    std::string          fname;
};

struct ZipperArchive
{
    boost::property_tree::ptree profile, config;
    std::vector<EntryBuffer>    entries;
};

const constexpr char *CONFIG_FNAME  = "config.ini";
const constexpr char *PROFILE_FNAME = "prusaslicer.ini";

ZipperArchive read_zipper_archive(const std::string &zipfname,
                                  const std::vector<std::string> &includes,
                                  const std::vector<std::string> &excludes);

// Extract the print profile form the archive onto 'out'.
// Returns a profile that has correct parameters to use for model reconstruction
// even if the needed parameters were not fully found in the archive's metadata.
// The inout argument shall be a usable fallback profile if the archive
// has totally corrupted metadata.
std::pair<DynamicPrintConfig, ConfigSubstitutions> extract_profile(
    const ZipperArchive &arch, DynamicPrintConfig &inout);

} // namespace Slic3r

#endif // ZIPPERARCHIVEIMPORT_HPP
