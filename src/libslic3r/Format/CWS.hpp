///|/ Copyright (c) SuperSlicer Joseph Lenox <lenox.joseph@gmail.com>
///|/
///|/ SuperSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_format_CWS_HPP
#define slic3r_format_CWS_HPP

#include "libslic3r/Format/SL1.hpp"
#include "libslic3r/I18N.hpp"
#include "SLAArchiveWriter.hpp"
#include "SLAArchiveFormatRegistry.hpp"

namespace Slic3r {
// "Masked" CWS as used by Malyan S100
class MaskedCWSArchive : public SL1Archive {
public:
    MaskedCWSArchive() = default;
    explicit MaskedCWSArchive(const SLAPrinterConfig &cfg): SL1Archive(cfg) {}
    explicit MaskedCWSArchive(SLAPrinterConfig &&cfg): SL1Archive(std::move(cfg)) {}
    // Export the print into an archive using the provided filename.
    void export_print(const std::string     fname,
                              const SLAPrint       &print,
                              const ThumbnailsList &thumbnails,
                              const std::string    &projectname = "") override;
};


inline Slic3r::ArchiveEntry masked_cws_format()
{
    Slic3r::ArchiveEntry entry(OutputFormat::ofMaskedCWS ,"CWS");
    entry.desc = L("Masked CWS / Malyan S100");
    entry.ext  = "cws";
    entry.ext_aliases = {"zip"};
    entry.wrfactoryfn = [] (const SLAPrinterConfig &cfg) { return std::make_unique<MaskedCWSArchive>(cfg); };

    return entry;
}
} // namespace Slic3r

#endif // slic3r_format_CWS_HPP
