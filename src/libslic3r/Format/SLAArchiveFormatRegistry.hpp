///|/ Copyright (c) Prusa Research 2023 Tomáš Mészáros @tamasmeszaros
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef SLA_ARCHIVE_FORMAT_REGISTRY_HPP
#define SLA_ARCHIVE_FORMAT_REGISTRY_HPP

#include "SLAArchiveWriter.hpp"
#include "SLAArchiveReader.hpp"
#include <cstring>

namespace Slic3r {

// Factory function that returns an implementation of SLAArchiveWriter given
// a printer configuration.
using ArchiveWriterFactory = std::function<
    std::unique_ptr<SLAArchiveWriter>(const SLAPrinterConfig &)
>;

// Factory function that returns an implementation of SLAArchiveReader
using ArchiveReaderFactory = std::function<
    std::unique_ptr<SLAArchiveReader>(const std::string       &fname,
                                      SLAImportQuality         quality,
                                      const ProgrFn & progr)
>;

struct ArchiveEntry {
    OutputFormat format;

    // Main ID for the format, for internal unique identification
    const char *id;

    // Generic description (usable in GUI) about an archive format. Should only
    // be marked for localization (macro L).
    const char *desc = "";

    // Main extension of the format.
    const char *ext = "zip";

    ArchiveWriterFactory wrfactoryfn;
    ArchiveReaderFactory rdfactoryfn;

    // Secondary, alias extensions
    std::vector<const char *> ext_aliases;

    explicit ArchiveEntry(OutputFormat format, const char *formatid) : format(format), id{formatid} {}

    ArchiveEntry(OutputFormat format,
                 const char *formatid,
                 const char *description,
                 const char *extension,
                 std::initializer_list<const char *> extaliases,
                 const ArchiveWriterFactory &wrfn,
                 const ArchiveReaderFactory &rdfn)
        : format(format)
        , id{formatid}
        , desc{description}
        , ext{extension}
        , wrfactoryfn{wrfn}
        , rdfactoryfn{rdfn}
        , ext_aliases{extaliases}
    {}

    bool operator <(const ArchiveEntry &other) const
    {
        return std::strcmp(id, other.id) < 0;
    }
};

std::vector<std::string> get_extensions(const ArchiveEntry &entry);

const std::map<OutputFormat, ArchiveEntry>& registered_sla_archives();

const ArchiveEntry * get_archive_entry(OutputFormat formatid);
const char * get_default_extension(OutputFormat formatid);
ArchiveWriterFactory get_writer_factory(OutputFormat formatid);
ArchiveReaderFactory get_reader_factory(OutputFormat formatid);

} // namespace Slic3r

#endif // ARCHIVEREGISTRY_HPP
