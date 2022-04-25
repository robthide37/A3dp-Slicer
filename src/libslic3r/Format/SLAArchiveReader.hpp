#ifndef SLAARCHIVEREADER_HPP
#define SLAARCHIVEREADER_HPP

#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ExPolygon.hpp"

struct indexed_triangle_set;

namespace Slic3r {

enum class SLAImportQuality { Accurate, Balanced, Fast };

class MissingProfileError : public RuntimeError { using RuntimeError::RuntimeError; };

using ProgrFn = std::function<bool(int)>;

class SLAArchiveReader {
public:

    virtual ~SLAArchiveReader() = default;

    virtual ConfigSubstitutions read(std::vector<ExPolygons> &slices,
                                     DynamicPrintConfig      &profile) = 0;

    virtual ConfigSubstitutions read(DynamicPrintConfig &profile) = 0;

    static std::unique_ptr<SLAArchiveReader> create(
        const std::string &fname,
        SLAImportQuality   quality = SLAImportQuality::Balanced,
        const ProgrFn     &progr   = [](int) { return false; });

    // Get the names of currently known archive reader implementations
    static const std::vector<const char *> & registered_archives();

    // Get the default file extensions belonging to an archive format
    static std::vector<const char *> get_extensions(const char *archtype);

    static const char * get_description(const char *archtype);
};

class ReaderUnimplementedError : public RuntimeError { using RuntimeError::RuntimeError; };

ConfigSubstitutions import_sla_archive(const std::string  &zipfname,
                                       DynamicPrintConfig &out);

ConfigSubstitutions import_sla_archive(
    const std::string       &zipfname,
    indexed_triangle_set    &out,
    DynamicPrintConfig      &profile,
    SLAImportQuality         quality = SLAImportQuality::Balanced,
    const ProgrFn &progr = [](int) { return true; });

} // namespace Slic3r

#endif // SLAARCHIVEREADER_HPP
