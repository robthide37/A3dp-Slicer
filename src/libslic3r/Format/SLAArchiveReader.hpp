#ifndef SLAARCHIVEREADER_HPP
#define SLAARCHIVEREADER_HPP

#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ExPolygon.hpp"

struct indexed_triangle_set;

namespace Slic3r {

enum class SLAImportQuality { Accurate, Balanced, Fast };

class MissingProfileError : public RuntimeError { using RuntimeError::RuntimeError; };

class SLAArchiveReader {
public:

    virtual ~SLAArchiveReader() = default;

    virtual ConfigSubstitutions read(std::vector<ExPolygons> &slices,
                                     DynamicPrintConfig      &profile) = 0;

    virtual ConfigSubstitutions read(DynamicPrintConfig &profile) = 0;

    static std::unique_ptr<SLAArchiveReader> create(
        const std::string       &fname,
        SLAImportQuality         quality = SLAImportQuality::Balanced,
        std::function<bool(int)> progr = [](int){ return false; });
};

class ReaderUnimplementedError : public RuntimeError { using RuntimeError::RuntimeError; };

ConfigSubstitutions import_sla_archive(const std::string  &zipfname,
                                       DynamicPrintConfig &out);

ConfigSubstitutions import_sla_archive(
    const std::string       &zipfname,
    indexed_triangle_set    &out,
    DynamicPrintConfig      &profile,
    SLAImportQuality         quality = SLAImportQuality::Balanced,
    std::function<bool(int)> progr = [](int) { return true; });

} // namespace Slic3r

#endif // SLAARCHIVEREADER_HPP
