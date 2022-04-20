#ifndef SLAARCHIVEREADER_HPP
#define SLAARCHIVEREADER_HPP

#include "libslic3r/PrintConfig.hpp"

struct indexed_triangle_set;

namespace Slic3r {

ConfigSubstitutions import_sla_archive(const std::string &zipfname, DynamicPrintConfig &out);

enum class SLAImportQuality { Accurate, Balanced, Fast };

ConfigSubstitutions import_sla_archive(
    const std::string       &zipfname,
    indexed_triangle_set    &out,
    DynamicPrintConfig      &profile,
    SLAImportQuality         quality = SLAImportQuality::Balanced,
    std::function<bool(int)> progr = [](int) { return true; });

class MissingProfileError : public RuntimeError { using RuntimeError::RuntimeError; };

} // namespace Slic3r

#endif // SLAARCHIVEREADER_HPP
