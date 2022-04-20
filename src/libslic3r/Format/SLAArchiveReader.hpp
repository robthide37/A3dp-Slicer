#ifndef SLAARCHIVEREADER_HPP
#define SLAARCHIVEREADER_HPP

#include "libslic3r/PrintConfig.hpp"

struct indexed_triangle_set;

namespace Slic3r {

ConfigSubstitutions import_sla_archive(const std::string &zipfname, DynamicPrintConfig &out);

ConfigSubstitutions import_sla_archive(
    const std::string &      zipfname,
    Vec2i                    windowsize,
    indexed_triangle_set &   out,
    DynamicPrintConfig &     profile,
    std::function<bool(int)> progr = [](int) { return true; });

inline ConfigSubstitutions import_sla_archive(
    const std::string &      zipfname,
    Vec2i                    windowsize,
    indexed_triangle_set &   out,
    std::function<bool(int)> progr = [](int) { return true; })
{
    DynamicPrintConfig profile;
    return import_sla_archive(zipfname, windowsize, out, profile, progr);
}

class MissingProfileError : public RuntimeError { using RuntimeError::RuntimeError; };

} // namespace Slic3r

#endif // SLAARCHIVEREADER_HPP
