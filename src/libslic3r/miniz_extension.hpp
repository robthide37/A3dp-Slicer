///|/ Copyright (c) Prusa Research 2019 - 2020 Tomáš Mészáros @tamasmeszaros
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef MINIZ_EXTENSION_HPP
#define MINIZ_EXTENSION_HPP

#include <string>
#include <miniz.h>

namespace Slic3r {

bool open_zip_reader(mz_zip_archive *zip, const std::string &fname_utf8);
bool open_zip_writer(mz_zip_archive *zip, const std::string &fname_utf8);
bool close_zip_reader(mz_zip_archive *zip);
bool close_zip_writer(mz_zip_archive *zip);

/***
* RAII wrapper for open_zip_reader & close_zip_reader
*/
class ZipReader
{
    /*const*/ bool m_success;
public:
    mz_zip_archive archive;

    ZipReader(const std::string &fname_utf8) {
        mz_zip_zero_struct(&archive);
        m_success = open_zip_reader(&archive, fname_utf8);
    }
    ~ZipReader() {
        if (m_success) {
            close_zip_reader(&archive);
        }
    }

    bool success() { return m_success; }
};

class MZ_Archive {
public:
    mz_zip_archive arch;
    
    MZ_Archive();
    
    static std::string get_errorstr(mz_zip_error mz_err);
    
    std::string get_errorstr() const
    {
        return get_errorstr(arch.m_last_error) + "!";
    }

    bool is_alive() const
    {
        return arch.m_zip_mode != MZ_ZIP_MODE_WRITING_HAS_BEEN_FINALIZED;
    }
};

} // namespace Slic3r

#endif // MINIZ_EXTENSION_HPP
