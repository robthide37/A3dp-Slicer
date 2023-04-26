#ifndef _SLIC3R_FORMAT_PWMX_HPP_
#define _SLIC3R_FORMAT_PWMX_HPP_

#include <string>

#include "SLAArchiveWriter.hpp"

#include "libslic3r/PrintConfig.hpp"

#define ANYCUBIC_SLA_FORMAT_VERSION_1     1
#define ANYCUBIC_SLA_FORMAT_VERSION_515 515
#define ANYCUBIC_SLA_FORMAT_VERSION_516 516
#define ANYCUBIC_SLA_FORMAT_VERSION_517 517

#define ANYCUBIC_SLA_FORMAT_VERSIONED(FILEFORMAT, NAME, VERSION) \
    { FILEFORMAT, { FILEFORMAT, [] (const auto &cfg) { return std::make_unique<AnycubicSLAArchive>(cfg, VERSION); } } }

#define ANYCUBIC_SLA_FORMAT(FILEFORMAT, NAME) \
    ANYCUBIC_SLA_FORMAT_VERSIONED(FILEFORMAT, NAME, ANYCUBIC_SLA_FORMAT_VERSION_1)

/**
    // Supports only ANYCUBIC_SLA_VERSION_1
    ANYCUBIC_SLA_FORMAT_VERSIONED("pws", "Photon / Photon S", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pw0", "Photon Zero", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwx", "Photon X", ANYCUBIC_SLA_VERSION_1),

    // Supports ANYCUBIC_SLA_VERSION_1 and ANYCUBIC_SLA_VERSION_515
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwmo", "Photon Mono", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwms", "Photon Mono SE", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("dlp", "Photon Ultra", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwmx", "Photon Mono X", ANYCUBIC_SLA_VERSION_1),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pmsq", "Photon Mono SQ", ANYCUBIC_SLA_VERSION_1),

    // Supports ANYCUBIC_SLA_VERSION_515 and ANYCUBIC_SLA_VERSION_516
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwma", "Photon Mono 4K", ANYCUBIC_SLA_VERSION_515),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pm3",  "Photon M3", ANYCUBIC_SLA_VERSION_515),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pm3m", "Photon M3 Max", ANYCUBIC_SLA_VERSION_515),

    // Supports NYCUBIC_SLA_VERSION_515 and ANYCUBIC_SLA_VERSION_516 and ANYCUBIC_SLA_VERSION_517
    ANYCUBIC_SLA_FORMAT_VERSIONED("pwmb", "Photon Mono X 6K / Photon M3 Plus", ANYCUBIC_SLA_VERSION_515),
    ANYCUBIC_SLA_FORMAT_VERSIONED("dl2p", "Photon Photon D2", ANYCUBIC_SLA_VERSION_515),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pmx2", "Photon Mono X2", ANYCUBIC_SLA_VERSION_515),
    ANYCUBIC_SLA_FORMAT_VERSIONED("pm3r", "Photon M3 Premium", ANYCUBIC_SLA_VERSION_515),
*/

namespace Slic3r {

class AnycubicSLAArchive: public SLAArchiveWriter {
    SLAPrinterConfig m_cfg;
    uint16_t m_version;

protected:
    std::unique_ptr<sla::RasterBase> create_raster() const override;
    sla::RasterEncoder get_encoder() const override;

    SLAPrinterConfig & cfg() { return m_cfg; }
    const SLAPrinterConfig & cfg() const { return m_cfg; }

public:
    
    AnycubicSLAArchive() = default;
    explicit AnycubicSLAArchive(const SLAPrinterConfig &cfg):
        m_cfg(cfg), m_version(ANYCUBIC_SLA_FORMAT_VERSION_1) {}
    explicit AnycubicSLAArchive(SLAPrinterConfig &&cfg):
        m_cfg(std::move(cfg)), m_version(ANYCUBIC_SLA_FORMAT_VERSION_1) {}

    explicit AnycubicSLAArchive(const SLAPrinterConfig &cfg, uint16_t version):
        m_cfg(cfg), m_version(version) {}
    explicit AnycubicSLAArchive(SLAPrinterConfig &&cfg, uint16_t version):
        m_cfg(std::move(cfg)), m_version(version) {}

    void export_print(const std::string     fname,
                      const SLAPrint       &print,
                      const ThumbnailsList &thumbnails,
                      const std::string    &projectname = "") override;
};


} // namespace Slic3r::sla

#endif // _SLIC3R_FORMAT_PWMX_HPP_
