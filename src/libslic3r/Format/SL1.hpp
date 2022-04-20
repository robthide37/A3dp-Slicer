#ifndef ARCHIVETRAITS_HPP
#define ARCHIVETRAITS_HPP

#include <string>

#include "SLAArchiveWriter.hpp"

#include "libslic3r/Zipper.hpp"
#include "libslic3r/PrintConfig.hpp"

namespace Slic3r {

class SL1Archive: public SLAArchiveWriter {
    SLAPrinterConfig m_cfg;
    
protected:
    std::unique_ptr<sla::RasterBase> create_raster() const override;
    sla::RasterEncoder get_encoder() const override;

    SLAPrinterConfig & cfg() { return m_cfg; }
    const SLAPrinterConfig & cfg() const { return m_cfg; }

    void export_print(Zipper &,
                      const SLAPrint       &print,
                      const ThumbnailsList &thumbnails,
                      const std::string    &projectname);

public:

    SL1Archive() = default;
    explicit SL1Archive(const SLAPrinterConfig &cfg): m_cfg(cfg) {}
    explicit SL1Archive(SLAPrinterConfig &&cfg): m_cfg(std::move(cfg)) {}

    void export_print(const std::string     fname,
                      const SLAPrint       &print,
                      const ThumbnailsList &thumbnails,
                      const std::string    &projectname = "") override;
};
    


} // namespace Slic3r::sla

#endif // ARCHIVETRAITS_HPP
