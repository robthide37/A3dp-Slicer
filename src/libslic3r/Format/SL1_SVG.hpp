#ifndef SL1_SVG_HPP
#define SL1_SVG_HPP

#include "SL1.hpp"

namespace Slic3r {

class SL1_SVGArchive: public SL1Archive {
protected:

    // Override the factory methods to produce svg instead of a real raster.
    std::unique_ptr<sla::RasterBase> create_raster() const override;
    sla::RasterEncoder get_encoder() const override;

public:

    using SL1Archive::SL1Archive;
};

} // namespace Slic3r

#endif // SL1_SVG_HPP
