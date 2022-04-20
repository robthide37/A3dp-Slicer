#include "SLAArchiveReader.hpp"
#include "SL1.hpp"
#include "SL1_SVG.hpp"

#include "libslic3r/SlicesToTriangleMesh.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>

#include <array>

namespace Slic3r {

std::unique_ptr<SLAArchiveReader> SLAArchiveReader::create(
    const std::string       &fname,
    SLAImportQuality         quality,
    std::function<bool(int)> progr)
{
    std::string ext = boost::filesystem::path(fname).extension().string();
    boost::algorithm::to_lower(ext);

    std::unique_ptr<SLAArchiveReader> ret;

    const char *SL1_ext[] = {
        SLAArchiveWriter::get_extension("SL1"),
        "sl1s",
        // ...
    };

    const char *SL2_ext[] = {
        SLAArchiveWriter::get_extension("SL2"),
        "sl1_svg",
        // ...
    };

    if (!ext.empty()) {
        if (ext.front() == '.')
            ext.erase(ext.begin());

        auto extcmp = [&ext](const auto &e) { return e == ext; };

        if (std::any_of(std::begin(SL1_ext), std::end(SL1_ext), extcmp)) {
            ret = std::make_unique<SL1Reader>(fname, quality, progr);
        } else if (std::any_of(std::begin(SL2_ext), std::end(SL2_ext), extcmp)) {
            ret = std::make_unique<SL1_SVGReader>(fname, quality, progr);
        }
    }

    return ret;
}

struct SliceParams { double layerh = 0., initial_layerh = 0.; };

static SliceParams get_slice_params(const DynamicPrintConfig &cfg)
{
    auto *opt_layerh = cfg.option<ConfigOptionFloat>("layer_height");
    auto *opt_init_layerh = cfg.option<ConfigOptionFloat>("initial_layer_height");

    if (!opt_layerh || !opt_init_layerh)
        throw MissingProfileError("Invalid SL1 / SL1S file");

    return SliceParams{opt_layerh->getFloat(), opt_init_layerh->getFloat()};
}

ConfigSubstitutions import_sla_archive(const std::string       &zipfname,
                                       indexed_triangle_set    &out,
                                       DynamicPrintConfig      &profile,
                                       SLAImportQuality         quality,
                                       std::function<bool(int)> progr)
{
    ConfigSubstitutions ret;

    if (auto reader = SLAArchiveReader::create(zipfname, quality, progr)) {
        std::vector<ExPolygons> slices;
        ret = reader->read(slices, profile);

        SliceParams slicp = get_slice_params(profile);

        if (!slices.empty())
            out = slices_to_mesh(slices, 0, slicp.layerh, slicp.initial_layerh);

    } else {
        throw ReaderUnimplementedError("Reader unimplemented");
    }

    return ret;
}

ConfigSubstitutions import_sla_archive(const std::string  &zipfname,
                                       DynamicPrintConfig &out)
{
    ConfigSubstitutions ret;

    if (auto reader = SLAArchiveReader::create(zipfname)) {
        ret = reader->read(out);
    } else {
        throw ReaderUnimplementedError("Reader unimplemented");
    }

    return ret;
}

} // namespace Slic3r
