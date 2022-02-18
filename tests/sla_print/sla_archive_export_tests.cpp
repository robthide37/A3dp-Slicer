#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include "libslic3r/SLAPrint.hpp"
#include "libslic3r/Format/SLAArchive.hpp"

#include <boost/filesystem.hpp>

using namespace Slic3r;

TEST_CASE("Archive export test", "[sla_archives]") {
    constexpr const char *PNAME = "20mm_cube";

    for (auto &archname : SLAArchive::registered_archives()) {
        INFO(std::string("Testing archive type: ") + archname);
        SLAPrint print;
        SLAFullPrintConfig fullcfg;

        auto m = Model::read_from_file(TEST_DATA_DIR PATH_SEPARATOR + std::string(PNAME) + ".obj", nullptr);

        fullcfg.set("sla_archive_format", archname);
        fullcfg.set("supports_enable", false);
        fullcfg.set("pad_enable", false);

        DynamicPrintConfig cfg;
        cfg.apply(fullcfg);

        print.set_status_callback([](const PrintBase::SlicingStatus&) {});
        print.apply(m, cfg);
        print.process();

        ThumbnailsList thumbnails;
        auto outputfname = std::string("output.") + SLAArchive::get_extension(archname);

        print.export_print(outputfname, thumbnails, PNAME);

        // Not much can be checked about the archives...
        REQUIRE(boost::filesystem::exists(outputfname));
    }
}
