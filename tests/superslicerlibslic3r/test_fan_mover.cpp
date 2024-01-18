
//#define CATCH_CONFIG_DISABLE

#include <catch_main.hpp>
#include "test_data.hpp"
#include <libslic3r/GCode/FanMover.hpp>
#include <libslic3r/ExtrusionEntity.hpp>

#include <boost/algorithm/string.hpp>

using namespace Slic3r;
using namespace Slic3r::Geometry;
using namespace Slic3r::Test;

// note: F600: 10mm/s

TEST_CASE("Simple gcode") {

    GCodeWriter writer;
    //what's used from the writer:
    writer.config.gcode_flavor.value = gcfMarlinFirmware;
    writer.config.gcode_comments.value = false;
    writer.config.fan_percentage.value = false; //0 -> 255
    //writer.tool[0] = nullptr;
    assert(writer.tool() == nullptr);
    assert(writer.get_tool(0) == nullptr);
    
    const std::string gcode = []() {
            std::string gcode;
            gcode += "M107\n";
            gcode += "G1 X0 Y0 F600\n"; // shouldn't move
            gcode += "G1 X20 Y0 E20\n"; // 2 sec
            gcode += "M106 S153\n";     // activate at 60%
            gcode += "G1 X40 Y0 E20\n"; // 2sec
            gcode += "M107\n";          // disable
            gcode += "G1 X60 Y0 E20\n"; // 2sec
            return gcode;
        }();
    
    SECTION("disable evrything")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        REQUIRE(gcode == processed_gcode);
    }

    SECTION("only kickstart 1sec")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   1      // fan_kickstart.value));
        );
        const std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        std::string wanted_gcode;
        wanted_gcode += "M107\n";
        wanted_gcode += "G1 X0 Y0 F600\n"; //shouldn't move
        wanted_gcode += "G1 X20 Y0 E20\n"; // 2 sec
        wanted_gcode += "M106 S255\n"; //activate with kickstart
        wanted_gcode += "G1 X26 Y0 E6\n"; // 1*0.6 = 0.6sec
        wanted_gcode += "M106 S153\n"; //go to good speed 60% of 1 second
        wanted_gcode += "G1 X40 Y0 E14\n"; // 2 - 1*0.6 = 1.4sec
        wanted_gcode += "M107\n"; //disable
        wanted_gcode += "G1 X60 Y0 E20\n"; // 2sec
        REQUIRE(wanted_gcode == processed_gcode);
    }
    
    SECTION("only speedup 1sec")
    {
        Slic3r::FanMover fan_mover(writer,
                                   1,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        std::string wanted_gcode;
        wanted_gcode += "M107\n";
        wanted_gcode += "G1 X0 Y0 F600\n"; //shouldn't move
        wanted_gcode += "G1 X10 Y0 E10\n"; // 1 sec
        wanted_gcode += "M106 S153\n"; //activate
        wanted_gcode += "G1 X20 Y0 E10\n"; // 1 sec
        wanted_gcode += "G1 X40 Y0 E20\n"; // 2sec
        wanted_gcode += "M107\n"; //disable
        wanted_gcode += "G1 X60 Y0 E20\n"; // 2sec
        REQUIRE(wanted_gcode == processed_gcode);
    }
    
    SECTION("speedup deactivated on not-overhangs")
    {
        Slic3r::FanMover fan_mover(writer,
                                   1,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   true, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string type_not_overhang;
        SECTION("erNone") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erNone) + "\n"; }
        SECTION("erPerimeter") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erPerimeter) + "\n"; }
        SECTION("erExternalPerimeter") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erExternalPerimeter) + "\n"; }
        SECTION("erInternalInfill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erInternalInfill) + "\n"; }
        SECTION("erSolidInfill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erSolidInfill) + "\n"; }
        SECTION("erTopSolidInfill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erTopSolidInfill) + "\n"; }
        SECTION("erBridgeInfill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erBridgeInfill) + "\n"; } //note: bridge infill isn't overhang
        SECTION("erInternalBridgeInfill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erInternalBridgeInfill) + "\n"; }
        SECTION("erThinWall") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erThinWall) + "\n"; }
        SECTION("erGapFill") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erGapFill) + "\n"; }
        SECTION("erSkirt") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erSkirt) + "\n"; }
        SECTION("erSupportMaterial") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erSupportMaterial) + "\n"; }
        SECTION("erSupportMaterialInterface") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erSupportMaterialInterface) + "\n"; }
        SECTION("erWipeTower") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erWipeTower) + "\n"; }
        SECTION("erMilling") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erMilling) + "\n"; }
        SECTION("erCustom") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erCustom) + "\n"; }
        SECTION("erMixed") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erMixed) + "\n"; }
        SECTION("erTravel") { type_not_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erTravel) + "\n"; }
        std::string processed_gcode = fan_mover.process_gcode(type_not_overhang + gcode, true);
        REQUIRE(type_not_overhang + gcode == processed_gcode);
    }
    SECTION("speedup activated on overhangs")
    {
        Slic3r::FanMover fan_mover(writer,
                                   1,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   true, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string type_overhang;
        SECTION("erOverhangPerimeter") { type_overhang = std::string(";TYPE:") + ExtrusionEntity::role_to_string(ExtrusionRole::erOverhangPerimeter) + "\n"; }
        std::string processed_gcode = fan_mover.process_gcode(type_overhang + gcode, true);
        std::string wanted_gcode;
        wanted_gcode += "M107\n";
        wanted_gcode += "G1 X0 Y0 F600\n"; //shouldn't move
        wanted_gcode += "G1 X10 Y0 E10\n"; // 1 sec
        wanted_gcode += "M106 S153\n"; //activate
        wanted_gcode += "G1 X20 Y0 E10\n"; // 1 sec
        wanted_gcode += "G1 X40 Y0 E20\n"; // 2sec
        wanted_gcode += "M107\n"; //disable
        wanted_gcode += "G1 X60 Y0 E20\n"; // 2sec
        REQUIRE(type_overhang + wanted_gcode == processed_gcode);
    }
    
    SECTION("erase double M107")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string m107 = "M107\n";
        std::string processed_gcode = fan_mover.process_gcode(m107+gcode+m107, false);
        processed_gcode += fan_mover.process_gcode(gcode+m107, true);
        REQUIRE(gcode+gcode.substr(5) == processed_gcode); // substr to remove first m107
    }
    SECTION("erase M107-like M106 -> M107")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   true,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string m106 = "M106 S0\n";
        std::string processed_gcode = fan_mover.process_gcode(gcode+m106, false);
        processed_gcode += fan_mover.process_gcode(gcode, true);
        REQUIRE(gcode+gcode.substr(5) == processed_gcode); // substr to remove first m107
    }
    //TODO?
    //SECTION("erase M106 -> M107")
    //{
    //    Slic3r::FanMover fan_mover(writer,
    //                               0,     // fan_speedup_time.value,
    //                               false, // with_D_option
    //                               true,  // use_relative_e_distances.value,
    //                               false, // fan_speedup_overhangs.value,
    //                               0      // fan_kickstart.value));
    //    );
    //    std::string m106 = "M106 S125\n";
    //    std::string processed_gcode = fan_mover.process_gcode(m106+gcode+m106, false);
    //    processed_gcode += fan_mover.process_gcode(gcode, true);
    //    REQUIRE(gcode+m106 == processed_gcode); // substr to remove first m107
    //}
}


TEST_CASE("Simple gcode (absoute)")
{
    GCodeWriter writer;
    // what's used from the writer:
    writer.config.gcode_flavor.value   = gcfMarlinFirmware;
    writer.config.gcode_comments.value = false;
    writer.config.fan_percentage.value = false; // 0 -> 255
    // writer.tool[0] = nullptr;
    assert(writer.tool() == nullptr);
    assert(writer.get_tool(0) == nullptr);

    const std::string gcode = []() {
        std::string gcode;
        gcode += "M107\n";
        gcode += "G1 X0 Y0 F600\n"; // shouldn't move
        gcode += "G1 X20 Y0 E20\n"; // 2 sec
        gcode += "M106 S153\n";     // activate at 60%
        gcode += "G1 X40 Y0 E40\n"; // 2sec
        gcode += "M107\n";          // disable
        gcode += "G1 X60 Y0 E60\n"; // 2sec
        return gcode;
    }();

    
    SECTION("disable evrything")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   false,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        REQUIRE(gcode == processed_gcode);
    }

    SECTION("only kickstart 1sec")
    {
        Slic3r::FanMover fan_mover(writer,
                                   0,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   false,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   1      // fan_kickstart.value));
        );
        const std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        std::string wanted_gcode;
        wanted_gcode += "M107\n";
        wanted_gcode += "G1 X0 Y0 F600\n"; //shouldn't move
        wanted_gcode += "G1 X20 Y0 E20\n"; // 2 sec
        wanted_gcode += "M106 S255\n"; //activate with kickstart
        wanted_gcode += "G1 X26 Y0 E26\n"; // 1*0.6 = 0.6sec
        wanted_gcode += "M106 S153\n"; //go to good speed 60% of 1 second
        wanted_gcode += "G1 X40 Y0 E40\n"; // 2 - 1*0.6 = 1.4sec
        wanted_gcode += "M107\n"; //disable
        wanted_gcode += "G1 X60 Y0 E60\n"; // 2sec
        REQUIRE(wanted_gcode == processed_gcode);
    }
    
    SECTION("only speedup 1sec")
    {
        Slic3r::FanMover fan_mover(writer,
                                   1,     // fan_speedup_time.value,
                                   false, // with_D_option
                                   false,  // use_relative_e_distances.value,
                                   false, // fan_speedup_overhangs.value,
                                   0      // fan_kickstart.value));
        );
        std::string processed_gcode = fan_mover.process_gcode(gcode, true);
        std::string wanted_gcode;
        wanted_gcode += "M107\n";
        wanted_gcode += "G1 X0 Y0 F600\n"; //shouldn't move
        wanted_gcode += "G1 X10 Y0 E10\n"; // 1 sec
        wanted_gcode += "M106 S153\n"; //activate
        wanted_gcode += "G1 X20 Y0 E20\n"; // 1 sec
        wanted_gcode += "G1 X40 Y0 E40\n"; // 2sec
        wanted_gcode += "M107\n"; //disable
        wanted_gcode += "G1 X60 Y0 E60\n"; // 2sec
        REQUIRE(wanted_gcode == processed_gcode);
    }
}

//TODO: complex gcode
