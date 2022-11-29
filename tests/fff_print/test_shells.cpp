#include <catch2/catch.hpp>

#include "libslic3r/GCodeReader.hpp"

#include "test_data.hpp" // get access to init_print, etc

using namespace Slic3r::Test;
using namespace Slic3r;

SCENARIO("Shells", "[Shells]") {
    GIVEN("20mm box") {
        auto test = [](const DynamicPrintConfig &config){
            std::string gcode = Slic3r::Test::slice({ Slic3r::Test::TestMesh::cube_20x20x20 }, config);
            
            std::vector<coord_t> zs;
            std::set<coord_t> layers_with_solid_infill;
            std::set<coord_t> layers_with_bridge_infill;
            const double solid_infill_speed = config.opt_float("solid_infill_speed") * 60;
            const double bridge_speed       = config.opt_float("bridge_speed") * 60;
            GCodeReader parser;
            parser.parse_buffer(gcode,
                [&zs, &layers_with_solid_infill, &layers_with_bridge_infill, solid_infill_speed, bridge_speed]
                    (Slic3r::GCodeReader &self, const Slic3r::GCodeReader::GCodeLine &line)
            {
                double z = line.new_Z(self);
                REQUIRE(z >= 0);
                if (z > 0) {
                    coord_t scaled_z = scaled<float>(z);
                    zs.emplace_back(scaled_z);
                    if (line.extruding(self) && line.dist_XY(self) > 0) {
                        double f = line.new_F(self);
                        if (std::abs(f - solid_infill_speed) < EPSILON)
                            layers_with_solid_infill.insert(scaled_z);
                        if (std::abs(f - bridge_speed) < EPSILON)
                            layers_with_bridge_infill.insert(scaled_z);
                    }
                }
            });
            sort_remove_duplicates(zs);

            auto has_solid_infill  = [&layers_with_solid_infill](coord_t z) { return layers_with_solid_infill.find(z) != layers_with_solid_infill.end(); };
            auto has_bridge_infill = [&layers_with_bridge_infill](coord_t z) { return layers_with_bridge_infill.find(z) != layers_with_bridge_infill.end(); };
            auto has_shells        = [&has_solid_infill, &has_bridge_infill, &zs](int layer_idx) { coord_t z = zs[layer_idx]; return has_solid_infill(z) || has_bridge_infill(z); };
            const int bottom_solid_layers = config.opt_int("bottom_solid_layers");
            const int top_solid_layers    = config.opt_int("top_solid_layers");
            THEN("correct number of bottom solid layers") {
                for (int i = 0; i < bottom_solid_layers; ++ i)
                    REQUIRE(has_shells(i));
                for (int i = bottom_solid_layers; i < int(zs.size() / 2); ++ i)
                    REQUIRE(! has_shells(i));
            }
            THEN("correct number of top solid layers") {
                for (int i = 0; i < top_solid_layers; ++ i)
                    REQUIRE(has_shells(int(zs.size()) - i - 1));
                for (int i = top_solid_layers; i < int(zs.size() / 2); ++ i)
                    REQUIRE(! has_shells(int(zs.size()) - i - 1));
            }
            if (top_solid_layers > 0) {
                THEN("solid infill speed is used on solid infill") {
                    for (int i = 0; i < top_solid_layers - 1; ++ i) {
                        auto z = zs[int(zs.size()) - i - 1];
                        REQUIRE(has_solid_infill(z));
                        REQUIRE(! has_bridge_infill(z));
                    }
                }
                THEN("bridge used in first solid layer over sparse infill") {
                    auto z = zs[int(zs.size()) - top_solid_layers];
                    REQUIRE(! has_solid_infill(z));
                    REQUIRE(has_bridge_infill(z));
                }
            }
        };

        auto config = Slic3r::DynamicPrintConfig::full_print_config_with({
                { "skirts",                 0 },
                { "perimeters",             0 },
                { "solid_infill_speed",     99 },
                { "top_solid_infill_speed", 99 },
                { "bridge_speed",           72 },
                { "first_layer_speed",      "100%" },
                { "cooling",                "0" }
            });

        WHEN("three top and bottom layers") {
            // proper number of shells is applied
            config.set_deserialize_strict({
                { "top_solid_layers",       3 },
                { "bottom_solid_layers",    3 }
            });
            test(config);
        }

        WHEN("zero top and bottom layers") {
            // no shells are applied when both top and bottom are set to zero
            config.set_deserialize_strict({
                { "top_solid_layers",       0 },
                { "bottom_solid_layers",    0 }
            });
            test(config);
        }

        WHEN("three top and bottom layers, zero infill") {
            // proper number of shells is applied even when fill density is none
            config.set_deserialize_strict({
                { "perimeters",             1 },
                { "top_solid_layers",       3 },
                { "bottom_solid_layers",    3 }
            });
            test(config);
        }
    }
}
