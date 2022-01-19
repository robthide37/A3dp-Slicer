#include <catch2/catch.hpp>

#include <memory>

#include "libslic3r/GCode/FindReplace.hpp"

using namespace Slic3r;

SCENARIO("Find/Replace with plain text", "[GCodeFindReplace]") {
    GIVEN("G-code") {
        const std::string gcode =
            "G1 Z0; home\n"
            "G1 Z1; move up\n"
            "G1 X0 Y1 Z1; perimeter\n"
            "G1 X13 Y32 Z1; infill\n"
            "G1 X13 Y32 Z1; wipe\n";
        WHEN("Replace \"move up\" with \"move down\", case sensitive") {
            GCodeFindReplace find_replace({ "move up", "move down", "" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move up\" with \"move down\", case insensitive") {
            GCodeFindReplace find_replace({ "move up", "move down", "i" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move UP\" with \"move down\", case insensitive") {
            GCodeFindReplace find_replace({ "move UP", "move down", "i" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move up\" with \"move down\", case sensitive") {
            GCodeFindReplace find_replace({ "move UP", "move down", "" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }

        // Whole word
        WHEN("Replace \"move up\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move up", "move down", "w" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move u\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move u", "move down", "w" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }
        WHEN("Replace \"ove up\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move u", "move down", "w" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }

        // Multi-line replace
        WHEN("Replace \"move up\\nG1 X0 \" with \"move down\\nG0 X1 \"") {
            GCodeFindReplace find_replace({ "move up\\nG1 X0 ", "move down\\nG0 X1 ", "" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G0 X1 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        // Multi-line replace, whole word.
        WHEN("Replace \"move up\\nG1 X0\" with \"move down\\nG0 X1\", whole word") {
            GCodeFindReplace find_replace({ "move up\\nG1 X0", "move down\\nG0 X1", "w" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G0 X1 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        // Multi-line replace, whole word, fails.
        WHEN("Replace \"move up\\nG1 X\" with \"move down\\nG0 X\", whole word") {
            GCodeFindReplace find_replace({ "move up\\nG1 X", "move down\\nG0 X", "w" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }
    }

    GIVEN("G-code with decimals") {
        const std::string gcode =
            "G1 Z0.123; home\n"
            "G1 Z1.21; move up\n"
            "G1 X0 Y.33 Z.431 E1.2; perimeter\n";
        WHEN("Regular expression NOT processed in non-regex mode") {
            GCodeFindReplace find_replace({ "( [XYZEF]-?)\\.([0-9]+)", "\\10.\\2", "" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }
    }
}

SCENARIO("Find/Replace with regexp", "[GCodeFindReplace]") {
    GIVEN("G-code") {
        const std::string gcode =
            "G1 Z0; home\n"
            "G1 Z1; move up\n"
            "G1 X0 Y1 Z1; perimeter\n"
            "G1 X13 Y32 Z1; infill\n"
            "G1 X13 Y32 Z1; wipe\n";
        WHEN("Replace \"move up\" with \"move down\", case sensitive") {
            GCodeFindReplace find_replace({ "move up", "move down", "r" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move up\" with \"move down\", case insensitive") {
            GCodeFindReplace find_replace({ "move up", "move down", "ri" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move UP\" with \"move down\", case insensitive") {
            GCodeFindReplace find_replace({ "move UP", "move down", "ri" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move up\" with \"move down\", case sensitive") {
            GCodeFindReplace find_replace({ "move UP", "move down", "r" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }

        // Whole word
        WHEN("Replace \"move up\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move up", "move down", "rw" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G1 X0 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        WHEN("Replace \"move u\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move u", "move down", "rw" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }
        WHEN("Replace \"ove up\" with \"move down\", whole word") {
            GCodeFindReplace find_replace({ "move u", "move down", "rw" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }

        // Multi-line replace
        WHEN("Replace \"move up\\nG1 X0 \" with \"move down\\nG0 X1 \"") {
            GCodeFindReplace find_replace({ "move up\\nG1 X0 ", "move down\\nG0 X1 ", "r" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G0 X1 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        // Multi-line replace, whole word.
        WHEN("Replace \"move up\\nG1 X0\" with \"move down\\nG0 X1\", whole word") {
            GCodeFindReplace find_replace({ "move up\\nG1 X0", "move down\\nG0 X1", "rw" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0; home\n"
                // substituted
                "G1 Z1; move down\n"
                "G0 X1 Y1 Z1; perimeter\n"
                "G1 X13 Y32 Z1; infill\n"
                "G1 X13 Y32 Z1; wipe\n");
        }
        // Multi-line replace, whole word, fails.
        WHEN("Replace \"move up\\nG1 X\" with \"move down\\nG0 X\", whole word") {
            GCodeFindReplace find_replace({ "move up\\nG1 X", "move down\\nG0 X", "rw" });
            REQUIRE(find_replace.process_layer(gcode) == gcode);
        }
    }

    GIVEN("G-code with decimals") {
        const std::string gcode =
            "G1 Z0.123; home\n"
            "G1 Z1.21; move up\n"
            "G1 X0 Y.33 Z.431 E1.2; perimeter\n";
        WHEN("Missing zeros before dot filled in") {
            GCodeFindReplace find_replace({ "( [XYZEF]-?)\\.([0-9]+)", "\\10.\\2", "r" });
            REQUIRE(find_replace.process_layer(gcode) ==
                "G1 Z0.123; home\n"
                "G1 Z1.21; move up\n"
                "G1 X0 Y0.33 Z0.431 E1.2; perimeter\n");
        }
    }
}
