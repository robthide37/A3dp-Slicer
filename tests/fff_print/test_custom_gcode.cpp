#include <catch2/catch.hpp>

#include <numeric>
#include <sstream>

#include "libslic3r/Config.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/libslic3r.h"

#include "test_data.hpp"

using namespace Slic3r;

#if 0
SCENARIO("Output file format", "[CustomGCode]")
{
    WHEN("output_file_format set") {
        auto config = Slic3r::DynamicPrintConfig::full_print_config_with({
            { "output_filename_format",         "ts_[travel_speed]_lh_[layer_height].gcode" },
            { "start_gcode",                    "TRAVEL:[travel_speed] HEIGHT:[layer_height]\n" }
        });

        Print print;
        Model model;
        Test::init_print({ Test::TestMesh::cube_2x20x10 }, print, model, config);
    
        std::string output_file = print.output_filepath();
        THEN("print config options are replaced in output filename") {
            output_file.find(std::string("ts_") + )
        }
        my ($t, $h) = map $config->$_, qw(travel_speed layer_height);
        ok $output_file =~ /ts_${t}_/, '';
        ok $output_file =~ /lh_$h\./, 'region config options are replaced in output filename';
    
        std::string gcode = print.gcode(print);
        THEN("print config options are replaced in custom G-code") {
            ok $gcode =~ /TRAVEL:$t/, '';    
        }
        THEN("region config options are replaced in custom G-code") {
            ok $gcode =~ /HEIGHT:$h/, '';
        }
    }
}

SCENARIO("Custom G-code", "[CustomGCode]")
{
    WHEN("start_gcode and layer_gcode set") {
        auto config = Slic3r::DynamicPrintConfig::full_print_config_with({
            { "start_gcode", "_MY_CUSTOM_START_GCODE_" },  // to avoid dealing with the nozzle lift in start G-code
            { "layer_gcode", "_MY_CUSTOM_LAYER_GCODE_" }
        });
        GCodeReader parser;
        bool        last_move_was_z_change = false;
        int         num_layer_changes_not_applied = 0;
        parser.parse_buffer(Slic3r::Test::slice({ Test::TestMesh::cube_2x20x10 }, config), 
            [&last_move_was_z_change, &num_layer_changes_not_applied](Slic3r::GCodeReader &self, const Slic3r::GCodeReader::GCodeLine &line)
        {
            if (line.extruding(self)) {
                if (! was_extruding)
                    seam_points.emplace_back(self.xy_scaled());
                was_extruding = true;
            } else if (! line.cmd_is("M73")) {
                // skips remaining time lines (M73)
                was_extruding = false;
            }
            if (last_move_was_z_change != line.cmd_is("_MY_CUSTOM_LAYER_GCODE_"))
                ++ num_layer_changes_not_applied;
            last_move_was_z_change = line.dist_Z(self) > 0;
        });
        THEN("custom layer G-code is applied after Z move and before other moves");
    };

{
    my $config = Slic3r::Config->new;
    $config->set('nozzle_diameter', [0.6,0.6,0.6,0.6]);
    $config->set('extruder', 2);
    $config->set('first_layer_temperature', [200,205]);
    
    {
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        my $gcode = Slic3r::Test::gcode($print);
        ok $gcode =~ /M104 S205 T1/, 'temperature set correctly for non-zero yet single extruder';
        ok $gcode !~ /M104 S\d+ T0/, 'unused extruder correctly ignored';
    }
    
    $config->set('infill_extruder', 1);
    {
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        my $gcode = Slic3r::Test::gcode($print);
        ok $gcode =~ /M104 S200 T0/, 'temperature set correctly for first extruder';
        ok $gcode =~ /M104 S205 T1/, 'temperature set correctly for second extruder';
    }
    
    my @start_gcode = (qq!
;__temp0:[first_layer_temperature_0]__
;__temp1:[first_layer_temperature_1]__
;__temp2:[first_layer_temperature_2]__
    !, qq!
;__temp0:{first_layer_temperature[0]}__
;__temp1:{first_layer_temperature[1]}__
;__temp2:{first_layer_temperature[2]}__
    !);
    my @syntax_description = (' (legacy syntax)', ' (new syntax)');
    for my $i (0, 1) {
        $config->set('start_gcode', $start_gcode[$i]);
        {
            my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
            my $gcode = Slic3r::Test::gcode($print);
            # we use the [infill_extruder] placeholder to make sure this test doesn't
            # catch a false positive caused by the unparsed start G-code option itself
            # being embedded in the G-code
            ok $gcode =~ /temp0:200/, 'temperature placeholder for first extruder correctly populated' . $syntax_description[$i];
            ok $gcode =~ /temp1:205/, 'temperature placeholder for second extruder correctly populated' . $syntax_description[$i];
            ok $gcode =~ /temp2:200/, 'temperature placeholder for unused extruder populated with first value' . $syntax_description[$i];
        }
    }

    $config->set('start_gcode', qq!
;substitution:{if infill_extruder==1}extruder1
         {elsif infill_extruder==2}extruder2
         {else}extruder3{endif}
    !);
    {
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        my $gcode = Slic3r::Test::gcode($print);
        ok $gcode =~ /substitution:extruder1/, 'if / else / endif - first block returned';
    }
}

{
    my $config = Slic3r::Config::new_from_defaults;
    $config->set('before_layer_gcode', ';BEFORE [layer_num]');
    $config->set('layer_gcode', ';CHANGE [layer_num]');
    $config->set('support_material', 1);
    $config->set('layer_height', 0.2);
    my $print = Slic3r::Test::init_print('overhang', config => $config);
    my $gcode = Slic3r::Test::gcode($print);
    
    my @before = ();
    my @change = ();
    foreach my $line (split /\R+/, $gcode) {
        if ($line =~ /;BEFORE (\d+)/) {
            push @before, $1;
        } elsif ($line =~ /;CHANGE (\d+)/) {
            push @change, $1;
            fail 'inconsistent layer_num before and after layer change'
                if $1 != $before[-1];
        }
    }
    is_deeply \@before, \@change, 'layer_num is consistent before and after layer changes';
    ok !defined(first { $change[$_] != $change[$_-1]+1 } 1..$#change),
        'layer_num grows continously';  # i.e. no duplicates or regressions
}

{
    my $config = Slic3r::Config->new;
    $config->set('nozzle_diameter', [0.6,0.6,0.6,0.6,0.6]);
    $config->set('start_gcode', qq!
;substitution:{if infill_extruder==1}if block
         {elsif infill_extruder==2}elsif block 1
         {elsif infill_extruder==3}elsif block 2
         {elsif infill_extruder==4}elsif block 3
         {else}endif block{endif}
    !);
    my @returned = ('', 'if block', 'elsif block 1', 'elsif block 2', 'elsif block 3', 'endif block');
    for my $i (1,2,3,4,5) {
        $config->set('infill_extruder', $i);
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        my $gcode = Slic3r::Test::gcode($print);
        my $found_other = 0;
        for my $j (1,2,3,4,5) {
            next if $i == $j;
            $found_other = 1 if $gcode =~ /substitution:$returned[$j]/;
        }
        ok $gcode =~ /substitution:$returned[$i]/, 'if / else / endif - ' . $returned[$i] . ' returned';
        ok !$found_other, 'if / else / endif - only ' . $returned[$i] . ' returned';
    }
}

{
    my $config = Slic3r::Config->new;
    $config->set('nozzle_diameter', [0.6,0.6,0.6,0.6]);
    $config->set('start_gcode', 
        ';substitution:{if infill_extruder==1}{if perimeter_extruder==1}block11{else}block12{endif}' .
        '{elsif infill_extruder==2}{if perimeter_extruder==1}block21{else}block22{endif}' .
        '{else}{if perimeter_extruder==1}block31{else}block32{endif}{endif}:end');
    for my $i (1,2,3) {
        $config->set('infill_extruder', $i);
        for my $j (1,2) {
            $config->set('perimeter_extruder', $j);
            my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
            my $gcode = Slic3r::Test::gcode($print);
            ok $gcode =~ /substitution:block$i$j:end/, "two level if / else / endif - block$i$j returned";
        }
    }
}

{
    my $config = Slic3r::Config->new;
    $config->set('start_gcode', 
        ';substitution:{if notes=="MK2"}MK2{elsif notes=="MK3"}MK3{else}MK1{endif}:end');
    for my $printer_name ("MK2", "MK3", "MK1") {
        $config->set('notes', $printer_name);
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        my $gcode = Slic3r::Test::gcode($print);
        ok $gcode =~ /substitution:$printer_name:end/, "printer name $printer_name matched";
    }
}

{
    my $config = Slic3r::Config::new_from_defaults;
    $config->set('complete_objects', 1);
    $config->set('between_objects_gcode', '_MY_CUSTOM_GCODE_');
    my $print = Slic3r::Test::init_print('20mm_cube', config => $config, duplicate => 3);
    my $gcode = Slic3r::Test::gcode($print);
    is scalar(() = $gcode =~ /^_MY_CUSTOM_GCODE_/gm), 2, 'between_objects_gcode is applied correctly';
}

#endif