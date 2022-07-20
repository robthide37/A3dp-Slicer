use Test::More;
use strict;
use warnings;

BEGIN {
    use FindBin;
    use lib "$FindBin::Bin/../lib";
    use local::lib "$FindBin::Bin/../local-lib";
}

use List::Util qw(first);
use Slic3r;
use Slic3r::Surface ':types';
use Slic3r::Test;

plan tests => 8;

{
    my $test = sub {
        my ($config) = @_;
        
        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
        ok my $gcode = Slic3r::Test::gcode($print), "infill_every_layers does not crash";
        
        my $tool = undef;
        my %layers = ();        # layer_z => 1
        my %layer_infill = ();  # layer_z => has_infill
        Slic3r::GCode::Reader->new->parse($gcode, sub {
            my ($self, $cmd, $args, $info) = @_;
            
            if ($cmd =~ /^T(\d+)/) {
                $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0 && $tool != $config->support_material_extruder-1) {
                $layer_infill{$self->Z} //= 0;
                if ($tool == $config->infill_extruder-1) {
                    $layer_infill{$self->Z} = 1;
                }
            }
            # Previously, all G-code commands had a fixed number of decimal points with means with redundant zeros after decimal points.
            # We changed this behavior and got rid of these redundant padding zeros, which caused this test to fail
            # because the position in Z-axis is compared as a string, and previously, G-code contained the following two commands:
            # "G1 Z5 F5000 ; lift nozzle"
            # "G1 Z5.000 F7800.000"
            # That has a different Z-axis position from the view of string comparisons of floating-point numbers.
            # To correct the computation of the number of printed layers, even in the case of string comparisons of floating-point numbers,
            # we filtered out the G-code command with the commend 'lift nozzle'.
            $layers{$args->{Z}} = 1 if $cmd eq 'G1' && $info->{dist_Z} && index($info->{comment}, 'lift nozzle') == -1;
        });
        
        my $layers_with_perimeters = scalar(keys %layer_infill);
        my $layers_with_infill = grep $_ > 0,  values %layer_infill;
        is scalar(keys %layers), $layers_with_perimeters+$config->raft_layers, 'expected number of layers';
        
        if ($config->raft_layers == 0) {
            # first infill layer printed directly on print bed is not combined, so we don't consider it.
            $layers_with_infill--;
            $layers_with_perimeters--;
        }
        
        # we expect that infill is generated for half the number of combined layers
        # plus for each single layer that was not combined (remainder)
        is $layers_with_infill,
            int($layers_with_perimeters/$config->infill_every_layers) + ($layers_with_perimeters % $config->infill_every_layers),
            'infill is only present in correct number of layers';
    };
    
    my $config = Slic3r::Config::new_from_defaults;
    $config->set('layer_height', 0.2);
    $config->set('first_layer_height', 0.2);
    $config->set('nozzle_diameter', [0.5,0.5,0.5,0.5]);
    $config->set('infill_every_layers', 2);
    $config->set('perimeter_extruder', 1);
    $config->set('infill_extruder', 2);
    $config->set('wipe_into_infill', 0);
    $config->set('support_material_extruder', 3);
    $config->set('support_material_interface_extruder', 3);
    $config->set('top_solid_layers', 0);
    $config->set('bottom_solid_layers', 0);
    $test->($config);
    
    $config->set('skirts', 0);  # prevent usage of perimeter_extruder in raft layers
    $config->set('raft_layers', 5);
    $test->($config);
}

{
    my $config = Slic3r::Config::new_from_defaults;
    $config->set('layer_height', 0.2);
    $config->set('first_layer_height', 0.2);
    $config->set('nozzle_diameter', [0.5]);
    $config->set('infill_every_layers', 2);
    
    my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
    $print->process;
    
    ok defined(first { @{$_->get_region(0)->fill_surfaces->filter_by_type(S_TYPE_INTERNALVOID)} > 0 }
        @{$print->print->get_object(0)->layers}),
        'infill combination produces internal void surfaces';
    
    # we disable combination after infill has been generated
    $config->set('infill_every_layers', 1);
    $print->apply($print->print->model->clone, $config);
    $print->process;
    
    ok !(defined first { @{$_->get_region(0)->fill_surfaces} == 0 }
        @{$print->print->get_object(0)->layers}),
            'infill combination is idempotent';
}

__END__
