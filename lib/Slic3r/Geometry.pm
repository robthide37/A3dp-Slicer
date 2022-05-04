package Slic3r::Geometry;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);

# Exported by this module. The last section starting with convex_hull is exported by Geometry.xsp
our @EXPORT_OK = qw(
    PI epsilon 

    scale
    unscale
    scaled_epsilon

    X Y Z
    convex_hull
    chained_path_from
    deg2rad
    rad2deg
);

use constant PI => 4 * atan2(1, 1);
use constant A => 0;
use constant B => 1;
use constant X1 => 0;
use constant Y1 => 1;
use constant X2 => 2;
use constant Y2 => 3;

sub epsilon () { 1E-4 }
sub scaled_epsilon () { epsilon / &Slic3r::SCALING_FACTOR }

sub scale   ($) { $_[0] / &Slic3r::SCALING_FACTOR }
sub unscale ($) { $_[0] * &Slic3r::SCALING_FACTOR }

# 2D
sub bounding_box {
    my ($points) = @_;
    
    my @x = map $_->x, @$points;
    my @y = map $_->y, @$points;    #,,
    my @bb = (undef, undef, undef, undef);
    for (0..$#x) {
        $bb[X1] = $x[$_] if !defined $bb[X1] || $x[$_] < $bb[X1];
        $bb[X2] = $x[$_] if !defined $bb[X2] || $x[$_] > $bb[X2];
        $bb[Y1] = $y[$_] if !defined $bb[Y1] || $y[$_] < $bb[Y1];
        $bb[Y2] = $y[$_] if !defined $bb[Y2] || $y[$_] > $bb[Y2];
    }
    
    return @bb[X1,Y1,X2,Y2];
}

1;
