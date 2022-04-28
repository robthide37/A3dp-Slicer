package Slic3r::Geometry::Clipper;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
	offset 
	offset2_ex
    diff_ex diff union_ex 
    union);

1;
