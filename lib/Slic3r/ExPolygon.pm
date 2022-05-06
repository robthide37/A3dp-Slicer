package Slic3r::ExPolygon;
use strict;
use warnings;

# an ExPolygon is a polygon with holes

sub bounding_box {
    my $self = shift;
    return $self->contour->bounding_box;
}

1;
