package Slic3r::Line;
use strict;
use warnings;

# a line is a two-points line
use parent 'Slic3r::Polyline';

sub grow {
    my $self = shift;
    return Slic3r::Polyline->new(@$self)->grow(@_);
}

1;
