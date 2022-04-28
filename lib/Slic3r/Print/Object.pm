package Slic3r::Print::Object;
# extends c++ class Slic3r::PrintObject (Print.xsp)
use strict;
use warnings;

use List::Util qw(min max sum first);
use Slic3r::Flow ':roles';
use Slic3r::Print::State ':steps';
use Slic3r::Surface ':types';

sub layers {
    my $self = shift;
    return [ map $self->get_layer($_), 0..($self->layer_count - 1) ];
}

sub support_layers {
    my $self = shift;
    return [ map $self->get_support_layer($_), 0..($self->support_layer_count - 1) ];
}

1;
