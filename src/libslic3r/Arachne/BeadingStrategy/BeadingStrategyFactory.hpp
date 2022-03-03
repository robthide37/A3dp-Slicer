// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef BEADING_STRATEGY_FACTORY_H
#define BEADING_STRATEGY_FACTORY_H

#include "BeadingStrategy.hpp"
#include "../../Point.hpp"

namespace Slic3r::Arachne
{

class BeadingStrategyFactory
{
public:
    static BeadingStrategyPtr makeStrategy
    (
        const BeadingStrategyType type,
        const coord_t preferred_bead_width_outer = scaled<coord_t>(0.0005),
        const coord_t preferred_bead_width_inner = scaled<coord_t>(0.0005),
        const coord_t preferred_transition_length = scaled<coord_t>(0.0004),
        const float transitioning_angle = M_PI / 4.0,
        const bool print_thin_walls = false,
        const coord_t min_bead_width = 0,
        const coord_t min_feature_size = 0,
        const double wall_split_middle_threshold = 0.5,
        const double wall_add_middle_threshold = 0.5,
        const coord_t max_bead_count = 0,
        const coord_t outer_wall_offset = 0,
        const int inward_distributed_center_wall_count = 2,
        const double minimum_variable_line_width = 0.5
    );
};

} // namespace Slic3r::Arachne
#endif // BEADING_STRATEGY_FACTORY_H
