//Copyright (c) 2021 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include "BeadingStrategyFactory.hpp"

#include "LimitedBeadingStrategy.hpp"
#include "CenterDeviationBeadingStrategy.hpp"
#include "WideningBeadingStrategy.hpp"
#include "DistributedBeadingStrategy.hpp"
#include "RedistributeBeadingStrategy.hpp"
#include "OuterWallInsetBeadingStrategy.hpp"

#include <limits>
#include <boost/log/trivial.hpp>

namespace Slic3r::Arachne
{

coord_t getWeightedAverage(const coord_t preferred_bead_width_outer, const coord_t preferred_bead_width_inner, const coord_t max_bead_count)
{
    if(max_bead_count > preferred_bead_width_outer - preferred_bead_width_inner)
    {
        //The difference between outer and inner bead width would be spread out across so many lines that rounding errors would destroy the difference.
        //Also catches the case of max_bead_count being "infinite" (max integer).
        return (preferred_bead_width_outer + preferred_bead_width_inner) / 2;
    }
    if (max_bead_count > 2)
    {
        return ((preferred_bead_width_outer * 2) + preferred_bead_width_inner * (max_bead_count - 2)) / max_bead_count;
    }
    if (max_bead_count <= 0)
    {
        return preferred_bead_width_inner;
    }
    return preferred_bead_width_outer;
}

BeadingStrategyPtr BeadingStrategyFactory::makeStrategy
(
    const BeadingStrategyType type,
    const coord_t preferred_bead_width_outer,
    const coord_t preferred_bead_width_inner,
    const coord_t preferred_transition_length,
    const float transitioning_angle,
    const bool print_thin_walls,
    const coord_t min_bead_width,
    const coord_t min_feature_size,
    const double wall_split_middle_threshold,
    const double wall_add_middle_threshold,
    const coord_t max_bead_count,
    const coord_t outer_wall_offset,
    const int inward_distributed_center_wall_count,
    const double minimum_variable_line_width
)
{
    using std::make_unique;
    using std::move;
    const coord_t bar_preferred_wall_width = getWeightedAverage(preferred_bead_width_outer, preferred_bead_width_inner, max_bead_count);
    BeadingStrategyPtr ret;
    switch (type)
    {
        case BeadingStrategyType::Center:
            ret = make_unique<CenterDeviationBeadingStrategy>(bar_preferred_wall_width, transitioning_angle, wall_split_middle_threshold, wall_add_middle_threshold);
            break;
        case BeadingStrategyType::Distributed:
            ret = make_unique<DistributedBeadingStrategy>(bar_preferred_wall_width, preferred_transition_length, transitioning_angle, wall_split_middle_threshold, wall_add_middle_threshold, std::numeric_limits<int>::max());
            break;
        case BeadingStrategyType::InwardDistributed:
            ret = make_unique<DistributedBeadingStrategy>(bar_preferred_wall_width, preferred_transition_length, transitioning_angle, wall_split_middle_threshold, wall_add_middle_threshold, inward_distributed_center_wall_count);
            break;
        default:
            BOOST_LOG_TRIVIAL(error) << "Cannot make strategy!";
            return nullptr;
    }
    
    if(print_thin_walls)
    {
        BOOST_LOG_TRIVIAL(debug) << "Applying the Widening Beading meta-strategy with minimum input width " << min_feature_size << " and minimum output width " << min_bead_width << ".";
        ret = make_unique<WideningBeadingStrategy>(move(ret), min_feature_size, min_bead_width);
    }
    if (max_bead_count > 0)
    {
        BOOST_LOG_TRIVIAL(debug) << "Applying the Redistribute meta-strategy with outer-wall width = " << preferred_bead_width_outer << ", inner-wall width = " << preferred_bead_width_inner;
        ret = make_unique<RedistributeBeadingStrategy>(preferred_bead_width_outer, preferred_bead_width_inner, minimum_variable_line_width, move(ret));
        //Apply the LimitedBeadingStrategy last, since that adds a 0-width marker wall which other beading strategies shouldn't touch.
        BOOST_LOG_TRIVIAL(debug) << "Applying the Limited Beading meta-strategy with maximum bead count = " << max_bead_count << ".";
        ret = make_unique<LimitedBeadingStrategy>(max_bead_count, move(ret));
    }
    
    if (outer_wall_offset > 0)
    {
        BOOST_LOG_TRIVIAL(debug) << "Applying the OuterWallOffset meta-strategy with offset = " << outer_wall_offset << ".";
        ret = make_unique<OuterWallInsetBeadingStrategy>(outer_wall_offset, move(ret));
    }
    return ret;
}
} // namespace Slic3r::Arachne
