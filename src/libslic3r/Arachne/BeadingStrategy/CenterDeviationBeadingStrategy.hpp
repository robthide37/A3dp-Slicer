// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.


#ifndef CENTER_DEVIATION_BEADING_STRATEGY_H
#define CENTER_DEVIATION_BEADING_STRATEGY_H

#include "BeadingStrategy.hpp"

namespace Slic3r::Arachne
{

/*!
 * This beading strategy makes the deviation in the thickness of the part
 * entirely compensated by the innermost wall.
 *
 * The outermost walls all use the ideal width, as far as possible.
 */
class CenterDeviationBeadingStrategy : public BeadingStrategy
{
  private:
    // For uneven numbers of lines: Minimum line width for which the middle line will be split into two lines.
    coord_t minimum_line_width_split;

    // For even numbers of lines: Minimum line width for which a new middle line will be added between the two innermost lines.
    coord_t minimum_line_width_add;

  public:
    CenterDeviationBeadingStrategy(coord_t pref_bead_width,
                                   double transitioning_angle,
                                   double wall_split_middle_threshold,
                                   double wall_add_middle_threshold);

    ~CenterDeviationBeadingStrategy() override{};
    Beading compute(coord_t thickness, coord_t bead_count) const override;
    coord_t getOptimalThickness(coord_t bead_count) const override;
    coord_t getTransitionThickness(coord_t lower_bead_count) const override;
    coord_t getOptimalBeadCount(coord_t thickness) const override;
};

} // namespace Slic3r::Arachne
#endif // CENTER_DEVIATION_BEADING_STRATEGY_H
