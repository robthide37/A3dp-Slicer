///|/ Copyright (c) Prusa Research 2021 - 2023 Lukáš Matěna @lukasmatena, Vojtěch Bubník @bubnikv, Lukáš Hejl @hejllukas
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_FillLightning_hpp_
#define slic3r_FillLightning_hpp_

#include "FillBase.hpp"

namespace Slic3r {

class PrintObject;

namespace FillLightning {

class Generator;
// To keep the definition of Octree opaque, we have to define a custom deleter.
struct GeneratorDeleter { void operator()(Generator *p); };
using  GeneratorPtr = std::unique_ptr<Generator, GeneratorDeleter>;

GeneratorPtr build_generator(const PrintObject &print_object, const coordf_t fill_density, const std::function<void()> &throw_on_cancel_callback);

class Filler : public Slic3r::Fill
{
public:
    ~Filler() override = default;

    Generator   *generator { nullptr };
protected:
    Fill* clone() const override { return new Filler(*this); }

    void _fill_surface_single(const FillParams              &params,
                              unsigned int                   thickness_layers,
                              const std::pair<float, Point> &direction,
                              ExPolygon                      expolygon,
                              Polylines &polylines_out) const override;

    // Let the G-code export reoder the infill lines.
	bool no_sort() const override { return false; }
};

} // namespace FillAdaptive
} // namespace Slic3r

#endif // slic3r_FillLightning_hpp_
