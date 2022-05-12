#include "../Print.hpp"

#include "FillLightning.hpp"
#include "Lightning/Generator.hpp"
#include "../Surface.hpp"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace Slic3r::FillLightning {

Polylines Filler::fill_surface(const Surface *surface, const FillParams &params)
{
    const Layer &layer = generator->getTreesForLayer(this->layer_id);
    return layer.convertToLines(to_polygons(surface->expolygon), generator->infilll_extrusion_width());
}

void GeneratorDeleter::operator()(Generator *p) {
    delete p;
}

GeneratorPtr build_generator(const PrintObject &print_object, const std::function<void()> &throw_on_cancel_callback)
{
    return GeneratorPtr(new Generator(print_object, throw_on_cancel_callback));
}

} // namespace Slic3r::FillAdaptive
