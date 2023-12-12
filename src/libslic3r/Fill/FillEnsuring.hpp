///|/ Copyright (c) Prusa Research 2023 Pavel Mikuš @Godrak, Lukáš Hejl @hejllukas
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_FillEnsuring_hpp_
#define slic3r_FillEnsuring_hpp_

#include "FillBase.hpp"
#include "FillRectilinear.hpp"

namespace Slic3r {

class FillEnsuring : public Fill
{
protected:
    ThickPolylines make_fill_polylines(const Surface *surface, const FillParams &params, bool stop_vibrations, bool fill_gaps, bool connect_extrusions) const;

public:
    Fill *clone() const override { return new FillEnsuring(*this); }
    ~FillEnsuring() override = default;
    Polylines fill_surface(const Surface *surface, const FillParams &params) const override
    {
        throw new Slic3r::RuntimeError("error, trying to fillsurface a FillEnsuring");
        return {};
    };
    ThickPolylines fill_surface_arachne(const Surface *surface, const FillParams &params) const override
    {
        return this->make_fill_polylines(surface, params, true, true, true);
    };

protected:
    void fill_surface_single_arachne(const Surface &surface, const FillParams &params, ThickPolylines &thick_polylines_out);

    bool no_sort() const override { return true; }

    // PrintRegionConfig is used for computing overlap between boundary contour and inner Rectilinear infill.
    const PrintRegionConfig *print_region_config = nullptr;

    friend class Layer;
};

} // namespace Slic3r

#endif // slic3r_FillEnsuring_hpp_
