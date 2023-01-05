#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Arachne/WallToolPaths.hpp"

#include "FillEnsuring.hpp"

#include <boost/log/trivial.hpp>

namespace Slic3r {

ThickPolylines FillEnsuring::fill_surface_arachne(const Surface *surface, const FillParams &params)
{
    // Perform offset.
    Slic3r::ExPolygons expp = this->overlap != 0. ? offset_ex(surface->expolygon, scaled<float>(this->overlap)) : ExPolygons{surface->expolygon};
    // Create the infills for each of the regions.
    ThickPolylines thick_polylines_out;
    for (ExPolygon &ex_poly : expp)
        fill_surface_single_arachne(Surface(*surface, std::move(ex_poly)), params, thick_polylines_out);

    return thick_polylines_out;
}

void FillEnsuring::fill_surface_single_arachne(const Surface &surface, const FillParams &params, ThickPolylines &thick_polylines_out)
{
    assert(params.use_arachne);
    assert(this->print_config != nullptr && this->print_object_config != nullptr && this->print_region_config != nullptr);

    coord_t                scaled_spacing = scaled<coord_t>(this->spacing);
    Polygons               polygons       = to_polygons(surface.expolygon);
    Arachne::WallToolPaths wall_tool_paths(polygons, scaled_spacing, scaled_spacing, 1, 0, params.layer_height, *this->print_object_config, *this->print_config);
    if (std::vector<Arachne::VariableWidthLines> loop = wall_tool_paths.getToolPaths(); !loop.empty()) {
        assert(loop.size() == 1);

        size_t firts_poly_idx = thick_polylines_out.size();
        Point  last_pos(0, 0);
        for (const Arachne::ExtrusionLine &extrusion : loop.front()) {
            if (extrusion.empty())
                continue;

            ThickPolyline thick_polyline = Arachne::to_thick_polyline(extrusion);
            if (thick_polyline.length() == 0.)
                //FIXME this should not happen.
                continue;
            assert(thick_polyline.size() > 1);
            assert(thick_polyline.length() > 0.);
            //assert(thick_polyline.points.size() == thick_polyline.width.size());
            if (extrusion.is_closed)
                thick_polyline.start_at_index(nearest_point_index(thick_polyline.points, last_pos));

            assert(thick_polyline.size() > 1);
            //assert(thick_polyline.points.size() == thick_polyline.width.size());
            thick_polylines_out.emplace_back(std::move(thick_polyline));
            last_pos = thick_polylines_out.back().last_point();
        }

        // clip the paths to prevent the extruder from getting exactly on the first point of the loop
        // Keep valid paths only.
        size_t j = firts_poly_idx;
        for (size_t i = firts_poly_idx; i < thick_polylines_out.size(); ++i) {
            assert(thick_polylines_out[i].size() > 1);
            assert(thick_polylines_out[i].length() > 0.);
            //assert(thick_polylines_out[i].points.size() == thick_polylines_out[i].width.size());
            thick_polylines_out[i].clip_end(this->loop_clipping);
            assert(thick_polylines_out[i].size() > 1);
            if (thick_polylines_out[i].is_valid()) {
                if (j < i)
                    thick_polylines_out[j] = std::move(thick_polylines_out[i]);
                ++j;
            }
        }
        if (j < thick_polylines_out.size())
            thick_polylines_out.erase(thick_polylines_out.begin() + int(j), thick_polylines_out.end());
    }

    // Remaining infill area will be filled with classic Rectilinear infill.
    ExPolygons infill_contour = union_ex(wall_tool_paths.getInnerContour());
    if (offset_ex(infill_contour, -float(scaled_spacing / 2.)).empty())
        infill_contour.clear(); // Infill region is too small, so let's filter it out.

    Polygons pp;
    for (ExPolygon &ex : infill_contour)
        ex.simplify_p(scaled<double>(params.resolution), &pp);

    // Collapse too narrow infill areas and append them to thick_polylines_out.
    const auto min_perimeter_infill_spacing = coord_t(scaled_spacing * (1. - INSET_OVERLAP_TOLERANCE));
    const auto infill_overlap = coord_t(scale_(this->print_region_config->get_abs_value("infill_overlap", this->spacing)));
    for (ExPolygon &ex_poly : offset2_ex(union_ex(pp), float(-min_perimeter_infill_spacing / 2.), float(infill_overlap + min_perimeter_infill_spacing / 2.))) {
        Polylines polylines;
        if (Surface new_surface(surface, std::move(ex_poly)); !fill_surface_by_lines(&new_surface, params, 0.f, 0.f, polylines))
            BOOST_LOG_TRIVIAL(error) << "FillEnsuring::fill_surface() failed to fill a region.";
        append(thick_polylines_out, to_thick_polylines(std::move(polylines), scaled<coord_t>(this->spacing)));
    }
}

} // namespace Slic3r
