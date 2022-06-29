#include "PerimeterGenerator.hpp"
#include "ClipperUtils.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "ShortestPath.hpp"
#include "clipper/clipper_z.hpp"

#include "Arachne/WallToolPaths.hpp"
#include "Arachne/utils/ExtrusionLine.hpp"

#include <cmath>
#include <cassert>
#include <stack>

namespace Slic3r {

ExtrusionPaths thick_polyline_to_extrusion_paths(const ThickPolyline &thick_polyline, ExtrusionRole role, const Flow &flow, const float tolerance, const float merge_tolerance)
{
    ExtrusionPaths paths;
    ExtrusionPath path(role);
    ThickLines lines = thick_polyline.thicklines();
    
    for (int i = 0; i < (int)lines.size(); ++i) {
        const ThickLine& line = lines[i];
        
        const coordf_t line_len = line.length();
        if (line_len < SCALED_EPSILON) continue;
        
        double thickness_delta = fabs(line.a_width - line.b_width);
        if (thickness_delta > tolerance) {
            const auto segments = (unsigned int)ceil(thickness_delta / tolerance);
            const coordf_t seg_len = line_len / segments;
            Points pp;
            std::vector<coordf_t> width;
            {
                pp.push_back(line.a);
                width.push_back(line.a_width);
                for (size_t j = 1; j < segments; ++j) {
                    pp.push_back((line.a.cast<double>() + (line.b - line.a).cast<double>().normalized() * (j * seg_len)).cast<coord_t>());
                    
                    coordf_t w = line.a_width + (j*seg_len) * (line.b_width-line.a_width) / line_len;
                    width.push_back(w);
                    width.push_back(w);
                }
                pp.push_back(line.b);
                width.push_back(line.b_width);
                
                assert(pp.size() == segments + 1u);
                assert(width.size() == segments*2);
            }
            
            // delete this line and insert new ones
            lines.erase(lines.begin() + i);
            for (size_t j = 0; j < segments; ++j) {
                ThickLine new_line(pp[j], pp[j+1]);
                new_line.a_width = width[2*j];
                new_line.b_width = width[2*j+1];
                lines.insert(lines.begin() + i + j, new_line);
            }
            
            -- i;
            continue;
        }
        
        const double w = fmax(line.a_width, line.b_width);
        if (path.polyline.points.empty()) {
            path.polyline.append(line.a);
            path.polyline.append(line.b);
            // Convert from spacing to extrusion width based on the extrusion model
            // of a square extrusion ended with semi circles.
            Flow new_flow = (role == erOverhangPerimeter && flow.bridge()) ? flow : flow.with_width(unscale<float>(w) + flow.height() * float(1. - 0.25 * PI));
            #ifdef SLIC3R_DEBUG
            printf("  filling %f gap\n", flow.width);
            #endif
            path.mm3_per_mm  = new_flow.mm3_per_mm();
            path.width       = new_flow.width();
            path.height      = new_flow.height();
        } else {
            thickness_delta = fabs(scale_(flow.width()) - w);
            if (thickness_delta <= merge_tolerance) {
                // the width difference between this line and the current flow width is 
                // within the accepted tolerance
                path.polyline.append(line.b);
            } else {
                // we need to initialize a new line
                paths.emplace_back(std::move(path));
                path = ExtrusionPath(role);
                -- i;
            }
        }
    }
    if (path.polyline.is_valid())
        paths.emplace_back(std::move(path));
    return paths;
}

static void variable_width(const ThickPolylines& polylines, ExtrusionRole role, const Flow &flow, std::vector<ExtrusionEntity*> &out)
{
	// This value determines granularity of adaptive width, as G-code does not allow
	// variable extrusion within a single move; this value shall only affect the amount
	// of segments, and any pruning shall be performed before we apply this tolerance.
	const float tolerance = float(scale_(0.05));
	for (const ThickPolyline &p : polylines) {
		ExtrusionPaths paths = thick_polyline_to_extrusion_paths(p, role, flow, tolerance, tolerance);
		// Append paths to collection.
		if (! paths.empty()) {
			if (paths.front().first_point() == paths.back().last_point())
				out.emplace_back(new ExtrusionLoop(std::move(paths)));
			else {
				for (ExtrusionPath &path : paths)
					out.emplace_back(new ExtrusionPath(std::move(path)));
			}
		}
	}
}

// Hierarchy of perimeters.
class PerimeterGeneratorLoop {
public:
    // Polygon of this contour.
    Polygon                             polygon;
    // Is it a contour or a hole?
    // Contours are CCW oriented, holes are CW oriented.
    bool                                is_contour;
    // Depth in the hierarchy. External perimeter has depth = 0. An external perimeter could be both a contour and a hole.
    unsigned short                      depth;
    // Should this contur be fuzzyfied on path generation?
    bool                                fuzzify;
    // Children contour, may be both CCW and CW oriented (outer contours or holes).
    std::vector<PerimeterGeneratorLoop> children;
    
    PerimeterGeneratorLoop(const Polygon &polygon, unsigned short depth, bool is_contour, bool fuzzify) : 
        polygon(polygon), is_contour(is_contour), depth(depth), fuzzify(fuzzify) {}
    // External perimeter. It may be CCW or CW oriented (outer contour or hole contour).
    bool is_external() const { return this->depth == 0; }
    // An island, which may have holes, but it does not have another internal island.
    bool is_internal_contour() const;
};

// Thanks Cura developers for this function.
static void fuzzy_polygon(Polygon &poly, double fuzzy_skin_thickness, double fuzzy_skin_point_dist)
{
    const double min_dist_between_points = fuzzy_skin_point_dist * 3. / 4.; // hardcoded: the point distance may vary between 3/4 and 5/4 the supplied value
    const double range_random_point_dist = fuzzy_skin_point_dist / 2.;
    double dist_left_over = double(rand()) * (min_dist_between_points / 2) / double(RAND_MAX); // the distance to be traversed on the line before making the first new point
    Point* p0 = &poly.points.back();
    Points out;
    out.reserve(poly.points.size());
    for (Point &p1 : poly.points)
    { // 'a' is the (next) new point between p0 and p1
        Vec2d  p0p1      = (p1 - *p0).cast<double>();
        double p0p1_size = p0p1.norm();
        // so that p0p1_size - dist_last_point evaulates to dist_left_over - p0p1_size
        double dist_last_point = dist_left_over + p0p1_size * 2.;
        for (double p0pa_dist = dist_left_over; p0pa_dist < p0p1_size;
            p0pa_dist += min_dist_between_points + double(rand()) * range_random_point_dist / double(RAND_MAX))
        {
            double r = double(rand()) * (fuzzy_skin_thickness * 2.) / double(RAND_MAX) - fuzzy_skin_thickness;
            out.emplace_back(*p0 + (p0p1 * (p0pa_dist / p0p1_size) + perp(p0p1).cast<double>().normalized() * r).cast<coord_t>());
            dist_last_point = p0pa_dist;
        }
        dist_left_over = p0p1_size - dist_last_point;
        p0 = &p1;
    }
    while (out.size() < 3) {
        size_t point_idx = poly.size() - 2;
        out.emplace_back(poly[point_idx]);
        if (point_idx == 0)
            break;
        -- point_idx;
    }
    if (out.size() >= 3)
        poly.points = std::move(out);
}

// Thanks Cura developers for this function.
static void fuzzy_extrusion_line(Arachne::ExtrusionLine &ext_lines, double fuzzy_skin_thickness, double fuzzy_skin_point_dist)
{
    const double min_dist_between_points = fuzzy_skin_point_dist * 3. / 4.; // hardcoded: the point distance may vary between 3/4 and 5/4 the supplied value
    const double range_random_point_dist = fuzzy_skin_point_dist / 2.;
    double       dist_left_over          = double(rand()) * (min_dist_between_points / 2) / double(RAND_MAX); // the distance to be traversed on the line before making the first new point

    auto                                   *p0 = &ext_lines.front();
    std::vector<Arachne::ExtrusionJunction> out;
    out.reserve(ext_lines.size());
    for (auto &p1 : ext_lines) {
        if (p0->p == p1.p) { // Connect endpoints.
            out.emplace_back(p1.p, p1.w, p1.perimeter_index);
            continue;
        }

        // 'a' is the (next) new point between p0 and p1
        Vec2d  p0p1      = (p1.p - p0->p).cast<double>();
        double p0p1_size = p0p1.norm();
        // so that p0p1_size - dist_last_point evaulates to dist_left_over - p0p1_size
        double dist_last_point = dist_left_over + p0p1_size * 2.;
        for (double p0pa_dist = dist_left_over; p0pa_dist < p0p1_size; p0pa_dist += min_dist_between_points + double(rand()) * range_random_point_dist / double(RAND_MAX)) {
            double r = double(rand()) * (fuzzy_skin_thickness * 2.) / double(RAND_MAX) - fuzzy_skin_thickness;
            out.emplace_back(p0->p + (p0p1 * (p0pa_dist / p0p1_size) + perp(p0p1).cast<double>().normalized() * r).cast<coord_t>(), p1.w, p1.perimeter_index);
            dist_last_point = p0pa_dist;
        }
        dist_left_over = p0p1_size - dist_last_point;
        p0             = &p1;
    }

    while (out.size() < 3) {
        size_t point_idx = ext_lines.size() - 2;
        out.emplace_back(ext_lines[point_idx].p, ext_lines[point_idx].w, ext_lines[point_idx].perimeter_index);
        if (point_idx == 0)
            break;
        -- point_idx;
    }

    if (ext_lines.back().p == ext_lines.front().p) // Connect endpoints.
        out.back().p = out.front().p;

    if (out.size() >= 3)
        ext_lines.junctions = std::move(out);
}

using PerimeterGeneratorLoops = std::vector<PerimeterGeneratorLoop>;

static ExtrusionEntityCollection traverse_loops(const PerimeterGenerator &perimeter_generator, const PerimeterGeneratorLoops &loops, ThickPolylines &thin_walls)
{
    // loops is an arrayref of ::Loop objects
    // turn each one into an ExtrusionLoop object
    ExtrusionEntityCollection   coll;
    Polygon                     fuzzified;
    for (const PerimeterGeneratorLoop &loop : loops) {
        bool is_external = loop.is_external();
        
        ExtrusionRole role;
        ExtrusionLoopRole loop_role;
        role = is_external ? erExternalPerimeter : erPerimeter;
        if (loop.is_internal_contour()) {
            // Note that we set loop role to ContourInternalPerimeter
            // also when loop is both internal and external (i.e.
            // there's only one contour loop).
            loop_role = elrContourInternalPerimeter;
        } else {
            loop_role = elrDefault;
        }
        
        // detect overhanging/bridging perimeters
        ExtrusionPaths paths;
        const Polygon &polygon = loop.fuzzify ? fuzzified : loop.polygon;
        if (loop.fuzzify) {
            fuzzified = loop.polygon;
            fuzzy_polygon(fuzzified, scaled<float>(perimeter_generator.config->fuzzy_skin_thickness.value), scaled<float>(perimeter_generator.config->fuzzy_skin_point_dist.value));
        }
        if (perimeter_generator.config->overhangs && perimeter_generator.layer_id > perimeter_generator.object_config->raft_layers
            && ! ((perimeter_generator.object_config->support_material || perimeter_generator.object_config->support_material_enforce_layers > 0) && 
                  perimeter_generator.object_config->support_material_contact_distance.value == 0)) {
            // get non-overhang paths by intersecting this loop with the grown lower slices
            extrusion_paths_append(
                paths,
                intersection_pl({ polygon }, perimeter_generator.lower_slices_polygons()),
                role,
                is_external ? perimeter_generator.ext_mm3_per_mm()           : perimeter_generator.mm3_per_mm(),
                is_external ? perimeter_generator.ext_perimeter_flow.width() : perimeter_generator.perimeter_flow.width(),
                (float)perimeter_generator.layer_height);
            
            // get overhang paths by checking what parts of this loop fall 
            // outside the grown lower slices (thus where the distance between
            // the loop centerline and original lower slices is >= half nozzle diameter
            extrusion_paths_append(
                paths,
                diff_pl({ polygon }, perimeter_generator.lower_slices_polygons()),
                erOverhangPerimeter,
                perimeter_generator.mm3_per_mm_overhang(),
                perimeter_generator.overhang_flow.width(),
                perimeter_generator.overhang_flow.height());
            
            // Reapply the nearest point search for starting point.
            // We allow polyline reversal because Clipper may have randomly reversed polylines during clipping.
            chain_and_reorder_extrusion_paths(paths, &paths.front().first_point());
        } else {
            ExtrusionPath path(role);
            path.polyline   = polygon.split_at_first_point();
            path.mm3_per_mm = is_external ? perimeter_generator.ext_mm3_per_mm()           : perimeter_generator.mm3_per_mm();
            path.width      = is_external ? perimeter_generator.ext_perimeter_flow.width() : perimeter_generator.perimeter_flow.width();
            path.height     = (float)perimeter_generator.layer_height;
            paths.push_back(path);
        }
        
        coll.append(ExtrusionLoop(std::move(paths), loop_role));
    }
    
    // Append thin walls to the nearest-neighbor search (only for first iteration)
    if (! thin_walls.empty()) {
        variable_width(thin_walls, erExternalPerimeter, perimeter_generator.ext_perimeter_flow, coll.entities);
        thin_walls.clear();
    }
    
    // Traverse children and build the final collection.
	Point zero_point(0, 0);
	std::vector<std::pair<size_t, bool>> chain = chain_extrusion_entities(coll.entities, &zero_point);
    ExtrusionEntityCollection out;
    for (const std::pair<size_t, bool> &idx : chain) {
		assert(coll.entities[idx.first] != nullptr);
        if (idx.first >= loops.size()) {
            // This is a thin wall.
			out.entities.reserve(out.entities.size() + 1);
            out.entities.emplace_back(coll.entities[idx.first]);
			coll.entities[idx.first] = nullptr;
            if (idx.second)
				out.entities.back()->reverse();
        } else {
            const PerimeterGeneratorLoop &loop = loops[idx.first];
            assert(thin_walls.empty());
            ExtrusionEntityCollection children = traverse_loops(perimeter_generator, loop.children, thin_walls);
            out.entities.reserve(out.entities.size() + children.entities.size() + 1);
            ExtrusionLoop *eloop = static_cast<ExtrusionLoop*>(coll.entities[idx.first]);
            coll.entities[idx.first] = nullptr;
            if (loop.is_contour) {
                eloop->make_counter_clockwise();
                out.append(std::move(children.entities));
                out.entities.emplace_back(eloop);
            } else {
                eloop->make_clockwise();
                out.entities.emplace_back(eloop);
                out.append(std::move(children.entities));
            }
        }
    }
    return out;
}

static ClipperLib_Z::Paths clip_extrusion(const ClipperLib_Z::Path &subject, const ClipperLib_Z::Paths &clip, ClipperLib_Z::ClipType clipType)
{
    ClipperLib_Z::Clipper clipper;
    clipper.ZFillFunction([](const ClipperLib_Z::IntPoint &e1bot, const ClipperLib_Z::IntPoint &e1top, const ClipperLib_Z::IntPoint &e2bot,
                             const ClipperLib_Z::IntPoint &e2top, ClipperLib_Z::IntPoint &pt) {
        ClipperLib_Z::IntPoint start = e1bot;
        ClipperLib_Z::IntPoint end   = e1top;

        if (start.z() <= 0 && end.z() <= 0) {
            start = e2bot;
            end   = e2top;
        }

        assert(start.z() > 0 && end.z() > 0);

        // Interpolate extrusion line width.
        double length_sqr = (end - start).cast<double>().squaredNorm();
        double dist_sqr   = (pt - start).cast<double>().squaredNorm();
        double t          = std::sqrt(dist_sqr / length_sqr);

        pt.z() = start.z() + coord_t((end.z() - start.z()) * t);
    });

    clipper.AddPath(subject, ClipperLib_Z::ptSubject, false);
    clipper.AddPaths(clip, ClipperLib_Z::ptClip, true);

    ClipperLib_Z::PolyTree clipped_polytree;
    ClipperLib_Z::Paths    clipped_paths;
    clipper.Execute(clipType, clipped_polytree, ClipperLib_Z::pftNonZero, ClipperLib_Z::pftNonZero);
    ClipperLib_Z::PolyTreeToPaths(clipped_polytree, clipped_paths);

    return clipped_paths;
}

struct PerimeterGeneratorArachneExtrusion
{
    Arachne::ExtrusionLine *extrusion = nullptr;
    // Indicates if closed ExtrusionLine is a contour or a hole. Used it only when ExtrusionLine is a closed loop.
    bool is_contour = false;
    // Should this extrusion be fuzzyfied on path generation?
    bool fuzzify = false;
};

static ExtrusionEntityCollection traverse_extrusions(const PerimeterGenerator &perimeter_generator, std::vector<PerimeterGeneratorArachneExtrusion> &pg_extrusions)
{
    ExtrusionEntityCollection extrusion_coll;
    for (PerimeterGeneratorArachneExtrusion &pg_extrusion : pg_extrusions) {
        Arachne::ExtrusionLine *extrusion = pg_extrusion.extrusion;
        if (extrusion->empty())
            continue;

        const bool    is_external = extrusion->inset_idx == 0;
        ExtrusionRole role        = is_external ? erExternalPerimeter : erPerimeter;

        if (pg_extrusion.fuzzify)
            fuzzy_extrusion_line(*extrusion, scaled<float>(perimeter_generator.config->fuzzy_skin_thickness.value), scaled<float>(perimeter_generator.config->fuzzy_skin_point_dist.value));

        ExtrusionPaths paths;
        // detect overhanging/bridging perimeters
        if (perimeter_generator.config->overhangs && perimeter_generator.layer_id > perimeter_generator.object_config->raft_layers
            && ! ((perimeter_generator.object_config->support_material || perimeter_generator.object_config->support_material_enforce_layers > 0) &&
                 perimeter_generator.object_config->support_material_contact_distance.value == 0)) {

            ClipperLib_Z::Path extrusion_path;
            extrusion_path.reserve(extrusion->size());
            for (const Arachne::ExtrusionJunction &ej : extrusion->junctions)
                extrusion_path.emplace_back(ej.p.x(), ej.p.y(), ej.w);

            ClipperLib_Z::Paths lower_slices_paths;
            lower_slices_paths.reserve(perimeter_generator.lower_slices_polygons().size());
            for (const Polygon &poly : perimeter_generator.lower_slices_polygons()) {
                lower_slices_paths.emplace_back();
                ClipperLib_Z::Path &out = lower_slices_paths.back();
                out.reserve(poly.points.size());
                for (const Point &pt : poly.points)
                    out.emplace_back(pt.x(), pt.y(), 0);
            }

            // get non-overhang paths by intersecting this loop with the grown lower slices
            extrusion_paths_append(paths, clip_extrusion(extrusion_path, lower_slices_paths, ClipperLib_Z::ctIntersection), role,
                                   is_external ? perimeter_generator.ext_perimeter_flow : perimeter_generator.perimeter_flow);

            // get overhang paths by checking what parts of this loop fall
            // outside the grown lower slices (thus where the distance between
            // the loop centerline and original lower slices is >= half nozzle diameter
            extrusion_paths_append(paths, clip_extrusion(extrusion_path, lower_slices_paths, ClipperLib_Z::ctDifference), erOverhangPerimeter,
                                   perimeter_generator.overhang_flow);

            // Reapply the nearest point search for starting point.
            // We allow polyline reversal because Clipper may have randomly reversed polylines during clipping.
            // Arachne sometimes creates extrusion with zero-length (just two same endpoints);
            if (!paths.empty())
                chain_and_reorder_extrusion_paths(paths, &paths.front().first_point());
        } else {
            extrusion_paths_append(paths, *extrusion, role, is_external ? perimeter_generator.ext_perimeter_flow : perimeter_generator.perimeter_flow);
        }

        // Append paths to collection.
        if (!paths.empty()) {
            if (extrusion->is_closed) {
                ExtrusionLoop extrusion_loop(std::move(paths));
                // Restore the orientation of the extrusion loop.
                if (pg_extrusion.is_contour)
                    extrusion_loop.make_counter_clockwise();
                else
                    extrusion_loop.make_clockwise();

                extrusion_coll.append(std::move(extrusion_loop));
            } else
                for (ExtrusionPath &path : paths)
                    extrusion_coll.append(ExtrusionPath(std::move(path)));
        }
    }

    return extrusion_coll;
}

// Thanks, Cura developers, for implementing an algorithm for generating perimeters with variable width (Arachne) that is based on the paper
// "A framework for adaptive width control of dense contour-parallel toolpaths in fused deposition modeling"
void PerimeterGenerator::process_arachne()
{
    // other perimeters
    m_mm3_per_mm               	  = this->perimeter_flow.mm3_per_mm();
    coord_t perimeter_spacing     = this->perimeter_flow.scaled_spacing();

    // external perimeters
    m_ext_mm3_per_mm           	   = this->ext_perimeter_flow.mm3_per_mm();
    coord_t ext_perimeter_width    = this->ext_perimeter_flow.scaled_width();
    coord_t ext_perimeter_spacing  = this->ext_perimeter_flow.scaled_spacing();
    coord_t ext_perimeter_spacing2 = scaled<coord_t>(0.5f * (this->ext_perimeter_flow.spacing() + this->perimeter_flow.spacing()));

    // overhang perimeters
    m_mm3_per_mm_overhang         = this->overhang_flow.mm3_per_mm();

    // solid infill
    coord_t solid_infill_spacing  = this->solid_infill_flow.scaled_spacing();

    // prepare grown lower layer slices for overhang detection
    if (this->lower_slices != nullptr && this->config->overhangs) {
        // We consider overhang any part where the entire nozzle diameter is not supported by the
        // lower layer, so we take lower slices and offset them by half the nozzle diameter used
        // in the current layer
        double nozzle_diameter = this->print_config->nozzle_diameter.get_at(this->config->perimeter_extruder-1);
        m_lower_slices_polygons = offset(*this->lower_slices, float(scale_(+nozzle_diameter/2)));
    }

    // we need to process each island separately because we might have different
    // extra perimeters for each one
    for (const Surface &surface : this->slices->surfaces) {
        // detect how many perimeters must be generated for this island
        int        loop_number = this->config->perimeters + surface.extra_perimeters - 1; // 0-indexed loops
        ExPolygons last        = offset_ex(surface.expolygon.simplify_p(m_scaled_resolution), - float(ext_perimeter_width / 2. - ext_perimeter_spacing / 2.));
        Polygons   last_p      = to_polygons(last);

        Arachne::WallToolPaths wallToolPaths(last_p, ext_perimeter_spacing, perimeter_spacing, coord_t(loop_number + 1), 0, *this->object_config, *this->print_config);
        std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();
        loop_number = int(perimeters.size()) - 1;

        int start_perimeter = int(perimeters.size()) - 1;
        int end_perimeter   = -1;
        int direction       = -1;

        if (this->config->external_perimeters_first) {
            start_perimeter = 0;
            end_perimeter   = int(perimeters.size());
            direction       = 1;
        }

        std::vector<Arachne::ExtrusionLine *> all_extrusions;
        for (int perimeter_idx = start_perimeter; perimeter_idx != end_perimeter; perimeter_idx += direction) {
            if (perimeters[perimeter_idx].empty())
                continue;
            for (Arachne::ExtrusionLine &wall : perimeters[perimeter_idx])
                all_extrusions.emplace_back(&wall);
        }

        // Find topological order with constraints from extrusions_constrains.
        std::vector<size_t>              blocked(all_extrusions.size(), 0); // Value indicating how many extrusions it is blocking (preceding extrusions) an extrusion.
        std::vector<std::vector<size_t>> blocking(all_extrusions.size());   // Each extrusion contains a vector of extrusions that are blocked by this extrusion.
        std::unordered_map<const Arachne::ExtrusionLine *, size_t> map_extrusion_to_idx;
        for (size_t idx = 0; idx < all_extrusions.size(); idx++)
            map_extrusion_to_idx.emplace(all_extrusions[idx], idx);

        auto extrusions_constrains = Arachne::WallToolPaths::getRegionOrder(all_extrusions, this->config->external_perimeters_first);
        for (auto [before, after] : extrusions_constrains) {
            auto after_it = map_extrusion_to_idx.find(after);
            ++blocked[after_it->second];
            blocking[map_extrusion_to_idx.find(before)->second].emplace_back(after_it->second);
        }

        std::vector<bool> processed(all_extrusions.size(), false);          // Indicate that the extrusion was already processed.
        Point             current_position = all_extrusions.empty() ? Point::Zero() : all_extrusions.front()->junctions.front().p; // Some starting position.
        std::vector<PerimeterGeneratorArachneExtrusion> ordered_extrusions;         // To store our result in. At the end we'll std::swap.
        ordered_extrusions.reserve(all_extrusions.size());

        while (ordered_extrusions.size() < all_extrusions.size()) {
            size_t best_candidate    = 0;
            double best_distance_sqr = std::numeric_limits<double>::max();
            bool   is_best_closed    = false;

            std::vector<size_t> available_candidates;
            for (size_t candidate = 0; candidate < all_extrusions.size(); ++candidate) {
                if (processed[candidate] || blocked[candidate])
                    continue; // Not a valid candidate.
                available_candidates.push_back(candidate);
            }

            std::sort(available_candidates.begin(), available_candidates.end(), [&all_extrusions](const size_t a_idx, const size_t b_idx) -> bool {
                return all_extrusions[a_idx]->is_closed < all_extrusions[b_idx]->is_closed;
            });

            for (const size_t candidate_path_idx : available_candidates) {
                auto& path = all_extrusions[candidate_path_idx];

                if (path->junctions.empty()) { // No vertices in the path. Can't find the start position then or really plan it in. Put that at the end.
                    if (best_distance_sqr == std::numeric_limits<double>::max()) {
                        best_candidate = candidate_path_idx;
                        is_best_closed = path->is_closed;
                    }
                    continue;
                }

                const Point candidate_position = path->junctions.front().p;
                double      distance_sqr       = (current_position - candidate_position).cast<double>().norm();
                if (distance_sqr < best_distance_sqr) { // Closer than the best candidate so far.
                    if (path->is_closed || (!path->is_closed && best_distance_sqr != std::numeric_limits<double>::max()) || (!path->is_closed && !is_best_closed)) {
                        best_candidate    = candidate_path_idx;
                        best_distance_sqr = distance_sqr;
                        is_best_closed    = path->is_closed;
                    }
                }
            }

            auto &best_path = all_extrusions[best_candidate];
            ordered_extrusions.push_back({best_path, best_path->is_contour(), false});
            processed[best_candidate] = true;
            for (size_t unlocked_idx : blocking[best_candidate])
                blocked[unlocked_idx]--;

            if(!best_path->junctions.empty()) { //If all paths were empty, the best path is still empty. We don't upate the current position then.
                if(best_path->is_closed)
                    current_position = best_path->junctions[0].p; //We end where we started.
                else
                    current_position = best_path->junctions.back().p; //Pick the other end from where we started.
            }
        }

        if (this->layer_id > 0 && this->config->fuzzy_skin != FuzzySkinType::None) {
            std::vector<PerimeterGeneratorArachneExtrusion *> closed_loop_extrusions;
            for (PerimeterGeneratorArachneExtrusion &extrusion : ordered_extrusions)
                if (extrusion.extrusion->inset_idx == 0) {
                    if (extrusion.extrusion->is_closed && this->config->fuzzy_skin == FuzzySkinType::External) {
                        closed_loop_extrusions.emplace_back(&extrusion);
                    } else {
                        extrusion.fuzzify = true;
                    }
                }

            if (this->config->fuzzy_skin == FuzzySkinType::External) {
                ClipperLib_Z::Paths loops_paths;
                loops_paths.reserve(closed_loop_extrusions.size());
                for (const auto &cl_extrusion : closed_loop_extrusions) {
                    assert(cl_extrusion->extrusion->junctions.front() == cl_extrusion->extrusion->junctions.back());
                    size_t             loop_idx = &cl_extrusion - &closed_loop_extrusions.front();
                    ClipperLib_Z::Path loop_path;
                    loop_path.reserve(cl_extrusion->extrusion->junctions.size() - 1);
                    for (auto junction_it = cl_extrusion->extrusion->junctions.begin(); junction_it != std::prev(cl_extrusion->extrusion->junctions.end()); ++junction_it)
                        loop_path.emplace_back(junction_it->p.x(), junction_it->p.y(), loop_idx);
                    loops_paths.emplace_back(loop_path);
                }

                ClipperLib_Z::Clipper clipper;
                clipper.AddPaths(loops_paths, ClipperLib_Z::ptSubject, true);
                ClipperLib_Z::PolyTree loops_polytree;
                clipper.Execute(ClipperLib_Z::ctUnion, loops_polytree, ClipperLib_Z::pftEvenOdd, ClipperLib_Z::pftEvenOdd);

                for (const ClipperLib_Z::PolyNode *child_node : loops_polytree.Childs) {
                    // The whole contour must have the same index.
                    coord_t polygon_idx  = child_node->Contour.front().z();
                    bool    has_same_idx = std::all_of(child_node->Contour.begin(), child_node->Contour.end(),
                                                       [&polygon_idx](const ClipperLib_Z::IntPoint &point) -> bool { return polygon_idx == point.z(); });
                    if (has_same_idx)
                        closed_loop_extrusions[polygon_idx]->fuzzify = true;
                }
            }
        }

        if (ExtrusionEntityCollection extrusion_coll = traverse_extrusions(*this, ordered_extrusions); !extrusion_coll.empty())
            this->loops->append(extrusion_coll);

        ExPolygons    infill_contour = union_ex(wallToolPaths.getInnerContour());
        const coord_t spacing        = (perimeters.size() == 1) ? ext_perimeter_spacing2 : perimeter_spacing;
        if (offset_ex(infill_contour, -float(spacing / 2.)).empty())
            infill_contour.clear(); // Infill region is too small, so let's filter it out.

        // create one more offset to be used as boundary for fill
        // we offset by half the perimeter spacing (to get to the actual infill boundary)
        // and then we offset back and forth by half the infill spacing to only consider the
        // non-collapsing regions
        coord_t inset =
            (loop_number < 0) ? 0 :
            (loop_number == 0) ?
                                // one loop
                ext_perimeter_spacing:
                // two or more loops?
                perimeter_spacing;

        inset = coord_t(scale_(this->config->get_abs_value("infill_overlap", unscale<double>(inset))));
        Polygons pp;
        for (ExPolygon &ex : infill_contour)
            ex.simplify_p(m_scaled_resolution, &pp);
        // collapse too narrow infill areas
        const auto    min_perimeter_infill_spacing = coord_t(solid_infill_spacing * (1. - INSET_OVERLAP_TOLERANCE));
        // append infill areas to fill_surfaces
        this->fill_surfaces->append(
            offset2_ex(
                union_ex(pp),
                float(- min_perimeter_infill_spacing / 2.),
                float(inset + min_perimeter_infill_spacing / 2.)),
            stInternal);
    }
}

void PerimeterGenerator::process_classic()
{
    // other perimeters
    m_mm3_per_mm               		= this->perimeter_flow.mm3_per_mm();
    coord_t perimeter_width         = this->perimeter_flow.scaled_width();
    coord_t perimeter_spacing       = this->perimeter_flow.scaled_spacing();
    
    // external perimeters
    m_ext_mm3_per_mm           		= this->ext_perimeter_flow.mm3_per_mm();
    coord_t ext_perimeter_width     = this->ext_perimeter_flow.scaled_width();
    coord_t ext_perimeter_spacing   = this->ext_perimeter_flow.scaled_spacing();
    coord_t ext_perimeter_spacing2  = scaled<coord_t>(0.5f * (this->ext_perimeter_flow.spacing() + this->perimeter_flow.spacing()));
    
    // overhang perimeters
    m_mm3_per_mm_overhang      		= this->overhang_flow.mm3_per_mm();
    
    // solid infill
    coord_t solid_infill_spacing    = this->solid_infill_flow.scaled_spacing();
    
    // Calculate the minimum required spacing between two adjacent traces.
    // This should be equal to the nominal flow spacing but we experiment
    // with some tolerance in order to avoid triggering medial axis when
    // some squishing might work. Loops are still spaced by the entire
    // flow spacing; this only applies to collapsing parts.
    // For ext_min_spacing we use the ext_perimeter_spacing calculated for two adjacent
    // external loops (which is the correct way) instead of using ext_perimeter_spacing2
    // which is the spacing between external and internal, which is not correct
    // and would make the collapsing (thus the details resolution) dependent on 
    // internal flow which is unrelated.
    coord_t min_spacing         = coord_t(perimeter_spacing      * (1 - INSET_OVERLAP_TOLERANCE));
    coord_t ext_min_spacing     = coord_t(ext_perimeter_spacing  * (1 - INSET_OVERLAP_TOLERANCE));
    bool    has_gap_fill 		= this->config->gap_fill_enabled.value && this->config->gap_fill_speed.value > 0;

    // prepare grown lower layer slices for overhang detection
    if (this->lower_slices != NULL && this->config->overhangs) {
        // We consider overhang any part where the entire nozzle diameter is not supported by the
        // lower layer, so we take lower slices and offset them by half the nozzle diameter used 
        // in the current layer
        double nozzle_diameter = this->print_config->nozzle_diameter.get_at(this->config->perimeter_extruder-1);
        m_lower_slices_polygons = offset(*this->lower_slices, float(scale_(+nozzle_diameter/2)));
    }

    // we need to process each island separately because we might have different
    // extra perimeters for each one
    for (const Surface &surface : this->slices->surfaces) {
        // detect how many perimeters must be generated for this island
        int        loop_number = this->config->perimeters + surface.extra_perimeters - 1;  // 0-indexed loops
        ExPolygons last        = union_ex(surface.expolygon.simplify_p(m_scaled_resolution));
        ExPolygons gaps;
        if (loop_number >= 0) {
            // In case no perimeters are to be generated, loop_number will equal to -1.
            std::vector<PerimeterGeneratorLoops> contours(loop_number+1);    // depth => loops
            std::vector<PerimeterGeneratorLoops> holes(loop_number+1);       // depth => loops
            ThickPolylines thin_walls;
            // we loop one time more than needed in order to find gaps after the last perimeter was applied
            for (int i = 0;; ++ i) {  // outer loop is 0
                // Calculate next onion shell of perimeters.
                ExPolygons offsets;
                if (i == 0) {
                    // the minimum thickness of a single loop is:
                    // ext_width/2 + ext_spacing/2 + spacing/2 + width/2
                    offsets = this->config->thin_walls ? 
                        offset2_ex(
                            last,
                            - float(ext_perimeter_width / 2. + ext_min_spacing / 2. - 1),
                            + float(ext_min_spacing / 2. - 1)) :
                        offset_ex(last, - float(ext_perimeter_width / 2.));
                    // look for thin walls
                    if (this->config->thin_walls) {
                        // the following offset2 ensures almost nothing in @thin_walls is narrower than $min_width
                        // (actually, something larger than that still may exist due to mitering or other causes)
                        coord_t min_width = coord_t(scale_(this->ext_perimeter_flow.nozzle_diameter() / 3));
                        ExPolygons expp = opening_ex(
                            // medial axis requires non-overlapping geometry
                            diff_ex(last, offset(offsets, float(ext_perimeter_width / 2.) + ClipperSafetyOffset)),
                            float(min_width / 2.));
                        // the maximum thickness of our thin wall area is equal to the minimum thickness of a single loop
                        for (ExPolygon &ex : expp)
                            ex.medial_axis(ext_perimeter_width + ext_perimeter_spacing2, min_width, &thin_walls);
                    }
                    if (m_spiral_vase && offsets.size() > 1) {
                    	// Remove all but the largest area polygon.
                    	keep_largest_contour_only(offsets);
                    }
                } else {
                    //FIXME Is this offset correct if the line width of the inner perimeters differs
                    // from the line width of the infill?
                    coord_t distance = (i == 1) ? ext_perimeter_spacing2 : perimeter_spacing;
                    offsets = this->config->thin_walls ?
                        // This path will ensure, that the perimeters do not overfill, as in 
                        // prusa3d/Slic3r GH #32, but with the cost of rounding the perimeters
                        // excessively, creating gaps, which then need to be filled in by the not very 
                        // reliable gap fill algorithm.
                        // Also the offset2(perimeter, -x, x) may sometimes lead to a perimeter, which is larger than
                        // the original.
                        offset2_ex(last,
                                - float(distance + min_spacing / 2. - 1.),
                                float(min_spacing / 2. - 1.)) :
                        // If "detect thin walls" is not enabled, this paths will be entered, which 
                        // leads to overflows, as in prusa3d/Slic3r GH #32
                        offset_ex(last, - float(distance));
                    // look for gaps
                    if (has_gap_fill)
                        // not using safety offset here would "detect" very narrow gaps
                        // (but still long enough to escape the area threshold) that gap fill
                        // won't be able to fill but we'd still remove from infill area
                        append(gaps, diff_ex(
                            offset(last,    - float(0.5 * distance)),
                            offset(offsets,   float(0.5 * distance + 10))));  // safety offset
                }
                if (offsets.empty()) {
                    // Store the number of loops actually generated.
                    loop_number = i - 1;
                    // No region left to be filled in.
                    last.clear();
                    break;
                } else if (i > loop_number) {
                    // If i > loop_number, we were looking just for gaps.
                    break;
                }
                {
                    const bool fuzzify_contours = this->config->fuzzy_skin != FuzzySkinType::None && i == 0 && this->layer_id > 0;
                    const bool fuzzify_holes    = fuzzify_contours && this->config->fuzzy_skin == FuzzySkinType::All;
                    for (const ExPolygon &expolygon : offsets) {
    	                // Outer contour may overlap with an inner contour,
    	                // inner contour may overlap with another inner contour,
    	                // outer contour may overlap with itself.
    	                //FIXME evaluate the overlaps, annotate each point with an overlap depth,
                        // compensate for the depth of intersection.
                        contours[i].emplace_back(expolygon.contour, i, true, fuzzify_contours);

                        if (! expolygon.holes.empty()) {
                            holes[i].reserve(holes[i].size() + expolygon.holes.size());
                            for (const Polygon &hole : expolygon.holes)
                                holes[i].emplace_back(hole, i, false, fuzzify_holes);
                        }
                    }
                }
                last = std::move(offsets);
                if (i == loop_number && (! has_gap_fill || this->config->fill_density.value == 0)) {
                	// The last run of this loop is executed to collect gaps for gap fill.
                	// As the gap fill is either disabled or not 
                	break;
                }
            }

            // nest loops: holes first
            for (int d = 0; d <= loop_number; ++ d) {
                PerimeterGeneratorLoops &holes_d = holes[d];
                // loop through all holes having depth == d
                for (int i = 0; i < (int)holes_d.size(); ++ i) {
                    const PerimeterGeneratorLoop &loop = holes_d[i];
                    // find the hole loop that contains this one, if any
                    for (int t = d + 1; t <= loop_number; ++ t) {
                        for (int j = 0; j < (int)holes[t].size(); ++ j) {
                            PerimeterGeneratorLoop &candidate_parent = holes[t][j];
                            if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                                candidate_parent.children.push_back(loop);
                                holes_d.erase(holes_d.begin() + i);
                                -- i;
                                goto NEXT_LOOP;
                            }
                        }
                    }
                    // if no hole contains this hole, find the contour loop that contains it
                    for (int t = loop_number; t >= 0; -- t) {
                        for (int j = 0; j < (int)contours[t].size(); ++ j) {
                            PerimeterGeneratorLoop &candidate_parent = contours[t][j];
                            if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                                candidate_parent.children.push_back(loop);
                                holes_d.erase(holes_d.begin() + i);
                                -- i;
                                goto NEXT_LOOP;
                            }
                        }
                    }
                    NEXT_LOOP: ;
                }
            }
            // nest contour loops
            for (int d = loop_number; d >= 1; -- d) {
                PerimeterGeneratorLoops &contours_d = contours[d];
                // loop through all contours having depth == d
                for (int i = 0; i < (int)contours_d.size(); ++ i) {
                    const PerimeterGeneratorLoop &loop = contours_d[i];
                    // find the contour loop that contains it
                    for (int t = d - 1; t >= 0; -- t) {
                        for (size_t j = 0; j < contours[t].size(); ++ j) {
                            PerimeterGeneratorLoop &candidate_parent = contours[t][j];
                            if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                                candidate_parent.children.push_back(loop);
                                contours_d.erase(contours_d.begin() + i);
                                -- i;
                                goto NEXT_CONTOUR;
                            }
                        }
                    }
                    NEXT_CONTOUR: ;
                }
            }
            // at this point, all loops should be in contours[0]
            ExtrusionEntityCollection entities = traverse_loops(*this, contours.front(), thin_walls);
            // if brim will be printed, reverse the order of perimeters so that
            // we continue inwards after having finished the brim
            // TODO: add test for perimeter order
            if (this->config->external_perimeters_first || 
                (this->layer_id == 0 && this->object_config->brim_width.value > 0))
                entities.reverse();
            // append perimeters for this slice as a collection
            if (! entities.empty())
                this->loops->append(entities);
        } // for each loop of an island

        // fill gaps
        if (! gaps.empty()) {
            // collapse 
            double min = 0.2 * perimeter_width * (1 - INSET_OVERLAP_TOLERANCE);
            double max = 2. * perimeter_spacing;
            ExPolygons gaps_ex = diff_ex(
                //FIXME offset2 would be enough and cheaper.
                opening_ex(gaps, float(min / 2.)),
                offset2_ex(gaps, - float(max / 2.), float(max / 2. + ClipperSafetyOffset)));
            ThickPolylines polylines;
            for (const ExPolygon &ex : gaps_ex)
                ex.medial_axis(max, min, &polylines);
            if (! polylines.empty()) {
				ExtrusionEntityCollection gap_fill;
				variable_width(polylines, erGapFill, this->solid_infill_flow, gap_fill.entities);
                /*  Make sure we don't infill narrow parts that are already gap-filled
                    (we only consider this surface's gaps to reduce the diff() complexity).
                    Growing actual extrusions ensures that gaps not filled by medial axis
                    are not subtracted from fill surfaces (they might be too short gaps
                    that medial axis skips but infill might join with other infill regions
                    and use zigzag).  */
                //FIXME Vojtech: This grows by a rounded extrusion width, not by line spacing,
                // therefore it may cover the area, but no the volume.
                last = diff_ex(last, gap_fill.polygons_covered_by_width(10.f));
				this->gap_fill->append(std::move(gap_fill.entities));
			}
        }

        // create one more offset to be used as boundary for fill
        // we offset by half the perimeter spacing (to get to the actual infill boundary)
        // and then we offset back and forth by half the infill spacing to only consider the
        // non-collapsing regions
        coord_t inset = 
            (loop_number < 0) ? 0 :
            (loop_number == 0) ?
                // one loop
                ext_perimeter_spacing / 2 :
                // two or more loops?
                perimeter_spacing / 2;
        // only apply infill overlap if we actually have one perimeter
        if (inset > 0)
            inset -= coord_t(scale_(this->config->get_abs_value("infill_overlap", unscale<double>(inset + solid_infill_spacing / 2))));
        // simplify infill contours according to resolution
        Polygons pp;
        for (ExPolygon &ex : last)
            ex.simplify_p(m_scaled_resolution, &pp);
        // collapse too narrow infill areas
        coord_t min_perimeter_infill_spacing = coord_t(solid_infill_spacing * (1. - INSET_OVERLAP_TOLERANCE));
        // append infill areas to fill_surfaces
        this->fill_surfaces->append(
            offset2_ex(
                union_ex(pp),
                float(- inset - min_perimeter_infill_spacing / 2.),
                float(min_perimeter_infill_spacing / 2.)),
            stInternal);
    } // for each island
}

bool PerimeterGeneratorLoop::is_internal_contour() const
{
    // An internal contour is a contour containing no other contours
    if (! this->is_contour)
        return false;
    for (const PerimeterGeneratorLoop &loop : this->children)
        if (loop.is_contour)
            return false;
    return true;
}

}
