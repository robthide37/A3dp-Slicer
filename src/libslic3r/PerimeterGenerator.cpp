#include "clipper/clipper_z.hpp"

#include "PerimeterGenerator.hpp"

#include "BoundingBox.hpp"
#include "BridgeDetector.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Geometry.hpp"
#include "Geometry/MedialAxis.hpp"
#include "Geometry.hpp"
#include "Line.hpp"
#include "Milling/MillingPostProcess.hpp"
#include "Polygon.hpp"
#include "ShortestPath.hpp"
#include "SVG.hpp"

#include "Arachne/WallToolPaths.hpp"
#include "Arachne/utils/ExtrusionLine.hpp"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <list>
#include <stack>
#include <unordered_map>
#include <vector>

#include <boost/log/trivial.hpp>

//#define ARACHNE_DEBUG

#ifdef ARACHNE_DEBUG
#include "SVG.hpp"
#include "Utils.hpp"
#endif

namespace Slic3r {

//struct CollectionSimplifyVisitor : public ExtrusionVisitor
//{
//    ExtrusionEntityCollection *previous = nullptr;
//    virtual void default_use(ExtrusionEntity& entity) override {};
//    virtual void use(ExtrusionEntityCollection &coll) override
//    {
//        if (previous && coll.entities().size() == 1) {
//            previous->append_move_from(coll);
//        } else {
//            ExtrusionEntityCollection *old_previous = previous;
//            previous                                = &coll;
//            for (size_t i = 0; i < coll.entities().size(); ++i) { coll.set_entities()[i]->visit(*this); }
//            previous = old_previous;
//        }
//    }
//};

PerimeterGeneratorLoops get_all_Childs(PerimeterGeneratorLoop loop) {
    PerimeterGeneratorLoops ret;
    for (PerimeterGeneratorLoop &child : loop.children) {
        ret.push_back(child);
        PerimeterGeneratorLoops vals = get_all_Childs(child);
        ret.insert(ret.end(), vals.begin(), vals.end());
    }
    return ret;
}

// Thanks Cura developers for this function.
static void fuzzy_polygon(Polygon& poly, coordf_t fuzzy_skin_thickness, coordf_t fuzzy_skin_point_dist)
{
    const double min_dist_between_points = fuzzy_skin_point_dist * 3. / 4.; // hardcoded: the point distance may vary between 3/4 and 5/4 the supplied value
    const double range_random_point_dist = fuzzy_skin_point_dist / 2.;
    double dist_left_over = double(rand()) * (min_dist_between_points / 2) / double(RAND_MAX); // the distance to be traversed on the line before making the first new point
    Point* p0 = &poly.points.back();
    Points out;
    out.reserve(poly.points.size());
    for (Point& p1 : poly.points)
    { // 'a' is the (next) new point between p0 and p1
        Vec2d  p0p1 = (p1 - *p0).cast<double>();
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
    if (out.size() < 3) {
        size_t point_idx = std::max(size_t(0), poly.size() - 2);
        while (out.size() < 3) {
            out.emplace_back(poly[point_idx]);
            if (point_idx == 0)
                break;
            --point_idx;
        }
    }
    if (out.size() >= 3)
        poly.points = std::move(out);
}

// Thanks Cura developers for this function.
//supermerill: doesn't work
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
        out.front().p = out.back().p;

    if (out.size() >= 3)
        ext_lines.junctions = std::move(out);
}

void convert_to_clipperpath(const Polygons& source, ClipperLib_Z::Paths& dest) {
    dest.clear();
    dest.reserve(source.size());
    for (const Polygon& poly : source) {
        dest.emplace_back();
        ClipperLib_Z::Path& out = dest.back();
        out.reserve(poly.points.size());
        for (const Point& pt : poly.points)
            out.emplace_back(pt.x(), pt.y(), 0);
    }
}

#ifdef ARACHNE_DEBUG
static void export_perimeters_to_svg(const std::string &path, const Polygons &contours, const std::vector<Arachne::VariableWidthLines> &perimeters, const ExPolygons &infill_area)
{
    coordf_t    stroke_width = scale_(0.03);
    BoundingBox bbox         = get_extents(contours);
    bbox.offset(scale_(1.));
    ::Slic3r::SVG svg(path.c_str(), bbox);

    svg.draw(infill_area, "cyan");

    for (const Arachne::VariableWidthLines &perimeter : perimeters)
        for (const Arachne::ExtrusionLine &extrusion_line : perimeter) {
            ThickPolyline thick_polyline = to_thick_polyline(extrusion_line);
            svg.draw({thick_polyline}, "green", "blue", stroke_width);
        }

    for (const Line &line : to_lines(contours))
        svg.draw(line, "red", stroke_width);
}
#endif

// Thanks, Cura developers, for implementing an algorithm for generating perimeters with variable width (Arachne) that is based on the paper
// "A framework for adaptive width control of dense contour-parallel toolpaths in fused deposition modeling"
ProcessSurfaceResult PerimeterGenerator::process_arachne(int& loop_number, const Surface& surface) {

    ProcessSurfaceResult result;

    coord_t scaled_resolution = get_resolution(0, false, &surface);
    scaled_resolution = (scaled_resolution < SCALED_EPSILON ? SCALED_EPSILON : scaled_resolution);
    coord_t ext_displacement = (this->get_ext_perimeter_width() / 2. - this->get_ext_perimeter_spacing() / 2.);
    ExPolygons last = (ext_displacement != 0)
        ? offset_ex(surface.expolygon.simplify_p(scaled_resolution),  -ext_displacement)
        : union_ex(surface.expolygon.simplify_p(scaled_resolution));

    //increase surface for milling_post-process
    if (this->mill_extra_size > SCALED_EPSILON) {
        if (unmillable.empty())
            last = offset_ex(last, mill_extra_size);
        else {
            ExPolygons growth = diff_ex(offset_ex(last, mill_extra_size), unmillable, ApplySafetyOffset::Yes);
            last.insert(last.end(), growth.begin(), growth.end());
            last = union_ex(last);
        }
    }

    // only_one_perimeter_top, from orca
    std::vector<Arachne::VariableWidthLines> out_shell;
    if (loop_number > 0 && this->config->only_one_perimeter_top && !surface.has_mod_bridge() && upper_slices != nullptr) {
        this->throw_if_canceled();
        // Check if current layer has surfaces that are not covered by upper layer (i.e., top surfaces)
        ExPolygons non_top_polygons;
        ExPolygons fill_clip;

        //has to set the outer polygon to the centerline of the external perimeter
        split_top_surfaces(lower_slices, upper_slices, offset_ex(last, -this->ext_perimeter_spacing/2), result.top_fills, non_top_polygons, result.fill_clip);

        if (result.top_fills.empty()) {
            // No top surfaces, no special handling needed
        } else {
            // First we slice the outer shell
            const Polygons         last_p = to_polygons(last);
            Arachne::WallToolPaths wallToolPaths(last_p, this->get_ext_perimeter_spacing(),this->get_ext_perimeter_width(), 
                                                 this->get_perimeter_spacing(), this->get_perimeter_width(), 1, coord_t(0),
                                                 this->layer->height, *this->object_config, *this->print_config);
            out_shell = wallToolPaths.getToolPaths();
            // Make sure infill not overlap with wall
            // offset the InnerContour as arachne use bounds and not centerline
            result.top_fills = intersection_ex(result.top_fills, offset_ex(wallToolPaths.getInnerContour(), this->ext_perimeter_spacing/2));

            if (!result.top_fills.empty()) {
                // Then get the inner part that needs more walls
                // reduce the not-top fill to the bound for arachne (as arachne doesn't use the centerline but the boundary)
                // note: you can also diff_ex(offset_ex(result.top_fills, this->perimeter_spacing / 2), wallToolPaths.getInnerContour());  this should have similar results
                last = intersection_ex(offset_ex(non_top_polygons, -this->perimeter_spacing / 2), wallToolPaths.getInnerContour());
                //{
                //    static int i = 0;
                //    i++;
                //    std::stringstream stri;
                //    stri << this->layer->id() << "_M_" << i << "_only_one_peri"
                //         << ".svg";
                //    SVG svg(stri.str());
                //    //svg.draw(to_polylines(old_last), "green");
                //    //svg.draw(to_polylines(offset_ex(old_last, -this->ext_perimeter_spacing / 2)), "lime");
                //    //svg.draw(to_polylines(old_top), "blue");
                //    svg.draw(to_polylines(result.top_fills), "cyan");
                //    svg.draw(to_polylines(result.fill_clip), "pink");
                //    svg.draw(to_polylines(wallToolPaths.getInnerContour()), "orange");
                //    svg.draw(to_polylines(non_top_polygons), "red");
                //    svg.draw(to_polylines(last), "brown");
                //    svg.Close();
                //}
                loop_number--;
            } else {
                // Give up the outer shell because we don't have any meaningful top surface
                out_shell.clear();
            }
        }
    }

    const Polygons last_p = to_polygons(last);
    Arachne::WallToolPaths wallToolPaths(last_p, this->get_ext_perimeter_spacing(), this->get_ext_perimeter_width(), 
        this->get_perimeter_spacing(), this->get_perimeter_width(), loop_number, coord_t(0), 
        this->layer->height, *this->object_config, *this->print_config);
    std::vector<Arachne::VariableWidthLines> perimeters = wallToolPaths.getToolPaths();

#if _DEBUG
    for (auto perimeter : perimeters) {
        for (Arachne::ExtrusionLine &extrusion : perimeter) {
            if (extrusion.isZeroLength())
                continue;
            for (Slic3r::Arachne::ExtrusionJunction &junction : extrusion.junctions) {
                Point pt = junction.p;
                assert(unscaled(pt.x()) < 10000 && unscaled(pt.x()) > -10000);
                assert(unscaled(pt.y()) < 10000 && unscaled(pt.y()) > -10000);
            }
        }
    }
#endif
    // hack to fix points that go to the moon. https://github.com/supermerill/SuperSlicer/issues/4032
    // get max dist possible
    BoundingBox bb;
    for (ExPolygon &expo : last) bb.merge(expo.contour.points);
    const coordf_t max_dist = bb.min.distance_to_square(bb.max);
    //detect astray points and delete them
    for (Arachne::VariableWidthLines &perimeter : perimeters) {
        this->throw_if_canceled();
        for (auto it_extrusion = perimeter.begin(); it_extrusion != perimeter.end();) {
            Point last_point = bb.min;
            for (auto it_junction = it_extrusion->junctions.begin(); it_junction != it_extrusion->junctions.end();) {
                coordf_t dist = it_junction->p.distance_to_square(last_point);
                if (dist > max_dist) {
                    it_junction = it_extrusion->junctions.erase(it_junction);
                } else {
                    last_point = it_junction->p;
                    ++it_junction;
                }
            }
            if (it_extrusion->junctions.size() < 2) {
                it_extrusion = perimeter.erase(it_extrusion);
            } else {
                ++it_extrusion;
            }
        }
    }

    // only_one_perimeter_top, from orca
    if (!out_shell.empty()) {
        // Combine outer shells
        size_t inset_offset = 0;
        for (auto &p : out_shell) {
            for (auto &l : p) {
                if (l.inset_idx + 1 > inset_offset) {
                    inset_offset = l.inset_idx + 1;
                }
            }
        }
        for (auto &p : perimeters) {
            for (auto &l : p) { l.inset_idx += inset_offset; }
        }
        perimeters.insert(perimeters.begin(), out_shell.begin(), out_shell.end());
    }

    loop_number = int(perimeters.size());

#ifdef ARACHNE_DEBUG
        {
            static int iRun = 0;
            export_perimeters_to_svg(debug_out_path("arachne-perimeters-%d-%d.svg", layer_id, iRun++), to_polygons(last), perimeters, union_ex(wallToolPaths.getInnerContour()));
        }
#endif

    // All closed ExtrusionLine should have the same the first and the last point.
    // But in rare cases, Arachne produce ExtrusionLine marked as closed but without
    // equal the first and the last point.
    assert([&perimeters = std::as_const(perimeters)]() -> bool {
        for (const Arachne::VariableWidthLines& perimeter : perimeters)
            for (const Arachne::ExtrusionLine& el : perimeter)
                if (el.is_closed && el.junctions.front().p != el.junctions.back().p)
                    return false;
        return true;
    }());

    int start_perimeter = int(perimeters.size()) - 1;
    int end_perimeter = -1;
    int direction = -1;

    if (this->config->external_perimeters_first) {
        start_perimeter = 0;
        end_perimeter = int(perimeters.size());
        direction = 1;
    }

    std::vector<Arachne::ExtrusionLine*> all_extrusions;
    for (int perimeter_idx = start_perimeter; perimeter_idx != end_perimeter; perimeter_idx += direction) {
        if (perimeters[perimeter_idx].empty())
            continue;
        for (Arachne::ExtrusionLine& wall : perimeters[perimeter_idx])
            all_extrusions.emplace_back(&wall);
    }

    // Find topological order with constraints from extrusions_constrains.
    std::vector<size_t>              blocked(all_extrusions.size(), 0); // Value indicating how many extrusions it is blocking (preceding extrusions) an extrusion.
    std::vector<std::vector<size_t>> blocking(all_extrusions.size());   // Each extrusion contains a vector of extrusions that are blocked by this extrusion.
    std::unordered_map<const Arachne::ExtrusionLine*, size_t> map_extrusion_to_idx;
    for (size_t idx = 0; idx < all_extrusions.size(); idx++)
        map_extrusion_to_idx.emplace(all_extrusions[idx], idx);
    
    //TODO: order extrusion for contour/hole separatly
    bool reverse_order = this->config->external_perimeters_first
        || (this->object_config->brim_width.value > 0 && this->layer->id() == 0)
        || (this->object_config->brim_width_interior.value > 0 && this->layer->id() == 0);
    auto extrusions_constrains = Arachne::WallToolPaths::getRegionOrder(all_extrusions, reverse_order);
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
        this->throw_if_canceled();
        size_t best_candidate = 0;
        double best_distance_sqr = std::numeric_limits<double>::max();
        bool   is_best_closed = false;

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
            double      distance_sqr = (current_position - candidate_position).cast<double>().norm();
            if (distance_sqr < best_distance_sqr) { // Closer than the best candidate so far.
                if (path->is_closed || (!path->is_closed && best_distance_sqr != std::numeric_limits<double>::max()) || (!path->is_closed && !is_best_closed)) {
                    best_candidate = candidate_path_idx;
                    best_distance_sqr = distance_sqr;
                    is_best_closed = path->is_closed;
                }
            }
        }

        Arachne::ExtrusionLine* best_path = all_extrusions[best_candidate];
        ordered_extrusions.push_back({ best_path, best_path->is_contour(), false });
        processed[best_candidate] = true;
        for (size_t unlocked_idx : blocking[best_candidate])
            blocked[unlocked_idx]--;

        if (!best_path->junctions.empty()) { //If all paths were empty, the best path is still empty. We don't upate the current position then.
            if (best_path->is_closed)
                current_position = best_path->junctions[0].p; //We end where we started.
            else
                current_position = best_path->junctions.back().p; //Pick the other end from where we started.
        }
    }

    // fuzzify
    if (this->layer->id() > 0 && this->config->fuzzy_skin != FuzzySkinType::None) {
        std::vector<PerimeterGeneratorArachneExtrusion*> closed_loop_extrusions;
        for (PerimeterGeneratorArachneExtrusion& extrusion : ordered_extrusions)
            if (extrusion.extrusion->inset_idx == 0 || this->config->fuzzy_skin == FuzzySkinType::All) {
                if (extrusion.extrusion->is_closed && this->config->fuzzy_skin == FuzzySkinType::External) {
                    closed_loop_extrusions.emplace_back(&extrusion);
                } else {
                    extrusion.fuzzify = true;
                }
            }

        if (this->config->fuzzy_skin == FuzzySkinType::External) {
            ClipperLib_Z::Paths loops_paths;
            loops_paths.reserve(closed_loop_extrusions.size());
            for (const auto& cl_extrusion : closed_loop_extrusions) {
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

            for (const ClipperLib_Z::PolyNode* child_node : loops_polytree.Childs) {
                // The whole contour must have the same index.
                coord_t polygon_idx = child_node->Contour.front().z();
                bool    has_same_idx = std::all_of(child_node->Contour.begin(), child_node->Contour.end(),
                    [&polygon_idx](const ClipperLib_Z::IntPoint& point) -> bool { return polygon_idx == point.z(); });
                if (has_same_idx)
                    closed_loop_extrusions[polygon_idx]->fuzzify = true;
            }
        }
    }
    
        this->throw_if_canceled();
    if (ExtrusionEntityCollection extrusion_coll = _traverse_extrusions(ordered_extrusions); !extrusion_coll.empty()) {
        extrusion_coll.set_can_sort_reverse(false, false);
        this->loops->append(extrusion_coll);
    }

    ExPolygons    infill_contour = union_ex(wallToolPaths.getInnerContour());
    const coord_t spacing = (perimeters.size() == 1) ? ext_perimeter_spacing2 : perimeter_spacing;
    if (offset_ex(infill_contour, -float(spacing / 2.)).empty())
        infill_contour.clear(); // Infill region is too small, so let's filter it out.

    result.inner_perimeter = infill_contour;

    return result;
}

void PerimeterGenerator::split_top_surfaces(const ExPolygons *lower_slices,
                                            const ExPolygons *upper_slices,
                                            const ExPolygons &orig_polygons,
                                            ExPolygons &      top_fills,
                                            ExPolygons &      non_top_polygons,
                                            ExPolygons &      fill_clip)
{
    // other perimeters
    coord_t perimeter_width   = this->perimeter_flow.scaled_width();
    coord_t perimeter_spacing = this->perimeter_flow.scaled_spacing();

    // external perimeters
    coord_t ext_perimeter_width   = this->ext_perimeter_flow.scaled_width();
    coord_t ext_perimeter_spacing = this->ext_perimeter_flow.scaled_spacing();

    double  fill_nozzle_diameter = this->solid_infill_flow.nozzle_diameter();

    bool has_gap_fill = this->config->gap_fill_enabled &&
                        this->object_config->perimeter_generator.value == PerimeterGeneratorType::Classic;

    // split the polygons with top/not_top
    // get the offset from solid surface anchor*
    const int32_t peri_count = this->config->perimeters.value;
    const double max_perimeters_width = unscaled(double(ext_perimeter_width + perimeter_spacing * int(peri_count - int(1)))); 
    coord_t offset_top_surface = scale_t(config->external_infill_margin.get_abs_value(peri_count == 0 ? 0. : max_perimeters_width));
    // if possible, try to not push the extra perimeters inside the sparse infill
    if (offset_top_surface > 0.9 * (peri_count <= 1 ? 0. : (perimeter_spacing * (peri_count - 1))))
        offset_top_surface -= coord_t(0.9 * (peri_count <= 1 ? 0. : (perimeter_spacing * (peri_count - 1))));
    else
        offset_top_surface = 0;
    // don't takes into account too thin areas
    // skip if the exposed area is smaller than "min_width_top_surface"
    coordf_t min_width_top_surface = std::max(coordf_t(ext_perimeter_spacing / 2 + 10),
                                            scale_d(this->config->min_width_top_surface.get_abs_value(unscaled(perimeter_width))));

    Polygons grown_upper_slices;
    if (!this->config->only_one_perimeter_top_other_algo.value) {
        grown_upper_slices = offset(*upper_slices, min_width_top_surface);
    } else {
        ExPolygons grown_accumulator;
        // make thin upper surfaces disapear with -+offset_top_surface
        // do offset2 per island, to avoid big blob merging
        // remove polygon too thin (but don't mess with holes)
        for (const ExPolygon &expoly_to_grow : *this->upper_slices) {
            // only offset the contour, as it can merge holes
            Polygons contour = offset2(ExPolygons{ExPolygon{expoly_to_grow.contour}}, -offset_top_surface,
                                       offset_top_surface + min_width_top_surface +
                                           (this->mill_extra_size > SCALED_EPSILON ? (double) mill_extra_size : 0));
            if (!contour.empty()) {
                if (expoly_to_grow.holes.empty()) {
                    for (Polygon &p : contour) grown_accumulator.push_back(ExPolygon{p});
                } else {
                    Polygons holes = expoly_to_grow.holes;
                    for (Polygon &h : holes) h.reverse();
                    holes = offset(holes,
                                   -min_width_top_surface -
                                       ((this->mill_extra_size > SCALED_EPSILON) ? (double) mill_extra_size : 0));
                    for (ExPolygon p : diff_ex(contour, holes)) grown_accumulator.push_back(p);
                }
            }
        }
        grown_upper_slices = union_(grown_accumulator);
    }

    // get boungding box of last
    BoundingBox last_box = get_extents(orig_polygons);
    last_box.offset(SCALED_EPSILON);

    // get the Polygons upper the polygon this layer
    Polygons upper_polygons_series_clipped = ClipperUtils::clip_clipper_polygons_with_subject_bbox(grown_upper_slices, last_box);

    // set the clip to a virtual "second perimeter"
    fill_clip = offset_ex(orig_polygons, -coordf_t(ext_perimeter_spacing));
    // Check whether surface be bridge or not
    ExPolygons bridge_checker;
    if (lower_slices != nullptr) {
        // BBS: get the Polygons below the polygon this layer
        Polygons lower_polygons_series_clipped = ClipperUtils::clip_clipper_polygons_with_subject_bbox(*lower_slices, last_box);
        coordf_t bridge_offset = std::max(coordf_t(ext_perimeter_spacing), coordf_t(perimeter_width));
        // SoftFever: improve bridging
        const coordf_t bridge_margin = scale_d(this->config->bridged_infill_margin.get_abs_value(unscaled(perimeter_width)));
        bridge_checker = offset_ex(diff_ex(orig_polygons, lower_polygons_series_clipped, ApplySafetyOffset::Yes),
                                   1.5 * bridge_offset + bridge_margin + perimeter_spacing / 2);
    }
    ExPolygons delete_bridge = diff_ex(orig_polygons, bridge_checker, ApplySafetyOffset::Yes);
    // get the real top surface
    ExPolygons top_polygons;
    if (this->mill_extra_size < SCALED_EPSILON) {
        top_polygons = diff_ex(delete_bridge, upper_polygons_series_clipped, ApplySafetyOffset::Yes);
    } else if (this->unmillable.empty()) {
        top_polygons = diff_ex(delete_bridge, offset_ex(upper_polygons_series_clipped, (double) mill_extra_size), ApplySafetyOffset::Yes);
    } else {
        top_polygons = diff_ex(delete_bridge,
                               diff_ex(offset_ex(upper_polygons_series_clipped, (double) mill_extra_size), 
                                   unmillable, ApplySafetyOffset::Yes));
    }
    // save the top area for gap fill, or something. Made by BB/orca, but no comment.
    ExPolygons temp_gap = diff_ex(top_polygons, fill_clip);
    // get the not-top surface, from the "real top" but enlarged by external_infill_margin (and the
    // min_width_top_surface we removed a bit before)
    // also remove the ext_perimeter_spacing/2 width because we are faking the external perimeter, and we will remove ext_perimeter_spacing2
    ExPolygons inner_polygons = diff_ex(orig_polygons,
                                        offset_ex(top_polygons, offset_top_surface + min_width_top_surface - double(ext_perimeter_spacing / 2)),
                                        ApplySafetyOffset::Yes);
    // get the enlarged top surface, by using inner_polygons instead of upper_slices, and clip it for it to be exactly
    // the polygons to fill.
    top_polygons = diff_ex(fill_clip, inner_polygons, ApplySafetyOffset::Yes);
    // increase by half peri the inner space to fill the frontier between last and stored.
    top_fills = union_ex(top_fills, top_polygons);
    // set the clip to the external wall but go back inside by infill_extrusion_width/2 to be sure the extrusion won't
    // go outside even with a 100% overlap.
    double infill_spacing_unscaled = this->config->infill_extrusion_width.get_abs_value(fill_nozzle_diameter);
    if (infill_spacing_unscaled == 0)
        infill_spacing_unscaled = Flow::auto_extrusion_width(frInfill, fill_nozzle_diameter);
    fill_clip = offset_ex(orig_polygons, double(ext_perimeter_spacing / 2) - scale_(infill_spacing_unscaled / 2));

    non_top_polygons = intersection_ex(inner_polygons, orig_polygons);
    // Made by BB/orca, but no comment. Plz test it and report the usefullness.
    if (has_gap_fill)
        non_top_polygons = union_ex(non_top_polygons, temp_gap);

    //if (!top_fills.empty() && !non_top_polygons.empty()) {
    //    static int i = 0;
    //    i++;
    //    std::stringstream stri;
    //    stri << this->layer->id() << "_1_" << i << "_only_one_peri"
    //         << ".svg";
    //    SVG svg(stri.str());
    //    svg.draw(to_polylines(top_fills), "green");
    //    svg.draw(to_polylines(inner_polygons), "yellow");
    //    svg.draw(to_polylines(top_polygons), "cyan");
    //    svg.draw(to_polylines(orig_polygons), "orange");
    //    svg.draw(to_polylines(non_top_polygons), "red");
    //    svg.Close();
    //}
}

void PerimeterGenerator::process()
{
    // other perimeters
    this->m_mm3_per_mm = this->perimeter_flow.mm3_per_mm();
    this->perimeter_width = this->perimeter_flow.scaled_width();
    this->perimeter_spacing = this->perimeter_flow.scaled_spacing();

    // external perimeters
    this->m_ext_mm3_per_mm = this->ext_perimeter_flow.mm3_per_mm();
    this->ext_perimeter_width = this->ext_perimeter_flow.scaled_width();
    //spacing between two external perimeter (where you don't have the space to add other loops)
    this->ext_perimeter_spacing = this->ext_perimeter_flow.scaled_spacing();
    //spacing between external perimeter and the second
    this->ext_perimeter_spacing2 = ext_perimeter_spacing / 2 + perimeter_spacing / 2; //this->ext_perimeter_flow.scaled_spacing(this->perimeter_flow);

    // overhang perimeters
    this->m_mm3_per_mm_overhang = this->overhang_flow.mm3_per_mm();

    // solid infill
    this->solid_infill_spacing = this->solid_infill_flow.scaled_spacing();

    // infill gap to add vs perimeter (useful if using perimeter bonding)
    this->infill_gap = 0;

    this->round_peri = this->config->perimeter_round_corners.value;
    this->min_round_spacing = round_peri ? perimeter_width / 10 : 0;

    // perimeter bonding set.
    if (this->perimeter_flow.spacing_ratio() == 1
        && this->ext_perimeter_flow.spacing_ratio() == 1
        && this->config->external_perimeters_first
        && this->object_config->perimeter_bonding.value > 0) {
        this->infill_gap = (this->object_config->perimeter_bonding.get_abs_value(1)) * ext_perimeter_spacing;
        this->ext_perimeter_spacing2 -= infill_gap;
    }

    // Calculate the minimum required spacing between two adjacent traces.
    // This should be equal to the nominal flow spacing but we experiment
    // with some tolerance in order to avoid triggering medial axis when
    // some squishing might work. Loops are still spaced by the entire
    // flow spacing; this only applies to collapsing parts.
    // For ext_min_spacing we use the ext_perimeter_spacing calculated for two adjacent
    // external loops (which is the correct way) instead of using ext_perimeter_spacing2
    // which is the spacing between external and internal, which is not correct
    // and would make the collapsing (thus the details resolution) dependent on 
    // internal flow which is unrelated. <- i don't undertand, so revert to ext_perimeter_spacing2
    //const coord_t min_spacing     = (coord_t)( perimeter_spacing      * (1 - 0.05/*INSET_OVERLAP_TOLERANCE*/) );
    //const coord_t ext_min_spacing = (coord_t)( ext_perimeter_spacing2  * (1 - 0.05/*INSET_OVERLAP_TOLERANCE*/) );
    // now the tolerance is built in thin_perimeter settings

    // prepare grown lower layer slices for overhang detection
    if (this->lower_slices != NULL && (this->config->overhangs_width.value > 0 || this->config->overhangs_width_speed.value > 0)) {
        // We consider overhang any part where the entire nozzle diameter is not supported by the
        // lower layer, so we take lower slices and offset them by overhangs_width of the nozzle diameter used 
        // in the current layer

        //we use a range to avoid threshold issues.
        coord_t overhangs_width_flow = scale_t(config->overhangs_width.get_abs_value(this->overhang_flow.nozzle_diameter()));
        coord_t overhangs_width_speed = scale_t(config->overhangs_width_speed.get_abs_value(this->overhang_flow.nozzle_diameter()));
        coord_t min_feature = std::min(overhangs_width_flow, overhangs_width_speed) / 10;
        coord_t overhangs_width_flow_90 = coord_t(overhangs_width_flow * 0.99);
        coord_t overhangs_width_flow_110 = coord_t(overhangs_width_flow * 1.15);
        coord_t overhangs_width_speed_90 = coord_t(overhangs_width_speed * 0.99);
        coord_t overhangs_width_speed_110 = coord_t(overhangs_width_speed * 1.15);

        //flow offset should be greater than speed offset because the flow apply also the speed.
        //check if overhangs_width_speed is low enough to be relevant (if flow is activated)
        if (overhangs_width_flow > 0 && overhangs_width_speed + this->overhang_flow.nozzle_diameter() * 0.01 > overhangs_width_flow) {
            overhangs_width_speed = 0;
            overhangs_width_speed_90 = 0;
            overhangs_width_speed_110 = 0;
        }
        if (overhangs_width_flow > 0) {
            if (overhangs_width_flow_90 < overhangs_width_speed_110) {
                overhangs_width_speed_110 = overhangs_width_flow_90 = (overhangs_width_flow + overhangs_width_speed) / 2;
            }
        }

        if (overhangs_width_speed > 0 || overhangs_width_flow > 0) {
            ExPolygons simplified;
            //simplify the lower slices if too high (means low number) resolution (we can be very aggressive here)
            if (this->print_config->resolution < min_feature / 2) {
                for (const ExPolygon& expoly : *lower_slices) {
                    expoly.simplify(min_feature, &simplified);
                }
            }

            if (overhangs_width_speed > 0 && (overhangs_width_speed < overhangs_width_flow || overhangs_width_flow == 0)) {
                this->_lower_slices_bridge_speed_small = offset((simplified.empty() ? *this->lower_slices : simplified), (coordf_t)overhangs_width_speed_90 - (coordf_t)(ext_perimeter_width / 2));
                this->_lower_slices_bridge_speed_big = offset((simplified.empty() ? *this->lower_slices : simplified), (coordf_t)overhangs_width_speed_110 - (coordf_t)(ext_perimeter_width / 2));
                if (use_arachne) {
                    convert_to_clipperpath(this->_lower_slices_bridge_speed_small, this->_lower_slices_bridge_speed_small_clipperpaths);
                    convert_to_clipperpath(this->_lower_slices_bridge_speed_big, this->_lower_slices_bridge_speed_big_clipperpaths);
                }
            }
            if (overhangs_width_flow > 0) {
                if (overhangs_width_speed_110 == overhangs_width_flow_90 && overhangs_width_speed < overhangs_width_flow) {
                    this->_lower_slices_bridge_flow_small = this->_lower_slices_bridge_speed_big;
                    if (use_arachne) this->_lower_slices_bridge_flow_small_clipperpaths = this->_lower_slices_bridge_speed_big_clipperpaths;
                } else {
                    this->_lower_slices_bridge_flow_small = offset((simplified.empty() ? *this->lower_slices : simplified), (coordf_t)overhangs_width_flow_90 - (coordf_t)(ext_perimeter_width / 2));
                    if (use_arachne) convert_to_clipperpath(this->_lower_slices_bridge_flow_small, this->_lower_slices_bridge_flow_small_clipperpaths);
                }
                this->_lower_slices_bridge_flow_big = offset((simplified.empty() ? *this->lower_slices : simplified), (coordf_t)overhangs_width_flow_110 - (coordf_t)(ext_perimeter_width / 2));
                if (use_arachne) convert_to_clipperpath(this->_lower_slices_bridge_flow_big, this->_lower_slices_bridge_flow_big_clipperpaths);
            }
        }
    }
    this->throw_if_canceled();

    // have to grown the perimeters if mill post-process
    MillingPostProcess miller(this->slices, this->lower_slices, config, object_config, print_config);
    bool have_to_grow_for_miller = miller.can_be_milled(layer) && config->milling_extra_size.get_abs_value(1) > 0;
    this->mill_extra_size = 0;
    if (have_to_grow_for_miller) {
        this->unmillable = miller.get_unmillable_areas(layer);
        double spacing_vs_width = ext_perimeter_flow.width() - ext_perimeter_flow.spacing();
        this->mill_extra_size = scale_(config->milling_extra_size.get_abs_value(spacing_vs_width));
        have_to_grow_for_miller = this->mill_extra_size > SCALED_EPSILON;
    }

    // we need to process each island separately because we might have different
    // extra perimeters for each one
    Surfaces all_surfaces = this->slices->surfaces;

    processs_no_bridge(all_surfaces);

    int surface_idx = 0;
    const int extra_odd_perimeter = (config->extra_perimeters_odd_layers && layer->id() % 2 == 1 ? 1 : 0);
    for (const Surface& surface : all_surfaces) {
        surface_idx++;

        // detect how many perimeters must be generated for this island
        int nb_loop_contour = this->config->perimeters;
        if (nb_loop_contour > 0)
            nb_loop_contour += extra_odd_perimeter + surface.extra_perimeters;
        int nb_loop_holes = this->config->perimeters_hole;
        if (nb_loop_holes > 0)
            nb_loop_holes += extra_odd_perimeter + surface.extra_perimeters;
        
        if (nb_loop_contour < 0)
            nb_loop_contour = std::max(0, nb_loop_holes);
        if (nb_loop_holes < 0)
            nb_loop_holes = std::max(0, nb_loop_contour);

            if (print_config->spiral_vase) {
                if (layer->id() >= config->bottom_solid_layers) {
                    nb_loop_contour = 1;
                    nb_loop_holes = 0;
                }
            }

        if ((layer->id() == 0 && this->config->only_one_perimeter_first_layer) ||
            (this->config->only_one_perimeter_top && this->upper_slices == NULL)) {
            nb_loop_contour = std::min(nb_loop_contour, 1);
            nb_loop_holes = std::min(nb_loop_holes, 1);
        }

        ProcessSurfaceResult surface_process_result;
        //core generation
        if (use_arachne && nb_loop_holes == nb_loop_contour) {
            surface_process_result = process_arachne(nb_loop_contour, surface);
            nb_loop_holes = nb_loop_contour; // nb_loop_contour is in/out
        } else {
            surface_process_result = process_classic(nb_loop_contour, nb_loop_holes, surface);
        }
        this->throw_if_canceled();


        // create one more offset to be used as boundary for fill
        // we offset by half the perimeter spacing (to get to the actual infill boundary)
        // and then we offset back and forth by half the infill spacing to only consider the
        // non-collapsing regions
        coord_t infill_peri_overlap = 0;
        // only apply infill overlap if we actually have one perimeter
        if (nb_loop_contour > 0 || nb_loop_holes > 0) {
            //give the overlap size to let the infill do his overlap
            //add overlap if at least one perimeter
            coordf_t perimeter_spacing_for_encroach = 0;
            if(this->config->perimeters == 1)
                perimeter_spacing_for_encroach = this->ext_perimeter_flow.spacing();
            else if(this->config->only_one_perimeter_top.value)
                //note: use the min of the two to avoid overextrusion if only one perimeter top
                // TODO: only do that if there is a top & a not-top surface
                perimeter_spacing_for_encroach = std::min(this->perimeter_flow.spacing(), this->ext_perimeter_flow.spacing());
            else //if(layerm->region().config().perimeters > 1)
                perimeter_spacing_for_encroach = this->perimeter_flow.spacing();
            infill_peri_overlap = scale_t(this->config->get_abs_value("infill_overlap", perimeter_spacing_for_encroach));
        }

        // simplify infill contours according to resolution
        Polygons not_filled_p;
        for (const ExPolygon& ex : surface_process_result.inner_perimeter)
            ex.simplify_p(scale_t(std::max(this->print_config->resolution.value, print_config->resolution_internal / 4)), &not_filled_p);
        ExPolygons not_filled_exp = union_ex(not_filled_p);
        // collapse too narrow infill areas
        coord_t min_perimeter_infill_spacing = (coord_t)(this->get_solid_infill_spacing() * (1. - INSET_OVERLAP_TOLERANCE));
        ExPolygons infill_exp;
        infill_exp = offset2_ex(not_filled_exp,
            double(- min_perimeter_infill_spacing / 2 + infill_peri_overlap - this->get_infill_gap()),
            double(min_perimeter_infill_spacing / 2));
        //special branch if gap : don't inset away from gaps!
        ExPolygons gap_fill_exps;
        if (!surface_process_result.gap_srf.empty()) {
            //not_filled_exp = union_ex(not_filled_p);
            infill_exp = offset2_ex(not_filled_exp,
                double(- min_perimeter_infill_spacing / 2 + infill_peri_overlap - this->get_infill_gap()),
                double(min_perimeter_infill_spacing / 2));
            //remove gaps surfaces
            not_filled_p.clear();
            for (ExPolygon& ex : surface_process_result.gap_srf)
                ex.simplify_p(scale_t(std::max(this->print_config->resolution.value, print_config->resolution_internal / 4)), &not_filled_p);
            gap_fill_exps = union_ex(not_filled_p);
            gap_fill_exps = offset_ex(gap_fill_exps, -infill_peri_overlap);
            infill_exp = diff_ex(infill_exp, gap_fill_exps);
        }
        
        this->throw_if_canceled();
        //if any top_fills, grow them by ext_perimeter_spacing/2 to have the real un-anchored fill
        ExPolygons top_infill_exp = intersection_ex(surface_process_result.fill_clip, offset_ex(surface_process_result.top_fills, double(this->get_ext_perimeter_spacing() / 2)));
        if (!surface_process_result.top_fills.empty()) {
            infill_exp = union_ex(infill_exp, offset_ex(top_infill_exp, double(infill_peri_overlap)));
        }
        // append infill areas to fill_surfaces
        this->fill_surfaces->append(infill_exp, stPosInternal | stDensSparse);

        if (infill_peri_overlap != 0) {
            ExPolygons polyWithoutOverlap;
            if (min_perimeter_infill_spacing / 2 > infill_peri_overlap)
                polyWithoutOverlap = offset2_ex(
                    not_filled_exp,
                    double(- infill_gap - min_perimeter_infill_spacing / 2 + infill_peri_overlap),
                    double(min_perimeter_infill_spacing / 2 - infill_peri_overlap));
            else
                polyWithoutOverlap = offset_ex(
                    not_filled_exp,
                    double(- this->get_infill_gap()));
            if (!gap_fill_exps.empty()) {
                polyWithoutOverlap = diff_ex(polyWithoutOverlap, gap_fill_exps);
            }
            if (!surface_process_result.top_fills.empty()) {
                polyWithoutOverlap = union_ex(polyWithoutOverlap, top_infill_exp);
            }
            this->fill_no_overlap.insert(this->fill_no_overlap.end(), polyWithoutOverlap.begin(), polyWithoutOverlap.end());
            /*{
                static int isaqsdsdfsdfqzfn = 0;
                std::stringstream stri;
                stri << this->layer->id() << "_2_end_makeperimeter_" << isaqsdsdfsdfqzfn++ << ".svg";
                SVG svg(stri.str());
                svg.draw(to_polylines(infill_exp), "blue");
                svg.draw(to_polylines(fill_no_overlap), "cyan");
                svg.draw(to_polylines(not_filled_exp), "green");
                svg.draw(to_polylines(last_no_gaps), "yellow");
                //svg.draw(to_polylines(offset_ex(surface_process_result.fill_clip, ext_perimeter_spacing / 2)), "brown");
                svg.draw(to_polylines(top_infill_exp), "orange");
                svg.Close();
            }*/
        }
    }

}

void PerimeterGenerator::processs_no_bridge(Surfaces& all_surfaces) {
    //store surface for bridge infill to avoid unsupported perimeters (but the first one, this one is always good)
    if (this->config->no_perimeter_unsupported_algo != npuaNone
        && this->lower_slices != NULL && !this->lower_slices->empty()) {
        coordf_t bridged_infill_margin = scale_d(this->config->bridged_infill_margin.get_abs_value(unscaled(this->ext_perimeter_width)));

        for (size_t surface_idx = 0; surface_idx < all_surfaces.size(); surface_idx++) {
            Surface* surface = &all_surfaces[surface_idx];
            ExPolygons last = { surface->expolygon };
            //compute our unsupported surface
            ExPolygons unsupported = diff_ex(last, *this->lower_slices, ApplySafetyOffset::Yes);
            if (!unsupported.empty()) {
                //remove small overhangs
                ExPolygons unsupported_filtered = offset2_ex(unsupported, double(-this->get_perimeter_spacing()), double(this->get_perimeter_spacing()));
                if (!unsupported_filtered.empty()) {
                    //to_draw.insert(to_draw.end(), last.begin(), last.end());
                    //extract only the useful part of the lower layer. The safety offset is really needed here.
                    ExPolygons support = diff_ex(last, unsupported, ApplySafetyOffset::Yes);
                    if (!unsupported.empty()) {
                        //only consider the part that can be bridged (really, by the bridge algorithm)
                        //first, separate into islands (ie, each ExPlolygon)
                        int numploy = 0;
                        //only consider the bottom layer that intersect unsupported, to be sure it's only on our island.
                        ExPolygonCollection lower_island(support);
                        //a detector per island
                        ExPolygons bridgeable;
                        for (ExPolygon unsupported : unsupported_filtered) {
                            BridgeDetector detector( unsupported,
                                lower_island.expolygons,
                                this->overhang_flow.scaled_spacing(),
                                scale_t(this->print_config->bridge_precision.get_abs_value(this->overhang_flow.spacing())),
                                this->layer->id());
                            if (detector.detect_angle(Geometry::deg2rad(this->config->bridge_angle.value)))
                                expolygons_append(bridgeable, union_ex(detector.coverage()));
                        }
                        if (!bridgeable.empty()) {
                            //check if we get everything or just the bridgeable area
                            if (this->config->no_perimeter_unsupported_algo.value == npuaNoPeri || this->config->no_perimeter_unsupported_algo.value == npuaFilled) {
                                //we bridge everything, even the not-bridgeable bits
                                for (size_t i = 0; i < unsupported_filtered.size();) {
                                    ExPolygon& poly_unsupp = *(unsupported_filtered.begin() + i);
                                    Polygons contour_simplified = poly_unsupp.contour.simplify(perimeter_spacing);
                                    ExPolygon poly_unsupp_bigger = poly_unsupp;
                                    Polygons contour_bigger = offset(poly_unsupp_bigger.contour, bridged_infill_margin);
                                    if (contour_bigger.size() == 1) poly_unsupp_bigger.contour = contour_bigger[0];

                                    //check convex, has some bridge, not overhang
                                    if (contour_simplified.size() == 1 && contour_bigger.size() == 1 && contour_simplified[0].concave_points().size() == 0
                                        && intersection_ex(bridgeable, ExPolygons{ poly_unsupp }).size() > 0
                                        && diff_ex(ExPolygons{ poly_unsupp_bigger }, union_ex(last, offset_ex(bridgeable, bridged_infill_margin + perimeter_spacing / 2)), ApplySafetyOffset::Yes).size() == 0
                                        ) {
                                        //ok, keep it
                                        i++;
                                    } else {
                                        unsupported_filtered.erase(unsupported_filtered.begin() + i);
                                    }
                                }
                                unsupported_filtered = intersection_ex(last,
                                    offset2_ex(unsupported_filtered, double(-perimeter_spacing / 2), double(bridged_infill_margin + perimeter_spacing / 2)));
                                if (this->config->no_perimeter_unsupported_algo.value == npuaFilled) {
                                    for (ExPolygon& expol : unsupported_filtered) {
                                        //check if the holes won't be covered by the upper layer
                                        //TODO: if we want to do that, we must modify the geometry before making perimeters.
                                        //if (this->upper_slices != nullptr && !this->upper_slices->expolygons.empty()) {
                                        //    for (Polygon &poly : expol.holes) poly.make_counter_clockwise();
                                        //    float perimeterwidth = this->config->perimeters == 0 ? 0 : (this->ext_perimeter_flow.scaled_width() + (this->config->perimeters - 1) + this->perimeter_flow.scaled_spacing());
                                        //    std::cout << "test upper slices with perimeterwidth=" << perimeterwidth << "=>" << offset_ex(this->upper_slices->expolygons, -perimeterwidth).size();
                                        //    if (intersection(Polygons() = { expol.holes }, to_polygons(offset_ex(this->upper_slices->expolygons, -this->ext_perimeter_flow.scaled_width() / 2))).empty()) {
                                        //        std::cout << " EMPTY";
                                        //        expol.holes.clear();
                                        //    } else {
                                        //    }
                                        //    std::cout << "\n";
                                        //} else {
                                        expol.holes.clear();
                                        //}

                                        //detect inside volume
                                        for (size_t surface_idx_other = 0; surface_idx_other < all_surfaces.size(); surface_idx_other++) {
                                            if (surface_idx == surface_idx_other) continue;
                                            if (intersection_ex(ExPolygons() = { expol }, ExPolygons() = { all_surfaces[surface_idx_other].expolygon }).size() > 0) {
                                                //this means that other_surf was inside an expol holes
                                                //as we removed them, we need to add a new one
                                                ExPolygons new_poly = offset2_ex(ExPolygons{ all_surfaces[surface_idx_other].expolygon }, double(-bridged_infill_margin - perimeter_spacing), double(perimeter_spacing));
                                                if (new_poly.size() == 1) {
                                                    all_surfaces[surface_idx_other].expolygon = new_poly[0];
                                                    expol.holes.push_back(new_poly[0].contour);
                                                    expol.holes.back().make_clockwise();
                                                } else {
                                                    for (size_t idx = 0; idx < new_poly.size(); idx++) {
                                                        Surface new_surf = all_surfaces[surface_idx_other];
                                                        new_surf.expolygon = new_poly[idx];
                                                        all_surfaces.push_back(new_surf);
                                                        expol.holes.push_back(new_poly[idx].contour);
                                                        expol.holes.back().make_clockwise();
                                                    }
                                                    all_surfaces.erase(all_surfaces.begin() + surface_idx_other);
                                                    if (surface_idx_other < surface_idx) {
                                                        surface_idx--;
                                                        surface = &all_surfaces[surface_idx];
                                                    }
                                                    surface_idx_other--;
                                                }
                                            }
                                        }
                                    }

                                }
                                //TODO: add other polys as holes inside this one (-margin)
                            } else if (this->config->no_perimeter_unsupported_algo.value == npuaBridgesOverhangs || this->config->no_perimeter_unsupported_algo.value == npuaBridges) {
                                //simplify to avoid most of artefacts from printing lines.
                                ExPolygons bridgeable_simplified;
                                for (ExPolygon& poly : bridgeable) {
                                    poly.simplify(perimeter_spacing, &bridgeable_simplified);
                                }
                                bridgeable_simplified = offset2_ex(bridgeable_simplified, -ext_perimeter_width, ext_perimeter_width);
                                //bridgeable_simplified = intersection_ex(bridgeable_simplified, unsupported_filtered);
                                //offset by perimeter spacing because the simplify may have reduced it a bit.
                                //it's not dangerous as it will be intersected by 'unsupported' later
                                //FIXME: add overlap in this->fill_surfaces->append
                                //FIXME: it overlap inside unsuppported not-bridgeable area!

                                //bridgeable_simplified = offset2_ex(bridgeable_simplified, (double)-perimeter_spacing, (double)perimeter_spacing * 2);
                                //ExPolygons unbridgeable = offset_ex(diff_ex(unsupported, bridgeable_simplified), perimeter_spacing * 3 / 2);
                                //ExPolygons unbridgeable = intersection_ex(unsupported, diff_ex(unsupported_filtered, offset_ex(bridgeable_simplified, ext_perimeter_width / 2)));
                                //unbridgeable = offset2_ex(unbridgeable, -ext_perimeter_width, ext_perimeter_width);


                                if (this->config->no_perimeter_unsupported_algo.value == npuaBridges) {
                                    ExPolygons unbridgeable = unsupported_filtered;
                                    for (ExPolygon& expol : unbridgeable)
                                        expol.holes.clear();
                                    unbridgeable = diff_ex(unbridgeable, bridgeable_simplified);
                                    unbridgeable = offset2_ex(unbridgeable, -ext_perimeter_width * 2, ext_perimeter_width * 2);
                                    ExPolygons bridges_temp = offset2_ex(intersection_ex(last, diff_ex(unsupported_filtered, unbridgeable), ApplySafetyOffset::Yes), -ext_perimeter_width / 4, ext_perimeter_width / 4);
                                    //remove the overhangs section from the surface polygons
                                    ExPolygons reference = last;
                                    last = diff_ex(last, unsupported_filtered);
                                    //ExPolygons no_bridge = diff_ex(offset_ex(unbridgeable, ext_perimeter_width * 3 / 2), last);
                                    //bridges_temp = diff_ex(bridges_temp, no_bridge);
                                    coordf_t offset_to_do = bridged_infill_margin;
                                    bool first = true;
                                    unbridgeable = diff_ex(unbridgeable, offset_ex(bridges_temp, ext_perimeter_width));
                                    while (offset_to_do > ext_perimeter_width * 1.5) {
                                        unbridgeable = offset2_ex(unbridgeable, -ext_perimeter_width / 4, ext_perimeter_width * 2.25, ClipperLib::jtSquare);
                                        bridges_temp = diff_ex(bridges_temp, unbridgeable);
                                        bridges_temp = offset_ex(bridges_temp, ext_perimeter_width, ClipperLib::jtMiter, 6.);
                                        unbridgeable = diff_ex(unbridgeable, offset_ex(bridges_temp, ext_perimeter_width));
                                        offset_to_do -= ext_perimeter_width;
                                        first = false;
                                    }
                                    unbridgeable = offset_ex(unbridgeable, ext_perimeter_width + offset_to_do, ClipperLib::jtSquare);
                                    bridges_temp = diff_ex(bridges_temp, unbridgeable);
                                    unsupported_filtered = offset_ex(bridges_temp, offset_to_do);
                                    unsupported_filtered = intersection_ex(unsupported_filtered, reference);
                                } else {
                                    ExPolygons unbridgeable = intersection_ex(unsupported, diff_ex(unsupported_filtered, offset_ex(bridgeable_simplified, ext_perimeter_width / 2)));
                                    unbridgeable = offset2_ex(unbridgeable, -ext_perimeter_width, ext_perimeter_width);
                                    unsupported_filtered = unbridgeable;

                                    ////put the bridge area inside the unsupported_filtered variable
                                    //unsupported_filtered = intersection_ex(last,
                                    //    diff_ex(
                                    //    offset_ex(bridgeable_simplified, (double)perimeter_spacing / 2),
                                    //    unbridgeable
                                    //    )
                                    //    );
                                }
                            } else {
                                unsupported_filtered.clear();
                            }
                        } else {
                            unsupported_filtered.clear();
                        }
                    }

                    if (!unsupported_filtered.empty()) {

                        //add this directly to the infill list.
                        // this will avoid to throw wrong offsets into a good polygons
                        this->fill_surfaces->append(
                            unsupported_filtered,
                            stPosInternal | stDensSparse);

                        // store the results
                        last = diff_ex(last, unsupported_filtered, ApplySafetyOffset::Yes);
                        //remove "thin air" polygons (note: it assumes that all polygons below will be extruded)
                        for (int i = 0; i < last.size(); i++) {
                            if (intersection_ex(support, ExPolygons() = { last[i] }).empty()) {
                                this->fill_surfaces->append(
                                    ExPolygons() = { last[i] },
                                    stPosInternal | stDensSparse);
                                last.erase(last.begin() + i);
                                i--;
                            }
                        }
                    }
                }
            }
            if (last.size() == 0) {
                all_surfaces.erase(all_surfaces.begin() + surface_idx);
                surface_idx--;
            } else {
                surface->expolygon = last[0];
                for (size_t idx = 1; idx < last.size(); idx++) {
                    all_surfaces.emplace_back(*surface, last[idx]);
                }
            }
        }
    }
}

Polygons get_contours(const ExPolygons &expolys) {
    Polygons polys;
    for (const ExPolygon &expoly : expolys) {
        assert(expoly.contour.is_counter_clockwise());
        polys.push_back(expoly.contour);
    }
    return polys;
}

Polygons as_contour(const Polygons &holes) {
    Polygons out;
    for (const Polygon &hole : holes) {
        assert(hole.is_clockwise());
        out.push_back(hole);
        out.back().make_counter_clockwise();
    }
    return out;
}

Polygons get_holes_as_contour(const ExPolygon &expoly) {
    Polygons polys;
    for(const Polygon hole : expoly.holes){
        assert(hole.is_clockwise());
        polys.push_back(hole);
        polys.back().make_counter_clockwise();
    }
    return polys;
}

Polygons get_holes_as_contour(const ExPolygons &expolys) {
    Polygons polys;
    for(const ExPolygon &expoly : expolys) 
        for(const Polygon hole : expoly.holes){
            assert(hole.is_clockwise());
            polys.push_back(hole);
            polys.back().make_counter_clockwise();
        }
    return polys;
}

// expolygon representing the perimeter path
struct ExPolygonAsynch
{
    enum ExPolygonAsynchType {
        epatGrowHole,
        epatShrinkContour
    };
    ExPolygonAsynchType type;
    ExPolygon expoly;
    // shrink the contour by this value to get the end of the spacing (should be negative, to shrink from centerline or edge)
    coordf_t  offset_contour_inner;
    // shrink the contour by this value to get the external shell (the spacing position) (can be negative to grow from centreline, and be positive to shrink from surface polygon)
    coordf_t  offset_contour_outer;
    // grow the holes by this value to get the end of the spacing (should be negative, to grow from centerline or edge)
    coordf_t  offset_holes_inner;
    // grow the holes by this value to get the external shell (the spacing position) (should be the same value as offset_contour_outer)
    coordf_t  offset_holes_outer;
    
};

// next_onion can be partially filled
void grow_holes_only(std::vector<ExPolygonAsynch> &unmoveable_contours,
                     ExPolygons &                      next_onion,
                     coordf_t                          spacing,
                     coordf_t                          overlap_spacing,
                     bool                              round_peri,
                     float                             min_round_spacing = 3.f)
{
    assert(spacing > 0);
    assert(overlap_spacing >= 0);
    Polygons new_contours;
    for (size_t idx_unmoveable = 0; idx_unmoveable < unmoveable_contours.size(); ++idx_unmoveable) {
        ExPolygonAsynch & unmoveable_contour = unmoveable_contours[idx_unmoveable];
        assert(unmoveable_contour.type == ExPolygonAsynch::ExPolygonAsynchType::epatGrowHole);
        ExPolygon &expoly = unmoveable_contour.expoly;
        assert(unmoveable_contour.offset_holes_inner <= 0);
        // grow fake contours, can now have fake holes and/or less fake contours.
        Polygons ok_holes = offset(get_holes_as_contour(expoly),
                                   -unmoveable_contour.offset_holes_inner + spacing / 2 + overlap_spacing,
                                        (round_peri ? ClipperLib::JoinType::jtRound :
                                                        ClipperLib::JoinType::jtMiter),
                                        (round_peri ? min_round_spacing : 3));
        for (size_t i = 0; i < ok_holes.size(); ++i) {
            if (ok_holes[i].is_clockwise()) {
                // hole, it's a new peri, move it.
                new_contours.push_back(std::move(ok_holes[i]));
                ok_holes.erase(ok_holes.begin() + i);
                new_contours.back().make_counter_clockwise();
                i--;
            }
        }
        ok_holes = union_(ok_holes);
        for (const Polygon &p : ok_holes) assert(p.is_counter_clockwise());
        //shrink contour, can now be multiple contour.
        coordf_t computed_offset = unmoveable_contour.offset_contour_inner;
        computed_offset -= spacing / 2;
        computed_offset -= overlap_spacing;
        Polygons ex_contour_offset = offset(Polygons{expoly.contour}, computed_offset);
        bool ex_contour_offset_now_fake_hole = false;
        for (size_t idx_hole = 0; idx_hole < ok_holes.size(); ++idx_hole) {
            const Polygon &hole = ok_holes[idx_hole];
            assert(hole.is_counter_clockwise());
            // Check if it can fuse with contour
            // TODO: bounding box for quicker cut search
            auto it_contour_candidate_for_fuse = ex_contour_offset.begin();
            Polygons fused_contour;
            while (it_contour_candidate_for_fuse != ex_contour_offset.end()) {
                ExPolygons result = diff_ex(Polygons{*it_contour_candidate_for_fuse}, Polygons{hole});
                // Only two options here, it can fuse and then there is 1 or more contour, no holes.
                // Or it don't touch the contour and so nothing happen. (the hole can be inside or outside)
                // SO, we can check it it slip or if the contour has been modified
                if (result.size() > 1 || (result.size() == 1 && result.front().contour != *it_contour_candidate_for_fuse)) {
                    for (ExPolygon &expoly : result) assert(expoly.holes.empty());
                    // now use this one.
                    append(fused_contour, to_polygons(result));
                    ex_contour_offset_now_fake_hole = true;
                    // remove from useful holes
                    ok_holes.erase(ok_holes.begin() + idx_hole);
                    idx_hole--;
                    it_contour_candidate_for_fuse = ex_contour_offset.erase(it_contour_candidate_for_fuse);
                } else {
                    ++it_contour_candidate_for_fuse;
                }
            }
            if(!fused_contour.empty())
                append(ex_contour_offset, std::move(fused_contour));
        }
        // if moved from unmoveable_contours to growing_contours, then move the expoly
        if (ex_contour_offset_now_fake_hole) {
            // add useful holes to the contours, and push them
            if (overlap_spacing != 0)
                append(next_onion, offset_ex(diff_ex(ex_contour_offset, ok_holes), overlap_spacing));
            else
                append(next_onion, diff_ex(ex_contour_offset, ok_holes));
            // remove from unmoveable
            unmoveable_contours.erase(unmoveable_contours.begin() + idx_unmoveable);
            idx_unmoveable--;
        } else {
            // update holes
            // shrink to centerline
            if (overlap_spacing != 0)
                ok_holes = offset(ok_holes, -overlap_spacing);
            /*test*/ for (const Polygon &p : ok_holes) assert(p.is_counter_clockwise());
            ExPolygons new_unmoveable = diff_ex(Polygons{unmoveable_contour.expoly.contour}, ok_holes);
            // check if size is good. It's not possible to split the peri: it isn't srhunk, and holes intersect are alreedy detected (not unmoveable anymore)
            assert(new_unmoveable.size() <= 1);
            if (new_unmoveable.empty()) {
                // remove from unmoveable
                unmoveable_contours.erase(unmoveable_contours.begin() + idx_unmoveable);
                idx_unmoveable--;
            } else if(new_unmoveable.size() == 1){
                // update
                unmoveable_contour.expoly = new_unmoveable.front();
                unmoveable_contour.offset_holes_inner = -spacing / 2;
                unmoveable_contour.offset_holes_outer = spacing / 2;
            } else {
                assert(false);
                //add all
                for(ExPolygon & new_expoly : new_unmoveable)
                    unmoveable_contours.push_back({unmoveable_contour.type, new_expoly, unmoveable_contour.offset_contour_inner,
                                                   unmoveable_contour.offset_contour_outer, -spacing / 2,
                                                   spacing / 2});
                // remove from unmoveable
                unmoveable_contours.erase(unmoveable_contours.begin() + idx_unmoveable);
                idx_unmoveable--;
            }
        }
    }
}


// next_onion can be partially filled
void grow_contour_only(std::vector<ExPolygonAsynch> &unmoveable_holes, coordf_t spacing, coordf_t overlap_spacing, bool round_peri, float min_round_spacing = 3.f) {
    assert(spacing > 0);
    assert(overlap_spacing >= 0);
    Polygons new_contours;
    // mutable size to allow insert at the same time.
    size_t unmoveable_holes_size = unmoveable_holes.size();
    for (size_t idx_unmoveable = 0; idx_unmoveable < unmoveable_holes_size; ++idx_unmoveable) {
        ExPolygonAsynch & unmoveable_hole = unmoveable_holes[idx_unmoveable];
        assert(unmoveable_hole.type == ExPolygonAsynch::ExPolygonAsynchType::epatShrinkContour);
        ExPolygon &expoly = unmoveable_hole.expoly;
        // shrink contour, can now have more contours.
        assert(unmoveable_hole.offset_contour_inner <=0);
        Polygons ok_contours = offset(expoly.contour, unmoveable_hole.offset_contour_inner - spacing/2 - overlap_spacing,
                                        (round_peri ? ClipperLib::JoinType::jtRound :
                                                        ClipperLib::JoinType::jtMiter),
                                        (round_peri ? min_round_spacing : 3));
        //we shrunk -> new peri can appear, holes can disapear, but there is already none.
        for (const Polygon &p : ok_contours) assert(p.is_counter_clockwise());
        //grow holes to right size
        assert(-unmoveable_hole.offset_holes_inner + spacing/2 - overlap_spacing > 0);
        Polygons original_holes = get_holes_as_contour(expoly);
        Polygons offsetted_holes = offset(original_holes, -unmoveable_hole.offset_holes_inner + spacing/2 + overlap_spacing);
        // remove fake periemter, i don't want them.
        for (size_t i = 0; i < offsetted_holes.size(); ++i) {
            if (offsetted_holes[i].is_clockwise()) {
                offsetted_holes.erase(offsetted_holes.begin() + i);
                new_contours.back().make_counter_clockwise();
                i--;
            }
        }
        offsetted_holes = union_(offsetted_holes);
        for (const Polygon &p : offsetted_holes) assert(p.is_counter_clockwise());

        for (Polygon simple_contour : ok_contours) {
            // remove holes
            ExPolygons test_expoly = diff_ex(Polygons{simple_contour}, offsetted_holes);
            if (overlap_spacing != 0) {
                test_expoly = offset_ex(test_expoly, overlap_spacing);
            }
            if (test_expoly.size() == 1) {
                // no merge, then i can use the right hole size
                ExPolygons new_unmoveable_hole = diff_ex(Polygons{test_expoly[0].contour}, original_holes);
                // diff with smaller holes, so it has to be only one contour.
                assert(new_unmoveable_hole.size() == 1);
                expoly                               = new_unmoveable_hole[0];
                unmoveable_hole.offset_contour_inner = -spacing / 2;
                unmoveable_hole.offset_contour_outer = spacing / 2;
            } else {
                // a hole cut it, or clear it.
                for (ExPolygon &new_expoly : test_expoly) {
                    ExPolygons new_unmoveable_holes = diff_ex(Polygons{new_expoly.contour}, original_holes);
                    for(ExPolygon & new_unmoveable_hole : new_unmoveable_holes)
                        unmoveable_holes.push_back({unmoveable_hole.type, new_unmoveable_hole, -spacing / 2, spacing / 2,
                                                unmoveable_hole.offset_holes_inner, unmoveable_hole.offset_holes_outer});
                }
                unmoveable_holes.erase(unmoveable_holes.begin() + idx_unmoveable);
                idx_unmoveable--;
                unmoveable_holes_size--;
            }
        }
        //we shrink periemter, so it doesn't create holes, so we don't have anythign to add to next_onion.
    }
}

ProcessSurfaceResult PerimeterGenerator::process_classic(int& contour_count, int& holes_count, const Surface& surface)
{
    ProcessSurfaceResult results;
    ExPolygons gaps;
    //this var store infill surface removed from last to not add any more perimeters to it.
    // simplification already done at slicing
    //simplify the loop to avoid artifacts when shrinking almost-0 segments
    coord_t resolution = get_resolution(0, false, &surface);
    ExPolygons last    = union_ex(surface.expolygon.simplify_p((resolution < SCALED_EPSILON ? SCALED_EPSILON : resolution)));
    double last_area   = -1;

    // list of Expolygons where contour or holes aren't growing.
    std::vector<ExPolygonAsynch> last_asynch;
    bool last_asynch_initialized = false;

    if (contour_count > 0 || holes_count > 0) {

        //increase surface for milling_post-process
        if (this->mill_extra_size > SCALED_EPSILON) {
            if (unmillable.empty())
                last = offset_ex(last, mill_extra_size);
            else {
                //FIXME only work if mill_extra_size < mill_nozzle/2 (becasue it's the extra offset from unmillable)
                //FIXME overhangs if mill_extra_size is too big
                //FIXME merge with process_arachne?
                ExPolygons growth = diff_ex(offset_ex(last, mill_extra_size), unmillable, ApplySafetyOffset::Yes);
                last.insert(last.end(), growth.begin(), growth.end());
                last = union_ex(last);
            }
        }

        this->throw_if_canceled();
        // Add perimeters on overhangs : initialization
        ExPolygons overhangs_unsupported;
        if ((this->config->extra_perimeters_overhangs || (this->config->overhangs_reverse && this->layer->id() % 2 == 1))
            && !last.empty() && this->lower_slices != NULL && !this->lower_slices->empty()) {
            //remove holes from lower layer, we only ant that for overhangs, not bridges!
            ExPolygons lower_without_holes;
            for (const ExPolygon& exp : *this->lower_slices)
                lower_without_holes.emplace_back(to_expolygon(exp.contour));
            // opening is offset2-+
            overhangs_unsupported = opening_ex(diff_ex(last, lower_without_holes, ApplySafetyOffset::Yes), scale_t(print_config->resolution_internal));
            if (!overhangs_unsupported.empty()) {
                //only consider overhangs and let bridges alone
                //only consider the part that can be bridged (really, by the bridge algorithm)
                //first, separate into islands (ie, each ExPlolygon)
                //only consider the bottom layer that intersect unsupported, to be sure it's only on our island.
                const ExPolygonCollection lower_island(diff_ex(last, overhangs_unsupported));
                ExPolygons bridgeable;
                for (ExPolygon unsupported : overhangs_unsupported) {
                    BridgeDetector detector( unsupported,
                        lower_island.expolygons,
                        this->overhang_flow.scaled_spacing(),
                        scale_t(this->print_config->bridge_precision.get_abs_value(this->overhang_flow.spacing())),
                        this->layer->id());
                    if (detector.detect_angle(Geometry::deg2rad(this->config->bridge_angle.value)))
                        expolygons_append(bridgeable, union_ex(detector.coverage()));
                }
                if (!bridgeable.empty()) {
                    //simplify to avoid most of artefacts from printing lines.
                    ExPolygons bridgeable_simplified;
                    for (const ExPolygon& poly : bridgeable) {
                        poly.simplify(perimeter_spacing / 2, &bridgeable_simplified);
                    }

                    //offset by perimeter spacing because the simplify may have reduced it a bit.
                    if (!bridgeable_simplified.empty()) {
                        bridgeable_simplified = offset_ex(bridgeable_simplified, double(perimeter_spacing));
                        overhangs_unsupported = diff_ex(overhangs_unsupported, bridgeable_simplified, ApplySafetyOffset::Yes);
                    }
                }
            }
        }
        bool has_steep_overhang = false;
        if (this->layer->id() % 2 == 1 && this->config->overhangs_reverse //check if my option is set and good layer
            && !last.empty() && this->lower_slices != NULL && !this->lower_slices->empty() //has something to work with 
            ) {
            ExPolygons overhangs = diff_ex(last, *lower_slices);
            coord_t offset = scale_t(config->overhangs_reverse_threshold.get_abs_value(this->perimeter_flow.width()));
            //version with: scale_(std::tan(PI * (0.5f / 90) * config->overhangs_reverse_threshold.value ) * this->layer->height)

            if (offset_ex(overhangs, -offset / 2.).size() > 0) {
                //allow this loop to be printed in reverse
                has_steep_overhang = true;
            }
        }

        // In case no perimeters are to be generated, contour_count / holes_count will equal to 0.            
        std::vector<PerimeterGeneratorLoops> contours(contour_count);    // depth => loops
        std::vector<PerimeterGeneratorLoops> holes(holes_count);       // depth => loops
        ThickPolylines thin_walls_thickpolys;
        ExPolygons no_last_gapfill;
        // we loop one time more than needed in order to find gaps after the last perimeter was applied
        for (int perimeter_idx = 0;; ++perimeter_idx) {  // outer loop is 0
            this->throw_if_canceled();

            // We can add more perimeters if there are uncovered overhangs
            // improvement for future: find a way to add perimeters only where it's needed.
            bool has_overhang = false;
            if (this->config->extra_perimeters_overhangs && !last.empty() && !overhangs_unsupported.empty()) {
                overhangs_unsupported = intersection_ex(overhangs_unsupported, last, ApplySafetyOffset::Yes);
                if (overhangs_unsupported.size() > 0) {
                    //please don't stop adding perimeter yet.
                    has_overhang = true;
                }
            }

            // allow this perimeter to overlap itself?
            float thin_perimeter = this->config->thin_perimeters.get_abs_value(1);
            if (perimeter_idx > 0 && thin_perimeter != 0) {
                thin_perimeter = this->config->thin_perimeters_all.get_abs_value(1);
            }
            bool allow_perimeter_anti_hysteresis = thin_perimeter >= 0;
            if (thin_perimeter < 0) {
                thin_perimeter = -thin_perimeter;
            }
            if (thin_perimeter < 0.02) { // can create artifacts
                thin_perimeter = 0;
            }

            // Calculate next onion shell of perimeters.
            // this variable stored the next onion
            ExPolygons next_onion;
            // like next_onion, but with all polygons, even ones that didn't grow and so won't be added as periemter
            ExPolygons area_used;
            ExPolygons* all_next_onion = &next_onion;

            if (perimeter_idx == 0) {
                // compute next onion
                    // the minimum thickness of a single loop is:
                    // ext_width/2 + ext_spacing/2 + spacing/2 + width/2
                coordf_t good_spacing    = ext_perimeter_width / 2;
                coordf_t overlap_spacing = (1 - thin_perimeter) * ext_perimeter_spacing / 2;
                if (holes_count == 0 || contour_count == 0) {

                    if (holes_count == 0) {
                        for (ExPolygon &expoly : last) { 
                            last_asynch.push_back(ExPolygonAsynch{ExPolygonAsynch::ExPolygonAsynchType::epatShrinkContour, expoly,
                                // inner_offset                                                 outer_offset (go spacing limit)
                                -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2, -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2,
                                -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2, -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2});
                        }
                        last_asynch_initialized = true;
                        grow_contour_only(last_asynch, get_perimeter_spacing(), 0 /*no overlap for external*/, false /*no round peri for external*/);
                    } else {
                        for (ExPolygon &expoly : last) { 
                            last_asynch.push_back(ExPolygonAsynch{ExPolygonAsynch::ExPolygonAsynchType::epatGrowHole, expoly,
                                // inner_offset                                                 outer_offset (go spacing limit)
                                -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2, -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2,
                                -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2, -coordf_t(this->get_perimeter_width() - get_perimeter_spacing())/2});
                        }
                        last_asynch_initialized = true;
                        grow_holes_only(last_asynch, next_onion, get_perimeter_spacing(), 0 /*no overlap for external*/, false /*no round peri for external*/);
                    }
                } else {
                    if (thin_perimeter > 0.98) {
                        next_onion = offset_ex(last, -(float) (ext_perimeter_width / 2),
                                               ClipperLib::JoinType::jtMiter, 3);
                    } else {
                        coordf_t good_spacing    = ext_perimeter_width / 2;
                        coordf_t overlap_spacing = (1 - thin_perimeter) * ext_perimeter_spacing / 2;
                        next_onion               = offset2_ex(last, -(float) (good_spacing + overlap_spacing - 1),
                                                +(float) (overlap_spacing + 1), ClipperLib::JoinType::jtMiter, 3);
                    }
                    if (thin_perimeter < 0.7) {
                        // offset2_ex can create artifacts, if too big. see superslicer#2428
                        next_onion = intersection_ex(next_onion, offset_ex(last, -(float) (ext_perimeter_width / 2),
                                                                           ClipperLib::JoinType::jtMiter, 3));
                    }
                }
                
                bool special_area = contour_count == 0 || holes_count == 0;
                if (special_area && (this->config->thin_walls || m_spiral_vase)) {
                    area_used = next_onion;
                    for(auto& expolycontainer : last_asynch)
                        area_used.push_back(expolycontainer.expoly);
                    all_next_onion = &area_used;
                }
                // look for thin walls
                if (this->config->thin_walls) {

                    // detect edge case where a curve can be split in multiple small chunks.
                    if (allow_perimeter_anti_hysteresis && !special_area) {
                        std::vector<float> divs = { 2.1f, 1.9f, 2.2f, 1.75f, 1.5f }; //don't go too far, it's not possible to print thin wall after that
                        size_t idx_div = 0;
                        while (next_onion.size() > last.size() && idx_div < divs.size()) {
                            float div = divs[idx_div];
                            //use a sightly bigger spacing to try to drastically improve the split, that can lead to very thick gapfill
                            ExPolygons next_onion_secondTry = offset2_ex(
                                last,
                                -(float)((ext_perimeter_width / 2) + (ext_perimeter_spacing / div) - 1),
                                +(float)((ext_perimeter_spacing / div) - 1));
                            if (next_onion.size() > next_onion_secondTry.size() * 1.2 && next_onion.size() > next_onion_secondTry.size() + 2) {
                                next_onion = next_onion_secondTry;
                            }
                            idx_div++;
                        }
                    }

                    // the following offset2 ensures almost nothing in @thin_walls is narrower than $min_width
                    // (actually, something larger than that still may exist due to mitering or other causes)
                    coord_t min_width = scale_t(this->config->thin_walls_min_width.get_abs_value(this->ext_perimeter_flow.nozzle_diameter()));

                    ExPolygons no_thin_zone = offset_ex(*all_next_onion, double(ext_perimeter_width / 2), jtSquare);
                    // medial axis requires non-overlapping geometry
                    ExPolygons thin_zones = diff_ex(last, no_thin_zone, ApplySafetyOffset::Yes);
                    //don't use offset2_ex, because we don't want to merge the zones that have been separated.
                        //a very little bit of overlap can be created here with other thin polygons, but it's more useful than worisome.
                    ExPolygons half_thins = offset_ex(thin_zones, double(-min_width / 2));
                    //simplify them
                    for (ExPolygon& half_thin : half_thins) {
                        half_thin.remove_point_too_near(ext_perimeter_width/20);
                    }
                    //we push the bits removed and put them into what we will use as our anchor
                    if (half_thins.size() > 0) {
                        no_thin_zone = diff_ex(last, offset_ex(half_thins, double(min_width / 2 - SCALED_EPSILON)), ApplySafetyOffset::Yes);
                    }
                    ExPolygons thins;
                    // compute a bit of overlap to anchor thin walls inside the print.
                    for (ExPolygon& half_thin : half_thins) {
                        //growing back the polygon
                        ExPolygons thin = offset_ex(half_thin, double(min_width / 2));
                        assert(thin.size() <= 1);
                        if (thin.empty()) continue;
                        coord_t thin_walls_overlap = scale_t(this->config->thin_walls_overlap.get_abs_value(this->ext_perimeter_flow.nozzle_diameter()));
                        ExPolygons anchor = intersection_ex(offset_ex(half_thin, double(min_width / 2) +
                            (float)(thin_walls_overlap), jtSquare), no_thin_zone, ApplySafetyOffset::Yes);
                        ExPolygons bounds = union_ex(thin, anchor, ApplySafetyOffset::Yes);
                        for (ExPolygon& bound : bounds) {
                            if (!intersection_ex(thin[0], bound).empty()) {
                                //be sure it's not too small to extrude reliably
                                thin[0].remove_point_too_near(ext_perimeter_width / 10);
                                if (thin[0].area() > min_width * (ext_perimeter_width + ext_perimeter_spacing)) {
                                    thins.push_back(thin[0]);
                                    bound.remove_point_too_near(ext_perimeter_width / 10);
                                    // the maximum thickness of our thin wall area is equal to the minimum thickness of a single loop (*1.2 because of circles approx. and enlrgment from 'div')
                                    Slic3r::Geometry::MedialAxis ma{ thin[0], (coord_t)((ext_perimeter_width + ext_perimeter_spacing) * 1.2),
                                        min_width, scale_t(this->layer->height) };
                                    ma.use_bounds(bound)
                                        .use_min_real_width(scale_t(this->ext_perimeter_flow.nozzle_diameter()))
                                        .use_tapers(thin_walls_overlap)
                                        .set_min_length(ext_perimeter_width + ext_perimeter_spacing)
                                        .build(thin_walls_thickpolys);
                                }
                                break;
                            }
                        }
                    }
                    // use perimeters to extrude area that can't be printed by thin walls
                    // it's a bit like re-add thin area into perimeter area.
                    // it can over-extrude a bit, but it's for a better good.
                    if(!special_area) {
                        if (thin_perimeter > 0.98)
                            next_onion = union_ex(next_onion, offset_ex(diff_ex(last, thins, ApplySafetyOffset::Yes),
                                -(float)(ext_perimeter_width / 2),
                                ClipperLib::JoinType::jtMiter,
                                3));
                        else if (thin_perimeter > 0.01)
                            next_onion = union_ex(next_onion, offset2_ex(diff_ex(last, thins, ApplySafetyOffset::Yes),
                                -(float)((ext_perimeter_width / 2) + ((1 - thin_perimeter) * ext_perimeter_spacing / 4)),
                                (float)((1 - thin_perimeter) * ext_perimeter_spacing / 4),
                                ClipperLib::JoinType::jtMiter,
                                3));
                        else
                            next_onion = union_ex(next_onion, offset2_ex(diff_ex(last, thins, ApplySafetyOffset::Yes),
                                -(float)((ext_perimeter_width / 2) + (ext_perimeter_spacing / 4)),
                                (float)(ext_perimeter_spacing / 4),
                                ClipperLib::JoinType::jtMiter,
                                3));
                        //simplify the loop to avoid almost-0 segments
                        resolution = get_resolution(1, false, &surface);
                        ExPolygons next_onion_temp;
                        for (ExPolygon& exp : next_onion)
                            exp.simplify((resolution < SCALED_EPSILON ? SCALED_EPSILON : resolution), &next_onion_temp);
                        //mask
                        next_onion = intersection_ex(next_onion_temp, last);
                    }
                }
                if (m_spiral_vase && all_next_onion->size() > 1) {
                    assert(contour_count > 0);
                    // Remove all but the largest area polygon.
                    keep_largest_contour_only(*all_next_onion);
                }
            } else {
                //FIXME Is this offset correct if the line width of the inner perimeters differs
                // from the line width of the infill?
                coord_t good_spacing = (perimeter_idx == 1) ? ext_perimeter_spacing2 : perimeter_spacing;
                if (thin_perimeter <= 0.98) {
                    // This path will ensure, that the perimeters do not overfill, as in 
                    // prusa3d/Slic3r GH #32, but with the cost of rounding the perimeters
                    // excessively, creating gaps, which then need to be filled in by the not very
                    // reliable gap fill algorithm.
                    // Also the offset2(perimeter, -x, x) may sometimes lead to a perimeter, which is larger than
                    // the original.
                    next_onion = offset2_ex(last,
                        -(float)(good_spacing + (1 - thin_perimeter) * perimeter_spacing / 2 - 1),
                        +(float)((1 - thin_perimeter) * perimeter_spacing / 2 - 1),
                        (round_peri ? ClipperLib::JoinType::jtRound : ClipperLib::JoinType::jtMiter),
                        (round_peri ? min_round_spacing : 3));
                    if (allow_perimeter_anti_hysteresis) {
                        // now try with different min spacing if we fear some hysteresis
                        // TODO, do that for each polygon from last, instead to do for all of them in one go.
                        ExPolygons no_thin_onion = offset_ex(last, double(-good_spacing));
                        if (last_area < 0) {
                            last_area = 0;
                            for (const ExPolygon &expoly : last) { last_area += expoly.area(); }
                        }
                        double new_area = 0;
                        for (const ExPolygon &expoly : next_onion) { new_area += expoly.area(); }

                        std::vector<float> divs{1.8f, 1.6f}; // don't over-extrude, so don't use divider >2
                        size_t             idx_div = 0;
                        while ((next_onion.size() > no_thin_onion.size() ||
                                (new_area != 0 && last_area > new_area * 100)) &&
                               idx_div < divs.size()) {
                            float div = divs[idx_div];
                            //use a sightly bigger spacing to try to drastically improve the split, that can lead to very thick gapfill
                            ExPolygons next_onion_secondTry = offset2_ex(
                                last,
                                -(float)(good_spacing + (1 - thin_perimeter) * (perimeter_spacing / div) - 1),
                                +(float)((1 - thin_perimeter) * (perimeter_spacing / div) - 1));
                            if (next_onion.size() > next_onion_secondTry.size() * 1.2 && next_onion.size() > next_onion_secondTry.size() + 2) {
                                // don't get it if it creates too many
                                next_onion = next_onion_secondTry;
                            } else if (next_onion.size() > next_onion_secondTry.size() || last_area > new_area * 100) {
                                // don't get it if it's too small
                                double area_new = 0;
                                for (const ExPolygon &expoly : next_onion_secondTry) { area_new += expoly.area(); }
                                if (last_area > area_new * 100 || new_area == 0) {
                                    next_onion = next_onion_secondTry;
                                }
                            }
                            idx_div++;
                        }
                        last_area = new_area;
                    }
                } else {
                    // If "overlapping_perimeters" is enabled, this paths will be entered, which 
                    // leads to overflows, as in prusa3d/Slic3r GH #32
                    next_onion = offset_ex(last, double(-good_spacing),
                        (round_peri ? ClipperLib::JoinType::jtRound : ClipperLib::JoinType::jtMiter),
                        (round_peri ? min_round_spacing : 3));
                }
                
                std::vector<ExPolygonAsynch> *touse;
                std::vector<ExPolygonAsynch> copy;
                if (perimeter_idx < std::max(contour_count, holes_count)) {
                    touse = &last_asynch;
                } else {
                    // for gap fill only: use a copy
                    copy = last_asynch;
                    touse = &copy;
                }
                if (contour_count > perimeter_idx && holes_count <= perimeter_idx) {
                    grow_contour_only(*touse, good_spacing, (1 - thin_perimeter) * perimeter_spacing / 2,
                                        round_peri, min_round_spacing);
                }
                if (holes_count > perimeter_idx && contour_count <= perimeter_idx) {
                    grow_holes_only(*touse, next_onion, good_spacing,
                                    (1 - thin_perimeter) * perimeter_spacing / 2, round_peri, min_round_spacing);
                }

                bool special_area = contour_count == 0 || holes_count == 0;
                if (special_area && (this->config->thin_walls || m_spiral_vase)) {
                    area_used = next_onion;
                    for (auto &expolycontainer : *touse) area_used.push_back(expolycontainer.expoly);
                    all_next_onion = &area_used;
                }


                // look for gaps
                if (this->config->gap_fill_enabled.value
                    //check if we are going to have an other perimeter
                    && (perimeter_idx < std::max(contour_count, holes_count) || has_overhang || all_next_onion->empty() ||
                        (this->config->gap_fill_last.value && perimeter_idx == std::max(contour_count, holes_count)))) {
                    // not using safety offset here would "detect" very narrow gaps
                    // (but still long enough to escape the area threshold) that gap fill
                    // won't be able to fill but we'd still remove from infill area
                    no_last_gapfill = offset_ex(*all_next_onion, 0.5f * good_spacing + 10,
                        (round_peri ? ClipperLib::JoinType::jtRound : ClipperLib::JoinType::jtMiter),
                        (round_peri ? min_round_spacing : 3));
                    if (perimeter_idx == 1) {
                        append(gaps, diff_ex(
                            offset_ex(last, -0.5f * this->get_ext_perimeter_spacing()),
                            no_last_gapfill));  // safety offset
                    } else {
                        append(gaps, diff_ex(
                            offset_ex(last, -0.5f * this->get_perimeter_spacing()),
                            no_last_gapfill));  // safety offset
                    }
                }
            }
            //{
            //    static int aodfjiaz = 0;
            //    std::stringstream stri;
            //    stri << this->layer->id() << "_perimeter_loop_" << (aodfjiaz++) << ".svg";
            //    SVG svg(stri.str());
            //    svg.draw(surface.expolygon, "grey");
            //    svg.draw(to_polylines(last), "yellow");
            //    svg.draw(to_polylines(next_onion), "green");
            //    svg.Close();
            //}

            if (next_onion.empty() && last_asynch.empty()) {
                // Store the number of loops actually generated.
                if (perimeter_idx < contour_count) {
                    assert(contours.size() == contour_count);
                    for(size_t i = perimeter_idx; i<contours.size(); i++)
                        assert(contours[perimeter_idx].empty());
                    contour_count = perimeter_idx;
                    contours.resize(contour_count);
                }
                if (perimeter_idx < holes_count) {
                    assert(holes.size() == holes_count);
                    for(size_t i = perimeter_idx; i<holes.size(); i++)
                        assert(holes[perimeter_idx].empty());
                    holes_count = perimeter_idx;
                    holes.resize(holes_count);
                }
                // No region left to be filled in.
                last.clear();
                break;
            } else if (perimeter_idx >= std::max(contour_count, holes_count)) {
                if (has_overhang) {
                    contour_count++;
                    holes_count++; //TODO: only increase the ones that are needed (or just use 2.7)
                    contours.emplace_back();
                    holes.emplace_back();
                } else {
                    // If perimeter_idx > loop_number, we were looking just for gaps.
                    break;
                }
            }
            if (contour_count <= perimeter_idx && !next_onion.empty()) {
                assert(contour_count <= perimeter_idx);
                assert(holes_count > perimeter_idx);
                //assert(contours.size() == perimeter_idx);
                contour_count = perimeter_idx + 1;
                while (contours.size() < contour_count) {
                    contours.emplace_back();
                }
            }
            
            assert(contours.size() == contour_count);
            assert(holes.size() == holes_count);

            // fuzzify params
            const bool fuzzify_contours = this->config->fuzzy_skin != FuzzySkinType::None && perimeter_idx == 0 && this->layer->id() > 0;
            const bool fuzzify_holes = this->config->fuzzy_skin == FuzzySkinType::Shell && perimeter_idx == 0 && this->layer->id() > 0 ;
            const bool fuzzify_all = this->config->fuzzy_skin == FuzzySkinType::All && this->layer->id() > 0 ;

            //push last_asynch or next_onion into contours & holes
            if (!last_asynch.empty()) {
                // we already put the last hole, now add contours.
                for (auto &exp : last_asynch) {
                    if (exp.type == ExPolygonAsynch::epatShrinkContour) {
                        assert(next_onion.empty());
                        assert(exp.expoly.contour.is_counter_clockwise());
                        if (exp.expoly.contour.length() > SCALED_EPSILON) // TODO: atleastLength
                            contours[perimeter_idx].emplace_back(exp.expoly.contour, perimeter_idx, true,
                                                                 has_steep_overhang, fuzzify_contours || fuzzify_all);
                    } else {
                        // we already put the last contour, now add holes
                        // contours from hole collapse is added via next_onion
                        assert(exp.type == ExPolygonAsynch::epatGrowHole);
                        for (auto &hole : exp.expoly.holes) {
                            assert(hole.is_clockwise());
                            if(hole.length() > SCALED_EPSILON) // TODO: atleastLength
                                holes[perimeter_idx].emplace_back(hole, perimeter_idx, false, has_steep_overhang,
                                                              fuzzify_contours || fuzzify_all);
                        }
                    }
                }
            }
            for (const ExPolygon& expolygon : next_onion) {
                //TODO: add width here to allow variable width (if we want to extrude a sightly bigger perimeter, see thin wall)
                if(contour_count > perimeter_idx && expolygon.contour.length() > SCALED_EPSILON) // TODO: atleastLength
                    contours[perimeter_idx].emplace_back(expolygon.contour, perimeter_idx, true, has_steep_overhang, fuzzify_contours || fuzzify_all);
                if (!expolygon.holes.empty() && holes_count > perimeter_idx) {
                    holes[perimeter_idx].reserve(holes[perimeter_idx].size() + expolygon.holes.size());
                    for (const Polygon& hole : expolygon.holes)
                        if(hole.length() > SCALED_EPSILON) // TODO: atleastLength
                            holes[perimeter_idx].emplace_back(hole, perimeter_idx, false, has_steep_overhang, fuzzify_holes || fuzzify_all);
                }
            }

            //simplify the loop to avoid artifacts when shrinking almost-0 segments
            resolution = get_resolution(perimeter_idx + 1, false, &surface);
            last.clear();
            for (ExPolygon &exp : next_onion)
                exp.simplify((resolution < SCALED_EPSILON ? SCALED_EPSILON : resolution), &last);

            // store surface for top infill if only_one_perimeter_top
            if (perimeter_idx == 0 && (config->only_one_perimeter_top && this->upper_slices != NULL)
                && contour_count > 1 && holes_count > 1) {
                ExPolygons next;
                split_top_surfaces(this->lower_slices, this->upper_slices, last, results.top_fills, next,
                                   results.fill_clip);
                last = next;
            }
            
            //if next turn we are in asynch mode, move from last to last_asynch
            if ( !last_asynch_initialized && (
                (holes_count == perimeter_idx + 1 && contour_count > perimeter_idx + 1) ||
                (contour_count == perimeter_idx + 1 && holes_count > perimeter_idx + 1))) {
                coordf_t last_spacing = perimeter_idx == 0 ? get_ext_perimeter_spacing() / 2 :
                                                             get_perimeter_spacing() / 2;
                // populate last_asynch from last
                for (ExPolygon &expoly : last) {
                    last_asynch.push_back(
                        {holes_count == perimeter_idx + 1 ? ExPolygonAsynch::ExPolygonAsynchType::epatShrinkContour :
                                                            ExPolygonAsynch::ExPolygonAsynchType::epatGrowHole,
                         std::move(expoly), -last_spacing, last_spacing, -last_spacing, last_spacing});
                }
                last.clear();
                last_asynch_initialized = true;
            }
        }

        // fuzzify
        const bool fuzzify_gapfill = this->config->fuzzy_skin == FuzzySkinType::All && this->layer->id() > 0;
        // check for extracting extra perimeters from gapfill
        if (!gaps.empty()) {
            // if needed, add it to the first empty contour list
            const size_t contours_size = contour_count;
            assert(contours.size() == contour_count);
            //first, find loops and try to extract a perimeter from them.
            for (size_t gap_idx = 0; gap_idx < gaps.size(); gap_idx++) {
                ExPolygon& expoly = gaps[gap_idx];
                if (!expoly.holes.empty()) {
                    //this is a a sort of a loop
                    //try to see if it's possible to add a "perimeter"
                    ExPolygons contour_expolygon = offset_ex(expoly, -(float)(perimeter_spacing / 2), ClipperLib::jtMiter, 3);
                    if (contour_expolygon.size() == 1 && !contour_expolygon.front().holes.empty()) {
                        //OK
                        // update list & variable to let the new perimeter be taken into account
                        contour_count = contours_size + 1;
                        if (contours_size >= contours.size()) {
                            contours.emplace_back();
                            holes.emplace_back();
                        }
                        assert(contours.size() == contour_count);
                        //Add the new perimeter
                        contours[contours_size].emplace_back(contour_expolygon.front().contour, contours_size, true, has_steep_overhang, fuzzify_gapfill);
                        //create the new gapfills
                        ExPolygons gapfill_area = offset_ex(Polygons{ expoly.contour }, -(float)(perimeter_spacing));
                        ExPolygons to_add = intersection_ex(ExPolygons{ expoly }, gapfill_area);
                        //add the new gapfill
                        if (to_add.size() == 0)
                            expoly.clear();
                        else
                            expoly = to_add.front();
                        for (size_t j = 1; j < to_add.size(); j++)
                            gaps.push_back(to_add[j]);
                    }
                }
            }
        }
        assert(contours.size() == contour_count);
        //assert(holes.size() == holes_count);
        // nest loops: holes first
        for (int d = 0; d < holes_count; ++d) {
            PerimeterGeneratorLoops& holes_d = holes[d];
            // loop through all holes having depth == d
            for (int hole_idx = 0; hole_idx < (int)holes_d.size(); ++hole_idx) {
                const PerimeterGeneratorLoop& loop = holes_d[hole_idx];
                assert(loop.polygon.length() > SCALED_EPSILON);
                // find the hole loop that contains this one, if any
                for (int t = d + 1; t < holes_count; ++t) {
                    for (int j = 0; j < (int)holes[t].size(); ++j) {
                        PerimeterGeneratorLoop& candidate_parent = holes[t][j];
                        if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                            candidate_parent.children.push_back(loop);
                            holes_d.erase(holes_d.begin() + hole_idx);
                            --hole_idx;
                            goto NEXT_LOOP;
                        }
                    }
                }
                // if no hole contains this hole, find the contour loop that contains it
                for (int t = contours.size() - 1; t >= 0; --t) {
                    for (int j = 0; j < (int)contours[t].size(); ++j) {
                        PerimeterGeneratorLoop& candidate_parent = contours[t][j];
                        if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                            candidate_parent.children.push_back(loop);
                            holes_d.erase(holes_d.begin() + hole_idx);
                            --hole_idx;
                            goto NEXT_LOOP;
                        }
                    }
                }
                // no perimeter, then add the hole like a perimeter.
                while(d >= contours.size())
                    contours.emplace_back();
                contours[d].push_back(loop);
                holes_d.erase(holes_d.begin() + hole_idx);
                --hole_idx;
            NEXT_LOOP:;
            }
        }
        // nest contour loops
        for (int d = contours.size() - 1; d >= 1; --d) {
            PerimeterGeneratorLoops& contours_d = contours[d];
            // loop through all contours having depth == d
            for (int contour_idx = 0; contour_idx < (int)contours_d.size(); ++contour_idx) {
                const PerimeterGeneratorLoop& loop = contours_d[contour_idx];
                assert(loop.polygon.length() > SCALED_EPSILON);
                // find the contour loop that contains it
                for (int t = d - 1; t >= 0; --t) {
                    for (size_t j = 0; j < contours[t].size(); ++j) {
                        PerimeterGeneratorLoop& candidate_parent = contours[t][j];
                        if (candidate_parent.polygon.contains(loop.polygon.first_point())) {
                            candidate_parent.children.push_back(loop);
                            contours_d.erase(contours_d.begin() + contour_idx);
                            --contour_idx;
                            goto NEXT_CONTOUR;
                        }
                    }
                }
                //can't find one, put in front
                if (contours.front().empty()) {
                    contours.front().push_back(loop);
                } else {
                    contours.front().front().children.push_back(loop);
                }
                contours_d.erase(contours_d.begin() + contour_idx);
                --contour_idx;
            NEXT_CONTOUR:;
            }
        }
        //remove all empty perimeters
        while(contours.size() > 1 && contours.front().empty())
            contours.erase(contours.begin());
        // fuse all unfused 
        // at this point, all loops should be in contours[0] (= contours.front() )
        // or no perimeters nor holes have been generated, too small area.

        assert(contours.empty() || contours.front().size() >= 1);
        // collection of loops to add into loops
        ExtrusionEntityCollection peri_entities;
        if (!contours.empty()) {
            if (config->perimeter_loop.value) {
                // onlyone_perimeter = >fusion all perimeterLoops
                for (PerimeterGeneratorLoop &loop : contours.front()) {
                    ExtrusionLoop extr_loop = this->_traverse_and_join_loops(loop, get_all_Childs(loop),
                                                                             loop.polygon.points.front());
                    // ExtrusionLoop extr_loop = this->_traverse_and_join_loops_old(loop, loop.polygon.points.front(), true);
                    if (extr_loop.paths.back().polyline.back() != extr_loop.paths.front().polyline.front()) {
                        extr_loop.paths.back().polyline.append(extr_loop.paths.front().polyline.front());
                        assert(false);
                    }
                    peri_entities.append(extr_loop);
                }

                // append thin walls
                if (!thin_walls_thickpolys.empty()) {
                    if (this->object_config->thin_walls_merge) {
                        _merge_thin_walls(peri_entities, thin_walls_thickpolys);
                    } else {
                        peri_entities.append(
                            Geometry::thin_variable_width(thin_walls_thickpolys, erThinWall, this->ext_perimeter_flow,
                                                          std::max(ext_perimeter_width / 4,
                                                                   scale_t(this->print_config->resolution)),
                                                          false));
                    }
                    thin_walls_thickpolys.clear();
                }
            } else {
                if (this->object_config->thin_walls_merge) {
                    ThickPolylines no_thin_walls;
                    peri_entities = this->_traverse_loops(contours.front(), no_thin_walls);
                    _merge_thin_walls(peri_entities, thin_walls_thickpolys);
                } else {
                    peri_entities = this->_traverse_loops(contours.front(), thin_walls_thickpolys);
                }
            }
        }
#if _DEBUG
        peri_entities.visit(LoopAssertVisitor{});
#endif
        // remove the un-needed top collection if only one child.
        //peri_entities.visit(CollectionSimplifyVisitor{});
        if (peri_entities.entities().size() == 1) {
            if (ExtrusionEntityCollection *coll_child = dynamic_cast<ExtrusionEntityCollection *>(
                    peri_entities.set_entities().front());
                coll_child != nullptr) {
                peri_entities.set_can_sort_reverse(coll_child->can_sort(), coll_child->can_reverse());
                peri_entities.append_move_from(*coll_child);
                peri_entities.remove(0);
            }
        }


        //{
        //    static int aodfjiaqsdz = 0;
        //    std::stringstream stri;
        //    stri << this->layer->id() << "_perimeter_loops_" << (aodfjiaqsdz++) << ".svg";
        //    SVG svg(stri.str());
        //    svg.draw(surface.expolygon, "grey");
        //    struct TempVisitor : public ExtrusionVisitorRecursiveConst {
        //        SVG* svg;
        //        virtual void use(const ExtrusionPath& path) override { svg->draw(path.polyline, "green"); }
        //    } bbvisitor;
        //    bbvisitor.svg = &svg;
        //    peri_entities.visit(bbvisitor);
        //    svg.Close();
        //}

        // append perimeters for this slice as a collection
        if (!peri_entities.empty()) {
            //move it, to avoid to clone evrything and then delete it
            this->loops->append(peri_entities);
        }
    } // for each loop of an island
#if _DEBUG
    this->loops->visit(LoopAssertVisitor{});
#endif

    // fill gaps
    ExPolygons gaps_ex;
    if (!gaps.empty()) {
        // collapse 
        coordf_t min = 0.2 * perimeter_width * (1 - INSET_OVERLAP_TOLERANCE);
        //be sure we don't gapfill where the perimeters are already touching each other (negative spacing).
        min = std::max(min, double(Flow::new_from_spacing((float)EPSILON, (float)this->perimeter_flow.nozzle_diameter(), (float)this->layer->height, (float)this->perimeter_flow.spacing_ratio(), false).scaled_width()));
        coordf_t real_max = 2.5 * perimeter_spacing;
        const coordf_t minwidth = scale_d(this->config->gap_fill_min_width.get_abs_value(unscaled((double)perimeter_width)));
        const coordf_t maxwidth = scale_d(this->config->gap_fill_max_width.get_abs_value(unscaled((double)perimeter_width)));
        const coord_t minlength = scale_t(this->config->gap_fill_min_length.get_abs_value(unscaled((double)perimeter_width)));
        if (minwidth > 0) {
            min = std::max(min, minwidth);
        }
        coordf_t max = real_max;
        if (maxwidth > 0) {
            max = std::min(max, maxwidth);
        }
        const coord_t gapfill_extension = scale_t(this->config->gap_fill_extension.get_abs_value(unscaled((double)perimeter_width)));
        //remove areas that are too big (shouldn't occur...)
        ExPolygons too_big = offset2_ex(gaps, double(-max / 2), double(+max / 2));
        ExPolygons gaps_ex_to_test = too_big.empty() ? gaps : diff_ex(gaps, too_big, ApplySafetyOffset::Yes);
        const double minarea = scale_d(scale_d(this->config->gap_fill_min_area.get_abs_value(unscaled((double)perimeter_width) * unscaled((double)perimeter_width))));
        // check each gapfill area to see if it's printable.
        for (const ExPolygon& expoly : gaps_ex_to_test) {
            this->throw_if_canceled();
            //remove too small gaps that are too hard to fill.
            //ie one that are smaller than an extrusion with width of min and a length of max.
            if (expoly.area() > minarea) {
                const coordf_t offset_test = min * 0.5;
                ExPolygons expoly_after_shrink_test = offset_ex(ExPolygons{expoly}, -offset_test);
                //if the shrink split the area in multipe bits
                if (expoly_after_shrink_test.size() > 1) {
                    //remove too small bits
                    for (int exp_idx = 0; exp_idx < expoly_after_shrink_test.size(); exp_idx++) {
                        if (expoly_after_shrink_test[exp_idx].area() < (SCALED_EPSILON * SCALED_EPSILON * 4)) {
                            expoly_after_shrink_test.erase(expoly_after_shrink_test.begin() + exp_idx);
                            exp_idx--;
                        } else {
                            ExPolygons wider = offset_ex(ExPolygons{ expoly_after_shrink_test[exp_idx] }, offset_test);
                            if (wider.empty() || wider[0].area() < minarea) {
                                expoly_after_shrink_test.erase(expoly_after_shrink_test.begin() + exp_idx);
                                exp_idx--;
                            }
                        }
                    }
                    //maybe some areas are a just bit too thin, try with just a little more offset to remove them.
                    const coordf_t offset_test_2 = min * 0.8;
                    ExPolygons expoly_after_shrink_test2 = offset_ex(ExPolygons{expoly}, -offset_test_2);
                    for (int exp_idx = 0; exp_idx < expoly_after_shrink_test2.size(); exp_idx++) {
                        if (expoly_after_shrink_test2[exp_idx].area() < (SCALED_EPSILON * SCALED_EPSILON * 4)) {
                            expoly_after_shrink_test2.erase(expoly_after_shrink_test2.begin() + exp_idx);
                            exp_idx--;

                        } else {
                            ExPolygons wider = offset_ex(ExPolygons{ expoly_after_shrink_test2[exp_idx] }, offset_test_2);
                            if (wider.empty() || wider[0].area() < minarea) {
                                expoly_after_shrink_test2.erase(expoly_after_shrink_test2.begin() + exp_idx);
                                exp_idx--;
                            }
                        }
                    }
                    //it's better if there are significantly less extrusions
                    if (expoly_after_shrink_test.size() / 1.42 > expoly_after_shrink_test2.size()) {
                        expoly_after_shrink_test2 = offset_ex(expoly_after_shrink_test2, offset_test_2);
                        //insert with move instead of copy
                        std::move(expoly_after_shrink_test2.begin(), expoly_after_shrink_test2.end(), std::back_inserter(gaps_ex));
                    } else {
                        expoly_after_shrink_test = offset_ex(expoly_after_shrink_test, offset_test);
                        std::move(expoly_after_shrink_test.begin(), expoly_after_shrink_test.end(), std::back_inserter(gaps_ex));
                    }
                } else {
                    expoly_after_shrink_test = offset_ex(expoly_after_shrink_test, offset_test);
                    std::move(expoly_after_shrink_test.begin(), expoly_after_shrink_test.end(), std::back_inserter(gaps_ex));
                }
            }
        }
        // create lines from the area
        ThickPolylines polylines;
        for (const ExPolygon& ex : gaps_ex) {
            Geometry::MedialAxis md{ ex, coord_t(real_max), coord_t(min), coord_t(this->layer->height) };
            if (minlength > 0) {
                md.set_min_length(minlength);
            }
            if (gapfill_extension > 0) {
                md.set_extension_length(gapfill_extension);
            }
            md.set_biggest_width(max);
            md.build(polylines);
        }
        // create extrusion from lines
        Flow gap_fill_flow = Flow::new_from_width(this->perimeter_flow.width(), this->perimeter_flow.nozzle_diameter(),this->perimeter_flow.height(), this->config->gap_fill_overlap.get_abs_value(1.), false);
        if (!polylines.empty()) {
            this->gap_fill->append(Geometry::thin_variable_width(
                polylines,
                erGapFill, 
                gap_fill_flow, 
                scale_t(this->print_config->resolution_internal),
                true));
            /*  Make sure we don't infill narrow parts that are already gap-filled
                (we only consider this surface's gaps to reduce the diff() complexity).
                Growing actual extrusions ensures that gaps not filled by medial axis
                are not subtracted from fill surfaces (they might be too short gaps
                that medial axis skips but infill might join with other infill regions
                and use zigzag).  */
                // get clean surface of gap
            results.gap_srf = union_ex(offset(this->gap_fill->polygons_covered_by_width(float(SCALED_EPSILON) / 10), float(SCALED_EPSILON / 2)));
            // intersection to ignore the bits of gapfill tha may be over infill, as it's epsilon and there may be some voids here anyway.
            results.gap_srf = intersection_ex(results.gap_srf, gaps_ex);
            // the diff(last, gap) will be done after, as we have to keep the last un-gapped to avoid unneeded gap/infill offset
        }
    }

    if (contour_count == 0 && holes_count == 0) {
        // for the infill shell, move it a little bit inside so the extrusion tip don't go over the sides.
        results.inner_perimeter = offset_ex(last, -(this->get_perimeter_width() - get_perimeter_spacing()) / 2);
    } else {
        coordf_t last_spacing = std::max(contour_count, holes_count) == 1 ?
                                                      get_ext_perimeter_spacing() / 2 :
                                                      get_perimeter_spacing() / 2;
        results.inner_perimeter = offset_ex(last, -last_spacing);
        if (!last_asynch.empty()) {
            // merge with last_async
            for (auto &exp : last_asynch) {
                if (exp.offset_contour_inner == exp.offset_holes_inner) {
                    append(results.inner_perimeter, offset_ex(exp.expoly, exp.offset_contour_inner));
                } else {
                    // offset contour & holes separatly
                    // first holes:
                    assert(exp.offset_holes_inner <= 0);
                    Polygons holes = offset(get_holes_as_contour(exp.expoly), -exp.offset_holes_inner);
                    // we are growing (fake) perimeter, so it can creates holes.
                    for (size_t i = 0; i < holes.size(); ++i) {
                        Polygon &fakeperi = holes[i];
                        if (fakeperi.is_clockwise()) {
                            // put real perimeters in results.inner_perimeter
                            fakeperi.make_counter_clockwise();
                            results.inner_perimeter.push_back(ExPolygon(fakeperi));
                            holes.erase(holes.begin() + i);
                            i--;
                        }
                    }
                    // now shrink perimeter
                    Polygons perimeters = offset(exp.expoly.contour, exp.offset_contour_inner);
                    // as it shrink, it can creates more perimeter, not a big deal.
                    for (auto &p : perimeters) assert(p.is_counter_clockwise());
                    
                    // now diff and add
                    append(results.inner_perimeter, diff_ex(perimeters, holes));
                }
            }
        }
    }

    return results;
}


ExtrusionPaths PerimeterGenerator::create_overhangs(const Polyline& loop_polygons, ExtrusionRole role, bool is_external) const {
    ExtrusionPaths paths;
    const double overhangs_width = this->config->overhangs_width.get_abs_value(this->overhang_flow.nozzle_diameter());
    const double overhangs_width_speed = this->config->overhangs_width_speed.get_abs_value(this->overhang_flow.nozzle_diameter());
    if ( 0 == overhangs_width && 0 == overhangs_width_speed) {
        //error
        ExtrusionPath path(role, false);
        path.polyline = loop_polygons;
        path.mm3_per_mm = is_external ? this->ext_mm3_per_mm()           : this->mm3_per_mm();
        path.width = is_external ?      this->ext_perimeter_flow.width() : this->perimeter_flow.width();
        path.height = (float)this->layer->height;
        assert(path.mm3_per_mm == path.mm3_per_mm);
        assert(path.width == path.width);
        assert(path.height == path.height);
        paths.push_back(path);
        return paths;
    
    }
    //set the fan & speed before the flow
    Polylines ok_polylines = { loop_polygons };

    Polylines small_speed;
    Polylines big_speed;
    bool no_small_flow = _lower_slices_bridge_speed_big == _lower_slices_bridge_flow_small;
    Polylines small_flow;
    Polylines big_flow;
#ifdef _DEBUG
    for (Polyline& poly : ok_polylines)
        for (int i = 0; i < poly.points.size() - 1; i++)
            assert(poly.points[i] != poly.points[i + 1]);
#endif

    Polylines* previous = &ok_polylines;
    if (overhangs_width_speed > 0 && (overhangs_width_speed < overhangs_width || overhangs_width == 0)) {
        if (!this->_lower_slices_bridge_speed_small.empty()) {
            small_speed = diff_pl(*previous, this->_lower_slices_bridge_speed_small);
#ifdef _DEBUG
            for (Polyline& poly : small_speed) //                       assert small_speed
                for (int i = 0; i < poly.points.size() - 1; i++) //     assert small_speed
                    assert(poly.points[i] != poly.points[i + 1]); //    assert small_speed
#endif
            if (!small_speed.empty()) {
                *previous = intersection_pl(*previous, this->_lower_slices_bridge_speed_small);
#ifdef _DEBUG
                for (Polyline& poly : *previous) //                         assert previous
                    for (int i = 0; i < poly.points.size() - 1; i++) //     assert previous
                        assert(poly.points[i] != poly.points[i + 1]); //    assert previous
#endif
                previous = &small_speed;
            }
        }
        if (!this->_lower_slices_bridge_speed_big.empty()) {
            big_speed = diff_pl(*previous, this->_lower_slices_bridge_speed_big);
#ifdef _DEBUG
            for (Polyline& poly : big_speed) //                         assert big_speed
                for (int i = 0; i < poly.points.size() - 1; i++) //     assert big_speed
                    assert(poly.points[i] != poly.points[i + 1]); //    assert big_speed
#endif
            if (!big_speed.empty()) {
                *previous = intersection_pl(*previous, this->_lower_slices_bridge_speed_big);
#ifdef _DEBUG
                for (Polyline& poly : *previous) //                         assert previous
                    for (int i = 0; i < poly.points.size() - 1; i++) //     assert previous
                        assert(poly.points[i] != poly.points[i + 1]); //    assert previous
#endif
                previous = &big_speed;
            }
        }
    }
    if (overhangs_width > 0) {
        if (!this->_lower_slices_bridge_flow_small.empty()) {
            small_flow = diff_pl(*previous, this->_lower_slices_bridge_flow_small);
#ifdef _DEBUG
            for (Polyline& poly : small_flow) //                        assert small_flow
                for (int i = 0; i < poly.points.size() - 1; i++) //     assert small_flow
                    assert(poly.points[i] != poly.points[i + 1]); //    assert small_flow
#endif
            if (!small_flow.empty()) {
                *previous = intersection_pl(*previous, this->_lower_slices_bridge_flow_small);
#ifdef _DEBUG
                for (Polyline& poly : *previous) //                         assert previous
                    for (int i = 0; i < poly.points.size() - 1; i++) //     assert previous
                        assert(poly.points[i] != poly.points[i + 1]); //    assert previous
#endif
                previous = &small_flow;
            }
        }
        if (!this->_lower_slices_bridge_flow_big.empty()) {
            big_flow = diff_pl(*previous, this->_lower_slices_bridge_flow_big);
#ifdef _DEBUG
            for (Polyline& poly : big_flow) //                          assert big_flow
                for (int i = 0; i < poly.points.size() - 1; i++) //     assert big_flow
                    assert(poly.points[i] != poly.points[i + 1]); //    assert big_flow
#endif
            if (!big_flow.empty()) {
                *previous = intersection_pl(*previous, this->_lower_slices_bridge_flow_big);
#ifdef _DEBUG
                for (Polyline& poly : *previous) //                         assert previous
                    for (int i = 0; i < poly.points.size() - 1; i++) //     assert previous
                        assert(poly.points[i] != poly.points[i + 1]); //    assert previous
#endif
                previous = &big_flow;
            }
        }
    }

    //note: layer height is used to identify the path type
    if (!ok_polylines.empty()) {
        //fast track
        if (small_speed.empty() && big_speed.empty() && small_flow.empty() && big_flow.empty()) {
            ExtrusionPath path(role, false);
            path.polyline = loop_polygons;
            path.mm3_per_mm = is_external ? this->ext_mm3_per_mm() : this->mm3_per_mm();
            path.width = is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width();
            path.height = (float)this->layer->height;
            return { path };
        }
        extrusion_paths_append(
            paths,
            ok_polylines,
            role,
            is_external ? this->ext_mm3_per_mm() : this->mm3_per_mm(),
            is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width(),
            0);
    }
    if (!small_speed.empty()) {
        extrusion_paths_append(
            paths,
            small_speed,
            erOverhangPerimeter,
            is_external ? this->ext_mm3_per_mm() : this->mm3_per_mm(),
            is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width(),
            no_small_flow ? 2 : 1);
    }
    if (!big_speed.empty()) {
        extrusion_paths_append(
            paths,
            big_speed,
            erOverhangPerimeter,
            is_external ? this->ext_mm3_per_mm() : this->mm3_per_mm(),
            is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width(),
            no_small_flow ? 3 : 2);
    }
    if (!small_flow.empty()) {
        extrusion_paths_append(
            paths,
            small_flow,
            erOverhangPerimeter,
            this->mm3_per_mm_overhang(),
            this->overhang_flow.width(),
            3);
    }
    if (!big_flow.empty()) {
        extrusion_paths_append(
            paths,
            big_flow,
            erOverhangPerimeter,
            this->mm3_per_mm_overhang(),
            this->overhang_flow.width(),
            4);
    }

    // reapply the nearest point search for starting point
    // We allow polyline reversal because Clipper may have randomly reversed polylines during clipping.
    if(!paths.empty())
        chain_and_reorder_extrusion_paths(paths, &paths.front().first_point());

    bool has_normal = !ok_polylines.empty();
    bool has_speed = !small_speed.empty() || !big_speed.empty();
    bool has_flow = !small_flow.empty() || !big_flow.empty();

    std::function<void(ExtrusionPaths&, const std::function<bool(ExtrusionPath&, ExtrusionPath&, ExtrusionPath&)>&)> foreach = [](ExtrusionPaths &paths, const std::function<bool(ExtrusionPath&, ExtrusionPath&, ExtrusionPath&)>& doforeach) {
        if (paths.size() > 2)
            for (int i = 1; i < paths.size() - 1; i++) {
                if (doforeach(paths[i - 1], paths[i], paths[i + 1])) {
                    paths.erase(paths.begin() + i);
                    i--;
                    if (paths[i].height == paths[i + 1].height) {
                        paths[i].polyline.append(paths[i + 1].polyline);
                        paths.erase(paths.begin() + i + 1);
                    }
                }
            }
        if (paths.size() > 2)
            if (doforeach(paths[paths.size() - 2], paths.back(), paths.front())) {
                paths.erase(paths.end() - 1);
                if (paths.back().height == paths.front().height) {
                    paths.back().polyline.append(paths.front().polyline);
                    paths.front().polyline.swap(paths.back().polyline);
                    paths.erase(paths.end() - 1);
                }
            }
        if (paths.size() > 2)
            if (doforeach(paths.back(), paths.front(), paths[1])) {
                paths.erase(paths.begin());
                if (paths.back().height == paths.front().height) {
                    paths.back().polyline.append(paths.front().polyline);
                    paths.front().polyline.swap(paths.back().polyline);
                    paths.erase(paths.end() - 1);
                }
            }
    };

    if (paths.size() > 2) {
        double min_length = this->perimeter_flow.scaled_width() * 2;
        double ok_length = this->perimeter_flow.scaled_width() * 20;

        //curr will be deleted by 'foreach' (our caller, see above) if the return value is true. So its points need to be merged in prev or next.
        foreach(paths, [min_length, ok_length](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            if (curr.length() < min_length) {
                float diff_height = std::abs(prev.height - curr.height) - std::abs(next.height - curr.height);
                //have to choose the rigth path
                if (diff_height < 0 || (diff_height == 0 && prev.length() > next.length())) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                return true;
            } else if(((int)curr.height) % 2 == 1 && curr.length() > ok_length){
                curr.height++;
                if (prev.height == curr.height) {
                    prev.polyline.append(curr.polyline);
                    return true;
                } else if (next.height == curr.height) {
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        });

        foreach(paths, [](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            if (curr.height == 3) {
                //have to choose the rigth path
                if (prev.height == 4 || (prev.height == 2 && next.height < 2)) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                return true;
            }
            return false;
        });
        foreach(paths, [](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            if (curr.height == 1) {
                //have to choose the rigth path
                if (prev.height == 2 || (prev.height == 0 && next.height > 2)) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                return true;
            }
            return false;
        });
    }
    if(paths.size() == 2){
        double min_length = this->perimeter_flow.scaled_width() * 2;
        if (paths.front().length() < min_length) {
            paths.front().polyline.append(paths.back().polyline);
            paths.back().polyline.swap(paths.front().polyline);
            paths.erase(paths.begin());
        }else if (paths.back().length() < min_length) {
            paths.front().polyline.append(paths.back().polyline);
            paths.erase(paths.begin() + 1);
        }
    }
    //set correct height
    for (ExtrusionPath& path : paths) {
        path.height = path.height < 3 ? (float)this->layer->height : this->overhang_flow.height();
    }

    return paths;
}

#ifdef _DEBUG
void test_overhangs(const std::vector<ClipperLib_Z::Path>& path1, const std::vector<ClipperLib_Z::Path>& path2, Points &outer_points)
{
    for (const ClipperLib_Z::Path &poly : path1)
        for (int i = 0; i < poly.size() - 1; i++)
            assert(poly[i] != poly[i + 1]);
    for (const ClipperLib_Z::Path &poly : path2)
        for (int i = 0; i < poly.size() - 1; i++)
            assert(poly[i] != poly[i + 1]);
    // check if points are equal
    //Points path2_points;
    //Points path1_points;
    //for (auto &line : path2)
    //    for (auto &pt : line) path2_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
    //for (auto &line : path1)
    //    for (auto &pt : line) path1_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
    //for (Point &pt : path2_points) {
    //    bool found    = false;
    //    bool in_outer = false;
    //    Point pt2_almost;
    //    for (Point &pt2 : path1_points)
    //        if (pt.coincides_with_epsilon(pt2)) {
    //            found = true;
    //            pt2_almost = pt2;
    //            break;
    //        }
    //    for (Point &pt2 : outer_points)
    //        if (pt.coincides_with_epsilon(pt2)) {
    //            found    = true;
    //            in_outer = true;
    //            pt2_almost = pt2;
    //            break;
    //        }
    //    assert(found);
    //    found    = false;
    //    in_outer = false;
    //    for (Point &pt2 : path1_points)
    //        if (pt.coincides_with(pt2)) {
    //            found = true;
    //            break;
    //        }
    //    for (Point &pt2 : outer_points)
    //        if (pt.coincides_with(pt2)) {
    //            found    = true;
    //            in_outer = true;
    //            break;
    //        }
    //    assert(found);
    //}
    //  => points can be different from diff & intersect
    // TODO: create a new operation that create the diff & intersect at the same time
}
#endif

//TODO: transform to ExtrusionMultiPath instead of ExtrusionPaths
ExtrusionPaths PerimeterGenerator::create_overhangs(const ClipperLib_Z::Path& arachne_path, ExtrusionRole role, bool is_external) const {
    ExtrusionPaths paths;
    const bool is_loop = Point{ arachne_path.front().x(), arachne_path.front().y() }.coincides_with_epsilon(Point{ arachne_path.back().x(), arachne_path.back().y() });
    const double overhangs_width = this->config->overhangs_width.get_abs_value(this->overhang_flow.nozzle_diameter());
    const double overhangs_width_speed = this->config->overhangs_width_speed.get_abs_value(this->overhang_flow.nozzle_diameter());
    if (0 == overhangs_width && 0 == overhangs_width_speed) {
        //error
        //assert(path.mm3_per_mm == path.mm3_per_mm);
        //assert(path.width == path.width);
        //assert(path.height == path.height);
        append(paths, Geometry::unsafe_variable_width(Arachne::to_thick_polyline(arachne_path),
            role,
            is_external ? this->ext_perimeter_flow : this->perimeter_flow,
            std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
            (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10));
        //(const ThickPolyline& polyline, const ExtrusionRole role, const Flow& flow, const coord_t resolution_internal, const coord_t tolerance)
        for (ExtrusionPath& path : paths) {
            //these variable_width paths aren't gapfill, they are proper perimeters
            path.set_can_reverse(false);
        }
        return paths;

    }
    //set the fan & speed before the flow
    ClipperLib_Z::Paths ok_polylines = { arachne_path };

    ClipperLib_Z::Paths small_speed;
    ClipperLib_Z::Paths big_speed;
    bool no_small_flow = _lower_slices_bridge_speed_big_clipperpaths == _lower_slices_bridge_flow_small_clipperpaths;
    ClipperLib_Z::Paths small_flow;
    ClipperLib_Z::Paths big_flow;
#ifdef _DEBUG
    for (ClipperLib_Z::Path& poly : ok_polylines)
        for (int i = 0; i < poly.size() - 1; i++)
            assert(poly[i] != poly[i + 1]);
#endif

    std::vector<ClipperLib_Z::Path>* previous = &ok_polylines;
    if (overhangs_width_speed > 0 && (overhangs_width_speed < overhangs_width || overhangs_width == 0)) {
        if (!this->_lower_slices_bridge_speed_small_clipperpaths.empty()) {
            //small_speed = diff_pl(*previous, this->_lower_slices_bridge_speed_small);
#ifdef _DEBUG
            Points outer_points;
            for(auto & line: *previous) for(auto &pt : line) outer_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
#endif
            small_speed = clip_extrusion(*previous, this->_lower_slices_bridge_speed_small_clipperpaths, ClipperLib_Z::ctDifference);
#ifdef _DEBUG
            for (ClipperLib_Z::Path& poly : small_speed) //                       assert small_speed
                for (int i = 0; i < poly.size() - 1; i++) //     assert small_speed
                    assert(poly[i] != poly[i + 1]); //    assert small_speed
#endif
            if (!small_speed.empty()) {
                *previous = clip_extrusion(*previous, this->_lower_slices_bridge_speed_small_clipperpaths, ClipperLib_Z::ctIntersection);
#ifdef _DEBUG
                test_overhangs(small_speed, *previous, outer_points);
                test_overhangs(*previous, small_speed, outer_points);
#endif
                previous = &small_speed;
            }
        }
        if (!this->_lower_slices_bridge_speed_big.empty()) {
#ifdef _DEBUG
            Points outer_points;
            for(auto & line: *previous) for(auto &pt : line) outer_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
#endif
            big_speed = clip_extrusion(*previous, this->_lower_slices_bridge_speed_big_clipperpaths, ClipperLib_Z::ctDifference);
#ifdef _DEBUG
            for (ClipperLib_Z::Path& poly : big_speed) //                         assert big_speed
                for (int i = 0; i < poly.size() - 1; i++) //     assert big_speed
                    assert(poly[i] != poly[i + 1]); //    assert big_speed
#endif
            if (!big_speed.empty()) {
                *previous = clip_extrusion(*previous, this->_lower_slices_bridge_speed_big_clipperpaths, ClipperLib_Z::ctIntersection);
#ifdef _DEBUG
                test_overhangs(big_speed, *previous, outer_points);
                test_overhangs(*previous, big_speed, outer_points);
#endif
                previous = &big_speed;
            }
        }
    }
    if (overhangs_width > 0) {
        if (!this->_lower_slices_bridge_flow_small.empty()) {
#ifdef _DEBUG
            Points outer_points;
            for(auto & line: *previous) for(auto &pt : line) outer_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
#endif
            small_flow = clip_extrusion(*previous, this->_lower_slices_bridge_flow_small_clipperpaths, ClipperLib_Z::ctDifference);
#ifdef _DEBUG
            for (ClipperLib_Z::Path& poly : small_flow) //                        assert small_flow
                for (int i = 0; i < poly.size() - 1; i++) //     assert small_flow
                    assert(poly[i] != poly[i + 1]); //    assert small_flow
#endif
            if (!small_flow.empty()) {
                *previous = clip_extrusion(*previous, this->_lower_slices_bridge_flow_small_clipperpaths, ClipperLib_Z::ctIntersection);
#ifdef _DEBUG
                test_overhangs(small_flow, *previous, outer_points);
                test_overhangs(*previous, small_flow, outer_points);
#endif
                previous = &small_flow;
            }
        }
        if (!this->_lower_slices_bridge_flow_big.empty()) {
#ifdef _DEBUG
            Points outer_points;
            for(auto & line: *previous) for(auto &pt : line) outer_points.emplace_back(coord_t(pt.x()), coord_t(pt.y()));
#endif
            big_flow = clip_extrusion(*previous, this->_lower_slices_bridge_flow_big_clipperpaths, ClipperLib_Z::ctDifference);
#ifdef _DEBUG
            for (ClipperLib_Z::Path& poly : big_flow) //                          assert big_flow
                for (int i = 0; i < poly.size() - 1; i++) //     assert big_flow
                    assert(poly[i] != poly[i + 1]); //    assert big_flow
#endif
            if (!big_flow.empty()) {
                *previous = clip_extrusion(*previous, this->_lower_slices_bridge_flow_big_clipperpaths, ClipperLib_Z::ctIntersection);
#ifdef _DEBUG
                test_overhangs(big_flow, *previous, outer_points);
                test_overhangs(*previous, big_flow, outer_points);
#endif
                previous = &big_flow;
            }
        }
    }

    //note: layer height is used to identify the path type
    if (!ok_polylines.empty()) {
        //fast track
        if (small_speed.empty() && big_speed.empty() && small_flow.empty() && big_flow.empty()) {
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(Arachne::to_thick_polyline(arachne_path),
                    role,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            // thickpaths can be empty if extrusion_path is too short
            assert(thickpaths.empty() || thickpaths.front().first_point().x() == arachne_path.front().x());
            assert(thickpaths.empty() || thickpaths.front().first_point().y() == arachne_path.front().y());
            assert(thickpaths.empty() || thickpaths.back().last_point().x() == arachne_path.back().x());
            assert(thickpaths.empty() || thickpaths.back().last_point().y() == arachne_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                path.set_can_reverse(!is_loop);
                paths.push_back(std::move(path));
            }
            return paths;
        }
        for (const ClipperLib_Z::Path& extrusion_path : ok_polylines) {
            auto thick_poly = Arachne::to_thick_polyline(extrusion_path);
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(thick_poly,
                    role,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            // thickpaths can be empty if extrusion_path is too short
            assert(thickpaths.empty() || thickpaths.front().first_point().x() == extrusion_path.front().x());
            assert(thickpaths.empty() || thickpaths.front().first_point().y() == extrusion_path.front().y());
            assert(thickpaths.empty() || thickpaths.back().last_point().x() == extrusion_path.back().x());
            assert(thickpaths.empty() || thickpaths.back().last_point().y() == extrusion_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                path.set_can_reverse(!is_loop);
                path.height = 0;
                paths.push_back(std::move(path));
            }
        }
    }
    if (!small_speed.empty()) {
        for (const ClipperLib_Z::Path& extrusion_path : small_speed) {
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(Arachne::to_thick_polyline(extrusion_path),
                    erOverhangPerimeter,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            // thickpaths can be empty if extrusion_path is too short
            assert(thickpaths.empty() || thickpaths.front().first_point().x() == extrusion_path.front().x());
            assert(thickpaths.empty() || thickpaths.front().first_point().y() == extrusion_path.front().y());
            assert(thickpaths.empty() || thickpaths.back().last_point().x() == extrusion_path.back().x());
            assert(thickpaths.empty() || thickpaths.back().last_point().y() == extrusion_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                path.set_can_reverse(!is_loop);
                path.height = no_small_flow ? 2 : 1;
                paths.push_back(std::move(path));
            }
        }
    }
    if (!big_speed.empty()) {
        for (const ClipperLib_Z::Path& extrusion_path : big_speed) {
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(Arachne::to_thick_polyline(extrusion_path),
                    erOverhangPerimeter,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            // thickpaths can be empty if extrusion_path is too short
            assert(thickpaths.empty() || thickpaths.front().first_point().x() == extrusion_path.front().x());
            assert(thickpaths.empty() || thickpaths.front().first_point().y() == extrusion_path.front().y());
            assert(thickpaths.empty() || thickpaths.back().last_point().x() == extrusion_path.back().x());
            assert(thickpaths.empty() || thickpaths.back().last_point().y() == extrusion_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                path.set_can_reverse(!is_loop);
                path.height = no_small_flow ? 3 : 2;
                paths.push_back(std::move(path));
            }
        }
    }
    if (!small_flow.empty()) {
        for (const ClipperLib_Z::Path& extrusion_path : small_flow) {
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(Arachne::to_thick_polyline(extrusion_path),
                    erOverhangPerimeter,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            // thickpaths can be empty if extrusion_path is too short
            assert(thickpaths.empty() || thickpaths.front().first_point().x() == extrusion_path.front().x());
            assert(thickpaths.empty() || thickpaths.front().first_point().y() == extrusion_path.front().y());
            assert(thickpaths.empty() || thickpaths.back().last_point().x() == extrusion_path.back().x());
            assert(thickpaths.empty() || thickpaths.back().last_point().y() == extrusion_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                // change flow to overhang one if too much.
                if (path.mm3_per_mm > this->overhang_flow.mm3_per_mm() ){
                    path.mm3_per_mm = this->overhang_flow.mm3_per_mm();
                    path.height = this->overhang_flow.height();
                    path.width = this->overhang_flow.width();
                }
                path.set_can_reverse(!is_loop);
                path.height = 3;
                paths.push_back(std::move(path));
            }
        }
    }
    if (!big_flow.empty()) {
        for (const ClipperLib_Z::Path& extrusion_path : big_flow) {
            ExtrusionPaths thickpaths = Geometry::unsafe_variable_width(Arachne::to_thick_polyline(extrusion_path),
                    erOverhangPerimeter,
                    is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                    std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                    (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10);
            if (thickpaths.empty()) {
                // Note: can create problem with chain_and_reorder_extrusion_paths
                assert(extrusion_path.size() < 2 ||
                       Point(extrusion_path.front().x(), extrusion_path.front().y())
                           .coincides_with_epsilon(Point(extrusion_path.back().x(), extrusion_path.back().y())));
                continue;
            }
#ifdef _DEBUG
            for (int i = 1; i < thickpaths.size(); i++) {
                assert(thickpaths[i - 1].last_point() == thickpaths[i].first_point());
            }
#endif
            assert(thickpaths.front().first_point().x() == extrusion_path.front().x());
            assert(thickpaths.front().first_point().y() == extrusion_path.front().y());
            assert(thickpaths.back().last_point().x() == extrusion_path.back().x());
            assert(thickpaths.back().last_point().y() == extrusion_path.back().y());
            for (ExtrusionPath& path : thickpaths) {
                // change flow to overhang one if too much.
                if (path.mm3_per_mm > this->overhang_flow.mm3_per_mm()) {
                    path.mm3_per_mm = this->overhang_flow.mm3_per_mm();
                    path.height = this->overhang_flow.height();
                    path.width = this->overhang_flow.width();
                }
                path.set_can_reverse(!is_loop);
                path.height = 4;
                paths.push_back(std::move(path));
            }
        }
    }

    //FIXME from here, it's exactly the same as the other create_overhangs, please merge that into a function.
    
    //(or not)
    Point first_point(arachne_path.front().x(), arachne_path.front().y());
    // reapply the nearest point search for starting point
    // We allow polyline reversal because Clipper may have randomly reversed polylines during clipping.
    if (!paths.empty())
        chain_and_reorder_extrusion_paths(paths, &first_point);

    //check if evvrything is okay (it can fail)
    bool not_sorted_enough = false;
    for (int i = 1; i < paths.size(); i++) {
        if (!paths[i - 1].last_point().coincides_with_epsilon(paths[i].first_point())) {
            not_sorted_enough = true;
            break;
        }
    }
    if (not_sorted_enough) {
        Point other_point = paths[1].first_point();
        chain_and_reorder_extrusion_paths(paths, &other_point);
        auto path = paths.back();
        paths.erase(paths.end()-1);
        paths.insert(paths.begin(), path);
        bool not_sorted_enough = false;
        for (int i = 1; i < paths.size(); i++) {
            if (paths[i - 1].last_point().coincides_with_epsilon(paths[i].first_point())) {
                not_sorted_enough = true;
                break;
            }
        }
        if (not_sorted_enough) {
            // do it manually by brute-force
            // TODO
            chain_and_reorder_extrusion_paths(paths, &first_point);
        }
    }

    for (int i = 1; i < paths.size(); i++) {
        // diff/inter can generate points with ~3-5 unit of diff.
        if (paths[i - 1].last_point() != paths[i].first_point()) {
            assert(paths[i - 1].last_point().coincides_with_epsilon(paths[i].first_point()));
            Point middle = (paths[i - 1].last_point() + paths[i].first_point()) / 2;
            paths[i - 1].polyline.set_points().back() = middle;
            paths[i].polyline.set_points().front() = middle;
        }
    }
#ifdef _DEBUG
    for (int i = 1; i < paths.size(); i++) {
        assert(paths[i - 1].last_point() == paths[i].first_point());
    }
#endif

    bool has_normal = !ok_polylines.empty();
    bool has_speed = !small_speed.empty() || !big_speed.empty();
    bool has_flow = !small_flow.empty() || !big_flow.empty();

    std::function<void(ExtrusionPaths&, const std::function<bool(ExtrusionPath&, ExtrusionPath&, ExtrusionPath&)>&)> foreach = [is_loop](ExtrusionPaths& paths, const std::function<bool(ExtrusionPath&, ExtrusionPath&, ExtrusionPath&)>& doforeach) {
        if (paths.size() > 2)
            for (int i = 1; i < paths.size() - 1; i++) {
                if (doforeach(paths[i - 1], paths[i], paths[i + 1])) {
                    paths.erase(paths.begin() + i);
                    i--;
                    if (paths[i].height == paths[i + 1].height) {
                        paths[i].polyline.append(paths[i + 1].polyline);
                        paths.erase(paths.begin() + i + 1);
                    }
                }
            }
        if (is_loop) {
            if (paths.size() > 2) {
                if (doforeach(paths[paths.size() - 2], paths.back(), paths.front())) {
                    paths.erase(paths.end() - 1);
                    if (paths.back().height == paths.front().height) {
                        paths.back().polyline.append(paths.front().polyline);
                        paths.front().polyline.swap(paths.back().polyline);
                        paths.erase(paths.end() - 1);
                    }
                }
            }
            if (paths.size() > 2) {
                if (doforeach(paths.back(), paths.front(), paths[1])) {
                    paths.erase(paths.begin());
                    if (paths.back().height == paths.front().height) {
                        paths.back().polyline.append(paths.front().polyline);
                        paths.front().polyline.swap(paths.back().polyline);
                        paths.erase(paths.end() - 1);
                    }
                }
            }
        }
    };

    if (paths.size() > 2) {
        double min_length = this->perimeter_flow.scaled_width() * 2;
        double ok_length = this->perimeter_flow.scaled_width() * 20;
        foreach(paths, [min_length, ok_length](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            assert(prev.last_point() == curr.first_point());
            if (curr.length() < min_length) {
                // if too small, merge it with prev or next, whatever is best.
                float diff_height = std::abs(prev.height - curr.height) - std::abs(next.height - curr.height);
                //have to choose the rigth path
                if (diff_height < 0 || (diff_height == 0 && prev.length() > next.length())) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                assert(prev.last_point() == next.first_point());
                return true;
            } else if (((int)curr.height) % 2 == 1 && curr.length() > ok_length) {
                // not too small, and we are still undecided (%2==1) of the type of this section
                // so choose to add it to the prev, or next, or keep it in limbo until another foreach.
                curr.height++;
                if (prev.height == curr.height) {
                    //merge to prev
                    prev.polyline.append(curr.polyline);
                    return true;
                } else if (next.height == curr.height) {
                    //merge to next
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                    return true;
                } else {
                    //keep it in limbo
                    return false;
                }
            }
            return false;
            });

        foreach(paths, [](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            if (curr.height == 3) {
                //have to choose the rigth path
                if (prev.height == 4 || (prev.height == 2 && next.height < 2)) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                return true;
            }
            return false;
            });
        foreach(paths, [](ExtrusionPath& prev, ExtrusionPath& curr, ExtrusionPath& next) {
            if (curr.height == 1) {
                //have to choose the rigth path
                if (prev.height == 2 || (prev.height == 0 && next.height > 2)) {
                    //merge to previous
                    assert(prev.last_point() == curr.first_point());
                    assert(curr.polyline.size() > 1);
                    prev.polyline.append(curr.polyline);
                } else {
                    //merge to next
                    assert(curr.last_point() == next.first_point());
                    assert(curr.polyline.size() > 1);
                    curr.polyline.append(next.polyline);
                    next.polyline.swap(curr.polyline);
                }
                return true;
            }
            return false;
            });
    }
    if (paths.size() == 2) {
        double min_length = this->perimeter_flow.scaled_width() * 2;
        if (paths.front().length() < min_length) {
            paths.front().polyline.append(paths.back().polyline);
            paths.back().polyline.swap(paths.front().polyline);
            paths.erase(paths.begin());
        } else if (paths.back().length() < min_length) {
            paths.front().polyline.append(paths.back().polyline);
            paths.erase(paths.begin() + 1);
        }
    }
    //set correct height
    for (ExtrusionPath& path : paths) {
        path.height = path.height < 3 ? (float)this->layer->height : this->overhang_flow.height();
    }

    return paths;
}

// Thanks Cura developers for this function.
static void fuzzy_paths(ExtrusionPaths& paths, coordf_t fuzzy_skin_thickness, coordf_t fuzzy_skin_point_dist)
{
    const coordf_t min_dist_between_points = fuzzy_skin_point_dist * 3. / 4.; // hardcoded: the point distance may vary between 3/4 and 5/4 the supplied value
    const coordf_t range_random_point_dist = fuzzy_skin_point_dist / 2.;
    coordf_t dist_next_point = //min_dist_between_points / 4 + (coordf_t(rand()) * range_random_point_dist / double(RAND_MAX)); // the distance to be traversed on the line before making the first new point
        coordf_t(rand()) * (min_dist_between_points / 2) / double(RAND_MAX); // the distance to be traversed on the line before making the first new point

    // check if the paths length is enough for at least 3 points, or return.
    {
        coordf_t min_dist = min_dist_between_points * 3;
        for (const ExtrusionPath& path : paths) {
            min_dist -= path.length();
            if (min_dist < 0)
                break;
        }
        if (min_dist > 0) {
            // Too small, can't fuzzy.
            return;
        }
    }

    const Point last_point = paths.back().last_point();
    //not always a loop, with arachne
    bool is_loop = paths.front().first_point() == last_point;
#ifdef _DEBUG
    if (is_loop)
        assert(paths.back().last_point() == paths.front().first_point());
    for (int i = 1; i < paths.size(); i++) {
        assert(paths[i - 1].last_point() == paths[i].first_point());
    }
#endif
    Point p0 = is_loop ? last_point : paths.front().first_point();
    const Point* previous_point = is_loop ? &last_point : &paths.front().first_point();
    for (size_t idx_path = 0; idx_path < paths.size(); idx_path++) {
        ExtrusionPath& path = paths[idx_path];
        Points out;
        size_t next_idx = 1;
        assert(path.size() > 1);
        // it always follow
        assert(p0 == path.polyline.front());
        out.reserve(path.polyline.size());
        out.push_back(*previous_point);
        for (; next_idx < path.polyline.size(); next_idx++)
        {
            const Point& p1 = path.polyline.get_points()[next_idx];
            // 'a' is the (next) new point between p0 and p1
            Vec2d  p0p1 = (p1 - p0).cast<double>();
            coordf_t p0p1_size = p0p1.norm();
            //skip points too close to each other.
            if (dist_next_point < p0p1_size) {
                coordf_t p0pa_dist;
                for (p0pa_dist = dist_next_point; p0pa_dist < p0p1_size;
                    p0pa_dist += min_dist_between_points + coordf_t(rand()) * range_random_point_dist / double(RAND_MAX))
                {
                    coordf_t r = coordf_t(rand()) * (fuzzy_skin_thickness * 2.) / double(RAND_MAX) - fuzzy_skin_thickness;
                    out.emplace_back(p0 + (p0p1 * (p0pa_dist / p0p1_size) + perp(p0p1).cast<double>().normalized() * r).cast<coord_t>());
                }
                dist_next_point = p0pa_dist - p0p1_size;
                p0 = p1;
            } else {
                dist_next_point -= p0p1_size;
            }
        }
        if (out.size() <= 1) {
            //too small, erase
            path.polyline.clear();
            paths.erase(paths.begin() + idx_path);
            idx_path--;
        } else {
            p0 = path.polyline.back();
            path.polyline = out;
            previous_point = &path.polyline.back();
        }
    }
    assert(!paths.empty());
    if (is_loop) {
        assert(paths.front().polyline.front() != paths.back().polyline.back());
        //the first point is the old one. remove it and try to make another point if needed.
        if (paths.front().size() > 2 && fuzzy_skin_point_dist * 2 > paths.back().last_point().distance_to(paths.front().polyline.get_points()[1])) {
            //distance small enough and enough points to delete the first, just erase
            paths.front().polyline.clip_first_point();
        }//TODO: else
        //loop -> last point is the same as the first
        paths.back().polyline.append(paths.front().polyline.front());
        assert(paths.front().polyline.front() == paths.back().polyline.back());
    } else {
        //line -> ensure you end with the same last point
        paths.back().polyline.append(last_point);
    }
#ifdef _DEBUG
    if (is_loop)
        assert(paths.back().last_point() == paths.front().first_point());
    for (int i = 1; i < paths.size(); i++) {
        assert(paths[i - 1].last_point() == paths[i].first_point());
    }
#endif
}

ExtrusionEntityCollection PerimeterGenerator::_traverse_loops(
    const PerimeterGeneratorLoops &loops, ThickPolylines &thin_walls, int count_since_overhang /*= -1*/) const
{
    // loops is an arrayref of ::Loop objects
    // turn each one into an ExtrusionLoop object
    ExtrusionEntitiesPtr coll;
    for (const PerimeterGeneratorLoop &loop : loops) {
        bool is_external = loop.is_external();
        
        ExtrusionRole role;
        ExtrusionLoopRole loop_role = ExtrusionLoopRole::elrDefault;
        role = is_external ? erExternalPerimeter : erPerimeter;
        if (loop.is_internal_contour()) {
            // Note that we set loop role to ContourInternalPerimeter
            // also when loop is both internal and external (i.e.
            // there's only one contour loop).
            loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrInternal);
        }
        if (!loop.is_contour) {
            loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrHole);
        }
        if (loop.children.empty()) {
            loop_role = ExtrusionLoopRole(loop_role | ExtrusionLoopRole::elrFirstLoop);
        }
        if (this->config->external_perimeters_vase.value && this->config->external_perimeters_first.value && is_external) {
            if ((loop.is_contour && this->config->external_perimeters_nothole.value) || (!loop.is_contour && this->config->external_perimeters_hole.value)) {
                loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrVase);
            }
        }
        
        // detect overhanging/bridging perimeters
        ExtrusionPaths paths;

        bool can_overhang = this->config->overhangs_width_speed.value > 0
            && this->layer->id() > object_config->raft_layers;
        if(this->object_config->support_material && this->object_config->support_material_contact_distance_type.value == zdNone)
            can_overhang = false;
        if (can_overhang) {
            paths = this->create_overhangs(loop.polygon.split_at_first_point(), role, is_external);
        } else {
            ExtrusionPath path(role, false);
            path.polyline   = loop.polygon.split_at_first_point();
            path.mm3_per_mm = is_external ? this->ext_mm3_per_mm()           : this->mm3_per_mm();
            path.width      = is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width();
            path.height     = (float) this->layer->height;
            assert(path.mm3_per_mm == path.mm3_per_mm);
            assert(path.width == path.width);
            assert(path.height == path.height);
            paths.push_back(path);
        }
        if (loop.fuzzify) {
            double nozle_diameter = is_external ? this->ext_perimeter_flow.nozzle_diameter() : this->perimeter_flow.nozzle_diameter();
            double fuzzy_skin_thickness = config->fuzzy_skin_thickness.get_abs_value(nozle_diameter);
            double fuzzy_skin_point_dist = config->fuzzy_skin_point_dist.get_abs_value(nozle_diameter);
            fuzzy_paths(paths, scale_d(fuzzy_skin_thickness), scale_d(fuzzy_skin_point_dist));
        }

        coll.push_back(new ExtrusionLoop(paths, loop_role));
    }
    
    // append thin walls to the nearest-neighbor search (only for first iteration)
    if (!thin_walls.empty()) {
        ExtrusionEntitiesPtr tw = Geometry::thin_variable_width(thin_walls, erThinWall, this->ext_perimeter_flow, std::max(ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)), false);
        coll.insert(coll.end(), tw.begin(), tw.end());
        //don't add again
        thin_walls.clear();
    }
    // traverse children and build the final collection
    Point zero_point(0, 0);
    //result is  [idx, needReverse] ?
    std::vector<std::pair<size_t, bool>> chain = chain_extrusion_entities(coll, &zero_point);
    ExtrusionEntityCollection coll_out;
    if (chain.empty()) return coll_out;

    //little check: if you have external holes with only one extrusion and internal things, please draw the internal first, just in case it can help print the hole better.
    std::vector<std::pair<size_t, bool>> better_chain;
    {
        std::vector<std::pair<size_t, bool>> alone_holes;
        std::vector<std::pair<size_t, bool>> keep_ordering;
        std::vector<std::pair<size_t, bool>> thin_walls;
        for (const std::pair<size_t, bool> &idx : chain) {
            if (idx.first < loops.size())
                if (!loops[idx.first].is_external() ||
                    (!loops[idx.first].is_contour && !loops[idx.first].children.empty()))
                    alone_holes.push_back(idx);
                else
                    keep_ordering.push_back(idx);
            else
                thin_walls.push_back(idx);
        }
        append(better_chain, std::move(alone_holes));
        append(better_chain, std::move(keep_ordering));
        append(better_chain, std::move(thin_walls));
    }
    assert(better_chain.size() == chain.size());
    
    // if brim will be printed, reverse the order of perimeters so that
    // we continue inwards after having finished the brim
    const bool reverse_contour  = (this->layer->id() == 0 && this->object_config->brim_width.value > 0) ||
                           (this->config->external_perimeters_first.value && this->config->external_perimeters_nothole.value);
    const bool reverse_hole = (this->layer->id() == 0 && this->object_config->brim_width_interior.value > 0) || 
                           (this->config->external_perimeters_first.value && this->config->external_perimeters_hole.value);

    //move from coll to coll_out and getting children of each in the same time. (deep first)
    for (const std::pair<size_t, bool> &idx : better_chain) {
        
        if (idx.first >= loops.size()) {
            // this is a thin wall
            // let's get it from the sorted collection as it might have been reversed
            coll_out.set_entities().reserve(coll_out.entities().size() + 1);
            coll_out.set_entities().emplace_back(coll[idx.first]);
            coll[idx.first] = nullptr;
            if (idx.second)
                coll_out.entities().back()->reverse();
            //if thin extrusion is a loop, make it ccw like a normal contour.
            if (ExtrusionLoop* loop = dynamic_cast<ExtrusionLoop*>(coll_out.entities().back())) {
                loop->make_counter_clockwise();
            }
        } else {
            const PerimeterGeneratorLoop &loop = loops[idx.first];
            ExtrusionLoop* eloop = static_cast<ExtrusionLoop*>(coll[idx.first]);
            bool has_overhang = false;
            if (this->config->overhangs_speed_enforce.value > 0) {
                for (const ExtrusionPath& path : eloop->paths) {
                    if (path.role() == erOverhangPerimeter) {
                        has_overhang = true;
                        break;
                    }
                }
                if (has_overhang || ( count_since_overhang>=0 && this->config->overhangs_speed_enforce.value > count_since_overhang)) {
                    //enforce
                    for (ExtrusionPath& path : eloop->paths) {
                        if (path.role() == erPerimeter || path.role() == erExternalPerimeter) {
                            path.set_role(erOverhangPerimeter);
                        }
                    }
                }

            }
            assert(thin_walls.empty());
            ExtrusionEntityCollection children = this->_traverse_loops(loop.children, thin_walls, has_overhang ? 1 : count_since_overhang < 0 ? -1 : (count_since_overhang+1));
            coll[idx.first] = nullptr;
            bool has_steep_overhangs_this_loop = false;
            if (loop.is_steep_overhang && this->layer->id() % 2 == 1 && !config->perimeter_reverse) {
                has_steep_overhangs_this_loop = HasRoleVisitor::search(*eloop, HasThisRoleVisitor{erOverhangPerimeter});
            }
            if ((loop.is_contour && !reverse_contour) || (!loop.is_contour && reverse_hole)) {
                //note: this->layer->id() % 2 == 1 already taken into account in the is_steep_overhang compute (to save time).
                // if contour: reverse if steep_overhang & odd. if hole: the opposite
                bool clockwise = ((config->perimeter_reverse || has_steep_overhangs_this_loop) && this->layer->id() % 2 == 1) == loop.is_contour;
                if (clockwise)
                    eloop->make_clockwise();
                else
                    eloop->make_counter_clockwise();
                //ensure that our children are printed before us
                if (!children.empty()) {
                    ExtrusionEntityCollection print_child_beforeplz;
                    print_child_beforeplz.set_can_sort_reverse(false, false);
                    if (children.entities().size() > 1 && (children.can_reverse() || children.can_sort())) {
                        print_child_beforeplz.append(children);
                    } else {
                        print_child_beforeplz.append_move_from(children);
                    }
                    print_child_beforeplz.append(*eloop);
                    coll_out.append(std::move(print_child_beforeplz));
                } else {
                    coll_out.append(*eloop);
                }
            } else {
                bool counter_clockwise = ((config->perimeter_reverse || has_steep_overhangs_this_loop) && this->layer->id() % 2 == 1) != loop.is_contour;
                // if hole: reverse if steep_overhang & odd. if contour: the opposite
                if (counter_clockwise)
                    eloop->make_counter_clockwise();
                else
                    eloop->make_clockwise();
                // ensure that our children are printed after us
                if (!children.empty()) {
                    ExtrusionEntityCollection print_child_beforeplz;
                    print_child_beforeplz.set_can_sort_reverse(false, false);
                    print_child_beforeplz.append(*eloop);
                    if (children.entities().size() > 1 && (children.can_reverse() || children.can_sort())) {
                        print_child_beforeplz.append(children);
                    } else {
                        print_child_beforeplz.append_move_from(children);
                    }
                    coll_out.append(std::move(print_child_beforeplz));
                } else {
                    coll_out.append(*eloop);
                }
            }
        }
    }
    return coll_out;
}

ExtrusionEntityCollection PerimeterGenerator::_traverse_extrusions(std::vector<PerimeterGeneratorArachneExtrusion>& pg_extrusions)
{
    ExtrusionEntityCollection extrusion_coll;
    size_t biggest_inset_idx = 0;
    for (PerimeterGeneratorArachneExtrusion& pg_extrusion : pg_extrusions) {
        biggest_inset_idx = std::max(biggest_inset_idx, pg_extrusion.extrusion->inset_idx);
    }
    for (PerimeterGeneratorArachneExtrusion& pg_extrusion : pg_extrusions) {
        Arachne::ExtrusionLine* extrusion = pg_extrusion.extrusion;
        if (extrusion->isZeroLength())
            continue;

        const bool    is_external = extrusion->inset_idx == 0;
        ExtrusionLoopRole loop_role = ExtrusionLoopRole::elrDefault;
        ExtrusionRole role = is_external ? erExternalPerimeter : erPerimeter;
        if (biggest_inset_idx == extrusion->inset_idx) {
            // Note that we set loop role to ContourInternalPerimeter
            // also when loop is both internal and external (i.e.
            // there's only one contour loop).
            loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrInternal | ExtrusionLoopRole::elrFirstLoop);
        }
        if (!pg_extrusion.is_contour) {
            loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrHole);
        }
        if (this->config->external_perimeters_vase.value && this->config->external_perimeters_first.value && is_external) {
            if ((pg_extrusion.is_contour && this->config->external_perimeters_nothole.value) || (!pg_extrusion.is_contour && this->config->external_perimeters_hole.value)) {
                loop_role = (ExtrusionLoopRole)(loop_role | ExtrusionLoopRole::elrVase);
            }
        }

        // fuzzy_extrusion_line() don't work. I can use fuzzy_paths() anyway, not a big deal.
        //if (pg_extrusion.fuzzify)
        //    fuzzy_extrusion_line(*extrusion, scale_d(this->config->fuzzy_skin_thickness.value), scale_d(this->config->fuzzy_skin_point_dist.value));

        ExtrusionPaths paths;
        // detect overhanging/bridging perimeters
        if (this->config->overhangs_width_speed > 0 && this->layer->id() > this->object_config->raft_layers
            && !((this->object_config->support_material || this->object_config->support_material_enforce_layers > 0) &&
                this->object_config->support_material_contact_distance.value == 0)) {
            ClipperLib_Z::Path extrusion_path;
            extrusion_path.reserve(extrusion->size());
            for (const Arachne::ExtrusionJunction& ej : extrusion->junctions) {
                //remove duplicate points from arachne
                if(extrusion_path.empty() || 
                    (ej.p.x() != extrusion_path.back().x() || ej.p.y() != extrusion_path.back().y()))
                    extrusion_path.emplace_back(ej.p.x(), ej.p.y(), ej.w);
            }
            if(extrusion->is_closed)
                assert((extrusion_path.front() - extrusion_path.back()).norm() < SCALED_EPSILON);
            paths = this->create_overhangs(extrusion_path, role, is_external);
            
            // Reapply the nearest point search for starting point.
            // We allow polyline reversal because Clipper may have randomly reversed polylines during clipping.
            // Arachne sometimes creates extrusion with zero-length (just two same endpoints);
            if (!paths.empty()) {
                Point start_point = paths.front().first_point();
                if (!extrusion->is_closed) {
                    // Especially for open extrusion, we need to select a starting point that is at the start
                    // or the end of the extrusions to make one continuous line. Also, we prefer a non-overhang
                    // starting point.
                    struct PointInfo
                    {
                        size_t occurrence  = 0;
                        bool   is_overhang = false;
                    };
                    std::unordered_map<Point, PointInfo, PointHash> point_occurrence;
                    for (const ExtrusionPath &path : paths) {
                        ++point_occurrence[path.first_point()].occurrence;
                        ++point_occurrence[path.last_point()].occurrence;
                        if (path.role() == erOverhangPerimeter) {
                            point_occurrence[path.first_point()].is_overhang = true;
                            point_occurrence[path.last_point()].is_overhang  = true;
                        }
                    }

                    // Prefer non-overhang point as a starting point.
                    for (const std::pair<Point, PointInfo> pt : point_occurrence)
                        if (pt.second.occurrence == 1) {
                            start_point = pt.first;
                            if (!pt.second.is_overhang) {
                                start_point = pt.first;
                                break;
                            }
                        }
                }

                chain_and_reorder_extrusion_paths(paths, &start_point);
                for(size_t i = 1; i< paths.size(); ++i)
                    assert(paths[i-1].last_point().coincides_with_epsilon(paths[i].first_point()));
                if(extrusion->is_closed)
                    assert(paths.back().last_point().coincides_with_epsilon(paths.front().first_point()));
            }
        } else {
            append(paths, Geometry::unsafe_variable_width(Arachne::to_thick_polyline(*extrusion),
                role,
                is_external ? this->ext_perimeter_flow : this->perimeter_flow,
                std::max(this->ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)),
                (is_external ? this->ext_perimeter_flow : this->perimeter_flow).scaled_width() / 10));
        }

        // Apply fuzzify
        if (!paths.empty() && pg_extrusion.fuzzify) {
            double nozle_diameter = is_external ? this->ext_perimeter_flow.nozzle_diameter() : this->perimeter_flow.nozzle_diameter();
            double fuzzy_skin_thickness = config->fuzzy_skin_thickness.get_abs_value(nozle_diameter);
            double fuzzy_skin_point_dist = config->fuzzy_skin_point_dist.get_abs_value(nozle_diameter);
           fuzzy_paths(paths, scale_d(fuzzy_skin_thickness), scale_d(fuzzy_skin_point_dist));
        }

        //set to overhang speed if any chunk is overhang
        bool has_overhang = false;
        if (this->config->overhangs_speed_enforce.value > 0) {
            for (const ExtrusionPath& path : paths) {
                if (path.role() == erOverhangPerimeter) {
                    has_overhang = true;
                    break;
                }
            }
            if (has_overhang) {
                //enforce
                for (ExtrusionPath& path : paths) {
                    if (path.role() == erPerimeter || path.role() == erExternalPerimeter) {
                        path.set_role(erOverhangPerimeter);
                    }
                }
            }
        }

        // Append paths to collection.
        if (!paths.empty()) {
            for (size_t idx_path = 0; idx_path < paths.size(); ++idx_path) {
                if (idx_path > 0)
                    assert(paths[idx_path - 1].last_point().coincides_with_epsilon(paths[idx_path].first_point()));
                for (size_t idx_pt = 1; idx_pt < paths[idx_path].size(); ++idx_pt) {
                    assert(!paths[idx_path].polyline.get_points()[idx_pt - 1].coincides_with_epsilon(paths[idx_path].polyline.get_points()[idx_pt]));
                }
            }
            if (extrusion->is_closed) {
                assert(paths.back().last_point().coincides_with_epsilon(paths.front().first_point()));
                ExtrusionLoop extrusion_loop(std::move(paths), loop_role);
                // Restore the orientation of the extrusion loop.
                //TODO: use if (loop.is_steep_overhang && this->layer->id() % 2 == 1) to make_clockwise => need to detect is_steep_overhang on the arachne path
                if ((config->perimeter_reverse /*|| pg_extrusion.is_steep_overhang*/ && this->layer->id() % 2 == 1) == pg_extrusion.is_contour)
                    extrusion_loop.make_counter_clockwise();
                else
                    extrusion_loop.make_clockwise();
#if _DEBUG
                for (auto it = std::next(extrusion_loop.paths.begin()); it != extrusion_loop.paths.end(); ++it) {
                    assert(it->polyline.size() >= 2);
                    assert(std::prev(it)->last_point() == it->first_point());
                }
                // first & last points can be very near each other but sometimes not exactly.
                assert(extrusion_loop.paths.front().first_point().coincides_with_epsilon(extrusion_loop.paths.back().last_point()));
#endif
                //ensure the start & end points are the same.
                extrusion_loop.paths.front().polyline.set_points().front() = extrusion_loop.paths.back().last_point();
                assert(extrusion_loop.paths.front().first_point() == (extrusion_loop.paths.back().last_point()));

                extrusion_coll.append(std::move(extrusion_loop));
            } else {

                // Because we are processing one ExtrusionLine all ExtrusionPaths should form one connected path.
                // But there is possibility that due to numerical issue there is poss
                assert([&paths = std::as_const(paths)]() -> bool {
                    for (auto it = std::next(paths.begin()); it != paths.end(); ++it)
                        if (std::prev(it)->last_point() != it->first_point())
                            return false;
                    return true;
                }());
                ExtrusionMultiPath multi_path;
                multi_path.paths.push_back(std::move(paths.front()));
                multi_path.set_can_reverse(true);

                for (auto it_path = std::next(paths.begin()); it_path != paths.end(); ++it_path) {
                    if (!multi_path.paths.back().last_point().coincides_with_epsilon(it_path->first_point())) {
                        extrusion_coll.append(std::move(multi_path));
                        multi_path = ExtrusionMultiPath();
                        multi_path.set_can_reverse(true);
                    }
                    it_path->set_can_reverse(false);
                    multi_path.paths.push_back(std::move(*it_path));
                }

                extrusion_coll.append(ExtrusionMultiPath(std::move(multi_path)));
            }
        }
    }

    return extrusion_coll;
}

void PerimeterGenerator::_merge_thin_walls(ExtrusionEntityCollection &extrusions, ThickPolylines &thin_walls) const {
    //TODO: find a way to avoid double copy (from EntityCollection to ChangeFlow to searcher.search_result.loop
    class ChangeFlow : public ExtrusionVisitor {
    public:
        ChangeFlow(coordf_t resolution) : resolution_sqr(resolution * resolution) {}
        float percent_extrusion;
        std::vector<ExtrusionPath> paths;
        const Point* first_point = nullptr;
        coordf_t resolution_sqr;
        //TODO real travel with role & width
        void ensure_travel_to(const Point &pt) {
            assert(!paths.empty());
            Point last_point = paths.back().last_point();
            if (last_point != pt) {
                if (last_point.distance_to_square(pt) < resolution_sqr) {
                    paths.back().polyline.set_points().back() = pt;
                } else {
                    //add travel
                    ExtrusionPath travel(paths.back().role(), false);
                    travel.width = paths.back().width;
                    travel.height = paths.back().height;
                    travel.mm3_per_mm = 0;
                    travel.polyline.append(last_point);
                    travel.polyline.append(pt);
                    paths.push_back(travel);
                }
            }
        }
        virtual void use(ExtrusionPath &path) override {
            //ensure the loop is continue.
            if (first_point != nullptr) {
                if (*first_point != path.first_point()) {
                    if (first_point->distance_to_square(path.first_point()) < resolution_sqr) {
                        path.polyline.set_points()[0] = *first_point;
                    } else {
                        //add travel
                        ExtrusionPath travel(path.role(), false);
                        travel.width = path.width;
                        travel.height = path.height;
                        travel.mm3_per_mm = 0;
                        travel.polyline.append(*first_point);
                        travel.polyline.append(path.first_point());
                        paths.push_back(travel);
                    }
                }
                first_point = nullptr;
            }
            path.mm3_per_mm *= percent_extrusion;
            path.width *= percent_extrusion;
            paths.push_back(path);
        }
        virtual void use(ExtrusionPath3D &path3D) override { /*shouldn't happen*/ }
        virtual void use(ExtrusionMultiPath &multipath) override { /*shouldn't happen*/ }
        virtual void use(ExtrusionMultiPath3D &multipath) { /*shouldn't happen*/ }
        virtual void use(ExtrusionLoop &loop) override {
            for (ExtrusionPath &path : loop.paths)
                this->use(path);
        }
        virtual void use(ExtrusionEntityCollection &collection) override {
            for (ExtrusionEntity *entity : collection.entities())
                entity->visit(*this);
        }
    };
    struct BestPoint {
        //Point p;
        ExtrusionPath *path;
        size_t idx_path;
        ExtrusionLoop *loop;
        size_t idx_line;
        Line line;
        double dist;
        bool from_start;
    };
    //use a visitor to found the best point.
    class SearchBestPoint : public ExtrusionVisitor {
    public:
        ThickPolyline* thin_wall;
        BestPoint search_result;
        size_t idx_path;
        ExtrusionLoop *current_loop = nullptr;
        virtual void use(ExtrusionPath &path) override {
            //don't consider other thin walls.
            if (path.role() == erThinWall) return;
            //for each segment
            Lines lines = path.polyline.lines();
            for (size_t idx_line = 0; idx_line < lines.size(); idx_line++) {
                //look for nearest point
                double dist = lines[idx_line].distance_to_squared(thin_wall->front());
                if (dist < search_result.dist) {
                    search_result.path = &path;
                    search_result.idx_path = idx_path;
                    search_result.idx_line = idx_line;
                    search_result.line = lines[idx_line];
                    search_result.dist = dist;
                    search_result.from_start = true;
                    search_result.loop = current_loop;
                }
                dist = lines[idx_line].distance_to_squared(thin_wall->back());
                if (dist < search_result.dist) {
                    search_result.path = &path;
                    search_result.idx_path = idx_path;
                    search_result.idx_line = idx_line;
                    search_result.line = lines[idx_line];
                    search_result.dist = dist;
                    search_result.from_start = false;
                    search_result.loop = current_loop;
                }
            }
        }
        virtual void use(ExtrusionPath3D &path3D) override { /*shouldn't happen*/ }
        virtual void use(ExtrusionMultiPath &multipath) override { /*shouldn't happen*/ }
        virtual void use(ExtrusionMultiPath3D &multipath) { /*shouldn't happen*/ }
        virtual void use(ExtrusionLoop &loop) override {
            ExtrusionLoop * last_loop = current_loop;
            current_loop = &loop;
            //for each extrusion path
            idx_path = 0;
            for (ExtrusionPath &path : loop.paths) {
                this->use(path);
                idx_path++;
            }
            current_loop = last_loop;
        }
        virtual void use(ExtrusionEntityCollection &collection) override {
            //for each loop? (or other collections)
            for (ExtrusionEntity *entity : collection.entities())
                entity->visit(*this);
        }
    };
    //max dist to branch: ~half external perimeter width
    coord_t max_width = this->ext_perimeter_flow.scaled_width();
    SearchBestPoint searcher;
    ThickPolylines not_added;
    //search the best extusion/point to branch into
     //for each thin wall
    int idx = 0;
    for (ThickPolyline &tw : thin_walls) {
        searcher.thin_wall = &tw;
        searcher.search_result.dist = double(max_width);
        searcher.search_result.dist *= searcher.search_result.dist;
        searcher.search_result.path = nullptr;
        searcher.use(extrusions);
        idx++;
        //now insert thin wall if it has a point
        //it found a segment
        if (searcher.search_result.path != nullptr) {
#if _DEBUG
            searcher.search_result.loop->visit(LoopAssertVisitor{});
            ExtrusionLoop orig_loop = *searcher.search_result.loop;
#endif
            if (!searcher.search_result.from_start)
                tw.reverse();
            //save old path, as it may be destroyed before being re-created and we want to keep its parameters.
            ExtrusionPath path_to_split = *searcher.search_result.path; // TODO: 2.7: just save hte pathsettigns
            //get the point
            Point point = tw.front().projection_onto(searcher.search_result.line);
            //we have to create 3 paths: 1: thinwall extusion, 2: thinwall return, 3: end of the path
            //create new path : end of the path
            Polyline poly_after;
            poly_after.points.push_back(point);
            poly_after.points.insert(poly_after.points.end(),
                searcher.search_result.path->polyline.get_points().begin() + searcher.search_result.idx_line + 1,
                searcher.search_result.path->polyline.get_points().end());
            searcher.search_result.path->polyline.set_points().erase(
                searcher.search_result.path->polyline.set_points().begin() + searcher.search_result.idx_line + 1,
                searcher.search_result.path->polyline.set_points().end());
            size_t idx_path_before = searcher.search_result.idx_path;
            size_t idx_path_to_add = idx_path_before + 1;
            //check if the first part of the split polyline is long enough.
            assert(!searcher.search_result.loop->paths[idx_path_before].empty());
            if (searcher.search_result.loop->paths[idx_path_before].length() 
                + searcher.search_result.loop->paths[idx_path_before].polyline.back().distance_to(point) < SCALED_EPSILON) {
                //not long enough, move point to first point and destroy it
                assert(!searcher.search_result.loop->paths[idx_path_before].empty());
                point = searcher.search_result.loop->paths[idx_path_before].first_point();
                searcher.search_result.loop->paths.erase(searcher.search_result.loop->paths.begin() + idx_path_before);
                idx_path_to_add = idx_path_before;
                idx_path_before--;
            } else {
                //long enough
                if (!searcher.search_result.loop->paths[idx_path_before].polyline.back().coincides_with_epsilon(point)) {
                    // add the point to previous path
                    searcher.search_result.loop->paths[idx_path_before].polyline.append(point);
                } else {
                    //point too near of the start, just move the point
                    point                     = searcher.search_result.loop->paths[idx_path_before].polyline.back();
                    poly_after.points.front() = point;
                }
            }
            // remove next point if too near to point for the poly_after
            if (poly_after.size() > 1 && poly_after.points[0].coincides_with_epsilon(poly_after.points[1])) {
                poly_after.points.erase(poly_after.begin() + 1);
            }
            assert(idx_path_before > searcher.search_result.loop->paths.size() || searcher.search_result.loop->paths[idx_path_before].size() >= 2);
            assert(idx_path_before > searcher.search_result.loop->paths.size() || searcher.search_result.loop->paths[idx_path_before].length() > SCALED_EPSILON);
            // check if poly_after is big enough to be added
            if (poly_after.points.size() <= 1 || poly_after.length()  < SCALED_EPSILON) {
                //if there is a path after it, then move a little bit the first point
                if (searcher.search_result.loop->paths.size() > idx_path_to_add) {
                    assert(searcher.search_result.loop->paths[idx_path_to_add].first_point().coincides_with_epsilon(point));
                    searcher.search_result.loop->paths[idx_path_to_add].polyline.set_points().front() = point;
                    // set return point, to move our return point
                    poly_after.clear();
                    poly_after.points.push_back(point);
                } else {
                    //use last point as the end pos
                    Point end_pt = poly_after.back();
                    poly_after.clear();
                    poly_after.points.push_back(end_pt);
                }
            } else {
                assert(poly_after.length() > SCALED_EPSILON);
                searcher.search_result.loop->paths.insert(searcher.search_result.loop->paths.begin() + idx_path_to_add, 
                    ExtrusionPath(poly_after, path_to_split));
            }
            assert(idx_path_before > searcher.search_result.loop->paths.size() || searcher.search_result.loop->paths[idx_path_before].polyline.size() > 1);
            assert(poly_after.size() > 0);
            
            //create thin wall path exttrusion
            ExtrusionEntityCollection tws;
            tws.append(Geometry::thin_variable_width({ tw }, erThinWall, this->ext_perimeter_flow, std::max(ext_perimeter_flow.scaled_width() / 10, scale_t(this->print_config->resolution)), false));
            assert(!tws.entities().empty());
#if _DEBUG
                tws.visit(LoopAssertVisitor{});
#endif
            ChangeFlow change_flow(std::max(scale_t(this->print_config->resolution), SCALED_EPSILON));
            if (tws.entities().size() == 1 && tws.entities()[0]->is_loop()) {
                //loop, just add it 
                change_flow.first_point = &point;
                change_flow.percent_extrusion = 1;
                change_flow.use(tws);
                // ChangeFlow added the first move if needed, now add the second
                change_flow.ensure_travel_to(point);
                //add move around
                searcher.search_result.loop->paths.insert(searcher.search_result.loop->paths.begin() + idx_path_to_add,
                    change_flow.paths.begin(), change_flow.paths.end());
                ////add move to -> ??? i don't remember why i wrote that, so here it's removed.
                assert(poly_after.first_point() == point);
                //if (poly_after.first_point() != point) {
                //    assert(poly_after.first_point().coincides_with_epsilon(point));
                //    assert(searcher.search_result.loop->paths.size() > idx_path_to_add);
                //    assert(poly_after.first_point().coincides_with_epsilon(searcher.search_result.loop->paths[idx_path_to_add].polyline.set_points().front()));
                //    searcher.search_result.loop->paths[idx_path_to_add].polyline.set_points().front() = poly_after.first_point();
                //}
#if _DEBUG
                searcher.search_result.loop->visit(LoopAssertVisitor{});
#endif
            } else {
                //first add the return path
                //ExtrusionEntityCollection tws_second = tws; // this does a deep copy
                change_flow.first_point = &poly_after.first_point(); // end at the start of the next path
                change_flow.percent_extrusion = 0.1f;
                change_flow.use(tws); // tws_second); //does not need the deep copy if the change_flow copy the content instead of re-using it.
                //force reverse
                for (ExtrusionPath &path : change_flow.paths)
                    path.reverse();
                std::reverse(change_flow.paths.begin(), change_flow.paths.end());
                //std::reverse(change_flow.paths.begin(), change_flow.paths.end());
                searcher.search_result.loop->paths.insert(searcher.search_result.loop->paths.begin() + idx_path_to_add,
                    change_flow.paths.begin(), change_flow.paths.end()); //TODO 2.7:change role by a kind of thinwalltravel that won't be considered for seam
                //add the real extrusion path
                change_flow.first_point = &point; // start at the end of previous extrusion
                change_flow.percent_extrusion = 9.f; // 0.9 but as we modified it by 0.1 just before, has to multiply by 10
                change_flow.paths = std::vector<ExtrusionPath>();
                change_flow.use(tws);
#if _DEBUG
                for (ExtrusionPath &path : change_flow.paths)
                    path.visit(LoopAssertVisitor{});
#endif
                searcher.search_result.loop->paths.insert(searcher.search_result.loop->paths.begin() + idx_path_to_add,
                    change_flow.paths.begin(), change_flow.paths.end());
#if _DEBUG
                searcher.search_result.loop->visit(LoopAssertVisitor{});
#endif
            }
        } else {
            not_added.push_back(tw);
        }
    }

    //now add thinwalls that have no anchor (make them reversable)
    extrusions.append(Geometry::thin_variable_width(not_added, erThinWall, this->ext_perimeter_flow, std::max(ext_perimeter_flow.scaled_width() / 4, scale_t(this->print_config->resolution)), true));
#if _DEBUG
        extrusions.visit(LoopAssertVisitor{});
#endif
}

PerimeterIntersectionPoint
PerimeterGenerator::_get_nearest_point(const PerimeterGeneratorLoops &children, ExtrusionLoop &myPolylines, const coord_t dist_cut, const coord_t max_dist) const {
    //find best points of intersections
    PerimeterIntersectionPoint intersect;
    intersect.distance = 0x7FFFFFFF; // ! assumption on intersect type & max value
    intersect.idx_polyline_outter = -1;
    intersect.idx_children = -1;
    for (size_t idx_child = 0; idx_child < children.size(); idx_child++) {
        const PerimeterGeneratorLoop &child = children[idx_child];
        for (size_t idx_poly = 0; idx_poly < myPolylines.paths.size(); idx_poly++) {
            //if (myPolylines.paths[idx_poly].extruder_id == (unsigned int)-1) continue;
            if (myPolylines.paths[idx_poly].length() < dist_cut + perimeter_flow.scaled_width()/20) continue;

            if ((myPolylines.paths[idx_poly].role() == erExternalPerimeter || child.is_external() )
                && (this->object_config->seam_position.value != SeamPosition::spRandom && this->object_config->seam_position.value != SeamPosition::spAllRandom)) {
                //first, try to find 2 point near enough //TODO: use seam placer or at least an equivalent.
                for (size_t idx_point = 0; idx_point < myPolylines.paths[idx_poly].polyline.size(); idx_point++) {
                    const Point &p = myPolylines.paths[idx_poly].polyline.get_points()[idx_point];
                    const Point &nearest_p = *child.polygon.closest_point(p);
                    const double dist = nearest_p.distance_to(p);
                    //Try to find a point in the far side, aligning them
                    if (dist + dist_cut / 20 < intersect.distance || 
                        (config->perimeter_loop_seam.value == spRear && (intersect.idx_polyline_outter <0 || p.y() > intersect.outter_best.y())
                            && dist <= max_dist && intersect.distance + dist_cut / 20)) {
                        //ok, copy the idx
                        intersect.distance = (coord_t)nearest_p.distance_to(p);
                        intersect.idx_children = idx_child;
                        intersect.idx_polyline_outter = idx_poly;
                        intersect.outter_best = p;
                        intersect.child_best = nearest_p;
                    }
                }
            } else {
                //first, try to find 2 point near enough
                for (size_t idx_point = 0; idx_point < myPolylines.paths[idx_poly].polyline.size(); idx_point++) {
                    const Point &p = myPolylines.paths[idx_poly].polyline.get_points()[idx_point];
                    const Point &nearest_p = *child.polygon.closest_point(p);
                    const double dist = nearest_p.distance_to(p);
                    if (dist + SCALED_EPSILON < intersect.distance || 
                        (config->perimeter_loop_seam.value == spRear && (intersect.idx_polyline_outter<0 || p.y() < intersect.outter_best.y())
                            && dist <= max_dist && intersect.distance + dist_cut / 20)) {
                        //ok, copy the idx
                        intersect.distance = (coord_t)nearest_p.distance_to(p);
                        intersect.idx_children = idx_child;
                        intersect.idx_polyline_outter = idx_poly;
                        intersect.outter_best = p;
                        intersect.child_best = nearest_p;
                    }
                }
            }
        }
    }
    if (intersect.distance <= max_dist) {
        return intersect;
    }

    for (size_t idx_child = 0; idx_child < children.size(); idx_child++) {
        const PerimeterGeneratorLoop &child = children[idx_child];
        for (size_t idx_poly = 0; idx_poly < myPolylines.paths.size(); idx_poly++) {
            //if (myPolylines.paths[idx_poly].extruder_id == (unsigned int)-1) continue;
            if (myPolylines.paths[idx_poly].length() < dist_cut + perimeter_flow.scaled_width() / 20) continue;

            //second, try to check from one of my points
            //don't check the last point, as it's used to go outter, can't use it to go inner.
            for (size_t idx_point = 1; idx_point < myPolylines.paths[idx_poly].polyline.size()-1; idx_point++) {
                const Point &p = myPolylines.paths[idx_poly].polyline.get_points()[idx_point];
                Point nearest_p = child.polygon.point_projection(p).first;
                coord_t dist = (coord_t)nearest_p.distance_to(p);
                //if no projection, go to next
                if (dist == 0) continue;
                if (dist + SCALED_EPSILON / 2 < intersect.distance) {
                    //ok, copy the idx
                    intersect.distance = dist;
                    intersect.idx_children = idx_child;
                    intersect.idx_polyline_outter = idx_poly;
                    intersect.outter_best = p;
                    intersect.child_best = nearest_p;
                }
            }
        }
    }
    if (intersect.distance <= max_dist) {
        return intersect;
    }

    for (size_t idx_child = 0; idx_child < children.size(); idx_child++) {
        const PerimeterGeneratorLoop &child = children[idx_child];
        for (size_t idx_poly = 0; idx_poly < myPolylines.paths.size(); idx_poly++) {
            //if (myPolylines.paths[idx_poly].extruder_id == (unsigned int)-1) continue;
            if (myPolylines.paths[idx_poly].length() < dist_cut + perimeter_flow.scaled_width() / 20) continue;

            //lastly, try to check from one of his points
            for (size_t idx_point = 0; idx_point < child.polygon.size(); idx_point++) {
                const Point &p = child.polygon.points[idx_point];
                Point nearest_p = myPolylines.paths[idx_poly].polyline.point_projection(p).first;
                coord_t dist = (coord_t)nearest_p.distance_to(p);
                //if no projection, go to next
                if (dist == 0) continue;
                if (dist + SCALED_EPSILON / 2 < intersect.distance) {
                    //ok, copy the idx
                    intersect.distance = dist;
                    intersect.idx_children = idx_child;
                    intersect.idx_polyline_outter = idx_poly;
                    intersect.outter_best = nearest_p;
                    intersect.child_best = p;
                }
            }
        }
    }
    return intersect;
}


ExtrusionLoop
PerimeterGenerator::_extrude_and_cut_loop(const PerimeterGeneratorLoop &loop, const Point entry_point, const Line &direction, bool enforce_loop) const
{

    bool need_to_reverse = false;
    Polyline initial_polyline;
    coord_t dist_cut = (coord_t)scale_(this->print_config->nozzle_diameter.get_at(this->config->perimeter_extruder - 1));

    //fuzzify first in this case, as it's a bit complicated to do it after.
    Polygon fuzzy_poly;
    if (loop.fuzzify) {
        fuzzy_poly = loop.polygon;
        double nozle_diameter = loop.is_external() ? this->ext_perimeter_flow.nozzle_diameter() : this->perimeter_flow.nozzle_diameter();
        double fuzzy_skin_thickness = config->fuzzy_skin_thickness.get_abs_value(nozle_diameter);
        double fuzzy_skin_point_dist = config->fuzzy_skin_point_dist.get_abs_value(nozle_diameter);
        fuzzy_polygon(fuzzy_poly, scale_d(fuzzy_skin_thickness), scale_d(fuzzy_skin_point_dist));
    }
    const Polygon& poly_to_use = loop.fuzzify ? fuzzy_poly : loop.polygon;

    if (poly_to_use.size() < 3) return ExtrusionLoop(elrDefault);
    if (poly_to_use.length() < dist_cut * 2) {
        if (enforce_loop) {
            //do something to still use it
            dist_cut = poly_to_use.length() / 4;
        } else {
            //reduce it ot a single-point loop that will eb emrged inside the complex path
            ExtrusionLoop single_point(elrDefault);
            Polyline poly_point;
            poly_point.append(poly_to_use.centroid());
            single_point.paths.emplace_back(
                loop.is_external() ? erExternalPerimeter : erPerimeter,
                (double)(loop.is_external() ? this->ext_mm3_per_mm() : this->mm3_per_mm()),
                (float)(loop.is_external() ? this->ext_perimeter_flow.width() : this->perimeter_flow.width()),
                (float)(this->layer->height),
                false/*can't reverse*/);
            single_point.paths.back().polyline = poly_point;
            return single_point;
        }
    }
    const size_t idx_closest_from_entry_point = poly_to_use.closest_point_index(entry_point);
    if (poly_to_use.points[idx_closest_from_entry_point].distance_to(entry_point) > SCALED_EPSILON * 2) {
        //create new Point
        //get first point
        size_t idx_before = -1;
        for (size_t idx_p_a = 0; idx_p_a < poly_to_use.points.size(); ++idx_p_a) {
            Line l(poly_to_use.points[idx_p_a], poly_to_use.points[(idx_p_a + 1 == poly_to_use.points.size()) ? 0 : (idx_p_a + 1)]);
            if (entry_point.distance_to(l) < SCALED_EPSILON) {
                idx_before = idx_p_a;
                break;
            }
        }
        if (idx_before == (size_t)-1) std::cerr << "ERROR: _traverse_and_join_loops : idx_before can't be finded to create new point\n";
        initial_polyline = poly_to_use.split_at_index(idx_before);
        initial_polyline.points.push_back(entry_point);
        initial_polyline.points[0] = entry_point;
    } else {
        initial_polyline = poly_to_use.split_at_index(idx_closest_from_entry_point);
    }
   

    //std::vector<PerimeterPolylineNode> myPolylines;
    ExtrusionLoop my_loop;

    //overhang / notoverhang
    {
        bool is_external = loop.is_external();

        ExtrusionRole role;
        ExtrusionLoopRole loop_role;
        role = is_external ? erExternalPerimeter : erPerimeter;
        if (loop.is_internal_contour()) {
            // Note that we set loop role to ContourInternalPerimeter
            // also when loop is both internal and external (i.e.
            // there's only one contour loop).
            loop_role = elrInternal;
        } else {
            loop_role = elrDefault;
        }
        if (!loop.is_contour) {
            loop_role = (ExtrusionLoopRole)(loop_role | elrHole);
        }

        // detect overhanging/bridging perimeters
        if ( this->config->overhangs_width_speed.value > 0 && this->layer->id() > 0
            && !(this->object_config->support_material && this->object_config->support_material_contact_distance_type.value == zdNone)) {
            ExtrusionPaths paths = this->create_overhangs(initial_polyline, role, is_external);
            
            if (direction.length() > 0) {
                Polyline direction_polyline;
                for (ExtrusionPath &path : paths) {
                    if(direction_polyline.size() == 0 || direction_polyline.points.back() != path.first_point())
                        direction_polyline.points.insert(direction_polyline.points.end(), path.polyline.get_points().begin(), path.polyline.get_points().end());
                }
                for (int i = 0; i < direction_polyline.points.size() - 1; i++)
                    assert(direction_polyline.points[i] != direction_polyline.points[i + 1]);
                if (direction_polyline.length() > perimeter_flow.scaled_width() / 8) {
                    direction_polyline.clip_start(perimeter_flow.scaled_width() / 20);
                    direction_polyline.clip_end(perimeter_flow.scaled_width() / 20);
                }
                coord_t dot = direction.dot(Line(direction_polyline.points.back(), direction_polyline.points.front()));
                need_to_reverse = dot>0;
            }
            if (need_to_reverse) {
                std::reverse(paths.begin(), paths.end());
            }
            //search for the first path
            size_t good_idx = 0; 
            for (size_t idx_path = 0; idx_path < paths.size(); idx_path++) {
                const ExtrusionPath &path = paths[idx_path];
                if (need_to_reverse) {
                    if (path.polyline.back().coincides_with_epsilon(initial_polyline.front())) {
                        good_idx = idx_path;
                        break;
                    }
                } else {
                    if (path.polyline.front().coincides_with_epsilon(initial_polyline.front())) {
                        good_idx = idx_path;
                        break;
                    }
                }
            }
            for (size_t idx_path = good_idx; idx_path < paths.size(); idx_path++) {
                ExtrusionPath &path = paths[idx_path];
                if (need_to_reverse) path.reverse();
                my_loop.paths.push_back(path);
            }
            for (size_t idx_path = 0; idx_path < good_idx; idx_path++) {
                ExtrusionPath &path = paths[idx_path];
                if (need_to_reverse) path.reverse();
                my_loop.paths.push_back(path);
            }
        } else {

            if (direction.length() > 0) {
                Polyline direction_polyline = initial_polyline;
                direction_polyline.clip_start(perimeter_flow.scaled_width() / 20);
                direction_polyline.clip_end(perimeter_flow.scaled_width() / 20);
                coord_t dot = direction.dot(Line(direction_polyline.back(), direction_polyline.front()));
                need_to_reverse = dot>0;
            }

            ExtrusionPath path(role, false);
            path.polyline = initial_polyline;
            if (need_to_reverse) path.polyline.reverse();
            path.mm3_per_mm = is_external ? this->ext_mm3_per_mm() : this->mm3_per_mm();
            path.width = is_external ? this->ext_perimeter_flow.width() : this->perimeter_flow.width();
            path.height = (float)(this->layer->height);
            my_loop.paths.push_back(path);
        }

    }
    
    return my_loop;
}

ExtrusionLoop
PerimeterGenerator::_traverse_and_join_loops(const PerimeterGeneratorLoop &loop, const PerimeterGeneratorLoops &children, const Point entry_point) const
{
    //std::cout << " === ==== _traverse_and_join_loops ==== ===\n";
    // other perimeters
    //this->_mm3_per_mm = this->perimeter_flow.mm3_per_mm();
    //coord_t perimeter_width = this->perimeter_flow.scaled_width();
    const coord_t perimeter_spacing = this->perimeter_flow.scaled_spacing();

    //// external perimeters
    //this->_ext_mm3_per_mm = this->ext_perimeter_flow.mm3_per_mm();
    //coord_t ext_perimeter_width = this->ext_perimeter_flow.scaled_width();
    const coord_t ext_perimeter_spacing = this->ext_perimeter_flow.scaled_spacing();
    //coord_t ext_perimeter_spacing2 = this->ext_perimeter_flow.scaled_spacing(this->perimeter_flow);

    //const coord_t dist_cut = (coord_t)scale_(this->print_config->nozzle_diameter.get_at(this->config->perimeter_extruder - 1));
    //TODO change this->external_perimeter_flow.scaled_width() if it's the first one!
    const coord_t max_width_extrusion = this->perimeter_flow.scaled_width();
    ExtrusionLoop my_loop = _extrude_and_cut_loop(loop, entry_point, Line{ {0,0},{0,0} }, true);

    int child_idx = 0;
    //Polylines myPolylines = { myPolyline };
    //iterate on each point ot find the best place to go into the child
    PerimeterGeneratorLoops childs = children;
    while (!childs.empty()) {
        child_idx++;
        PerimeterIntersectionPoint nearest = this->_get_nearest_point(childs, my_loop, coord_t(this->perimeter_flow.scaled_width()), coord_t(this->perimeter_flow.scaled_width()* 1.42));
        if (nearest.idx_children == (size_t)-1) {
            //return ExtrusionEntityCollection();
            break;
        } else {
            const PerimeterGeneratorLoop &child = childs[nearest.idx_children];
            //std::cout << "c." << child_idx << " === i have " << my_loop.paths.size() << " paths" << " == cut_path_is_ccw size " << path_is_ccw.size() << "\n";
            //std::cout << "change to child " << nearest.idx_children << " @ " << unscale(nearest.outter_best.x) << ":" << unscale(nearest.outter_best.y)
            //    << ", idxpolyline = " << nearest.idx_polyline_outter << "\n";
            //PerimeterGeneratorLoops less_childs = childs;
            //less_childs.erase(less_childs.begin() + nearest.idx_children);
            //create new node with recursive ask for the inner perimeter & COPY of the points, ready to be cut
            my_loop.paths.insert(my_loop.paths.begin() + nearest.idx_polyline_outter + 1, my_loop.paths[nearest.idx_polyline_outter]);

            // outer_start == outer_end
            ExtrusionPath *outer_start = &my_loop.paths[nearest.idx_polyline_outter];
            ExtrusionPath *outer_end = &my_loop.paths[nearest.idx_polyline_outter + 1];
            Line deletedSection;

            //cut our polyline, so outer_start has no common point with outer_end
            //separate them
            size_t nearest_idx_outter = outer_start->polyline.closest_point_index(nearest.outter_best);
            if (outer_start->polyline.get_points()[nearest_idx_outter].coincides_with_epsilon(nearest.outter_best)) {
                if (nearest_idx_outter < outer_start->polyline.size() - 1) {
                    outer_start->polyline.set_points().erase(outer_start->polyline.set_points().begin() + nearest_idx_outter + 1, outer_start->polyline.set_points().end());
                }
                if (nearest_idx_outter > 0) {
                    outer_end->polyline.set_points().erase(outer_end->polyline.set_points().begin(), outer_end->polyline.set_points().begin() + nearest_idx_outter);
                }
            } else {
                //get first point
                size_t idx_before = -1;
                for (size_t idx_p_a = 0; idx_p_a < outer_start->polyline.size() - 1; ++idx_p_a) {
                    Line l(outer_start->polyline.get_points()[idx_p_a], outer_start->polyline.get_points()[idx_p_a + 1]);
                    if (nearest.outter_best.distance_to(l) < SCALED_EPSILON) {
                        idx_before = idx_p_a;
                        break;
                    }
                }
                if (idx_before == (size_t)-1) {
                    std::cerr << "ERROR: idx_before can't be finded\n";
                    continue;
                }

                outer_start->polyline.set_points().erase(outer_start->polyline.set_points().begin() + idx_before + 1, outer_start->polyline.set_points().end());
                outer_start->polyline.append(nearest.outter_best);

                if (idx_before < outer_end->polyline.size()-1)
                    outer_end->polyline.set_points().erase(outer_end->polyline.set_points().begin(), outer_end->polyline.set_points().begin() + idx_before + 1);
                else
                    outer_end->polyline.set_points().erase(outer_end->polyline.set_points().begin()+1, outer_end->polyline.set_points().end());
                outer_end->polyline.append_before(nearest.outter_best);
            }
            Polyline to_reduce = outer_start->polyline.as_polyline();
            if (to_reduce.size()>1 && to_reduce.length() > (perimeter_flow.scaled_width() / 10)) to_reduce.clip_end(perimeter_flow.scaled_width() / 20);
            deletedSection.a = to_reduce.back();
            to_reduce = outer_end->polyline.as_polyline();
            if (to_reduce.size()>1 && to_reduce.length() > (perimeter_flow.scaled_width() / 10)) to_reduce.clip_start(perimeter_flow.scaled_width() / 20);
            deletedSection.b = to_reduce.front();
            
            //get the inner loop to connect to us.
            ExtrusionLoop child_loop = _extrude_and_cut_loop(child, nearest.child_best, deletedSection);

            const coord_t inner_child_spacing = child.is_external() ? ext_perimeter_spacing : perimeter_spacing;
            const coord_t outer_start_spacing = scale_(outer_start->width - outer_start->height * (1. - 0.25 * PI));
            const coord_t outer_end_spacing = scale_(outer_end->width - outer_end->height * (1. - 0.25 * PI));

            //FIXME: if child_loop has no point or 1 point or not enough space !!!!!!!
            const size_t child_paths_size = child_loop.paths.size();
            if (child_paths_size == 0) continue;
            my_loop.paths.insert(my_loop.paths.begin() + nearest.idx_polyline_outter + 1, child_loop.paths.begin(), child_loop.paths.end());
            
            //add paths into my_loop => need to re-get the refs
            outer_start = &my_loop.paths[nearest.idx_polyline_outter];
            outer_end = &my_loop.paths[nearest.idx_polyline_outter + child_paths_size + 1];
            ExtrusionPath *inner_start = &my_loop.paths[nearest.idx_polyline_outter+1];
            ExtrusionPath *inner_end = &my_loop.paths[nearest.idx_polyline_outter + child_paths_size];
            //TRIM
            //choose trim direction
            if (outer_start->polyline.size() == 1 && outer_end->polyline.size() == 1) {
                //do nothing
            } else if (outer_start->polyline.size() == 1) {
                outer_end->polyline.clip_start(double(outer_end_spacing));
                if (inner_end->polyline.length() > inner_child_spacing)
                    inner_end->polyline.clip_end(double(inner_child_spacing));
                else
                    inner_end->polyline.clip_end(inner_end->polyline.length() / 2);
            } else if (outer_end->polyline.size() == 1) {
                outer_start->polyline.clip_end(double(outer_start_spacing));
                if (inner_start->polyline.length() > inner_child_spacing)
                    inner_start->polyline.clip_start(double(inner_child_spacing));
                else
                    inner_start->polyline.clip_start(inner_start->polyline.length()/2);
            } else {
                coord_t length_poly_1 = (coord_t)outer_start->polyline.length();
                coord_t length_poly_2 = (coord_t)outer_end->polyline.length();
                coord_t length_trim_1 = outer_start_spacing / 2;
                coord_t length_trim_2 = outer_end_spacing / 2;
                if (length_poly_1 < length_trim_1) {
                    length_trim_2 = length_trim_1 + length_trim_2 - length_poly_1;
                }
                if (length_poly_2 < length_trim_1) {
                    length_trim_1 = length_trim_1 + length_trim_2 - length_poly_2;
                }
                if (length_poly_1 > length_trim_1) {
                    outer_start->polyline.clip_end(double(length_trim_1));
                } else {
                    outer_start->polyline.set_points().erase(outer_start->polyline.set_points().begin() + 1, outer_start->polyline.set_points().end());
                }
                if (length_poly_2 > length_trim_2) {
                    outer_end->polyline.clip_start(double(length_trim_2));
                } else {
                    outer_end->polyline.set_points().erase(outer_end->polyline.set_points().begin(), outer_end->polyline.set_points().end() - 1);
                }
                
                length_poly_1 = coord_t(inner_start->polyline.length());
                length_poly_2 = coord_t(inner_end->polyline.length());
                length_trim_1 = inner_child_spacing / 2;
                length_trim_2 = inner_child_spacing / 2;
                if (length_poly_1 < length_trim_1) {
                    length_trim_2 = length_trim_1 + length_trim_2 - length_poly_1;
                }
                if (length_poly_2 < length_trim_1) {
                    length_trim_1 = length_trim_1 + length_trim_2 - length_poly_2;
                }
                if (length_poly_1 > length_trim_1) {
                    inner_start->polyline.clip_start(double(length_trim_1));
                } else {
                    inner_start->polyline.set_points().erase(
                        inner_start->polyline.set_points().begin(),
                        inner_start->polyline.set_points().end() - 1);
                }
                if (length_poly_2 > length_trim_2) {
                    inner_end->polyline.clip_end(double(length_trim_2));
                } else {
                    inner_end->polyline.set_points().erase(
                        inner_end->polyline.set_points().begin() + 1,
                        inner_end->polyline.set_points().end());
                }
            }

            //last check to see if we need a reverse
            {
                Line l1(outer_start->polyline.back(), inner_start->polyline.front());
                Line l2(inner_end->polyline.back(), outer_end->polyline.front());
                Point p_inter(0, 0);
                bool is_interect = l1.intersection(l2, &p_inter);
                if (is_interect && p_inter.distance_to(l1) < SCALED_EPSILON && p_inter.distance_to(l2) < SCALED_EPSILON) {
                    //intersection! need to reverse!
                    std::reverse(my_loop.paths.begin() + nearest.idx_polyline_outter + 1, my_loop.paths.begin() + nearest.idx_polyline_outter + child_paths_size + 1);
                    for (size_t idx = nearest.idx_polyline_outter + 1; idx < nearest.idx_polyline_outter + child_paths_size + 1; idx++) {
                        my_loop.paths[idx].reverse();
                    }
                    outer_start = &my_loop.paths[nearest.idx_polyline_outter];
                    inner_start = &my_loop.paths[nearest.idx_polyline_outter + 1];
                    inner_end = &my_loop.paths[nearest.idx_polyline_outter + child_paths_size];
                    outer_end = &my_loop.paths[nearest.idx_polyline_outter + child_paths_size + 1];
                }

            }

            //now add extrusionPAths to connect the two loops
            ExtrusionPaths travel_path_begin;// (ExtrusionRole::erTravel, 0, outer_start->width, outer_start->height);
            //travel_path_begin.extruder_id = -1;
            ExtrusionPaths travel_path_end;// (ExtrusionRole::erTravel, 0, outer_end->width, outer_end->height);
            //travel_path_end.extruder_id = -1;
            double dist_travel = outer_start->polyline.back().distance_to(inner_start->polyline.front());
            if (dist_travel > max_width_extrusion*1.5 && this->config->fill_density.value > 0) {
                travel_path_begin.emplace_back(ExtrusionRole::erPerimeter, outer_start->mm3_per_mm, outer_start->width, outer_start->height, false);
                travel_path_begin.emplace_back(ExtrusionRole::erTravel, 0, outer_start->width/10, outer_start->height, false);
                travel_path_begin.emplace_back(ExtrusionRole::erPerimeter, outer_start->mm3_per_mm, outer_start->width, outer_start->height, false);
                //travel_path_begin[0].extruder_id = -1;
                //travel_path_begin[1].extruder_id = -1;
                //travel_path_begin[2].extruder_id = -1;
                Line line(outer_start->polyline.back(), inner_start->polyline.front());
                Point p_dist_cut_extrude = (line.b - line.a);
                p_dist_cut_extrude.x() = (coord_t)(p_dist_cut_extrude.x() * ((double)max_width_extrusion) / (line.length() * 2));
                p_dist_cut_extrude.y() = (coord_t)(p_dist_cut_extrude.y() * ((double)max_width_extrusion) / (line.length() * 2));
                //extrude a bit after the turn, to close the loop
                Point p_start_travel = line.a;
                p_start_travel += p_dist_cut_extrude;
                travel_path_begin[0].polyline.append(outer_start->polyline.back());
                travel_path_begin[0].polyline.append(p_start_travel);
                //extrude a bit before the final turn, to close the loop
                Point p_end_travel = line.b;
                p_end_travel -= p_dist_cut_extrude;
                travel_path_begin[2].polyline.append(p_end_travel);
                travel_path_begin[2].polyline.append(inner_start->polyline.front());
                //fake travel in the middle
                travel_path_begin[1].polyline.append(p_start_travel);
                travel_path_begin[1].polyline.append(p_end_travel);
            } else {
                // the path is small enough to extrude all along.
                double flow_mult = 1;
                if (dist_travel > max_width_extrusion && this->config->fill_density.value > 0) {
                    // the path is a bit too long, reduce the extrusion flow.
                    flow_mult = max_width_extrusion / dist_travel;
                }
                travel_path_begin.emplace_back(ExtrusionRole::erPerimeter, outer_start->mm3_per_mm * flow_mult, (float)(outer_start->width * flow_mult), outer_start->height, false);
                //travel_path_begin[0].extruder_id = -1;
                travel_path_begin[0].polyline.append(outer_start->polyline.back());
                travel_path_begin[0].polyline.append(inner_start->polyline.front());
            }
            dist_travel = inner_end->polyline.back().distance_to(outer_end->polyline.front());
            if (dist_travel > max_width_extrusion*1.5 && this->config->fill_density.value > 0) {
                travel_path_end.emplace_back(ExtrusionRole::erPerimeter, outer_end->mm3_per_mm, outer_end->width, outer_end->height, false);
                travel_path_end.emplace_back(ExtrusionRole::erTravel, 0, outer_end->width/10, outer_end->height, false);
                travel_path_end.emplace_back(ExtrusionRole::erPerimeter, outer_end->mm3_per_mm, outer_end->width, outer_end->height, false);
                //travel_path_end[0].extruder_id = -1;
                //travel_path_end[1].extruder_id = -1;
                //travel_path_end[2].extruder_id = -1;
                Line line(inner_end->polyline.back(), outer_end->polyline.front());
                Point p_dist_cut_extrude = (line.b - line.a);
                p_dist_cut_extrude.x() = (coord_t)(p_dist_cut_extrude.x() * ((double)max_width_extrusion) / (line.length() * 2));
                p_dist_cut_extrude.y() = (coord_t)(p_dist_cut_extrude.y() * ((double)max_width_extrusion) / (line.length() * 2));
                //extrude a bit after the turn, to close the loop
                Point p_start_travel_2 = line.a;
                p_start_travel_2 += p_dist_cut_extrude;
                travel_path_end[0].polyline.append(inner_end->polyline.back());
                travel_path_end[0].polyline.append(p_start_travel_2);
                //extrude a bit before the final turn, to close the loop
                Point p_end_travel_2 = line.b;
                p_end_travel_2 -= p_dist_cut_extrude;
                travel_path_end[2].polyline.append(p_end_travel_2);
                travel_path_end[2].polyline.append(outer_end->polyline.front());
                //fake travel in the middle
                travel_path_end[1].polyline.append(p_start_travel_2);
                travel_path_end[1].polyline.append(p_end_travel_2);
            } else {
                // the path is small enough to extrude all along.
                double flow_mult = 1;
                if (dist_travel > max_width_extrusion && this->config->fill_density.value > 0) {
                    // the path is a bit too long, reduce the extrusion flow.
                    flow_mult = max_width_extrusion / dist_travel;
                }
                travel_path_end.emplace_back(ExtrusionRole::erPerimeter, outer_end->mm3_per_mm * flow_mult, (float)(outer_end->width * flow_mult), outer_end->height, false);
                //travel_path_end[0].extruder_id = -1;
                travel_path_end[0].polyline.append(inner_end->polyline.back());
                travel_path_end[0].polyline.append(outer_end->polyline.front());
            }
            //check if we add path or reuse bits
            //FIXME
            /*if (outer_start->polyline.points.size() == 1) {
                outer_start->polyline = travel_path_begin.front().polyline;
                travel_path_begin.erase(travel_path_begin.begin());
                outer_start->extruder_id = -1;
            } else if (outer_end->polyline.points.size() == 1) {
                outer_end->polyline = travel_path_end.back().polyline;
                travel_path_end.erase(travel_path_end.end() - 1);
                outer_end->extruder_id = -1;
            }*/
            //add paths into my_loop => after that all ref are wrong!
            for (size_t i = travel_path_end.size() - 1; i < travel_path_end.size(); i--) {
                my_loop.paths.insert(my_loop.paths.begin() + nearest.idx_polyline_outter + child_paths_size + 1, travel_path_end[i]);
            }
            for (size_t i = travel_path_begin.size() - 1; i < travel_path_begin.size(); i--) {
                my_loop.paths.insert(my_loop.paths.begin() + nearest.idx_polyline_outter + 1, travel_path_begin[i]);
            }
        }
        //remove one-point extrusion
        //FIXME prevent this instead of patching here?
        for (size_t i = 0; i < my_loop.paths.size(); i++) {
            if (my_loop.paths[i].polyline.size() < 2) {
                if (my_loop.paths[i].polyline.size() == 1)
                    BOOST_LOG_TRIVIAL(warning) << "erase one-point extrusion : layer " << this->layer->id() << " " << my_loop.paths[i].polyline.front().x() << ":" << my_loop.paths[i].polyline.front().y() << "\n";
                my_loop.paths.erase(my_loop.paths.begin() + i);
                i--;
            }
        }

        //update for next loop
        childs.erase(childs.begin() + nearest.idx_children);
    }

    return my_loop;
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

coord_t PerimeterGenerator::get_resolution(size_t perimeter_id, bool is_overhang, const Surface* srf) const
{
    coord_t reso = scale_t(this->print_config->resolution.value);
    if (reso == 0) reso = SCALED_EPSILON;
    return reso;
    //deactivated because with full perimeter on tube, the innermost perimeter can be very rough, and not a circle anymore.
    //if on top or bottom, use external resolution.
    //if (is_overhang || perimeter_id == 0)
    //    return reso;
    //if(srf && srf->has_pos_top())
    //    return reso;
    //// for each perimeter, reduce the precision by a factor 3
    //int mult = (int)std::pow(2, perimeter_id);
    //coord_t reso_internal = scale_t(this->print_config->resolution_internal.value);
    //if(reso_internal < reso * mult)
    //    return reso_internal;
    //return reso * mult;
}

}
