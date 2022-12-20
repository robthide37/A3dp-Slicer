#include "RegionExpansion.hpp"

#include <libslic3r/ClipperZUtils.hpp>
#include <libslic3r/ClipperUtils.hpp>

namespace Slic3r {
namespace Algorithm {

// similar to expolygons_to_zpaths(), but each contour is expanded before converted to zpath.
// The expanded contours are then opened (the first point is repeated at the end).
static ClipperLib_Z::Paths expolygons_to_zpaths_expanded_opened(
    const ExPolygons &src, const float expansion, coord_t base_idx)
{
    ClipperLib_Z::Paths out;
    out.reserve(2 * std::accumulate(src.begin(), src.end(), size_t(0),
        [](const size_t acc, const ExPolygon &expoly) { return acc + expoly.num_contours(); }));
    coord_t z = base_idx;
    ClipperLib::ClipperOffset offsetter;
    offsetter.ShortestEdgeLength = expansion * ClipperOffsetShortestEdgeFactor;
    ClipperLib::Paths expansion_cache;
    for (const ExPolygon &expoly : src) {
        for (size_t icontour = 0; icontour < expoly.num_contours(); ++ icontour) {
            // Execute reorients the contours so that the outer most contour has a positive area. Thus the output
            // contours will be CCW oriented even though the input paths are CW oriented.
            // Offset is applied after contour reorientation, thus the signum of the offset value is reversed.
            offsetter.Clear();
            offsetter.AddPath(expoly.contour_or_hole(icontour).points, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
            expansion_cache.clear();
            offsetter.Execute(expansion_cache, icontour == 0 ? expansion : -expansion);
            append(out, ClipperZUtils::to_zpaths<true>(expansion_cache, z));
        }
        ++ z;
    }
    return out;
}

// Paths were created by splitting closed polygons into open paths and then by clipping them.
// Thus some pieces of the clipped polygons may now become split at the ends of the source polygons.
// Those ends are sorted lexicographically in "splits".
// Reconnect those split pieces.
static inline void merge_splits(ClipperLib_Z::Paths &paths, std::vector<std::pair<ClipperLib_Z::IntPoint, int>> &splits)
{
    for (auto it_path = paths.begin(); it_path != paths.end(); ) {
        ClipperLib_Z::Path &path = *it_path;
        assert(path.size() >= 2);
        bool merged = false;
        if (path.size() >= 2) {
            const ClipperLib_Z::IntPoint &front = path.front();
            const ClipperLib_Z::IntPoint &back  = path.back();
            // The path before clipping was supposed to cross the clipping boundary or be fully out of it.
            // Thus the clipped contour is supposed to become open.
            assert(front.x() != back.x() || front.y() != back.y());
            if (front.x() != back.x() || front.y() != back.y()) {
                // Look up the ends in "splits", possibly join the contours.
                // "splits" maps into the other piece connected to the same end point.
                auto find_end = [&splits](const ClipperLib_Z::IntPoint &pt) -> std::pair<ClipperLib_Z::IntPoint, int>* {
                    auto it = std::lower_bound(splits.begin(), splits.end(), pt,
                        [](const auto &l, const auto &r){ return ClipperZUtils::zpoint_lower(l.first, r); });
                    return it != splits.end() && it->first == pt ? &(*it) : nullptr;
                };
                auto *end = find_end(front);
                bool  end_front = true;
                if (! end) {
                    end_front = false;
                    end = find_end(back);
                }
                if (end) {
                    // This segment ends at a split point of the source closed contour before clipping.
                    if (end->second == -1) {
                        // Open end was found, not matched yet.
                        end->second = int(it_path - paths.begin());
                    } else {
                        // Open end was found and matched with end->second
                        ClipperLib_Z::Path &other_path = paths[end->second];
                        polylines_merge(other_path, other_path.front() == end->first, std::move(path), end_front);
                        if (std::next(it_path) == paths.end()) {
                            paths.pop_back();
                            break;
                        }
                        path = std::move(paths.back());
                        paths.pop_back();
                        merged = true;
                    }
                }
            }
        }
        if (! merged)
            ++ it_path;
    }
}

static ClipperLib::Paths wavefront_initial(ClipperLib::ClipperOffset &co, const ClipperLib::Paths &polylines, float offset)
{
    ClipperLib::Paths out;
    out.reserve(polylines.size());
    ClipperLib::Paths out_this;
    for (const ClipperLib::Path &path : polylines) {
        co.Clear();
        co.AddPath(path, jtRound, ClipperLib::etOpenRound);
        co.Execute(out_this, offset);
        append(out, std::move(out_this));
    }
    return out;
}

// Input polygons may consist of multiple expolygons, even nested expolygons.
// After inflation some polygons may thus overlap, however the overlap is being resolved during the successive
// clipping operation, thus it is not being done here.
static ClipperLib::Paths wavefront_step(ClipperLib::ClipperOffset &co, const ClipperLib::Paths &polygons, float offset)
{
    ClipperLib::Paths out;
    out.reserve(polygons.size());
    ClipperLib::Paths out_this;
    for (const ClipperLib::Path &polygon : polygons) {
        co.Clear();
        // Execute reorients the contours so that the outer most contour has a positive area. Thus the output
        // contours will be CCW oriented even though the input paths are CW oriented.
        // Offset is applied after contour reorientation, thus the signum of the offset value is reversed.
        co.AddPath(polygon, jtRound, ClipperLib::etClosedPolygon);
        bool ccw = ClipperLib::Orientation(polygon);
        co.Execute(out_this, ccw ? offset : - offset);
        if (! ccw) {
            // Reverse the resulting contours.
            for (ClipperLib::Path &path : out_this)
                std::reverse(path.begin(), path.end());
        }
        append(out, std::move(out_this));
    }
    return out;
}

static ClipperLib::Paths wavefront_clip(const ClipperLib::Paths &wavefront, const Polygons &clipping)
{
    ClipperLib::Clipper clipper;
    clipper.AddPaths(wavefront, ClipperLib::ptSubject, true);
    clipper.AddPaths(ClipperUtils::PolygonsProvider(clipping),  ClipperLib::ptClip, true);
    ClipperLib::Paths out;
    clipper.Execute(ClipperLib::ctIntersection, out, ClipperLib::pftPositive, ClipperLib::pftPositive);
    return out;
}

static Polygons propagate_wave_from_boundary(
    ClipperLib::ClipperOffset   &co,
    // Seed of the wave: Open polylines very close to the boundary.
    const ClipperLib::Paths     &seed,
    // Boundary inside which the waveform will propagate.
    const ExPolygon             &boundary,
    // How much to inflate the seed lines to produce the first wave area.
    const float                  initial_step,
    // How much to inflate the first wave area and the successive wave areas in each step.
    const float                  other_step,
    // Number of inflate steps after the initial step.
    const size_t                 num_other_steps,
    // Maximum inflation of seed contours over the boundary. Used to trim boundary to speed up
    // clipping during wave propagation.
    const float                  max_inflation)
{
    assert(! seed.empty() && seed.front().size() >= 2);
    Polygons clipping = ClipperUtils::clip_clipper_polygons_with_subject_bbox(boundary, get_extents<true>(seed).inflated(max_inflation));
    ClipperLib::Paths polygons = wavefront_clip(wavefront_initial(co, seed, initial_step), clipping);
    // Now offset the remaining 
    for (size_t ioffset = 0; ioffset < num_other_steps; ++ ioffset)
        polygons = wavefront_clip(wavefront_step(co, polygons, other_step), clipping);
    return to_polygons(polygons);
}

// Calculating radius discretization according to ClipperLib offsetter code, see void ClipperOffset::DoOffset(double delta)
inline double clipper_round_offset_error(double offset, double arc_tolerance)
{
    static constexpr const double def_arc_tolerance = 0.25;
    const double y =
        arc_tolerance <= 0 ?
            def_arc_tolerance :
            arc_tolerance > offset * def_arc_tolerance ?
                offset * def_arc_tolerance :
                arc_tolerance;
    double steps = std::min(M_PI / std::acos(1. - y / offset), offset * M_PI);
    return offset * (1. - cos(M_PI / steps));
}

// Returns regions per source ExPolygon expanded into boundary.
std::vector<Polygons> expand_expolygons(
    // Source regions that are supposed to touch the boundary.
    // Boundaries of source regions touching the "boundary" regions will be expanded into the "boundary" region.
    const ExPolygons    &src,
    const ExPolygons    &boundary,
    // Scaled expansion value
    float                full_expansion,
    // Expand by waves of expansion_step size (expansion_step is scaled).
    float                expansion_step,
    // Don't take more than max_nr_steps for small expansion_step.
    size_t               max_nr_expansion_steps)
{
    assert(full_expansion > 0);
    assert(expansion_step > 0);
    assert(max_nr_expansion_steps > 0);

    // Initial expansion of src to make the source regions intersect with boundary regions just a bit.
    float                  tiny_expansion;
    // How much to inflate the seed lines to produce the first wave area.
    float                  initial_step;
    // How much to inflate the first wave area and the successive wave areas in each step.
    float                  other_step;
    // Number of inflate steps after the initial step.
    size_t                 num_other_steps;
    // Maximum inflation of seed contours over the boundary. Used to trim boundary to speed up
    // clipping during wave propagation.
    float                  max_inflation;

    // Offsetter to be applied for all inflation waves. Its accuracy is set with the block below.
    ClipperLib::ClipperOffset co;

    {
        // Initial expansion of src to make the source regions intersect with boundary regions just a bit.
        // The expansion should not be too tiny, but also small enough, so the following expansion will
        // compensate for tiny_expansion and bring the wave back to the boundary without producing
        // ugly cusps where it touches the boundary.
        tiny_expansion = std::min(0.25f * full_expansion, scaled<float>(0.05f));
        size_t nsteps = size_t(ceil((full_expansion - tiny_expansion) / expansion_step));
        if (max_nr_expansion_steps > 0)
            nsteps = std::min(nsteps, max_nr_expansion_steps);
        assert(nsteps > 0);
        initial_step = (full_expansion - tiny_expansion) / nsteps;
        if (nsteps > 1 && 0.25 * initial_step < tiny_expansion) {
            // Decrease the step size by lowering number of steps.
            nsteps       = std::max<size_t>(1, (floor((full_expansion - tiny_expansion) / (4. * tiny_expansion))));
            initial_step = (full_expansion - tiny_expansion) / nsteps;
        }
        if (0.25 * initial_step < tiny_expansion || nsteps == 1) {
            tiny_expansion = 0.2f * full_expansion;
            initial_step   = 0.8f * full_expansion;
        }
        other_step = initial_step;
        num_other_steps = nsteps - 1;

        // Accuracy of the offsetter for wave propagation.
        co.ArcTolerance       = float(scale_(0.1));
        co.ShortestEdgeLength = std::abs(initial_step * ClipperOffsetShortestEdgeFactor);

        // Maximum inflation of seed contours over the boundary. Used to trim boundary to speed up
        // clipping during wave propagation. Needs to be in sync with the offsetter accuracy.
        // Clipper positive round offset should rather offset less than more.
        // Still a little bit of additional offset was added.
        max_inflation = (tiny_expansion + nsteps * initial_step) * 1.1;
//                (clipper_round_offset_error(tiny_expansion, co.ArcTolerance) + nsteps * clipper_round_offset_error(initial_step, co.ArcTolerance) * 1.5; // Account for uncertainty 
    }

    using Intersection  = ClipperZUtils::ClipperZIntersectionVisitor::Intersection;
    using Intersections = ClipperZUtils::ClipperZIntersectionVisitor::Intersections;

    ClipperLib_Z::Paths expansion_seeds;
    Intersections       intersections;

    coord_t             idx_boundary_begin = 1;
    coord_t             idx_boundary_end;
    coord_t             idx_src_end;

    {
        ClipperLib_Z::Clipper zclipper;
        ClipperZUtils::ClipperZIntersectionVisitor visitor(intersections);
        zclipper.ZFillFunction(visitor.clipper_callback());
        // as closed contours
        {
            ClipperLib_Z::Paths zboundary = ClipperZUtils::expolygons_to_zpaths(boundary, idx_boundary_begin);
            idx_boundary_end = idx_boundary_begin + coord_t(zboundary.size());
            zclipper.AddPaths(zboundary, ClipperLib_Z::ptClip, true);
        }
        // as open contours
        std::vector<std::pair<ClipperLib_Z::IntPoint, int>> zsrc_splits;
        {
            ClipperLib_Z::Paths zsrc = expolygons_to_zpaths_expanded_opened(src, tiny_expansion, idx_boundary_end);
            zclipper.AddPaths(zsrc, ClipperLib_Z::ptSubject, false);
            idx_src_end = idx_boundary_end + coord_t(zsrc.size());
            zsrc_splits.reserve(zsrc.size());
            for (const ClipperLib_Z::Path &path : zsrc) {
                assert(path.size() >= 2);
                assert(path.front() == path.back());
                zsrc_splits.emplace_back(path.front(), -1);
            }
            std::sort(zsrc_splits.begin(), zsrc_splits.end(), [](const auto &l, const auto &r){ return ClipperZUtils::zpoint_lower(l.first, r.first); });
        }
        ClipperLib_Z::PolyTree polytree;
        zclipper.Execute(ClipperLib_Z::ctIntersection, polytree, ClipperLib_Z::pftNonZero, ClipperLib_Z::pftNonZero);
        ClipperLib_Z::PolyTreeToPaths(std::move(polytree), expansion_seeds);
        merge_splits(expansion_seeds, zsrc_splits);
    }

    // Sort paths into their respective islands.
    // Each src x boundary will be processed (wave expanded) independently.
    // Multiple pieces of a single src may intersect the same boundary.
    struct SeedOrigin {
        int     src;
        int     boundary;
        int     seed;
    };
    std::vector<SeedOrigin> map_seeds;
    map_seeds.reserve(expansion_seeds.size());
    int iseed = 0;
    for (const ClipperLib_Z::Path &path : expansion_seeds) {
        assert(path.size() >= 2);
        const ClipperLib_Z::IntPoint &front = path.front();
        const ClipperLib_Z::IntPoint &back  = path.back();
        // Both ends of a seed segment are supposed to be inside a single boundary expolygon.
        assert(front.z() < 0);
        assert(back.z() < 0);
        const Intersection *intersection = nullptr;
        auto intersection_point_valid = [idx_boundary_end, idx_src_end](const Intersection &is) {
            return is.first >= 1 && is.first < idx_boundary_end &&
                   is.second >= idx_boundary_end && is.second < idx_src_end;
        };
        if (front.z() < 0) {
            const Intersection &is = intersections[- front.z() - 1];
            assert(intersection_point_valid(is));
            if (intersection_point_valid(is))
                intersection = &is;
        }
        if (! intersection && back.z() < 0) {
            const Intersection &is = intersections[- back.z() - 1];
            assert(intersection_point_valid(is));
            if (intersection_point_valid(is))
                intersection = &is;
        }
        if (intersection) {
            // The path intersects the boundary contour at least at one side. 
            map_seeds.push_back({ intersection->second - idx_boundary_end, intersection->first - 1, iseed });
        }
        ++ iseed;
    }
    // Sort the seeds by their intersection boundary and source contour.
    std::sort(map_seeds.begin(), map_seeds.end(), [](const auto &l, const auto &r){ 
        return l.boundary < r.boundary || (l.boundary == r.boundary && l.src < r.src); 
    });

    std::vector<Polygons> out(src.size(), Polygons{});
    ClipperLib::Paths     paths;
    for (auto it_seed = map_seeds.begin(); it_seed != map_seeds.end();) {
        auto it = it_seed;
        paths.clear();
        for (; it != map_seeds.end() && it->boundary == it_seed->boundary && it->src == it_seed->src; ++ it)
            paths.emplace_back(ClipperZUtils::from_zpath(expansion_seeds[it->seed]));
        // Propagate the wavefront while clipping it with the trimmed boundary.
        // Collect the expanded polygons, merge them with the source polygons.
        append(out[it_seed->src], propagate_wave_from_boundary(co, paths, boundary[it_seed->boundary], initial_step, other_step, num_other_steps, max_inflation));
        it_seed = it;
    }

    return out;
}

} // Algorithm
} // Slic3r
