#include "Layer.hpp"
#include <clipper/clipper_z.hpp>
#include "ClipperUtils.hpp"
#include "Print.hpp"
#include "Fill/Fill.hpp"
#include "ShortestPath.hpp"
#include "SVG.hpp"
#include "BoundingBox.hpp"

#include <boost/log/trivial.hpp>

namespace Slic3r {

Layer::~Layer()
{
    this->lower_layer = this->upper_layer = nullptr;
    for (LayerRegion *region : m_regions)
        delete region;
    m_regions.clear();
}

// Test whether whether there are any slices assigned to this layer.
bool Layer::empty() const
{
	for (const LayerRegion *layerm : m_regions)
        if (layerm != nullptr && ! layerm->slices().empty())
            // Non empty layer.
            return false;
    return true;
}

LayerRegion* Layer::add_region(const PrintRegion *print_region)
{
    m_regions.emplace_back(new LayerRegion(this, print_region));
    return m_regions.back();
}

// merge all regions' slices to get islands
void Layer::make_slices()
{
    ExPolygons slices;
    if (m_regions.size() == 1) {
        // optimization: if we only have one region, take its slices
        slices = to_expolygons(m_regions.front()->slices().surfaces);
    } else {
        Polygons slices_p;
        for (LayerRegion *layerm : m_regions)
            polygons_append(slices_p, to_polygons(layerm->slices().surfaces));
        slices = union_safety_offset_ex(slices_p);
    }
    
    this->lslices.clear();
    this->lslices.reserve(slices.size());
    
    // prepare ordering points
    Points ordering_points;
    ordering_points.reserve(slices.size());
    for (const ExPolygon &ex : slices)
        ordering_points.push_back(ex.contour.first_point());
    
    // sort slices
    std::vector<Points::size_type> order = chain_points(ordering_points);
    
    // populate slices vector
    for (size_t i : order)
        this->lslices.emplace_back(std::move(slices[i]));
}

// used by Layer::build_up_down_graph()
[[nodiscard]] static ClipperLib_Z::Paths expolygons_to_zpaths(const ExPolygons &expolygons, coord_t isrc)
{
    size_t num_paths = 0;
    for (const ExPolygon &expolygon : expolygons)
        num_paths += expolygon.num_contours();

    ClipperLib_Z::Paths out;
    out.reserve(num_paths);

    for (const ExPolygon &expolygon : expolygons) {
        for (size_t icontour = 0; icontour < expolygon.num_contours(); ++ icontour) {
            const Polygon &contour = expolygon.contour_or_hole(icontour);
            out.emplace_back();
            ClipperLib_Z::Path &path = out.back();
            path.reserve(contour.size());
            for (const Point &p : contour.points)
                path.push_back({ p.x(), p.y(), isrc });
        }
        ++ isrc;
    }

    return out;
}

// used by Layer::build_up_down_graph()
static void connect_layer_slices(
    Layer                                           &below,
    Layer                                           &above,
    const ClipperLib_Z::PolyTree                    &polytree,
    const std::vector<std::pair<coord_t, coord_t>>  &intersections,
    const coord_t                                    offset_below,
    const coord_t                                    offset_above
#ifndef NDEBUG
    , const coord_t                                    offset_end
#endif // NDEBUG
    )
{
    class Visitor {
    public:
        Visitor(const std::vector<std::pair<coord_t, coord_t>> &intersections, 
            Layer &below, Layer &above, const coord_t offset_below, const coord_t offset_above
#ifndef NDEBUG
            , const coord_t offset_end
#endif // NDEBUG
            ) :
            m_intersections(intersections), m_below(below), m_above(above), m_offset_below(offset_below), m_offset_above(offset_above)
#ifndef NDEBUG
            , m_offset_end(offset_end) 
#endif // NDEBUG
            {}

        void visit(const ClipperLib_Z::PolyNode &polynode)
        {
#ifndef NDEBUG
            auto assert_intersection_valid = [this](int i, int j) {
                assert(i != j);
                if (i > j)
                    std::swap(i, j);
                assert(i >= m_offset_below);
                assert(i < m_offset_above);
                assert(j >= m_offset_above);
                assert(j < m_offset_end);
                return true;
            };
#endif // NDEBUG
            if (polynode.Contour.size() >= 3) {
                // If there is an intersection point, it should indicate which contours (one from layer below, the other from layer above) intersect.
                // Otherwise the contour is fully inside another contour.
                int32_t i = 0, j = 0;
                for (int icontour = 0; icontour <= polynode.ChildCount(); ++ icontour) {
                    const bool                first   = icontour == 0;
                    const ClipperLib_Z::Path &contour = first ? polynode.Contour : polynode.Childs[icontour - 1]->Contour;
                    if (contour.size() >= 3) {
                        if (first) {
                            i = contour.front().z();
                            j = i;
                            if (i < 0) {
                                std::tie(i, j) = m_intersections[-i - 1];
                                assert(assert_intersection_valid(i, j));
                                goto end;
                            }
                        }
                        for (const ClipperLib_Z::IntPoint& pt : contour) {
                            j = pt.z();
                            if (j < 0) {
                                std::tie(i, j) = m_intersections[-j - 1];
                                assert(assert_intersection_valid(i, j));
                                goto end;
                            }
                            else if (i != j)
                                goto end;
                        }
                    }
                }
            end:
                bool found = false;
                if (i == j) {
                    // The contour is completely inside another contour.
                    Point pt(polynode.Contour.front().x(), polynode.Contour.front().y());
                    if (i < m_offset_above) {
                        // Index of an island below. Look-it up in the island above.
                        assert(i >= m_offset_below);
                        i -= m_offset_below;
                        for (int l = int(m_above.lslices_ex.size()) - 1; l >= 0; -- l) {
                            LayerSlice &lslice = m_above.lslices_ex[l];
                            if (lslice.bbox.contains(pt) && m_above.lslices[l].contains(pt)) {
                                found = true;
                                j = l;
                                assert(i >= 0 && i < m_below.lslices_ex.size());
                                assert(j >= 0 && j < m_above.lslices_ex.size());
                                break;
                            }
                        }
                    } else {
                        // Index of an island above. Look-it up in the island below.
                        assert(j < m_offset_end);
                        j -= m_offset_above;
                        for (int l = int(m_below.lslices_ex.size()) - 1; l >= 0; -- l) {
                            LayerSlice &lslice = m_below.lslices_ex[l];
                            if (lslice.bbox.contains(pt) && m_below.lslices[l].contains(pt)) {
                                found = true;
                                i = l;
                                assert(i >= 0 && i < m_below.lslices_ex.size());
                                assert(j >= 0 && j < m_above.lslices_ex.size());
                                break;
                            }
                        }
                    }
                } else {
                    assert(assert_intersection_valid(i, j));
                    if (i > j)
                        std::swap(i, j);
                    i -= m_offset_below;
                    j -= m_offset_above;
                    assert(i >= 0 && i < m_below.lslices_ex.size());
                    assert(j >= 0 && j < m_above.lslices_ex.size());
                    found = true;
                }
                if (found) {
                    // Subtract area of holes from the area of outer contour.
                    double area = ClipperLib_Z::Area(polynode.Contour);
                    for (int icontour = 0; icontour < polynode.ChildCount(); ++ icontour)
                        area -= ClipperLib_Z::Area(polynode.Childs[icontour]->Contour);
                    // Store the links and area into the contours.
                    LayerSlice::Links &links_below = m_below.lslices_ex[i].overlaps_above;
                    LayerSlice::Links &links_above = m_above.lslices_ex[j].overlaps_below;
                    LayerSlice::Link key{ j };
                    auto it_below = std::lower_bound(links_below.begin(), links_below.end(), key, [](auto &l, auto &r){ return l.slice_idx < r.slice_idx; });
                    if (it_below != links_below.end() && it_below->slice_idx == j) {
                        it_below->area += area;
                    } else {
                        auto it_above = std::lower_bound(links_above.begin(), links_above.end(), key, [](auto &l, auto &r){ return l.slice_idx < r.slice_idx; });
                        if (it_above != links_above.end() && it_above->slice_idx == i) {
                            it_above->area += area;
                        } else {
                            // Insert into one of the two vectors.
                            bool take_below = false;
                            if (links_below.size() < LayerSlice::LinksStaticSize)
                                take_below = false;
                            else if (links_above.size() >= LayerSlice::LinksStaticSize) {
                                size_t shift_below = links_below.end() - it_below;
                                size_t shift_above = links_above.end() - it_above;
                                take_below = shift_below < shift_above;
                            }
                            if (take_below)
                                links_below.insert(it_below, { j, float(area) });
                            else
                                links_above.insert(it_above, { i, float(area) });
                        }
                    }
                }
            }
            for (int i = 0; i < polynode.ChildCount(); ++ i)
                for (int j = 0; j < polynode.Childs[i]->ChildCount(); ++ j)
                    this->visit(*polynode.Childs[i]->Childs[j]);
        }

    private:
        const std::vector<std::pair<coord_t, coord_t>> &m_intersections;
        Layer                                          &m_below;
        Layer                                          &m_above;
        const coord_t                                   m_offset_below;
        const coord_t                                   m_offset_above;
#ifndef NDEBUG
        const coord_t                                   m_offset_end;
#endif // NDEBUG
    } visitor(intersections, below, above, offset_below, offset_above
#ifndef NDEBUG
        , offset_end
#endif // NDEBUG
    );

    for (int i = 0; i < polytree.ChildCount(); ++ i)
        visitor.visit(*polytree.Childs[i]);

#ifndef NDEBUG
    // Verify that only one directional link is stored: either from bottom slice up or from upper slice down.
    for (int32_t islice = 0; islice < below.lslices_ex.size(); ++ islice) {
        LayerSlice::Links &links1 = below.lslices_ex[islice].overlaps_above;
        for (LayerSlice::Link &link1 : links1) {
            LayerSlice::Links &links2 = above.lslices_ex[link1.slice_idx].overlaps_below;
            assert(! std::binary_search(links2.begin(), links2.end(), link1, [](auto &l, auto &r){ return l.slice_idx < r.slice_idx; }));
        }
    }
    for (int32_t islice = 0; islice < above.lslices_ex.size(); ++ islice) {
        LayerSlice::Links &links1 = above.lslices_ex[islice].overlaps_below;
        for (LayerSlice::Link &link1 : links1) {
            LayerSlice::Links &links2 = below.lslices_ex[link1.slice_idx].overlaps_above;
            assert(! std::binary_search(links2.begin(), links2.end(), link1, [](auto &l, auto &r){ return l.slice_idx < r.slice_idx; }));
        }
    }
#endif // NDEBUG

    // Scatter the links, but don't sort them yet.
    for (int32_t islice = 0; islice < below.lslices_ex.size(); ++ islice)
        for (LayerSlice::Link &link : below.lslices_ex[islice].overlaps_above)
            above.lslices_ex[link.slice_idx].overlaps_below.push_back({ islice, link.area });
    for (int32_t islice = 0; islice < above.lslices_ex.size(); ++ islice)
        for (LayerSlice::Link &link : above.lslices_ex[islice].overlaps_below)
            below.lslices_ex[link.slice_idx].overlaps_above.push_back({ islice, link.area });
    // Sort the links.
    for (LayerSlice &lslice : below.lslices_ex)
        std::sort(lslice.overlaps_above.begin(), lslice.overlaps_above.end(), [](const LayerSlice::Link &l, const LayerSlice::Link &r){ return l.slice_idx < r.slice_idx; });
    for (LayerSlice &lslice : above.lslices_ex)
        std::sort(lslice.overlaps_below.begin(), lslice.overlaps_below.end(), [](const LayerSlice::Link &l, const LayerSlice::Link &r){ return l.slice_idx < r.slice_idx; });
}

void Layer::build_up_down_graph(Layer& below, Layer& above)
{
    coord_t             paths_below_offset = 0;
    ClipperLib_Z::Paths paths_below = expolygons_to_zpaths(below.lslices, paths_below_offset);
    coord_t             paths_above_offset = paths_below_offset + coord_t(below.lslices.size());
    ClipperLib_Z::Paths paths_above = expolygons_to_zpaths(above.lslices, paths_above_offset);
#ifndef NDEBUG
    coord_t             paths_end = paths_above_offset + coord_t(above.lslices.size());
#endif // NDEBUG

    class ZFill {
    public:
        ZFill() = default;
        void reset() { m_intersections.clear(); }
        void operator()(
            const ClipperLib_Z::IntPoint& e1bot, const ClipperLib_Z::IntPoint& e1top,
            const ClipperLib_Z::IntPoint& e2bot, const ClipperLib_Z::IntPoint& e2top,
            ClipperLib_Z::IntPoint& pt) {
            coord_t srcs[4]{ e1bot.z(), e1top.z(), e2bot.z(), e2top.z() };
            coord_t* begin = srcs;
            coord_t* end = srcs + 4;
            std::sort(begin, end);
            end = std::unique(begin, end);
            assert(begin + 2 == end);
            if (begin + 1 == end)
                pt.z() = *begin;
            else if (begin + 2 <= end) {
                // store a -1 based negative index into the "intersections" vector here.
                m_intersections.emplace_back(srcs[0], srcs[1]);
                pt.z() = -coord_t(m_intersections.size());
            }
        }
        const std::vector<std::pair<coord_t, coord_t>>& intersections() const { return m_intersections; }

    private:
        std::vector<std::pair<coord_t, coord_t>> m_intersections;
    } zfill;

    ClipperLib_Z::Clipper  clipper;
    ClipperLib_Z::PolyTree result;
    clipper.ZFillFunction(
        [&zfill](const ClipperLib_Z::IntPoint &e1bot, const ClipperLib_Z::IntPoint &e1top, 
                 const ClipperLib_Z::IntPoint &e2bot, const ClipperLib_Z::IntPoint &e2top, ClipperLib_Z::IntPoint &pt)
        { return zfill(e1bot, e1top, e2bot, e2top, pt); });
    clipper.AddPaths(paths_below, ClipperLib_Z::ptSubject, true);
    clipper.AddPaths(paths_above, ClipperLib_Z::ptClip, true);
    clipper.Execute(ClipperLib_Z::ctIntersection, result, ClipperLib_Z::pftNonZero, ClipperLib_Z::pftNonZero);

    connect_layer_slices(below, above, result, zfill.intersections(), paths_below_offset, paths_above_offset
#ifndef NDEBUG
        , paths_end
#endif // NDEBUG
        );
}

static inline bool layer_needs_raw_backup(const Layer *layer)
{
    return ! (layer->regions().size() == 1 && (layer->id() > 0 || layer->object()->config().elefant_foot_compensation.value == 0));
}

void Layer::backup_untyped_slices()
{
    if (layer_needs_raw_backup(this)) {
        for (LayerRegion *layerm : m_regions)
            layerm->m_raw_slices = to_expolygons(layerm->slices().surfaces);
    } else {
        assert(m_regions.size() == 1);
        m_regions.front()->m_raw_slices.clear();
    }
}

void Layer::restore_untyped_slices()
{
    if (layer_needs_raw_backup(this)) {
        for (LayerRegion *layerm : m_regions)
            layerm->m_slices.set(layerm->m_raw_slices, stInternal);
    } else {
        assert(m_regions.size() == 1);
        m_regions.front()->m_slices.set(this->lslices, stInternal);
    }
}

// Similar to Layer::restore_untyped_slices()
// To improve robustness of detect_surfaces_type() when reslicing (working with typed slices), see GH issue #7442.
// Only resetting layerm->slices if Slice::extra_perimeters is always zero or it will not be used anymore
// after the perimeter generator.
void Layer::restore_untyped_slices_no_extra_perimeters()
{
    if (layer_needs_raw_backup(this)) {
        for (LayerRegion *layerm : m_regions)
        	if (! layerm->region().config().extra_perimeters.value)
            	layerm->m_slices.set(layerm->m_raw_slices, stInternal);
    } else {
    	assert(m_regions.size() == 1);
    	LayerRegion *layerm = m_regions.front();
    	// This optimization is correct, as extra_perimeters are only reused by prepare_infill() with multi-regions.
        //if (! layerm->region().config().extra_perimeters.value)
        	layerm->m_slices.set(this->lslices, stInternal);
    }
}

ExPolygons Layer::merged(float offset_scaled) const
{
	assert(offset_scaled >= 0.f);
    // If no offset is set, apply EPSILON offset before union, and revert it afterwards.
	float offset_scaled2 = 0;
	if (offset_scaled == 0.f) {
		offset_scaled  = float(  EPSILON);
		offset_scaled2 = float(- EPSILON);
    }
    Polygons polygons;
	for (LayerRegion *layerm : m_regions) {
		const PrintRegionConfig &config = layerm->region().config();
		// Our users learned to bend Slic3r to produce empty volumes to act as subtracters. Only add the region if it is non-empty.
		if (config.bottom_solid_layers > 0 || config.top_solid_layers > 0 || config.fill_density > 0. || config.perimeters > 0)
			append(polygons, offset(layerm->slices().surfaces, offset_scaled));
	}
    ExPolygons out = union_ex(polygons);
	if (offset_scaled2 != 0.f)
		out = offset_ex(out, offset_scaled2);
    return out;
}

// Here the perimeters are created cummulatively for all layer regions sharing the same parameters influencing the perimeters.
// The perimeter paths and the thin fills (ExtrusionEntityCollection) are assigned to the first compatible layer region.
// The resulting fill surface is split back among the originating regions.
void Layer::make_perimeters()
{
    BOOST_LOG_TRIVIAL(trace) << "Generating perimeters for layer " << this->id();
    
    // keep track of regions whose perimeters we have already generated
    std::vector<unsigned char>                              done(m_regions.size(), false);
    std::vector<uint32_t>                                   layer_region_ids;
    std::vector<std::pair<ExtrusionRange, ExtrusionRange>>  perimeter_and_gapfill_ranges;
    ExPolygons                                              fill_expolygons;
    std::vector<ExPolygonRange>                             fill_expolygons_ranges;
    SurfacesPtr                                             surfaces_to_merge;
    SurfacesPtr                                             surfaces_to_merge_temp;

    auto layer_region_reset_perimeters = [](LayerRegion &layerm) {
        layerm.m_perimeters.clear();
        layerm.m_fills.clear();
        layerm.m_thin_fills.clear();
        layerm.m_fill_expolygons.clear();
        layerm.m_fill_expolygons_bboxes.clear();
        layerm.m_fill_expolygons_composite.clear();
        layerm.m_fill_expolygons_composite_bboxes.clear();
    };

    // Remove layer islands, remove references to perimeters and fills from these layer islands to LayerRegion ExtrusionEntities.
    for (LayerSlice &lslice : this->lslices_ex)
        lslice.islands.clear();

    for (LayerRegionPtrs::iterator layerm = m_regions.begin(); layerm != m_regions.end(); ++ layerm)
        if (size_t region_id = layerm - m_regions.begin(); ! done[region_id]) {
            layer_region_reset_perimeters(**layerm);
            if (! (*layerm)->slices().empty()) {
    	        BOOST_LOG_TRIVIAL(trace) << "Generating perimeters for layer " << this->id() << ", region " << region_id;
    	        done[region_id] = true;
    	        const PrintRegionConfig &config = (*layerm)->region().config();
    	        
                perimeter_and_gapfill_ranges.clear();
                fill_expolygons.clear();
                fill_expolygons_ranges.clear();
                surfaces_to_merge.clear();

    	        // find compatible regions
                layer_region_ids.clear();
    	        layer_region_ids.push_back(region_id);
    	        for (LayerRegionPtrs::const_iterator it = layerm + 1; it != m_regions.end(); ++it)
    	            if (! (*it)->slices().empty()) {
    		            LayerRegion* other_layerm = *it;
    		            const PrintRegionConfig &other_config = other_layerm->region().config();
    		            if (config.perimeter_extruder             == other_config.perimeter_extruder
    		                && config.perimeters                  == other_config.perimeters
    		                && config.perimeter_speed             == other_config.perimeter_speed
    		                && config.external_perimeter_speed    == other_config.external_perimeter_speed
    		                && (config.gap_fill_enabled ? config.gap_fill_speed.value : 0.) == 
                               (other_config.gap_fill_enabled ? other_config.gap_fill_speed.value : 0.)
    		                && config.overhangs                   == other_config.overhangs
    		                && config.opt_serialize("perimeter_extrusion_width") == other_config.opt_serialize("perimeter_extrusion_width")
    		                && config.thin_walls                  == other_config.thin_walls
    		                && config.external_perimeters_first   == other_config.external_perimeters_first
    		                && config.infill_overlap              == other_config.infill_overlap
                            && config.fuzzy_skin                  == other_config.fuzzy_skin
                            && config.fuzzy_skin_thickness        == other_config.fuzzy_skin_thickness
                            && config.fuzzy_skin_point_dist       == other_config.fuzzy_skin_point_dist)
    		            {
                            layer_region_reset_perimeters(*other_layerm);
    		                layer_region_ids.push_back(it - m_regions.begin());
    		                done[it - m_regions.begin()] = true;
    		            }
    		        }

    	        if (layer_region_ids.size() == 1) {  // optimization
    	            (*layerm)->make_perimeters((*layerm)->slices(), perimeter_and_gapfill_ranges, fill_expolygons, fill_expolygons_ranges);
                    this->sort_perimeters_into_islands((*layerm)->slices(), region_id, perimeter_and_gapfill_ranges, std::move(fill_expolygons), fill_expolygons_ranges, layer_region_ids);
    	        } else {
    	            SurfaceCollection new_slices;
    	            // Use the region with highest infill rate, as the make_perimeters() function below decides on the gap fill based on the infill existence.
    	            LayerRegion *layerm_config = m_regions[layer_region_ids.front()];
    	            {
    	                // Merge slices (surfaces) according to number of extra perimeters.
    	                for (uint32_t region_id : layer_region_ids) {
                            LayerRegion &layerm = *m_regions[region_id];
    	                    for (const Surface &surface : layerm.slices())
                                surfaces_to_merge.emplace_back(&surface);
    	                    if (layerm.region().config().fill_density > layerm_config->region().config().fill_density)
    	                    	layerm_config = &layerm;
    	                }
                        std::sort(surfaces_to_merge.begin(), surfaces_to_merge.end(), [](const Surface *l, const Surface *r){ return l->extra_perimeters < r->extra_perimeters; });
                        for (size_t i = 0; i < surfaces_to_merge.size();) {
                            size_t j = i;
                            const Surface &first = *surfaces_to_merge[i];
                            size_t extra_perimeters = first.extra_perimeters;
                            for (; j < surfaces_to_merge.size() && surfaces_to_merge[j]->extra_perimeters == extra_perimeters; ++ j) ;
                            if (i + 1 == j)
                                // Nothing to merge, just copy.
                                new_slices.surfaces.emplace_back(*surfaces_to_merge[i]);
                            else {
                                surfaces_to_merge_temp.assign(surfaces_to_merge.begin() + i, surfaces_to_merge.begin() + j);
                                new_slices.append(offset_ex(surfaces_to_merge_temp, ClipperSafetyOffset), first);
                            }
                            i = j;
                        }
    	            }
    	            // make perimeters
    	            layerm_config->make_perimeters(new_slices, perimeter_and_gapfill_ranges, fill_expolygons, fill_expolygons_ranges);
                    this->sort_perimeters_into_islands(new_slices, region_id, perimeter_and_gapfill_ranges, std::move(fill_expolygons), fill_expolygons_ranges, layer_region_ids);
    	        }
    	    }
        }
    BOOST_LOG_TRIVIAL(trace) << "Generating perimeters for layer " << this->id() << " - Done";
}

void Layer::sort_perimeters_into_islands(
    // Slices for which perimeters and fill_expolygons were just created.
    // The slices may have been created by merging multiple source slices with the same perimeter parameters.
    const SurfaceCollection                                         &slices,
    // Region where the perimeters, gap fills and fill expolygons are stored.
    const uint32_t                                                   region_id,
    // Perimeters and gap fills produced by the perimeter generator for the slices,
    // sorted by the source slices.
    const std::vector<std::pair<ExtrusionRange, ExtrusionRange>>    &perimeter_and_gapfill_ranges,
    // Fill expolygons produced for all source slices above.
    ExPolygons                                                      &&fill_expolygons,
    // Fill expolygon ranges sorted by the source slices.
    const std::vector<ExPolygonRange>                               &fill_expolygons_ranges,
    // If the current layer consists of multiple regions, then the fill_expolygons above are split by the source LayerRegion surfaces.
    const std::vector<uint32_t>                                     &layer_region_ids)
{
    LayerRegion &this_layer_region = *m_regions[region_id];

    // Bounding boxes of fill_expolygons.
    BoundingBoxes fill_expolygons_bboxes;
    fill_expolygons_bboxes.reserve(fill_expolygons.size());
    for (const ExPolygon &expolygon : fill_expolygons)
        fill_expolygons_bboxes.emplace_back(get_extents(expolygon));


    // Take one sample point for each source slice, to be used to sort source slices into layer slices.
    // source slice index + its sample.
    std::vector<std::pair<uint32_t, Point>> perimeter_slices_queue;
    perimeter_slices_queue.reserve(slices.size());
    for (uint32_t islice = 0; islice < uint32_t(slices.size()); ++ islice) {
        const std::pair<ExtrusionRange, ExtrusionRange> &extrusions = perimeter_and_gapfill_ranges[islice];
        Point sample;
        bool  sample_set = false;
        if (! extrusions.first.empty()) {
            sample     = this_layer_region.perimeters().entities[*extrusions.first.begin()]->first_point();
            sample_set = true;
        } else if (! extrusions.second.empty()) {
            sample     = this_layer_region.thin_fills().entities[*extrusions.second.begin()]->first_point();
            sample_set = true;
        } else {
            for (uint32_t iexpoly : fill_expolygons_ranges[islice])
                if (const ExPolygon &expoly = fill_expolygons[iexpoly]; ! expoly.empty()) {
                    sample     = expoly.contour.points.front();
                    sample_set = true;
                    break;
                }
        }
        // There may be a valid empty island.
        // assert(sample_set);
        if (sample_set)
            perimeter_slices_queue.emplace_back(islice, sample);
    }

    // Map of source fill_expolygon into region and fill_expolygon of that region.
    // -1: not set
    std::vector<std::pair<int, int>> map_expolygon_to_region_and_fill;
    const bool                       has_multiple_regions = layer_region_ids.size() > 1;
    assert(has_multiple_regions || layer_region_ids.size() == 1);
    // assign fill_surfaces to each layer
    if (! fill_expolygons.empty()) {
        if (has_multiple_regions) {
            // Sort the bounding boxes lexicographically.
            std::vector<uint32_t> fill_expolygons_bboxes_sorted(fill_expolygons_bboxes.size());
            std::iota(fill_expolygons_bboxes_sorted.begin(), fill_expolygons_bboxes_sorted.end(), 0);
            std::sort(fill_expolygons_bboxes_sorted.begin(), fill_expolygons_bboxes_sorted.end(), [&fill_expolygons_bboxes](uint32_t lhs, uint32_t rhs){
                const BoundingBox &bbl = fill_expolygons_bboxes[lhs];
                const BoundingBox &bbr = fill_expolygons_bboxes[rhs];
                return bbl.min < bbr.min || (bbl.min == bbr.min && bbl.max < bbr.max);
            });
            map_expolygon_to_region_and_fill.assign(fill_expolygons.size(), std::make_pair(-1, -1));
            for (uint32_t region_idx : layer_region_ids) {
                LayerRegion &l = *m_regions[region_idx];
                l.m_fill_expolygons = intersection_ex(l.slices().surfaces, fill_expolygons);
                l.m_fill_expolygons_bboxes.reserve(l.fill_expolygons().size());
                for (const ExPolygon &expolygon : l.fill_expolygons()) {
                    BoundingBox bbox = get_extents(expolygon);
                    l.m_fill_expolygons_bboxes.emplace_back(bbox);
                    auto it_bbox = std::lower_bound(fill_expolygons_bboxes_sorted.begin(), fill_expolygons_bboxes_sorted.end(), bbox, [&fill_expolygons_bboxes](uint32_t lhs, const BoundingBox &bbr){
                        const BoundingBox &bbl = fill_expolygons_bboxes[lhs];
                        return bbl.min < bbr.min || (bbl.min == bbr.min && bbl.max < bbr.max);
                    });
                    if (it_bbox != fill_expolygons_bboxes_sorted.end())
                        if (uint32_t fill_id = *it_bbox; fill_expolygons_bboxes[fill_id] == bbox) {
                            // With a very high probability the two expolygons match exactly. Confirm that.
                            if (expolygons_match(expolygon, fill_expolygons[fill_id])) {
                                std::pair<int, int> &ref = map_expolygon_to_region_and_fill[fill_id];
                                // Only one expolygon produced by intersection with LayerRegion surface may match an expolygon of fill_expolygons.
                                assert(ref.first == -1);
                                ref.first  = region_idx;
                                ref.second = int(&expolygon - l.fill_expolygons().data());
                            }
                        }
                }
            }
        } else {
            this_layer_region.m_fill_expolygons        = std::move(fill_expolygons);
            this_layer_region.m_fill_expolygons_bboxes = std::move(fill_expolygons_bboxes);
        }
    }

    // Sort perimeter extrusions, thin fill extrusions and fill expolygons into islands.
    std::vector<uint32_t> region_fill_sorted_last;
    auto insert_into_island = [
        // Region where the perimeters, gap fills and fill expolygons are stored.
        region_id, 
        // Whether there are infills with different regions generated for this LayerSlice.
        has_multiple_regions,
        // Perimeters and gap fills to be sorted into islands.
        &perimeter_and_gapfill_ranges,
        // Infill regions to be sorted into islands.
        &fill_expolygons, &fill_expolygons_bboxes, &fill_expolygons_ranges,
        // Mapping of fill_expolygon to region and its infill.
        &map_expolygon_to_region_and_fill,
        // Output
        &regions = m_regions, &lslices_ex = this->lslices_ex,
        // fill_expolygons and fill_expolygons_bboxes need to be sorted into contiguous sequence by island,
        // thus region_fill_sorted_last contains last fill_expolygon processed (meaning sorted).
        &region_fill_sorted_last]
        (int lslice_idx, int source_slice_idx) {
        lslices_ex[lslice_idx].islands.push_back({});
        LayerIsland &island = lslices_ex[lslice_idx].islands.back();
        island.perimeters = LayerExtrusionRange(region_id, perimeter_and_gapfill_ranges[source_slice_idx].first);
        island.thin_fills = perimeter_and_gapfill_ranges[source_slice_idx].second;
        if (ExPolygonRange fill_range = fill_expolygons_ranges[source_slice_idx]; ! fill_range.empty()) {
            if (has_multiple_regions) {
                // Check whether the fill expolygons of this island were split into multiple regions.
                island.fill_region_id = LayerIsland::fill_region_composite_id;
                for (uint32_t fill_idx : fill_range) {
                    const std::pair<int, int> &kvp = map_expolygon_to_region_and_fill[fill_idx];
                    if (kvp.first == -1 || (island.fill_region_id != -1 && island.fill_region_id != kvp.second)) {
                        island.fill_region_id = LayerIsland::fill_region_composite_id;
                        break;
                    } else
                        island.fill_region_id = kvp.second;
                }
                if (island.fill_expolygons_composite()) {
                    // They were split, thus store the unsplit "composite" expolygons into the region of perimeters.
                    LayerRegion &this_layer_region = *regions[region_id];
                    auto begin = uint32_t(this_layer_region.fill_expolygons_composite().size());
                    this_layer_region.m_fill_expolygons_composite.reserve(this_layer_region.fill_expolygons_composite().size() + fill_range.size());
                    std::move(fill_expolygons.begin() + *fill_range.begin(), fill_expolygons.begin() + *fill_range.end(), std::back_inserter(this_layer_region.m_fill_expolygons_composite));
                    this_layer_region.m_fill_expolygons_composite_bboxes.insert(this_layer_region.m_fill_expolygons_composite_bboxes.end(), 
                        fill_expolygons_bboxes.begin() + *fill_range.begin(), fill_expolygons_bboxes.begin() + *fill_range.end());
                    island.fill_expolygons = ExPolygonRange(begin, uint32_t(this_layer_region.fill_expolygons_composite().size()));
                } else {
                    if (region_fill_sorted_last.empty())
                        region_fill_sorted_last.assign(regions.size(), 0);
                    uint32_t &last = region_fill_sorted_last[island.fill_region_id];
                    // They were not split and they belong to the same region.
                    // Sort the region m_fill_expolygons to a continuous span.
                    uint32_t begin = last;
                    LayerRegion &layerm = *regions[island.fill_region_id];
                    for (uint32_t fill_id : fill_range) {
                        uint32_t region_fill_id = map_expolygon_to_region_and_fill[fill_id].second;
                        assert(region_fill_id >= last);
                        if (region_fill_id > last) {
                            std::swap(layerm.m_fill_expolygons[region_fill_id], layerm.m_fill_expolygons[last]);
                            std::swap(layerm.m_fill_expolygons_bboxes[region_fill_id], layerm.m_fill_expolygons_bboxes[last]);
                        }
                        ++ last;
                    }
                    island.fill_expolygons = ExPolygonRange(begin, last);
                }
            } else {
                // Layer island is made of one fill region only.
                island.fill_expolygons = fill_range;
                island.fill_region_id  = region_id;
            }
        }
    };

    // First sort into islands using exact fit.
    // Traverse the slices in an increasing order of bounding box size, so that the islands inside another islands are tested first,
    // so we can just test a point inside ExPolygon::contour and we may skip testing the holes.
    auto point_inside_surface = [&lslices = this->lslices, &lslices_ex = this->lslices_ex](size_t lslice_idx, const Point &point) {
        const BoundingBox &bbox = lslices_ex[lslice_idx].bbox;
        return point.x() >= bbox.min.x() && point.x() < bbox.max.x() &&
               point.y() >= bbox.min.y() && point.y() < bbox.max.y() &&
               // Exact match: Don't just test whether a point is inside the outer contour of an island,
               // test also whether the point is not inside some hole of the same expolygon.
               // This is unfortunatelly necessary because the point may be inside an expolygon of one of this expolygon's hole
               // and missed due to numerical issues.
               lslices[lslice_idx].contains(point);
    };
    for (int lslice_idx = int(lslices_ex.size()) - 1; lslice_idx >= 0 && ! perimeter_slices_queue.empty(); -- lslice_idx)
        for (auto it_source_slice = perimeter_slices_queue.begin(); it_source_slice != perimeter_slices_queue.end(); ++ it_source_slice)
            if (point_inside_surface(lslice_idx, it_source_slice->second)) {
                insert_into_island(lslice_idx, it_source_slice->first);
                if (std::next(it_source_slice) != perimeter_slices_queue.end())
                    // Remove the current slice & point pair from the queue.
                    *it_source_slice = perimeter_slices_queue.back();
                perimeter_slices_queue.pop_back();
                break;
            }
    // If anything fails to be sorted in using exact fit, try to find a closest island.
    auto point_inside_surface_dist2 =
        [&lslices = this->lslices, &lslices_ex = this->lslices_ex, bbox_eps = scaled<coord_t>(this->object()->print()->config().gcode_resolution.value) + SCALED_EPSILON]
        (const size_t lslice_idx, const Point &point) {
        const BoundingBox &bbox = lslices_ex[lslice_idx].bbox;
        return 
            point.x() < bbox.min.x() - bbox_eps || point.x() > bbox.max.x() + bbox_eps ||
            point.y() < bbox.min.y() - bbox_eps || point.y() > bbox.max.y() + bbox_eps ?
            std::numeric_limits<double>::max() :
            (lslices[lslice_idx].point_projection(point) - point).cast<double>().squaredNorm();
    };
    for (auto it_source_slice  = perimeter_slices_queue.begin(); it_source_slice != perimeter_slices_queue.end(); ++ it_source_slice) {
        double d2min = std::numeric_limits<double>::max();
        int    lslice_idx_min = -1;
        for (int lslice_idx = int(lslices_ex.size()) - 1; lslice_idx >= 0; -- lslice_idx)
            if (double d2 = point_inside_surface_dist2(lslice_idx, it_source_slice->second); d2 < d2min) {
                d2min = d2;
                lslice_idx_min = lslice_idx;
            }
        assert(lslice_idx_min != -1);
        insert_into_island(lslice_idx_min, it_source_slice->first);
    }
}

void Layer::export_region_slices_to_svg(const char *path) const
{
    BoundingBox bbox;
    for (const auto *region : m_regions)
        for (const auto &surface : region->slices())
            bbox.merge(get_extents(surface.expolygon));
    Point legend_size = export_surface_type_legend_to_svg_box_size();
    Point legend_pos(bbox.min(0), bbox.max(1));
    bbox.merge(Point(std::max(bbox.min(0) + legend_size(0), bbox.max(0)), bbox.max(1) + legend_size(1)));

    SVG svg(path, bbox);
    const float transparency = 0.5f;
    for (const auto *region : m_regions)
        for (const auto &surface : region->slices())
            svg.draw(surface.expolygon, surface_type_to_color_name(surface.surface_type), transparency);
    export_surface_type_legend_to_svg(svg, legend_pos);
    svg.Close(); 
}

// Export to "out/LayerRegion-name-%d.svg" with an increasing index with every export.
void Layer::export_region_slices_to_svg_debug(const char *name) const
{
    static size_t idx = 0;
    this->export_region_slices_to_svg(debug_out_path("Layer-slices-%s-%d.svg", name, idx ++).c_str());
}

void Layer::export_region_fill_surfaces_to_svg(const char *path) const
{
    BoundingBox bbox;
    for (const auto *region : m_regions)
        for (const auto &surface : region->slices())
            bbox.merge(get_extents(surface.expolygon));
    Point legend_size = export_surface_type_legend_to_svg_box_size();
    Point legend_pos(bbox.min(0), bbox.max(1));
    bbox.merge(Point(std::max(bbox.min(0) + legend_size(0), bbox.max(0)), bbox.max(1) + legend_size(1)));

    SVG svg(path, bbox);
    const float transparency = 0.5f;
    for (const auto *region : m_regions)
        for (const auto &surface : region->slices())
            svg.draw(surface.expolygon, surface_type_to_color_name(surface.surface_type), transparency);
    export_surface_type_legend_to_svg(svg, legend_pos);
    svg.Close();
}

// Export to "out/LayerRegion-name-%d.svg" with an increasing index with every export.
void Layer::export_region_fill_surfaces_to_svg_debug(const char *name) const
{
    static size_t idx = 0;
    this->export_region_fill_surfaces_to_svg(debug_out_path("Layer-fill_surfaces-%s-%d.svg", name, idx ++).c_str());
}

BoundingBox get_extents(const LayerRegion &layer_region)
{
    BoundingBox bbox;
    if (! layer_region.slices().empty()) {
        bbox = get_extents(layer_region.slices().surfaces.front());
        for (auto it = layer_region.slices().surfaces.cbegin() + 1; it != layer_region.slices().surfaces.cend(); ++ it)
            bbox.merge(get_extents(*it));
    }
    return bbox;
}

BoundingBox get_extents(const LayerRegionPtrs &layer_regions)
{
    BoundingBox bbox;
    if (!layer_regions.empty()) {
        bbox = get_extents(*layer_regions.front());
        for (auto it = layer_regions.begin() + 1; it != layer_regions.end(); ++it)
            bbox.merge(get_extents(**it));
    }
    return bbox;
}

}
