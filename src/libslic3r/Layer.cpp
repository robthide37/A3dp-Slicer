#include "Layer.hpp"
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
    const coord_t                                    offset_above,
    const coord_t                                    offset_end)
{
    class Visitor {
    public:
        Visitor(const std::vector<std::pair<coord_t, coord_t>> &intersections, 
            Layer &below, Layer &above, const coord_t offset_below, const coord_t offset_above, const coord_t offset_end) :
            m_intersections(intersections), m_below(below), m_above(above), m_offset_below(offset_below), m_offset_above(offset_above), m_offset_end(offset_end) {}

        void visit(const ClipperLib_Z::PolyNode &polynode)
        {
            if (polynode.Contour.size() >= 3) {
                int32_t i = 0, j = 0;
                double  area = 0;
                for (int icontour = 0; icontour <= polynode.ChildCount(); ++ icontour) {
                    const ClipperLib_Z::Path &contour = icontour == 0 ? polynode.Contour : polynode.Childs[icontour - 1]->Contour;
                    if (contour.size() >= 3) {
                        area = ClipperLib_Z::Area(contour);
                        int32_t i = contour.front().z();
                        int32_t j = i;
                        if (i < 0) {
                            std::tie(i, j) = m_intersections[-i - 1];
                        } else {
                            for (const ClipperLib_Z::IntPoint& pt : contour) {
                                j = pt.z();
                                if (j < 0) {
                                    std::tie(i, j) = m_intersections[-j - 1];
                                    goto end;
                                }
                                else if (i != j)
                                    goto end;
                            }
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
                                break;
                            }
                        }
                    } else {
                        // Index of an island above. Look-it up in the island below.
                        assert(j < m_offset_end);
                        j -= m_offset_below;
                        for (int l = int(m_below.lslices_ex.size()) - 1; l >= 0; -- l) {
                            LayerSlice &lslice = m_below.lslices_ex[l];
                            if (lslice.bbox.contains(pt) && m_below.lslices[l].contains(pt)) {
                                found = true;
                                i = l;
                                break;
                            }
                        }
                    }
                } else {
                    if (i > j)
                        std::swap(i, j);
                    assert(i >= m_offset_below);
                    assert(i < m_offset_above);
                    i -= m_offset_below;
                    assert(j >= m_offset_above);
                    assert(j < m_offset_end);
                    j -= m_offset_above;
                    found = true;
                }
                if (found) {
                    // Subtract area of holes from the area of outer contour.
                    for (int icontour = 0; icontour < polynode.ChildCount(); ++ icontour)
                        area -= ClipperLib_Z::Area(polynode.Childs[icontour]->Contour);
                    // Store the links and area into the contours.
                    LayerSlice::Links &links_below = m_below.lslices_ex[i].overlaps_above;
                    LayerSlice::Links &links_above = m_above.lslices_ex[i].overlaps_below;
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
        const coord_t                                   m_offset_end;
    } visitor(intersections, below, above, offset_below, offset_above, offset_end);

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
        LayerSlice::Links &links1 = above.lslices_ex[islice].overlaps_above;
        for (LayerSlice::Link &link1 : links1) {
            LayerSlice::Links &links2 = below.lslices_ex[link1.slice_idx].overlaps_below;
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
    coord_t             paths_end = paths_above_offset + coord_t(above.lslices.size());

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

    connect_layer_slices(below, above, result, zfill.intersections(), paths_below_offset, paths_above_offset, paths_end);
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
    std::vector<unsigned char> done(m_regions.size(), false);
    
    LayerRegionPtrs layerms;
    for (LayerRegionPtrs::iterator layerm = m_regions.begin(); layerm != m_regions.end(); ++ layerm) 
    	if ((*layerm)->slices().empty()) {
 			(*layerm)->m_perimeters.clear();
 			(*layerm)->m_fills.clear();
 			(*layerm)->m_thin_fills.clear();
    	} else {
	        size_t region_id = layerm - m_regions.begin();
	        if (done[region_id])
	            continue;
	        BOOST_LOG_TRIVIAL(trace) << "Generating perimeters for layer " << this->id() << ", region " << region_id;
	        done[region_id] = true;
	        const PrintRegionConfig &config = (*layerm)->region().config();
	        
	        // find compatible regions
            layerms.clear();
	        layerms.push_back(*layerm);
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
			 			other_layerm->m_perimeters.clear();
			 			other_layerm->m_fills.clear();
			 			other_layerm->m_thin_fills.clear();
		                layerms.push_back(other_layerm);
		                done[it - m_regions.begin()] = true;
		            }
		        }
	        
            SurfaceCollection fill_surfaces;
	        if (layerms.size() == 1) {  // optimization
                (*layerm)->m_fill_expolygons.clear();
	            (*layerm)->m_fill_surfaces.clear();
	            (*layerm)->make_perimeters((*layerm)->slices(), &fill_surfaces);
	            (*layerm)->m_fill_expolygons = to_expolygons(fill_surfaces.surfaces);
	        } else {
	            SurfaceCollection new_slices;
	            // Use the region with highest infill rate, as the make_perimeters() function below decides on the gap fill based on the infill existence.
	            LayerRegion *layerm_config = layerms.front();
	            {
	                // group slices (surfaces) according to number of extra perimeters
	                std::map<unsigned short, Surfaces> slices;  // extra_perimeters => [ surface, surface... ]
	                for (LayerRegion *layerm : layerms) {
	                    for (const Surface &surface : layerm->slices())
	                        slices[surface.extra_perimeters].emplace_back(surface);
	                    if (layerm->region().config().fill_density > layerm_config->region().config().fill_density)
	                    	layerm_config = layerm;
                        layerm->m_fill_surfaces.clear();
                        layerm->m_fill_expolygons.clear();
	                }
	                // merge the surfaces assigned to each group
	                for (std::pair<const unsigned short,Surfaces> &surfaces_with_extra_perimeters : slices)
	                    new_slices.append(offset_ex(surfaces_with_extra_perimeters.second, ClipperSafetyOffset), surfaces_with_extra_perimeters.second.front());
	            }
	            // make perimeters
	            layerm_config->make_perimeters(new_slices, &fill_surfaces);
	            // assign fill_surfaces to each layer
                if (! fill_surfaces.empty()) {
                    // Separate the fill surfaces.
                    for (LayerRegion *l : layerms)
                       l->m_fill_expolygons = intersection_ex(fill_surfaces.surfaces, l->slices().surfaces);
                }
	        }
	    }
    BOOST_LOG_TRIVIAL(trace) << "Generating perimeters for layer " << this->id() << " - Done";
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
