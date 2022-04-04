#include "SupportableIssuesSearch.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>

#include "libslic3r/Layer.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include "libslic3r/ClipperUtils.hpp"

namespace Slic3r {
namespace SupportableIssues {

struct Params {
    float bridge_distance = 5.0f;
    float printable_protrusion_distance = 1.0f;
};

namespace Impl {

struct LayerDescriptor {
    Vec2f centroid { 0.0f, 0.0f };
    size_t segments_count { 0 };
    float perimeter_length { 0.0f };
};

struct EdgeGridWrapper {
    EdgeGridWrapper(coord_t resolution, ExPolygons ex_polys) :
            ex_polys(ex_polys) {

        grid.create(this->ex_polys, resolution);
        grid.calculate_sdf();
    }
    EdgeGrid::Grid grid;
    ExPolygons ex_polys;
}
;

EdgeGridWrapper compute_layer_merged_edge_grid(const Layer *layer) {
    static const float eps = float(scale_(layer->object()->config().slice_closing_radius.value));
    // merge with offset
    ExPolygons merged = layer->merged(eps);
    // ofsset back
    ExPolygons layer_outline = offset_ex(merged, -eps);

    float min_region_flow_width { };
    for (const auto *region : layer->regions()) {
        min_region_flow_width = std::max(min_region_flow_width, region->flow(FlowRole::frExternalPerimeter).width());
    }
    std::cout << "min_region_flow_width: " << min_region_flow_width << std::endl;
    return EdgeGridWrapper(scale_(min_region_flow_width), layer_outline);
}

void check_extrusion_entity_stability(const ExtrusionEntity *entity, const EdgeGridWrapper &supported_grid,
        const Params &params) {
    if (entity->is_collection()){
        for (const auto* e: static_cast<ExtrusionEntityCollection>(entity).entities){
            check_extrusion_entity_stability(e, supported_grid, params);
        }
    } else { //single extrusion path, with possible varying parameters
        entity->as_polyline().points;
    }


}

void check_layer_stability(const PrintObject *po, size_t layer_idx, const Params &params) {
    if (layer_idx == 0) {
        // first layer is usually ok
        return;
    }
    const Layer *layer = po->get_layer(layer_idx);
    const Layer *prev_layer = layer->lower_layer;
    EdgeGridWrapper supported_grid = compute_layer_merged_edge_grid(prev_layer);

    for (const LayerRegion *layer_region : layer->regions()) {
        coordf_t flow_width = coordf_t(
                scale_(layer_region->flow(FlowRole::frExternalPerimeter).width()));
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                if (perimeter->role() == ExtrusionRole::erExternalPerimeter) {
                    check_extrusion_entity_stability(perimeter, supported_grid, params);
                } // ex_perimeter
            } // perimeter
        } // ex_entity
    }
}

} //Impl End

void quick_search(const PrintObject *po, const Params &params = Params { }) {
    using namespace Impl;
    std::vector<LayerDescriptor> descriptors(po->layer_count());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, po->layer_count()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    const Layer *layer = po->get_layer(layer_idx);

                    LayerDescriptor &descriptor = descriptors[layer_idx];
                    size_t point_count { 0 };

                    for (const LayerRegion *layer_region : layer->regions()) {
                        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
        if (perimeter->role() == ExtrusionRole::erExternalPerimeter) {
            assert(perimeter->is_loop());
            descriptor.segments_count++;
            const ExtrusionLoop *loop = static_cast<const ExtrusionLoop*>(perimeter);
            for (const ExtrusionPath& path : loop->paths) {
                Vec2f prev_pos = unscale(path.polyline.last_point()).cast<float>();
                for (size_t p_idx = 0; p_idx < path.polyline.points.size(); ++p_idx) {
                    point_count++;
                    Vec2f point_pos = unscale(path.polyline.points[p_idx]).cast<float>();
                    descriptor.centroid += point_pos;
                    descriptor.perimeter_length += (point_pos - prev_pos).norm();
                    prev_pos = point_pos;
                } //point
            } //path
        } // ex_perimeter
    } // perimeter
} // ex_entity
} // region

                    descriptor.centroid /= float(point_count);

                } // layer
            } // thread
            );

    for (size_t desc_idx = 0; desc_idx < descriptors.size(); ++desc_idx) {
        const LayerDescriptor &descriptor = descriptors[desc_idx];
        std::cout << "SIS layer idx: " << desc_idx << " reg count: " << descriptor.segments_count << " len: "
                << descriptor.perimeter_length <<
                " centroid: " << descriptor.centroid.x() << " | " << descriptor.centroid.y() << std::endl;
        if (desc_idx > 0) {
            const LayerDescriptor &prev = descriptors[desc_idx - 1];
            std::cout << "SIS diff: " << desc_idx << " reg count: "
                    << (int(descriptor.segments_count) - int(prev.segments_count)) <<
                    " len: " << (descriptor.perimeter_length - prev.perimeter_length) <<
                    " centroid: " << (descriptor.centroid - prev.centroid).norm() << std::endl;
        }
    }

}

}
}
