#include "SupportSpotsGenerator.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <unordered_set>
#include <stack>

#include "AABBTreeLines.hpp"
#include "KDTreeIndirect.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "Geometry/ConvexHull.hpp"

#define DEBUG_FILES

#ifdef DEBUG_FILES
#include <boost/nowide/cstdio.hpp>
#include "libslic3r/Color.hpp"
#endif

namespace Slic3r {

class ExtrusionLine
{
public:
    ExtrusionLine() :
            a(Vec2f::Zero()), b(Vec2f::Zero()), len(0.0f), origin_entity(nullptr) {
    }
    ExtrusionLine(const Vec2f &_a, const Vec2f &_b, const ExtrusionEntity *origin_entity) :
            a(_a), b(_b), len((_a - _b).norm()), origin_entity(origin_entity) {
    }

    float length() {
        return (a - b).norm();
    }

    bool is_external_perimeter() const {
        assert(origin_entity != nullptr);
        return origin_entity->role() == erExternalPerimeter;
    }

    Vec2f a;
    Vec2f b;
    float len;
    const ExtrusionEntity *origin_entity;

    bool support_point_generated = false;
    float malformation = 0.0f;

    static const constexpr int Dim = 2;
    using Scalar = Vec2f::Scalar;
};

auto get_a(ExtrusionLine &&l) {
    return l.a;
}
auto get_b(ExtrusionLine &&l) {
    return l.b;
}

namespace SupportSpotsGenerator {

SupportPoint::SupportPoint(const Vec3f &position, float force, const Vec3f &direction) :
        position(position), force(force), direction(direction) {
}

class LinesDistancer {
private:
    std::vector<ExtrusionLine> lines;
    AABBTreeIndirect::Tree<2, float> tree;

public:
    explicit LinesDistancer(std::vector<ExtrusionLine> &lines) :
            lines(lines) {
        tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(this->lines);
    }

    // negative sign means inside
    float signed_distance_from_lines(const Vec2f &point, size_t &nearest_line_index_out,
            Vec2f &nearest_point_out) const {
        auto distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, point, nearest_line_index_out,
                nearest_point_out);
        if (distance < 0)
            return std::numeric_limits<float>::infinity();

        distance = sqrt(distance);
        const ExtrusionLine &line = lines[nearest_line_index_out];
        Vec2f v1 = line.b - line.a;
        Vec2f v2 = point - line.a;
        if ((v1.x() * v2.y()) - (v1.y() * v2.x()) > 0.0) {
            distance *= -1;
        }
        return distance;
    }

    const ExtrusionLine& get_line(size_t line_idx) const {
        return lines[line_idx];
    }

    const std::vector<ExtrusionLine>& get_lines() const {
        return lines;
    }
};

static const size_t NULL_ISLAND = std::numeric_limits<size_t>::max();

class PixelGrid {
    Vec2f pixel_size;
    Vec2f origin;
    Vec2f size;
    Vec2i pixel_count;

    std::vector<size_t> pixels { };

public:
    PixelGrid(const PrintObject *po, float resolution) {
        pixel_size = Vec2f(resolution, resolution);

        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec2f min = unscale(Vec2crd(-size_half.x(), -size_half.y())).cast<float>();
        Vec2f max = unscale(Vec2crd(size_half.x(), size_half.y())).cast<float>();

        origin = min;
        size = max - min;
        pixel_count = size.cwiseQuotient(pixel_size).cast<int>() + Vec2i::Ones();

        pixels.resize(pixel_count.y() * pixel_count.x());
        clear();
    }

    void distribute_edge(const Vec2f &p1, const Vec2f &p2, size_t value) {
        Vec2f dir = (p2 - p1);
        float length = dir.norm();
        if (length < 0.1) {
            return;
        }
        float step_size = this->pixel_size.x() / 2.0;

        float distributed_length = 0;
        while (distributed_length < length) {
            float next_len = std::min(length, distributed_length + step_size);
            Vec2f location = p1 + ((next_len / length) * dir);
            this->access_pixel(location) = value;

            distributed_length = next_len;
        }
    }

    void clear() {
        for (size_t &val : pixels) {
            val = NULL_ISLAND;
        }
    }

    float pixel_area() const {
        return this->pixel_size.x() * this->pixel_size.y();
    }

    size_t get_pixel(const Vec2i &coords) const {
        return pixels[this->to_pixel_index(coords)];
    }

    Vec2i get_pixel_count() {
        return pixel_count;
    }

    Vec2f get_pixel_center(const Vec2i &coords) const {
        return origin + coords.cast<float>().cwiseProduct(this->pixel_size)
                + this->pixel_size.cwiseQuotient(Vec2f(2.0f, 2.0f));
    }

private:
    Vec2i to_pixel_coords(const Vec2f &position) const {
        Vec2i pixel_coords = (position - this->origin).cwiseQuotient(this->pixel_size).cast<int>();
        return pixel_coords;
    }

    size_t to_pixel_index(const Vec2i &pixel_coords) const {
        assert(pixel_coords.x() >= 0);
        assert(pixel_coords.x() < pixel_count.x());
        assert(pixel_coords.y() >= 0);
        assert(pixel_coords.y() < pixel_count.y());

        return pixel_coords.y() * pixel_count.x() + pixel_coords.x();
    }

    size_t& access_pixel(const Vec2f &position) {
        return pixels[this->to_pixel_index(this->to_pixel_coords(position))];
    }
};

struct SupportGridFilter {
private:
    Vec3f cell_size;
    Vec3f origin;
    Vec3f size;
    Vec3i cell_count;

    std::unordered_set<size_t> taken_cells { };

public:
    SupportGridFilter(const PrintObject *po, float voxel_size) {
        cell_size = Vec3f(voxel_size, voxel_size, voxel_size);

        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec3f min = unscale(Vec3crd(-size_half.x(), -size_half.y(), 0)).cast<float>() - cell_size;
        Vec3f max = unscale(Vec3crd(size_half.x(), size_half.y(), po->height())).cast<float>() + cell_size;

        origin = min;
        size = max - min;
        cell_count = size.cwiseQuotient(cell_size).cast<int>() + Vec3i::Ones();
    }

    Vec3i to_cell_coords(const Vec3f &position) const {
        Vec3i cell_coords = (position - this->origin).cwiseQuotient(this->cell_size).cast<int>();
        return cell_coords;
    }

    size_t to_cell_index(const Vec3i &cell_coords) const {
        assert(cell_coords.x() >= 0);
        assert(cell_coords.x() < cell_count.x());
        assert(cell_coords.y() >= 0);
        assert(cell_coords.y() < cell_count.y());
        assert(cell_coords.z() >= 0);
        assert(cell_coords.z() < cell_count.z());

        return cell_coords.z() * cell_count.x() * cell_count.y()
                + cell_coords.y() * cell_count.x()
                + cell_coords.x();
    }

    Vec3f get_cell_center(const Vec3i &cell_coords) const {
        return origin + cell_coords.cast<float>().cwiseProduct(this->cell_size)
                + this->cell_size.cwiseQuotient(Vec3f(2.0f, 2.0f, 2.0));
    }

    void take_position(const Vec3f &position) {
        taken_cells.insert(to_cell_index(to_cell_coords(position)));
    }

    bool position_taken(const Vec3f &position) const {
        return taken_cells.find(to_cell_index(to_cell_coords(position))) != taken_cells.end();
    }

};

struct IslandConnection {
    float area;
    Vec3f centroid_accumulator;
};

struct Island {
    std::unordered_map<size_t, IslandConnection> connected_islands;
    std::vector<Vec3f> pivot_points; // for support points present on this layer (or bed extrusions)
    float volume;
    Vec3f volume_centroid_accumulator;
    float sticking_force; // for support points present on this layer (or bed extrusions)
    Vec3f sticking_centroid_accumulator;

    std::vector<ExtrusionLine> external_lines;
};

struct LayerIslands {
    std::vector<Island> islands;
    float layer_z;
};

float get_flow_width(const LayerRegion *region, ExtrusionRole role) {
    switch (role) {
        case ExtrusionRole::erBridgeInfill:
            return region->flow(FlowRole::frExternalPerimeter).width();
        case ExtrusionRole::erExternalPerimeter:
            return region->flow(FlowRole::frExternalPerimeter).width();
        case ExtrusionRole::erGapFill:
            return region->flow(FlowRole::frInfill).width();
        case ExtrusionRole::erPerimeter:
            return region->flow(FlowRole::frPerimeter).width();
        case ExtrusionRole::erSolidInfill:
            return region->flow(FlowRole::frSolidInfill).width();
        case ExtrusionRole::erInternalInfill:
            return region->flow(FlowRole::frInfill).width();
        case ExtrusionRole::erTopSolidInfill:
            return region->flow(FlowRole::frTopSolidInfill).width();
        default:
            return region->flow(FlowRole::frPerimeter).width();
    }
}

// Accumulator of current extruion path properties
// It remembers unsuported distance and maximum accumulated curvature over that distance.
// Used to determine local stability issues (too long bridges, extrusion curves into air)
struct ExtrusionPropertiesAccumulator {
    float distance = 0; //accumulated distance
    float curvature = 0; //accumulated signed ccw angles
    float max_curvature = 0; //max absolute accumulated value

    void add_distance(float dist) {
        distance += dist;
    }

    void add_angle(float ccw_angle) {
        curvature += ccw_angle;
        max_curvature = std::max(max_curvature, std::abs(curvature));
    }

    void reset() {
        distance = 0;
        curvature = 0;
        max_curvature = 0;
    }
};

void check_extrusion_entity_stability(const ExtrusionEntity *entity,
        std::vector<ExtrusionLine> &checked_lines_out,
        float print_z,
        const LayerRegion *layer_region,
        const LinesDistancer &prev_layer_lines,
        Issues &issues,
        const Params &params) {

    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            check_extrusion_entity_stability(e, checked_lines_out, print_z, layer_region, prev_layer_lines,
                    issues, params);
        }
    } else { //single extrusion path, with possible varying parameters
        const auto to_vec3f = [print_z](const Vec2f &point) {
            return Vec3f(point.x(), point.y(), print_z);
        };
        Points points { };
        entity->collect_points(points);
        std::vector<ExtrusionLine> lines;
        lines.reserve(points.size() * 1.5);
        lines.emplace_back(unscaled(points[0]).cast<float>(), unscaled(points[0]).cast<float>(), entity);
        for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
            Vec2f start = unscaled(points[point_idx]).cast<float>();
            Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
            Vec2f v = next - start; // vector from next to current
            float dist_to_next = v.norm();
            v.normalize();
            int lines_count = int(std::ceil(dist_to_next / params.bridge_distance));
            float step_size = dist_to_next / lines_count;
            for (int i = 0; i < lines_count; ++i) {
                Vec2f a(start + v * (i * step_size));
                Vec2f b(start + v * ((i + 1) * step_size));
                lines.emplace_back(a, b, entity);
            }
        }

        ExtrusionPropertiesAccumulator bridging_acc { };
        ExtrusionPropertiesAccumulator malformation_acc { };
        bridging_acc.add_distance(params.bridge_distance + 1.0f); // Initialise unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter starts and short loops into air.
        const float flow_width = get_flow_width(layer_region, entity->role());

        for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx) {
            ExtrusionLine &current_line = lines[line_idx];
            float curr_angle = 0;
            if (line_idx + 1 < lines.size()) {
                const Vec2f v1 = current_line.b - current_line.a;
                const Vec2f v2 = lines[line_idx + 1].b - lines[line_idx + 1].a;
                curr_angle = angle(v1, v2);
            }
            bridging_acc.add_angle(curr_angle);
            malformation_acc.add_angle(std::max(0.0f, curr_angle));

            size_t nearest_line_idx;
            Vec2f nearest_point;
            float dist_from_prev_layer = prev_layer_lines.signed_distance_from_lines(current_line.b, nearest_line_idx,
                    nearest_point);

            if (fabs(dist_from_prev_layer) < flow_width) {
                bridging_acc.reset();
            } else {
                bridging_acc.add_distance(current_line.len);
                if (bridging_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                > params.bridge_distance
                        / (1.0f + (bridging_acc.max_curvature
                                * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    issues.support_points.emplace_back(to_vec3f(current_line.b), 0.0f, Vec3f(0.f, 0.0f, -1.0f));
                    current_line.support_point_generated = true;
                    bridging_acc.reset();
                }
            }

            //malformation
            if (fabs(dist_from_prev_layer) < flow_width * 2.0f) {
                const ExtrusionLine &nearest_line = prev_layer_lines.get_line(nearest_line_idx);
                current_line.malformation += 0.9 * nearest_line.malformation;
            }
            if (dist_from_prev_layer > flow_width * 0.3) {
                malformation_acc.add_distance(current_line.len);
                current_line.malformation += 0.15f
                        * (0.8f + 0.2f * malformation_acc.max_curvature / (1.0f + 0.5f * malformation_acc.distance));
            } else {
                malformation_acc.reset();
            }
        }
        checked_lines_out.insert(checked_lines_out.end(), lines.begin(), lines.end());
    }
}

std::tuple<LayerIslands, PixelGrid> reckon_islands(
        const Layer *layer, bool first_layer,
        size_t prev_layer_islands_count,
        const PixelGrid &prev_layer_grid,
        const std::vector<ExtrusionLine> &layer_lines,
        const Params &params) {

    //extract extrusions (connected paths from multiple lines) from the layer_lines. belonging to single polyline is determined by origin_entity ptr.
    // result is a vector of [start, end) index pairs into the layer_lines vector
    std::vector<std::pair<size_t, size_t>> extrusions; //start and end idx (one beyond last extrusion) [start,end)
    const ExtrusionEntity *current_ex = nullptr;
    for (size_t lidx = 0; lidx < layer_lines.size(); ++lidx) {
        const ExtrusionLine &line = layer_lines[lidx];
        if (line.origin_entity == current_ex) {
            extrusions.back().second = lidx + 1;
        } else {
            extrusions.emplace_back(lidx, lidx + 1);
            current_ex = line.origin_entity;
        }
    }

    std::vector<LinesDistancer> islands; // these search trees will be used to determine to which island does the extrusion begin
    std::vector<std::vector<size_t>> island_extrusions; //final assigment of each extrusion to an island
    // initliaze the search from external perimeters - at the beginning, there is island candidate for each external perimeter.
    // some of them will disappear (e.g. holes)
    for (size_t e = 0; e < extrusions.size(); ++e) {
        if (layer_lines[extrusions[e].first].is_external_perimeter()) {
            std::vector<ExtrusionLine> copy(extrusions[e].second - extrusions[e].first);
            for (size_t ex_line_idx = extrusions[e].first; ex_line_idx < extrusions[e].second; ++ex_line_idx) {
                copy[ex_line_idx - extrusions[e].first] = layer_lines[ex_line_idx];
            }
            islands.emplace_back(copy);
            island_extrusions.push_back( { e });
        }
    }
    // backup code if islands not found - this can currently happen, as external perimeters may be also pure overhang perimeters, and there is no
    // way to distinguish external extrusions with total certainty.
    // If that happens, just make the first extrusion into island - it may be wrong, but it won't crash.
    if (islands.empty() && !extrusions.empty()) {
        std::vector<ExtrusionLine> copy(extrusions[0].second - extrusions[0].first);
        for (size_t ex_line_idx = extrusions[0].first; ex_line_idx < extrusions[0].second; ++ex_line_idx) {
            copy[ex_line_idx - extrusions[0].first] = layer_lines[ex_line_idx];
        }
        islands.emplace_back(copy);
        island_extrusions.push_back( { 0 });
    }

    // assign non external extrusions to islands
    for (size_t e = 0; e < extrusions.size(); ++e) {
        if (!layer_lines[extrusions[e].first].is_external_perimeter()) {
            bool island_assigned = false;
            for (size_t i = 0; i < islands.size(); ++i) {
                size_t _idx;
                Vec2f _pt;
                if (islands[i].signed_distance_from_lines(layer_lines[extrusions[e].first].a, _idx, _pt) < 0) {
                    island_extrusions[i].push_back(e);
                    island_assigned = true;
                    break;
                }
            }
            if (!island_assigned) { // If extrusion is not assigned for some reason, push it into the first island. As with the previous backup code,
                // it may be wrong, but it won't crash
                island_extrusions[0].push_back(e);
            }
        }
    }
    // merge islands which are embedded within each other (mainly holes)
    for (size_t i = 0; i < islands.size(); ++i) {
        if (islands[i].get_lines().empty()) {
            continue;
        }
        for (size_t j = 0; j < islands.size(); ++j) {
            if (islands[j].get_lines().empty() || i == j) {
                continue;
            }
            size_t _idx;
            Vec2f _pt;
            if (islands[i].signed_distance_from_lines(islands[j].get_line(0).a, _idx, _pt) < 0) {
                island_extrusions[i].insert(island_extrusions[i].end(), island_extrusions[j].begin(),
                        island_extrusions[j].end());
                island_extrusions[j].clear();
            }
        }
    }

    float flow_width = get_flow_width(layer->regions()[0], erExternalPerimeter);
    // after filtering the layer lines into islands, build the result LayerIslands structure.
    LayerIslands result { };
    result.layer_z = layer->slice_z;
    std::vector<size_t> line_to_island_mapping(layer_lines.size(), NULL_ISLAND);
    for (const std::vector<size_t> &island_ex : island_extrusions) {
        if (island_ex.empty()) {
            continue;
        }

        Island island { };
        island.external_lines.insert(island.external_lines.end(),
                layer_lines.begin() + extrusions[island_ex[0]].first,
                layer_lines.begin() + extrusions[island_ex[0]].second);
        for (size_t extrusion_idx : island_ex) {
            for (size_t lidx = extrusions[extrusion_idx].first; lidx < extrusions[extrusion_idx].second; ++lidx) {
                line_to_island_mapping[lidx] = result.islands.size();
                const ExtrusionLine &line = layer_lines[lidx];
                float volume = line.origin_entity->min_mm3_per_mm() * line.len;
                island.volume += volume;
                island.volume_centroid_accumulator += to_3d(Vec2f((line.a + line.b) / 2.0f), float(layer->print_z))
                        * volume;

                if (first_layer) {
                    float sticking_force = line.len * flow_width * params.base_adhesion;
                    island.sticking_force += sticking_force;
                    island.sticking_centroid_accumulator += sticking_force
                            * to_3d(Vec2f((line.a + line.b) / 2.0f), float(layer->print_z));
                    if (line.is_external_perimeter()) {
                        island.pivot_points.push_back(to_3d(Vec2f(line.b), float(layer->print_z)));
                    }
                } else if (layer_lines[lidx].support_point_generated) {
                    float support_interface_area = params.support_points_interface_radius
                            * params.support_points_interface_radius
                            * float(PI);
                    float sticking_force = support_interface_area * params.support_adhesion;
                    island.sticking_force += sticking_force;
                    island.sticking_centroid_accumulator += sticking_force
                            * to_3d(Vec2f(line.b), float(layer->print_z));
                    island.pivot_points.push_back(to_3d(Vec2f(line.b), float(layer->print_z)));
                }
            }
        }
        result.islands.push_back(island);
    }

    //LayerIslands structure built. Now determine connections and their areas to the previous layer using raterization.
    PixelGrid current_layer_grid = prev_layer_grid;
    current_layer_grid.clear();
    // build index image of current layer
    tbb::parallel_for(tbb::blocked_range<size_t>(0, layer_lines.size()),
            [&layer_lines, &current_layer_grid, &line_to_island_mapping](
                    tbb::blocked_range<size_t> r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    size_t island = line_to_island_mapping[i];
                    const ExtrusionLine &line = layer_lines[i];
                    current_layer_grid.distribute_edge(line.a, line.b, island);
                }
            });

    //compare the image of previous layer with the current layer. For each pair of overlapping valid pixels, add pixel area to the respective island connection
    for (size_t x = 0; x < size_t(current_layer_grid.get_pixel_count().x()); ++x) {
        for (size_t y = 0; y < size_t(current_layer_grid.get_pixel_count().y()); ++y) {
            Vec2i coords = Vec2i(x, y);
            if (current_layer_grid.get_pixel(coords) != NULL_ISLAND
                    && prev_layer_grid.get_pixel(coords) != NULL_ISLAND) {
                IslandConnection& connection = result.islands[current_layer_grid.get_pixel(coords)]
                                                              .connected_islands[prev_layer_grid.get_pixel(coords)];
                connection.area += current_layer_grid.pixel_area();
                connection.centroid_accumulator += to_3d(current_layer_grid.get_pixel_center(coords), result.layer_z) * current_layer_grid.pixel_area();
            }
        }
    }

    return {result, current_layer_grid};
}

struct ObjectPart {
    float volume { };
    Vec3f volume_centroid_accumulator = Vec3f::Zero();
    float sticking_force { };
    Vec3f sticking_centroid_accumulator = Vec3f::Zero();
    std::vector<Vec3f> pivot_points { };

    void add(const ObjectPart &other) {
        this->volume_centroid_accumulator += other.volume_centroid_accumulator;
        this->volume += other.volume;
        this->sticking_force += other.sticking_force;
        this->sticking_centroid_accumulator += other.sticking_centroid_accumulator;
        this->pivot_points.insert(this->pivot_points.end(), other.pivot_points.begin(), other.pivot_points.end());
    }

    ObjectPart(const Island &island) {
        this->volume = island.volume;
        this->volume_centroid_accumulator = island.volume_centroid_accumulator;
        this->sticking_force = island.sticking_force;
        this->sticking_centroid_accumulator = island.sticking_centroid_accumulator;
        this->pivot_points = island.pivot_points;
    }

    ObjectPart() = default;
};

struct WeakestConnection {
    float area = 0.0f;
    Vec3f centroid_accumulator = Vec3f::Zero();

    void add(const WeakestConnection& other) {
        this->area += other.area;
        this->centroid_accumulator += other.centroid_accumulator;
    }
};

Issues check_global_stability(const std::vector<LayerIslands> &islands_graph, const Params &params) {
    Issues issues { };
    size_t next_part_idx = 0;
    std::unordered_map<size_t, ObjectPart> active_object_parts;
    std::unordered_map<size_t, size_t> prev_island_to_object_part_mapping;
    std::unordered_map<size_t, size_t> next_island_to_object_part_mapping;

    std::unordered_map<size_t, WeakestConnection> prev_island_to_weakest_connection;
    std::unordered_map<size_t, WeakestConnection> next_island_to_weakest_connection;

    for (size_t layer_idx = 0; layer_idx < islands_graph.size(); ++layer_idx) {
        float layer_z = islands_graph[layer_idx].layer_z;
        std::unordered_set<size_t> layer_active_parts;
        std::cout << "at layer: " << layer_idx << "  the following island to object mapping is used:" << std::endl;
        for (const auto &m : prev_island_to_object_part_mapping) {
            std::cout << "island " << m.first << " maps to part " << m.second << std::endl;
            Vec3f connection_center = prev_island_to_weakest_connection[m.first].centroid_accumulator /  prev_island_to_weakest_connection[m.first].area;
            std::cout << " island has weak point with connection area: " <<
                    prev_island_to_weakest_connection[m.first].area << " and center: " <<
                    connection_center.x() << " " << connection_center.y() << " " << connection_center.z() << std::endl;
        }

        for (size_t island_idx = 0; island_idx < islands_graph[layer_idx].islands.size(); ++island_idx) {
            const Island &island = islands_graph[layer_idx].islands[island_idx];
            if (island.connected_islands.empty()) { //new object part emerging
                size_t part_idx = next_part_idx;
                next_part_idx++;
                active_object_parts.emplace(part_idx, ObjectPart(island));
                next_island_to_object_part_mapping.emplace(island_idx, part_idx);
                next_island_to_weakest_connection.emplace(island_idx,
                        WeakestConnection{INFINITY, Vec3f::Zero()});
                layer_active_parts.insert(part_idx);
            } else {
                size_t final_part_idx{};
                WeakestConnection transfered_weakest_connection{};
                WeakestConnection new_weakest_connection{};
                // MERGE parts
                {
                    std::unordered_set<size_t> part_indices;
                    for (const auto &connection : island.connected_islands) {
                        part_indices.insert(prev_island_to_object_part_mapping.at(connection.first));
                        transfered_weakest_connection.add(prev_island_to_weakest_connection.at(connection.first));
                        new_weakest_connection.area += connection.second.area;
                        new_weakest_connection.centroid_accumulator += connection.second.centroid_accumulator;
                    }
                    final_part_idx = *part_indices.begin();
                    for (size_t part_idx : part_indices) {
                        if (final_part_idx != part_idx) {
                            std::cout << "at layer: " << layer_idx << "  merging object part: " << part_idx
                                    << " into final part: " << final_part_idx << std::endl;
                            active_object_parts.at(final_part_idx).add(active_object_parts.at(part_idx));
                            active_object_parts.erase(part_idx);
                        }
                    }
                }
                auto estimate_strength = [layer_z](const WeakestConnection& conn){
                    float radius = fsqrt(conn.area / PI);
                    float arm_len_estimate = std::max(0.001f, layer_z - (conn.centroid_accumulator.z() / conn.area));
                    return radius * conn.area / arm_len_estimate;
                };

                if (estimate_strength(transfered_weakest_connection) < estimate_strength(new_weakest_connection)) {
                    new_weakest_connection = transfered_weakest_connection;
                }
                next_island_to_weakest_connection.emplace(island_idx, new_weakest_connection);
                next_island_to_object_part_mapping.emplace(island_idx, final_part_idx);
                ObjectPart &part = active_object_parts[final_part_idx];
                part.add(ObjectPart(island));
                layer_active_parts.insert(final_part_idx);
            }
        }

        std::unordered_set<size_t> parts_to_delete;
        for (const auto &part : active_object_parts) {
            if (layer_active_parts.find(part.first) == layer_active_parts.end()) {
                parts_to_delete.insert(part.first);
            } else {
                std::cout << "at layer " << layer_idx << " part is still active: " << part.first << std::endl;
            }
        }
        for (size_t part_id : parts_to_delete) {
            active_object_parts.erase(part_id);
            std::cout << " at layer: " << layer_idx << " removing object part " << part_id << std::endl;
        }
        prev_island_to_object_part_mapping = next_island_to_object_part_mapping;
        next_island_to_object_part_mapping.clear();
        prev_island_to_weakest_connection = next_island_to_weakest_connection;
        next_island_to_weakest_connection.clear();

    }
    return issues;
}

/*

 // islands_graph.back() refers to the top most (current) layer
 for (size_t island_idx = 0; island_idx < islands_graph.back().islands.size(); ++island_idx) {
 Island &island = islands_graph.back().islands[island_idx];

 std::vector<ExtrusionLine> island_external_lines;
 for (size_t lidx : islands_lines[island_idx]) {
 island_external_lines.push_back(layer_lines[lidx]);
 }
 LinesDistancer island_lines_dist(island_external_lines);
 Accumulator acc = island; // in acc, we accumulate the mass and other properties of the object part as we traverse the islands down to bed
 // There is one object part for each island at the top most layer, and each one is computed individually -
 // Some of the calculations will be done multiple times
 int layer_idx = islands_graph.size() - 1;
 // traverse the islands graph down, and for each connection area, calculate if it holds or breaks
 while (acc.connected_islands.size() > 0) {
 //test for break between layer_idx and layer_idx -1;
 LayerIslands below = islands_graph[layer_idx - 1]; // must exist, see while condition
 layer_idx--;
 // initialize variables that we will accumulate over all islands, which are connected to the current object part
 std::vector<Vec2f> pivot_points;
 Vec2f sticking_centroid;
 float connection_area = 0;
 for (const auto &pair : acc.connected_islands) {
 const Island &below_i = below.islands[pair.first];
 Vec2f centroid = (below_i.volume_centroid_accumulator / below_i.volume).head<2>(); // centroid of the island 'below_i'; TODO it should be centroid of the connection area
 pivot_points.push_back(centroid); // for object parts, we also consider breaking pivots in the centroids of the islands
 sticking_centroid += centroid * pair.second; // pair.second is connection area in mm^2
 connection_area += pair.second;
 }
 sticking_centroid /= connection_area; //normalize to get final sticking centroid
 for (const Vec3f &p_point : acc.pivot_points) {
 pivot_points.push_back(p_point.head<2>());
 }
 // Now we have accumulated pivot points, connection area and sticking centroid of the whole layer to the current object part

 // create KD tree over current pivot points
 auto coord_fn = [&pivot_points](size_t idx, size_t dim) {
 return pivot_points[idx][dim];
 };
 KDTreeIndirect<2, float, decltype(coord_fn)> pivot_points_tree(coord_fn, pivot_points.size());

 // iterate over extrusions at top layer island, check each for stability
 for (const ExtrusionLine &line : island_external_lines) {
 Vec2f line_dir = (line.b - line.a).normalized();
 Vec2f pivot_site_search_point = line.b + line_dir * 300.0f;
 size_t pivot_idx = find_closest_point(pivot_points_tree, pivot_site_search_point);
 const Vec2f &pivot = pivot_points[pivot_idx];

 float sticking_arm = (pivot - sticking_centroid).norm();
 float sticking_torque = sticking_arm * connection_area * params.tensile_strength; // For breakage in between layers, we compute with tensile strength, not bed adhesion

 float mass = acc.volume * params.filament_density;
 const Vec3f &mass_centorid = acc.volume_centroid_accumulator / acc.volume;
 float weight = mass * params.gravity_constant;
 float weight_arm = (pivot - mass_centorid.head<2>()).norm();
 float weight_torque = weight_arm * weight;

 float bed_movement_arm = mass_centorid.z();
 float bed_movement_force = params.max_acceleration * mass;
 float bed_movement_torque = bed_movement_force * bed_movement_arm;

 Vec3f extruder_pressure_direction = to_3d(line_dir, 0.0f);
 extruder_pressure_direction.z() = -0.2 - line.malformation * 0.5;
 extruder_pressure_direction.normalize();
 float conflict_torque_arm = (to_3d(Vec2f(pivot - line.b), print_z).cross(
 extruder_pressure_direction)).norm();
 float extruder_conflict_force = params.tolerable_extruder_conflict_force +
 std::min(line.malformation, 1.0f) * params.malformations_additive_conflict_extruder_force;
 float extruder_conflict_torque = extruder_conflict_force * conflict_torque_arm;

 float total_torque = bed_movement_torque + extruder_conflict_torque - weight_torque - sticking_torque;

 if (total_torque > 0) {
 Vec2f target_point { };
 size_t _idx { };
 island_lines_dist.signed_distance_from_lines(pivot_site_search_point, _idx, target_point);
 if (!supports_presence_grid.position_taken(to_3d(target_point, print_z))) {
 float area = params.support_points_interface_radius * params.support_points_interface_radius
 * float(PI);
 float sticking_force = area * params.support_adhesion;
 Vec3f support_point = to_3d(target_point, print_z);
 island.pivot_points.push_back(support_point);
 island.sticking_force += sticking_force;
 island.sticking_centroid_accumulator += sticking_force * support_point;
 issues.support_points.emplace_back(support_point,
 extruder_conflict_torque - sticking_torque, extruder_pressure_direction);
 supports_presence_grid.take_position(to_3d(target_point, print_z));
 }
 }
 #if 0
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: sticking_arm: " << sticking_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: sticking_torque: " << sticking_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: weight_arm: " << sticking_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: weight_torque: " << weight_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: bed_movement_arm: " << bed_movement_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: bed_movement_torque: " << bed_movement_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: conflict_torque_arm: " << conflict_torque_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: extruder_conflict_torque: " << extruder_conflict_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: total_torque: " << total_torque << "   printz: " << print_z;
 #endif
 }

 std::unordered_map<size_t, float> tmp = acc.connected_islands;
 acc.connected_islands.clear();
 // finally, add gathered islands to accumulator, and continue down to next layer
 for (const auto &pair : tmp) {
 const Island &below_i = below.islands[pair.first];
 for (const auto &below_islands : below_i.connected_islands) {
 acc.connected_islands[below_islands.first] += below_islands.second;
 }
 for (const Vec3f &pivot_p : below_i.pivot_points) {
 acc.pivot_points.push_back(pivot_p);
 }
 acc.sticking_centroid_accumulator += below_i.sticking_centroid_accumulator;
 acc.sticking_force += below_i.sticking_force;
 acc.volume += below_i.volume;
 acc.volume_centroid_accumulator += below_i.volume_centroid_accumulator;
 }
 }

 // We have arrived to the bed level, now check for stability of the object part on the bed
 std::vector<Vec2f> pivot_points;
 for (const Vec3f &p_point : acc.pivot_points) {
 pivot_points.push_back(p_point.head<2>());
 }
 auto coord_fn = [&pivot_points](size_t idx, size_t dim) {
 return pivot_points[idx][dim];
 };
 KDTreeIndirect<2, float, decltype(coord_fn)> pivot_points_tree(coord_fn, pivot_points.size());

 for (const ExtrusionLine &line : island_external_lines) {
 Vec2f line_dir = (line.b - line.a).normalized();
 Vec2f pivot_site_search_point = line.b + line_dir * 300.0f;
 size_t pivot_idx = find_closest_point(pivot_points_tree, pivot_site_search_point);
 const Vec2f &pivot = pivot_points[pivot_idx];

 const Vec2f &sticking_centroid = acc.sticking_centroid_accumulator.head<2>() / acc.sticking_force;
 float sticking_arm = (pivot - sticking_centroid).norm();
 float sticking_torque = sticking_arm * acc.sticking_force;

 float mass = acc.volume * params.filament_density;
 const Vec3f &mass_centorid = acc.volume_centroid_accumulator / acc.volume;
 float weight = mass * params.gravity_constant;
 float weight_arm = (pivot - mass_centorid.head<2>()).norm();
 float weight_torque = weight_arm * weight;

 float bed_movement_arm = mass_centorid.z();
 float bed_movement_force = params.max_acceleration * mass;
 float bed_movement_torque = bed_movement_force * bed_movement_arm;

 Vec3f extruder_pressure_direction = to_3d(line_dir, 0.0f);
 extruder_pressure_direction.z() = -0.2 - line.malformation * 0.5;
 extruder_pressure_direction.normalize();
 float conflict_torque_arm = (to_3d(Vec2f(pivot - line.b), print_z).cross(
 extruder_pressure_direction)).norm();
 float extruder_conflict_force = params.tolerable_extruder_conflict_force +
 std::min(line.malformation, 1.0f) * params.malformations_additive_conflict_extruder_force;
 float extruder_conflict_torque = extruder_conflict_force * conflict_torque_arm;

 float total_torque = bed_movement_torque + extruder_conflict_torque - weight_torque - sticking_torque;

 if (total_torque > 0) {
 Vec2f target_point;
 size_t _idx;
 island_lines_dist.signed_distance_from_lines(pivot_site_search_point, _idx, target_point);
 if (!supports_presence_grid.position_taken(to_3d(target_point, print_z))) {
 float area = params.support_points_interface_radius * params.support_points_interface_radius
 * float(PI);
 float sticking_force = area * params.support_adhesion;
 Vec3f support_point = to_3d(target_point, print_z);
 island.pivot_points.push_back(support_point);
 island.sticking_force += sticking_force;
 island.sticking_centroid_accumulator += sticking_force * support_point;
 issues.support_points.emplace_back(support_point,
 extruder_conflict_torque - sticking_torque, extruder_pressure_direction);
 supports_presence_grid.take_position(to_3d(target_point, print_z));
 }
 }
 #if 0
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: sticking_arm: " << sticking_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: sticking_torque: " << sticking_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: weight_arm: " << sticking_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: weight_torque: " << weight_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: bed_movement_arm: " << bed_movement_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: bed_movement_torque: " << bed_movement_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: conflict_torque_arm: " << conflict_torque_arm;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: extruder_conflict_torque: " << extruder_conflict_torque;
 BOOST_LOG_TRIVIAL(debug)
 << "SSG: total_torque: " << total_torque << "   printz: " << print_z;
 #endif
 }
 }

 return issues;
 */

std::tuple<Issues, std::vector<LayerIslands>> check_extrusions_and_build_graph(const PrintObject *po,
        const Params &params) {
#ifdef DEBUG_FILES
    FILE *segmentation_f = boost::nowide::fopen(debug_out_path("segmentation.obj").c_str(), "w");
    FILE *malform_f = boost::nowide::fopen(debug_out_path("malformations.obj").c_str(), "w");
#endif

    Issues issues { };
    std::vector<LayerIslands> islands_graph;
    std::vector<ExtrusionLine> layer_lines;
    float flow_width = get_flow_width(po->layers()[po->layer_count() - 1]->regions()[0], erExternalPerimeter);
    PixelGrid prev_layer_grid(po, flow_width);

    // PREPARE BASE LAYER
    const Layer *layer = po->layers()[0];
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                Points points { };
                perimeter->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next, perimeter };
                    layer_lines.push_back(line);
                }
                if (perimeter->is_loop()) {
                    Vec2f start = unscaled(points[points.size() - 1]).cast<float>();
                    Vec2f next = unscaled(points[0]).cast<float>();
                    ExtrusionLine line { start, next, perimeter };
                    layer_lines.push_back(line);
                }
            } // perimeter
        } // ex_entity
        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                Points points { };
                fill->collect_points(points);
                for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                    Vec2f start = unscaled(points[point_idx]).cast<float>();
                    Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                    ExtrusionLine line { start, next, fill };
                    layer_lines.push_back(line);
                }
            } // fill
        } // ex_entity
    } // region

    auto [layer_islands, layer_grid] = reckon_islands(layer, true, 0, prev_layer_grid,
            layer_lines, params);
    islands_graph.push_back(std::move(layer_islands));
#ifdef DEBUG_FILES
    for (size_t x = 0; x < size_t(layer_grid.get_pixel_count().x()); ++x) {
        for (size_t y = 0; y < size_t(layer_grid.get_pixel_count().y()); ++y) {
            Vec2i coords = Vec2i(x, y);
            size_t island_idx = layer_grid.get_pixel(coords);
            if (layer_grid.get_pixel(coords) != NULL_ISLAND) {
                Vec2f pos = layer_grid.get_pixel_center(coords);
                size_t pseudornd = ((island_idx + 127) * 33331 + 6907) % 23;
                Vec3f color = value_to_rgbf(0.0f, float(23), float(pseudornd));
                fprintf(segmentation_f, "v %f %f %f  %f %f %f\n", pos[0],
                        pos[1], layer->print_z, color[0], color[1], color[2]);
            }
        }
    }
    for (const auto &line : layer_lines) {
        if (line.malformation > 0.0f) {
            Vec3f color = value_to_rgbf(0, 1.0f, line.malformation);
            fprintf(malform_f, "v %f %f %f  %f %f %f\n", line.b[0],
                    line.b[1], layer->print_z, color[0], color[1], color[2]);
        }
    }
#endif
    LinesDistancer external_lines(layer_lines);
    layer_lines.clear();
    prev_layer_grid = layer_grid;

    for (size_t layer_idx = 1; layer_idx < po->layer_count(); ++layer_idx) {
        const Layer *layer = po->layers()[layer_idx];
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    check_extrusion_entity_stability(perimeter, layer_lines, layer->print_z, layer_region,
                            external_lines, issues, params);
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill
                            || fill->role() == ExtrusionRole::erBridgeInfill) {
                        check_extrusion_entity_stability(fill, layer_lines, layer->print_z, layer_region,
                                external_lines, issues, params);
                    } else {
                        Points points { };
                        fill->collect_points(points);
                        for (int point_idx = 0; point_idx < int(points.size() - 1); ++point_idx) {
                            Vec2f start = unscaled(points[point_idx]).cast<float>();
                            Vec2f next = unscaled(points[point_idx + 1]).cast<float>();
                            ExtrusionLine line { start, next, fill };
                            layer_lines.push_back(line);
                        }
                    }
                } // fill
            } // ex_entity
        } // region

        auto [layer_islands, layer_grid] = reckon_islands(layer, true, 0, prev_layer_grid,
                layer_lines, params);
        islands_graph.push_back(std::move(layer_islands));

#ifdef DEBUG_FILES
        for (size_t x = 0; x < size_t(layer_grid.get_pixel_count().x()); ++x) {
            for (size_t y = 0; y < size_t(layer_grid.get_pixel_count().y()); ++y) {
                Vec2i coords = Vec2i(x, y);
                size_t island_idx = layer_grid.get_pixel(coords);
                if (layer_grid.get_pixel(coords) != NULL_ISLAND) {
                    Vec2f pos = layer_grid.get_pixel_center(coords);
                    size_t pseudornd = ((island_idx + 127) * 33331 + 6907) % 23;
                    Vec3f color = value_to_rgbf(0.0f, float(23), float(pseudornd));
                    fprintf(segmentation_f, "v %f %f %f  %f %f %f\n", pos[0],
                            pos[1], layer->print_z, color[0], color[1], color[2]);
                }
            }
        }
        for (const auto &line : layer_lines) {
            if (line.malformation > 0.0f) {
                Vec3f color = value_to_rgbf(0, 1.0f, line.malformation);
                fprintf(malform_f, "v %f %f %f  %f %f %f\n", line.b[0],
                        line.b[1], layer->print_z, color[0], color[1], color[2]);
            }
        }
#endif
        external_lines = LinesDistancer(layer_lines);
        layer_lines.clear();
        prev_layer_grid = layer_grid;
    }

#ifdef DEBUG_FILES
    fclose(segmentation_f);
    fclose(malform_f);
#endif

    return {issues, islands_graph};
}

#ifdef DEBUG_FILES
void debug_export(Issues issues, std::string file_name) {
    Slic3r::CNumericLocalesSetter locales_setter;
    {
        FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_supports.obj").c_str()).c_str(), "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Debug files: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < issues.support_points.size(); ++i) {
            fprintf(fp, "v %f %f %f  %f %f %f\n", issues.support_points[i].position(0),
                    issues.support_points[i].position(1),
                    issues.support_points[i].position(2), 1.0, 0.0, 1.0);
        }

        fclose(fp);
    }
}
#endif

std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
    return {};
}

Issues full_search(const PrintObject *po, const Params &params) {
    auto [local_issues, graph] = check_extrusions_and_build_graph(po, params);
    Issues global_issues = check_global_stability(graph, params);
#ifdef DEBUG_FILES
    debug_export(local_issues, "local_issues");
    debug_export(global_issues, "global_issues");
#endif

    global_issues.support_points.insert(global_issues.support_points.end(),
            local_issues.support_points.begin(), local_issues.support_points.end());

    return global_issues;
}

} //SupportableIssues End
}

