#include "SupportSpotsGenerator.hpp"

#include "ExtrusionEntity.hpp"
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

// #define DETAILED_DEBUG_LOGS
// #define DEBUG_FILES

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
    ExtrusionLine(const Vec2f &a, const Vec2f &b, const ExtrusionEntity *origin_entity) :
            a(a), b(b), len((a - b).norm()), origin_entity(origin_entity) {
    }

    float length() {
        return (a - b).norm();
    }

    bool is_external_perimeter() const {
        assert(origin_entity != nullptr);
        return origin_entity->role() == erExternalPerimeter || origin_entity->role() == erOverhangPerimeter;
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

SupportPoint::SupportPoint(const Vec3f &position, float force, float spot_radius, const Vec3f &direction) :
        position(position), force(force), spot_radius(spot_radius), direction(direction) {
}

class LinesDistancer {
private:
    std::vector<ExtrusionLine> lines;
    AABBTreeIndirect::Tree<2, float> tree;

public:
    explicit LinesDistancer(const std::vector<ExtrusionLine> &lines) :
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
    float area { };
    Vec3f centroid_accumulator = Vec3f::Zero();
    Vec2f second_moment_of_area_accumulator = Vec2f::Zero();
    float second_moment_of_area_covariance_accumulator { };

    void add(const IslandConnection &other) {
        this->area += other.area;
        this->centroid_accumulator += other.centroid_accumulator;
        this->second_moment_of_area_accumulator += other.second_moment_of_area_accumulator;
        this->second_moment_of_area_covariance_accumulator += other.second_moment_of_area_covariance_accumulator;
    }

    void print_info(const std::string &tag) {
        Vec3f centroid = centroid_accumulator / area;
        Vec2f variance =
                (second_moment_of_area_accumulator / area - centroid.head<2>().cwiseProduct(centroid.head<2>()));
        float covariance = second_moment_of_area_covariance_accumulator / area - centroid.x() * centroid.y();
        std::cout << tag << std::endl;
        std::cout << "area: " << area << std::endl;
        std::cout << "centroid: " << centroid.x() << " " << centroid.y() << " " << centroid.z() << std::endl;
        std::cout << "variance: " << variance.x() << " " << variance.y() << std::endl;
        std::cout << "covariance: " << covariance << std::endl;
    }
};

struct Island {
    std::unordered_map<size_t, IslandConnection> connected_islands { };
    float volume { };
    Vec3f volume_centroid_accumulator = Vec3f::Zero();
    float sticking_area { }; // for support points present on this layer (or bed extrusions)
    Vec3f sticking_centroid_accumulator = Vec3f::Zero();
    Vec2f sticking_second_moment_of_area_accumulator = Vec2f::Zero();
    float sticking_second_moment_of_area_covariance_accumulator { };

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

// Accumulator of current extrusion path properties
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

// base function: ((e^(((1)/(x^(2)+1)))-1)/(e-1))
// checkout e.g. here: https://www.geogebra.org/calculator
float gauss(float value, float mean_x_coord, float mean_value, float falloff_speed) {
    float shifted = value - mean_x_coord;
    float denominator = falloff_speed * shifted * shifted + 1.0f;
    float exponent = 1.0f / denominator;
    return mean_value * (std::exp(exponent) - 1.0f) / (std::exp(1.0f) - 1.0f);
}

void push_lines(const ExtrusionEntity *e, std::vector<ExtrusionLine>& destination)
{
    assert(!e->is_collection());
    Polyline pl = e->as_polyline();
    for (int point_idx = 0; point_idx < int(pl.points.size() - 1); ++point_idx) {
        Vec2f         start = unscaled(pl.points[point_idx]).cast<float>();
        Vec2f         next  = unscaled(pl.points[point_idx + 1]).cast<float>();
        ExtrusionLine line{start, next, e};
        destination.push_back(line);
    }
}

std::vector<ExtrusionLine> to_short_lines(const ExtrusionEntity *e, float length_limit)
{
    assert(!e->is_collection());
    Polyline pl = e->as_polyline();
    std::vector<ExtrusionLine> lines;
    lines.reserve(pl.points.size() * 1.5f);
    for (int point_idx = 0; point_idx < int(pl.points.size() - 1); ++point_idx) {
        Vec2f start        = unscaled(pl.points[point_idx]).cast<float>();
        Vec2f next         = unscaled(pl.points[point_idx + 1]).cast<float>();
        Vec2f v            = next - start; // vector from next to current
        float dist_to_next = v.norm();
        v.normalize();
        int   lines_count = int(std::ceil(dist_to_next / length_limit));
        float step_size   = dist_to_next / lines_count;
        for (int i = 0; i < lines_count; ++i) {
            Vec2f a(start + v * (i * step_size));
            Vec2f b(start + v * ((i + 1) * step_size));
            lines.emplace_back(a, b, e);
        }
    }
    return lines;
}

void check_extrusion_entity_stability(const ExtrusionEntity *entity,
        std::vector<ExtrusionLine> &checked_lines_out,
        float layer_z,
        const LayerRegion *layer_region,
        const LinesDistancer &prev_layer_lines,
        Issues &issues,
        const Params &params) {

    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            check_extrusion_entity_stability(e, checked_lines_out, layer_z, layer_region, prev_layer_lines,
                    issues, params);
        }
    } else { //single extrusion path, with possible varying parameters
        const auto to_vec3f = [layer_z](const Vec2f &point) {
            return Vec3f(point.x(), point.y(), layer_z);
        };
        float overhang_dist = tan(params.overhang_angle_deg * PI / 180.0f) * layer_region->layer()->height;
        float min_malformation_dist = tan(params.malformation_angle_span_deg.first * PI / 180.0f)
                * layer_region->layer()->height;
        float max_malformation_dist = tan(params.malformation_angle_span_deg.second * PI / 180.0f)
                * layer_region->layer()->height;

        std::vector<ExtrusionLine> lines = to_short_lines(entity, params.bridge_distance);
        if (lines.empty()) return;

        ExtrusionPropertiesAccumulator bridging_acc { };
        ExtrusionPropertiesAccumulator malformation_acc { };
        bridging_acc.add_distance(params.bridge_distance + 1.0f);
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
            // malformation in concave angles does not happen
            malformation_acc.add_angle(std::max(0.0f, curr_angle));

            size_t nearest_line_idx;
            Vec2f nearest_point;
            float dist_from_prev_layer = prev_layer_lines.signed_distance_from_lines(current_line.b, nearest_line_idx,
                    nearest_point);

            if (fabs(dist_from_prev_layer) < overhang_dist) {
                bridging_acc.reset();
            } else {
                bridging_acc.add_distance(current_line.len);
                // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                bool in_layer_dist_condition = bridging_acc.distance
                        > params.bridge_distance / (1.0f + (bridging_acc.max_curvature
                                * params.bridge_distance_decrease_by_curvature_factor / PI));
                bool between_layers_condition = fabs(dist_from_prev_layer) > flow_width ||
                        prev_layer_lines.get_line(nearest_line_idx).malformation > 3.0f * layer_region->layer()->height;

                if (in_layer_dist_condition && between_layers_condition) {
                    issues.support_points.emplace_back(to_vec3f(current_line.b), 0.0f, params.support_points_interface_radius, Vec3f(0.f, 0.0f, -1.0f));
                    current_line.support_point_generated = true;
                    bridging_acc.reset();
                }
            }

            //malformation
            if (fabs(dist_from_prev_layer) < 3.0f * flow_width) {
                const ExtrusionLine &nearest_line = prev_layer_lines.get_line(nearest_line_idx);
                current_line.malformation += 0.9 * nearest_line.malformation;
            }
            if (dist_from_prev_layer > min_malformation_dist && dist_from_prev_layer < max_malformation_dist) {
                malformation_acc.add_distance(current_line.len);
                current_line.malformation += layer_region->layer()->height *
                                             (0.5f + 1.5f * (malformation_acc.max_curvature / PI) *
                                                         gauss(malformation_acc.distance, 5.0f, 1.0f, 0.2f));
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

    //extract extrusions (connected paths from multiple lines) from the layer_lines. Grouping by the same polyline is determined by common origin_entity ptr.
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

    std::vector<LinesDistancer> islands; // these search trees will be used to determine to which island does the extrusion belong.
    std::vector<std::vector<size_t>> island_extrusions; //final assigment of each extrusion to an island.
    // initliaze the search from external perimeters - at the beginning, there is island candidate for each external perimeter.
    // some of them will disappear (e.g. holes)
    for (size_t e = 0; e < extrusions.size(); ++e) {
        if (layer_lines[extrusions[e].first].origin_entity->is_loop() &&
                layer_lines[extrusions[e].first].is_external_perimeter()) {
            std::vector<ExtrusionLine> copy(extrusions[e].second - extrusions[e].first);
            for (size_t ex_line_idx = extrusions[e].first; ex_line_idx < extrusions[e].second; ++ex_line_idx) {
                copy[ex_line_idx - extrusions[e].first] = layer_lines[ex_line_idx];
            }
            islands.emplace_back(copy);
            island_extrusions.push_back( { e });
        }
    }

    // backup code if islands not found
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
        if (!layer_lines[extrusions[e].first].origin_entity->is_loop() ||
                !layer_lines[extrusions[e].first].is_external_perimeter()) {
            bool island_assigned = false;
            for (size_t i = 0; i < islands.size(); ++i) {
                if (island_extrusions[i].empty()) {
                    continue;
                }
                size_t idx = 0;
                Vec2f pt = Vec2f::Zero();
                if (islands[i].signed_distance_from_lines(layer_lines[extrusions[e].first].a, idx, pt) < 0) {
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
                float volume = line.len * layer->height * flow_width * PI / 4.0f;
                island.volume += volume;
                island.volume_centroid_accumulator += to_3d(Vec2f((line.a + line.b) / 2.0f), float(layer->slice_z))
                        * volume;

                if (first_layer) {
                    float sticking_area = line.len * flow_width;
                    island.sticking_area += sticking_area;
                    Vec2f middle = Vec2f((line.a + line.b) / 2.0f);
                    island.sticking_centroid_accumulator += sticking_area * to_3d(middle, float(layer->slice_z));
                    // Bottom infill lines can be quite long, and algined, so the middle approximaton used above does not work
                    Vec2f dir = (line.b - line.a).normalized();
                    float segment_length = flow_width; // segments of size flow_width
                    for (float segment_middle_dist = std::min(line.len, segment_length * 0.5f);
                            segment_middle_dist < line.len;
                            segment_middle_dist += segment_length) {
                        Vec2f segment_middle = line.a + segment_middle_dist * dir;
                        island.sticking_second_moment_of_area_accumulator += segment_length * flow_width
                                * segment_middle.cwiseProduct(segment_middle);
                        island.sticking_second_moment_of_area_covariance_accumulator += segment_length * flow_width
                                * segment_middle.x()
                                * segment_middle.y();
                    }
                } else if (layer_lines[lidx].support_point_generated) {
                    float sticking_area = line.len * flow_width;
                    island.sticking_area += sticking_area;
                    island.sticking_centroid_accumulator += sticking_area * to_3d(line.b, float(layer->slice_z));
                    island.sticking_second_moment_of_area_accumulator += sticking_area * line.b.cwiseProduct(line.b);
                    island.sticking_second_moment_of_area_covariance_accumulator += sticking_area * line.b.x()
                            * line.b.y();
                }
            }
        }
        result.islands.push_back(island);
    }

    //LayerIslands structure built. Now determine connections and their areas to the previous layer using rasterization.
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
                IslandConnection &connection = result.islands[current_layer_grid.get_pixel(coords)]
                        .connected_islands[prev_layer_grid.get_pixel(coords)];
                Vec2f current_coords = current_layer_grid.get_pixel_center(coords);
                connection.area += current_layer_grid.pixel_area();
                connection.centroid_accumulator += to_3d(current_coords, result.layer_z)
                        * current_layer_grid.pixel_area();
                connection.second_moment_of_area_accumulator += current_coords.cwiseProduct(current_coords)
                        * current_layer_grid.pixel_area();
                connection.second_moment_of_area_covariance_accumulator += current_coords.x() * current_coords.y()
                        * current_layer_grid.pixel_area();
            }
        }
    }

    // filter out very small connection areas, they brake the graph building
    for (Island &island : result.islands) {
        std::vector<size_t> conns_to_remove;
        for (const auto &conn : island.connected_islands) {
            if (conn.second.area < params.connections_min_considerable_area) { conns_to_remove.push_back(conn.first); }
        }
        for (size_t conn : conns_to_remove) { island.connected_islands.erase(conn); }
    }

    return {result, current_layer_grid};
}

struct CoordinateFunctor {
    const std::vector<Vec3f> *coordinates;
    CoordinateFunctor(const std::vector<Vec3f> *coords) :
            coordinates(coords) {
    }
    CoordinateFunctor() :
            coordinates(nullptr) {
    }

    const float& operator()(size_t idx, size_t dim) const {
        return coordinates->operator [](idx)[dim];
    }
};

class ObjectPart {
    float volume { };
    Vec3f volume_centroid_accumulator = Vec3f::Zero();
    float sticking_area { };
    Vec3f sticking_centroid_accumulator = Vec3f::Zero();
    Vec2f sticking_second_moment_of_area_accumulator = Vec2f::Zero();
    float sticking_second_moment_of_area_covariance_accumulator { };

public:
    ObjectPart() = default;

    ObjectPart(const Island &island) {
        this->volume = island.volume;
        this->volume_centroid_accumulator = island.volume_centroid_accumulator;
        this->sticking_area = island.sticking_area;
        this->sticking_centroid_accumulator = island.sticking_centroid_accumulator;
        this->sticking_second_moment_of_area_accumulator = island.sticking_second_moment_of_area_accumulator;
        this->sticking_second_moment_of_area_covariance_accumulator =
                island.sticking_second_moment_of_area_covariance_accumulator;
    }

    float get_volume() const {
        return volume;
    }

    void add(const ObjectPart &other) {
        this->volume_centroid_accumulator += other.volume_centroid_accumulator;
        this->volume += other.volume;
        this->sticking_area += other.sticking_area;
        this->sticking_centroid_accumulator += other.sticking_centroid_accumulator;
        this->sticking_second_moment_of_area_accumulator += other.sticking_second_moment_of_area_accumulator;
        this->sticking_second_moment_of_area_covariance_accumulator +=
                other.sticking_second_moment_of_area_covariance_accumulator;
    }

    void add_support_point(const Vec3f &position, float sticking_area) {
        this->sticking_area += sticking_area;
        this->sticking_centroid_accumulator += sticking_area * position;
        this->sticking_second_moment_of_area_accumulator += sticking_area
                * position.head<2>().cwiseProduct(position.head<2>());
        this->sticking_second_moment_of_area_covariance_accumulator += sticking_area
                * position.x() * position.y();
    }

    float compute_directional_xy_variance(
            const Vec2f &line_dir,
            const Vec3f &centroid_accumulator,
            const Vec2f &second_moment_of_area_accumulator,
            const float &second_moment_of_area_covariance_accumulator,
            const float &area) const {
        assert(area > 0);
        Vec3f centroid = centroid_accumulator / area;
        Vec2f variance = (second_moment_of_area_accumulator / area
                - centroid.head<2>().cwiseProduct(centroid.head<2>()));
        float covariance = second_moment_of_area_covariance_accumulator / area - centroid.x() * centroid.y();
        // Var(aX+bY)=a^2*Var(X)+b^2*Var(Y)+2*a*b*Cov(X,Y)
        float directional_xy_variance = line_dir.x() * line_dir.x() * variance.x()
                + line_dir.y() * line_dir.y() * variance.y() +
                2.0f * line_dir.x() * line_dir.y() * covariance;
#ifdef DETAILED_DEBUG_LOGS
        BOOST_LOG_TRIVIAL(debug)
        << "centroid: " << centroid.x() << "  " << centroid.y() << "  " << centroid.z();
        BOOST_LOG_TRIVIAL(debug)
        << "variance: " << variance.x() << "  " << variance.y();
        BOOST_LOG_TRIVIAL(debug)
        << "covariance: " << covariance;
        BOOST_LOG_TRIVIAL(debug)
        << "directional_xy_variance: " << directional_xy_variance;
#endif
        return directional_xy_variance;
    }

    float compute_elastic_section_modulus(
            const Vec2f &line_dir,
            const Vec3f &extreme_point,
            const Vec3f &centroid_accumulator,
            const Vec2f &second_moment_of_area_accumulator,
            const float &second_moment_of_area_covariance_accumulator,
            const float &area) const {

        float directional_xy_variance = compute_directional_xy_variance(
                line_dir,
                centroid_accumulator,
                second_moment_of_area_accumulator,
                second_moment_of_area_covariance_accumulator,
                area);
        if (directional_xy_variance < EPSILON) {
            return 0.0f;
        }
        Vec3f centroid = centroid_accumulator / area;
        float extreme_fiber_dist = line_alg::distance_to(
                Linef(centroid.head<2>().cast<double>(),
                        (centroid.head<2>() + Vec2f(line_dir.y(), -line_dir.x())).cast<double>()),
                extreme_point.head<2>().cast<double>());
        float elastic_section_modulus = area * directional_xy_variance / extreme_fiber_dist;

#ifdef DETAILED_DEBUG_LOGS
        BOOST_LOG_TRIVIAL(debug)
        << "extreme_fiber_dist: " << extreme_fiber_dist;
        BOOST_LOG_TRIVIAL(debug)
        << "elastic_section_modulus: " << elastic_section_modulus;
#endif

        return elastic_section_modulus;
    }

    float is_stable_while_extruding(
            const IslandConnection &connection,
            const ExtrusionLine &extruded_line,
            const Vec3f &extreme_point,
            float layer_z,
            const Params &params) const {

        Vec2f line_dir = (extruded_line.b - extruded_line.a).normalized();
        const Vec3f &mass_centroid = this->volume_centroid_accumulator / this->volume;
        float mass = this->volume * params.filament_density;
        float weight = mass * params.gravity_constant;

        float movement_force = params.max_acceleration * mass;

        float extruder_conflict_force = params.standard_extruder_conflict_force +
                std::min(extruded_line.malformation, 1.0f) * params.malformations_additive_conflict_extruder_force;

        // section for bed calculations
        {
            if (this->sticking_area < EPSILON)
                return 1.0f;

            Vec3f bed_centroid = this->sticking_centroid_accumulator / this->sticking_area;
            float bed_yield_torque = -compute_elastic_section_modulus(
                    line_dir,
                    extreme_point,
                    this->sticking_centroid_accumulator,
                    this->sticking_second_moment_of_area_accumulator,
                    this->sticking_second_moment_of_area_covariance_accumulator,
                    this->sticking_area)
                    * params.get_bed_adhesion_yield_strength();

            Vec2f bed_weight_arm = (mass_centroid.head<2>() - bed_centroid.head<2>());
            float bed_weight_arm_len = bed_weight_arm.norm();
            float bed_weight_dir_xy_variance = compute_directional_xy_variance(bed_weight_arm,
                    this->sticking_centroid_accumulator,
                    this->sticking_second_moment_of_area_accumulator,
                    this->sticking_second_moment_of_area_covariance_accumulator,
                    this->sticking_area);
            float bed_weight_sign = bed_weight_arm_len < 2.0f * sqrt(bed_weight_dir_xy_variance) ? -1.0f : 1.0f;
            float bed_weight_torque = bed_weight_sign * bed_weight_arm_len * weight;

            float bed_movement_arm = std::max(0.0f, mass_centroid.z() - bed_centroid.z());
            float bed_movement_torque = movement_force * bed_movement_arm;

            float bed_conflict_torque_arm = layer_z - bed_centroid.z();
            float bed_extruder_conflict_torque = extruder_conflict_force * bed_conflict_torque_arm;

            float bed_total_torque = bed_movement_torque + bed_extruder_conflict_torque + bed_weight_torque
                    + bed_yield_torque;

#ifdef DETAILED_DEBUG_LOGS
            BOOST_LOG_TRIVIAL(debug)
            << "bed_centroid: " << bed_centroid.x() << "  " << bed_centroid.y() << "  " << bed_centroid.z();
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_yield_torque: " << bed_yield_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_weight_arm: " << bed_weight_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_weight_torque: " << bed_weight_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_movement_arm: " << bed_movement_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_movement_torque: " << bed_movement_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_conflict_torque_arm: " << bed_conflict_torque_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: extruded_line.malformation: " << extruded_line.malformation;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: extruder_conflict_force: " << extruder_conflict_force;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: bed_extruder_conflict_torque: " << bed_extruder_conflict_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: total_torque: " << bed_total_torque << "   layer_z: " << layer_z;
#endif

            if (bed_total_torque > 0)
                return bed_total_torque / bed_conflict_torque_arm;
        }

        //section for weak connection calculations
        {
            if (connection.area < EPSILON)
                return 1.0f;

            Vec3f conn_centroid = connection.centroid_accumulator / connection.area;

            if (layer_z - conn_centroid.z() < 3.0f) {
                return -1.0f;
            }
            float conn_yield_torque = compute_elastic_section_modulus(
                    line_dir,
                    extreme_point,
                    connection.centroid_accumulator,
                    connection.second_moment_of_area_accumulator,
                    connection.second_moment_of_area_covariance_accumulator,
                    connection.area) * params.material_yield_strength;

            float conn_weight_arm = (conn_centroid.head<2>() - mass_centroid.head<2>()).norm();
            float conn_weight_torque = conn_weight_arm * weight * (conn_centroid.z() / layer_z);

            float conn_movement_arm = std::max(0.0f, mass_centroid.z() - conn_centroid.z());
            float conn_movement_torque = movement_force * conn_movement_arm;

            float conn_conflict_torque_arm = layer_z - conn_centroid.z();
            float conn_extruder_conflict_torque = extruder_conflict_force * conn_conflict_torque_arm;

            float conn_total_torque = conn_movement_torque + conn_extruder_conflict_torque + conn_weight_torque
                    - conn_yield_torque;

#ifdef DETAILED_DEBUG_LOGS
            BOOST_LOG_TRIVIAL(debug)
            << "bed_centroid: " << conn_centroid.x() << "  " << conn_centroid.y() << "  " << conn_centroid.z();
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_yield_torque: " << conn_yield_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_weight_arm: " << conn_weight_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_weight_torque: " << conn_weight_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_movement_arm: " << conn_movement_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_movement_torque: " << conn_movement_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_conflict_torque_arm: " << conn_conflict_torque_arm;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: conn_extruder_conflict_torque: " << conn_extruder_conflict_torque;
            BOOST_LOG_TRIVIAL(debug)
            << "SSG: total_torque: " << conn_total_torque << "   layer_z: " << layer_z;
#endif

            return conn_total_torque / conn_conflict_torque_arm;
        }
    }
};

#ifdef DETAILED_DEBUG_LOGS
void debug_print_graph(const std::vector<LayerIslands> &islands_graph) {
    std::cout << "BUILT ISLANDS GRAPH:" << std::endl;
    for (size_t layer_idx = 0; layer_idx < islands_graph.size(); ++layer_idx) {
        std::cout << "ISLANDS AT LAYER: " << layer_idx << "  AT HEIGHT: " << islands_graph[layer_idx].layer_z
                << std::endl;
        for (size_t island_idx = 0; island_idx < islands_graph[layer_idx].islands.size(); ++island_idx) {
            const Island &island = islands_graph[layer_idx].islands[island_idx];
            std::cout << "        ISLAND " << island_idx << std::endl;
            std::cout << "              volume: " << island.volume << std::endl;
            std::cout << "              sticking_area: " << island.sticking_area << std::endl;
            std::cout << "              connected_islands count: " << island.connected_islands.size() << std::endl;
            std::cout << "              lines count: " << island.external_lines.size() << std::endl;
        }
    }
    std::cout << "END OF GRAPH" << std::endl;
}
#endif

class ActiveObjectParts {
    size_t next_part_idx = 0;
    std::unordered_map<size_t, ObjectPart> active_object_parts;
    std::unordered_map<size_t, size_t> active_object_parts_id_mapping;

public:
    size_t get_flat_id(size_t id) {
        size_t index = active_object_parts_id_mapping.at(id);
        while (index != active_object_parts_id_mapping.at(index)) {
            index = active_object_parts_id_mapping.at(index);
        }
        size_t i = id;
        while (index != active_object_parts_id_mapping.at(i)) {
            size_t next = active_object_parts_id_mapping[i];
            active_object_parts_id_mapping[i] = index;
            i = next;
        }
        return index;
    }

    ObjectPart& access(size_t id) {
        return this->active_object_parts.at(this->get_flat_id(id));
    }

    size_t insert(const Island &island) {
        this->active_object_parts.emplace(next_part_idx, ObjectPart(island));
        this->active_object_parts_id_mapping.emplace(next_part_idx, next_part_idx);
        return next_part_idx++;
    }

    void merge(size_t from, size_t to) {
        size_t to_flat = this->get_flat_id(to);
        size_t from_flat = this->get_flat_id(from);
        active_object_parts.at(to_flat).add(active_object_parts.at(from_flat));
        active_object_parts.erase(from_flat);
        active_object_parts_id_mapping[from] = to_flat;
    }
};

Issues check_global_stability(SupportGridFilter supports_presence_grid,
        const std::vector<LayerIslands> &islands_graph, const Params &params) {
#ifdef DETAILED_DEBUG_LOGS
    debug_print_graph(islands_graph);
#endif

    Issues issues { };
    ActiveObjectParts active_object_parts { };
    std::unordered_map<size_t, size_t> prev_island_to_object_part_mapping;
    std::unordered_map<size_t, size_t> next_island_to_object_part_mapping;

    std::unordered_map<size_t, IslandConnection> prev_island_weakest_connection;
    std::unordered_map<size_t, IslandConnection> next_island_weakest_connection;

    for (size_t layer_idx = 0; layer_idx < islands_graph.size(); ++layer_idx) {
        float layer_z = islands_graph[layer_idx].layer_z;

#ifdef DETAILED_DEBUG_LOGS
        for (const auto &m : prev_island_to_object_part_mapping) {
            std::cout << "island " << m.first << " maps to part " << m.second << std::endl;
            prev_island_weakest_connection[m.first].print_info("connection info:");
        }
#endif

        for (size_t island_idx = 0; island_idx < islands_graph[layer_idx].islands.size(); ++island_idx) {
            const Island &island = islands_graph[layer_idx].islands[island_idx];
            if (island.connected_islands.empty()) { //new object part emerging
                size_t part_id = active_object_parts.insert(island);
                next_island_to_object_part_mapping.emplace(island_idx, part_id);
                next_island_weakest_connection.emplace(island_idx,
                        IslandConnection { 1.0f, Vec3f::Zero(), Vec2f { INFINITY, INFINITY } });
            } else {
                size_t final_part_id { };
                IslandConnection transfered_weakest_connection { };
                IslandConnection new_weakest_connection { };
                // MERGE parts
                {
                    std::unordered_set<size_t> parts_ids;
                    for (const auto &connection : island.connected_islands) {
                        size_t part_id = active_object_parts.get_flat_id(
                                prev_island_to_object_part_mapping.at(connection.first));
                        parts_ids.insert(part_id);
                        transfered_weakest_connection.add(prev_island_weakest_connection.at(connection.first));
                        new_weakest_connection.add(connection.second);
                    }
                    final_part_id = *parts_ids.begin();
                    for (size_t part_id : parts_ids) {
                        if (final_part_id != part_id) {
                            active_object_parts.merge(part_id, final_part_id);
                        }
                    }
                }
                auto estimate_conn_strength = [layer_z](const IslandConnection &conn) {
                    Vec3f centroid = conn.centroid_accumulator / conn.area;
                    Vec2f variance = (conn.second_moment_of_area_accumulator / conn.area
                            - centroid.head<2>().cwiseProduct(centroid.head<2>()));
                    float xy_variance = variance.x() + variance.y();
                    float arm_len_estimate = std::max(1.0f, layer_z - (conn.centroid_accumulator.z() / conn.area));
                    return conn.area * sqrt(xy_variance) / arm_len_estimate;
                };

#ifdef DETAILED_DEBUG_LOGS
                new_weakest_connection.print_info("new_weakest_connection");
                transfered_weakest_connection.print_info("transfered_weakest_connection");
#endif

                if (estimate_conn_strength(transfered_weakest_connection)
                        > estimate_conn_strength(new_weakest_connection)) {
                    transfered_weakest_connection = new_weakest_connection;
                }
                next_island_weakest_connection.emplace(island_idx, transfered_weakest_connection);
                next_island_to_object_part_mapping.emplace(island_idx, final_part_id);
                ObjectPart &part = active_object_parts.access(final_part_id);
                part.add(ObjectPart(island));
            }
        }

        prev_island_to_object_part_mapping = next_island_to_object_part_mapping;
        next_island_to_object_part_mapping.clear();
        prev_island_weakest_connection = next_island_weakest_connection;
        next_island_weakest_connection.clear();

        // All object parts updated, inactive parts removed and weakest point of each island updated as well.
        // Now compute the stability of each active object part, adding supports where necessary, and also
        // check each island whether the weakest point is strong enough. If not, add supports as well.

        for (size_t island_idx = 0; island_idx < islands_graph[layer_idx].islands.size(); ++island_idx) {
            const Island &island = islands_graph[layer_idx].islands[island_idx];
            ObjectPart &part = active_object_parts.access(prev_island_to_object_part_mapping[island_idx]);

            IslandConnection &weakest_conn = prev_island_weakest_connection[island_idx];
#ifdef DETAILED_DEBUG_LOGS
            weakest_conn.print_info("weakest connection info: ");
#endif
            LinesDistancer island_lines_dist(island.external_lines);
            float unchecked_dist = params.min_distance_between_support_points + 1.0f;

            for (const ExtrusionLine &line : island.external_lines) {
                if ((unchecked_dist + line.len < params.min_distance_between_support_points
                        && line.malformation < 0.3f) || line.len == 0) {
                    unchecked_dist += line.len;
                } else {
                    unchecked_dist = line.len;
                    Vec2f target_point;
                    size_t idx;
                    Vec3f pivot_site_search_point = to_3d(Vec2f(line.b + (line.b - line.a).normalized() * 300.0f),
                            layer_z);
                    island_lines_dist.signed_distance_from_lines(pivot_site_search_point.head<2>(), idx,
                            target_point);
                    Vec3f support_point = to_3d(target_point, layer_z);
                    auto force = part.is_stable_while_extruding(weakest_conn, line, support_point, layer_z, params);
                    if (force > 0) {
                        if (!supports_presence_grid.position_taken(support_point)) {
                            float orig_area = params.support_points_interface_radius * params.support_points_interface_radius * float(PI);
                            // artifically lower the area for materials that have strong bed adhesion, as this adhesion does not apply for support points
                            float       altered_area = orig_area * params.get_support_spots_adhesion_strength() / params.get_bed_adhesion_yield_strength();
                            part.add_support_point(support_point, altered_area);

                            float radius = part.get_volume() < params.small_parts_threshold ? params.small_parts_support_points_interface_radius : params.support_points_interface_radius;
                            issues.support_points.emplace_back(support_point, force, radius, to_3d(Vec2f(line.b - line.a).normalized(), 0.0f));
                            supports_presence_grid.take_position(support_point);

                            weakest_conn.area += altered_area;
                            weakest_conn.centroid_accumulator += support_point * altered_area;
                            weakest_conn.second_moment_of_area_accumulator += altered_area *
                                                                              support_point.head<2>().cwiseProduct(support_point.head<2>());
                            weakest_conn.second_moment_of_area_covariance_accumulator += altered_area * support_point.x() *
                                                                                         support_point.y();
                        }
                    }
                }
            }
        }
        //end of iteration over layer
    }
    return issues;
}

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
    PixelGrid prev_layer_grid(po, flow_width*2.0f);

// PREPARE BASE LAYER
    const Layer *layer = po->layers()[0];
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                push_lines(perimeter, layer_lines);
            } // perimeter
        } // ex_entity
        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                push_lines(fill, layer_lines);
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
                        pos[1], layer->slice_z, color[0], color[1], color[2]);
            }
        }
    }
    for (const auto &line : layer_lines) {
        if (line.malformation > 0.0f) {
            Vec3f color = value_to_rgbf(-EPSILON, 1.0f, line.malformation);
            fprintf(malform_f, "v %f %f %f  %f %f %f\n", line.b[0],
                    line.b[1], layer->slice_z, color[0], color[1], color[2]);
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
                    check_extrusion_entity_stability(perimeter, layer_lines, layer->slice_z, layer_region,
                            external_lines, issues, params);
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill
                            || fill->role() == ExtrusionRole::erBridgeInfill) {
                        check_extrusion_entity_stability(fill, layer_lines, layer->slice_z, layer_region,
                                external_lines, issues, params);
                    } else {
                        push_lines(fill, layer_lines);
                    }
                } // fill
            } // ex_entity
        } // region

        auto [layer_islands, layer_grid] = reckon_islands(layer, false, 0, prev_layer_grid,
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
                            pos[1], layer->slice_z, color[0], color[1], color[2]);
                }
            }
        }
        for (const auto &line : layer_lines) {
            if (line.malformation > 0.0f) {
                Vec3f color = value_to_rgbf(0, 1.0f, line.malformation);
                fprintf(malform_f, "v %f %f %f  %f %f %f\n", line.b[0],
                        line.b[1], layer->slice_z, color[0], color[1], color[2]);
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

// std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
//     return {};
// }

Issues full_search(const PrintObject *po, const Params &params) {
    auto [local_issues, graph] = check_extrusions_and_build_graph(po, params);
    Issues global_issues = check_global_stability( { po, params.min_distance_between_support_points }, graph, params);
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

