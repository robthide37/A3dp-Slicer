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
            a(Vec2f::Zero()), b(Vec2f::Zero()), len(0.0f), external_perimeter(false) {
    }
    ExtrusionLine(const Vec2f &_a, const Vec2f &_b, bool external_perimeter) :
            a(_a), b(_b), len((_a - _b).norm()), external_perimeter(external_perimeter) {
    }

    float length() {
        return (a - b).norm();
    }

    Vec2f a;
    Vec2f b;
    float len;
    bool external_perimeter;

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

struct Island {
    std::unordered_map<size_t, float> islands_under_with_connection_area;
    std::vector<Vec3f> pivot_points;
    float volume;
    Vec3f volume_centroid;
    float sticking_force; // for support points present on this layer (or bed extrusions)
    Vec3f sticking_centroid;
};

struct LayerIslands {
    std::vector<Island> islands;
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
        bool is_ex_perimeter = entity->role() == erExternalPerimeter;
        lines.emplace_back(unscaled(points[0]).cast<float>(), unscaled(points[0]).cast<float>(), is_ex_perimeter);
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
                lines.emplace_back(a, b, is_ex_perimeter);
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

            if (dist_from_prev_layer < flow_width) {
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
                current_line.malformation += 0.15
                        * (0.8 + 0.2 * malformation_acc.max_curvature / (1.0f + 0.5f * malformation_acc.distance));
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

    BOOST_LOG_TRIVIAL(debug) << "SSG: reckon islands on printz: " << layer->print_z;

    std::vector<std::pair<size_t, size_t>> extrusions; //start and end idx (one beyond last extrusion) [start,end)
    Vec2f current_pt = layer_lines[0].a;
    std::pair<size_t, size_t> current_ext(0, 1);
    for (size_t lidx = 0; lidx < layer_lines.size(); ++lidx) {
        const ExtrusionLine &line = layer_lines[lidx];
        if (line.a == current_pt) {
            current_ext.second = lidx + 1;
        } else {
            extrusions.push_back(current_ext);
            current_ext.first = lidx;
            current_ext.second = lidx + 1;
        }
        current_pt = line.b;
    }

    BOOST_LOG_TRIVIAL(debug) << "SSG: layer_lines size: " << layer_lines.size();

    std::vector<LinesDistancer> islands;
    std::vector<std::vector<size_t>> island_extrusions;
    for (size_t e = 0; e < extrusions.size(); ++e) {
        if (layer_lines[extrusions[e].first].external_perimeter) {
            std::vector<ExtrusionLine> copy(extrusions[e].second - extrusions[e].first);
            for (size_t ex_line_idx = extrusions[e].first; ex_line_idx < extrusions[e].second; ++ex_line_idx) {
                copy[ex_line_idx - extrusions[e].first] = layer_lines[ex_line_idx];
            }
            islands.emplace_back(copy);
            island_extrusions.push_back( { e });
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "SSG: external perims: " << islands.size();

    for (size_t i = 0; i < islands.size(); ++i) {
        for (size_t e = 0; e < extrusions.size(); ++e) {
            if (!layer_lines[extrusions[e].first].external_perimeter) {
                size_t _idx;
                Vec2f _pt;
                if (islands[i].signed_distance_from_lines(layer_lines[extrusions[e].first].a, _idx, _pt) < 0) {
                    island_extrusions[i].push_back(e);
                }
            }
        }
    }

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

    BOOST_LOG_TRIVIAL(debug) << "SSG: filter islands";

    float flow_width = get_flow_width(layer->regions()[0], erExternalPerimeter);

    LayerIslands result { };
    std::vector<size_t> line_to_island_mapping(layer_lines.size(), NULL_ISLAND);
    for (const std::vector<size_t> &island_ex : island_extrusions) {
        if (island_ex.empty()) {
            continue;
        }

        Island island { };
        for (size_t extrusion_idx : island_ex) {
            for (size_t lidx = extrusions[extrusion_idx].first; lidx < extrusions[extrusion_idx].second; ++lidx) {
                line_to_island_mapping[lidx] = result.islands.size();
                const ExtrusionLine &line = layer_lines[lidx];
                float volume = line.len * flow_width * layer->height * 0.7; // 1/sqrt(2) compensation for cylindrical shape
                island.volume += volume;
                island.volume_centroid += to_3d(Vec2f((line.a + line.b) / 2.0f), float(layer->print_z)) * volume;

                if (first_layer) {
                    float sticking_force = line.len * flow_width * params.base_adhesion;
                    island.sticking_force += sticking_force;
                    island.sticking_centroid += sticking_force
                            * to_3d(Vec2f((line.a + line.b) / 2.0f), float(layer->print_z));
                    if (line.external_perimeter) {
                        island.pivot_points.push_back(to_3d(Vec2f(line.b), float(layer->print_z)));
                    }
                } else if (layer_lines[lidx].support_point_generated) {
                    float support_interface_area = params.support_points_interface_radius
                            * params.support_points_interface_radius
                            * float(PI);
                    float sticking_force = support_interface_area * params.support_adhesion;
                    island.sticking_force += sticking_force;
                    island.sticking_centroid += sticking_force * to_3d(Vec2f(line.b), float(layer->print_z));
                    island.pivot_points.push_back(to_3d(Vec2f(line.b), float(layer->print_z)));
                }
            }
        }
        result.islands.push_back(island);
    }

    BOOST_LOG_TRIVIAL(debug)
    << "SSG: There are " << result.islands.size() << " islands on printz: " << layer->print_z;

    PixelGrid current_layer_grid = prev_layer_grid;
    current_layer_grid.clear();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, layer_lines.size()),
            [&layer_lines, &current_layer_grid, &line_to_island_mapping](
                    tbb::blocked_range<size_t> r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    size_t island = line_to_island_mapping[i];
                    const ExtrusionLine &line = layer_lines[i];
                    current_layer_grid.distribute_edge(line.a, line.b, island);
                }
            });

    BOOST_LOG_TRIVIAL(debug) << "SSG: rasterized";

    for (size_t x = 0; x < size_t(current_layer_grid.get_pixel_count().x()); ++x) {
        for (size_t y = 0; y < size_t(current_layer_grid.get_pixel_count().y()); ++y) {
            Vec2i coords = Vec2i(x, y);
            if (current_layer_grid.get_pixel(coords) != NULL_ISLAND
                    && prev_layer_grid.get_pixel(coords) != NULL_ISLAND) {
                result.islands[current_layer_grid.get_pixel(coords)].islands_under_with_connection_area[prev_layer_grid.get_pixel(coords)] +=
                        current_layer_grid.pixel_area();
            }
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "SSG: connection area computed";

    return {result, current_layer_grid};
}

Issues check_object_stability(const PrintObject *po, const Params &params) {
    Issues issues { };
    std::vector<ExtrusionLine> layer_lines;
    float flow_width = get_flow_width(po->layers()[po->layer_count()-1]->regions()[0], erExternalPerimeter);
    PixelGrid prev_layer_grid(po, flow_width);
    BOOST_LOG_TRIVIAL(debug) << "SSG: flow width: " << flow_width;

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
                    ExtrusionLine line { start, next, perimeter->role() == erExternalPerimeter };
                    layer_lines.push_back(line);
                }
                if (perimeter->is_loop()) {
                    Vec2f start = unscaled(points[points.size() - 1]).cast<float>();
                    Vec2f next = unscaled(points[0]).cast<float>();
                    ExtrusionLine line { start, next, perimeter->role() == erExternalPerimeter };
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
                    ExtrusionLine line { start, next, false };
                    layer_lines.push_back(line);
                }
            } // fill
        } // ex_entity
    } // region

    auto [layer_islands, layer_grid] = reckon_islands(layer, true, 0, prev_layer_grid, layer_lines, params);
    std::remove_if(layer_lines.begin(), layer_lines.end(), [](const ExtrusionLine &line) {
        return !line.external_perimeter;
    });
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
                            ExtrusionLine line { start, next, false };
                            layer_lines.push_back(line);
                        }
                    }
                } // fill
            } // ex_entity
        } // region

        auto [layer_islands, layer_grid] = reckon_islands(layer, true, 0, prev_layer_grid, layer_lines, params);
        std::remove_if(layer_lines.begin(), layer_lines.end(), [](const ExtrusionLine &line) {
            return !line.external_perimeter;
        });
        layer_lines = std::vector<ExtrusionLine>();
        LinesDistancer external_lines(layer_lines);
        layer_lines.clear();
        prev_layer_grid = layer_grid;
    }

    return issues;
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
    check_object_stability(po, params);
    return {};
}

Issues
full_search(const PrintObject *po, const Params &params) {
    auto issues = check_object_stability(po, params);
#ifdef DEBUG_FILES
    debug_export(issues, "issues");
#endif
    return issues;

}
} //SupportableIssues End
}

