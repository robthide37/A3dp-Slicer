#include "SupportableIssuesSearch.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <unordered_set>
#include <stack>

#include "libslic3r/Layer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "Geometry/ConvexHull.hpp"
#include "PolygonPointTest.hpp"

#define DEBUG_FILES

#ifdef DEBUG_FILES
#include <boost/nowide/cstdio.hpp>
#include "libslic3r/Color.hpp"
#endif

namespace Slic3r {
namespace SupportableIssues {

void Issues::add(const Issues &layer_issues) {
    supports_nedded.insert(supports_nedded.end(), layer_issues.supports_nedded.begin(),
            layer_issues.supports_nedded.end());
    curling_up.insert(curling_up.end(), layer_issues.curling_up.begin(), layer_issues.curling_up.end());
}

bool Issues::empty() const {
    return supports_nedded.empty() && curling_up.empty();
}

SupportPoint::SupportPoint(const Vec3f &position, float weight) :
        position(position), weight(weight) {
}

CurledFilament::CurledFilament(const Vec3f &position, float estimated_height) :
        position(position), estimated_height(estimated_height) {
}

CurledFilament::CurledFilament(const Vec3f &position) :
        position(position), estimated_height(0.0f) {
}

struct Cell {
    float volume;
    float curled_height;
    int island_id = std::numeric_limits<int>::max();
};

struct CentroidAccumulator {
    Polygon convex_hull { };
    Points points { };
    Vec3f accumulated_value = Vec3f::Zero();
    float accumulated_volume { };
    float base_area { };
    float additional_supports_adhesion { };
    const float base_height { };

    explicit CentroidAccumulator(float base_height) :
            base_height(base_height) {
    }

    void calculate_base_hull() {
        convex_hull = Geometry::convex_hull(points);
        assert(convex_hull.is_counter_clockwise());
    }
};

struct CentroidAccumulators {
    std::unordered_map<int, size_t> mapping;
    std::vector<CentroidAccumulator> acccumulators;

    explicit CentroidAccumulators(size_t reserve_count) {
        acccumulators.reserve(reserve_count);
    }

    CentroidAccumulator& create_accumulator(int id, float base_height) {
        mapping[id] = acccumulators.size();
        acccumulators.push_back(CentroidAccumulator { base_height });
        return this->access(id);
    }

    CentroidAccumulator& access(int id) {
        return acccumulators[mapping[id]];
    }

    void merge_to(int from_id, int to_id) {
        CentroidAccumulator &from_acc = this->access(from_id);
        CentroidAccumulator &to_acc = this->access(to_id);
        if (&from_acc == &to_acc) {
            return;
        }
        to_acc.accumulated_value += from_acc.accumulated_value;
        to_acc.accumulated_volume += from_acc.accumulated_volume;
        to_acc.points.insert(to_acc.points.end(), from_acc.points.begin(), from_acc.points.end());
        to_acc.calculate_base_hull();
        to_acc.base_area += from_acc.base_area;
        mapping[from_id] = mapping[to_id];
    }
};

struct BalanceDistributionGrid {
    BalanceDistributionGrid() = default;

    void init(const PrintObject *po, size_t layer_idx_begin, size_t layer_idx_end) {
        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec3crd min = Vec3crd(-size_half.x(), -size_half.y(), 0);
        Vec3crd max = Vec3crd(size_half.x(), size_half.y(), po->height());

        cell_size = Vec3crd { int(cell_height * 2), int(cell_height * 2), int(cell_height) };
        assert(cell_size.x() == cell_size.y());

        global_origin = min;
        global_size = max - min;
        global_cell_count = global_size.cwiseQuotient(cell_size) + Vec3i::Ones();

        coord_t local_min_z = scale_(po->layers()[layer_idx_begin]->print_z);
        coord_t local_max_z = scale_(po->layers()[layer_idx_end > 0 ? layer_idx_end - 1 : 0]->print_z);
        int local_min_z_index = local_min_z / cell_size.z();
        int local_max_z_index = local_max_z / cell_size.z() + 1;

        local_z_index_offset = local_min_z_index;
        local_z_cell_count = local_max_z_index + 1 - local_min_z_index;

        cells.resize(local_z_cell_count * global_cell_count.y() * global_cell_count.x());
    }

    Vec3i to_global_cell_coords(const Vec3i &local_cell_coords) const {
        return local_cell_coords + local_z_index_offset * Vec3i::UnitZ();
    }

    Vec3i to_local_cell_coords(const Vec3i &global_cell_coords) const {
        return global_cell_coords - local_z_index_offset * Vec3i::UnitZ();
    }

    Vec3i to_global_cell_coords(const Point &p, float print_z) const {
        Vec3crd position = Vec3crd { p.x(), p.y(), int(scale_(print_z)) };
        Vec3i cell_coords = (position - this->global_origin).cwiseQuotient(this->cell_size);
        return cell_coords;
    }

    Vec3i to_global_cell_coords(const Vec3f &position) const {
        Vec3crd scaled_position = scaled(position);
        Vec3i cell_coords = (scaled_position - this->global_origin).cwiseQuotient(this->cell_size);
        return cell_coords;
    }

    Vec3i to_local_cell_coords(const Point &p, float print_z) const {
        Vec3i cell_coords = this->to_global_cell_coords(p, print_z);
        return this->to_local_cell_coords(cell_coords);
    }

    size_t to_cell_index(const Vec3i &local_cell_coords) const {
        assert(local_cell_coords.x() >= 0);
        assert(local_cell_coords.x() < global_cell_count.x());
        assert(local_cell_coords.y() >= 0);
        assert(local_cell_coords.y() < global_cell_count.y());
        assert(local_cell_coords.z() >= 0);
        assert(local_cell_coords.z() < local_z_cell_count);

        return local_cell_coords.z() * global_cell_count.x() * global_cell_count.y()
                + local_cell_coords.y() * global_cell_count.x()
                + local_cell_coords.x();
    }

    Vec3crd get_cell_center(const Vec3i &global_cell_coords) const {
        return global_origin + global_cell_coords.cwiseProduct(this->cell_size)
                + this->cell_size.cwiseQuotient(Vec3crd(2, 2, 2));
    }

    Cell& access_cell(const Point &p, float print_z) {
        return cells[this->to_cell_index(to_local_cell_coords(p, print_z))];
    }

    Cell& access_cell(const Vec3f &unscaled_position) {
        return cells[this->to_cell_index(this->to_local_cell_coords(this->to_global_cell_coords(unscaled_position)))];
    }

    Cell& access_cell(const Vec3i &local_cell_coords) {
        return cells[this->to_cell_index(local_cell_coords)];
    }

    const Cell& access_cell(const Vec3i &local_cell_coords) const {
        return cells[this->to_cell_index(local_cell_coords)];
    }

    void distribute_edge(const Point &p1, const Point &p2, float print_z, float unscaled_width, float unscaled_height) {
        Vec2d dir = (p2 - p1).cast<double>();
        double length = dir.norm();
        if (length < 0.1) {
            return;
        }
        double step_size = this->cell_size.x() / 2.0;

        float diameter = unscaled_height * unscaled_width * 0.7071f; // constant to simulate somewhat elliptical shape (1/sqrt(2))

        double distributed_length = 0;
        while (distributed_length < length) {
            double next_len = std::min(length, distributed_length + step_size);
            double current_dist_payload = next_len - distributed_length;

            Point location = p1 + ((next_len / length) * dir).cast<coord_t>();
            float payload = unscale<float>(current_dist_payload) * diameter;
            this->access_cell(location, print_z).volume += payload;

            distributed_length = next_len;
        }
    }

    void merge(const BalanceDistributionGrid &other) {
        int z_start = std::max(local_z_index_offset, other.local_z_index_offset);
        int z_end = std::min(local_z_index_offset + local_z_cell_count,
                other.local_z_index_offset + other.local_z_cell_count);

        for (int x = 0; x < global_cell_count.x(); ++x) {
            for (int y = 0; y < global_cell_count.y(); ++y) {
                for (int z = z_start; z < z_end; ++z) {
                    Vec3i global_coords { x, y, z };
                    Vec3i local_coords = this->to_local_cell_coords(global_coords);
                    Vec3i other_local_coords = other.to_local_cell_coords(global_coords);
                    this->access_cell(local_coords).volume += other.access_cell(other_local_coords).volume;
                }
            }
        }
    }

    void analyze(Issues &issues, const Params &params) {
        const auto validate_xy_coords = [&](const Vec2i &local_coords) {
            return local_coords.x() >= 0 && local_coords.y() >= 0 && local_coords.x() < this->global_cell_count.x()
                    && local_coords.y() < this->global_cell_count.y();
        };
        CentroidAccumulators accumulators(issues.supports_nedded.size() + 4);
        auto custom_comparator = [](const Vec2i &left, const Vec2i &right) {
            if (left.x() == right.x()) {
                return left.y() < right.y();
            }
            return left.x() < right.x();
        };

        int next_island_id = -1;
        std::set<Vec2i, decltype(custom_comparator)> coords_to_check(custom_comparator);
        for (int y = 0; y < global_cell_count.y(); ++y) {
            for (int x = 0; x < global_cell_count.x(); ++x) {
                Cell &origin_cell = this->access_cell(Vec3i(x, y, 0));
                if (origin_cell.volume > 0 && origin_cell.island_id == std::numeric_limits<int>::max()) {
                    CentroidAccumulator &acc = accumulators.create_accumulator(next_island_id, 0);
                    coords_to_check.clear();
                    coords_to_check.insert(Vec2i(x, y));
                    while (!coords_to_check.empty()) {
                        Vec2i current_coords = *coords_to_check.begin();
                        coords_to_check.erase(coords_to_check.begin());
                        if (!validate_xy_coords(current_coords)) {
                            continue;
                        }
                        Cell &cell = this->access_cell(Vec3i(current_coords.x(), current_coords.y(), 0));
                        if (cell.volume <= 0 || cell.island_id != std::numeric_limits<int>::max()) {
                            continue;
                        }
                        cell.island_id = next_island_id;
                        Vec3crd cell_center = this->get_cell_center(
                                Vec3i(current_coords.x(), current_coords.y(), local_z_index_offset));
                        acc.points.push_back(Point(cell_center.head<2>()));
                        acc.accumulated_value += unscale(cell_center).cast<float>() * cell.volume;
                        acc.accumulated_volume += cell.volume;

                        for (int y_offset = -1; y_offset <= 1; ++y_offset) {
                            for (int x_offset = -1; x_offset <= 1; ++x_offset) {
                                if (y_offset != 0 || x_offset != 0) {
                                    coords_to_check.insert(
                                            Vec2i(current_coords.x() + x_offset, current_coords.y() + y_offset));
                                }
                            }
                        }
                    }
                    next_island_id--;
                    acc.calculate_base_hull();
                    acc.base_area = unscale<float>(unscale<float>(acc.convex_hull.area())); //apply unscale 2x, it has units of area
                }
            }
        }

        std::sort(issues.supports_nedded.begin(), issues.supports_nedded.end(),
                [](const SupportPoint &left, const SupportPoint &right) {
                    return left.position.z() < right.position.z();
                });
        for (int index = 0; index < int(issues.supports_nedded.size()); ++index) {
            Vec3i local_coords = this->to_local_cell_coords(
                    this->to_global_cell_coords(issues.supports_nedded[index].position));
            this->access_cell(local_coords).island_id = index;
            CentroidAccumulator &acc = accumulators.create_accumulator(index,
                    issues.supports_nedded[index].position.z());
            acc.points.push_back(Point(scaled(Vec2f(issues.supports_nedded[index].position.head<2>()))));
            acc.calculate_base_hull();
            acc.base_area = params.support_points_interface_area;
        }

        for (const CurledFilament &curling : issues.curling_up) {
            this->access_cell(curling.position).curled_height += curling.estimated_height;
        }

        std::unordered_set<int> modified_acc_ids;
        modified_acc_ids.reserve(issues.supports_nedded.size() + 1);
        for (int z = 1; z < local_z_cell_count; ++z) {
            std::cout << "current z: " << z << std::endl;

            modified_acc_ids.clear();

            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    Cell &current = this->access_cell(Vec3i(x, y, z));

                    // distribute curling
                    if (current.volume > 0) {
                        float curled_height = 0;
                        for (int y_offset = -2; y_offset <= 2; ++y_offset) {
                            for (int x_offset = -2; x_offset <= 2; ++x_offset) {
                                if (validate_xy_coords(Vec2i(x + x_offset, y + y_offset))) {
                                    Cell &under = this->access_cell(Vec3i(x + x_offset, y + y_offset, z - 1));
                                    curled_height = std::max(curled_height, under.curled_height);
                                }
                            }
                        }
                        bool curled = current.curled_height > 0;
                        current.curled_height += std::max(0.0f, float(curled_height - unscaled(this->cell_size.z())));
                        if (!curled) {
                            current.curled_height /= 4.0f;
                        }
                    }

                    // distribute islands info
                    if (current.volume > 0 && current.island_id == std::numeric_limits<int>::max()) {
                        int min_island_id_found = std::numeric_limits<int>::max();
                        std::unordered_set<int> ids_to_merge { };
                        for (int y_offset = -2; y_offset <= 2; ++y_offset) {
                            for (int x_offset = -2; x_offset <= 2; ++x_offset) {
                                if (validate_xy_coords(Vec2i(x + x_offset, y + y_offset))) {
                                    Cell &under = this->access_cell(Vec3i(x + x_offset, y + y_offset, z - 1));
                                    if (under.island_id < min_island_id_found) {
                                        min_island_id_found = under.island_id;
                                    }
                                    ids_to_merge.insert(under.island_id);
                                }
                            }
                        }
                        // assign island and update its info
                        if (min_island_id_found < std::numeric_limits<int>::max()) {
                            ids_to_merge.erase(std::numeric_limits<int>::max());
                            ids_to_merge.erase(min_island_id_found);
                            current.island_id = min_island_id_found;
                            for (auto id : ids_to_merge) {
                                accumulators.merge_to(id, min_island_id_found);
                            }

                            CentroidAccumulator &acc = accumulators.access(min_island_id_found);
                            acc.accumulated_value += current.volume
                                    * unscale(this->get_cell_center(this->to_global_cell_coords(Vec3i(x, y, z)))).cast<
                                            float>();
                            acc.accumulated_volume += current.volume;
                            modified_acc_ids.insert(min_island_id_found);
                        }
                    }
                }
            }

            std::cout << " check all active accumulators " << std::endl;

            for (int acc_id : modified_acc_ids) {

                std::cout << "Z:  " << z << "   controlling acc id: " << acc_id << std::endl;

                CentroidAccumulator &acc = accumulators.access(acc_id);
                Vec3f centroid = acc.accumulated_value / acc.accumulated_volume;

                std::cout << "acc.accumulated_value :   " << acc.accumulated_value.x() << "  "
                        << acc.accumulated_value.y() << " " << acc.accumulated_value.z() << std::endl;
                std::cout << "acc.accumulated_volume :   " << acc.accumulated_volume << std::endl;
                std::cout << "centroid:   " << centroid.x() << "  " << centroid.y() << " " << centroid.z() << std::endl;

                //determine signed shortest distance to the convex hull
                Point centroid_base_projection = Point(scaled(Vec2f(centroid.head<2>())));
                Point pivot;
                double distance_scaled_sq = std::numeric_limits<double>::max();
                bool inside = true;
                if (acc.convex_hull.points.size() == 1) {
                    pivot = acc.convex_hull.points[0];
                    distance_scaled_sq = (pivot - centroid_base_projection).squaredNorm();
                    inside = true;
                } else {
                    for (Line line : acc.convex_hull.lines()) {
                        Point closest_point;
                        double dist_sq = line.distance_to_squared(centroid_base_projection, &closest_point);
                        if (dist_sq < distance_scaled_sq) {
                            pivot = closest_point;
                            distance_scaled_sq = dist_sq;
                        }
                        if ((centroid_base_projection - closest_point).cast<double>().dot(line.normal().cast<double>())
                                > 0) {
                            inside = false;
                        }
                    }
                }

                Vec3f pivot_3d;
                pivot_3d << unscale(pivot).cast<float>(), acc.base_height;
                float embedded_distance = unscaled(sqrt(distance_scaled_sq));
                float centroid_pivot_distance = (centroid - pivot_3d).norm();
                float base_center_pivot_distance = float(unscale(Vec2crd(acc.convex_hull.centroid() - pivot)).norm());

                std::cout << "centroid inside ? " << inside << "  and embedded distance is: " << embedded_distance
                        << std::endl;

                bool additional_supports_needed = false;
                float sticking_force = acc.base_area
                        * (acc.base_height == 0 ? params.base_adhesion : params.support_adhesion);
                float sticking_torque = base_center_pivot_distance * sticking_force;

                std::cout << "sticking force: " << sticking_force << "   sticking torque: " << sticking_torque
                        << std::endl;

                float xy_movement_force = acc.accumulated_volume * params.filament_density * params.max_acceleration;
                float xy_movement_torque = xy_movement_force * centroid_pivot_distance;

                std::cout << "xy_movement_force: " << xy_movement_force << "  xy_movement_torque: "
                        << xy_movement_torque << std::endl;

                float weight = acc.accumulated_volume * params.filament_density * params.gravity_constant;
                float weight_torque = embedded_distance * weight;
                if (!inside) {
                    weight_torque *= -1;
                }
                std::cout << "weight: " << weight << "  weight_torque: " << weight_torque << std::endl;

                float extruder_conflict_torque = params.tolerable_extruder_conflict_force * 2.0f
                        * centroid_pivot_distance;
                std::cout << "extruder_conflict_torque: " << extruder_conflict_torque << std::endl;

                float total_momentum = sticking_torque + weight_torque - xy_movement_torque - extruder_conflict_torque;
                additional_supports_needed = total_momentum < 0;

                std::cout << "total_momentum: " << total_momentum << std::endl;
                std::cout << "additional supports needed: " << additional_supports_needed << std::endl;

                if (additional_supports_needed) {
                    Vec2f attractor_dir =
                            unscale(Vec2crd(inside ?
                                                     pivot - centroid_base_projection :
                                                     centroid_base_projection - pivot)).cast<float>().normalized();
                    Vec2f attractor = unscale(centroid_base_projection).cast<float>() + 10000 * attractor_dir;

                    std::cout << " attractor:  " << attractor.x() << " | " << attractor.y() << std::endl;

                    double min_dist = std::numeric_limits<double>::max();
                    Vec3f support_point = centroid;
                    Vec2i coords = Vec2i(0, 0);
                    for (int y = 0; y < global_cell_count.y(); ++y) {
                        for (int x = 0; x < global_cell_count.x(); ++x) {
                            Cell &cell = this->access_cell(Vec3i(x, y, z));
                            if (cell.island_id != std::numeric_limits<int>::max() &&
                                    &accumulators.access(cell.island_id) == &acc) {
                                Vec3f cell_center =
                                        unscale(this->get_cell_center(this->to_global_cell_coords(Vec3i(x, y, z)))).cast<
                                                float>();
                                float dist_sq = (cell_center.head<2>() - attractor).squaredNorm();
                                if (dist_sq < min_dist) {
                                    min_dist = dist_sq;
                                    support_point = cell_center;
                                    coords = Vec2i(x, y);
                                }
                            }
                        }
                    }

                    int final_height_coords = z;
                    while (final_height_coords > 0
                            && this->access_cell(Vec3i(coords.x(), coords.y(), final_height_coords)).volume > 0) {
                        final_height_coords--;
                    }
                    support_point.z() = unscaled(
                            (final_height_coords + this->local_z_index_offset) * this->cell_size.z());
                    float expected_force = total_momentum / (support_point - pivot_3d).norm();

                    std::cout << " new support point:  " << support_point.x() << " | " << support_point.y() << " | "
                            << support_point.z() << std::endl;
                    std::cout << " expected_force:  " << expected_force << std::endl;

                    issues.supports_nedded.emplace_back(support_point, expected_force);
                    acc.points.push_back(Point::new_scale(Vec2f(support_point.head<2>())));
                    acc.base_area += params.support_points_interface_area;
                    acc.calculate_base_hull();
                }

            }
        }
    }

#ifdef DEBUG_FILES
    void debug_export() const {
        Slic3r::CNumericLocalesSetter locales_setter;
        {
            FILE *volume_grid_file = boost::nowide::fopen(debug_out_path("volume_grid.obj").c_str(), "w");
            FILE *islands_grid_file = boost::nowide::fopen(debug_out_path("islands_grid.obj").c_str(), "w");
            FILE *curling_grid_file = boost::nowide::fopen(debug_out_path("curling_grid.obj").c_str(), "w");

            if (volume_grid_file == nullptr || islands_grid_file == nullptr || curling_grid_file == nullptr) {
                BOOST_LOG_TRIVIAL(error)
                << "Debug files: Couldn't open debug file for writing, destination: " << debug_out_path("");
                return;
            }

            float max_volume = 0;
            int min_island_id = 0;
            int max_island_id = 0;
            float max_curling_height = 0;

            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    for (int z = 0; z < local_z_cell_count; ++z) {
                        const Cell &cell = access_cell(Vec3i(x, y, z));
                        max_volume = std::max(max_volume, cell.volume);
                        if (cell.island_id != std::numeric_limits<int>::max()) {
                            min_island_id = std::min(min_island_id, cell.island_id);
                            max_island_id = std::max(max_island_id, cell.island_id);
                        }
                        max_curling_height = std::max(max_curling_height, cell.curled_height);
                    }
                }
            }

            for (int x = 0; x < global_cell_count.x(); ++x) {
                for (int y = 0; y < global_cell_count.y(); ++y) {
                    for (int z = 0; z < local_z_cell_count; ++z) {
                        Vec3f center = unscale(get_cell_center(to_global_cell_coords(Vec3i { x, y, z }))).cast<float>();
                        const Cell &cell = access_cell(Vec3i(x, y, z));
                        if (cell.volume != 0) {
                            auto volume_color = value_to_rgbf(0, cell.volume, cell.volume);
                            fprintf(volume_grid_file, "v %f %f %f  %f %f %f\n", center(0), center(1), center(2),
                                    volume_color.x(), volume_color.y(), volume_color.z());
                        }
                        if (cell.island_id != std::numeric_limits<int>::max()) {
                            auto island_color = value_to_rgbf(min_island_id, max_island_id + 1, cell.island_id);
                            fprintf(islands_grid_file, "v %f %f %f  %f %f %f\n", center(0), center(1), center(2),
                                    island_color.x(), island_color.y(), island_color.z());
                        }
                        if (cell.curled_height > 0) {
                            auto curling_color = value_to_rgbf(0, max_curling_height, cell.curled_height);
                            fprintf(curling_grid_file, "v %f %f %f  %f %f %f\n", center(0), center(1), center(2),
                                    curling_color.x(), curling_color.y(), curling_color.z());
                        }
                    }
                }
            }

            fclose(volume_grid_file);
            fclose(islands_grid_file);
            fclose(curling_grid_file);
        }
    }
#endif

    static constexpr float cell_height = scale_(0.3f);

    Vec3crd cell_size { };

    Vec3crd global_origin { };
    Vec3crd global_size { };
    Vec3i global_cell_count { };

    int local_z_index_offset { };
    int local_z_cell_count { };
    std::vector<Cell> cells { };

}
;

namespace Impl {

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

        for (size_t i = 0; i < issues.supports_nedded.size(); ++i) {
            fprintf(fp, "v %f %f %f  %f %f %f\n", issues.supports_nedded[i].position(0),
                    issues.supports_nedded[i].position(1),
                    issues.supports_nedded[i].position(2), 1.0, 0.0, 1.0);
        }

        fclose(fp);
    }
    {
        FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_curling.obj").c_str()).c_str(), "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error)
            << "Debug files: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < issues.curling_up.size(); ++i) {
            fprintf(fp, "v %f %f %f  %f %f %f\n", issues.curling_up[i].position(0),
                    issues.curling_up[i].position(1),
                    issues.curling_up[i].position(2), 0.0, 1.0, 0.0);
        }
        fclose(fp);
    }

}
#endif

EdgeGridWrapper compute_layer_edge_grid(const Layer *layer) {
    float min_region_flow_width { 1.0f };
    for (const auto *region : layer->regions()) {
        min_region_flow_width = std::min(min_region_flow_width,
                region->flow(FlowRole::frExternalPerimeter).width());
    }
    std::vector<Points> lines;
    for (const LayerRegion *layer_region : layer->regions()) {
        for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
            lines.push_back(Points { });
            ex_entity->collect_points(lines.back());
        } // ex_entity

        for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
            lines.push_back(Points { });
            ex_entity->collect_points(lines.back());
        } // ex_entity
    }

    return EdgeGridWrapper(scale_(min_region_flow_width), lines);
}

//TODO needs revision
coordf_t get_flow_width(const LayerRegion *region, ExtrusionRole role) {
    switch (role) {
        case ExtrusionRole::erBridgeInfill:
            return region->flow(FlowRole::frExternalPerimeter).scaled_width();
        case ExtrusionRole::erExternalPerimeter:
            return region->flow(FlowRole::frExternalPerimeter).scaled_width();
        case ExtrusionRole::erGapFill:
            return region->flow(FlowRole::frInfill).scaled_width();
        case ExtrusionRole::erPerimeter:
            return region->flow(FlowRole::frPerimeter).scaled_width();
        case ExtrusionRole::erSolidInfill:
            return region->flow(FlowRole::frSolidInfill).scaled_width();
        case ExtrusionRole::erInternalInfill:
            return region->flow(FlowRole::frInfill).scaled_width();
        case ExtrusionRole::erTopSolidInfill:
            return region->flow(FlowRole::frTopSolidInfill).scaled_width();
        default:
            return region->flow(FlowRole::frPerimeter).scaled_width();
    }
}

coordf_t get_max_allowed_distance(ExtrusionRole role, coordf_t flow_width, bool external_perimeters_first,
        const Params &params) { // <= distance / flow_width (can be larger for perimeter, if not external perimeter first)
    if ((role == ExtrusionRole::erExternalPerimeter || role == ExtrusionRole::erOverhangPerimeter)
            && (external_perimeters_first)) {
        return params.max_first_ex_perim_unsupported_distance_factor * flow_width;
    } else {
        return params.max_unsupported_distance_factor * flow_width;
    }
}

struct SegmentAccumulator {
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

Issues check_extrusion_entity_stability(const ExtrusionEntity *entity, float print_z, const LayerRegion *layer_region,
        const EdgeGridWrapper &supported_grid, const Params &params) {

    Issues issues { };
    if (entity->is_collection()) {
        for (const auto *e : static_cast<const ExtrusionEntityCollection*>(entity)->entities) {
            issues.add(
                    check_extrusion_entity_stability(e, print_z, layer_region, supported_grid, params));
        }
    } else { //single extrusion path, with possible varying parameters
        //prepare stack of points on the extrusion path. If there are long segments, additional points might be pushed onto the stack during the algorithm.
        std::stack<Point> points { };
        for (const auto &p : entity->as_polyline().points) {
            points.push(p);
        }

        SegmentAccumulator supports_acc { };
        supports_acc.add_distance(params.bridge_distance + 1.0f); // initialize unsupported distance with larger than tolerable distance ->
        // -> it prevents extruding perimeter start and short loops into air.

        const auto to_vec3f = [print_z](const Point &point) {
            Vec2f tmp = unscale(point).cast<float>();
            return Vec3f(tmp.x(), tmp.y(), print_z);
        };
        float region_height = layer_region->layer()->height;

        Point prev_point = points.top(); // prev point of the path. Initialize with first point.
        Vec3f prev_fpoint = to_vec3f(prev_point);
        coordf_t flow_width = get_flow_width(layer_region, entity->role());
        bool external_perimters_first = layer_region->region().config().external_perimeters_first;
        const coordf_t max_allowed_dist_from_prev_layer = get_max_allowed_distance(entity->role(), flow_width,
                external_perimters_first, params);

        while (!points.empty()) {
            Point point = points.top();
            points.pop();
            Vec3f fpoint = to_vec3f(point);
            float edge_len = (fpoint - prev_fpoint).norm();

            coordf_t dist_from_prev_layer { 0 };
            if (!supported_grid.signed_distance(point, flow_width, dist_from_prev_layer)) { // dist from prev layer not found, assume empty layer
                issues.supports_nedded.push_back(SupportPoint(fpoint, 1.0f));
                supports_acc.reset();
            }

            float angle = 0;
            if (!points.empty()) {
                const Vec2f v1 = (fpoint - prev_fpoint).head<2>();
                const Vec2f v2 = unscale(points.top()).cast<float>() - fpoint.head<2>();
                float dot = v1(0) * v2(0) + v1(1) * v2(1);
                float cross = v1(0) * v2(1) - v1(1) * v2(0);
                angle = float(atan2(float(cross), float(dot))); // ccw angle, TODO replace with angle func, once it gets into master
            }

            supports_acc.add_angle(angle);

            if (dist_from_prev_layer > max_allowed_dist_from_prev_layer) { //extrusion point is unsupported
                supports_acc.add_distance(edge_len); // for algorithm simplicity, expect that the whole line between prev and current point is unsupported

                if (supports_acc.distance // if unsupported distance is larger than bridge distance linearly decreased by curvature, enforce supports.
                > params.bridge_distance
                        / (1.0f + (supports_acc.max_curvature
                                * params.bridge_distance_decrease_by_curvature_factor / PI))) {
                    issues.supports_nedded.push_back(SupportPoint(fpoint, 1.0f));
                    supports_acc.reset();
                }
            } else {
                supports_acc.reset();
            }

            // Estimation of short curvy segments which are not supported -> problems with curling
            if (dist_from_prev_layer > -max_allowed_dist_from_prev_layer * 0.7071) { //extrusion point is unsupported or poorly supported
                issues.curling_up.push_back(
                        CurledFilament(fpoint, 2.0f * region_height + region_height * 6.0f * std::abs(angle) / PI));
            }

            prev_point = point;
            prev_fpoint = fpoint;

            if (!points.empty()) { //oversampling if necessary
                Vec2f next = unscale(points.top()).cast<float>();
                Vec2f reverse_v = fpoint.head<2>() - next; // vector from next to current
                float dist_to_next = reverse_v.norm();
                reverse_v.normalize();
                int new_points_count = dist_to_next / params.bridge_distance;
                float step_size = dist_to_next / (new_points_count + 1);
                for (int i = 1; i <= new_points_count; ++i) {
                    points.push(Point::new_scale(Vec2f(next + reverse_v * (i * step_size))));
                }
            }

        }
    }
    return issues;
}

void distribute_layer_volume(const PrintObject *po, size_t layer_idx,
        BalanceDistributionGrid &balance_grid) {
    const Layer *layer = po->get_layer(layer_idx);
    for (const LayerRegion *region : layer->regions()) {
        for (const ExtrusionEntity *collections : region->fills.entities) {
            for (const ExtrusionEntity *entity : static_cast<const ExtrusionEntityCollection*>(collections)->entities) {
                for (const Line &line : entity->as_polyline().lines()) {
                    balance_grid.distribute_edge(line.a, line.b, layer->print_z,
                            unscale<float>(get_flow_width(region, entity->role())), layer->height);
                }
            }
        }
        for (const ExtrusionEntity *collections : region->perimeters.entities) {
            for (const ExtrusionEntity *entity : static_cast<const ExtrusionEntityCollection*>(collections)->entities) {
                for (const Line &line : entity->as_polyline().lines()) {
                    balance_grid.distribute_edge(line.a, line.b, layer->print_z,
                            unscale<float>(get_flow_width(region, entity->role())), layer->height);
                }
            }
        }
    }
}

Issues check_layer_stability(const PrintObject *po, size_t layer_idx, bool full_check,
        const Params &params) {
    const Layer *layer = po->get_layer(layer_idx);
    //Prepare edge grid of previous layer, will be used to check if the extrusion path is supported
    EdgeGridWrapper supported_grid = compute_layer_edge_grid(layer->lower_layer);

    Issues issues { };
    if (full_check) { // If full check; check stability of perimeters, gap fills, and bridges.
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    issues.add(
                            check_extrusion_entity_stability(perimeter, layer->print_z, layer_region,
                                    supported_grid, params));
                } // perimeter
            } // ex_entity
            for (const ExtrusionEntity *ex_entity : layer_region->fills.entities) {
                for (const ExtrusionEntity *fill : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (fill->role() == ExtrusionRole::erGapFill
                            || fill->role() == ExtrusionRole::erBridgeInfill) {
                        issues.add(
                                check_extrusion_entity_stability(fill, layer->print_z, layer_region,
                                        supported_grid,
                                        params));
                    }
                } // fill
            } // ex_entity
        } // region

    } else { // If NOT full check, check only external perimeters
        for (const LayerRegion *layer_region : layer->regions()) {
            for (const ExtrusionEntity *ex_entity : layer_region->perimeters.entities) {
                for (const ExtrusionEntity *perimeter : static_cast<const ExtrusionEntityCollection*>(ex_entity)->entities) {
                    if (perimeter->role() == ExtrusionRole::erExternalPerimeter
                            || perimeter->role() == ExtrusionRole::erOverhangPerimeter) {
                        issues.add(
                                check_extrusion_entity_stability(perimeter, layer->print_z, layer_region,
                                        supported_grid, params));
                    }; // ex_perimeter
                } // perimeter
            } // ex_entity
        } //region
    }

    return issues;
}

} //Impl End

std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
    using namespace Impl;

    BalanceDistributionGrid grid { };
    grid.init(po, 0, po->layers().size());
    distribute_layer_volume(po, 0, grid);
    std::mutex grid_mutex;

    size_t layer_count = po->layer_count();
    std::vector<bool> layer_needs_supports(layer_count, false);
    tbb::parallel_for(tbb::blocked_range<size_t>(1, layer_count), [&](tbb::blocked_range<size_t> r) {
        BalanceDistributionGrid balance_grid { };
        balance_grid.init(po, r.begin(), r.end());

        for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
            distribute_layer_volume(po, layer_idx, balance_grid);
            auto layer_issues = check_layer_stability(po, layer_idx, false, params);
            if (!layer_issues.supports_nedded.empty()) {
                layer_needs_supports[layer_idx] = true;
            }
        }

        grid_mutex.lock();
        grid.merge(balance_grid);
        grid_mutex.unlock();
    });

    std::vector<size_t> problematic_layers;
    for (size_t index = 0; index < layer_needs_supports.size(); ++index) {
        if (layer_needs_supports[index]) {
            problematic_layers.push_back(index);
        }
    }
    return problematic_layers;
}

Issues full_search(const PrintObject *po, const Params &params) {
    using namespace Impl;

    BalanceDistributionGrid grid { };
    grid.init(po, 0, po->layers().size());
    distribute_layer_volume(po, 0, grid);
    std::mutex grid_mutex;

    size_t layer_count = po->layer_count();
    Issues found_issues = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, layer_count), Issues { },
            [&](tbb::blocked_range<size_t> r, const Issues &init) {
                BalanceDistributionGrid balance_grid { };
                balance_grid.init(po, r.begin(), r.end());
                Issues issues = init;
                for (size_t layer_idx = r.begin(); layer_idx < r.end(); ++layer_idx) {
                    distribute_layer_volume(po, layer_idx, balance_grid);
                    auto layer_issues = check_layer_stability(po, layer_idx, true, params);
                    if (!layer_issues.empty()) {
                        issues.add(layer_issues);
                    }
                }

                grid_mutex.lock();
                grid.merge(balance_grid);
                grid_mutex.unlock();

                return issues;
            },
            [](Issues left, const Issues &right) {
                left.add(right);
                return left;
            }
    );
#ifdef DEBUG_FILES
    Impl::debug_export(found_issues, "pre_issues");
#endif

    grid.analyze(found_issues, params);

#ifdef DEBUG_FILES
    grid.debug_export();
    Impl::debug_export(found_issues, "issues");
#endif

    return found_issues;
}

}
}
