#include "SupportSpotsGenerator.hpp"

#include "ExPolygon.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "GCode/ExtrusionProcessor.hpp"
#include "Line.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Print.hpp"
#include "PrintBase.hpp"
#include "Tesselate.hpp"
#include "libslic3r.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/parallel_reduce.h"
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <utility>
#include <vector>

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
    ExtrusionLine() : a(Vec2f::Zero()), b(Vec2f::Zero()), len(0.0), origin_entity(nullptr) {}
    ExtrusionLine(const Vec2f &a, const Vec2f &b, float len, const ExtrusionEntity *origin_entity)
        : a(a), b(b), len(len), origin_entity(origin_entity)
    {}

    ExtrusionLine(const Vec2f &a, const Vec2f &b)
        : a(a), b(b), len((a-b).norm()), origin_entity(nullptr)
    {}

    bool is_external_perimeter() const
    {
        assert(origin_entity != nullptr);
        return origin_entity->role().is_external_perimeter();
    }

    Vec2f                  a;
    Vec2f                  b;
    float                  len;
    const ExtrusionEntity *origin_entity;

    bool  support_point_generated = false;
    float form_quality            = 1.0f;
    float curled_up_height        = 0.0f;

    static const constexpr int Dim = 2;
    using Scalar                   = Vec2f::Scalar;
};

auto get_a(ExtrusionLine &&l) { return l.a; }
auto get_b(ExtrusionLine &&l) { return l.b; }

namespace SupportSpotsGenerator {

SupportPoint::SupportPoint(const Vec3f &position, float force, float spot_radius, const Vec2f &direction)
    : position(position), force(force), spot_radius(spot_radius), direction(direction)
{}

using LD = AABBTreeLines::LinesDistancer<ExtrusionLine>;

struct SupportGridFilter
{
private:
    Vec3f cell_size;
    Vec3f origin;
    Vec3f size;
    Vec3i cell_count;

    std::unordered_set<size_t> taken_cells{};

public:
    SupportGridFilter(const PrintObject *po, float voxel_size)
    {
        cell_size = Vec3f(voxel_size, voxel_size, voxel_size);

        Vec2crd size_half = po->size().head<2>().cwiseQuotient(Vec2crd(2, 2)) + Vec2crd::Ones();
        Vec3f   min       = unscale(Vec3crd(-size_half.x(), -size_half.y(), 0)).cast<float>() - cell_size;
        Vec3f   max       = unscale(Vec3crd(size_half.x(), size_half.y(), po->height())).cast<float>() + cell_size;

        origin     = min;
        size       = max - min;
        cell_count = size.cwiseQuotient(cell_size).cast<int>() + Vec3i::Ones();
    }

    Vec3i to_cell_coords(const Vec3f &position) const
    {
        Vec3i cell_coords = (position - this->origin).cwiseQuotient(this->cell_size).cast<int>();
        return cell_coords;
    }

    size_t to_cell_index(const Vec3i &cell_coords) const
    {
        assert(cell_coords.x() >= 0);
        assert(cell_coords.x() < cell_count.x());
        assert(cell_coords.y() >= 0);
        assert(cell_coords.y() < cell_count.y());
        assert(cell_coords.z() >= 0);
        assert(cell_coords.z() < cell_count.z());

        return cell_coords.z() * cell_count.x() * cell_count.y() + cell_coords.y() * cell_count.x() + cell_coords.x();
    }

    Vec3f get_cell_center(const Vec3i &cell_coords) const
    {
        return origin + cell_coords.cast<float>().cwiseProduct(this->cell_size) + this->cell_size.cwiseQuotient(Vec3f(2.0f, 2.0f, 2.0f));
    }

    void take_position(const Vec3f &position) { taken_cells.insert(to_cell_index(to_cell_coords(position))); }

    bool position_taken(const Vec3f &position) const
    {
        return taken_cells.find(to_cell_index(to_cell_coords(position))) != taken_cells.end();
    }
};

struct SliceConnection
{
    float area{};
    Vec3f centroid_accumulator              = Vec3f::Zero();
    Vec2f second_moment_of_area_accumulator = Vec2f::Zero();
    float second_moment_of_area_covariance_accumulator{};

    void add(const SliceConnection &other)
    {
        this->area += other.area;
        this->centroid_accumulator += other.centroid_accumulator;
        this->second_moment_of_area_accumulator += other.second_moment_of_area_accumulator;
        this->second_moment_of_area_covariance_accumulator += other.second_moment_of_area_covariance_accumulator;
    }

    void print_info(const std::string &tag)
    {
        Vec3f centroid   = centroid_accumulator / area;
        Vec2f variance   = (second_moment_of_area_accumulator / area - centroid.head<2>().cwiseProduct(centroid.head<2>()));
        float covariance = second_moment_of_area_covariance_accumulator / area - centroid.x() * centroid.y();
        std::cout << tag << std::endl;
        std::cout << "area: " << area << std::endl;
        std::cout << "centroid: " << centroid.x() << " " << centroid.y() << " " << centroid.z() << std::endl;
        std::cout << "variance: " << variance.x() << " " << variance.y() << std::endl;
        std::cout << "covariance: " << covariance << std::endl;
    }
};

float get_flow_width(const LayerRegion *region, ExtrusionRole role)
{
    if (role == ExtrusionRole::BridgeInfill) return region->flow(FlowRole::frExternalPerimeter).width();
    if (role == ExtrusionRole::ExternalPerimeter) return region->flow(FlowRole::frExternalPerimeter).width();
    if (role == ExtrusionRole::GapFill) return region->flow(FlowRole::frInfill).width();
    if (role == ExtrusionRole::Perimeter) return region->flow(FlowRole::frPerimeter).width();
    if (role == ExtrusionRole::SolidInfill) return region->flow(FlowRole::frSolidInfill).width();
    if (role == ExtrusionRole::InternalInfill) return region->flow(FlowRole::frInfill).width();
    if (role == ExtrusionRole::TopSolidInfill) return region->flow(FlowRole::frTopSolidInfill).width();
    // default
    return region->flow(FlowRole::frPerimeter).width();
}

std::vector<ExtrusionLine> to_short_lines(const ExtrusionEntity *e, float length_limit)
{
    assert(!e->is_collection());
    Polyline                   pl = e->as_polyline();
    std::vector<ExtrusionLine> lines;
    lines.reserve(pl.points.size() * 1.5f);
    for (int point_idx = 0; point_idx < int(pl.points.size()) - 1; ++point_idx) {
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
            lines.emplace_back(a, b, (a-b).norm(), e);
        }
    }
    return lines;
}

float estimate_curled_up_height(
    const ExtendedPoint &point, float layer_height, float flow_width, float prev_line_curled_height, Params params)
{
    float curled_up_height = 0.0f;
    if (fabs(point.distance) < 1.5 * flow_width) {
        curled_up_height = 0.85 * prev_line_curled_height;
    }
    if (point.distance > params.malformation_distance_factors.first * flow_width &&
        point.distance < params.malformation_distance_factors.second * flow_width && point.curvature > -0.1f) {
        float dist_factor = std::max(point.distance - params.malformation_distance_factors.first * flow_width, 0.01f) /
                            ((params.malformation_distance_factors.second - params.malformation_distance_factors.first) * flow_width);

        curled_up_height = layer_height * 2.0f * sqrt(sqrt(dist_factor)) * std::clamp(6.0f * point.curvature, 1.0f, 6.0f);
        curled_up_height = std::min(curled_up_height, params.max_curled_height_factor * layer_height);
    }

    return curled_up_height;
}

std::vector<ExtrusionLine> check_extrusion_entity_stability(const ExtrusionEntity                      *entity,
                                                            const LayerRegion                          *layer_region,
                                                            const LD                                   &prev_layer_lines,
                                                            const AABBTreeLines::LinesDistancer<Linef> &prev_layer_boundary,
                                                            const Params                               &params)
{
    if (entity->is_collection()) {
        std::vector<ExtrusionLine> checked_lines_out;
        checked_lines_out.reserve(prev_layer_lines.get_lines().size() / 3);
        for (const auto *e : static_cast<const ExtrusionEntityCollection *>(entity)->entities) {
            auto tmp = check_extrusion_entity_stability(e, layer_region, prev_layer_lines, prev_layer_boundary, params);
            checked_lines_out.insert(checked_lines_out.end(), tmp.begin(), tmp.end());
        }
        return checked_lines_out;
    } else { // single extrusion path, with possible varying parameters
        if (entity->length() < scale_(params.min_distance_to_allow_local_supports)) {
            return {};
        }

        const float flow_width = get_flow_width(layer_region, entity->role());

        std::vector<ExtendedPoint> annotated_points = estimate_points_properties<true, true, false, false>(entity->as_polyline().points,
                                                                                                           prev_layer_lines, flow_width,
                                                                                                           params.bridge_distance);

        std::vector<ExtrusionLine> lines_out;
        lines_out.reserve(annotated_points.size());
        float bridged_distance = annotated_points.front().position != annotated_points.back().position ? (params.bridge_distance + 1.0f) :
                                                                                                         0.0f;
        for (size_t i = 0; i < annotated_points.size(); ++i) {
            ExtendedPoint &curr_point = annotated_points[i];
            float          line_len   = i > 0 ? ((annotated_points[i - 1].position - curr_point.position).norm()) : 0.0f;
            ExtrusionLine  line_out{i > 0 ? annotated_points[i - 1].position.cast<float>() : curr_point.position.cast<float>(),
                                   curr_point.position.cast<float>(), line_len, entity};

            const ExtrusionLine nearest_prev_layer_line = prev_layer_lines.get_lines().size() > 0 ?
                                                        prev_layer_lines.get_line(curr_point.nearest_prev_layer_line) :
                                                        ExtrusionLine{};

            float sign = (prev_layer_boundary.distance_from_lines<true>(curr_point.position) + 0.5f * flow_width) < 0.0f ? -1.0f : 1.0f;
            curr_point.distance *= sign;

            float max_bridge_len = params.bridge_distance /
                                   ((1.0f + std::abs(curr_point.curvature)) * (1.0f + std::abs(curr_point.curvature)));

            if (curr_point.distance > 2.0f * flow_width) {
                line_out.form_quality = 0.8f;
                bridged_distance += line_len;
                if (bridged_distance > max_bridge_len) {
                    line_out.support_point_generated = true;
                    bridged_distance                 = 0.0f;
                }
            } else if (curr_point.distance > flow_width * (1.0 + std::clamp(curr_point.curvature, -0.30f, 0.20f))) {
                bridged_distance += line_len;
                line_out.form_quality = nearest_prev_layer_line.form_quality - 0.3f;
                if (line_out.form_quality < 0 && bridged_distance > max_bridge_len) {
                    line_out.support_point_generated = true;
                    line_out.form_quality            = 0.5f;
                    bridged_distance                 = 0.0f;
                }
            } else {
                bridged_distance = 0.0f;
            }

            line_out.curled_up_height = estimate_curled_up_height(curr_point, layer_region->layer()->height, flow_width,
                                                                  nearest_prev_layer_line.curled_up_height, params);

            lines_out.push_back(line_out);
        }

        return lines_out;
    }
}

// returns triangle area, first_moment_of_area_xy, second_moment_of_area_xy, second_moment_of_area_covariance
// none of the values is divided/normalized by area.
// The function computes intgeral over the area of the triangle, with function f(x,y) = x for first moments of area (y is analogous)
// f(x,y) = x^2 for second moment of area
// and f(x,y) = x*y for second moment of area covariance
std::tuple<float, Vec2f, Vec2f, float> compute_triangle_moments_of_area(const Vec2f &a, const Vec2f &b, const Vec2f &c)
{
    // based on the following guide:
    // Denote the vertices of S by a, b, c. Then the map
    //  g:(u,v)↦a+u(b−a)+v(c−a) ,
    //  which in coordinates appears as
    //  g:(u,v)↦{x(u,v)y(u,v)=a1+u(b1−a1)+v(c1−a1)=a2+u(b2−a2)+v(c2−a2) ,(1)
    //  obviously maps S′ bijectively onto S. Therefore the transformation formula for multiple integrals steps into action, and we obtain
    //  ∫Sf(x,y)d(x,y)=∫S′f(x(u,v),y(u,v))∣∣Jg(u,v)∣∣ d(u,v) .
    //  In the case at hand the Jacobian determinant is a constant: From (1) we obtain
    //  Jg(u,v)=det[xuyuxvyv]=(b1−a1)(c2−a2)−(c1−a1)(b2−a2) .
    //  Therefore we can write
    //  ∫Sf(x,y)d(x,y)=∣∣Jg∣∣∫10∫1−u0f~(u,v) dv du ,
    //  where f~ denotes the pullback of f to S′:
    //  f~(u,v):=f(x(u,v),y(u,v)) .
    //  Don't forget taking the absolute value of Jg!

    float jacobian_determinant_abs = std::abs((b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()));

    // coordinate transform: gx(u,v) = a.x + u * (b.x - a.x) + v * (c.x - a.x)
    // coordinate transform: gy(u,v) = a.y + u * (b.y - a.y) + v * (c.y - a.y)
    // second moment of area for x: f(x, y) = x^2;
    //              f(gx(u,v), gy(u,v)) = gx(u,v)^2 = ... (long expanded form)

    // result is Int_T func = jacobian_determinant_abs * Int_0^1 Int_0^1-u func(gx(u,v), gy(u,v)) dv du
    // integral_0^1 integral_0^(1 - u) (a + u (b - a) + v (c - a))^2 dv du = 1/12 (a^2 + a (b + c) + b^2 + b c + c^2)

    Vec2f second_moment_of_area_xy = jacobian_determinant_abs *
                                     (a.cwiseProduct(a) + b.cwiseProduct(b) + b.cwiseProduct(c) + c.cwiseProduct(c) +
                                      a.cwiseProduct(b + c)) /
                                     12.0f;
    // second moment of area covariance : f(x, y) = x*y;
    //              f(gx(u,v), gy(u,v)) = gx(u,v)*gy(u,v) = ... (long expanded form)
    //(a_1 + u * (b_1 - a_1) + v * (c_1 - a_1)) * (a_2 + u * (b_2 - a_2) + v * (c_2 - a_2))
    // ==    (a_1 + u (b_1 - a_1) + v (c_1 - a_1)) (a_2 + u (b_2 - a_2) + v (c_2 - a_2))

    // intermediate result: integral_0^(1 - u) (a_1 + u (b_1 - a_1) + v (c_1 - a_1)) (a_2 + u (b_2 - a_2) + v (c_2 - a_2)) dv =
    //  1/6 (u - 1) (-c_1 (u - 1) (a_2 (u - 1) - 3 b_2 u) - c_2 (u - 1) (a_1 (u - 1) - 3 b_1 u + 2 c_1 (u - 1)) + 3 b_1 u (a_2 (u - 1) - 2
    //  b_2 u) + a_1 (u - 1) (3 b_2 u - 2 a_2 (u - 1))) result = integral_0^1 1/6 (u - 1) (-c_1 (u - 1) (a_2 (u - 1) - 3 b_2 u) - c_2 (u -
    //  1) (a_1 (u - 1) - 3 b_1 u + 2 c_1 (u - 1)) + 3 b_1 u (a_2 (u - 1) - 2 b_2 u) + a_1 (u - 1) (3 b_2 u - 2 a_2 (u - 1))) du =
    //   1/24 (a_2 (b_1 + c_1) + a_1 (2 a_2 + b_2 + c_2) + b_2 c_1 + b_1 c_2 + 2 b_1 b_2 + 2 c_1 c_2)
    //  result is Int_T func = jacobian_determinant_abs * Int_0^1 Int_0^1-u func(gx(u,v), gy(u,v)) dv du
    float second_moment_of_area_covariance = jacobian_determinant_abs * (1.0f / 24.0f) *
                                             (a.y() * (b.x() + c.x()) + a.x() * (2.0f * a.y() + b.y() + c.y()) + b.y() * c.x() +
                                              b.x() * c.y() + 2.0f * b.x() * b.y() + 2.0f * c.x() * c.y());

    float area = jacobian_determinant_abs * 0.5f;

    Vec2f first_moment_of_area_xy = jacobian_determinant_abs * (a + b + c) / 6.0f;

    return {area, first_moment_of_area_xy, second_moment_of_area_xy, second_moment_of_area_covariance};
};

SliceConnection estimate_slice_connection(size_t slice_idx, const Layer *layer)
{
    SliceConnection connection;

    const LayerSlice   &slice       = layer->lslices_ex[slice_idx];
    ExPolygon           slice_poly  = layer->lslices[slice_idx];
    const Layer        *lower_layer = layer->lower_layer;

    ExPolygons below_polys{};
    for (const auto &link : slice.overlaps_below) { below_polys.push_back(lower_layer->lslices[link.slice_idx]); }
    ExPolygons overlap = intersection_ex({slice_poly}, below_polys);

    std::vector<Vec2f> triangles = triangulate_expolygons_2f(overlap);
    for (size_t idx = 0; idx < triangles.size(); idx += 3) {
        auto [area, first_moment_of_area, second_moment_area,
              second_moment_of_area_covariance] = compute_triangle_moments_of_area(triangles[idx], triangles[idx + 1], triangles[idx + 2]);
        connection.area += area;
        connection.centroid_accumulator += Vec3f(first_moment_of_area.x(), first_moment_of_area.y(), layer->print_z * area);
        connection.second_moment_of_area_accumulator += second_moment_area;
        connection.second_moment_of_area_covariance_accumulator += second_moment_of_area_covariance;
    }

    return connection;
};

class ObjectPart
{
public:
    float volume{};
    Vec3f volume_centroid_accumulator = Vec3f::Zero();
    float sticking_area{};
    Vec3f sticking_centroid_accumulator              = Vec3f::Zero();
    Vec2f sticking_second_moment_of_area_accumulator = Vec2f::Zero();
    float sticking_second_moment_of_area_covariance_accumulator{};

    ObjectPart() = default;

    void add(const ObjectPart &other)
    {
        this->volume_centroid_accumulator += other.volume_centroid_accumulator;
        this->volume += other.volume;
        this->sticking_area += other.sticking_area;
        this->sticking_centroid_accumulator += other.sticking_centroid_accumulator;
        this->sticking_second_moment_of_area_accumulator += other.sticking_second_moment_of_area_accumulator;
        this->sticking_second_moment_of_area_covariance_accumulator += other.sticking_second_moment_of_area_covariance_accumulator;
    }

    void add_support_point(const Vec3f &position, float sticking_area)
    {
        this->sticking_area += sticking_area;
        this->sticking_centroid_accumulator += sticking_area * position;
        this->sticking_second_moment_of_area_accumulator += sticking_area * position.head<2>().cwiseProduct(position.head<2>());
        this->sticking_second_moment_of_area_covariance_accumulator += sticking_area * position.x() * position.y();
    }

    float compute_directional_xy_variance(const Vec2f &line_dir,
                                          const Vec3f &centroid_accumulator,
                                          const Vec2f &second_moment_of_area_accumulator,
                                          const float &second_moment_of_area_covariance_accumulator,
                                          const float &area) const
    {
        assert(area > 0);
        Vec3f centroid   = centroid_accumulator / area;
        Vec2f variance   = (second_moment_of_area_accumulator / area - centroid.head<2>().cwiseProduct(centroid.head<2>()));
        float covariance = second_moment_of_area_covariance_accumulator / area - centroid.x() * centroid.y();
        // Var(aX+bY)=a^2*Var(X)+b^2*Var(Y)+2*a*b*Cov(X,Y)
        float directional_xy_variance = line_dir.x() * line_dir.x() * variance.x() + line_dir.y() * line_dir.y() * variance.y() +
                                        2.0f * line_dir.x() * line_dir.y() * covariance;
#ifdef DETAILED_DEBUG_LOGS
        BOOST_LOG_TRIVIAL(debug) << "centroid: " << centroid.x() << "  " << centroid.y() << "  " << centroid.z();
        BOOST_LOG_TRIVIAL(debug) << "variance: " << variance.x() << "  " << variance.y();
        BOOST_LOG_TRIVIAL(debug) << "covariance: " << covariance;
        BOOST_LOG_TRIVIAL(debug) << "directional_xy_variance: " << directional_xy_variance;
#endif
        return directional_xy_variance;
    }

    float compute_elastic_section_modulus(const Vec2f &line_dir,
                                          const Vec3f &extreme_point,
                                          const Vec3f &centroid_accumulator,
                                          const Vec2f &second_moment_of_area_accumulator,
                                          const float &second_moment_of_area_covariance_accumulator,
                                          const float &area) const
    {
        float directional_xy_variance = compute_directional_xy_variance(line_dir, centroid_accumulator, second_moment_of_area_accumulator,
                                                                        second_moment_of_area_covariance_accumulator, area);
        if (directional_xy_variance < EPSILON) { return 0.0f; }
        Vec3f centroid                = centroid_accumulator / area;
        float extreme_fiber_dist      = line_alg::distance_to(Linef(centroid.head<2>().cast<double>(),
                                                                    (centroid.head<2>() + Vec2f(line_dir.y(), -line_dir.x())).cast<double>()),
                                                              extreme_point.head<2>().cast<double>());
        float elastic_section_modulus = area * directional_xy_variance / extreme_fiber_dist;

#ifdef DETAILED_DEBUG_LOGS
        BOOST_LOG_TRIVIAL(debug) << "extreme_fiber_dist: " << extreme_fiber_dist;
        BOOST_LOG_TRIVIAL(debug) << "elastic_section_modulus: " << elastic_section_modulus;
#endif

        return elastic_section_modulus;
    }

    float is_stable_while_extruding(const SliceConnection &connection,
                                    const ExtrusionLine   &extruded_line,
                                    const Vec3f           &extreme_point,
                                    float                  layer_z,
                                    const Params          &params) const
    {
        Vec2f        line_dir      = (extruded_line.b - extruded_line.a).normalized();
        const Vec3f &mass_centroid = this->volume_centroid_accumulator / this->volume;
        float        mass          = this->volume * params.filament_density;
        float        weight        = mass * params.gravity_constant;

        float movement_force = params.max_acceleration * mass;

        float extruder_conflict_force = params.standard_extruder_conflict_force +
                                        std::min(extruded_line.curled_up_height, 1.0f) * params.malformations_additive_conflict_extruder_force;

        // section for bed calculations
        {
            if (this->sticking_area < EPSILON) return 1.0f;

            Vec3f bed_centroid     = this->sticking_centroid_accumulator / this->sticking_area;
            float bed_yield_torque = -compute_elastic_section_modulus(line_dir, extreme_point, this->sticking_centroid_accumulator,
                                                                      this->sticking_second_moment_of_area_accumulator,
                                                                      this->sticking_second_moment_of_area_covariance_accumulator,
                                                                      this->sticking_area) *
                                     params.get_bed_adhesion_yield_strength();

            Vec2f bed_weight_arm             = (mass_centroid.head<2>() - bed_centroid.head<2>());
            float bed_weight_arm_len         = bed_weight_arm.norm();
            float bed_weight_dir_xy_variance = compute_directional_xy_variance(bed_weight_arm, this->sticking_centroid_accumulator,
                                                                               this->sticking_second_moment_of_area_accumulator,
                                                                               this->sticking_second_moment_of_area_covariance_accumulator,
                                                                               this->sticking_area);
            float bed_weight_sign            = bed_weight_arm_len < 2.0f * sqrt(bed_weight_dir_xy_variance) ? -1.0f : 1.0f;
            float bed_weight_torque          = bed_weight_sign * bed_weight_arm_len * weight;

            float bed_movement_arm    = std::max(0.0f, mass_centroid.z() - bed_centroid.z());
            float bed_movement_torque = movement_force * bed_movement_arm;

            float bed_conflict_torque_arm      = layer_z - bed_centroid.z();
            float bed_extruder_conflict_torque = extruder_conflict_force * bed_conflict_torque_arm;

            float bed_total_torque = bed_movement_torque + bed_extruder_conflict_torque + bed_weight_torque + bed_yield_torque;

#ifdef DETAILED_DEBUG_LOGS
            BOOST_LOG_TRIVIAL(debug) << "bed_centroid: " << bed_centroid.x() << "  " << bed_centroid.y() << "  " << bed_centroid.z();
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_yield_torque: " << bed_yield_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_weight_arm: " << bed_weight_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_weight_torque: " << bed_weight_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_movement_arm: " << bed_movement_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_movement_torque: " << bed_movement_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_conflict_torque_arm: " << bed_conflict_torque_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: extruded_line.curled_up_height: " << extruded_line.curled_up_height;
            BOOST_LOG_TRIVIAL(debug) << "SSG: extruded_line.form_quality: " << extruded_line.form_quality;
            BOOST_LOG_TRIVIAL(debug) << "SSG: extruder_conflict_force: " << extruder_conflict_force;
            BOOST_LOG_TRIVIAL(debug) << "SSG: bed_extruder_conflict_torque: " << bed_extruder_conflict_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: total_torque: " << bed_total_torque << "   layer_z: " << layer_z;
#endif

            if (bed_total_torque > 0) return bed_total_torque / bed_conflict_torque_arm;
        }

        // section for weak connection calculations
        {
            if (connection.area < EPSILON) return 1.0f;

            Vec3f conn_centroid = connection.centroid_accumulator / connection.area;

            if (layer_z - conn_centroid.z() < 3.0f) { return -1.0f; }
            float conn_yield_torque = compute_elastic_section_modulus(line_dir, extreme_point, connection.centroid_accumulator,
                                                                      connection.second_moment_of_area_accumulator,
                                                                      connection.second_moment_of_area_covariance_accumulator,
                                                                      connection.area) *
                                      params.material_yield_strength;

            float conn_weight_arm    = (conn_centroid.head<2>() - mass_centroid.head<2>()).norm();
            float conn_weight_torque = conn_weight_arm * weight * (conn_centroid.z() / layer_z);

            float conn_movement_arm    = std::max(0.0f, mass_centroid.z() - conn_centroid.z());
            float conn_movement_torque = movement_force * conn_movement_arm;

            float conn_conflict_torque_arm      = layer_z - conn_centroid.z();
            float conn_extruder_conflict_torque = extruder_conflict_force * conn_conflict_torque_arm;

            float conn_total_torque = conn_movement_torque + conn_extruder_conflict_torque + conn_weight_torque - conn_yield_torque;

#ifdef DETAILED_DEBUG_LOGS
            BOOST_LOG_TRIVIAL(debug) << "conn_centroid: " << conn_centroid.x() << "  " << conn_centroid.y() << "  " << conn_centroid.z();
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_yield_torque: " << conn_yield_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_weight_arm: " << conn_weight_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_weight_torque: " << conn_weight_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_movement_arm: " << conn_movement_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_movement_torque: " << conn_movement_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_conflict_torque_arm: " << conn_conflict_torque_arm;
            BOOST_LOG_TRIVIAL(debug) << "SSG: conn_extruder_conflict_torque: " << conn_extruder_conflict_torque;
            BOOST_LOG_TRIVIAL(debug) << "SSG: total_torque: " << conn_total_torque << "   layer_z: " << layer_z;
#endif

            return conn_total_torque / conn_conflict_torque_arm;
        }
    }
};

// return new object part and actual area covered by extrusions
std::tuple<ObjectPart, float> build_object_part_from_slice(const LayerSlice &slice, const Layer *layer)
{
    ObjectPart new_object_part;
    float      area_covered_by_extrusions = 0;

    auto add_extrusions_to_object = [&new_object_part, &area_covered_by_extrusions](const ExtrusionEntity *e, const LayerRegion *region) {
        float                      flow_width = get_flow_width(region, e->role());
        const Layer               *l          = region->layer();
        float                      slice_z    = l->slice_z;
        float                      height     = l->height;
        std::vector<ExtrusionLine> lines      = to_short_lines(e, 5.0);
        for (const ExtrusionLine &line : lines) {
            float volume = line.len * height * flow_width * PI / 4.0f;
            area_covered_by_extrusions += line.len * flow_width;
            new_object_part.volume += volume;
            new_object_part.volume_centroid_accumulator += to_3d(Vec2f((line.a + line.b) / 2.0f), slice_z) * volume;

            if (l->bottom_z() < EPSILON) { // layer attached on bed
                float sticking_area = line.len * flow_width;
                new_object_part.sticking_area += sticking_area;
                Vec2f middle = Vec2f((line.a + line.b) / 2.0f);
                new_object_part.sticking_centroid_accumulator += sticking_area * to_3d(middle, slice_z);
                // Bottom infill lines can be quite long, and algined, so the middle approximaton used above does not work
                Vec2f dir            = (line.b - line.a).normalized();
                float segment_length = flow_width; // segments of size flow_width
                for (float segment_middle_dist = std::min(line.len, segment_length * 0.5f); segment_middle_dist < line.len;
                     segment_middle_dist += segment_length) {
                    Vec2f segment_middle = line.a + segment_middle_dist * dir;
                    new_object_part.sticking_second_moment_of_area_accumulator += segment_length * flow_width *
                                                                                  segment_middle.cwiseProduct(segment_middle);
                    new_object_part.sticking_second_moment_of_area_covariance_accumulator += segment_length * flow_width *
                                                                                             segment_middle.x() * segment_middle.y();
                }
            }
        }
    };

    for (const auto &island : slice.islands) {
        const LayerRegion *perimeter_region = layer->get_region(island.perimeters.region());
        for (const auto &perimeter_idx : island.perimeters) {
            for (const ExtrusionEntity *perimeter :
                 static_cast<const ExtrusionEntityCollection *>(perimeter_region->perimeters().entities[perimeter_idx])->entities) {
                add_extrusions_to_object(perimeter, perimeter_region);
            }
        }
        for (const LayerExtrusionRange &fill_range : island.fills) {
            const LayerRegion *fill_region = layer->get_region(fill_range.region());
            for (const auto &fill_idx : fill_range) {
                for (const ExtrusionEntity *fill :
                     static_cast<const ExtrusionEntityCollection *>(fill_region->fills().entities[fill_idx])->entities) {
                    add_extrusions_to_object(fill, fill_region);
                }
            }
        }
        for (const auto &thin_fill_idx : island.thin_fills) {
            add_extrusions_to_object(perimeter_region->thin_fills().entities[thin_fill_idx], perimeter_region);
        }
    }

    return {new_object_part, area_covered_by_extrusions};
}

class ActiveObjectParts
{
    size_t                                 next_part_idx = 0;
    std::unordered_map<size_t, ObjectPart> active_object_parts;
    std::unordered_map<size_t, size_t>     active_object_parts_id_mapping;

public:
    size_t get_flat_id(size_t id)
    {
        size_t index = active_object_parts_id_mapping.at(id);
        while (index != active_object_parts_id_mapping.at(index)) { index = active_object_parts_id_mapping.at(index); }
        size_t i = id;
        while (index != active_object_parts_id_mapping.at(i)) {
            size_t next                       = active_object_parts_id_mapping[i];
            active_object_parts_id_mapping[i] = index;
            i                                 = next;
        }
        return index;
    }

    ObjectPart &access(size_t id) { return this->active_object_parts.at(this->get_flat_id(id)); }

    size_t insert(const ObjectPart &new_part)
    {
        this->active_object_parts.emplace(next_part_idx, new_part);
        this->active_object_parts_id_mapping.emplace(next_part_idx, next_part_idx);
        return next_part_idx++;
    }

    void merge(size_t from, size_t to)
    {
        size_t to_flat   = this->get_flat_id(to);
        size_t from_flat = this->get_flat_id(from);
        active_object_parts.at(to_flat).add(active_object_parts.at(from_flat));
        active_object_parts.erase(from_flat);
        active_object_parts_id_mapping[from] = to_flat;
    }
};

SupportPoints check_stability(const PrintObject *po, const PrintTryCancel& cancel_func, const Params &params)
{
    SupportPoints            supp_points{};
    SupportGridFilter supports_presence_grid(po, params.min_distance_between_support_points);
    ActiveObjectParts active_object_parts{};
    LD                prev_layer_ext_perim_lines;

    std::unordered_map<size_t, size_t>          prev_slice_idx_to_object_part_mapping;
    std::unordered_map<size_t, size_t>          next_slice_idx_to_object_part_mapping;
    std::unordered_map<size_t, SliceConnection> prev_slice_idx_to_weakest_connection;
    std::unordered_map<size_t, SliceConnection> next_slice_idx_to_weakest_connection;

    for (size_t layer_idx = 0; layer_idx < po->layer_count(); ++layer_idx) {
        cancel_func();
        const Layer *layer                 = po->get_layer(layer_idx);
        float        bottom_z              = layer->bottom_z();
        auto create_support_point_position = [bottom_z](const Vec2f &layer_pos) { return Vec3f{layer_pos.x(), layer_pos.y(), bottom_z}; };

        for (size_t slice_idx = 0; slice_idx < layer->lslices_ex.size(); ++slice_idx) {
            const LayerSlice &slice             = layer->lslices_ex.at(slice_idx);
            auto [new_part, covered_area]       = build_object_part_from_slice(slice, layer);
            SliceConnection connection_to_below = estimate_slice_connection(slice_idx, layer);

#ifdef DETAILED_DEBUG_LOGS
            std::cout << "SLICE IDX: " << slice_idx << std::endl;
            for (const auto &link : slice.overlaps_below) {
                std::cout << "connected to slice below: " << link.slice_idx << "  by area : " << link.area << std::endl;
            }
            connection_to_below.print_info("CONNECTION TO BELOW");
#endif

            if (connection_to_below.area < EPSILON) { // new object part emerging
                size_t part_id = active_object_parts.insert(new_part);
                next_slice_idx_to_object_part_mapping.emplace(slice_idx, part_id);
                next_slice_idx_to_weakest_connection.emplace(slice_idx, connection_to_below);
            } else {
                size_t          final_part_id{};
                SliceConnection transfered_weakest_connection{};
                // MERGE parts
                {
                    std::unordered_set<size_t> parts_ids;
                    for (const auto &link : slice.overlaps_below) {
                        size_t part_id = active_object_parts.get_flat_id(prev_slice_idx_to_object_part_mapping.at(link.slice_idx));
                        parts_ids.insert(part_id);
                        transfered_weakest_connection.add(prev_slice_idx_to_weakest_connection.at(link.slice_idx));
                    }

                    final_part_id = *parts_ids.begin();
                    for (size_t part_id : parts_ids) {
                        if (final_part_id != part_id) { active_object_parts.merge(part_id, final_part_id); }
                    }
                }
                auto estimate_conn_strength = [bottom_z](const SliceConnection &conn) {
                    if (conn.area < EPSILON) { // connection is empty, does not exists. Return max strength so that it is not picked as the
                                               // weakest connection.
                        return INFINITY;
                    }
                    Vec3f centroid         = conn.centroid_accumulator / conn.area;
                    Vec2f variance         = (conn.second_moment_of_area_accumulator / conn.area -
                                      centroid.head<2>().cwiseProduct(centroid.head<2>()));
                    float xy_variance      = variance.x() + variance.y();
                    float arm_len_estimate = std::max(1.0f, bottom_z - (conn.centroid_accumulator.z() / conn.area));
                    return conn.area * sqrt(xy_variance) / arm_len_estimate;
                };

#ifdef DETAILED_DEBUG_LOGS
                connection_to_below.print_info("new_weakest_connection");
                transfered_weakest_connection.print_info("transfered_weakest_connection");
#endif

                if (estimate_conn_strength(transfered_weakest_connection) > estimate_conn_strength(connection_to_below)) {
                    transfered_weakest_connection = connection_to_below;
                }
                next_slice_idx_to_weakest_connection.emplace(slice_idx, transfered_weakest_connection);
                next_slice_idx_to_object_part_mapping.emplace(slice_idx, final_part_id);
                ObjectPart &part = active_object_parts.access(final_part_id);
                part.add(new_part);
            }
        }

        prev_slice_idx_to_object_part_mapping = next_slice_idx_to_object_part_mapping;
        next_slice_idx_to_object_part_mapping.clear();
        prev_slice_idx_to_weakest_connection = next_slice_idx_to_weakest_connection;
        next_slice_idx_to_weakest_connection.clear();

        std::vector<ExtrusionLine> current_layer_ext_perims_lines{};
        current_layer_ext_perims_lines.reserve(prev_layer_ext_perim_lines.get_lines().size());
        // All object parts updated, and for each slice we have coresponding weakest connection.
        // We can now check each slice and its corresponding weakest connection and object part for stability.
        for (size_t slice_idx = 0; slice_idx < layer->lslices_ex.size(); ++slice_idx) {
            const LayerSlice          &slice        = layer->lslices_ex.at(slice_idx);
            ObjectPart                &part         = active_object_parts.access(prev_slice_idx_to_object_part_mapping[slice_idx]);
            SliceConnection           &weakest_conn = prev_slice_idx_to_weakest_connection[slice_idx];

            std::vector<Linef> boundary_lines;
            for (const auto &link : slice.overlaps_below) { 
                auto ls = to_unscaled_linesf({layer->lower_layer->lslices[link.slice_idx]});
                boundary_lines.insert(boundary_lines.end(), ls.begin(), ls.end());
            }
            AABBTreeLines::LinesDistancer<Linef> prev_layer_boundary{std::move(boundary_lines)};


            std::vector<ExtrusionLine> current_slice_ext_perims_lines{};
            current_slice_ext_perims_lines.reserve(prev_layer_ext_perim_lines.get_lines().size() / layer->lslices_ex.size());
#ifdef DETAILED_DEBUG_LOGS
            weakest_conn.print_info("weakest connection info: ");
#endif
            // Function that is used when new support point is generated. It will update the ObjectPart stability, weakest conneciton info,
            // and the support presence grid and add the point to the issues.
            auto reckon_new_support_point = [&part, &weakest_conn, &supp_points, &supports_presence_grid, &params,
                                             &layer_idx](const Vec3f &support_point, float force, const Vec2f &dir) {
                if ((supports_presence_grid.position_taken(support_point) && force > 0) || layer_idx <= 1) {
                    return;
                }
                float area = params.support_points_interface_radius * params.support_points_interface_radius * float(PI);
                part.add_support_point(support_point, area);

                float radius = params.support_points_interface_radius;
                supp_points.emplace_back(support_point, force, radius, dir);
                if (force > 0) {
                    supports_presence_grid.take_position(support_point);
                }
                if (weakest_conn.area > EPSILON) { // Do not add it to the weakest connection if it is not valid - does not exist
                    weakest_conn.area += area;
                    weakest_conn.centroid_accumulator += support_point * area;
                    weakest_conn.second_moment_of_area_accumulator += area * support_point.head<2>().cwiseProduct(support_point.head<2>());
                    weakest_conn.second_moment_of_area_covariance_accumulator += area * support_point.x() * support_point.y();
                }
            };

            // first we will check local extrusion stability of bridges, then of perimeters. Perimeters are more important, they
            // account for most of the curling and possible crashes, so on them we will run also global stability check
            for (const auto &island : slice.islands) {
                // Support bridges where needed.
                for (const LayerExtrusionRange &fill_range : island.fills) {
                    const LayerRegion *fill_region = layer->get_region(fill_range.region());
                    for (const auto &fill_idx : fill_range) {
                        const ExtrusionEntity *entity = fill_region->fills().entities[fill_idx];
                        if (entity->role() == ExtrusionRole::BridgeInfill) {
                            for (const ExtrusionLine &bridge :
                                 check_extrusion_entity_stability(entity, fill_region, prev_layer_ext_perim_lines,prev_layer_boundary, params)) {
                                if (bridge.support_point_generated) {
                                    reckon_new_support_point(create_support_point_position(bridge.b), -EPSILON, Vec2f::Zero());
                                }
                            }
                        }
                    }
                }

                const LayerRegion *perimeter_region = layer->get_region(island.perimeters.region());
                for (const auto &perimeter_idx : island.perimeters) {
                    const ExtrusionEntity     *entity = perimeter_region->perimeters().entities[perimeter_idx];
                    std::vector<ExtrusionLine> perims = check_extrusion_entity_stability(entity, perimeter_region,
                                                                                         prev_layer_ext_perim_lines,prev_layer_boundary, params);
                    for (const ExtrusionLine &perim : perims) {
                        if (perim.support_point_generated) {
                            reckon_new_support_point(create_support_point_position(perim.b), -EPSILON, Vec2f::Zero());
                        }
                        if (perim.is_external_perimeter()) { current_slice_ext_perims_lines.push_back(perim); }
                    }
                }
            }

            LD    current_slice_lines_distancer(current_slice_ext_perims_lines);
            float unchecked_dist = params.min_distance_between_support_points + 1.0f;

            for (const ExtrusionLine &line : current_slice_ext_perims_lines) {
                if ((unchecked_dist + line.len < params.min_distance_between_support_points && line.curled_up_height < 0.3f) || line.len < EPSILON) {
                    unchecked_dist += line.len;
                } else {
                    unchecked_dist                = line.len;
                    Vec2f pivot_site_search_point = Vec2f(line.b + (line.b - line.a).normalized() * 300.0f);
                    auto [dist, nidx,
                          nearest_point]          = current_slice_lines_distancer.distance_from_lines_extra<false>(pivot_site_search_point);
                    Vec3f support_point           = create_support_point_position(nearest_point);
                    auto  force                   = part.is_stable_while_extruding(weakest_conn, line, support_point, bottom_z, params);
                    if (force > 0) { reckon_new_support_point(support_point, force, (line.b - line.a).normalized()); }
                }
            }
            current_layer_ext_perims_lines.insert(current_layer_ext_perims_lines.end(), current_slice_ext_perims_lines.begin(),
                                                  current_slice_ext_perims_lines.end());
        } // slice iterations
        prev_layer_ext_perim_lines = LD(current_layer_ext_perims_lines);
    } // layer iterations
    return supp_points;
}

#ifdef DEBUG_FILES
void debug_export(SupportPoints support_points, std::string file_name)
{
    Slic3r::CNumericLocalesSetter locales_setter;
    {
        FILE *fp = boost::nowide::fopen(debug_out_path((file_name + "_supports.obj").c_str()).c_str(), "w");
        if (fp == nullptr) {
            BOOST_LOG_TRIVIAL(error) << "Debug files: Couldn't open " << file_name << " for writing";
            return;
        }

        for (size_t i = 0; i < support_points.size(); ++i) {
            if (support_points[i].force <= 0) {
                fprintf(fp, "v %f %f %f  %f %f %f\n", support_points[i].position(0), support_points[i].position(1),
                        support_points[i].position(2), 0.0, 1.0, 0.0);
            } else {
                fprintf(fp, "v %f %f %f  %f %f %f\n", support_points[i].position(0), support_points[i].position(1),
                        support_points[i].position(2), 1.0, 0.0, 0.0);
            }
        }

        fclose(fp);
    }
}
#endif

// std::vector<size_t> quick_search(const PrintObject *po, const Params &params) {
//     return {};
// }
SupportPoints full_search(const PrintObject *po, const PrintTryCancel& cancel_func, const Params &params)
{
    SupportPoints supp_points = check_stability(po, cancel_func, params);
#ifdef DEBUG_FILES
    debug_export(supp_points, "issues");
#endif

    return supp_points;
}

void estimate_supports_malformations(SupportLayerPtrs &layers, float flow_width, const Params &params)
{
#ifdef DEBUG_FILES
    FILE *debug_file = boost::nowide::fopen(debug_out_path("supports_malformations.obj").c_str(), "w");
    FILE *full_file = boost::nowide::fopen(debug_out_path("supports_full.obj").c_str(), "w");
#endif

    AABBTreeLines::LinesDistancer<ExtrusionLine> prev_layer_lines{};

    for (SupportLayer *l : layers) {
        std::vector<ExtrusionLine> current_layer_lines;

        for (const ExtrusionEntity *extrusion : l->support_fills.flatten().entities) {
            Polyline pl = extrusion->as_polyline();
            Polygon  pol(pl.points);
            pol.make_counter_clockwise();

            auto annotated_points = estimate_points_properties<true, true, false, false>(pol.points, prev_layer_lines, flow_width);

            for (size_t i = 0; i < annotated_points.size(); ++i) {
                ExtendedPoint &curr_point = annotated_points[i];
                float          line_len   = i > 0 ? ((annotated_points[i - 1].position - curr_point.position).norm()) : 0.0f;
                ExtrusionLine  line_out{i > 0 ? annotated_points[i - 1].position.cast<float>() : curr_point.position.cast<float>(),
                                       curr_point.position.cast<float>(), line_len, extrusion};

                const ExtrusionLine nearest_prev_layer_line = prev_layer_lines.get_lines().size() > 0 ?
                                                                  prev_layer_lines.get_line(curr_point.nearest_prev_layer_line) :
                                                                  ExtrusionLine{};

                Vec2f v1 = (nearest_prev_layer_line.b - nearest_prev_layer_line.a);
                Vec2f v2 = (curr_point.position.cast<float>() - nearest_prev_layer_line.a);
                auto  d  = (v1.x() * v2.y()) - (v1.y() * v2.x());
                if (d > 0) {
                    curr_point.distance *= -1.0f;
                }

                line_out.curled_up_height = estimate_curled_up_height(curr_point, l->height, flow_width,
                                                                      nearest_prev_layer_line.curled_up_height, params);

                current_layer_lines.push_back(line_out);
            }
        }

        for (const ExtrusionLine &line : current_layer_lines) {
            if (line.curled_up_height > 0.3f) {
                l->malformed_lines.push_back(Line{Point::new_scale(line.a), Point::new_scale(line.b)});
            }
        }

#ifdef DEBUG_FILES
        for (const ExtrusionLine &line : current_layer_lines) {
            if (line.curled_up_height > 0.3f) {
                Vec3f color = value_to_rgbf(-EPSILON, l->height * params.max_curled_height_factor, line.curled_up_height);
                fprintf(debug_file, "v %f %f %f  %f %f %f\n", line.b[0], line.b[1], l->print_z, color[0], color[1], color[2]);
            }
        }
        for (const ExtrusionLine &line : current_layer_lines) {
            Vec3f color = value_to_rgbf(-EPSILON, l->height * params.max_curled_height_factor, line.curled_up_height);
            fprintf(full_file, "v %f %f %f  %f %f %f\n", line.b[0], line.b[1], l->print_z, color[0], color[1], color[2]);
        }
#endif

        prev_layer_lines = LD{current_layer_lines};
    }

#ifdef DEBUG_FILES
    fclose(debug_file);
    fclose(full_file);
#endif
}

void estimate_malformations(LayerPtrs &layers, const Params &params)
{
#ifdef DEBUG_FILES
    FILE *debug_file = boost::nowide::fopen(debug_out_path("object_malformations.obj").c_str(), "w");
    FILE *full_file = boost::nowide::fopen(debug_out_path("object_full.obj").c_str(), "w");
#endif

    LD prev_layer_lines{};

    for (Layer *l : layers) {
        std::vector<Linef> boundary_lines = l->lower_layer != nullptr ? to_unscaled_linesf(l->lower_layer->lslices) : std::vector<Linef>();
        AABBTreeLines::LinesDistancer<Linef> prev_layer_boundary{std::move(boundary_lines)};
        std::vector<ExtrusionLine>           current_layer_lines;
        for (const LayerRegion *layer_region : l->regions()) {
            for (const ExtrusionEntity *extrusion : layer_region->perimeters().flatten().entities) {
                Points extrusion_pts;
                extrusion->collect_points(extrusion_pts);
                float flow_width       = get_flow_width(layer_region, extrusion->role());
                auto  annotated_points = estimate_points_properties<true, false, false, false>(extrusion_pts, prev_layer_lines, flow_width,
                                                                                              params.bridge_distance);
                for (size_t i = 0; i < annotated_points.size(); ++i) {
                    ExtendedPoint &curr_point = annotated_points[i];
                    float          line_len   = i > 0 ? ((annotated_points[i - 1].position - curr_point.position).norm()) : 0.0f;
                    ExtrusionLine  line_out{i > 0 ? annotated_points[i - 1].position.cast<float>() : curr_point.position.cast<float>(),
                                           curr_point.position.cast<float>(), line_len, extrusion};

                    const ExtrusionLine nearest_prev_layer_line = prev_layer_lines.get_lines().size() > 0 ?
                                                                      prev_layer_lines.get_line(curr_point.nearest_prev_layer_line) :
                                                                      ExtrusionLine{};

                    float sign = (prev_layer_boundary.distance_from_lines<true>(curr_point.position) + 0.5f * flow_width) < 0.0f ? -1.0f :
                                                                                                                                   1.0f;
                    curr_point.distance *= sign;

                    line_out.curled_up_height = estimate_curled_up_height(curr_point, layer_region->layer()->height, flow_width,
                                                                          nearest_prev_layer_line.curled_up_height, params);

                    current_layer_lines.push_back(line_out);
                }
            }
        }

        for (const ExtrusionLine &line : current_layer_lines) {
            if (line.curled_up_height > 0.3f) {
                l->malformed_lines.push_back(Line{Point::new_scale(line.a), Point::new_scale(line.b)});
            }
        }

#ifdef DEBUG_FILES
        for (const ExtrusionLine &line : current_layer_lines) {
            if (line.curled_up_height > 0.3f) {
                Vec3f color = value_to_rgbf(-EPSILON, l->height * params.max_curled_height_factor, line.curled_up_height);
                fprintf(debug_file, "v %f %f %f  %f %f %f\n", line.b[0], line.b[1], l->print_z, color[0], color[1], color[2]);
            }
        }
         for (const ExtrusionLine &line : current_layer_lines) {
                Vec3f color = value_to_rgbf(-EPSILON, l->height * params.max_curled_height_factor, line.curled_up_height);
                fprintf(full_file, "v %f %f %f  %f %f %f\n", line.b[0], line.b[1], l->print_z, color[0], color[1], color[2]);
        }
#endif

        prev_layer_lines = LD{current_layer_lines};
    }

#ifdef DEBUG_FILES
    fclose(debug_file);
    fclose(full_file);
#endif
}

} // namespace SupportSpotsGenerator
} // namespace Slic3r
