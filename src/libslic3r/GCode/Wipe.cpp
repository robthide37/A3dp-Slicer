#include "Wipe.hpp"
#include "../GCode.hpp"

#include <string_view>

#include <Eigen/Geometry>

using namespace std::string_view_literals;

namespace Slic3r::GCode {

void Wipe::init(const PrintConfig &config, const std::vector<unsigned int> &extruders)
{
    this->reset_path();

    // Calculate maximum wipe length to accumulate by the wipe cache.
    // Paths longer than wipe_xy should never be needed for the wipe move.
    double wipe_xy = 0;
    const bool multimaterial = extruders.size() > 1;
    for (auto id : extruders)
        if (config.wipe.get_at(id)) {
            // Wipe length to extrusion ratio.
            const double xy_to_e = this->calc_xy_to_e_ratio(config, id);
            wipe_xy = std::max(wipe_xy, xy_to_e * config.retract_length.get_at(id));
            if (multimaterial)
                wipe_xy = std::max(wipe_xy, xy_to_e * config.retract_length_toolchange.get_at(id));
        }

    if (wipe_xy == 0)
        this->disable();
    else
        this->enable(wipe_xy);
}

void Wipe::set_path(SmoothPath &&path, bool reversed)
{
    this->reset_path();

    if (this->enabled() && ! path.empty()) {
        if (reversed) {
            m_path = std::move(path.back().path);
            Geometry::ArcWelder::reverse(m_path);
            int64_t len = Geometry::ArcWelder::estimate_path_length(m_path);
            for (auto it = std::next(path.rbegin()); len < m_wipe_len_max && it != path.rend(); ++ it) {
                if (it->path_attributes.role.is_bridge())
                    break; // Do not perform a wipe on bridges.
                assert(it->path.size() >= 2);
                assert(m_path.back().point == it->path.back().point);
                if (m_path.back().point != it->path.back().point)
                    // ExtrusionMultiPath is interrupted in some place. This should not really happen.
                    break;
                len += Geometry::ArcWelder::estimate_path_length(it->path);
                m_path.insert(m_path.end(), it->path.rbegin() + 1, it->path.rend());
            }
        } else {
            m_path = std::move(path.front().path);
            int64_t len = Geometry::ArcWelder::estimate_path_length(m_path);
            for (auto it = std::next(path.begin()); len < m_wipe_len_max && it != path.end(); ++ it) {
                if (it->path_attributes.role.is_bridge())
                    break; // Do not perform a wipe on bridges.
                assert(it->path.size() >= 2);
                assert(m_path.back().point == it->path.front().point);
                if (m_path.back().point != it->path.front().point)
                    // ExtrusionMultiPath is interrupted in some place. This should not really happen.
                    break;
                len += Geometry::ArcWelder::estimate_path_length(it->path);
                m_path.insert(m_path.end(), it->path.begin() + 1, it->path.end());
            }
        }
    }

    assert(m_path.empty() || m_path.size() > 1);
}

std::string Wipe::wipe(GCodeGenerator &gcodegen, bool toolchange)
{
    std::string     gcode;
    const Extruder &extruder = *gcodegen.writer().extruder();

    // Remaining quantized retraction length.
    if (double retract_length = extruder.retract_to_go(toolchange ? extruder.retract_length_toolchange() : extruder.retract_length()); 
        retract_length > 0 && this->has_path()) {
        // Delayed emitting of a wipe start tag.
        bool wiped = false;
        const double wipe_speed = this->calc_wipe_speed(gcodegen.writer().config);
        auto start_wipe = [&wiped, &gcode, &gcodegen, wipe_speed](){
            if (! wiped) {
                wiped = true;
                gcode += ";" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Wipe_Start) + "\n";
                gcode += gcodegen.writer().set_speed(wipe_speed * 60, {}, gcodegen.enable_cooling_markers() ? ";_WIPE"sv : ""sv);
            }
        };
        const double xy_to_e    = this->calc_xy_to_e_ratio(gcodegen.writer().config, extruder.id());
        auto         wipe_linear = [&gcode, &gcodegen, &retract_length, xy_to_e](const Vec2d &prev, Vec2d &p) {
            double segment_length = (p - prev).norm();
            // Quantize E axis as it is to be 
            double dE = GCodeFormatter::quantize_e(xy_to_e * segment_length);
            bool   done = false;
            if (dE > retract_length - EPSILON) {
                if (dE > retract_length + EPSILON)
                    // Shorten the segment.
                    p = gcodegen.point_to_gcode_quantized(prev + (p - prev) * (retract_length / dE));
                dE   = retract_length;
                done = true;
            }
            gcode += gcodegen.writer().extrude_to_xy(p, -dE, "wipe and retract"sv);
            retract_length -= dE;
            return done;
        };
        auto         wipe_arc = [&gcode, &gcodegen, &retract_length, xy_to_e](const Vec2d &prev, Vec2d &p, float radius_in, const bool ccw) {
            double radius = GCodeFormatter::quantize_xyzf(radius_in);
            Vec2f  center = Geometry::ArcWelder::arc_center(prev.cast<float>(), p.cast<float>(), float(radius), ccw);
            float  angle  = Geometry::ArcWelder::arc_angle(prev.cast<float>(), p.cast<float>(), float(radius));
            double segment_length = angle * std::abs(radius);
            double dE = GCodeFormatter::quantize_e(xy_to_e * segment_length);
            bool   done = false;
            if (dE > retract_length - EPSILON) {
                if (dE > retract_length + EPSILON)
                    // Shorten the segment.
                    p = gcodegen.point_to_gcode_quantized(Point(prev).rotated((ccw ? angle : -angle) * (retract_length / dE), center.cast<coord_t>()));
                dE   = retract_length;
                done = true;
            }
            gcode += gcodegen.writer().extrude_to_xy_G2G3R(p, radius, ccw, -dE, "wipe and retract"sv);
            retract_length -= dE;
            return done;
        };
        // Start with the current position, which may be different from the wipe path start in case of loop clipping.
        Vec2d prev = gcodegen.point_to_gcode_quantized(gcodegen.last_pos());
        auto  it   = this->path().begin();
        Vec2d p    = gcodegen.point_to_gcode_quantized(it->point + m_offset);
        ++ it;
        bool done = false;
        if (p != prev) {
            start_wipe();
            done = wipe_linear(prev, p);
        }
        if (! done) {
            prev = p;
            auto end = this->path().end();
            for (; it != end && ! done; ++ it) {
                p = gcodegen.point_to_gcode_quantized(it->point + m_offset);
                if (p != prev) {
                    start_wipe();
                    if (it->linear() ?
                        wipe_linear(prev, p) :
                        wipe_arc(prev, p, it->radius, it->ccw()))
                        break;
                    prev = p;
                }
            }
        }
        if (wiped) {
            // add tag for processor
            gcode += ";" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Wipe_End) + "\n";
            gcodegen.set_last_pos(gcodegen.gcode_to_point(p));
        }
    }

    // Prevent wiping again on the same path.
    this->reset_path();
    return gcode;
}

// Make a little move inwards before leaving loop after path was extruded,
// thus the current extruder position is at the end of a path and the path
// may not be closed in case the loop was clipped to hide a seam.
std::optional<Point> wipe_hide_seam(const SmoothPath &path, bool is_hole, double wipe_length)
{
    assert(! path.empty());
    assert(path.front().path.size() >= 2);
    assert(path.back().path.size() >= 2);

    // Heuristics for estimating whether there is a chance that the wipe move will fit inside a small perimeter
    // or that the wipe move direction could be calculated with reasonable accuracy.
    if (longer_than(path, 2.5 * wipe_length)) {
        // The print head will be moved away from path end inside the island.
        Point p_current = path.back().path.back().point;
        Point p_next = path.front().path.front().point;
        Point p_prev;
        {
            // Is the seam hiding gap large enough already?
            double l = wipe_length - (p_next - p_current).cast<double>().norm();
            if (l > 0) {
                // Not yet.
                std::optional<Point> n = sample_path_point_at_distance_from_start(path, l);
                assert(n);
                if (! n)
                    // Wipe move cannot be calculated, the loop is not long enough. This should not happen due to the longer_than() test above.
                    return {};
            }
            if (std::optional<Point> p = sample_path_point_at_distance_from_end(path, wipe_length); p)
                p_prev = *p;
            else
                // Wipe move cannot be calculated, the loop is not long enough. This should not happen due to the longer_than() test above.
                return {};
        }
        // Detect angle between last and first segment.
        // The side depends on the original winding order of the polygon (left for contours, right for holes).
        double angle_inside = angle(p_next - p_current, p_prev - p_current);
        assert(angle_inside >= -M_PI && angle_inside <= M_PI);
        // 3rd of this angle will be taken, thus make the angle monotonic before interpolation.
        if (is_hole) {
            if (angle_inside > 0)
                angle_inside -= 2.0 * M_PI;
        } else {
            if (angle_inside < 0)
                angle_inside += 2.0 * M_PI;
        }
        // Rotate the forward segment inside by 1/3 of the wedge angle.
        auto v_rotated = Eigen::Rotation2D(angle_inside) * (p_next - p_current).cast<double>().normalized();
        return std::make_optional<Point>(p_current + (v_rotated * wipe_length).cast<coord_t>());
    }

    return {};
}

} // namespace Slic3r::GCode
