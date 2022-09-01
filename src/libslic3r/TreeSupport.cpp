// Tree supports by Thomas Rahm, losely based on Tree Supports by CuraEngine.
// Original source of Thomas Rahm's tree supports:
// https://github.com/ThomasRahm/CuraEngine
//
// Original CuraEngine copyright:
// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "TreeSupport.hpp"
#include "BuildVolume.hpp"
#include "ClipperUtils.hpp"
#include "EdgeGrid.hpp"
#include "Fill/Fill.hpp"
#include "Layer.hpp"
#include "Print.hpp"
#include "MultiPoint.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "MutablePolygon.hpp"
#include "SupportMaterial.hpp"

#include <cassert>
#include <chrono>
#include <fstream>
#include <optional>
#include <stdio.h>
#include <string>
#ifdef _WIN32
    #include <windows.h> //todo Remove!  ONLY FOR PUBLIC BETA!!
#endif // _WIN32

#include <boost/log/trivial.hpp>

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

namespace Slic3r
{

enum class LineStatus
{
    INVALID,
    TO_MODEL,
    TO_MODEL_GRACIOUS,
    TO_MODEL_GRACIOUS_SAFE,
    TO_BP,
    TO_BP_SAFE
};

using LineInformation = std::vector<std::pair<Point, LineStatus>>;
using LineInformations = std::vector<LineInformation>;

static inline void validate_range(const Point &pt)
{
    static constexpr const int32_t hi = 65536 * 16384;
    if (pt.x() > hi || pt.y() > hi || -pt.x() > hi || -pt.y() > hi) 
      throw ClipperLib::clipperException("Coordinate outside allowed range");    
}

static inline void validate_range(const Points &points) 
{
    for (const Point &p : points)
        validate_range(p);
}

static inline void validate_range(const MultiPoint &mp) 
{
    validate_range(mp.points);
}

static inline void validate_range(const Polygons &polygons) 
{
    for (const Polygon &p : polygons)
        validate_range(p);
}

static inline void validate_range(const Polylines &polylines) 
{
    for (const Polyline &p : polylines)
        validate_range(p);
}

static inline void validate_range(const LineInformation &lines)
{
    for (const auto& p : lines)
        validate_range(p.first);
}

static inline void validate_range(const LineInformations &lines)
{
    for (const LineInformation &l : lines)
        validate_range(l);
}

static inline void clip_for_diff(const Polygon &src, const BoundingBox &bbox, Polygon &out)
{
    out.clear();
    const size_t cnt = src.points.size();
    if (cnt < 3)
        return;

    enum class Side {
        Left   = 1,
        Right  = 2,
        Top    = 4,
        Bottom = 8
    };

    auto sides = [bbox](const Point &p) {
        return  int(p.x() < bbox.min.x()) * int(Side::Left) +
                int(p.x() > bbox.max.x()) * int(Side::Right) +
                int(p.y() < bbox.min.y()) * int(Side::Bottom) +
                int(p.y() > bbox.max.y()) * int(Side::Top);
    };

    int sides_prev = sides(src.points.back());
    int sides_this = sides(src.points.front());
    const size_t last = cnt - 1;
    for (size_t i = 0; i < last; ++ i) {
        int sides_next = sides(src.points[i + 1]);
        if (// This point is inside. Take it.
            sides_this == 0 ||
            // Either this point is outside and previous or next is inside, or
            // the edge possibly cuts corner of the bounding box.
            (sides_prev & sides_this & sides_next) == 0) {
            out.points.emplace_back(src.points[i]);
            sides_prev = sides_this;
        } else {
            // All the three points (this, prev, next) are outside at the same side.
            // Ignore this point.
        }
        sides_this = sides_next;
    }
    // For the last point, if src is completely outside bbox, then out.points will be empty. Just use the first point instead.
    int sides_next = sides(out.points.empty() ? src.points.front() : out.points.front());
    if (// The last point is inside. Take it.
        sides_this == 0 ||
        // Either this point is outside and previous or next is inside, or
        // the edge possibly cuts corner of the bounding box.
        (sides_prev & sides_this & sides_next) == 0)
        out.points.emplace_back(src.points.back());
}

[[nodiscard]] static inline Polygon clip_for_diff(const Polygon &src, const BoundingBox &bbox)
{
    Polygon out;
    clip_for_diff(src, bbox, out);
    return out;
}

[[nodiscard]] static inline Polygons clip_for_diff(const Polygons &src, const BoundingBox &bbox)
{
    Polygons out;
    out.reserve(src.size());
    for (const Polygon &p : src)
        out.emplace_back(clip_for_diff(p, bbox));
    return out;
}

[[nodiscard]] static inline Polygons diff_clipped(const Polygons &src, const Polygons &clipping)
{
    return diff(src, clip_for_diff(clipping, get_extents(src).inflated(SCALED_EPSILON)));
}

static constexpr const auto tiny_area_threshold = sqr(scaled<double>(0.001));

static std::vector<std::pair<TreeSupport::TreeSupportSettings, std::vector<size_t>>> group_meshes(const Print &print, const std::vector<size_t> &print_object_ids)
{
    std::vector<std::pair<TreeSupport::TreeSupportSettings, std::vector<size_t>>> grouped_meshes;

    //FIXME this is ugly, it does not belong here.
    for (size_t object_id = 0; object_id < print_object_ids.size(); ++ object_id) {
        const PrintObject       &print_object  = *print.get_object(object_id);
        const PrintObjectConfig &object_config = print_object.config();
        if (object_config.support_material_interface_layers >= 2)
            TreeSupport::TreeSupportSettings::some_model_contains_thick_roof = true;
        if (object_config.support_material_contact_distance < EPSILON)
            // || min_feature_size < scaled<coord_t>(0.1) that is the minimum line width
            TreeSupport::TreeSupportSettings::has_to_rely_on_min_xy_dist_only = true;
    }

    size_t largest_printed_mesh_idx = 0;

    // Group all meshes that can be processed together. NOTE this is different from mesh-groups! Only one setting object is needed per group, 
    // as different settings in the same group may only occur in the tip, which uses the original settings objects from the meshes.
    for (size_t object_id = 0; object_id < print_object_ids.size(); ++ object_id) {
        const PrintObject       &print_object  = *print.get_object(object_id);
#ifndef NDEBUG
        const PrintObjectConfig &object_config = print_object.config();
#endif // NDEBUG
        // Support must be enabled and set to Tree style.
        assert(object_config.support_material);
        assert(object_config.support_material_style == smsTree);

        bool found_existing_group = false;
        TreeSupport::TreeSupportSettings next_settings{ TreeSupportMeshGroupSettings{ print_object } };
        //FIXME for now only a single object per group is enabled.
#if 0
        for (size_t idx = 0; idx < grouped_meshes.size(); ++ idx)
            if (next_settings == grouped_meshes[idx].first) {
                found_existing_group = true;
                grouped_meshes[idx].second.emplace_back(object_id);
                // handle some settings that are only used for performance reasons. This ensures that a horrible set setting intended to improve performance can not reduce it drastically.
                grouped_meshes[idx].first.performance_interface_skip_layers = std::min(grouped_meshes[idx].first.performance_interface_skip_layers, next_settings.performance_interface_skip_layers);
            }
#endif
        if (! found_existing_group)
            grouped_meshes.emplace_back(next_settings, std::vector<size_t>{ object_id });

        // no need to do this per mesh group as adaptive layers and raft setting are not setable per mesh.
        if (print.get_object(largest_printed_mesh_idx)->layers().back()->print_z < print_object.layers().back()->print_z)
            largest_printed_mesh_idx = object_id;
    }

#if 0
    {
        std::vector<coord_t> known_z(storage.meshes[largest_printed_mesh_idx].layers.size());
        for (size_t z = 0; z < storage.meshes[largest_printed_mesh_idx].layers.size(); z++)
            known_z[z] = storage.meshes[largest_printed_mesh_idx].layers[z].printZ;
        for (size_t idx = 0; idx < grouped_meshes.size(); ++ idx)
            grouped_meshes[idx].first.setActualZ(known_z);
    }
#endif

    return grouped_meshes;
}

#if 0
// todo remove as only for debugging relevant
[[nodiscard]] static std::string getPolygonAsString(const Polygons& poly)
{
    std::string ret;
    for (auto path : poly)
        for (Point p : path) {
            if (ret != "")
                ret += ",";
            ret += "(" + std::to_string(p.x()) + "," + std::to_string(p.y()) + ")";
        }
    return ret;
}
#endif

//todo Remove! Only relevant for public BETA!
static bool inline g_showed_critical_error = false;
static bool inline g_showed_performance_warning = false;
void TreeSupport::showError(std::string message, bool critical)
{ // todo Remove!  ONLY FOR PUBLIC BETA!!

#if defined(_WIN32) && defined(TREE_SUPPORT_SHOW_ERRORS)
    auto bugtype = std::string(critical ? " This is a critical bug. It may cause missing or malformed branches.\n" : "This bug should only decrease performance.\n");
    bool show    = (critical && !g_showed_critical_error) || (!critical && !g_showed_performance_warning);
    (critical ? g_showed_critical_error : g_showed_performance_warning) = true;

    if (show)
        MessageBoxA(nullptr, std::string("TreeSupport_2 MOD detected an error while generating the tree support.\nPlease report this back to me with profile and model.\nRevision 5.0\n" + message + "\n" + bugtype).c_str(), 
            "Bug detected!", MB_OK | MB_SYSTEMMODAL | MB_SETFOREGROUND | MB_ICONWARNING);
#endif // WIN32
}

[[nodiscard]] static const std::vector<Polygons> generate_overhangs(const PrintObject &print_object)
{
    std::vector<Polygons> out(print_object.layer_count(), Polygons{});

    const PrintObjectConfig &config                 = print_object.config();
    const bool               support_auto           = config.support_material_auto.value;
    const int                support_enforce_layers = config.support_material_enforce_layers.value;
    std::vector<Polygons>    enforcers_layers{ print_object.slice_support_enforcers() };
    std::vector<Polygons>    blockers_layers{ print_object.slice_support_blockers() };
    print_object.project_and_append_custom_facets(false, EnforcerBlockerType::ENFORCER, enforcers_layers);
    print_object.project_and_append_custom_facets(false, EnforcerBlockerType::BLOCKER, blockers_layers);
    const int                support_threshold      = config.support_material_threshold.value;
    const bool               support_threshold_auto = support_threshold == 0;
    // +1 makes the threshold inclusive
    double                   tan_threshold          = support_threshold_auto ? 0. : tan(M_PI * double(support_threshold + 1) / 180.);

    tbb::parallel_for(tbb::blocked_range<LayerIndex>(1, out.size()),
        [&print_object, &enforcers_layers, &blockers_layers, support_auto, support_enforce_layers, support_threshold_auto, tan_threshold, &out]
        (const tbb::blocked_range<LayerIndex> &range) {
        for (LayerIndex layer_id = range.begin(); layer_id < range.end(); ++ layer_id) {
            const Layer   &current_layer  = *print_object.get_layer(layer_id);
            const Layer   &lower_layer    = *print_object.get_layer(layer_id - 1);
            // Full overhangs with zero lower_layer_offset and no blockers applied.
            Polygons       raw_overhangs;
            bool           raw_overhangs_calculated = false;
            // Final overhangs.
            Polygons       overhangs;
            // For how many layers full overhangs shall be supported.
            const bool     enforced_layer = layer_id < support_enforce_layers;
            if (support_auto || enforced_layer) {
                float lower_layer_offset;
                if (enforced_layer)
                    lower_layer_offset = 0;
                else if (support_threshold_auto) {
                    float external_perimeter_width = 0;
                    for (const LayerRegion *layerm : lower_layer.regions())
                        external_perimeter_width += layerm->flow(frExternalPerimeter).scaled_width();
                    external_perimeter_width /= lower_layer.region_count();
                    lower_layer_offset = float(0.5 * external_perimeter_width);
                } else
                    lower_layer_offset = scaled<float>(lower_layer.height / tan_threshold);
                overhangs = lower_layer_offset == 0 ?
                    diff(current_layer.lslices, lower_layer.lslices) :
                    diff(current_layer.lslices, offset(lower_layer.lslices, lower_layer_offset));
                if (lower_layer_offset == 0) {
                    raw_overhangs = overhangs;
                    raw_overhangs_calculated = true;
                }
                if (! (enforced_layer || blockers_layers.empty() || blockers_layers[layer_id].empty()))
                    overhangs = diff(overhangs, blockers_layers[layer_id]);
            }
            if (! enforcers_layers.empty() && ! enforcers_layers[layer_id].empty())
                // Has some support enforcers at this layer, apply them to the overhangs, don't apply the support threshold angle.
                if (Polygons enforced_overhangs = intersection(raw_overhangs_calculated ? raw_overhangs : diff(current_layer.lslices, lower_layer.lslices), enforcers_layers[layer_id]);
                    ! enforced_overhangs.empty()) {
                    //FIXME this is a hack to make enforcers work on steep overhangs.
                    enforced_overhangs = diff(offset(enforced_overhangs, 
                        //FIXME this is a fudge constant!
                        scaled<float>(0.4)), 
                        lower_layer.lslices);
                    overhangs = overhangs.empty() ? std::move(enforced_overhangs) : union_(overhangs, enforced_overhangs);
                }
            out[layer_id] = std::move(overhangs);
        }
    });

    return out;
}

/*!
 * \brief Precalculates all avoidances, that could be required.
 *
 * \param storage[in] Background storage to access meshes.
 * \param currently_processing_meshes[in] Indexes of all meshes that are processed in this iteration
 */
[[nodiscard]] static LayerIndex precalculate(const Print &print, const std::vector<Polygons> &overhangs, const TreeSupport::TreeSupportSettings &config, const std::vector<size_t> &object_ids, TreeModelVolumes &volumes)
{
    // calculate top most layer that is relevant for support
    LayerIndex max_layer = 0;
    for (size_t object_id : object_ids) {
        const PrintObject &print_object         = *print.get_object(object_id);
        int                max_support_layer_id = 0;
        for (int layer_id = 1; layer_id < int(print_object.layer_count()); ++ layer_id)
            if (! overhangs[layer_id].empty())
                max_support_layer_id = layer_id;
        max_layer = std::max(max_support_layer_id - int(config.z_distance_top_layers), 0);
    }
    if (max_layer > 0)
        // The actual precalculation happens in TreeModelVolumes.
        volumes.precalculate(max_layer);
    return max_layer;
}

//FIXME this is an ugly wrapper interface for a single print object and a phony build volume.
void TreeSupport::generateSupportAreas(PrintObject& print_object)
{
    size_t idx = 0;
    for (PrintObject *po : print_object.print()->objects()) {
        if (po == &print_object)
            break;
        ++ idx;
    }
    this->generateSupportAreas(*print_object.print(), BuildVolume(Pointfs{ Vec2d{ -300., -300. }, Vec2d{ -300., +300. }, Vec2d{ +300., +300. }, Vec2d{ +300., -300. } }, 0.), { idx });
}

void TreeSupport::generateSupportAreas(Print &print, const BuildVolume &build_volume, const std::vector<size_t> &print_object_ids)
{
    g_showed_critical_error = false;
    g_showed_performance_warning = false;

    std::vector<std::pair<TreeSupport::TreeSupportSettings, std::vector<size_t>>> grouped_meshes = group_meshes(print, print_object_ids);
    if (grouped_meshes.empty())
        return;

    size_t counter = 0;

    // Process every mesh group. These groups can not be processed parallel, as the processing in each group is parallelized, and nested parallelization is disables and slow.
    for (std::pair<TreeSupportSettings, std::vector<size_t>> &processing : grouped_meshes)
    {
        // process each combination of meshes
        // this struct is used to easy retrieve setting. No other function except those in TreeModelVolumes and generateInitialAreas have knowledge of the existence of multiple meshes being processed.
        //FIXME this is a copy
        m_config = processing.first;
        BOOST_LOG_TRIVIAL(info) << "Processing support tree mesh group " << counter + 1 << " of " << grouped_meshes.size() << " containing " << grouped_meshes[counter].second.size() << " meshes.";
        auto t_start = std::chrono::high_resolution_clock::now();
#if 0
        std::vector<Polygons> exclude(num_support_layers);
        // get all already existing support areas and exclude them
        tbb::parallel_for(tbb::blocked_range<size_t>(0, num_support_layers),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx) {
                Polygons exlude_at_layer;
                append(exlude_at_layer, storage.support.supportLayers[layer_idx].support_bottom);
                append(exlude_at_layer, storage.support.supportLayers[layer_idx].support_roof);
                for (auto part : storage.support.supportLayers[layer_idx].support_infill_parts)
                    append(exlude_at_layer, part.outline);
                exclude[layer_idx] = union_(exlude_at_layer);
            }
        });
#endif
#ifdef SLIC3R_TREESUPPORTS_PROGRESS
        m_progress_multiplier = 1.0 / double(m_grouped_meshes.size());
        m_progress_offset = counter == 0 ? 0 : TREE_PROGRESS_TOTAL * (double(counter) * m_progress_multiplier);
#endif // SLIC3R_TREESUPPORT_PROGRESS
        PrintObject &print_object = *print.get_object(processing.second.front());
        m_volumes = TreeModelVolumes(print_object, build_volume, m_config.maximum_move_distance, m_config.maximum_move_distance_slow, processing.second.front(), m_progress_multiplier, m_progress_offset, /* additional_excluded_areas */{});

        //FIXME generating overhangs just for the furst mesh of the group.
        assert(processing.second.size() == 1);
        std::vector<Polygons>        overhangs = generate_overhangs(*print.get_object(processing.second.front()));

        // ### Precalculate avoidances, collision etc.
        size_t num_support_layers = precalculate(print, overhangs, processing.first, processing.second, m_volumes);
        if (num_support_layers == 0)
            continue;

        auto t_precalc = std::chrono::high_resolution_clock::now();

        // value is the area where support may be placed. As this is calculated in CreateLayerPathing it is saved and reused in drawAreas
        std::vector<std::set<SupportElement*>> move_bounds(num_support_layers);

        // ### Place tips of the support tree
        SupportGeneratorLayersPtr    bottom_contacts(num_support_layers, nullptr);
        SupportGeneratorLayersPtr    top_contacts(num_support_layers, nullptr);
        SupportGeneratorLayersPtr    top_interface_layers(num_support_layers, nullptr);
        SupportGeneratorLayersPtr    intermediate_layers(num_support_layers, nullptr);
        SupportGeneratorLayerStorage layer_storage;

        for (size_t mesh_idx : processing.second)
            generateInitialAreas(*print.get_object(mesh_idx), overhangs, move_bounds, top_contacts, top_interface_layers, layer_storage);
        auto t_gen = std::chrono::high_resolution_clock::now();

        // ### Propagate the influence areas downwards.
        createLayerPathing(move_bounds);
        auto t_path = std::chrono::high_resolution_clock::now();

        // ### Set a point in each influence area
        createNodesFromArea(move_bounds);
        auto t_place = std::chrono::high_resolution_clock::now();

        // ### draw these points as circles
        drawAreas(*print.get_object(processing.second.front()), overhangs, move_bounds, 
            bottom_contacts, top_contacts, intermediate_layers, layer_storage);

        auto t_draw = std::chrono::high_resolution_clock::now();
        auto dur_pre_gen = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_precalc - t_start).count();
        auto dur_gen = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_gen - t_precalc).count();
        auto dur_path = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_path - t_gen).count();
        auto dur_place = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_place - t_path).count();
        auto dur_draw = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_draw - t_place).count();
        auto dur_total = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_draw - t_start).count();
        BOOST_LOG_TRIVIAL(info) <<
            "Total time used creating Tree support for the currently grouped meshes: " << dur_total << " ms. "
            "Different subtasks:\nCalculating Avoidance: " << dur_pre_gen << " ms "
            "Creating inital influence areas: " << dur_gen << " ms "
            "Influence area creation: " << dur_path << "ms "
            "Placement of Points in InfluenceAreas: " << dur_place << "ms "
            "Drawing result as support " << dur_draw << " ms";
//        if (m_config.branch_radius==2121)
//            BOOST_LOG_TRIVIAL(error) << "Why ask questions when you already know the answer twice.\n (This is not a real bug, please dont report it.)";
        
        for (auto &layer : move_bounds) {
            for (auto elem : layer) {
                delete elem->area;
                delete elem;
            }
        }

        auto remove_undefined_layers = [](SupportGeneratorLayersPtr &layers) {
            layers.erase(std::remove_if(layers.begin(), layers.end(), [](const SupportGeneratorLayer* ptr) { return ptr == nullptr; }), layers.end());
        };
        remove_undefined_layers(bottom_contacts);
        remove_undefined_layers(top_contacts);
        remove_undefined_layers(intermediate_layers);

        // Produce the support G-code.
        // Used by both classic and tree supports.
        SupportParameters support_params(print_object);
        support_params.with_sheath = true;
        support_params.support_density = 0;
        SupportGeneratorLayersPtr interface_layers, base_interface_layers;
        SupportGeneratorLayersPtr raft_layers = generate_raft_base(print_object, support_params, print_object.slicing_parameters(), top_contacts, interface_layers, base_interface_layers, intermediate_layers, layer_storage);
#if 1 //#ifdef SLIC3R_DEBUG
        SupportGeneratorLayersPtr layers_sorted =
#endif // SLIC3R_DEBUG
            generate_support_layers(print_object, raft_layers, bottom_contacts, top_contacts, intermediate_layers, interface_layers, base_interface_layers);
        // Don't fill in the tree supports, make them hollow with just a single sheath line.
        generate_support_toolpaths(print_object.support_layers(), print_object.config(), support_params, print_object.slicing_parameters(),
            raft_layers, bottom_contacts, top_contacts, intermediate_layers, interface_layers, base_interface_layers);

 #if 0
//#ifdef SLIC3R_DEBUG
        {
            static int iRun = 0;
            ++ iRun;
            size_t layer_id = 0;
            for (int i = 0; i < int(layers_sorted.size());) {
                // Find the last layer with roughly the same print_z, find the minimum layer height of all.
                // Due to the floating point inaccuracies, the print_z may not be the same even if in theory they should.
                int j = i + 1;
                coordf_t zmax = layers_sorted[i]->print_z + EPSILON;
                bool empty = layers_sorted[i]->polygons.empty();
                for (; j < layers_sorted.size() && layers_sorted[j]->print_z <= zmax; ++j)
                    if (!layers_sorted[j]->polygons.empty())
                        empty = false;
                if (!empty) {
                    export_print_z_polygons_to_svg(
                        debug_out_path("support-%d-%lf.svg", iRun, layers_sorted[i]->print_z).c_str(),
                        layers_sorted.data() + i, j - i);
                    export_print_z_polygons_and_extrusions_to_svg(
                        debug_out_path("support-w-fills-%d-%lf.svg", iRun, layers_sorted[i]->print_z).c_str(),
                        layers_sorted.data() + i, j - i,
                        *print_object.support_layers()[layer_id]);
                    ++layer_id;
                }
                i = j;
            }
        }
#endif /* SLIC3R_DEBUG */

        ++ counter;
    }

//   storage.support.generated = true;
}

/*!
 * \brief Converts a Polygons object representing a line into the internal format.
 *
 * \param polylines[in] The Polyline that will be converted.
 * \param layer_idx[in] The current layer.
 * \return All lines of the \p polylines object, with information for each point regarding in which avoidance it is currently valid in.
 */
[[nodiscard]] static LineInformations convertLinesToInternal(
    const TreeModelVolumes &volumes, const TreeSupport::TreeSupportSettings &config,
    const Polylines &polylines, LayerIndex layer_idx)
{
    const bool xy_overrides_z = config.support_xy_overrides_z;

    LineInformations result;
    // Also checks if the position is valid, if it is NOT, it deletes that point
    for (const Polyline &line : polylines) {
        LineInformation res_line;
        for (Point p : line) {
            if (! contains(volumes.getAvoidance(config.getRadius(0), layer_idx, TreeModelVolumes::AvoidanceType::FastSafe, false, !xy_overrides_z), p))
                res_line.emplace_back(p, LineStatus::TO_BP_SAFE);
            else if (! contains(volumes.getAvoidance(config.getRadius(0), layer_idx, TreeModelVolumes::AvoidanceType::Fast, false, !xy_overrides_z), p))
                res_line.emplace_back(p, LineStatus::TO_BP);
            else if (config.support_rests_on_model && ! contains(volumes.getAvoidance(config.getRadius(0), layer_idx, TreeModelVolumes::AvoidanceType::FastSafe, true, !xy_overrides_z), p))
                res_line.emplace_back(p, LineStatus::TO_MODEL_GRACIOUS_SAFE);
            else if (config.support_rests_on_model && ! contains(volumes.getAvoidance(config.getRadius(0), layer_idx, TreeModelVolumes::AvoidanceType::Fast, true, !xy_overrides_z), p))
                res_line.emplace_back(p, LineStatus::TO_MODEL_GRACIOUS);
            else if (config.support_rests_on_model && ! contains(volumes.getCollision(config.getRadius(0), layer_idx, !xy_overrides_z), p))
                res_line.emplace_back(p, LineStatus::TO_MODEL);
            else if (!res_line.empty()) {
                result.emplace_back(res_line);
                res_line.clear();
            }
        }
        if (!res_line.empty()) {
            result.emplace_back(res_line);
            res_line.clear();
        }
    }

    validate_range(result);
    return result;
}

/*!
 * \brief Converts lines in internal format into a Polygons object representing these lines.
 *
 * \param lines[in] The lines that will be converted.
 * \return All lines of the \p lines object as a Polygons object.
 */
[[nodiscard]] static Polylines convertInternalToLines(LineInformations lines)
{
    Polylines result;
    for (LineInformation line : lines) {
        Polyline path;
        for (auto point_data : line)
            path.points.emplace_back(point_data.first);
        result.emplace_back(std::move(path));
    }
    validate_range(result);
    return result;
}

/*!
 * \brief Evaluates if a point has to be added now. Required for a splitLines call in generateInitialAreas.
 *
 * \param current_layer[in] The layer on which the point lies, point and its status.
 * \return whether the point is valid.
 */
[[nodiscard]] static bool evaluatePointForNextLayerFunction(
    const TreeModelVolumes &volumes, const TreeSupport::TreeSupportSettings &config,
    size_t current_layer, std::pair<Point, LineStatus> &p)
{
    using AvoidanceType = TreeSupport::AvoidanceType;
    if (! contains(volumes.getAvoidance(config.getRadius(0), current_layer - 1, p.second == LineStatus::TO_BP_SAFE ? AvoidanceType::FastSafe : AvoidanceType::Fast, false, !config.support_xy_overrides_z), p.first))
        return true;
    if (config.support_rests_on_model && (p.second != LineStatus::TO_BP && p.second != LineStatus::TO_BP_SAFE))
        return ! contains(
            p.second == LineStatus::TO_MODEL_GRACIOUS || p.second == LineStatus::TO_MODEL_GRACIOUS_SAFE ? 
                volumes.getAvoidance(config.getRadius(0), current_layer - 1, p.second == LineStatus::TO_MODEL_GRACIOUS_SAFE ? AvoidanceType::FastSafe : AvoidanceType::Fast, true, !config.support_xy_overrides_z) :
                volumes.getCollision(config.getRadius(0), current_layer - 1, !config.support_xy_overrides_z),
            p.first);
    return false;
}

/*!
 * \brief Evaluates which points of some lines are not valid one layer below and which are. Assumes all points are valid on the current layer. Validity is evaluated using supplied lambda.
 *
 * \param lines[in] The lines that have to be evaluated.
 * \param evaluatePoint[in] The function used to evaluate the points.
 * \return A pair with which points are still valid in the first slot and which are not in the second slot.
 */
template<typename EvaluatePointFn>
[[nodiscard]] static std::pair<LineInformations, LineInformations> splitLines(LineInformations lines, EvaluatePointFn evaluatePoint)
{
    // assumes all Points on the current line are valid

    LineInformations keep(1);
    LineInformations set_free(1);
    enum STATE
    {
        keeping,
        freeing
    };
    for (std::vector<std::pair<Point, LineStatus>> line : lines) {
        STATE current = keeping;
        LineInformation resulting_line;
        for (std::pair<Point, LineStatus> me : line) {
            if (evaluatePoint(me)) {
                if (keeping != current) {
                    if (!resulting_line.empty()) {
                        set_free.emplace_back(resulting_line);
                        resulting_line.clear();
                    }
                    current = keeping;
                }
                resulting_line.emplace_back(me);
            } else {
                if (freeing != current) {
                    if (!resulting_line.empty()) {
                        keep.emplace_back(resulting_line);
                        resulting_line.clear();
                    }
                    current = freeing;
                }
                resulting_line.emplace_back(me);
            }
        }
        if (!resulting_line.empty()) {
            if (current == keeping)
                keep.emplace_back(resulting_line);
            else
                set_free.emplace_back(resulting_line);
        }
    }
    validate_range(keep);
    validate_range(set_free);
    return std::pair<std::vector<std::vector<std::pair<Point, LineStatus>>>, std::vector<std::vector<std::pair<Point, LineStatus>>>>(keep, set_free);
}

// Ported from CURA's PolygonUtils::getNextPointWithDistance()
// Sample a next point at distance "dist" from start_pt on polyline segment (start_idx, start_idx + 1).
// Returns sample point and start index of its segment on polyline if such sample exists.
static std::optional<std::pair<Point, size_t>> polyline_sample_next_point_at_distance(const Points &polyline, const Point &start_pt, size_t start_idx, double dist)
{
    const double                dist2  = sqr(dist);
    const auto                  dist2i = int64_t(dist2);
    static constexpr const auto eps    = scaled<double>(0.01);

    for (size_t i = start_idx + 1; i < polyline.size(); ++ i) {
        const Point p1 = polyline[i];
        if ((p1 - start_pt).cast<int64_t>().squaredNorm() >= dist2i) {
            // The end point is outside the circle with center "start_pt" and radius "dist".
            const Point p0  = polyline[i - 1];
            Vec2d       v   = (p1 - p0).cast<double>();
            double      l2v = v.squaredNorm();
            if (l2v < sqr(eps)) {
                // Very short segment.
                Point c = (p0 + p1) / 2;
                if (std::abs((start_pt - c).cast<double>().norm() - dist) < eps)
                    return std::pair<Point, size_t>{ c, i - 1 };
                else
                    continue;
            }
            Vec2d p0f = (start_pt - p0).cast<double>();
            // Foot point of start_pt into v.
            Vec2d foot_pt = v * (p0f.dot(v) / l2v);
            // Vector from foot point of "start_pt" to "start_pt".
            Vec2d xf = p0f - foot_pt;
            // Squared distance of "start_pt" from the ray (p0, p1).
            double l2_from_line = xf.squaredNorm();
            double det = dist2 - l2_from_line;

            if (det > - SCALED_EPSILON) {
                // The ray (p0, p1) touches or intersects a circle centered at "start_pt" with radius "dist".
                // Distance of the circle intersection point from the foot point.
                double dist_circle_intersection = std::sqrt(std::max(0., det));
                if ((v - foot_pt).cast<double>().norm() > dist_circle_intersection) {
                    // Intersection of the circle with the segment (p0, p1) is on the right side (close to p1) from the foot point.
                    Point p = p0 + (foot_pt + v * (dist_circle_intersection / sqrt(l2v))).cast<coord_t>();
                    validate_range(p);
                    return std::pair<Point, size_t>{ p, i - 1 };
                }
            }
        }
    }
    return {};
}

/*!
 * \brief Eensures that every line segment is about distance in length. The resulting lines may differ from the original but all points are on the original
 *
 * \param input[in] The lines on which evenly spaced points should be placed.
 * \param distance[in] The distance the points should be from each other.
 * \param min_points[in] The amount of points that have to be placed. If not enough can be placed the distance will be reduced to place this many points.
 * \return A Polygons object containing the evenly spaced points. Does not represent an area, more a collection of points on lines.
 */
[[nodiscard]] static Polylines ensureMaximumDistancePolyline(const Polylines &input, double distance, size_t min_points)
{
    Polylines result;
    for (Polyline part : input) {
        if (part.empty())
            continue;

        double len = length(part.points);
        Polyline line;
        double current_distance = std::max(distance, scaled<double>(0.1));
        if (len < 2 * distance && min_points <= 1)
        {
            // Insert the opposite point of the first one.
            //FIXME pretty expensive
            Polyline pl(part);
            pl.clip_end(len / 2);
            line.points.emplace_back(pl.points.back());
        }
        else
        {
            size_t optimal_end_index = part.size() - 1;

            if (part.front() == part.back()) {
                size_t optimal_start_index = 0;
                // If the polyline was a polygon, there is a high chance it was an overhang. Overhangs that are <60� tend to be very thin areas, so lets get the beginning and end of them and ensure that they are supported.
                // The first point of the line will always be supported, so rotate the order of points in this polyline that one of the two corresponding points that are furthest from each other is in the beginning.
                // The other will be manually added (optimal_end_index)
                coord_t max_dist2_between_vertecies = 0;
                for (size_t idx = 0; idx < part.size() - 1; ++ idx) {
                    for (size_t inner_idx = 0; inner_idx < part.size() - 1; inner_idx++) {
                        if ((part[idx] - part[inner_idx]).cast<double>().squaredNorm() > max_dist2_between_vertecies) {
                            optimal_start_index = idx;
                            optimal_end_index = inner_idx;
                            max_dist2_between_vertecies = (part[idx] - part[inner_idx]).cast<double>().squaredNorm();
                        }
                    }
                }
                std::rotate(part.begin(), part.begin() + optimal_start_index, part.end() - 1);
                part[part.size() - 1] = part[0]; // restore that property that this polyline ends where it started.
                optimal_end_index = (part.size() + optimal_end_index - optimal_start_index - 1) % (part.size() - 1);
            }

            while (line.size() < min_points && current_distance >= scaled<double>(0.1))
            {
                line.clear();
                Point current_point = part[0];
                line.points.emplace_back(part[0]);
                if (min_points > 1 || (part[0] - part[optimal_end_index]).cast<double>().norm() > current_distance)
                    line.points.emplace_back(part[optimal_end_index]);
                size_t current_index = 0;
                std::optional<std::pair<Point, size_t>> next_point;
                double next_distance = current_distance;
                // Get points so that at least min_points are added and they each are current_distance away from each other. If that is impossible, decrease current_distance a bit.
                // The input are lines, that means that the line from the last to the first vertex does not have to exist, so exclude all points that are on this line!
                while ((next_point = polyline_sample_next_point_at_distance(part.points, current_point, current_index, next_distance))) {
                    // Not every point that is distance away, is valid, as it may be much closer to another point. This is especially the case when the overhang is very thin.
                    // So this ensures that the points are actually a certain distance from each other.
                    // This assurance is only made on a per polygon basis, as different but close polygon may not be able to use support below the other polygon.
                    double min_distance_to_existing_point = std::numeric_limits<double>::max();
                    for (Point p : line)
                        min_distance_to_existing_point = std::min(min_distance_to_existing_point, (p - next_point->first).cast<double>().norm());
                    if (min_distance_to_existing_point >= current_distance) {
                        // viable point was found. Add to possible result.
                        line.points.emplace_back(next_point->first);
                        current_point = next_point->first;
                        current_index = next_point->second;
                        next_distance = current_distance;
                    } else {
                        if (current_point == next_point->first) {
                            // In case a fixpoint is encountered, better aggressively overcompensate so the code does not become stuck here...
                            BOOST_LOG_TRIVIAL(warning) << "Tree Support: Encountered a fixpoint in polyline_sample_next_point_at_distance. This is expected to happen if the distance (currently " << next_distance << 
                                ") is smaller than 100";
                            TreeSupport::showError("Encountered issue while placing tips. Some tips may be missing.", true);
                            if (next_distance > 2 * current_distance)
                                // This case should never happen, but better safe than sorry.
                                break;
                            next_distance += current_distance;
                            continue;
                        }
                        // if the point was too close, the next possible viable point is at least distance-min_distance_to_existing_point away from the one that was just checked.
                        next_distance = std::max(current_distance - min_distance_to_existing_point, scaled<double>(0.1));
                        current_point = next_point->first;
                        current_index = next_point->second;
                    }
                }
                current_distance *= 0.9;
            }
        }
        result.emplace_back(std::move(line));
    }
    validate_range(result);
    return result;
}

/*!
 * \brief Returns Polylines representing the (infill) lines that will result in slicing the given area
 *
 * \param area[in] The area that has to be filled with infill.
 * \param roof[in] Whether the roofing or regular support settings should be used.
 * \param layer_idx[in] The current layer index.
 * \param support_infill_distance[in] The distance that should be between the infill lines.
 *
 * \return A Polygons object that represents the resulting infill lines.
 */
[[nodiscard]] static Polylines generateSupportInfillLines(
    const Polygons &polygon, const SupportParameters &support_params,
    bool roof, LayerIndex layer_idx, coord_t support_infill_distance)
{
#if 0
    Polygons gaps;
    // as we effectivly use lines to place our supportPoints we may use the Infill class for it, while not made for it it works perfect

    const EFillMethod pattern = roof ? config.roof_pattern : config.support_pattern;

//    const bool zig_zaggify_infill = roof ? pattern == EFillMethod::ZIG_ZAG : config.zig_zaggify_support;
    const bool connect_polygons = false;
    constexpr coord_t support_roof_overlap = 0;
    constexpr size_t infill_multiplier = 1;
    constexpr coord_t outline_offset = 0;
    const int support_shift = roof ? 0 : support_infill_distance / 2;
    const size_t wall_line_count = include_walls && !roof ? config.support_wall_count : 0;
    const Point infill_origin;
    constexpr Polygons* perimeter_gaps = nullptr;
    constexpr bool use_endpieces = true;
    const bool connected_zigzags = roof ? false : config.connect_zigzags;
    const size_t zag_skip_count = roof ? 0 : config.zag_skip_count;
    constexpr coord_t pocket_size = 0;
    std::vector<AngleRadians> angles = roof ? config.support_roof_angles : config.support_infill_angles;
    std::vector<VariableWidthLines> toolpaths;

    const coord_t z = config.getActualZ(layer_idx);
    int divisor = static_cast<int>(angles.size());
    int index = ((layer_idx % divisor) + divisor) % divisor;
    const AngleRadians fill_angle = angles[index];
    Infill roof_computation(pattern, true /* zig_zaggify_infill */, connect_polygons, polygon, 
        roof ? config.support_roof_line_width : config.support_line_width, support_infill_distance, support_roof_overlap, infill_multiplier, 
        fill_angle, z, support_shift, config.resolution, wall_line_count, infill_origin, 
        perimeter_gaps, connected_zigzags, use_endpieces, false /* skip_some_zags */, zag_skip_count, pocket_size);
    Polygons polygons;
    Polygons lines;
    roof_computation.generate(toolpaths, polygons, lines, config.settings);
    append(lines, to_polylines(polygons));
    return lines;
#else
#ifdef _WIN32
    if (! BoundingBox(Point::new_scale(-170., -170.), Point::new_scale(170., 170.)).contains(get_extents(polygon)))
        ::MessageBoxA(nullptr, "TreeSupport infill kravsky", "Bug detected!", MB_OK | MB_SYSTEMMODAL | MB_SETFOREGROUND | MB_ICONWARNING);
#endif // _WIN32

    const Flow            &flow   = roof ? support_params.support_material_interface_flow : support_params.support_material_flow;
    std::unique_ptr<Fill>  filler = std::unique_ptr<Fill>(Fill::new_from_type(roof ? support_params.interface_fill_pattern : support_params.base_fill_pattern));
    FillParams             fill_params;

    filler->layer_id = layer_idx;
    filler->spacing  = flow.spacing();
    filler->angle = roof ? 
        //fixme support_layer.interface_id() instead of layer_idx
        (support_params.interface_angle + (layer_idx & 1) ? float(- M_PI / 4.) : float(+ M_PI / 4.)) :
        support_params.base_angle;

    fill_params.density     = float(roof ? support_params.interface_density : scaled<float>(filler->spacing) / (scaled<float>(filler->spacing) + float(support_infill_distance)));
    fill_params.dont_adjust = true;

    Polylines out;
    for (ExPolygon &expoly : union_ex(polygon)) {
        // The surface type does not matter.
        assert(area(expoly) > 0.);
#ifdef _WIN32
        if (area(expoly) <= 0.)
            ::MessageBoxA(nullptr, "TreeSupport infill negative area", "Bug detected!", MB_OK | MB_SYSTEMMODAL | MB_SETFOREGROUND | MB_ICONWARNING);
#endif // _WIN32
        assert(intersecting_edges(to_polygons(expoly)).empty());
#ifdef _WIN32
        if (! intersecting_edges(to_polygons(expoly)).empty())
            ::MessageBoxA(nullptr, "TreeSupport infill self intersections", "Bug detected!", MB_OK | MB_SYSTEMMODAL | MB_SETFOREGROUND | MB_ICONWARNING);
#endif // _WIN32
        Surface surface(stInternal, std::move(expoly));
        try {
            Polylines pl = filler->fill_surface(&surface, fill_params);
            assert(pl.empty() || get_extents(surface.expolygon).inflated(SCALED_EPSILON).contains(get_extents(pl)));
#ifdef _WIN32
            if (! pl.empty() && ! get_extents(surface.expolygon).inflated(SCALED_EPSILON).contains(get_extents(pl)))
                ::MessageBoxA(nullptr, "TreeSupport infill failure", "Bug detected!", MB_OK | MB_SYSTEMMODAL | MB_SETFOREGROUND | MB_ICONWARNING);
#endif // _WIN32
            append(out, std::move(pl));
        } catch (InfillFailedException &) {
        }
    }
    validate_range(out);
    return out;
#endif
}

/*!
 * \brief Unions two Polygons. Ensures that if the input is non empty that the output also will be non empty.
 * \param first[in] The first Polygon.
 * \param second[in] The second Polygon.
 * \return The union of both Polygons
 */
[[nodiscard]] static Polygons safeUnion(const Polygons first, const Polygons second = Polygons())
{
    // unionPolygons can slowly remove Polygons under certain circumstances, because of rounding issues (Polygons that have a thin area).
    // This does not cause a problem when actually using it on large areas, but as influence areas (representing centerpoints) can be very thin, this does occur so this ugly workaround is needed
    // Here is an example of a Polygons object that will loose vertices when unioning, and will be gone after a few times unionPolygons was called:
    /*
    Polygons example;
    Polygon exampleInner;
    exampleInner.add(Point(120410,83599));//A
    exampleInner.add(Point(120384,83643));//B
    exampleInner.add(Point(120399,83618));//C
    exampleInner.add(Point(120414,83591));//D
    exampleInner.add(Point(120423,83570));//E
    exampleInner.add(Point(120419,83580));//F
    example.add(exampleInner);
    for(int i=0;i<10;i++){
         log("Iteration %d Example area: %f\n",i,area(example));
         example=example.unionPolygons();
    }
*/

    Polygons result;
    if (! first.empty() || ! second.empty()) {
        result = union_(first, second);
        if (result.empty()) {
            BOOST_LOG_TRIVIAL(debug) << "Caught an area destroying union, enlarging areas a bit.";
            // just take the few lines we have, and offset them a tiny bit. Needs to be offsetPolylines, as offset may aleady have problems with the area.
            result = union_(offset(to_polylines(first), scaled<float>(0.002), jtMiter, 1.2), offset(to_polylines(second), scaled<float>(0.002), jtMiter, 1.2));
        }
    }
    
    return result;
}

/*!
 * \brief Offsets (increases the area of) a polygons object in multiple steps to ensure that it does not lag through over a given obstacle.
 * \param me[in] Polygons object that has to be offset.
 * \param distance[in] The distance by which me should be offset. Expects values >=0.
 * \param collision[in] The area representing obstacles.
 * \param last_step_offset_without_check[in] The most it is allowed to offset in one step.
 * \param min_amount_offset[in] How many steps have to be done at least. As this uses round offset this increases the amount of vertices, which may be required if Polygons get very small. Required as arcTolerance is not exposed in offset, which should result with a similar result.
 * \return The resulting Polygons object.
 */
[[nodiscard]] static Polygons safeOffsetInc(const Polygons& me, coord_t distance, const Polygons& collision, coord_t safe_step_size, coord_t last_step_offset_without_check, size_t min_amount_offset)
{
    bool do_final_difference = last_step_offset_without_check == 0;
    Polygons ret = safeUnion(me); // ensure sane input
    
    // Trim the collision polygons with the region of interest for diff() efficiency.
    Polygons collision_trimmed_buffer;
    auto collision_trimmed = [&collision_trimmed_buffer, &collision, &ret, distance]() -> const Polygons& {
        if (collision_trimmed_buffer.empty() && ! collision.empty())
            collision_trimmed_buffer = clip_for_diff(collision, get_extents(ret).inflated(std::max(0, distance) + SCALED_EPSILON));
        return collision_trimmed_buffer;
    };

    if (distance == 0)
        return do_final_difference ? diff(ret, collision_trimmed()) : union_(ret);
    if (safe_step_size < 0 || last_step_offset_without_check < 0) {
        BOOST_LOG_TRIVIAL(error) << "Offset increase got invalid parameter!";
        TreeSupport::showError("Negative offset distance... How did you manage this ?", true);
        return do_final_difference ? diff(ret, collision_trimmed()) : union_(ret);
    }

    coord_t step_size = safe_step_size;
    int     steps = distance > last_step_offset_without_check ? (distance - last_step_offset_without_check) / step_size : 0;
    if (distance - steps * step_size > last_step_offset_without_check) {
        if ((steps + 1) * step_size <= distance)
            // This will be the case when last_step_offset_without_check >= safe_step_size
            ++ steps;
        else
            do_final_difference = true;
    }
    if (steps + (distance < last_step_offset_without_check || distance % step_size != 0) < min_amount_offset && min_amount_offset > 1) {
        // yes one can add a bool as the standard specifies that a result from compare operators has to be 0 or 1
        // reduce the stepsize to ensure it is offset the required amount of times
        step_size = distance / min_amount_offset;
        if (step_size >= safe_step_size) {
            // effectivly reduce last_step_offset_without_check
            step_size = safe_step_size;
            steps = min_amount_offset;
        } else
            steps = distance / step_size;
    }
    // offset in steps
    for (int i = 0; i < steps; ++ i) {
        ret = diff(offset(ret, step_size, ClipperLib::jtRound, scaled<float>(0.01)), collision_trimmed());
        // ensure that if many offsets are done the performance does not suffer extremely by the new vertices of jtRound.
        if (i % 10 == 7)
            ret = polygons_simplify(ret, scaled<double>(0.015));
    }
    // offset the remainder
    float last_offset = distance - steps * step_size;
    if (last_offset > SCALED_EPSILON)
        ret = offset(ret, distance - steps * step_size, ClipperLib::jtRound, scaled<float>(0.01));
    ret = polygons_simplify(ret, scaled<double>(0.015));

    if (do_final_difference)
        ret = diff(ret, collision_trimmed());
    return union_(ret);
}

static inline SupportGeneratorLayer& layer_initialize(
    SupportGeneratorLayer   &layer_new,
    const SupporLayerType    layer_type,
    const SlicingParameters &slicing_params,
    const size_t             layer_idx)
{
    layer_new.layer_type = layer_type;
    layer_new.print_z  = slicing_params.object_print_z_min + slicing_params.first_object_layer_height + layer_idx * slicing_params.layer_height;
    layer_new.height   = layer_idx == 0 ? slicing_params.first_object_layer_height : slicing_params.layer_height;
    layer_new.bottom_z = layer_idx == 0 ? slicing_params.object_print_z_min : layer_new.print_z - layer_new.height;
    return layer_new;
}

// Using the std::deque as an allocator.
inline SupportGeneratorLayer& layer_allocate(
    std::deque<SupportGeneratorLayer> &layer_storage,
    SupporLayerType                    layer_type,
    const SlicingParameters           &slicing_params,
    size_t                             layer_idx)
{
    //FIXME take raft into account.
    layer_storage.push_back(SupportGeneratorLayer());
    return layer_initialize(layer_storage.back(), layer_type, slicing_params, layer_idx);
}

inline SupportGeneratorLayer& layer_allocate(
    std::deque<SupportGeneratorLayer> &layer_storage,
    tbb::spin_mutex&                   layer_storage_mutex,
    SupporLayerType                    layer_type,
    const SlicingParameters           &slicing_params,
    size_t                             layer_idx)
{
    tbb::spin_mutex::scoped_lock lock(layer_storage_mutex);
    layer_storage.push_back(SupportGeneratorLayer());
    return layer_initialize(layer_storage.back(), layer_type, slicing_params, layer_idx);
}

void TreeSupport::generateInitialAreas(
    const PrintObject                       &print_object,
    const std::vector<Polygons>             &overhangs,
    std::vector<std::set<SupportElement*>>  &move_bounds,
    SupportGeneratorLayersPtr               &top_contacts,
    SupportGeneratorLayersPtr               &top_interface_layers,
    SupportGeneratorLayerStorage            &layer_storage)
{
    Polygon base_circle;
    const auto base_radius = scaled<int>(0.01);
    for (unsigned int i = 0; i < SUPPORT_TREE_CIRCLE_RESOLUTION; ++ i) {
        const double angle = static_cast<double>(i) / SUPPORT_TREE_CIRCLE_RESOLUTION * (2.0 * M_PI);
        base_circle.points.emplace_back(coord_t(cos(angle) * base_radius), coord_t(sin(angle) * base_radius));
    }
    TreeSupportMeshGroupSettings    mesh_group_settings(print_object);
    TreeSupportSettings             mesh_config{ mesh_group_settings };
    SupportParameters               support_params(print_object);
    support_params.with_sheath = true;
    support_params.support_density = 0;

    const size_t z_distance_delta = mesh_config.z_distance_top_layers + 1; // To ensure z_distance_top_layers are left empty between the overhang (zeroth empty layer), the support has to be added z_distance_top_layers+1 layers below

    const bool xy_overrides_z = mesh_config.support_xy_overrides_z;
#if 0
    if (mesh.overhang_areas.size() <= z_distance_delta)
        return;
#endif

    const coord_t connect_length = (mesh_config.support_line_width * 100. / mesh_group_settings.support_tree_top_rate) + std::max(2. * mesh_config.min_radius - 1.0 * mesh_config.support_line_width, 0.0);
    // As r*r=x*x+y*y (circle equation): If a circle with center at (0,0) the top most point is at (0,r) as in y=r.
    // This calculates how far one has to move on the x-axis so that y=r-support_line_width/2. 
    // In other words how far does one need to move on the x-axis to be support_line_width/2 away from the circle line.
    // As a circle is round this length is identical for every axis as long as the 90 degrees angle between both remains.
    const coord_t circle_length_to_half_linewidth_change = mesh_config.min_radius < mesh_config.support_line_width ? mesh_config.min_radius / 2 : sqrt(sqr(mesh_config.min_radius) - sqr(mesh_config.min_radius - mesh_config.support_line_width / 2));
    // Extra support offset to compensate for larger tip radiis. Also outset a bit more when z overwrites xy, because supporting something with a part of a support line is better than not supporting it at all.
    //FIXME Vojtech: This is not sufficient for support enforcers to work.
    //FIXME There is no account for the support overhang angle.
    //FIXME There is no account for the width of the collision regions.
    const coord_t extra_outset = std::max(coord_t(0), mesh_config.min_radius - mesh_config.support_line_width) + (xy_overrides_z ? 0 : mesh_config.support_line_width / 2)
        //FIXME this is a heuristic value for support enforcers to work.
//        + 10 * mesh_config.support_line_width;
        ;
    const size_t support_roof_layers = mesh_group_settings.support_roof_enable ? (mesh_group_settings.support_roof_height + mesh_config.layer_height / 2) / mesh_config.layer_height : 0;
    const bool roof_enabled = support_roof_layers != 0;
    const bool force_tip_to_roof = sqr<double>(mesh_config.min_radius) * M_PI > mesh_group_settings.minimum_roof_area && roof_enabled;
    //FIXME mesh_group_settings.support_angle does not apply to enforcers and also it does not apply to automatic support angle (by half the external perimeter width).
    const coord_t max_overhang_speed = (mesh_group_settings.support_angle < 0.5 * M_PI) ? (coord_t)(tan(mesh_group_settings.support_angle) * mesh_config.layer_height) : std::numeric_limits<coord_t>::max();
    const size_t max_overhang_insert_lag = std::max((size_t)round_up_divide(mesh_config.xy_distance, max_overhang_speed / 2), 2 * mesh_config.z_distance_top_layers); // cap for how much layer below the overhang a new support point may be added, as other than with regular support every new inserted point may cause extra material and time cost.  Could also be an user setting or differently calculated. Idea is that if an overhang does not turn valid in double the amount of layers a slope of support angle would take to travel xy_distance, nothing reasonable will come from it. The 2*z_distance_delta is only a catch for when the support angle is very high.

    //FIXME 
    size_t num_support_layers = print_object.layer_count();
    std::vector<std::unordered_set<Point, PointHash>> already_inserted(num_support_layers - z_distance_delta);

    std::mutex mutex_layer_storage, mutex_movebounds;
    tbb::parallel_for(tbb::blocked_range<size_t>(1, num_support_layers - z_distance_delta),
        [this, &print_object, &overhangs, &mesh_config, &mesh_group_settings, &support_params, 
         z_distance_delta, xy_overrides_z, force_tip_to_roof, roof_enabled, support_roof_layers, extra_outset, circle_length_to_half_linewidth_change, connect_length, max_overhang_insert_lag,
         &base_circle, &mutex_layer_storage, &mutex_movebounds, &top_contacts, &layer_storage, &already_inserted,
         &move_bounds, &base_radius](const tbb::blocked_range<size_t> &range) {
        for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx) {
            if (overhangs[layer_idx + z_distance_delta].empty())
                continue;
            // take the least restrictive avoidance possible
            Polygons relevant_forbidden;
            {
                const Polygons &relevant_forbidden_raw = (mesh_config.support_rests_on_model ?
                    (SUPPORT_TREE_ONLY_GRACIOUS_TO_MODEL ? m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx, AvoidanceType::Fast, true, !xy_overrides_z) :
                        m_volumes.getCollision(mesh_config.getRadius(0), layer_idx, !xy_overrides_z)) :
                    m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx, AvoidanceType::Fast, false, !xy_overrides_z));
                // prevent rounding errors down the line, points placed directly on the line of the forbidden area may not be added otherwise.
                relevant_forbidden = offset(union_ex(relevant_forbidden_raw), scaled<float>(0.005), jtMiter, 1.2);
            }

            auto generateLines = [&](const Polygons& area, bool roof, LayerIndex layer_idx) -> Polylines {
                const coord_t support_infill_distance = roof ? mesh_group_settings.support_roof_line_distance : mesh_group_settings.support_tree_branch_distance;
                return generateSupportInfillLines(area, support_params, roof, layer_idx, support_infill_distance);
            };

            auto addPointAsInfluenceArea = [&](std::pair<Point, LineStatus> p, size_t dtt, LayerIndex insert_layer, size_t dont_move_until, bool roof, bool skip_ovalisation)
            {
                bool to_bp = p.second == LineStatus::TO_BP || p.second == LineStatus::TO_BP_SAFE;
                bool gracious = to_bp || p.second == LineStatus::TO_MODEL_GRACIOUS || p.second == LineStatus::TO_MODEL_GRACIOUS_SAFE;
                bool safe_radius = p.second == LineStatus::TO_BP_SAFE || p.second == LineStatus::TO_MODEL_GRACIOUS_SAFE;
                if (!mesh_config.support_rests_on_model && !to_bp) {
                    BOOST_LOG_TRIVIAL(warning) << "Tried to add an invalid support point";
                    TreeSupport::showError("Unable to add tip. Some overhang may not be supported correctly.", true);
                    return;
                }
                Polygon circle;
                for (Point corner : base_circle)
                    circle.points.emplace_back(p.first + corner);
                {
                    std::lock_guard<std::mutex> critical_section_movebounds(mutex_movebounds);
                    if (! already_inserted[insert_layer].count(p.first / ((mesh_config.min_radius + 1) / 10))) {
                        // normalize the point a bit to also catch points which are so close that inserting it would achieve nothing
                        already_inserted[insert_layer].emplace(p.first / ((mesh_config.min_radius + 1) / 10));
                        SupportElement* elem = new SupportElement(dtt, insert_layer, p.first, to_bp, gracious, !xy_overrides_z, dont_move_until, roof, safe_radius, force_tip_to_roof, skip_ovalisation);
                        elem->area = new Polygons();
                        validate_range(circle);
                        elem->area->emplace_back(std::move(circle));
                        move_bounds[insert_layer].emplace(elem);
                    }
                }
            };

            auto addLinesAsInfluenceAreas = [&](LineInformations lines, size_t roof_tip_layers, LayerIndex insert_layer_idx, bool supports_roof, size_t dont_move_until)
            {
                validate_range(lines);
                // Add tip area as roof (happens when minimum roof area > minimum tip area) if possible
                size_t dtt_roof_tip;
                for (dtt_roof_tip = 0; dtt_roof_tip < roof_tip_layers && insert_layer_idx - dtt_roof_tip >= 1; dtt_roof_tip++)
                {
                    auto evaluateRoofWillGenerate = [&](std::pair<Point, LineStatus> p) {
                        //FIXME Vojtech: The circle is just shifted, it has a known size, the infill should fit all the time!
#if 0
                        Polygon roof_circle;
                        for (Point corner : base_circle)
                            roof_circle.points.emplace_back(p.first + corner * mesh_config.min_radius);
                        return !generateSupportInfillLines({ roof_circle }, mesh_config, true, insert_layer_idx - dtt_roof_tip, mesh_config.support_roof_line_distance).empty();
#else
                        return true;
#endif
                    };

                    std::pair<LineInformations, LineInformations> split = 
                        // keep all lines that are still valid on the next layer
                        splitLines(lines, [this, insert_layer_idx, dtt_roof_tip](std::pair<Point, LineStatus> &p){ return evaluatePointForNextLayerFunction(m_volumes, m_config, insert_layer_idx - dtt_roof_tip, p); });

                    for (LineInformation line : split.second) // add all points that would not be valid
                        for (std::pair<Point, LineStatus> point_data : line)
                            addPointAsInfluenceArea(point_data, 0, insert_layer_idx - dtt_roof_tip, roof_tip_layers - dtt_roof_tip, dtt_roof_tip != 0, false);

                    // not all roofs are guaranteed to actually generate lines, so filter these out and add them as points
                    split = splitLines(split.first, evaluateRoofWillGenerate);
                    lines = split.first;

                    for (LineInformation line : split.second)
                        for (std::pair<Point, LineStatus> point_data : line)
                            addPointAsInfluenceArea(point_data, 0, insert_layer_idx - dtt_roof_tip, roof_tip_layers - dtt_roof_tip, dtt_roof_tip != 0, false);

                    // add all tips as roof to the roof storage
                    Polygons added_roofs;
                    for (LineInformation line : lines)
                        for (std::pair<Point, LineStatus> p : line) {
                            Polygon roof_circle;
                            for (Point corner : base_circle)
                                roof_circle.points.emplace_back(p.first + corner * mesh_config.min_radius / base_radius);
                            added_roofs.emplace_back(roof_circle);
                        }
                    if (! added_roofs.empty()) {
                        added_roofs = union_(added_roofs);
                        {
                            std::lock_guard<std::mutex> lock(mutex_layer_storage);
                            SupportGeneratorLayer *&l = top_contacts[insert_layer_idx - dtt_roof_tip];
                            if (l == nullptr)
                                l = &layer_allocate(layer_storage, SupporLayerType::TopContact, print_object.slicing_parameters(), insert_layer_idx - dtt_roof_tip);
                            append(l->polygons, std::move(added_roofs));
                        }
                    }
                }

                for (LineInformation line : lines) {
                    bool disable_ovalistation = mesh_config.min_radius < 3 * mesh_config.support_line_width && roof_tip_layers == 0 && dtt_roof_tip == 0 && line.size() > 5; // If a line consists of enough tips, the assumption is that it is not a single tip, but part of a simulated support pattern. Ovalisation should be disabled for these to improve the quality of the lines when tip_diameter=line_width
                    for (auto point_data : line)
                        addPointAsInfluenceArea(point_data, 0, insert_layer_idx - dtt_roof_tip, dont_move_until > dtt_roof_tip ? dont_move_until - dtt_roof_tip : 0, dtt_roof_tip != 0 || supports_roof, disable_ovalistation);
                }
            };

            // every overhang has saved if a roof should be generated for it. This can NOT be done in the for loop as an area may NOT have a roof 
            // even if it is larger than the minimum_roof_area when it is only larger because of the support horizontal expansion and 
            // it would not have a roof if the overhang is offset by support roof horizontal expansion instead. (At least this is the current behavior of the regular support)
            Polygons overhang_regular;
            {
                const Polygons &overhang_raw = overhangs[layer_idx + z_distance_delta];
                overhang_regular = mesh_group_settings.support_offset == 0 ? 
                    overhang_raw :
                    safeOffsetInc(overhang_raw, mesh_group_settings.support_offset, relevant_forbidden, mesh_config.min_radius * 1.75 + mesh_config.xy_min_distance, 0, 1);
                // offset ensures that areas that could be supported by a part of a support line, are not considered unsupported overhang
                Polygons remaining_overhang = intersection(
                    diff(mesh_group_settings.support_offset == 0 ?
                            overhang_raw :
                            offset(union_ex(overhang_raw), mesh_group_settings.support_offset, jtMiter, 1.2),
                         offset(union_ex(overhang_regular), mesh_config.support_line_width * 0.5, jtMiter, 1.2)),
                    relevant_forbidden);

                // Offset the area to compensate for large tip radiis. Offset happens in multiple steps to ensure the tip is as close to the original overhang as possible.
                //+mesh_config.support_line_width / 80  to avoid calculating very small (useless) offsets because of rounding errors.
                //FIXME likely a better approach would be to find correspondences between the full overhang and the trimmed overhang
                // and if there is no correspondence, project the missing points to the clipping curve.
                for (coord_t extra_total_offset_acc = 0; ! remaining_overhang.empty() && extra_total_offset_acc + mesh_config.support_line_width / 8 < extra_outset; ) {
                    const coord_t offset_current_step = std::min(
                        extra_total_offset_acc + 2 * mesh_config.support_line_width > mesh_config.min_radius ?
                            mesh_config.support_line_width / 8 : 
                            circle_length_to_half_linewidth_change,
                        extra_outset - extra_total_offset_acc);
                    extra_total_offset_acc += offset_current_step;
                    const Polygons &raw_collision = m_volumes.getCollision(0, layer_idx, true);
                    const coord_t   offset_step   = mesh_config.xy_min_distance + mesh_config.support_line_width;
                    // Reducing the remaining overhang by the areas already supported.
                    //FIXME 1.5 * extra_total_offset_acc seems to be too much, it may remove some remaining overhang without being supported at all.
                    remaining_overhang = diff(remaining_overhang, safeOffsetInc(overhang_regular, 1.5 * extra_total_offset_acc, raw_collision, offset_step, 0, 1));
                    // Extending the overhangs by the inflated remaining overhangs.
                    overhang_regular   = union_(overhang_regular, diff(safeOffsetInc(remaining_overhang, extra_total_offset_acc, raw_collision, offset_step, 0, 1), relevant_forbidden));
                }
                // If the xy distance overrides the z distance, some support needs to be inserted further down.
                //=> Analyze which support points do not fit on this layer and check if they will fit a few layers down (while adding them an infinite amount of layers down would technically be closer the the setting description, it would not produce reasonable results. )
                if (xy_overrides_z)
                {
                    LineInformations overhang_lines;
                    {
                        //Vojtech: Generate support heads at support_tree_branch_distance spacing by producing a zig-zag infill at support_tree_branch_distance spacing,
                        // which is then resmapled 
                        // support_line_width to form a line here as otherwise most will be unsupported. Technically this violates branch distance, 
                        // mbut not only is this the only reasonable choice, but it ensures consistent behavior as some infill patterns generate 
                        // each line segment as its own polyline part causing a similar line forming behavior. Also it is assumed that 
                        // the area that is valid a layer below is to small for support roof.
                        Polylines polylines = ensureMaximumDistancePolyline(generateLines(remaining_overhang, false, layer_idx), mesh_config.min_radius, 1);
                        if (polylines.size() <= 3)
                            // add the outer wall to ensure it is correct supported instead
                            polylines = ensureMaximumDistancePolyline(to_polylines(remaining_overhang), connect_length, 3);
                        for (const auto &line : polylines) {
                            LineInformation res_line;
                            for (Point p : line)
                                res_line.emplace_back(p, LineStatus::INVALID);
                            overhang_lines.emplace_back(res_line);
                        }
                        validate_range(overhang_lines);
                    }

                    for (size_t lag_ctr = 1; lag_ctr <= max_overhang_insert_lag && !overhang_lines.empty() && layer_idx - coord_t(lag_ctr) >= 1; lag_ctr++) {
                        // get least restricted avoidance for layer_idx-lag_ctr
                        const Polygons &relevant_forbidden_below = (mesh_config.support_rests_on_model ? (SUPPORT_TREE_ONLY_GRACIOUS_TO_MODEL ? m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx - lag_ctr, AvoidanceType::Fast, true, !xy_overrides_z) : m_volumes.getCollision(mesh_config.getRadius(0), layer_idx - lag_ctr, !xy_overrides_z)) : m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx - lag_ctr, AvoidanceType::Fast, false, !xy_overrides_z));
                        // it is not required to offset the forbidden area here as the points wont change: If points here are not inside the forbidden area neither will they be later when placing these points, as these are the same points.
                        auto evaluatePoint = [&](std::pair<Point, LineStatus> p) { return contains(relevant_forbidden_below, p.first); };

                        std::pair<LineInformations, LineInformations> split = splitLines(overhang_lines, evaluatePoint); // keep all lines that are invalid
                        overhang_lines = split.first;
                        LineInformations fresh_valid_points = convertLinesToInternal(m_volumes, m_config, convertInternalToLines(split.second), layer_idx - lag_ctr); // set all now valid lines to their correct LineStatus. Easiest way is to just discard Avoidance information for each point and evaluate them again.
                        validate_range(fresh_valid_points);

                        addLinesAsInfluenceAreas(fresh_valid_points, (force_tip_to_roof && lag_ctr <= support_roof_layers) ? support_roof_layers : 0, layer_idx - lag_ctr, false, roof_enabled ? support_roof_layers : 0);
                    }
                }
            }

            Polygons overhang_roofs;
            std::vector<std::pair<ExPolygon, bool>> overhang_processing; 
            if (roof_enabled) {
                static constexpr const coord_t support_roof_offset = 0;
                overhang_roofs = safeOffsetInc(overhangs[layer_idx + z_distance_delta], support_roof_offset, relevant_forbidden, mesh_config.min_radius * 2 + mesh_config.xy_min_distance, 0, 1);
                if (mesh_group_settings.minimum_support_area > 0)
                    remove_small(overhang_roofs, mesh_group_settings.minimum_roof_area);
                overhang_regular = diff(overhang_regular, overhang_roofs, ApplySafetyOffset::Yes);
                for (ExPolygon &roof_part : union_ex(overhang_roofs))
                    overhang_processing.emplace_back(std::move(roof_part), true);
            }
            if (mesh_group_settings.minimum_support_area > 0)
                remove_small(overhang_regular, mesh_group_settings.minimum_support_area);

            for (ExPolygon &support_part : union_ex(overhang_regular))
                overhang_processing.emplace_back(std::move(support_part), false);

            for (const std::pair<ExPolygon, bool> &overhang_pair : overhang_processing) {
                const bool roof_allowed_for_this_part = overhang_pair.second;
                Polygons overhang_outset = to_polygons(overhang_pair.first);
                const size_t min_support_points = std::max(coord_t(1), std::min(coord_t(3), coord_t(total_length(overhang_outset) / connect_length)));
                LineInformations overhang_lines;
                Polygons last_overhang = overhang_outset;
                size_t dtt_roof = 0;
                // Sometimes roofs could be empty as the pattern does not generate lines if the area is narrow enough (i am looking at you, concentric infill).
                // To catch these cases the added roofs are saved to be evaluated later.
                std::vector<Polygons> added_roofs(support_roof_layers);

                // Assumption is that roof will support roof further up to avoid a lot of unnecessary branches. Each layer down it is checked whether the roof area 
                // is still large enough to be a roof and aborted as soon as it is not. This part was already reworked a few times, and there could be an argument 
                // made to change it again if there are actual issues encountered regarding supporting roofs.
                // Main problem is that some patterns change each layer, so just calculating points and checking if they are still valid an layer below is not useful, 
                // as the pattern may be different one layer below. Same with calculating which points are now no longer being generated as result from 
                // a decreasing roof, as there is no guarantee that a line will be above these points. Implementing a separate roof support behavior
                // for each pattern harms maintainability as it very well could be >100 LOC
                if (roof_allowed_for_this_part) {
                    for (dtt_roof = 0; dtt_roof < support_roof_layers && layer_idx - dtt_roof >= 1; dtt_roof++) {
                        // here the roof is handled. If roof can not be added the branches will try to not move instead
                        Polygons forbidden_next;
                        {
                            const Polygons &forbidden_next_raw = mesh_config.support_rests_on_model ? 
                                (SUPPORT_TREE_ONLY_GRACIOUS_TO_MODEL ? 
                                    m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx - (dtt_roof + 1), AvoidanceType::Fast, true, !xy_overrides_z) : 
                                    m_volumes.getCollision(mesh_config.getRadius(0), layer_idx - (dtt_roof + 1), !xy_overrides_z)) : 
                                m_volumes.getAvoidance(mesh_config.getRadius(0), layer_idx - (dtt_roof + 1), AvoidanceType::Fast, false, !xy_overrides_z);
                            // prevent rounding errors down the line
                            forbidden_next = offset(union_ex(forbidden_next_raw), scaled<float>(0.005), jtMiter, 1.2);
                        }
                        Polygons overhang_outset_next = diff(overhang_outset, forbidden_next);
                        if (area(overhang_outset_next) < mesh_group_settings.minimum_roof_area) {
                            // next layer down the roof area would be to small so we have to insert our roof support here. Also convert squaremicrons to squaremilimeter
                            size_t dtt_before = dtt_roof > 0 ? dtt_roof - 1 : 0;
                            if (dtt_roof != 0) {
                                // Produce support head points supporting an interface layer: First produce the interface lines, then sample them.
                                overhang_lines = convertLinesToInternal(m_volumes, m_config, 
                                    ensureMaximumDistancePolyline(generateLines(last_overhang, true, layer_idx - dtt_before), connect_length, 1), layer_idx - dtt_before);
                                overhang_lines = splitLines(overhang_lines, 
                                    [this, layer_idx, dtt_before](std::pair<Point, LineStatus> &p){ return evaluatePointForNextLayerFunction(m_volumes, m_config, layer_idx - dtt_before, p); }).first;
                            }
                            break;
                        }
                        added_roofs[dtt_roof] = overhang_outset;
                        last_overhang = overhang_outset;
                        overhang_outset = overhang_outset_next;
                    }
                }

                size_t layer_generation_dtt = std::max(dtt_roof, size_t(1)) - 1; // 1 inside max and -1 outside to avoid underflow. layer_generation_dtt=dtt_roof-1 if dtt_roof!=0;
                // if the roof should be valid, check that the area does generate lines. This is NOT guaranteed.
                if (overhang_lines.empty() && dtt_roof != 0 && generateLines(overhang_outset, true, layer_idx - layer_generation_dtt).empty())
                    for (size_t idx = 0; idx < dtt_roof; idx++) {
                        // check for every roof area that it has resulting lines. Remember idx 1 means the 2. layer of roof => higher idx == lower layer
                        if (generateLines(added_roofs[idx], true, layer_idx - idx).empty()) {
                            dtt_roof = idx;
                            layer_generation_dtt = std::max(dtt_roof, size_t(1)) - 1;
                            break;
                        }
                    }

                {
                    std::lock_guard<std::mutex> lock(mutex_layer_storage);
                    for (size_t idx = 0; idx < dtt_roof; ++ idx)
                        if (! added_roofs[idx].empty()) {
                            SupportGeneratorLayer *&l = top_contacts[layer_idx - idx];
                            if (l == nullptr)
                                l = &layer_allocate(layer_storage, SupporLayerType::TopContact, print_object.slicing_parameters(), layer_idx - idx);
                            // will be unioned in finalizeInterfaceAndSupportAreas
                            append(l->polygons, std::move(added_roofs[idx]));
                        }
                }

                if (overhang_lines.empty()) {
                    // support_line_width to form a line here as otherwise most will be unsupported. Technically this violates branch distance, but not only is this the only reasonable choice,
                    // but it ensures consistant behaviour as some infill patterns generate each line segment as its own polyline part causing a similar line forming behaviour. 
                    // This is not doen when a roof is above as the roof will support the model and the trees only need to support the roof
                    Polylines polylines = ensureMaximumDistancePolyline(generateLines(overhang_outset, dtt_roof != 0, layer_idx - layer_generation_dtt), dtt_roof == 0 ? mesh_config.min_radius / 2 : connect_length, 1); 
                    size_t point_count = 0;
                    for (const Polyline &poly : polylines)
                        point_count += poly.size();
                    if (point_count <= min_support_points) {
                        // add the outer wall (of the overhang) to ensure it is correct supported instead. Try placing the support points in a way that they fully support the outer wall, instead of just the with half of the the support line width.
                        // I assume that even small overhangs are over one line width wide, so lets try to place the support points in a way that the full support area generated from them 
                        // will support the overhang (if this is not done it may only be half). This WILL NOT be the case when supporting an angle of about < 60� so there is a fallback, 
                        // as some support is better than none.
                        Polygons reduced_overhang_outset = offset(union_ex(overhang_outset), -mesh_config.support_line_width / 2.2, jtMiter, 1.2);
                        polylines = ensureMaximumDistancePolyline(
                            to_polylines(!reduced_overhang_outset.empty() && area(offset(diff_ex(overhang_outset, reduced_overhang_outset), std::max(mesh_config.support_line_width, connect_length), jtMiter, 1.2)) < sqr(scaled<double>(0.001)) ?
                                reduced_overhang_outset :
                                overhang_outset),
                            connect_length, min_support_points);
                    }
                    LayerIndex last_insert_layer = layer_idx - dtt_roof;
                    overhang_lines = convertLinesToInternal(m_volumes, m_config, polylines, last_insert_layer);
                }

                if (int(dtt_roof) >= layer_idx && roof_allowed_for_this_part && ! overhang_outset.empty()) {
                    // reached buildplate
                    std::lock_guard<std::mutex> lock(mutex_layer_storage);
                    SupportGeneratorLayer*& l = top_contacts[0];
                    if (l == nullptr)
                        l = &layer_allocate(layer_storage, SupporLayerType::TopContact, print_object.slicing_parameters(), 0);
                    append(l->polygons, std::move(overhang_outset));
                } else // normal trees have to be generated
                    addLinesAsInfluenceAreas(overhang_lines, force_tip_to_roof ? support_roof_layers - dtt_roof : 0, layer_idx - dtt_roof, dtt_roof > 0, roof_enabled ? support_roof_layers - dtt_roof : 0);
            }
        }
    });
}

static unsigned int moveInside(const Polygons &polygons, Point &from, int distance = 0, int64_t maxDist2 = std::numeric_limits<int64_t>::max())
{
    Point  ret = from;
    double bestDist2 = std::numeric_limits<double>::max();
    auto   bestPoly = static_cast<unsigned int>(-1);
    bool   is_already_on_correct_side_of_boundary = false; // whether [from] is already on the right side of the boundary
    for (unsigned int poly_idx = 0; poly_idx < polygons.size(); ++ poly_idx) {
        const Polygon &poly = polygons[poly_idx];
        if (poly.size() < 2)
            continue;
        Point p0 = poly[poly.size() - 2];
        Point p1 = poly.back();
        // because we compare with vSize2 here (no division by zero), we also need to compare by vSize2 inside the loop
        // to avoid integer rounding edge cases
        bool projected_p_beyond_prev_segment = (p1 - p0).cast<int64_t>().dot((from - p0).cast<int64_t>()) >= (p1 - p0).cast<int64_t>().squaredNorm();
        for (const Point& p2 : poly) {
            // X = A + Normal(B-A) * (((B-A) dot (P-A)) / VSize(B-A));
            //   = A +       (B-A) *  ((B-A) dot (P-A)) / VSize2(B-A);
            // X = P projected on AB
            const Point& a = p1;
            const Point& b = p2;
            const Point& p = from;
            auto ab = (b - a).cast<int64_t>();
            auto ap = (p - a).cast<int64_t>();
            int64_t ab_length2 = ab.squaredNorm();
            if (ab_length2 <= 0) { //A = B, i.e. the input polygon had two adjacent points on top of each other.
                p1 = p2; //Skip only one of the points.
                continue;
            }
            int64_t dot_prod = ab.dot(ap);
            if (dot_prod <= 0) { // x is projected to before ab
                if (projected_p_beyond_prev_segment) { 
                    //  case which looks like:   > .
                    projected_p_beyond_prev_segment = false;
                    Point& x = p1;

                    auto dist2 = (x - p).cast<int64_t>().squaredNorm();
                    if (dist2 < bestDist2) {
                        bestDist2 = dist2;
                        bestPoly = poly_idx;
                        if (distance == 0)
                            ret = x;
                        else {
                            Vec2d  abd   = ab.cast<double>();
                            Vec2d  p1p2  = (p1 - p0).cast<double>();
                            double lab   = abd.norm();
                            double lp1p2 = p1p2.norm();
                            // inward direction irrespective of sign of [distance]
                            auto inward_dir = perp(abd * (scaled<double>(10.0) / lab) + p1p2 * (scaled<double>(10.0) / lp1p2));
                            // MM2INT(10.0) to retain precision for the eventual normalization
                            ret = x + (inward_dir * (distance / inward_dir.norm())).cast<coord_t>();
                            is_already_on_correct_side_of_boundary = inward_dir.dot((p - x).cast<double>()) * distance >= 0;
                        }
                    }
                } else {
                    projected_p_beyond_prev_segment = false;
                    p0 = p1;
                    p1 = p2;
                    continue;
                }
            } else if (dot_prod >= ab_length2) {
                // x is projected to beyond ab
                projected_p_beyond_prev_segment = true;
                p0 = p1;
                p1 = p2;
                continue;
            } else { 
                // x is projected to a point properly on the line segment (not onto a vertex). The case which looks like | .
                projected_p_beyond_prev_segment = false;
                Point x = a + (ab.cast<double>() * (double(dot_prod) / double(ab_length2))).cast<coord_t>();
                auto dist2 = (p - x).cast<int64_t>().squaredNorm();
                if (dist2 < bestDist2) {
                    bestDist2 = dist2;
                    bestPoly = poly_idx;
                    if (distance == 0)
                        ret = x;
                    else {
                        Vec2d abd = ab.cast<double>();
                        Vec2d inward_dir = perp(abd * (distance / abd.norm())); // inward or outward depending on the sign of [distance]
                        ret = x + inward_dir.cast<coord_t>();
                        is_already_on_correct_side_of_boundary = inward_dir.dot((p - x).cast<double>()) >= 0;
                    }
                }
            }
            p0 = p1;
            p1 = p2;
        }
    }
    // when the best point is already inside and we're moving inside, or when the best point is already outside and we're moving outside
    if (is_already_on_correct_side_of_boundary) {
        if (bestDist2 < distance * distance)
            from = ret;
        else {
            // from = from; // original point stays unaltered. It is already inside by enough distance
        }
        return bestPoly;
    } else if (bestDist2 < maxDist2) {
        from = ret;
        return bestPoly;
    }
    return -1;
}

/*!
 * \brief Merges Influence Areas if possible.
 *
 * Branches which do overlap have to be merged. This helper merges all elements in input with the elements into reduced_new_layer.
 * Elements in input_aabb are merged together if possible, while elements reduced_new_layer_aabb are not checked against each other.
 *
 * \param reduced_aabb[in,out] The already processed elements.
 * \param input_aabb[in] Not yet processed elements
 * \param to_bp_areas[in] The Elements of the current Layer that will reach the buildplate. Value is the influence area where the center of a circle of support may be placed.
 * \param to_model_areas[in] The Elements of the current Layer that do not have to reach the buildplate. Also contains main as every element that can reach the buildplate is not forced to.
 * Value is the influence area where the center of a circle of support may be placed.
 * \param influence_areas[in] The influence areas without avoidance removed.
 * \param insert_bp_areas[out] Elements to be inserted into the main dictionary after the Helper terminates.
 * \param insert_model_areas[out] Elements to be inserted into the secondary dictionary after the Helper terminates.
 * \param insert_influence[out] Elements to be inserted into the dictionary containing the largest possibly valid influence area (ignoring if the area may not be there because of avoidance)
 * \param erase[out] Elements that should be deleted from the above dictionaries.
 * \param layer_idx[in] The Index of the current Layer.
 */
static void mergeHelper(
    const TreeModelVolumes &volumes, const TreeSupport::TreeSupportSettings &config,
    std::map<TreeSupport::SupportElement, BoundingBox>& reduced_aabb, std::map<TreeSupport::SupportElement, BoundingBox>& input_aabb,
    const std::unordered_map<TreeSupport::SupportElement, Polygons>& to_bp_areas, const std::unordered_map<TreeSupport::SupportElement, Polygons>& to_model_areas,
    const std::map<TreeSupport::SupportElement, Polygons>& influence_areas,
    std::unordered_map<TreeSupport::SupportElement, Polygons>& insert_bp_areas, std::unordered_map<TreeSupport::SupportElement, Polygons>& insert_model_areas,
    std::unordered_map<TreeSupport::SupportElement, Polygons>& insert_influence, std::vector<TreeSupport::SupportElement>& erase, const LayerIndex layer_idx)
{
    using SupportElement = TreeSupport::SupportElement;

    const bool first_merge_iteration = reduced_aabb.empty(); // If this is the first iteration, all elements in input have to be merged with each other
    for (std::map<SupportElement, BoundingBox>::iterator influence_iter = input_aabb.begin(); influence_iter != input_aabb.end(); influence_iter++)
    {
        bool merged = false;
        BoundingBox influence_aabb = influence_iter->second;
        for (std::map<SupportElement, BoundingBox>::iterator reduced_check_iter = reduced_aabb.begin(); reduced_check_iter != reduced_aabb.end(); reduced_check_iter++)
        {
            // As every area has to be checked for overlaps with other areas, some fast heuristic is needed to abort early if clearly possible
            // This is so performance critical that using a map lookup instead of the direct access of the cached AABBs can have a surprisingly large performance impact
            BoundingBox aabb = reduced_check_iter->second;
            if (aabb.overlap(influence_aabb)) {
                if (!first_merge_iteration && input_aabb.count(reduced_check_iter->first))
                    break; // Do not try to merge elements that already should have been merged. Done for potential performance improvement.

                bool merging_gracious_and_non_gracious = reduced_check_iter->first.to_model_gracious != influence_iter->first.to_model_gracious; // we do not want to merge a gracious with a non gracious area as bad placement could negatively impact the dependability of the whole subtree
                bool merging_to_bp = reduced_check_iter->first.to_buildplate && influence_iter->first.to_buildplate;
                bool merging_min_and_regular_xy = reduced_check_iter->first.use_min_xy_dist != influence_iter->first.use_min_xy_dist; // could cause some issues with the increase of one area, as it is assumed that if the smaller is increased by the delta to the larger it is engulfed by it already. But because a different collision may be removed from the in drawArea generated circles, this assumption could be wrong.
                coord_t increased_to_model_radius = 0;
                size_t larger_to_model_dtt = 0;

                if (!merging_to_bp) {
                    coord_t infl_radius = config.getRadius(influence_iter->first); // get the real radius increase as the user does not care for the collision model.
                    coord_t redu_radius = config.getRadius(reduced_check_iter->first);
                    if (reduced_check_iter->first.to_buildplate != influence_iter->first.to_buildplate) {
                        if (reduced_check_iter->first.to_buildplate) {
                            if (infl_radius < redu_radius)
                                increased_to_model_radius = influence_iter->first.increased_to_model_radius + redu_radius - infl_radius;
                        } else {
                            if (infl_radius > redu_radius)
                                increased_to_model_radius = reduced_check_iter->first.increased_to_model_radius + infl_radius - redu_radius;
                        }
                    }
                    larger_to_model_dtt = std::max(influence_iter->first.distance_to_top, reduced_check_iter->first.distance_to_top);
                }

                // if a merge could place a stable branch on unstable ground, would be increasing the radius further than allowed to when merging to model and to_bp trees or would merge to model before it is known they will even been drawn the merge is skipped
                if (merging_min_and_regular_xy || merging_gracious_and_non_gracious || increased_to_model_radius > config.max_to_model_radius_increase || (!merging_to_bp && larger_to_model_dtt < config.min_dtt_to_model && !reduced_check_iter->first.supports_roof && !influence_iter->first.supports_roof))
                    continue;

                Polygons relevant_infl;
                Polygons relevant_redu;
                if (merging_to_bp) {
                    relevant_infl = to_bp_areas.count(influence_iter->first) ? to_bp_areas.at(influence_iter->first) : Polygons(); // influence_iter->first is a new element => not required to check if it was changed
                    relevant_redu = insert_bp_areas.count(reduced_check_iter->first) ? insert_bp_areas[reduced_check_iter->first] : (to_bp_areas.count(reduced_check_iter->first) ? to_bp_areas.at(reduced_check_iter->first) : Polygons());
                } else {
                    relevant_infl = to_model_areas.count(influence_iter->first) ? to_model_areas.at(influence_iter->first) : Polygons();
                    relevant_redu = insert_model_areas.count(reduced_check_iter->first) ? insert_model_areas[reduced_check_iter->first] : (to_model_areas.count(reduced_check_iter->first) ? to_model_areas.at(reduced_check_iter->first) : Polygons());
                }

                const bool red_bigger = config.getCollisionRadius(reduced_check_iter->first) > config.getCollisionRadius(influence_iter->first);
                std::pair<SupportElement, Polygons> smaller_rad = red_bigger ? std::pair<SupportElement, Polygons>(influence_iter->first, relevant_infl) : std::pair<SupportElement, Polygons>(reduced_check_iter->first, relevant_redu);
                std::pair<SupportElement, Polygons> bigger_rad = red_bigger ? std::pair<SupportElement, Polygons>(reduced_check_iter->first, relevant_redu) : std::pair<SupportElement, Polygons>(influence_iter->first, relevant_infl);
                const coord_t real_radius_delta = std::abs(config.getRadius(bigger_rad.first) - config.getRadius(smaller_rad.first));
                const coord_t smaller_collision_radius = config.getCollisionRadius(smaller_rad.first);

                // the area of the bigger radius is used to ensure correct placement regarding the relevant avoidance, so if that would change an invalid area may be created
                if (!bigger_rad.first.can_use_safe_radius && smaller_rad.first.can_use_safe_radius)
                    continue;

                // the bigger radius is used to verify that the area is still valid after the increase with the delta. If there were a point where the big influence area could be valid with can_use_safe_radius the element would already be can_use_safe_radius
                // the smaller radius, which gets increased by delta may reach into the area where use_min_xy_dist is no longer required.
                bool use_min_radius = bigger_rad.first.use_min_xy_dist && smaller_rad.first.use_min_xy_dist;

                // The idea is that the influence area with the smaller collision radius is increased by the radius difference.
                // If this area has any intersections with the influence area of the larger collision radius, a branch (of the larger collision radius) placed in this intersection, has already engulfed the branch of the smaller collision radius.
                // Because of this a merge may happen even if the influence areas (that represent possible center points of branches) do not intersect yet.
                // Remember that collision radius <= real radius as otherwise this assumption would be false.
                Polygons small_rad_increased_by_big_minus_small = safeOffsetInc(smaller_rad.second, real_radius_delta, volumes.getCollision(smaller_collision_radius, layer_idx - 1, use_min_radius), 2 * (config.xy_distance + smaller_collision_radius - 3), 0, 0); // -3 avoids possible rounding errors
                Polygons intersect = intersection(small_rad_increased_by_big_minus_small, bigger_rad.second);

                if (area(intersect) > tiny_area_threshold) { // dont use empty as a line is not empty, but for this use-case it very well may be (and would be one layer down as union does not keep lines)
                    if (area(offset(intersect, scaled<float>(-0.025), jtMiter, 1.2)) <= tiny_area_threshold) // check if the overlap is large enough (Small ares tend to attract rounding errors in clipper). While 25 was guessed as enough, i did not have reason to change it.
                        continue;

                    // Do the actual merge now that the branches are confirmed to be able to intersect.

                    // calculate which point is closest to the point of the last merge (or tip center if no merge above it has happened)
                    // used at the end to estimate where to best place the branch on the bottom most layer
                    // could be replaced with a random point inside the new area
                    Point new_pos = reduced_check_iter->first.next_position;
                    if (! contains(intersect, new_pos))
                        moveInside(intersect, new_pos);

                    if (increased_to_model_radius == 0)
                        increased_to_model_radius = std::max(reduced_check_iter->first.increased_to_model_radius, influence_iter->first.increased_to_model_radius);

                    SupportElement key(reduced_check_iter->first, influence_iter->first, layer_idx - 1, new_pos, increased_to_model_radius, config);

                    Polygons intersect_influence;
                    Polygons infl_small = insert_influence.count(smaller_rad.first) ? insert_influence[smaller_rad.first] : (influence_areas.count(smaller_rad.first) ? influence_areas.at(smaller_rad.first) : Polygons());
                    Polygons infl_big = insert_influence.count(bigger_rad.first) ? insert_influence[bigger_rad.first] : (influence_areas.count(bigger_rad.first) ? influence_areas.at(bigger_rad.first) : Polygons());
                    Polygons small_rad_increased_by_big_minus_small_infl = safeOffsetInc(infl_small, real_radius_delta, volumes.getCollision(smaller_collision_radius, layer_idx - 1, use_min_radius), 2 * (config.xy_distance + smaller_collision_radius - 3), 0, 0);
                    intersect_influence = intersection(small_rad_increased_by_big_minus_small_infl, infl_big); // if the one with the bigger radius with the lower radius removed overlaps we can merge
                    intersect_influence = safeUnion(intersect_influence, intersect); // Rounding errors again. Do not ask me where or why.

                    Polygons intersect_sec;
                    if (merging_to_bp && config.support_rests_on_model) {
                        if (key.to_model_gracious) {
                            Polygons sec_small = insert_model_areas.count(smaller_rad.first) ? insert_model_areas[smaller_rad.first] : (to_model_areas.count(smaller_rad.first) ? to_model_areas.at(smaller_rad.first) : Polygons());
                            Polygons sec_big = insert_model_areas.count(bigger_rad.first) ? insert_model_areas[bigger_rad.first] : (to_model_areas.count(bigger_rad.first) ? to_model_areas.at(bigger_rad.first) : Polygons());
                            Polygons small_rad_increased_by_big_minus_small_sec = safeOffsetInc(sec_small, real_radius_delta, volumes.getCollision(smaller_collision_radius, layer_idx - 1, use_min_radius), 2 * (config.xy_distance + smaller_collision_radius - 3), 0, 0);
                            intersect_sec = intersection(small_rad_increased_by_big_minus_small_sec, sec_big); // if the one with the bigger radius with the lower radius removed overlaps we can merge
                            intersect_influence = safeUnion(intersect_influence, intersect_sec); // still rounding errors
                        } else
                            intersect_sec = intersect_influence;
                    }

                    // remove the now merged elements from all buckets, as they do not exist anymore in their old form
                    insert_bp_areas.erase(reduced_check_iter->first);
                    insert_bp_areas.erase(influence_iter->first);
                    insert_model_areas.erase(reduced_check_iter->first);
                    insert_model_areas.erase(influence_iter->first);
                    insert_influence.erase(reduced_check_iter->first);
                    insert_influence.erase(influence_iter->first);

                    (merging_to_bp ? insert_bp_areas : insert_model_areas).emplace(key, intersect);
                    if (merging_to_bp && config.support_rests_on_model)
                        insert_model_areas.emplace(key, intersect_sec);
                    insert_influence.emplace(key, intersect_influence);

                    erase.emplace_back(reduced_check_iter->first);
                    erase.emplace_back(influence_iter->first);
                    Polygons merge = diff_clipped(offset(union_(intersect, intersect_sec), config.getRadius(key), ClipperLib::jtRound, scaled<float>(0.01)), volumes.getCollision(0, layer_idx - 1, false)); // regular union should be preferable here as Polygons tend to only become smaller through rounding errors (smaller!=has smaller area as holes have a negative area.). And if this area disappears because of rounding errors, the only downside is that it can not merge again on this layer.

                    reduced_aabb.erase(reduced_check_iter->first); // this invalidates reduced_check_iter
                    reduced_aabb.emplace(key, get_extents(merge));

                    merged = true;
                    break;
                }
            }
        }

        if (!merged)
            reduced_aabb[influence_iter->first] = influence_aabb;
    }
}

/*!
 * \brief Merges Influence Areas if possible.
 *
 * Branches which do overlap have to be merged. This manages the helper and uses a divide and conquer approach to parallelize this problem. This parallelization can at most accelerate the merging by a factor of 2.
 *
 * \param to_bp_areas[in] The Elements of the current Layer that will reach the buildplate.
 *  Value is the influence area where the center of a circle of support may be placed.
 * \param to_model_areas[in] The Elements of the current Layer that do not have to reach the buildplate. Also contains main as every element that can reach the buildplate is not forced to.
 *  Value is the influence area where the center of a circle of support may be placed.
 * \param influence_areas[in] The Elements of the current Layer without avoidances removed. This is the largest possible influence area for this layer.
 *  Value is the influence area where the center of a circle of support may be placed.
 * \param layer_idx[in] The current layer.
 */
static void mergeInfluenceAreas(
    const TreeModelVolumes &volumes, const TreeSupport::TreeSupportSettings &config,
    std::unordered_map<TreeSupport::SupportElement, Polygons>& to_bp_areas, std::unordered_map<TreeSupport::SupportElement, Polygons>& to_model_areas, std::map<TreeSupport::SupportElement, Polygons>& influence_areas, LayerIndex layer_idx)
{
    using SupportElement = TreeSupport::SupportElement;
    /*
     * Idea behind this is that the calculation of merges can be accelerated a bit using divide and conquer:
     * If two groups of areas are already merged, only all elements in group 2 have to be merged into group one.
     * This can only accelerate by factor 2 (as half the work is merging the last two groups).
     * The actual merge logic is found in mergeHelper. This function only manages parallelization of different mergeHelper calls.
     */

    const size_t input_size = influence_areas.size();
    if (input_size == 0)
        return;

    size_t num_threads = std::max(size_t(1), size_t(std::thread::hardware_concurrency())); // For some reason hardware concurrency can return 0;
    constexpr int min_elements_per_bucket = 2;

    // max_bucket_count is input_size/min_elements_per_bucket round down to the next 2^n.
    // The rounding to 2^n is to ensure improved performance, as every iteration two buckets will be merged, halving the amount of buckets.
    // If halving would cause an uneven count, e.g. 3 Then bucket 0 and 1 would have to be merged, and in the next iteration the last remaining buckets. This is assumed to not be optimal performance-wise.
    const size_t max_bucket_count = std::pow(2, std::floor(std::log(round_up_divide<int>(input_size, min_elements_per_bucket))));
    int bucket_count = std::min(max_bucket_count, num_threads); // do not use more buckets than available threads.

    // To achieve that every element in a bucket is already correctly merged with other elements in this bucket
    // an extra empty bucket is created for each bucket, and the elements are merged into the empty one.
    // Each thread will then process two buckets by merging all elements in the second bucket into the first one as mergeHelper will disable not trying to merge elements from the same bucket in this case.
    std::vector<std::map<SupportElement, Polygons>> buckets_area(2 * bucket_count);
    std::vector<std::map<SupportElement, BoundingBox>> buckets_aabb(2 * bucket_count);

    size_t position = 0, counter = 0;
    const size_t over_elements = input_size % bucket_count;
    const size_t elements_per_step = input_size / bucket_count;

    // split the data in x parts to be able to divide and conquer
    // the first "over_elements" of buckets gets elements_per_step+1 elements
    for (std::map<SupportElement, Polygons>::iterator iter = influence_areas.begin(); iter != influence_areas.end(); ++ iter) {
        buckets_area[position * 2 + 1].emplace(iter->first, iter->second); // only use every second bucket beginning with 1 as this makes the parallel call later easier as we assume everything in a bucket i%2==0 is already processed
        ++ counter;
        if ((counter == elements_per_step && position >= over_elements) || counter > elements_per_step) {
            position++;
            counter = 0;
        }
    }

    // precalculate the AABBs from the influence areas.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, bucket_count),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t idx = range.begin(); idx < range.end(); ++ idx) {
            // +=2 as in the beginning only uneven buckets will be filled
            size_t bucket_idx = 2 * idx + 1;
            for (const std::pair<const SupportElement, Polygons>& input_pair : buckets_area[bucket_idx])
                buckets_aabb[bucket_idx].emplace(input_pair.first, get_extents(input_pair.second).inflated(config.getRadius(input_pair.first)));
        }
    });

    while (buckets_area.size() > 1) {
        // Some temporary storage, of elements that have to be inserted or removed from the background storage. Only one per two buckets required
        std::vector<std::unordered_map<SupportElement, Polygons>> insert_main(buckets_area.size() / 2);
        std::vector<std::unordered_map<SupportElement, Polygons>> insert_secondary(buckets_area.size() / 2);
        std::vector<std::unordered_map<SupportElement, Polygons>> insert_influence(buckets_area.size() / 2);
        std::vector<std::vector<SupportElement>> erase(buckets_area.size() / 2);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, buckets_area.size() / 2),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t idx = range.begin(); idx < range.end(); ++ idx) {
                const size_t bucket_pair_idx = idx * 2;
                // Merge bucket_count adjacent to each other, merging uneven bucket numbers into even buckets
                mergeHelper(volumes, config, buckets_aabb[bucket_pair_idx], buckets_aabb[bucket_pair_idx + 1], to_bp_areas, to_model_areas, influence_areas, insert_main[bucket_pair_idx / 2], insert_secondary[bucket_pair_idx / 2], insert_influence[bucket_pair_idx / 2], erase[bucket_pair_idx / 2], layer_idx);
                // clear now irrelevant max_bucket_count, and delete them later
                buckets_area[bucket_pair_idx + 1].clear();
                buckets_aabb[bucket_pair_idx + 1].clear();
            }
        });

        for (size_t i = 0; i + 1 < buckets_area.size(); i += 2) {
            for (SupportElement &del : erase[i / 2]) {
                to_bp_areas.erase(del);
                to_model_areas.erase(del);
                influence_areas.erase(del);
            }
            for (const std::pair<const SupportElement, Polygons> &tup : insert_main[i / 2])
                to_bp_areas.emplace(std::move(tup));
            for (const std::pair<const SupportElement, Polygons> &tup : insert_secondary[i / 2])
                to_model_areas.emplace(std::move(tup));
            for (const std::pair<const SupportElement, Polygons> &tup : insert_influence[i / 2])
                influence_areas.emplace(std::move(tup));
        }

        buckets_area.erase(std::remove_if(buckets_area.begin(), buckets_area.end(), [&](const std::map<SupportElement, Polygons> &x) { return x.empty(); }), buckets_area.end());
        buckets_aabb.erase(std::remove_if(buckets_aabb.begin(), buckets_aabb.end(), [&](const std::map<SupportElement, BoundingBox> &x) { return x.empty(); }), buckets_aabb.end());
    }
}


std::optional<TreeSupport::SupportElement> TreeSupport::increaseSingleArea(AreaIncreaseSettings settings, LayerIndex layer_idx, SupportElement* parent, const Polygons& relevant_offset, Polygons& to_bp_data, Polygons& to_model_data, Polygons& increased, const coord_t overspeed, const bool mergelayer)
{
    SupportElement current_elem(parent); // also increases DTT by one
    Polygons check_layer_data;
    if (settings.increase_radius)
        current_elem.effective_radius_height += 1;
    coord_t radius = m_config.getCollisionRadius(current_elem);

    if (settings.move) {
        increased = relevant_offset;
        if (overspeed > 0) {
            const coord_t safe_movement_distance = (current_elem.use_min_xy_dist ? m_config.xy_min_distance : m_config.xy_distance) + (std::min(m_config.z_distance_top_layers, m_config.z_distance_bottom_layers) > 0 ? m_config.min_feature_size : 0);
            // The difference to ensure that the result not only conforms to wall_restriction, but collision/avoidance is done later. The higher last_safe_step_movement_distance comes exactly from the fact that the collision will be subtracted later.
            increased = safeOffsetInc(increased, overspeed, m_volumes.getWallRestriction(m_config.getCollisionRadius(*parent), layer_idx, parent->use_min_xy_dist), safe_movement_distance, safe_movement_distance + radius, 1);
        }
        if (settings.no_error && settings.move)
            // as ClipperLib::jtRound has to be used for offsets this simplify is VERY important for performance.
            polygons_simplify(increased, scaled<float>(0.025));
    } else 
        // if no movement is done the areas keep parent area as no move == offset(0)
        increased = *parent->area;

    if (mergelayer || current_elem.to_buildplate) {
        to_bp_data = safeUnion(diff_clipped(increased, m_volumes.getAvoidance(radius, layer_idx - 1, settings.type, false, settings.use_min_distance)));
        if (! current_elem.to_buildplate && area(to_bp_data) > tiny_area_threshold) {
            // mostly happening in the tip, but with merges one should check every time, just to be sure.
            current_elem.to_buildplate = true; // sometimes nodes that can reach the buildplate are marked as cant reach, tainting subtrees. This corrects it.
            BOOST_LOG_TRIVIAL(debug) << "Corrected taint leading to a wrong to model value on layer " << layer_idx - 1 << " targeting " << current_elem.target_height << " with radius " << radius;
        }
    }
    if (m_config.support_rests_on_model) {
        if (mergelayer || current_elem.to_model_gracious)
            to_model_data = safeUnion(diff_clipped(increased, m_volumes.getAvoidance(radius, layer_idx - 1, settings.type, true, settings.use_min_distance)));

        if (!current_elem.to_model_gracious) {
            if (mergelayer && area(to_model_data) >= tiny_area_threshold) {
                current_elem.to_model_gracious = true;
                BOOST_LOG_TRIVIAL(debug) << "Corrected taint leading to a wrong non gracious value on layer " << layer_idx - 1 << " targeting " << current_elem.target_height << " with radius " << radius;
            } else
                to_model_data = safeUnion(diff_clipped(increased, m_volumes.getCollision(radius, layer_idx - 1, settings.use_min_distance)));
        }
    }

    check_layer_data = current_elem.to_buildplate ? to_bp_data : to_model_data;

    if (settings.increase_radius && area(check_layer_data) > tiny_area_threshold) {
        auto validWithRadius = [&](coord_t next_radius) {
            if (m_volumes.ceilRadius(next_radius, settings.use_min_distance) <= m_volumes.ceilRadius(radius, settings.use_min_distance))
                return true;

            Polygons to_bp_data_2;
            if (current_elem.to_buildplate)
                to_bp_data_2 = diff_clipped(increased, m_volumes.getAvoidance(next_radius, layer_idx - 1, settings.type, false, settings.use_min_distance)); // regular union as output will not be used later => this area should always be a subset of the safeUnion one (i think)
            Polygons to_model_data_2;
            if (m_config.support_rests_on_model && !current_elem.to_buildplate)
                to_model_data_2 = diff_clipped(increased, 
                    current_elem.to_model_gracious ? 
                        m_volumes.getAvoidance(next_radius, layer_idx - 1, settings.type, true, settings.use_min_distance) :
                        m_volumes.getCollision(next_radius, layer_idx - 1, settings.use_min_distance));
            Polygons check_layer_data_2 = current_elem.to_buildplate ? to_bp_data_2 : to_model_data_2;
            return area(check_layer_data_2) > tiny_area_threshold;
        };
        coord_t ceil_radius_before = m_volumes.ceilRadius(radius, settings.use_min_distance);

        if (m_config.getCollisionRadius(current_elem) < m_config.increase_radius_until_radius && m_config.getCollisionRadius(current_elem) < m_config.getRadius(current_elem)) {
            coord_t target_radius = std::min(m_config.getRadius(current_elem), m_config.increase_radius_until_radius);
            coord_t current_ceil_radius = m_volumes.getRadiusNextCeil(radius, settings.use_min_distance);

            while (current_ceil_radius < target_radius && validWithRadius(m_volumes.getRadiusNextCeil(current_ceil_radius + 1, settings.use_min_distance)))
                current_ceil_radius = m_volumes.getRadiusNextCeil(current_ceil_radius + 1, settings.use_min_distance);
            size_t resulting_eff_dtt = current_elem.effective_radius_height;
            while (resulting_eff_dtt + 1 < current_elem.distance_to_top && m_config.getRadius(resulting_eff_dtt + 1, current_elem.elephant_foot_increases) <= current_ceil_radius && m_config.getRadius(resulting_eff_dtt + 1, current_elem.elephant_foot_increases) <= m_config.getRadius(current_elem))
                resulting_eff_dtt++;
            current_elem.effective_radius_height = resulting_eff_dtt;
        }
        radius = m_config.getCollisionRadius(current_elem);

        const coord_t foot_radius_increase = m_config.branch_radius * (std::max(m_config.diameter_scale_bp_radius - m_config.diameter_angle_scale_factor, 0.0));
        // Is nearly all of the time 1, but sometimes an increase of 1 could cause the radius to become bigger than recommendedMinRadius, which could cause the radius to become bigger than precalculated.
        double planned_foot_increase = std::min(1.0, double(m_config.recommendedMinRadius(layer_idx - 1) - m_config.getRadius(current_elem)) / foot_radius_increase);
//FIXME
        bool increase_bp_foot = planned_foot_increase > 0 && current_elem.to_buildplate;
//        bool increase_bp_foot = false;

        if (increase_bp_foot && m_config.getRadius(current_elem) >= m_config.branch_radius && m_config.getRadius(current_elem) >= m_config.increase_radius_until_radius)
            if (validWithRadius(m_config.getRadius(current_elem.effective_radius_height, current_elem.elephant_foot_increases + planned_foot_increase))) {
                current_elem.elephant_foot_increases += planned_foot_increase;
                radius = m_config.getCollisionRadius(current_elem);
            }

        if (ceil_radius_before != m_volumes.ceilRadius(radius, settings.use_min_distance)) {
            if (current_elem.to_buildplate)
                to_bp_data = safeUnion(diff_clipped(increased, m_volumes.getAvoidance(radius, layer_idx - 1, settings.type, false, settings.use_min_distance)));
            if (m_config.support_rests_on_model && (!current_elem.to_buildplate || mergelayer))
                to_model_data = safeUnion(diff_clipped(increased, 
                    current_elem.to_model_gracious ? 
                        m_volumes.getAvoidance(radius, layer_idx - 1, settings.type, true, settings.use_min_distance) :
                        m_volumes.getCollision(radius, layer_idx - 1, settings.use_min_distance)
                ));
            check_layer_data = current_elem.to_buildplate ? to_bp_data : to_model_data;
            if (area(check_layer_data) < tiny_area_threshold) {
                BOOST_LOG_TRIVIAL(error) << "Lost area by doing catch up from " << ceil_radius_before << " to radius " << m_volumes.ceilRadius(m_config.getCollisionRadius(current_elem), settings.use_min_distance);
                TreeSupport::showError("Area lost catching up radius. May not cause visible malformation.", true);
            }
        }
    }

    return area(check_layer_data) > tiny_area_threshold ? std::optional<SupportElement>(current_elem) : std::optional<SupportElement>();
}

void TreeSupport::increaseAreas(std::unordered_map<SupportElement, Polygons>& to_bp_areas, std::unordered_map<SupportElement, Polygons>& to_model_areas, std::map<SupportElement, Polygons>& influence_areas, std::vector<SupportElement*>& bypass_merge_areas, const std::vector<SupportElement*>& last_layer, const LayerIndex layer_idx, const bool mergelayer)
{
    std::mutex critical_sections;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, last_layer.size()),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t idx = range.begin(); idx < range.end(); ++ idx) {
            SupportElement* parent = last_layer[idx];

            SupportElement elem(parent); // also increases dtt

            const Polygons &wall_restriction = m_volumes.getWallRestriction(m_config.getCollisionRadius(*parent), layer_idx, parent->use_min_xy_dist); // Abstract representation of the model outline. If an influence area would move through it, it could teleport through a wall.

            Polygons to_bp_data, to_model_data;
            coord_t radius = m_config.getCollisionRadius(elem);

            // When the radius increases, the outer "support wall" of the branch will have been moved farther away from the center (as this is the definition of radius).
            // As it is not specified that the support_tree_angle has to be one of the center of the branch, it is here seen as the smaller angle of the outer wall of the branch, to the outer wall of the same branch one layer above.
            // As the branch may have become larger the distance between these 2 walls is smaller than the distance of the center points.
            // These extra distance is added to the movement distance possible for this layer.

            coord_t extra_speed = 5; // The extra speed is added to both movement distances. Also move 5 microns faster than allowed to avoid rounding errors, this may cause issues at VERY VERY small layer heights.
            coord_t extra_slow_speed = 0; // Only added to the slow movement distance.
            const coord_t ceiled_parent_radius = m_volumes.ceilRadius(m_config.getCollisionRadius(*parent), parent->use_min_xy_dist);
            coord_t projected_radius_increased = m_config.getRadius(parent->effective_radius_height + 1, parent->elephant_foot_increases);
            coord_t projected_radius_delta = projected_radius_increased - m_config.getCollisionRadius(*parent);

            // When z distance is more than one layer up and down the Collision used to calculate the wall restriction will always include the wall (and not just the xy_min_distance) of the layer above and below like this (d = blocked area because of z distance):
            /*
             *  layer z+1:dddddiiiiiioooo
             *  layer z+0:xxxxxdddddddddd
             *  layer z-1:dddddxxxxxxxxxx
             *  For more detailed visualisation see calculateWallRestrictions
             */
            const coord_t safe_movement_distance = (elem.use_min_xy_dist ? m_config.xy_min_distance : m_config.xy_distance) + (std::min(m_config.z_distance_top_layers, m_config.z_distance_bottom_layers) > 0 ? m_config.min_feature_size : 0);
            if (ceiled_parent_radius == m_volumes.ceilRadius(projected_radius_increased, parent->use_min_xy_dist) || projected_radius_increased < m_config.increase_radius_until_radius)
                // If it is guaranteed possible to increase the radius, the maximum movement speed can be increased, as it is assumed that the maximum movement speed is the one of the slower moving wall
                extra_speed += projected_radius_delta;
            else
                // if a guaranteed radius increase is not possible, only increase the slow speed
                extra_slow_speed += std::min(projected_radius_delta, (m_config.maximum_move_distance + extra_speed) - (m_config.maximum_move_distance_slow + extra_slow_speed)); // Ensure that the slow movement distance can not become larger than the fast one.

            if (m_config.layer_start_bp_radius > layer_idx && m_config.recommendedMinRadius(layer_idx - 1) < m_config.getRadius(elem.effective_radius_height + 1, elem.elephant_foot_increases)) {
                // can guarantee elephant foot radius increase
                if (ceiled_parent_radius == m_volumes.ceilRadius(m_config.getRadius(parent->effective_radius_height + 1, parent->elephant_foot_increases + 1), parent->use_min_xy_dist))
                    extra_speed += m_config.branch_radius * m_config.diameter_scale_bp_radius;
                else
                    extra_slow_speed += std::min(coord_t(m_config.branch_radius * m_config.diameter_scale_bp_radius), m_config.maximum_move_distance - (m_config.maximum_move_distance_slow + extra_slow_speed));
            }

            const coord_t fast_speed = m_config.maximum_move_distance + extra_speed;
            const coord_t slow_speed = m_config.maximum_move_distance_slow + extra_speed + extra_slow_speed;

            Polygons offset_slow, offset_fast;

            bool add = false;
            bool bypass_merge = false;
            constexpr bool increase_radius = true, no_error = true, use_min_radius = true, move = true; // aliases for better readability

            // Determine in which order configurations are checked if they result in a valid influence area. Check will stop if a valid area is found
            std::vector<AreaIncreaseSettings> order;
            auto insertSetting = [&](AreaIncreaseSettings settings, bool back) {
                if (std::find(order.begin(), order.end(), settings) == order.end()) {
                    if (back)
                        order.emplace_back(settings);
                    else
                        order.insert(order.begin(), settings);
                }
            };

            const bool parent_moved_slow = elem.last_area_increase.increase_speed < m_config.maximum_move_distance;
            const bool avoidance_speed_mismatch = parent_moved_slow && elem.last_area_increase.type != AvoidanceType::Slow;
            if (elem.last_area_increase.move && elem.last_area_increase.no_error && elem.can_use_safe_radius && !mergelayer && !avoidance_speed_mismatch && (elem.distance_to_top >= m_config.tip_layers || parent_moved_slow))
            {
                // assume that the avoidance type that was best for the parent is best for me. Makes this function about 7% faster.
                insertSetting({ elem.last_area_increase.type, elem.last_area_increase.increase_speed < m_config.maximum_move_distance ? slow_speed : fast_speed, increase_radius, elem.last_area_increase.no_error, !use_min_radius, elem.last_area_increase.move }, true);
                insertSetting({ elem.last_area_increase.type, elem.last_area_increase.increase_speed < m_config.maximum_move_distance ? slow_speed : fast_speed, !increase_radius, elem.last_area_increase.no_error, !use_min_radius, elem.last_area_increase.move }, true);
            }
            // branch may still go though a hole, so a check has to be done whether the hole was already passed, and the regular avoidance can be used.
            if (!elem.can_use_safe_radius)
            {
                // if the radius until which it is always increased can not be guaranteed, move fast. This is to avoid holes smaller than the real branch radius. This does not guarantee the avoidance of such holes, but ensures they are avoided if possible.
                // order.emplace_back(AvoidanceType::Slow,!increase_radius,no_error,!use_min_radius,move);
                insertSetting({ AvoidanceType::Slow, slow_speed, increase_radius, no_error, !use_min_radius, !move }, true); // did we go through the hole
                // in many cases the definition of hole is overly restrictive, so to avoid unnecessary fast movement in the tip, it is ignored there for a bit. This CAN cause a branch to go though a hole it otherwise may have avoided.
                if (elem.distance_to_top < round_up_divide(m_config.tip_layers, size_t(2)))
                    insertSetting({ AvoidanceType::Fast, slow_speed, increase_radius, no_error, !use_min_radius, !move }, true);
                insertSetting({ AvoidanceType::FastSafe, fast_speed, increase_radius, no_error, !use_min_radius, !move }, true); // did we manage to avoid the hole
                insertSetting({ AvoidanceType::FastSafe, fast_speed, !increase_radius, no_error, !use_min_radius, move }, true);
                insertSetting({ AvoidanceType::Fast, fast_speed, !increase_radius, no_error, !use_min_radius, move }, true);
            }
            else
            {
                insertSetting({ AvoidanceType::Slow, slow_speed, increase_radius, no_error, !use_min_radius, move }, true);
                // while moving fast to be able to increase the radius (b) may seems preferable (over a) this can cause the a sudden skip in movement, which looks similar to a layer shift and can reduce stability.
                // as such idx have chosen to only use the user setting for radius increases as a friendly recommendation.
                insertSetting({ AvoidanceType::Slow, slow_speed, !increase_radius, no_error, !use_min_radius, move }, true); // a
                if (elem.distance_to_top < m_config.tip_layers)
                {
                    insertSetting({ AvoidanceType::FastSafe, slow_speed, increase_radius, no_error, !use_min_radius, move }, true);
                }
                insertSetting({ AvoidanceType::FastSafe, fast_speed, increase_radius, no_error, !use_min_radius, move }, true); // b
                insertSetting({ AvoidanceType::FastSafe, fast_speed, !increase_radius, no_error, !use_min_radius, move }, true);
            }

            if (elem.use_min_xy_dist)
            {
                std::vector<AreaIncreaseSettings> new_order;
                // if the branch currently has to use min_xy_dist check if the configuration would also be valid with the regular xy_distance before checking with use_min_radius (Only happens when Support Distance priority is z overrides xy )
                for (AreaIncreaseSettings settings : order)
                {
                    new_order.emplace_back(settings);
                    new_order.push_back({ settings.type, settings.increase_speed, settings.increase_radius, settings.no_error, use_min_radius, settings.move });
                }
                order = new_order;
            }
            if (elem.to_buildplate || (elem.to_model_gracious && intersection(*parent->area, m_volumes.getPlaceableAreas(radius, layer_idx)).empty())) // error case
            {
                // it is normal that we wont be able to find a new area at some point in time if we wont be able to reach layer 0 aka have to connect with the model
                insertSetting({ AvoidanceType::Fast, fast_speed, !increase_radius, !no_error, elem.use_min_xy_dist, move }, true);
            }
            if (elem.distance_to_top < elem.dont_move_until && elem.can_use_safe_radius) // only do not move when holes would be avoided in every case.
                // Only do not move when already in a no hole avoidance with the regular xy distance.
                insertSetting({ AvoidanceType::Slow, 0, increase_radius, no_error, !use_min_radius, !move }, false);

            Polygons inc_wo_collision;
            // Check whether it is faster to calculate the area increased with the fast speed independently from the slow area, or time could be saved by reusing the slow area to calculate the fast one.
            // Calculated by comparing the steps saved when calcualting idependently with the saved steps when not.
            bool offset_independant_faster = (radius / safe_movement_distance - (((m_config.maximum_move_distance + extra_speed) < (radius + safe_movement_distance)) ? 1 : 0)) > (round_up_divide((extra_speed + extra_slow_speed + m_config.maximum_move_distance_slow), safe_movement_distance));
            for (AreaIncreaseSettings settings : order)
            {
                if (settings.move) {
                    if (offset_slow.empty() && (settings.increase_speed == slow_speed || !offset_independant_faster)) {
                        // offsetting in 2 steps makes our offsetted area rounder preventing (rounding) errors created by to pointy areas. At this point one can see that the Polygons class 
                        // was never made for precision in the single digit micron range.
                        offset_slow = safeOffsetInc(*parent->area, extra_speed + extra_slow_speed + m_config.maximum_move_distance_slow, wall_restriction, safe_movement_distance, offset_independant_faster ? safe_movement_distance + radius : 0, 2); 
                    }

                    if ((settings.increase_speed != slow_speed) && offset_fast.empty()) {
                        if (offset_independant_faster)
                            offset_fast = safeOffsetInc(*parent->area, extra_speed + m_config.maximum_move_distance, wall_restriction, safe_movement_distance, offset_independant_faster ? safe_movement_distance + radius : 0, 1);
                        else {
                            const coord_t delta_slow_fast = m_config.maximum_move_distance - (m_config.maximum_move_distance_slow + extra_slow_speed);
                            offset_fast = safeOffsetInc(offset_slow, delta_slow_fast, wall_restriction, safe_movement_distance, safe_movement_distance + radius, offset_independant_faster ? 2 : 1);
                        }
                    }
                }
                std::optional<SupportElement> result;
                if (!settings.no_error) { 
                    // ERROR CASE
                    // if the area becomes for whatever reason something that clipper sees as a line, offset would stop working, so ensure that even if if wrongly would be a line, it still actually has an area that can be increased
                    Polygons lines_offset = offset(to_polylines(*parent->area), scaled<float>(0.005), jtMiter, 1.2);
                    Polygons base_error_area = union_(*parent->area, lines_offset);
                    result = increaseSingleArea(settings, layer_idx, parent, base_error_area, to_bp_data, to_model_data, inc_wo_collision, (m_config.maximum_move_distance + extra_speed) * 1.5, mergelayer);
#ifdef TREE_SUPPORT_SHOW_ERRORS
                    BOOST_LOG_TRIVIAL(error)
#else // TREE_SUPPORT_SHOW_ERRORS
                    BOOST_LOG_TRIVIAL(warning)
#endif // TREE_SUPPORT_SHOW_ERRORS
                          << "Influence area could not be increased! Data about the Influence area: "
                             "Radius: " << radius << " at layer: " << layer_idx - 1 << " NextTarget: " << elem.next_height << " Distance to top: " << elem.distance_to_top << 
                             " Elephant foot increases " << elem.elephant_foot_increases << " use_min_xy_dist " << elem.use_min_xy_dist << " to buildplate " << elem.to_buildplate << 
                             " gracious " << elem.to_model_gracious << " safe " << elem.can_use_safe_radius << " until move " << elem.dont_move_until << " \n "
                             "Parent " << parent << ": Radius: " << m_config.getCollisionRadius(*parent) << " at layer: " << layer_idx << " NextTarget: " << parent->next_height << 
                             " Distance to top: " << parent->distance_to_top << " Elephant foot increases " << parent->elephant_foot_increases << "  use_min_xy_dist " << parent->use_min_xy_dist << 
                             " to buildplate " << parent->to_buildplate << " gracious " << parent->to_model_gracious << " safe " << parent->can_use_safe_radius << " until move " << parent->dont_move_until;
                    showError("Potentially lost branch!", true);
                } else
                    result = increaseSingleArea(settings, layer_idx, parent, settings.increase_speed == slow_speed ? offset_slow : offset_fast, to_bp_data, to_model_data, inc_wo_collision, 0, mergelayer);

                if (result)
                {
                    elem = *result;
                    radius = m_config.getCollisionRadius(elem);
                    elem.last_area_increase = settings;
                    add = true;
                    bypass_merge = !settings.move || (settings.use_min_distance && elem.distance_to_top < m_config.tip_layers); // do not merge if the branch should not move or the priority has to be to get farther away from the model.
                    if (settings.move)
                        elem.dont_move_until = 0;
                    else
                        elem.result_on_layer = parent->result_on_layer;

                    elem.can_use_safe_radius = settings.type != AvoidanceType::Fast;

                    if (!settings.use_min_distance)
                        elem.use_min_xy_dist = false;
                    if (!settings.no_error)
#ifdef TREE_SUPPORT_SHOW_ERRORS
                        BOOST_LOG_TRIVIAL(error) 
#else // TREE_SUPPORT_SHOW_ERRORS
                        BOOST_LOG_TRIVIAL(info)
#endif // TREE_SUPPORT_SHOW_ERRORS
                            << "Trying to keep area by moving faster than intended: Success";
                    break;
                }
                else if (!settings.no_error)
                    BOOST_LOG_TRIVIAL(error) << "Trying to keep area by moving faster than intended: FAILURE! WRONG BRANCHES LIKLY!";
            }

            if (add) {
                Polygons max_influence_area = safeUnion(diff_clipped(inc_wo_collision, m_volumes.getCollision(radius, layer_idx - 1, elem.use_min_xy_dist)), safeUnion(to_bp_data, to_model_data)); // union seems useless, but some rounding errors somewhere can cause to_bp_data to be slightly bigger than it should be
                {
                    std::lock_guard<std::mutex> critical_section_newLayer(critical_sections);
                    if (bypass_merge) {
                        validate_range(max_influence_area);
                        Polygons* new_area = new Polygons(max_influence_area);
                        SupportElement* next = new SupportElement(elem, new_area);
                        bypass_merge_areas.emplace_back(next);
                    } else {
                        influence_areas.emplace(elem, max_influence_area);
                        if (elem.to_buildplate)
                            to_bp_areas.emplace(elem, to_bp_data);
                        if (m_config.support_rests_on_model)
                            to_model_areas.emplace(elem, to_model_data);
                    }
                }
            }
            else
                parent->result_on_layer = Point(-1, -1); // If the bottom most point of a branch is set, later functions will assume that the position is valid, and ignore it. But as branches connecting with the model that are to small have to be culled, the bottom most point has to be not set. A point can be set on the top most tip layer (maybe more if it should not move for a few layers).
        }
    });
}


void TreeSupport::createLayerPathing(std::vector<std::set<SupportElement*>>& move_bounds)
{
#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    const double data_size_inverse = 1 / double(move_bounds.size());
    double progress_total = TREE_PROGRESS_PRECALC_AVO + TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_GENERATE_NODES;
#endif // SLIC3R_TREESUPPORTS_PROGRESS

    auto dur_inc = std::chrono::duration_values<std::chrono::nanoseconds>::zero();
    auto dur_merge = std::chrono::duration_values<std::chrono::nanoseconds>::zero();

    LayerIndex last_merge = move_bounds.size();
    bool new_element = false;

    size_t max_merge_every_x_layers = std::min(std::min(5000 / (std::max(m_config.maximum_move_distance, coord_t(100))), 1000 / std::max(m_config.maximum_move_distance_slow, coord_t(20))), 3000 / m_config.layer_height); // Ensures at least one merge operation per 3mm height, 50 layers, 1 mm movement of slow speed or 5mm movement of fast speed (whatever is lowest). Values were guessed.
    size_t merge_every_x_layers = 1;
    // Calculate the influence areas for each layer below (Top down)
    // This is done by first increasing the influence area by the allowed movement distance, and merging them with other influence areas if possible
    for (int layer_idx = int(move_bounds.size()) - 1; layer_idx > 0; -- layer_idx)
    {
        // merging is expensive and only parallelized to a max speedup of 2. As such it may be useful in some cases to only merge every few layers to improve performance.
        bool merge_this_layer = size_t(last_merge - layer_idx) >= merge_every_x_layers;

        if (new_element)
        {
            merge_this_layer = true;
            merge_every_x_layers = 1;
        }

        std::map<SupportElement, Polygons> influence_areas; // Over this map will be iterated when merging, as such it has to be ordered to ensure deterministic results.
        std::unordered_map<SupportElement, Polygons> to_bp_areas, to_model_areas; // The area of these SupportElement is not set, to avoid to much allocation and deallocation on the heap
        std::vector<SupportElement*> bypass_merge_areas; // Different to the other maps of SupportElements as these here have the area already set, as they are already to be inserted into move_bounds.

        auto ta = std::chrono::high_resolution_clock::now();

        std::vector<SupportElement*> last_layer;
        last_layer.insert(last_layer.begin(), move_bounds[layer_idx].begin(), move_bounds[layer_idx].end());

        // ### Increase the influence areas by the allowed movement distance
        increaseAreas(to_bp_areas, to_model_areas, influence_areas, bypass_merge_areas, last_layer, layer_idx, merge_this_layer);

        auto tb = std::chrono::high_resolution_clock::now();
        if (merge_this_layer)
        {
            bool reduced_by_merging = false;
            size_t count_before_merge = influence_areas.size();
            // ### Calculate which influence areas overlap, and merge them into a new influence area (simplified: an intersection of influence areas that have such an intersection)
            mergeInfluenceAreas(m_volumes, m_config, to_bp_areas, to_model_areas, influence_areas, layer_idx);

            last_merge = layer_idx;
            reduced_by_merging = count_before_merge > influence_areas.size();
            if (!reduced_by_merging && !new_element)
            {
                merge_every_x_layers = std::min(max_merge_every_x_layers, merge_every_x_layers + 1);
            }
        }
        auto tc = std::chrono::high_resolution_clock::now();

        dur_inc += tb - ta;
        dur_merge += tc - tb;

        new_element = !move_bounds[layer_idx - 1].empty();

        // Save calculated elements to output, and allocate Polygons on heap, as they will not be changed again.
        for (const std::pair<const SupportElement, Polygons> &tup : influence_areas) {
            const SupportElement &elem = tup.first;
            validate_range(tup.second);
            validate_range(safeUnion(tup.second));
            Polygons* new_area = new Polygons(safeUnion(tup.second));
            SupportElement* next = new SupportElement(elem, new_area);
            move_bounds[layer_idx - 1].emplace(next);

            if (area(*new_area) < tiny_area_threshold) {
                BOOST_LOG_TRIVIAL(error) << "Insert Error of Influence area on layer " << layer_idx - 1 << ". Origin of " << elem.parents.size() << " areas. Was to bp " << elem.to_buildplate;
                TreeSupport::showError("Insert error of area after merge.\n", true);
            }
        }

        // Place already fully constructed elements in the output.
        for (SupportElement* elem : bypass_merge_areas) {
            if (area(*elem->area) < tiny_area_threshold) {
                BOOST_LOG_TRIVIAL(error) << "Insert Error of Influence area bypass on layer " << layer_idx - 1;
                TreeSupport::showError("Insert error of area after bypassing merge.\n", true);
            }
            move_bounds[layer_idx - 1].emplace(elem);
        }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
        progress_total += data_size_inverse * TREE_PROGRESS_AREA_CALC;
        Progress::messageProgress(Progress::Stage::SUPPORT, progress_total * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
#endif
    }

    BOOST_LOG_TRIVIAL(info) << "Time spent with creating influence areas' subtasks: Increasing areas " << dur_inc.count() / 1000000 << " ms merging areas: " << dur_merge.count() / 1000000 << " ms";
}


void TreeSupport::setPointsOnAreas(const SupportElement* elem)
{
    // Based on the branch center point of the current layer, the point on the next (further up) layer is calculated.

    if (elem->result_on_layer == Point(-1, -1))
    {
        BOOST_LOG_TRIVIAL(error) << "Uninitialized support element";
        TreeSupport::showError("Uninitialized support element. A branch may be missing.\n", true);
        return;
    }

    for (SupportElement* next_elem : elem->parents)
    {
        if (next_elem->result_on_layer != Point(-1, -1)) // if the value was set somewhere else it it kept. This happens when a branch tries not to move after being unable to create a roof.
            continue;

        Point from = elem->result_on_layer;
        if (! contains(*next_elem->area, from)) {
            moveInside(*next_elem->area, from, 0); // Move inside has edgecases (see tests) so DONT use Polygons.inside to confirm correct move, Error with distance 0 is <= 1
            // it is not required to check if how far this move moved a point as is can be larger than maximum_movement_distance. While this seems like a problem it may for example occur after merges.
        }
        next_elem->result_on_layer = from;
        // do not call recursive because then amount of layers would be restricted by the stack size
    }
}

bool TreeSupport::setToModelContact(std::vector<std::set<SupportElement*>>& move_bounds, SupportElement* first_elem, const LayerIndex layer_idx)
{
    if (first_elem->to_model_gracious)
    {
        SupportElement* check = first_elem;

        std::vector<SupportElement*> checked;
        LayerIndex last_successfull_layer = layer_idx;
        bool set = false;

        // check for every layer upwards, up to the point where this influence area was created (either by initial insert or merge) if the branch could be placed on it, and highest up layer index.

        for (LayerIndex layer_check = layer_idx; check->next_height >= layer_check; layer_check++)
        {
            if (! intersection(*check->area, m_volumes.getPlaceableAreas(m_config.getCollisionRadius(*check), layer_check)).empty()) {
                set = true;
                last_successfull_layer = layer_check;
            }
            checked.emplace_back(check);
            if (check->parents.size() == 1)
            {
                check = check->parents[0];
            }
            else
            {
                break; // reached merge point
            }
        }

        // Could not find valid placement, even though it should exist => error handling
        if (!set)
        {
            if (SUPPORT_TREE_ONLY_GRACIOUS_TO_MODEL)
            {
                BOOST_LOG_TRIVIAL(warning) << "No valid placement found for to model gracious element on layer " << layer_idx << ": REMOVING BRANCH";
                TreeSupport::showError("Could not fine valid placement on model! Removing this branch...", true);
                for (LayerIndex layer = layer_idx; layer <= first_elem->next_height; layer++)
                {
                    move_bounds[layer].erase(checked[layer - layer_idx]);
                    delete checked[layer - layer_idx]->area;
                    delete checked[layer - layer_idx];
                }
            }
            else
            {
                BOOST_LOG_TRIVIAL(warning) << "No valid placement found for to model gracious element on layer " << layer_idx;
                TreeSupport::showError("Could not fine valid placement on model! Just placing it down anyway. Could cause floating branches.", true);
                first_elem->to_model_gracious = false;
                return setToModelContact(move_bounds, first_elem, layer_idx);
            }
        }

        for (LayerIndex layer = layer_idx + 1; layer < last_successfull_layer - 1; layer++)
        {
            move_bounds[layer].erase(checked[layer - layer_idx]);
            delete checked[layer - layer_idx]->area;
            delete checked[layer - layer_idx];
        }

        // Guess a point inside the influence area, in which the branch will be placed in.
        Point best = checked[last_successfull_layer - layer_idx]->next_position;
        if (! contains(*checked[last_successfull_layer - layer_idx]->area, best))
            moveInside(*checked[last_successfull_layer - layer_idx]->area, best);
        checked[last_successfull_layer - layer_idx]->result_on_layer = best;

        BOOST_LOG_TRIVIAL(debug) << "Added gracious Support On Model Point (" << best.x() << "," << best.y() << "). The current layer is " << last_successfull_layer;

        return last_successfull_layer != layer_idx;
    }
    else // can not add graceful => just place it here and hope for the best
    {
        Point best = first_elem->next_position;
        if (! contains(*first_elem->area, best))
            moveInside(*first_elem->area, best);
        first_elem->result_on_layer = best;
        first_elem->to_model_gracious = false;
        BOOST_LOG_TRIVIAL(debug) << "Added NON gracious Support On Model Point (" << best.x() << "," << best.y() << "). The current layer is " << layer_idx;
        return false;
    }
}

void TreeSupport::createNodesFromArea(std::vector<std::set<SupportElement*>>& move_bounds)
{
    // Initialize points on layer 0, with a "random" point in the influence area. Point is chosen based on an inaccurate estimate where the branches will split into two, but every point inside the influence area would produce a valid result.
    for (SupportElement* init : move_bounds[0]) {
        Point p = init->next_position;
        if (! contains(*init->area, p))
            moveInside(*init->area, p, 0);
        init->result_on_layer = p;
        setPointsOnAreas(init); // also set the parent nodes, as these will be required for the first iteration of the loop below
    }

    for (LayerIndex layer_idx = 1; layer_idx < LayerIndex(move_bounds.size()); ++ layer_idx) {
        std::unordered_set<SupportElement*> remove;
        for (SupportElement* elem : move_bounds[layer_idx]) {
            bool removed = false;
            // check if the resulting center point is not yet set
            if (elem->result_on_layer == Point(-1, -1)) {
                if (elem->to_buildplate || (!elem->to_buildplate && elem->distance_to_top < m_config.min_dtt_to_model && !elem->supports_roof)) {
                    if (elem->to_buildplate) {
                        BOOST_LOG_TRIVIAL(error) << "Uninitialized Influence area targeting " << elem->target_position.x() << "," << elem->target_position.y() << ") "
                            "at target_height: " << elem->target_height << " layer: " << layer_idx;
                        TreeSupport::showError("Uninitialized support element! A branch could be missing or exist partially.", true);
                    }
                    remove.emplace(elem); // we dont need to remove yet the parents as they will have a lower dtt and also no result_on_layer set
                    removed = true;
                    for (SupportElement* parent : elem->parents)
                        parent->result_on_layer = Point(-1, -1); // When the roof was not able to generate downwards enough, the top elements may have not moved, and have result_on_layer already set. As this branch needs to be removed => all parents result_on_layer have to be invalidated.
                    continue;
                } else {
                    // set the point where the branch will be placed on the model
                    removed = setToModelContact(move_bounds, elem, layer_idx);
                    if (removed)
                        remove.emplace(elem);
                }
            }

            if (!removed)
                setPointsOnAreas(elem); // element is valid now setting points in the layer above
        }

        // delete all not needed support elements
        for (SupportElement* del : remove) {
            move_bounds[layer_idx].erase(del);
            delete del->area;
            delete del;
        }
        remove.clear();
    }
}

void TreeSupport::generateBranchAreas(
    std::vector<std::pair<LayerIndex, SupportElement*>>             &linear_data, 
    std::vector<std::unordered_map<SupportElement*, Polygons>>      &layer_tree_polygons, 
    const std::map<SupportElement*, SupportElement*>                &inverse_tree_order)
{
#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    double progress_total = TREE_PROGRESS_PRECALC_AVO + TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_GENERATE_NODES + TREE_PROGRESS_AREA_CALC;
    constexpr int progress_report_steps = 10;
#endif // SLIC3R_TREESUPPORTS_PROGRESS

    Polygon branch_circle; // Pre-generate a circle with correct diameter so that we don't have to recompute those (co)sines every time.
    for (unsigned int i = 0; i < SUPPORT_TREE_CIRCLE_RESOLUTION; ++ i) {
        const double angle = static_cast<double>(i) / SUPPORT_TREE_CIRCLE_RESOLUTION * (2.0 * M_PI);
        branch_circle.points.emplace_back(coord_t(cos(angle) * m_config.branch_radius), coord_t(sin(angle) * m_config.branch_radius));
    }

    std::vector<Polygons> linear_inserts(linear_data.size());

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    const size_t progress_inserts_check_interval = linear_data.size() / progress_report_steps;
#endif // SLIC3R_TREESUPPORTS_PROGRESS

    std::mutex critical_sections;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, linear_data.size()),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t idx = range.begin(); idx < range.end(); ++ idx) {
            const LayerIndex      layer_idx  = linear_data[idx].first;
            const SupportElement *elem       = linear_data[idx].second;
            const auto            it_elem    = inverse_tree_order.find(const_cast<SupportElement*>(elem));
            const SupportElement* child_elem = it_elem == inverse_tree_order.end() ? nullptr : it_elem->second;
            const coord_t         radius     = m_config.getRadius(*elem);
            bool                  parent_uses_min = false;

            // Calculate multiple ovalized circles, to connect with every parent and child. Also generate regular circle for the current layer. Merge all these into one area.
            std::vector<std::pair<Point, coord_t>> movement_directions{ std::pair<Point, coord_t>(Point(0, 0), radius) };
            if (!elem->skip_ovalisation) {
                if (child_elem != nullptr) {
                    const Point movement = child_elem->result_on_layer - elem->result_on_layer;
                    movement_directions.emplace_back(movement, radius);
                }
                for (SupportElement *parent : elem->parents) {
                    const Point movement = parent->result_on_layer - elem->result_on_layer;
                    movement_directions.emplace_back(movement, std::max(m_config.getRadius(*parent), m_config.support_line_width));
                    parent_uses_min |= parent->use_min_xy_dist;
                }
            }

            double max_speed = 0;
            auto generateArea = [&volumes = m_volumes, layer_idx, elem, &branch_circle, branch_radius = m_config.branch_radius, support_line_width = m_config.support_line_width, &movement_directions, &max_speed, parent_uses_min](
                    coord_t aoffset) {
                Polygons poly;

                for (std::pair<Point, coord_t> movement : movement_directions) {
                    max_speed = std::max(max_speed, movement.first.cast<double>().norm());

                    // Visualization: https://jsfiddle.net/0zvcq39L/2/
                    // Ovalizes the circle to an ellipse, that contains both old center and new target position.
                    double used_scale = (movement.second + aoffset) / (1.0 * branch_radius);
                    Point center_position = elem->result_on_layer + movement.first / 2;
                    const double moveX = movement.first.x() / (used_scale * branch_radius);
                    const double moveY = movement.first.y() / (used_scale * branch_radius);
                    const double vsize_inv = 0.5 / (0.01 + std::sqrt(moveX * moveX + moveY * moveY));

                    double matrix[] = {
                        used_scale * (1 + moveX * moveX * vsize_inv),
                        used_scale * (0 + moveX * moveY * vsize_inv),
                        used_scale * (0 + moveX * moveY * vsize_inv),
                        used_scale * (1 + moveY * moveY * vsize_inv),
                    };
                    Polygon circle;
                    for (Point vertex : branch_circle)
                        circle.points.emplace_back(center_position + Point(matrix[0] * vertex.x() + matrix[1] * vertex.y(), matrix[2] * vertex.x() + matrix[3] * vertex.y()));
                    poly.emplace_back(std::move(circle));
                }

                poly = diff_clipped(offset(union_(poly), std::min(coord_t(50), support_line_width / 4), jtMiter, 1.2), 
                    volumes.getCollision(0, layer_idx, parent_uses_min || elem->use_min_xy_dist)); // There seem to be some rounding errors, causing a branch to be a tiny bit further away from the model that it has to be. This can cause the tip to be slightly further away front the overhang (x/y wise) than optimal. This fixes it, and for every other part, 0.05mm will not be noticed.
                return poly;
            };

            bool fast_relative_movement = max_speed > radius * 0.75;

            // ensure branch area will not overlap with model/collision. This can happen because of e.g. ovalization or increase_until_radius.
            linear_inserts[idx] = generateArea(0);

            if (fast_relative_movement || m_config.getRadius(*elem) - m_config.getCollisionRadius(*elem) > m_config.support_line_width) {
                // simulate the path the nozzle will take on the outermost wall
                // if multiple parts exist, the outer line will not go all around the support part potentially causing support material to be printed mid air
                ExPolygons nozzle_path = offset_ex(linear_inserts[idx], -m_config.support_line_width / 2);
                if (nozzle_path.size() > 1) {
                    // Just try to make the area a tiny bit larger.
                    linear_inserts[idx] = generateArea(m_config.support_line_width / 2);
                    nozzle_path = offset_ex(linear_inserts[idx], -m_config.support_line_width / 2);

                    // if larger area did not fix the problem, all parts off the nozzle path that do not contain the center point are removed, hoping for the best
                    if (nozzle_path.size() > 1) {
                        Polygons polygons_with_correct_center;
                        for (ExPolygon &part : nozzle_path) {
                            if (part.contains(elem->result_on_layer))
                                polygons_with_correct_center = union_(polygons_with_correct_center, part);
                            else {
                                // try a fuzzy inside as sometimes the point should be on the border, but is not because of rounding errors...
                                Point from = elem->result_on_layer;
                                Polygons to = to_polygons(std::move(part));
                                moveInside(to, from, 0);
                                if ((elem->result_on_layer - from).cast<double>().norm() < scaled<double>(0.025))
                                    polygons_with_correct_center = union_(polygons_with_correct_center, to);
                            }
                        }
                        // Increase the area again, to ensure the nozzle path when calculated later is very similar to the one assumed above.
                        linear_inserts[idx] = offset(polygons_with_correct_center, m_config.support_line_width / 2, jtMiter, 1.2);
                        linear_inserts[idx] = diff_clipped(linear_inserts[idx], m_volumes.getCollision(0, linear_data[idx].first, parent_uses_min || elem->use_min_xy_dist));
                    }
                }
            }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
            if (idx % progress_inserts_check_interval == 0) {
                std::lock_guard<std::mutex> critical_section_progress(critical_sections);
                progress_total += TREE_PROGRESS_GENERATE_BRANCH_AREAS / progress_report_steps;
                Progress::messageProgress(Progress::Stage::SUPPORT, progress_total * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
            }
#endif
        }
    });

    // single threaded combining all elements to the right layers. ONLY COPYS DATA!
    for (coord_t i = 0; i < static_cast<coord_t>(linear_data.size()); i++)
        layer_tree_polygons[linear_data[i].first].emplace(linear_data[i].second, linear_inserts[i]);
}

void TreeSupport::smoothBranchAreas(std::vector<std::unordered_map<SupportElement*, Polygons>>& layer_tree_polygons)
{
#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    double progress_total = TREE_PROGRESS_PRECALC_AVO + TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_GENERATE_NODES + TREE_PROGRESS_AREA_CALC + TREE_PROGRESS_GENERATE_BRANCH_AREAS;
#endif // SLIC3R_TREESUPPORTS_PROGRESS

    const coord_t max_radius_change_per_layer = 1 + m_config.support_line_width / 2; // this is the upper limit a radius may change per layer. +1 to avoid rounding errors

    // smooth upwards
    for (LayerIndex layer_idx = 0; layer_idx < LayerIndex(layer_tree_polygons.size()) - 1; ++ layer_idx) {
        std::vector<std::pair<SupportElement*, Polygons>> processing;
        processing.insert(processing.end(), layer_tree_polygons[layer_idx].begin(), layer_tree_polygons[layer_idx].end());
        std::vector<std::vector<std::pair<SupportElement*, Polygons>>> update_next(processing.size()); // with this a lock can be avoided

        tbb::parallel_for(tbb::blocked_range<size_t>(0, processing.size()),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t processing_idx = range.begin(); processing_idx < range.end(); ++ processing_idx) {
                std::pair<SupportElement*, Polygons> data_pair = processing[processing_idx];
                double max_outer_wall_distance = 0;
                bool do_something = false;
                for (SupportElement* parent : data_pair.first->parents)
                    if (m_config.getRadius(*parent) != m_config.getCollisionRadius(*parent)) {
                        do_something = true;
                        max_outer_wall_distance = std::max(max_outer_wall_distance, (data_pair.first->result_on_layer - parent->result_on_layer).cast<double>().norm() - (m_config.getRadius(*data_pair.first) - m_config.getRadius(*parent)));
                    }
                max_outer_wall_distance += max_radius_change_per_layer; // As this change is a bit larger than what usually appears, lost radius can be slowly reclaimed over the layers.
                if (do_something) {
                    Polygons max_allowed_area = offset(data_pair.second, float(max_outer_wall_distance), jtMiter, 1.2);
                    for (SupportElement* parent : data_pair.first->parents)
                        if (m_config.getRadius(*parent) != m_config.getCollisionRadius(*parent))
                            update_next[processing_idx].emplace_back(std::pair<SupportElement*, Polygons>(parent, intersection(layer_tree_polygons[layer_idx + 1][parent], max_allowed_area)));
                }
            }
        });

        for (std::vector<std::pair<SupportElement*, Polygons>> data_vector : update_next)
            for (std::pair<SupportElement*, Polygons> data_pair : data_vector)
                layer_tree_polygons[layer_idx + 1][data_pair.first] = data_pair.second;
    }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    progress_total += TREE_PROGRESS_SMOOTH_BRANCH_AREAS / 2;
    Progress::messageProgress(Progress::Stage::SUPPORT, progress_total * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL); // It is just assumed that both smoothing loops together are one third of the time spent in this function. This was guessed. As the whole function is only 10%, and the smoothing is hard to predict a progress report in the loop may be not useful.
#endif

    // smooth downwards
    std::unordered_set<SupportElement*> updated_last_iteration;
    for (int layer_idx = int(layer_tree_polygons.size()) - 2; layer_idx >= 0; -- layer_idx) {
        std::vector<std::pair<SupportElement*, Polygons>> processing;
        processing.insert(processing.end(), layer_tree_polygons[layer_idx].begin(), layer_tree_polygons[layer_idx].end());
        std::vector<std::pair<SupportElement*, Polygons>> update_next(processing.size(), std::pair<SupportElement*, Polygons>(nullptr, Polygons())); // with this a lock can be avoided

        tbb::parallel_for(tbb::blocked_range<size_t>(0, processing.size()),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t processing_idx = range.begin(); processing_idx < range.end(); ++ processing_idx) {
                std::pair<SupportElement*, Polygons> data_pair = processing[processing_idx];
                bool do_something = false;
                Polygons max_allowed_area;
                for (size_t idx = 0; idx < data_pair.first->parents.size(); ++ idx) {
                    SupportElement* parent = data_pair.first->parents[idx];
                    coord_t max_outer_line_increase = max_radius_change_per_layer;
                    Polygons result = offset(layer_tree_polygons[layer_idx + 1][parent], max_outer_line_increase, jtMiter, 1.2);
                    Point direction = data_pair.first->result_on_layer - parent->result_on_layer;
                    // move the polygons object
                    for (auto& outer : result)
                        for (Point& p : outer)
                            p += direction;
                    append(max_allowed_area, std::move(result));
                    do_something = do_something || updated_last_iteration.count(parent) || m_config.getCollisionRadius(*parent) != m_config.getRadius(*parent);
                }

                if (do_something) {
                    Polygons result = intersection(max_allowed_area, data_pair.second);
                    if (area(result) < area(data_pair.second))
                        update_next[processing_idx] = std::pair<SupportElement*, Polygons>(data_pair.first, result);
                }
            }
        });

        updated_last_iteration.clear();
        for (std::pair<Slic3r::TreeSupport::SupportElement*, Polygons> data_pair : update_next)
            if (data_pair.first != nullptr) {
                updated_last_iteration.emplace(data_pair.first);
                layer_tree_polygons[layer_idx][data_pair.first] = data_pair.second;
            }
    }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    progress_total += TREE_PROGRESS_SMOOTH_BRANCH_AREAS / 2;
    Progress::messageProgress(Progress::Stage::SUPPORT, progress_total * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
#endif
}

void TreeSupport::dropNonGraciousAreas(
    std::vector<std::unordered_map<SupportElement*, Polygons>>  &layer_tree_polygons,
    const std::vector<std::pair<LayerIndex, SupportElement*>>   &linear_data,
    std::vector<std::vector<std::pair<LayerIndex, Polygons>>>   &dropped_down_areas,
    const std::map<SupportElement*, SupportElement*>            &inverse_tree_order)
{
    tbb::parallel_for(tbb::blocked_range<size_t>(0, linear_data.size()),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t idx = range.begin(); idx < range.end(); ++ idx) {
            SupportElement* elem = linear_data[idx].second;
            bool non_gracious_model_contact = !elem->to_model_gracious && !inverse_tree_order.count(elem); // if a element has no child, it connects to whatever is below as no support further down for it will exist.
            if (non_gracious_model_contact) {
                Polygons rest_support = layer_tree_polygons[linear_data[idx].first][elem];
                LayerIndex counter = 1;
                while (area(rest_support) > tiny_area_threshold && counter < linear_data[idx].first) {
                    rest_support = diff_clipped(rest_support, m_volumes.getCollision(0, linear_data[idx].first - counter, false));
                    dropped_down_areas[idx].emplace_back(linear_data[idx].first - counter, rest_support);
                    counter++;
                }
            }
        }
    });
}

void TreeSupport::finalizeInterfaceAndSupportAreas(
    const PrintObject               &print_object,
    const std::vector<Polygons>     &overhangs,
    std::vector<Polygons>           &support_layer_storage,
    std::vector<Polygons>           &support_roof_storage,

    SupportGeneratorLayersPtr   	&bottom_contacts,
    SupportGeneratorLayersPtr   	&top_contacts,
    SupportGeneratorLayersPtr       &intermediate_layers,
    SupportGeneratorLayerStorage    &layer_storage)
{
    InterfacePreference interface_pref = m_config.interface_preference; // InterfacePreference::SUPPORT_LINES_OVERWRITE_INTERFACE;

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    double progress_total = TREE_PROGRESS_PRECALC_AVO + TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_GENERATE_NODES + TREE_PROGRESS_AREA_CALC + TREE_PROGRESS_GENERATE_BRANCH_AREAS + TREE_PROGRESS_SMOOTH_BRANCH_AREAS;
#endif // SLIC3R_TREESUPPORTS_PROGRESS

    // Iterate over the generated circles in parallel and clean them up. Also add support floor.
    tbb::spin_mutex layer_storage_mutex;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, support_layer_storage.size()),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx) {
            // Most of the time in this function is this union call. Can take 300+ ms when a lot of areas are to be unioned.
            support_layer_storage[layer_idx] = smooth_outward(union_(support_layer_storage[layer_idx]), m_config.support_line_width); //FIXME was .smooth(50);
            //smooth_outward(closing(std::move(bottom), closing_distance + minimum_island_radius, closing_distance, SUPPORT_SURFACES_OFFSET_PARAMETERS), smoothing_distance) :

            // simplify a bit, to ensure the output does not contain outrageous amounts of vertices. Should not be necessary, just a precaution.
            support_layer_storage[layer_idx] = polygons_simplify(support_layer_storage[layer_idx], std::min(scaled<double>(0.03), double(m_config.resolution)));
            // Subtract support lines of the branches from the roof
            SupportGeneratorLayer*& support_roof = top_contacts[layer_idx];
            if (! support_roof_storage[layer_idx].empty() || support_roof != nullptr) {
                if (support_roof == nullptr) {
                    support_roof = &layer_allocate(layer_storage, layer_storage_mutex, SupporLayerType::TopContact, print_object.slicing_parameters(), layer_idx);
                    support_roof->polygons = union_(support_roof_storage[layer_idx]);
                } else
                    support_roof->polygons = union_(support_roof->polygons, support_roof_storage[layer_idx]);

                if (! support_roof->polygons.empty() &&
                    area(intersection(support_layer_storage[layer_idx], support_roof->polygons)) > tiny_area_threshold) {
                    switch (interface_pref) {
                        case InterfacePreference::INTERFACE_AREA_OVERWRITES_SUPPORT:
                            support_layer_storage[layer_idx] = diff(support_layer_storage[layer_idx], support_roof->polygons);
                            break;
                        case InterfacePreference::SUPPORT_AREA_OVERWRITES_INTERFACE:
                            support_roof->polygons = diff(support_roof->polygons, support_layer_storage[layer_idx]);
                            break;
    //FIXME
    #if 1
                        case InterfacePreference::INTERFACE_LINES_OVERWRITE_SUPPORT:
                        case InterfacePreference::SUPPORT_LINES_OVERWRITE_INTERFACE:
                            assert(false);
                            [[fallthrough]];
    #else
                        case InterfacePreference::INTERFACE_LINES_OVERWRITE_SUPPORT:
                        {
                            // Hatch the support roof interfaces, offset them by their line width and subtract them from support base.
                            Polygons interface_lines = offset(to_polylines(
                                generateSupportInfillLines(support_roof->polygons, true, layer_idx, m_config.support_roof_line_distance)),
                                m_config.support_roof_line_width / 2);
                            support_layer_storage[layer_idx] = diff(support_layer_storage[layer_idx], interface_lines);
                            break;
                        }
                        case InterfacePreference::SUPPORT_LINES_OVERWRITE_INTERFACE:
                        {
                            // Hatch the support roof interfaces, offset them by their line width and subtract them from support base.
                            Polygons tree_lines = union_(offset(to_polylines(
                                generateSupportInfillLines(support_layer_storage[layer_idx], false, layer_idx, m_config.support_line_distance, true)),
                                m_config.support_line_width / 2));
                            // do not draw roof where the tree is. I prefer it this way as otherwise the roof may cut of a branch from its support below.
                            support_roof->polygons = diff(support_roof->polygons, tree_lines);
                            break;
                        }
    #endif
                        case InterfacePreference::NOTHING:
                            break;
                    }
                }
            }

            // Subtract support floors from the support area and add them to the support floor instead.
            if (m_config.support_bottom_layers > 0 && !support_layer_storage[layer_idx].empty()) {
                SupportGeneratorLayer*& support_bottom = bottom_contacts[layer_idx];
                Polygons layer_outset = diff_clipped(
                    m_config.support_bottom_offset > 0 ? offset(support_layer_storage[layer_idx], m_config.support_bottom_offset, jtMiter, 1.2) : support_layer_storage[layer_idx],
                    m_volumes.getCollision(0, layer_idx, false));
                Polygons floor_layer;
                size_t layers_below = 0;
                while (layers_below <= m_config.support_bottom_layers) {
                    // one sample at 0 layers below, another at m_config.support_bottom_layers. In-between samples at m_config.performance_interface_skip_layers distance from each other.
                    const size_t sample_layer = static_cast<size_t>(std::max(0, (static_cast<int>(layer_idx) - static_cast<int>(layers_below)) - static_cast<int>(m_config.z_distance_bottom_layers)));
                    //FIXME subtract the wipe tower 
                    append(floor_layer, intersection(layer_outset, overhangs[sample_layer]));
                    if (layers_below < m_config.support_bottom_layers)
                        layers_below = std::min(layers_below + m_config.performance_interface_skip_layers, m_config.support_bottom_layers);
                    else
                        break;
                }
                if (! floor_layer.empty()) {
                    if (support_bottom == nullptr)
                        support_bottom = &layer_allocate(layer_storage, layer_storage_mutex, SupporLayerType::BottomContact, print_object.slicing_parameters(), layer_idx);
                    support_bottom->polygons = union_(floor_layer, support_bottom->polygons);
                    support_layer_storage[layer_idx] = diff_clipped(support_layer_storage[layer_idx], offset(support_bottom->polygons, scaled<float>(0.01), jtMiter, 1.2)); // Subtract the support floor from the normal support.
                }
            }

            if (! support_layer_storage[layer_idx].empty()) {
                SupportGeneratorLayer *&l = intermediate_layers[layer_idx];
                if (l == nullptr)
                    l = &layer_allocate(layer_storage, layer_storage_mutex, SupporLayerType::Base, print_object.slicing_parameters(), layer_idx);
                append(l->polygons, union_(support_layer_storage[layer_idx]));
            }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
            {
                std::lock_guard<std::mutex> critical_section_progress(critical_sections);
                progress_total += TREE_PROGRESS_FINALIZE_BRANCH_AREAS / support_layer_storage.size();
                Progress::messageProgress(Progress::Stage::SUPPORT, progress_total * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
            }
#endif
#if 0
            {
                std::lock_guard<std::mutex> lock(critical_sections);
                if (!storage.support.supportLayers[layer_idx].support_infill_parts.empty() || !storage.support.supportLayers[layer_idx].support_roof.empty())
                    storage.support.layer_nr_max_filled_layer = std::max(storage.support.layer_nr_max_filled_layer, static_cast<int>(layer_idx));
            }
#endif
        }
    });
}

void TreeSupport::drawAreas(
    PrintObject                             &print_object,
    const std::vector<Polygons>             &overhangs,
    std::vector<std::set<SupportElement*>>  &move_bounds,

    SupportGeneratorLayersPtr            	&bottom_contacts,
    SupportGeneratorLayersPtr   	        &top_contacts,
    SupportGeneratorLayersPtr               &intermediate_layers,
    SupportGeneratorLayerStorage            &layer_storage)
{
    std::vector<Polygons> support_layer_storage(move_bounds.size());
    std::vector<Polygons> support_roof_storage(move_bounds.size());
    std::map<SupportElement*, SupportElement*> inverese_tree_order; // in the tree structure only the parents can be accessed. Inverse this to be able to access the children.
    std::vector<std::pair<LayerIndex, SupportElement*>> linear_data; // All SupportElements are put into a layer independent storage to improve parallelization. Was added at a point in time where this function had performance issues.
                                                                     // These were fixed by creating less initial points, but i do not see a good reason to remove a working performance optimization.
    for (LayerIndex layer_idx = 0; layer_idx < LayerIndex(move_bounds.size()); ++ layer_idx) {
        for (SupportElement* elem : move_bounds[layer_idx]) {
            if ((layer_idx > 0 && ((!inverese_tree_order.count(elem) && elem->target_height == layer_idx) || (inverese_tree_order.count(elem) && inverese_tree_order[elem]->result_on_layer == Point(-1, -1))))) // we either come from nowhere at the final layer or we had invalid parents 2. should never happen but just to be sure
                continue;
            for (SupportElement* par : elem->parents)
                if (par->result_on_layer != Point(-1, -1))
                    inverese_tree_order.emplace(par, elem);
            linear_data.emplace_back(layer_idx, elem);
        }
    }

    std::vector<std::unordered_map<SupportElement*, Polygons>> layer_tree_polygons(move_bounds.size()); // reorder the processed data by layers again. The map also could be a vector<pair<SupportElement*,Polygons>>.
    auto t_start = std::chrono::high_resolution_clock::now();
    // Generate the circles that will be the branches.
    generateBranchAreas(linear_data, layer_tree_polygons, inverese_tree_order);
    auto t_generate = std::chrono::high_resolution_clock::now();
    // In some edgecases a branch may go though a hole, where the regular radius does not fit. This can result in an apparent jump in branch radius. As such this cases need to be caught and smoothed out.
    smoothBranchAreas(layer_tree_polygons);
    auto t_smooth = std::chrono::high_resolution_clock::now();
    // drop down all trees that connect non gracefully with the model
    std::vector<std::vector<std::pair<LayerIndex, Polygons>>> dropped_down_areas(linear_data.size());
    dropNonGraciousAreas(layer_tree_polygons, linear_data, dropped_down_areas, inverese_tree_order);
    auto t_drop = std::chrono::high_resolution_clock::now();
    // single threaded combining all dropped down support areas to the right layers. ONLY COPYS DATA!
    for (coord_t i = 0; i < static_cast<coord_t>(dropped_down_areas.size()); i++)
        for (std::pair<LayerIndex, Polygons> &pair : dropped_down_areas[i])
            append(support_layer_storage[pair.first], std::move(pair.second));

    // single threaded combining all support areas to the right layers. ONLY COPYS DATA!
    for (LayerIndex layer_idx = 0; layer_idx < LayerIndex(layer_tree_polygons.size()); ++ layer_idx) {
        auto &this_layer_tree_polygons = layer_tree_polygons[layer_idx];
        auto &this_roofs               = support_roof_storage[layer_idx];
        auto &this_layers              = support_layer_storage[layer_idx];
        size_t cnt_roofs  = 0;
        size_t cnt_layers = 0;
        for (const std::pair<const SupportElement*, Polygons> &data_pair : this_layer_tree_polygons)
            ++ (data_pair.first->missing_roof_layers > data_pair.first->distance_to_top ? cnt_roofs : cnt_layers);
        this_roofs.reserve(this_roofs.size() + cnt_roofs);
        this_layers.reserve(this_layers.size() + cnt_layers);
        for (const std::pair<const SupportElement*, Polygons> &data_pair : this_layer_tree_polygons) {
            auto &src = const_cast<Polygons&>(data_pair.second);
            std::move(std::begin(src), std::end(src), std::back_inserter(data_pair.first->missing_roof_layers > data_pair.first->distance_to_top ? this_roofs : this_layers));
        }
    }

    finalizeInterfaceAndSupportAreas(print_object, overhangs, support_layer_storage, support_roof_storage,
        bottom_contacts, top_contacts, intermediate_layers, layer_storage);
    auto t_end = std::chrono::high_resolution_clock::now();

    auto dur_gen_tips = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_generate - t_start).count();
    auto dur_smooth = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_smooth - t_generate).count();
    auto dur_drop = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_drop - t_smooth).count();
    auto dur_finalize = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_drop).count();

    BOOST_LOG_TRIVIAL(info) << 
        "Time used for drawing subfuctions: generateBranchAreas: " << dur_gen_tips << " ms "
        "smoothBranchAreas: " << dur_smooth << " ms "
        "dropNonGraciousAreas: " << dur_drop << " ms "
        "finalizeInterfaceAndSupportAreas " << dur_finalize << " ms";
}

} // namespace Slic3r
