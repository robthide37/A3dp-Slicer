// Tree supports by Thomas Rahm, losely based on Tree Supports by CuraEngine.
// Original source of Thomas Rahm's tree supports:
// https://github.com/ThomasRahm/CuraEngine
//
// Original CuraEngine copyright:
// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "TreeModelVolumes.hpp"
#include "TreeSupport.hpp"

#include "BuildVolume.hpp"
#include "ClipperUtils.hpp"
#include "Flow.hpp"
#include "Layer.hpp"
#include "Point.hpp"
#include "Print.hpp"
#include "PrintConfig.hpp"

#include <boost/log/trivial.hpp>

#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

namespace Slic3r
{

TreeSupportMeshGroupSettings::TreeSupportMeshGroupSettings(const PrintObject &print_object)
{
    const PrintConfig       &print_config       = print_object.print()->config();
    const PrintObjectConfig &config             = print_object.config();
    const SlicingParameters &slicing_params     = print_object.slicing_parameters();
//    const std::vector<unsigned int>  printing_extruders = print_object.object_extruders();

    // Support must be enabled and set to Tree style.
    assert(config.support_material);
    assert(config.support_material_style == smsTree);

    coordf_t external_perimeter_width = 0.;
    for (size_t region_id = 0; region_id < print_object.num_printing_regions(); ++ region_id) {
        const PrintRegion &region = print_object.printing_region(region_id);
        external_perimeter_width = std::max(external_perimeter_width, coordf_t(region.flow(print_object, frExternalPerimeter, config.layer_height).width()));
    }

    this->layer_height              = scaled<coord_t>(config.layer_height.value);
    this->resolution                = scaled<coord_t>(print_config.resolution.value);
    this->min_feature_size          = scaled<coord_t>(config.min_feature_size.value);
    this->support_angle             = M_PI / 2. - config.support_material_angle * M_PI / 180.;
    this->support_line_width        = support_material_flow(&print_object, config.layer_height).scaled_width();
    this->support_roof_line_width   = support_material_interface_flow(&print_object, config.layer_height).scaled_width();
    //FIXME add it to SlicingParameters and reuse in both tree and normal supports?
    this->support_bottom_enable     = config.support_material_interface_layers.value > 0 && config.support_material_bottom_interface_layers.value != 0;
    this->support_bottom_height     = this->support_bottom_enable ?
        (config.support_material_bottom_interface_layers.value > 0 ?
            config.support_material_bottom_interface_layers.value :
            config.support_material_interface_layers.value) * this->layer_height :
        0;
    this->support_material_buildplate_only = config.support_material_buildplate_only;
//    this->support_xy_overrides_z    = 
    this->support_xy_distance       = scaled<coord_t>(config.support_material_xy_spacing.get_abs_value(external_perimeter_width));
    this->support_xy_distance_overhang = std::min(this->support_xy_distance, scaled<coord_t>(0.5 * external_perimeter_width));
    this->support_top_distance      = scaled<coord_t>(slicing_params.gap_support_object);
    this->support_bottom_distance   = scaled<coord_t>(slicing_params.gap_object_support);
//    this->support_interface_skip_height =
//    this->support_infill_angles     = 
    this->support_roof_enable       = config.support_material_interface_layers.value > 0;
    this->support_roof_height       = config.support_material_interface_layers.value * this->layer_height;
//    this->minimum_roof_area         = 
//    this->support_roof_angles       = 
    this->support_roof_pattern      = config.support_material_interface_pattern;
    this->support_pattern           = config.support_material_pattern;
    this->support_line_spacing      = scaled<coord_t>(config.support_material_spacing.value);
//    this->support_bottom_offset     = 
    this->support_wall_count        = config.support_material_with_sheath ? 1 : 0;
    this->support_roof_line_distance = scaled<coord_t>(config.support_material_interface_spacing.value) + this->support_roof_line_width;
//    this->minimum_support_area      = 
//    this->minimum_bottom_area       = 
//    this->support_offset            = 
}

static Polygons calculateMachineBorderCollision(Polygon machine_border)
{
    // Put a border of 1m around the print volume so that we don't collide.
#if 1
    //FIXME just returning no border will let tree support legs collide with print bed boundary
    return {};
#else
    //FIXME offsetting by 1000mm easily overflows int32_tr coordinate.
    Polygons out = offset(machine_border, scaled<float>(1000.), jtMiter, 1.2);
    machine_border.reverse(); // Makes the polygon negative so that we subtract the actual volume from the collision area.
    out.emplace_back(std::move(machine_border));
    return out;
#endif
}

TreeModelVolumes::TreeModelVolumes(
    const PrintObject &print_object,
    const BuildVolume &build_volume,
    const coord_t max_move, const coord_t max_move_slow, size_t current_mesh_idx, 
    double progress_multiplier, double progress_offset, const std::vector<Polygons>& additional_excluded_areas) :
    m_max_move{ std::max(max_move - 2, coord_t(0)) }, m_max_move_slow{ std::max(max_move_slow - 2, coord_t(0)) }, m_progress_multiplier{ progress_multiplier }, m_progress_offset{ progress_offset }, 
    // -2 to avoid rounding errors
    m_machine_border{ calculateMachineBorderCollision(build_volume.polygon()) }
{
    std::unordered_map<size_t, size_t> mesh_to_layeroutline_idx;

#if 0
    for (size_t mesh_idx = 0; mesh_idx < storage.meshes.size(); ++ mesh_idx) {
        SliceMeshStorage mesh = storage.meshes[mesh_idx];
        bool added = false;
        for (size_t idx = 0; idx < m_layer_outlines.size(); ++ idx)
            if (TreeSupport::TreeSupportSettings(m_layer_outlines[idx].first) == TreeSupport::TreeSupportSettings(mesh.settings)) {
                added = true;
                mesh_to_layeroutline_idx[mesh_idx] = idx;
            }
        if (! added) {
            mesh_to_layeroutline_idx[mesh_idx] = m_layer_outlines.size();
            m_layer_outlines.emplace_back(mesh.settings, std::vector<Polygons>(storage.support.supportLayers.size(), Polygons()));
        }
    }
    for (size_t idx = 0; idx < m_layer_outlines.size(); ++ idx) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_layer_outlines[idx].second.size()),
            [&](const tbb::blocked_range<size_t> &range) {
            for (const size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx)
                m_layer_outlines[idx].second[layer_idx] = union_(m_layer_outlines[idx].second[layer_idx]); 
        });
    }
    m_current_outline_idx = mesh_to_layeroutline_idx[current_mesh_idx];

    m_support_rests_on_model = false;
    m_min_resolution  = std::numeric_limits<coord_t>::max();
    for (auto data_pair : m_layer_outlines) {
        m_support_rests_on_model |= ! data_pair.first.support_material_buildplate_only;
        m_min_resolution = std::min(m_min_resolution, data_pair.first.resolution);
    }
#else
    {
        m_anti_overhang = print_object.slice_support_blockers();
        mesh_to_layeroutline_idx[0] = 0;
        TreeSupportMeshGroupSettings mesh_settings(print_object);
        m_layer_outlines.emplace_back(mesh_settings, std::vector<Polygons>{});
        m_current_outline_idx = 0;
        std::vector<Polygons> &outlines = m_layer_outlines.front().second;
        outlines.reserve(print_object.layer_count());
        for (const Layer *layer : print_object.layers())
            outlines.emplace_back(to_polygons(expolygons_simplify(layer->lslices, mesh_settings.resolution)));
    }
#endif

    const TreeSupport::TreeSupportSettings &config = m_layer_outlines[m_current_outline_idx].first;
    if (! config.support_xy_overrides_z) {
        m_current_min_xy_dist = config.xy_min_distance;
        if (TreeSupport::TreeSupportSettings::has_to_rely_on_min_xy_dist_only)
            m_current_min_xy_dist = std::max(m_current_min_xy_dist, coord_t(100));
        m_current_min_xy_dist_delta = std::max(config.xy_distance - m_current_min_xy_dist, coord_t(0));
    } else {
        m_current_min_xy_dist = config.xy_distance;
        m_current_min_xy_dist_delta = 0;
    }
    m_increase_until_radius = config.increase_radius_until_radius;
    m_radius_0 = config.getRadius(0);

#if 0
    for (size_t mesh_idx = 0; mesh_idx < storage.meshes.size(); mesh_idx++) {
        SliceMeshStorage mesh = storage.meshes[mesh_idx];
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_layer_outlines[mesh_to_layeroutline_idx[mesh_idx]].second.size()),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx)
                if (layer_idx < mesh.layer_nr_max_filled_layer) {
                    Polygons outline = extractOutlineFromMesh(mesh, layer_idx);
                    append(m_layer_outlines[mesh_to_layeroutline_idx[mesh_idx]].second[layer_idx], outline);
                }
        });
    }
    if (! additional_excluded_areas.empty()) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_anti_overhang.size()),
            [&](const tbb::blocked_range<size_t> &range) {
            for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx) {
                if (layer_idx < coord_t(additional_excluded_areas.size()))
                    append(m_anti_overhang[layer_idx], additional_excluded_areas[layer_idx]);
    //          if (SUPPORT_TREE_AVOID_SUPPORT_BLOCKER)
    //              append(m_anti_overhang[layer_idx], storage.support.supportLayers[layer_idx].anti_overhang);
    //FIXME block wipe tower
    //          if (storage.primeTower.enabled)
    //              append(m_anti_overhang[layer_idx], layer_idx == 0 ? storage.primeTower.outer_poly_first_layer : storage.primeTower.outer_poly);
                m_anti_overhang[layer_idx] = union_(m_anti_overhang[layer_idx]);
            }
        });
    }
#endif
}

void TreeModelVolumes::precalculate(coord_t max_layer)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    m_precalculated = true;

    // Get the config corresponding to one mesh that is in the current group. Which one has to be irrelevant. Not the prettiest way to do this, but it ensures some calculations that may be a bit more complex like inital layer diameter are only done in once.
    TreeSupport::TreeSupportSettings config(m_layer_outlines[m_current_outline_idx].first);

    // calculate which radius each layer in the tip may have.
    std::unordered_set<coord_t> possible_tip_radiis;
    for (size_t dtt = 0; dtt <= config.tip_layers; dtt++) {
        possible_tip_radiis.emplace(ceilRadius(config.getRadius(dtt)));
        possible_tip_radiis.emplace(ceilRadius(config.getRadius(dtt) + m_current_min_xy_dist_delta));
    }
    // It theoretically may happen in the tip, that the radius can change so much in-between 2 layers, that a ceil step is skipped (as in there is a radius r so that ceilRadius(radius(dtt))<ceilRadius(r)<ceilRadius(radius(dtt+1))). As such a radius will not reasonable happen in the tree and it will most likely not be requested, there is no need to calculate them. So just skip these.
    for (coord_t radius_eval = ceilRadius(1); radius_eval <= config.branch_radius; radius_eval = ceilRadius(radius_eval + 1))
        if (!possible_tip_radiis.count(radius_eval))
            m_ignorable_radii.emplace(radius_eval);

    // it may seem that the required avoidance can be of a smaller radius when going to model (no initial layer diameter for to model branches)
    // but as for every branch going towards the bp, the to model avoidance is required to check for possible merges with to model branches, this assumption is in-fact wrong.
    std::unordered_map<coord_t, LayerIndex> radius_until_layer;
    // while it is possible to calculate, up to which layer the avoidance should be calculated, this simulation is easier to understand, and does not need to be adjusted if something of the radius calculation is changed.
    // Overhead with an assumed worst case of 6600 layers was about 2ms
    for (LayerIndex simulated_dtt = 0; simulated_dtt <= max_layer; simulated_dtt++) {
        const LayerIndex current_layer = max_layer - simulated_dtt;
        const coord_t max_regular_radius = ceilRadius(config.getRadius(simulated_dtt, 0) + m_current_min_xy_dist_delta);
        const coord_t max_min_radius = ceilRadius(config.getRadius(simulated_dtt, 0)); // the maximal radius that the radius with the min_xy_dist can achieve
        const coord_t max_initial_layer_diameter_radius = ceilRadius(config.recommendedMinRadius(current_layer) + m_current_min_xy_dist_delta);
        if (!radius_until_layer.count(max_regular_radius))
            radius_until_layer[max_regular_radius] = current_layer;
        if (!radius_until_layer.count(max_min_radius))
            radius_until_layer[max_min_radius] = current_layer;
        if (!radius_until_layer.count(max_initial_layer_diameter_radius))
            radius_until_layer[max_initial_layer_diameter_radius] = current_layer;
    }

    // Copy to deque to use in parallel for later.
    std::deque<RadiusLayerPair> relevant_avoidance_radiis;
    std::deque<RadiusLayerPair> relevant_avoidance_radiis_to_model;
    relevant_avoidance_radiis.insert(relevant_avoidance_radiis.end(), radius_until_layer.begin(), radius_until_layer.end());
    relevant_avoidance_radiis_to_model.insert(relevant_avoidance_radiis_to_model.end(), radius_until_layer.begin(), radius_until_layer.end());

    // Append additional radiis needed for collision.

    radius_until_layer[ceilRadius(m_increase_until_radius, false)] = max_layer; // To calculate collision holefree for every radius, the collision of radius m_increase_until_radius will be required.
    // Collision for radius 0 needs to be calculated everywhere, as it will be used to ensure valid xy_distance in drawAreas.
    radius_until_layer[0] = max_layer;
    if (m_current_min_xy_dist_delta != 0)
        radius_until_layer[m_current_min_xy_dist_delta] = max_layer;

    std::deque<RadiusLayerPair> relevant_collision_radiis;
    relevant_collision_radiis.insert(relevant_collision_radiis.end(), radius_until_layer.begin(), radius_until_layer.end()); // Now that required_avoidance_limit contains the maximum of ild and regular required radius just copy.


    // ### Calculate the relevant collisions
    calculateCollision(relevant_collision_radiis);

    // calculate a separate Collisions with all holes removed. These are relevant for some avoidances that try to avoid holes (called safe)
    std::deque<RadiusLayerPair> relevant_hole_collision_radiis;
    for (RadiusLayerPair key : relevant_avoidance_radiis)
        if (key.first < m_increase_until_radius + m_current_min_xy_dist_delta)
            relevant_hole_collision_radiis.emplace_back(key);

    // ### Calculate collisions without holes, build from regular collision
    calculateCollisionHolefree(relevant_hole_collision_radiis);

    auto t_coll = std::chrono::high_resolution_clock::now();

    // ### Calculate the relevant avoidances in parallel as far as possible
    {
        tbb::task_group task_group;
        task_group.run([this, relevant_avoidance_radiis]{ calculateAvoidance(relevant_avoidance_radiis); });
        task_group.run([this, relevant_avoidance_radiis]{ calculateWallRestrictions(relevant_avoidance_radiis); });
        if (m_support_rests_on_model)
            task_group.run([this, relevant_avoidance_radiis_to_model]{ 
                calculatePlaceables(relevant_avoidance_radiis_to_model);
                calculateAvoidanceToModel(relevant_avoidance_radiis_to_model);
            });
        task_group.wait();
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    auto dur_col = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_coll - t_start).count();
    auto dur_avo = 0.001 * std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_coll).count();

    BOOST_LOG_TRIVIAL(info) << "Precalculating collision took" << dur_col << " ms. Precalculating avoidance took " << dur_avo << " ms.";
}

/*!
 * \brief Checks a cache for a given RadiusLayerPair and returns it if it is found
 * \param key RadiusLayerPair of the requested areas. The radius will be calculated up to the provided layer.
 * \return A wrapped optional reference of the requested area (if it was found, an empty optional if nothing was found)
 */
std::optional<std::reference_wrapper<const Polygons>> getArea(const TreeModelVolumes::RadiusLayerPolygonCache& cache, const TreeModelVolumes::RadiusLayerPair& key)
{
    const auto it = cache.find(key);
    return it == cache.end() ? std::optional<std::reference_wrapper<const Polygons>>{} : std::optional<std::reference_wrapper<const Polygons>>{ it->second };
}

const Polygons& TreeModelVolumes::getCollision(coord_t radius, LayerIndex layer_idx, bool min_xy_dist) const
{
    coord_t orig_radius = radius;
    std::optional<std::reference_wrapper<const Polygons>> result;
    if (!min_xy_dist)
        radius += m_current_min_xy_dist_delta;

    // special case as if a radius 0 is requested it could be to ensure correct xy distance. As such it is beneficial if the collision is as close to the configured values as possible.
    if (orig_radius != 0)
        radius = ceilRadius(radius);
    RadiusLayerPair key{ radius, layer_idx };

    {
        std::lock_guard<std::mutex> critical_section_support_max_layer_nr(*m_critical_avoidance_cache);
        result = getArea(m_collision_cache, key);
    }
    if (result)
        return result.value().get();
    if (m_precalculated) {
        BOOST_LOG_TRIVIAL(warning) << "Had to calculate collision at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
        TreeSupport::showError("Not precalculated Collision requested.", false);
    }
    const_cast<TreeModelVolumes*>(this)->calculateCollision(key);
    return getCollision(orig_radius, layer_idx, min_xy_dist);
}

const Polygons& TreeModelVolumes::getCollisionHolefree(coord_t radius, LayerIndex layer_idx, bool min_xy_dist) const
{
    coord_t orig_radius = radius;
    std::optional<std::reference_wrapper<const Polygons>> result;
    if (!min_xy_dist)
        radius += m_current_min_xy_dist_delta;
    if (radius >= m_increase_until_radius + m_current_min_xy_dist_delta)
        return getCollision(orig_radius, layer_idx, min_xy_dist);

    RadiusLayerPair key{ radius, layer_idx };

    {
        std::lock_guard<std::mutex> critical_section_support_max_layer_nr(*m_critical_collision_cache_holefree);
        result = getArea(m_collision_cache_holefree, key);
    }
    if (result)
        return result.value().get();
    if (m_precalculated) {
        BOOST_LOG_TRIVIAL(warning) << "Had to calculate collision holefree at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
        TreeSupport::showError("Not precalculated Holefree Collision requested.", false);
    }
    const_cast<TreeModelVolumes*>(this)->calculateCollisionHolefree(key);
    return getCollisionHolefree(orig_radius, layer_idx, min_xy_dist);
}

const Polygons& TreeModelVolumes::getAvoidance(coord_t radius, LayerIndex layer_idx, AvoidanceType type, bool to_model, bool min_xy_dist) const
{
    if (layer_idx == 0) // What on the layer directly above buildplate do i have to avoid to reach the buildplate ...
        return getCollision(radius, layer_idx, min_xy_dist);

    coord_t orig_radius = radius;

    std::optional<std::reference_wrapper<const Polygons>> result;

    if (!min_xy_dist)
        radius += m_current_min_xy_dist_delta;
    radius = ceilRadius(radius);

    if (radius >= m_increase_until_radius + m_current_min_xy_dist_delta && type == AvoidanceType::FAST_SAFE) // no holes anymore by definition at this request
        type = AvoidanceType::FAST;

    const RadiusLayerPair key{ radius, layer_idx };

    const RadiusLayerPolygonCache* cache_ptr = nullptr;
    std::mutex* mutex_ptr;
    if (!to_model && type == AvoidanceType::FAST) {
        cache_ptr = &m_avoidance_cache;
        mutex_ptr = m_critical_avoidance_cache.get();
    } else if (!to_model && type == AvoidanceType::SLOW) {
        cache_ptr = &m_avoidance_cache_slow;
        mutex_ptr = m_critical_avoidance_cache_slow.get();
    } else if (!to_model && type == AvoidanceType::FAST_SAFE) {
        cache_ptr = &m_avoidance_cache_hole;
        mutex_ptr = m_critical_avoidance_cache_holefree.get();
    } else if (to_model && type == AvoidanceType::FAST) {
        cache_ptr = &m_avoidance_cache_to_model;
        mutex_ptr = m_critical_avoidance_cache_to_model.get();
    } else if (to_model && type == AvoidanceType::SLOW) {
        cache_ptr = &m_avoidance_cache_to_model_slow;
        mutex_ptr = m_critical_avoidance_cache_to_model_slow.get();
    } else if (to_model && type == AvoidanceType::FAST_SAFE) {
        cache_ptr = &m_avoidance_cache_hole_to_model;
        mutex_ptr = m_critical_avoidance_cache_holefree_to_model.get();
    } else {
        BOOST_LOG_TRIVIAL(error) << "Invalid Avoidance Request";
        TreeSupport::showError("Invalid Avoidance Request.\n", true);
    }

    if (to_model) {
        {
            std::lock_guard<std::mutex> critical_section(*mutex_ptr);
            result = getArea(*cache_ptr, key);
        }
        if (result)
            return result.value().get();
        if (m_precalculated) {
            BOOST_LOG_TRIVIAL(warning) << "Had to calculate Avoidance to model at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
            TreeSupport::showError("Not precalculated Avoidance(to model) requested.", false);
        }
        const_cast<TreeModelVolumes*>(this)->calculateAvoidanceToModel(key);
    } else {
        {
            std::lock_guard<std::mutex> critical_section(*mutex_ptr);
            result = getArea(*cache_ptr, key);
        }
        if (result)
            return result.value().get();
        if (m_precalculated) {
            BOOST_LOG_TRIVIAL(warning) << "Had to calculate Avoidance at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
            TreeSupport::showError("Not precalculated Avoidance(to buildplate) requested.", false);
        }
        const_cast<TreeModelVolumes*>(this)->calculateAvoidance(key);
    }
    return getAvoidance(orig_radius, layer_idx, type, to_model, min_xy_dist); // retrive failed and correct result was calculated. Now it has to be retrived.
}

const Polygons& TreeModelVolumes::getPlaceableAreas(coord_t radius, LayerIndex layer_idx) const
{
    std::optional<std::reference_wrapper<const Polygons>> result;
    const coord_t orig_radius = radius;
    radius = ceilRadius(radius);
    RadiusLayerPair key{ radius, layer_idx };

    {
        std::lock_guard<std::mutex> critical_section(*m_critical_placeable_areas_cache);
        result = getArea(m_placeable_areas_cache, key);
    }
    if (result)
        return result.value().get();
    if (m_precalculated) {
        BOOST_LOG_TRIVIAL(warning) << "Had to calculate Placeable Areas at radius " << radius << " and layer " << layer_idx << ", but precalculate was called. Performance may suffer!";
        TreeSupport::showError("Not precalculated Placeable areas requested.", false);
    }
    if (radius != 0)
        const_cast<TreeModelVolumes*>(this)->calculatePlaceables(key);
    else
        getCollision(0, layer_idx, true);
    return getPlaceableAreas(orig_radius, layer_idx);
}

const Polygons& TreeModelVolumes::getWallRestriction(coord_t radius, LayerIndex layer_idx, bool min_xy_dist) const
{
    if (layer_idx == 0) // Should never be requested as there will be no going below layer 0 ..., but just to be sure some semi-sane catch. Alternative would be empty Polygon.
        return getCollision(radius, layer_idx, min_xy_dist);

    coord_t orig_radius = radius;
    min_xy_dist = min_xy_dist && m_current_min_xy_dist_delta > 0;

    std::optional<std::reference_wrapper<const Polygons>> result;

    radius = ceilRadius(radius);
    const RadiusLayerPair key{ radius, layer_idx };

    const RadiusLayerPolygonCache* cache_ptr = min_xy_dist ? &m_wall_restrictions_cache_min : &m_wall_restrictions_cache;

    if (min_xy_dist) {
        {
            std::lock_guard<std::mutex> critical_section(*m_critical_wall_restrictions_cache_min);
            result = getArea(*cache_ptr, key);
        }
        if (result)
            return result.value().get();
        if (m_precalculated) {
            BOOST_LOG_TRIVIAL(warning) << "Had to calculate Wall restricions at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
            TreeSupport::showError("Not precalculated Wall restriction of minimum xy distance requested ).", false);
        }
    } else {
        {
            std::lock_guard<std::mutex> critical_section(*m_critical_wall_restrictions_cache);
            result = getArea(*cache_ptr, key);
        }
        if (result)
            return result.value().get();
        if (m_precalculated) {
            BOOST_LOG_TRIVIAL(warning) << "Had to calculate Wall restricions at radius " << key.first << " and layer " << key.second << ", but precalculate was called. Performance may suffer!";
            TreeSupport::showError("Not precalculated Wall restriction requested ).", false);
        }
    }
    const_cast<TreeModelVolumes*>(this)->calculateWallRestrictions(key);
    return getWallRestriction(orig_radius, layer_idx, min_xy_dist); // Retrieve failed and correct result was calculated. Now it has to be retrieved.
}

coord_t TreeModelVolumes::ceilRadius(coord_t radius, bool min_xy_dist) const
{
    if (! min_xy_dist)
        radius += m_current_min_xy_dist_delta;
    return ceilRadius(radius);
}

coord_t TreeModelVolumes::getRadiusNextCeil(coord_t radius, bool min_xy_dist) const
{
    coord_t ceiled_radius = ceilRadius(radius, min_xy_dist);
    if (! min_xy_dist)
        ceiled_radius -= m_current_min_xy_dist_delta;
    return ceiled_radius;
}

[[nodiscard]] static inline Polygons simplify(const Polygons &polygons, coord_t resolution)
{
    //FIXME
    return polygons;
}

#if 0
Polygons TreeModelVolumes::extractOutlineFromMesh(const PrintObject &print_object, LayerIndex layer_idx) const
{
    constexpr bool external_polys_only = false;
    Polygons total;

    // similar to SliceDataStorage.getLayerOutlines but only for one mesh instead of for everyone

    if (mesh.settings.get<bool>("infill_mesh") || mesh.settings.get<bool>("anti_overhang_mesh"))
        return Polygons();

    const SliceLayer& layer = mesh.layers[layer_idx];

    layer.getOutlines(total, external_polys_only);
    if (mesh.settings.get<ESurfaceMode>("magic_mesh_surface_mode") != ESurfaceMode::NORMAL)
        total = union_(total, layer.openPolyLines.offsetPolyLine(100));
    coord_t resolution = mesh.settings.get<coord_t>("resolution");
    return simplify(total, resolution);
}
#endif

LayerIndex TreeModelVolumes::getMaxCalculatedLayer(coord_t radius, const RadiusLayerPolygonCache& map) const
{
    int max_layer = -1;
    // the placeable on model areas do not exist on layer 0, as there can not be model below it. As such it may be possible that layer 1 is available, but layer 0 does not exist.
    const RadiusLayerPair key_layer_1(radius, 1);
    if (getArea(map, key_layer_1))
        max_layer = 1;
    while (map.count(RadiusLayerPair(radius, max_layer + 1)))
        ++ max_layer;
    return max_layer;
}

void TreeModelVolumes::calculateCollision(std::deque<RadiusLayerPair> keys)
{
    tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size()),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t i = range.begin(); i < range.end(); ++ i) {
            coord_t radius = keys[i].first;
            RadiusLayerPair key(radius, 0);
            RadiusLayerPolygonCache data_outer;
            RadiusLayerPolygonCache data_placeable_outer;
            for (size_t outline_idx = 0; outline_idx < m_layer_outlines.size(); outline_idx++)
            {
                RadiusLayerPolygonCache data;
                RadiusLayerPolygonCache data_placeable;

                const coord_t layer_height = m_layer_outlines[outline_idx].first.layer_height;
                const bool support_rests_on_this_model = ! m_layer_outlines[outline_idx].first.support_material_buildplate_only;
                const coord_t z_distance_bottom = m_layer_outlines[outline_idx].first.support_bottom_distance;
                const size_t z_distance_bottom_layers = round_up_divide(z_distance_bottom, layer_height);
                const coord_t z_distance_top_layers = round_up_divide(m_layer_outlines[outline_idx].first.support_top_distance, layer_height);
                const LayerIndex max_required_layer = keys[i].second + std::max(coord_t(1), z_distance_top_layers);
                const coord_t xy_distance = outline_idx == m_current_outline_idx ? m_current_min_xy_dist : m_layer_outlines[outline_idx].first.support_xy_distance;
                // technically this causes collision for the normal xy_distance to be larger by m_current_min_xy_dist_delta for all not currently processing meshes as this delta will be added at request time.
                // avoiding this would require saving each collision for each outline_idx separately.
                // and later for each avoidance... But avoidance calculation has to be for the whole scene and can NOT be done for each outline_idx separately and combined later.
                // so avoiding this inaccuracy seems infeasible as it would require 2x the avoidance calculations => 0.5x the performance.
                coord_t min_layer_bottom;
                {
                    std::lock_guard<std::mutex> critical_section(*m_critical_collision_cache);
                    min_layer_bottom = getMaxCalculatedLayer(radius, m_collision_cache) - int(z_distance_bottom_layers);
                }

                if (min_layer_bottom < 0)
                    min_layer_bottom = 0;
                for (LayerIndex layer_idx = min_layer_bottom; layer_idx <= max_required_layer; layer_idx++) {
                    key.second = layer_idx;
                    Polygons collision_areas = m_machine_border;
                    if (size_t(layer_idx) < m_layer_outlines[outline_idx].second.size())
                        append(collision_areas, m_layer_outlines[outline_idx].second[layer_idx]);
                    // jtRound is not needed here, as the overshoot can not cause errors in the algorithm, because no assumptions are made about the model.
                    // if a key does not exist when it is accessed it is added!
                    append(data[key], offset(union_ex(collision_areas), radius + xy_distance, ClipperLib::jtMiter, 1.2));
                }

                // Add layers below, to ensure correct support_bottom_distance. Also save placeable areas of radius 0, if required for this mesh.
                for (int layer_idx = int(max_required_layer); layer_idx >= min_layer_bottom; -- layer_idx) {
                    key.second = layer_idx;
                    for (size_t layer_offset = 1; layer_offset <= z_distance_bottom_layers && layer_idx - coord_t(layer_offset) > min_layer_bottom; ++ layer_offset)
                        append(data[key], data[RadiusLayerPair(radius, layer_idx - layer_offset)]);
                    if (support_rests_on_this_model && radius == 0 && layer_idx < coord_t(1 + keys[i].second)) {
                        RadiusLayerPair key_next_layer(radius, layer_idx + 1);
                        //data[key] = union_(data[key]);
                        Polygons above = data[key_next_layer];
                        // just to be sure the area is correctly unioned as otherwise difference may behave unexpectedly.
                        //FIXME Vojtech: Why m_anti_overhang.size() > layer_idx + 1? Why +1?
                        above = m_anti_overhang.size() > layer_idx + 1 ? union_(above, m_anti_overhang[layer_idx]) : union_(above);
                        data_placeable[key_next_layer] = union_(data_placeable[key_next_layer], diff(data[key], above));
                    }
                }

                // Add collision layers above to ensure correct support_top_distance.
                for (LayerIndex layer_idx = min_layer_bottom; layer_idx <= max_required_layer; ++ layer_idx) {
                    key.second = layer_idx;
                    Polygons collisions = std::move(data[key]);
                    for (coord_t layer_offset = 1; layer_offset <= z_distance_top_layers && layer_offset + layer_idx <= max_required_layer; ++ layer_offset)
                        append(collisions, data[RadiusLayerPair(radius, layer_idx + layer_offset)]);
                    data[key] = m_anti_overhang.size() > layer_idx ? union_(collisions, offset(union_ex(m_anti_overhang[layer_idx]), radius, ClipperLib::jtMiter, 1.2)) : union_(collisions);
                }

                for (int layer_idx = int(max_required_layer); layer_idx > keys[i].second; -- layer_idx) {
                    // all these dont have the correct z_distance_top_layers as they can still have areas above them
                    auto it = data.find(RadiusLayerPair(radius, layer_idx));
                    if (it != data.end())
                        data.erase(it);
                }

                for (auto pair : data)
                    data_outer[pair.first] = union_(data_outer[pair.first], simplify(pair.second, m_min_resolution));
                if (radius == 0) {
                    for (auto pair : data_placeable)
                        data_placeable_outer[pair.first] = union_(data_placeable_outer[pair.first], simplify(pair.second, m_min_resolution));
                }
            }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_progress);
                if (m_precalculated && m_precalculation_progress < TREE_PROGRESS_PRECALC_COLL) {
                    m_precalculation_progress += TREE_PROGRESS_PRECALC_COLL / keys.size();
                    Progress::messageProgress(Progress::Stage::SUPPORT, m_precalculation_progress * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
                }
            }
#endif

            {
                std::lock_guard<std::mutex> critical_section(*m_critical_collision_cache);
                m_collision_cache.insert(data_outer.begin(), data_outer.end());
            }
            if (radius == 0) {
                std::lock_guard<std::mutex> critical_section(*m_critical_placeable_areas_cache);
                m_placeable_areas_cache.insert(data_placeable_outer.begin(), data_placeable_outer.end());
            }
        }
    });
}

void TreeModelVolumes::calculateCollisionHolefree(std::deque<RadiusLayerPair> keys)
{
    LayerIndex max_layer = 0;
    for (long long unsigned int i = 0; i < keys.size(); i++)
        max_layer = std::max(max_layer, keys[i].second);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, max_layer + 1),
        [&](const tbb::blocked_range<size_t> &range) {
        for (size_t layer_idx = range.begin(); layer_idx < range.end(); ++ layer_idx) {
            RadiusLayerPolygonCache data;
            for (RadiusLayerPair key : keys) {
                // Logically increase the collision by m_increase_until_radius
                coord_t radius = key.first;
                coord_t increase_radius_ceil = ceilRadius(m_increase_until_radius, false) - ceilRadius(radius, true);
                // this union is important as otherwise holes(in form of lines that will increase to holes in a later step) can get unioned onto the area.
                Polygons col = offset(union_ex(getCollision(m_increase_until_radius, layer_idx, false)), 5 - increase_radius_ceil, ClipperLib::jtRound, scaled<float>(0.01));
                col = simplify(col, m_min_resolution);
                data[RadiusLayerPair(radius, layer_idx)] = col;
            }

            std::lock_guard<std::mutex> critical_section(*m_critical_collision_cache_holefree);
            m_collision_cache_holefree.insert(data.begin(), data.end());
        }
    });
}

// ensures offsets are only done in sizes with a max step size per offset while adding the collision offset after each step, this ensures that areas cannot glitch through walls defined by the collision when offsetting to fast
static Polygons safeOffset(const Polygons& me, coord_t distance, ClipperLib::JoinType jt, coord_t max_safe_step_distance, const Polygons& collision)
{
    const size_t steps = std::abs(distance / max_safe_step_distance);
    assert(int64_t(distance) * int64_t(max_safe_step_distance) >= 0);
    ExPolygons ret = union_ex(me);
    for (size_t i = 0; i < steps; ++ i)
        ret = union_ex(union_(offset(ret, max_safe_step_distance, jt, jt == jtRound ? scaled<float>(0.01) : 1.2), collision));
    return union_(offset(ret, distance % max_safe_step_distance, jt, jt == jtRound ? scaled<float>(0.01) : 1.2), collision);
}

void TreeModelVolumes::calculateAvoidance(std::deque<RadiusLayerPair> keys)
{
    // For every RadiusLayer pair there are 3 avoidances that have to be calculate, calculated in the same paralell_for loop for better paralellisation.
    const std::vector<AvoidanceType> all_types = { AvoidanceType::SLOW, AvoidanceType::FAST_SAFE, AvoidanceType::FAST };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size() * 3),
        [&, keys, all_types](const tbb::blocked_range<size_t> &range) {
        for (size_t iter_idx = range.begin(); iter_idx < range.end(); ++ iter_idx) {
            size_t key_idx = iter_idx / 3;
            {
                size_t type_idx = iter_idx % all_types.size();
                AvoidanceType type = all_types[type_idx];
                const bool slow = type == AvoidanceType::SLOW;
                const bool holefree = type == AvoidanceType::FAST_SAFE;

                coord_t radius = keys[key_idx].first;
                LayerIndex max_required_layer = keys[key_idx].second;

                // do not calculate not needed safe avoidances
                if (holefree && radius >= m_increase_until_radius + m_current_min_xy_dist_delta)
                    continue;

                const coord_t offset_speed = slow ? m_max_move_slow : m_max_move;
                const coord_t max_step_move = std::max(1.9 * radius, m_current_min_xy_dist * 1.9);
                RadiusLayerPair key(radius, 0);
                Polygons latest_avoidance;
                LayerIndex start_layer;
                {
                    std::lock_guard<std::mutex> critical_section(*(slow ? m_critical_avoidance_cache_slow : holefree ? m_critical_avoidance_cache_holefree : m_critical_avoidance_cache));
                    start_layer = 1 + getMaxCalculatedLayer(radius, slow ? m_avoidance_cache_slow : holefree ? m_avoidance_cache_hole : m_avoidance_cache);
                }
                if (start_layer > max_required_layer) {
                    BOOST_LOG_TRIVIAL(debug) << "Requested calculation for value already calculated ?";
                    continue;
                }
                start_layer = std::max(start_layer, LayerIndex(1)); // Ensure StartLayer is at least 1 as if no avoidance was calculated getMaxCalculatedLayer returns -1
                std::vector<std::pair<RadiusLayerPair, Polygons>> data(max_required_layer + 1, std::pair<RadiusLayerPair, Polygons>(RadiusLayerPair(radius, -1), Polygons()));

                latest_avoidance = getAvoidance(radius, start_layer - 1, type, false, true); // minDist as the delta was already added, also avoidance for layer 0 will return the collision.

                // ### main loop doing the calculation
                for (LayerIndex layer = start_layer; layer <= max_required_layer; layer++) {
                    key.second = layer;
                    Polygons col = (slow && radius < m_increase_until_radius + m_current_min_xy_dist_delta) || holefree ? 
                        getCollisionHolefree(radius, layer, true) :
                        getCollision(radius, layer, true);
                    latest_avoidance = safeOffset(latest_avoidance, -offset_speed, ClipperLib::jtRound, -max_step_move, col);
                    latest_avoidance = simplify(latest_avoidance, m_min_resolution);
                    data[layer] = std::pair<RadiusLayerPair, Polygons>(key, latest_avoidance);
                }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
                {
                    std::lock_guard<std::mutex> critical_section(*m_critical_progress);
                    if (m_precalculated && m_precalculation_progress < TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_PRECALC_AVO) {
                        m_precalculation_progress += m_support_rests_on_model ? 0.4 : 1 * TREE_PROGRESS_PRECALC_AVO / (keys.size() * 3);
                        Progress::messageProgress(Progress::Stage::SUPPORT, m_precalculation_progress * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
                    }
                }
#endif

                {
                    std::lock_guard<std::mutex> critical_section(*(slow ? m_critical_avoidance_cache_slow : holefree ? m_critical_avoidance_cache_holefree : m_critical_avoidance_cache));
                    (slow ? m_avoidance_cache_slow : holefree ? m_avoidance_cache_hole : m_avoidance_cache).insert(data.begin(), data.end());
                }
            }
        }
    });
}

void TreeModelVolumes::calculatePlaceables(std::deque<RadiusLayerPair> keys)
{
    tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size()),
        [&, keys](const tbb::blocked_range<size_t> &range) {
        for (size_t key_idx = range.begin(); key_idx < range.end(); ++ key_idx) {
            const coord_t radius = keys[key_idx].first;
            const LayerIndex max_required_layer = keys[key_idx].second;
            std::vector<std::pair<RadiusLayerPair, Polygons>> data(max_required_layer + 1, std::pair<RadiusLayerPair, Polygons>(RadiusLayerPair(radius, -1), Polygons()));
            RadiusLayerPair key(radius, 0);

            LayerIndex start_layer;
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_placeable_areas_cache);
                start_layer = 1 + getMaxCalculatedLayer(radius, m_placeable_areas_cache);
            }
            if (start_layer > max_required_layer) {
                BOOST_LOG_TRIVIAL(debug) << "Requested calculation for value already calculated ?";
                continue;
            }

            if (start_layer == 0) {
                data[0] = std::pair<RadiusLayerPair, Polygons>(key, diff(m_machine_border, getCollision(radius, 0, true)));
                start_layer = 1;
            }

            for (LayerIndex layer = start_layer; layer <= max_required_layer; layer++) {
                key.second = layer;
                Polygons placeable = getPlaceableAreas(0, layer);
                placeable = simplify(placeable, m_min_resolution); // it is faster to do this here in each thread than once in calculateCollision.
                placeable = offset(union_ex(placeable), - radius, jtMiter, 1.2);
                data[layer] = std::pair<RadiusLayerPair, Polygons>(key, placeable);
            }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_progress);
                if (m_precalculated && m_precalculation_progress < TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_PRECALC_AVO) {
                    m_precalculation_progress += 0.2 * TREE_PROGRESS_PRECALC_AVO / (keys.size());
                    Progress::messageProgress(Progress::Stage::SUPPORT, m_precalculation_progress * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
                }
            }
#endif

            {
                std::lock_guard<std::mutex> critical_section(*m_critical_placeable_areas_cache);
                m_placeable_areas_cache.insert(data.begin(), data.end());
            }
        }
    });
}

void TreeModelVolumes::calculateAvoidanceToModel(std::deque<RadiusLayerPair> keys)
{
    // For every RadiusLayer pair there are 3 avoidances that have to be calculated, calculated in the same parallel_for loop for better parallelization.
    const std::vector<AvoidanceType> all_types = { AvoidanceType::SLOW, AvoidanceType::FAST_SAFE, AvoidanceType::FAST };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size() * 3),
        [&, keys, all_types](const tbb::blocked_range<size_t> &range) {
        for (size_t iter_idx = range.begin(); iter_idx < range.end(); ++ iter_idx) {
            size_t key_idx = iter_idx / 3;
            size_t type_idx = iter_idx % all_types.size();
            AvoidanceType type = all_types[type_idx];
            bool slow = type == AvoidanceType::SLOW;
            bool holefree = type == AvoidanceType::FAST_SAFE;
            coord_t radius = keys[key_idx].first;
            LayerIndex max_required_layer = keys[key_idx].second;

            // do not calculate not needed safe avoidances
            if (holefree && radius >= m_increase_until_radius + m_current_min_xy_dist_delta)
                continue;

            getPlaceableAreas(radius, max_required_layer); // ensuring Placeableareas are calculated
            const coord_t offset_speed = slow ? m_max_move_slow : m_max_move;
            const coord_t max_step_move = std::max(1.9 * radius, m_current_min_xy_dist * 1.9);
            Polygons latest_avoidance;
            std::vector<std::pair<RadiusLayerPair, Polygons>> data(max_required_layer + 1, std::pair<RadiusLayerPair, Polygons>(RadiusLayerPair(radius, -1), Polygons()));
            RadiusLayerPair key(radius, 0);

            LayerIndex start_layer;

            {
                std::lock_guard<std::mutex> critical_section(*(slow ? m_critical_avoidance_cache_to_model_slow : holefree ? m_critical_avoidance_cache_holefree_to_model : m_critical_avoidance_cache_to_model));
                start_layer = 1 + getMaxCalculatedLayer(radius, slow ? m_avoidance_cache_to_model_slow : holefree ? m_avoidance_cache_hole_to_model : m_avoidance_cache_to_model);
            }
            if (start_layer > max_required_layer) {
                BOOST_LOG_TRIVIAL(debug) << "Requested calculation for value already calculated ?";
                continue;
            }
            start_layer = std::max(start_layer, LayerIndex(1));
            latest_avoidance = getAvoidance(radius, start_layer - 1, type, true, true); // minDist as the delta was already added, also avoidance for layer 0 will return the collision.

            // ### main loop doing the calculation
            for (LayerIndex layer = start_layer; layer <= max_required_layer; layer++)
            {
                key.second = layer;
                Polygons col = getCollision(radius, layer, true);

                if ((slow && radius < m_increase_until_radius + m_current_min_xy_dist_delta) || holefree)
                {
                    col = getCollisionHolefree(radius, layer, true);
                }
                else
                {
                    col = getCollision(radius, layer, true);
                }

                latest_avoidance = diff(safeOffset(latest_avoidance, -offset_speed, ClipperLib::jtRound, -max_step_move, col), getPlaceableAreas(radius, layer));

                latest_avoidance = simplify(latest_avoidance, m_min_resolution);
                data[layer] = std::pair<RadiusLayerPair, Polygons>(key, latest_avoidance);
            }

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_progress);
                if (m_precalculated && m_precalculation_progress < TREE_PROGRESS_PRECALC_COLL + TREE_PROGRESS_PRECALC_AVO) {
                    m_precalculation_progress += 0.4 * TREE_PROGRESS_PRECALC_AVO / (keys.size() * 3);
                    Progress::messageProgress(Progress::Stage::SUPPORT, m_precalculation_progress * m_progress_multiplier + m_progress_offset, TREE_PROGRESS_TOTAL);
                }
            }
#endif

            {
                std::lock_guard<std::mutex> critical_section(*(slow ? m_critical_avoidance_cache_to_model_slow : holefree ? m_critical_avoidance_cache_holefree_to_model : m_critical_avoidance_cache_to_model));
                (slow ? m_avoidance_cache_to_model_slow : holefree ? m_avoidance_cache_hole_to_model : m_avoidance_cache_to_model).insert(data.begin(), data.end());
            }
        }
    });
}


void TreeModelVolumes::calculateWallRestrictions(std::deque<RadiusLayerPair> keys)
{
    // Wall restrictions are mainly important when they represent actual walls that are printed, and not "just" the configured z_distance, because technically valid placement is no excuse for moving through a wall.
    // As they exist to prevent accidentially moving though a wall at high speed between layers like thie (x = wall,i = influence area,o= empty space,d = blocked area because of z distance) Assume maximum movement distance is two characters and maximum safe movement distance of one character

    /* Potential issue addressed by the wall restrictions: Influence area may lag through a wall
     *  layer z+1:iiiiiiiiiiioooo
     *  layer z+0:xxxxxiiiiiiiooo
     *  layer z-1:ooooixxxxxxxxxx
     */

    // The radius for the upper collission has to be 0 as otherwise one may not enter areas that may be forbidden on layer_idx but not one below (c = not an influence area even though it should ):
    /*
     *  layer z+1:xxxxxiiiiiioo
     *  layer z+0:dddddiiiiiiio
     *  layer z-1:dddocdddddddd
     */
    // Also there can not just the collision of the lower layer be used because if it were:
    /*
     *  layer z+1:dddddiiiiiiiiiio
     *  layer z+0:xxxxxddddddddddc
     *  layer z-1:dddddxxxxxxxxxxc
     */
    // Or of the upper layer be used because if it were:
    /*
     *  layer z+1:dddddiiiiiiiiiio
     *  layer z+0:xxxxcddddddddddc
     *  layer z-1:ddddcxxxxxxxxxxc
     */

    // And just offseting with maximum movement distance (and not in multiple steps) could cause:
    /*
     *  layer z:   oxiiiiiiiiioo
     *  layer z-1: ixiiiiiiiiiii
     */

    tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size()),
        [&, keys](const tbb::blocked_range<size_t> &range) {
        for (size_t key_idx = range.begin(); key_idx < range.end(); ++ key_idx) {
            coord_t radius = keys[key_idx].first;
            RadiusLayerPair key(radius, 0);
            coord_t min_layer_bottom;
            RadiusLayerPolygonCache data;
            RadiusLayerPolygonCache data_min;

            {
                std::lock_guard<std::mutex> critical_section(*m_critical_wall_restrictions_cache);
                min_layer_bottom = getMaxCalculatedLayer(radius, m_wall_restrictions_cache);
            }

            if (min_layer_bottom < 1)
                min_layer_bottom = 1;

            for (LayerIndex layer_idx = min_layer_bottom; layer_idx <= keys[key_idx].second; layer_idx++) {
                key.second = layer_idx;
                LayerIndex layer_idx_below = layer_idx - 1;
                Polygons wall_restriction = intersection(getCollision(0, layer_idx, false), getCollision(radius, layer_idx_below, true)); // radius contains m_current_min_xy_dist_delta already if required
                wall_restriction = simplify(wall_restriction, m_min_resolution);
                data.emplace(key, wall_restriction);
                if (m_current_min_xy_dist_delta > 0)
                {
                    Polygons wall_restriction_min = intersection(getCollision(0, layer_idx, true), getCollision(radius, layer_idx_below, true));
                    wall_restriction = simplify(wall_restriction_min, m_min_resolution);
                    data_min.emplace(key, wall_restriction_min);
                }
            }
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_wall_restrictions_cache);
                m_wall_restrictions_cache.insert(data.begin(), data.end());
            }
            {
                std::lock_guard<std::mutex> critical_section(*m_critical_wall_restrictions_cache_min);
                m_wall_restrictions_cache_min.insert(data_min.begin(), data_min.end());
            }
        }
    });
}

coord_t TreeModelVolumes::ceilRadius(coord_t radius) const
{
    if (radius == 0)
        return 0;

    if (radius <= m_radius_0)
        return m_radius_0;

    if (SUPPORT_TREE_USE_EXPONENTIAL_COLLISION_RESOLUTION)
    {
        // generate SUPPORT_TREE_PRE_EXPONENTIAL_STEPS of radiis before starting to exponentially increase them.

        coord_t exponential_result = SUPPORT_TREE_EXPONENTIAL_THRESHOLD * SUPPORT_TREE_EXPONENTIAL_FACTOR;
        const coord_t stepsize = (exponential_result - m_radius_0) / (SUPPORT_TREE_PRE_EXPONENTIAL_STEPS + 1);
        coord_t result = m_radius_0;
        for (size_t step = 0; step < SUPPORT_TREE_PRE_EXPONENTIAL_STEPS; step++) {
            result += stepsize;
            if (result >= radius && !m_ignorable_radii.count(result))
                return result;
        }

        while (exponential_result < radius || m_ignorable_radii.count(exponential_result))
            exponential_result = std::max(coord_t(exponential_result * SUPPORT_TREE_EXPONENTIAL_FACTOR), exponential_result + SUPPORT_TREE_COLLISION_RESOLUTION);
        return exponential_result;
    }
    else
    { // generates equidistant steps of size SUPPORT_TREE_COLLISION_RESOLUTION starting from m_radius_0. If SUPPORT_TREE_USE_EXPONENTIAL_COLLISION_RESOLUTION then this code is dead, and can safely be removed.
        coord_t ceil_step_n = (radius - m_radius_0) / SUPPORT_TREE_COLLISION_RESOLUTION;
        coord_t resulting_ceil = m_radius_0 + (ceil_step_n + ((radius - m_radius_0) % SUPPORT_TREE_COLLISION_RESOLUTION != 0)) * SUPPORT_TREE_COLLISION_RESOLUTION;
        return 
            radius <= m_radius_0 && radius != 0 ? m_radius_0 :
            m_ignorable_radii.count(resulting_ceil) ? ceilRadius(resulting_ceil + 1) : resulting_ceil;
    }
}

}
