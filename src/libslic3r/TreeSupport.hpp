// Tree supports by Thomas Rahm, losely based on Tree Supports by CuraEngine.
// Original source of Thomas Rahm's tree supports:
// https://github.com/ThomasRahm/CuraEngine
//
// Original CuraEngine copyright:
// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef slic3r_TreeSupport_hpp
#define slic3r_TreeSupport_hpp

#include "TreeModelVolumes.hpp"
#include "Point.hpp"

#include <boost/functional/hash.hpp> // For combining hashes

#include "BoundingBox.hpp"

// #define TREE_SUPPORT_SHOW_ERRORS

#define SUPPORT_TREE_CIRCLE_RESOLUTION 25 // The number of vertices in each circle.

#ifdef SLIC3R_TREESUPPORTS_PROGRESS
    // The various stages of the process can be weighted differently in the progress bar.
    // These weights are obtained experimentally using a small sample size. Sensible weights can differ drastically based on the assumed default settings and model.
    #define TREE_PROGRESS_TOTAL 10000
    #define TREE_PROGRESS_PRECALC_COLL TREE_PROGRESS_TOTAL * 0.1
    #define TREE_PROGRESS_PRECALC_AVO TREE_PROGRESS_TOTAL * 0.4
    #define TREE_PROGRESS_GENERATE_NODES TREE_PROGRESS_TOTAL * 0.1
    #define TREE_PROGRESS_AREA_CALC TREE_PROGRESS_TOTAL * 0.3
    #define TREE_PROGRESS_DRAW_AREAS TREE_PROGRESS_TOTAL * 0.1
    #define TREE_PROGRESS_GENERATE_BRANCH_AREAS TREE_PROGRESS_DRAW_AREAS / 3
    #define TREE_PROGRESS_SMOOTH_BRANCH_AREAS TREE_PROGRESS_DRAW_AREAS / 3
    #define TREE_PROGRESS_FINALIZE_BRANCH_AREAS TREE_PROGRESS_DRAW_AREAS / 3
#endif // SLIC3R_TREESUPPORTS_PROGRESS

#define SUPPORT_TREE_ONLY_GRACIOUS_TO_MODEL false
#define SUPPORT_TREE_AVOID_SUPPORT_BLOCKER true
#define SUPPORT_TREE_USE_EXPONENTIAL_COLLISION_RESOLUTION true
#define SUPPORT_TREE_EXPONENTIAL_THRESHOLD 1000
#define SUPPORT_TREE_EXPONENTIAL_FACTOR 1.5
#define SUPPORT_TREE_PRE_EXPONENTIAL_STEPS 1
#define SUPPORT_TREE_COLLISION_RESOLUTION 500 // Only has an effect if SUPPORT_TREE_USE_EXPONENTIAL_COLLISION_RESOLUTION is false

namespace Slic3r
{

using LayerIndex = int;

//FIXME
class Print;
class PrintObject;
class SupportGeneratorLayer;
using SupportGeneratorLayerStorage	= std::deque<SupportGeneratorLayer>;
using SupportGeneratorLayersPtr		= std::vector<SupportGeneratorLayer*>;
/*!
 * \brief Generates a tree structure to support your models.
 */

class TreeSupport
{
public:
    using AvoidanceType = TreeModelVolumes::AvoidanceType;
    enum class InterfacePreference
    {
        INTERFACE_AREA_OVERWRITES_SUPPORT,
        SUPPORT_AREA_OVERWRITES_INTERFACE,
        INTERFACE_LINES_OVERWRITE_SUPPORT,
        SUPPORT_LINES_OVERWRITE_INTERFACE,
        NOTHING
    };

    /*!
     * \brief Creates an instance of the tree support generator.
     */
    TreeSupport() = default;

    /*!
     * \brief Create the areas that need support.
     *
     * These areas are stored inside the given SliceDataStorage object.
     * \param storage The data storage where the mesh data is gotten from and
     * where the resulting support areas are stored.
     */
    void generateSupportAreas(Print &print, const BuildVolume &build_volume, const std::vector<size_t>& print_object_ids);
    void generateSupportAreas(PrintObject &print_object);


    //todo Remove! Only relevant for public BETA!
    static bool inline showed_critical=false;
    static bool inline showed_performance=false;
    static void showError(std::string message,bool critical);

    struct TreeSupportSettings; // forward declaration as we need some config values in the merge case

    struct AreaIncreaseSettings
    {
        AvoidanceType   type                { AvoidanceType::FAST };
        coord_t         increase_speed      { 0 };
        bool            increase_radius     { false };
        bool            no_error            { false };
        bool            use_min_distance    { false };
        bool            move                { false };
        bool operator==(const AreaIncreaseSettings& other) const
        {
            return increase_radius == other.increase_radius && increase_speed == other.increase_speed && type == other.type && 
                   no_error == other.no_error && use_min_distance == other.use_min_distance && move == other.move;
        }
    };

    struct SupportElement
    {
        explicit SupportElement(
            coord_t distance_to_top, size_t target_height, Point target_position, bool to_buildplate, bool to_model_gracious, bool use_min_xy_dist, size_t dont_move_until, 
            bool supports_roof, bool can_use_safe_radius, bool force_tips_to_roof, bool skip_ovalisation) : 
            target_height(target_height), target_position(target_position), next_position(target_position), next_height(target_height), effective_radius_height(distance_to_top), 
            to_buildplate(to_buildplate), distance_to_top(distance_to_top), area(nullptr), result_on_layer(target_position), increased_to_model_radius(0), to_model_gracious(to_model_gracious), 
            elephant_foot_increases(0), use_min_xy_dist(use_min_xy_dist), supports_roof(supports_roof), dont_move_until(dont_move_until), can_use_safe_radius(can_use_safe_radius), 
            last_area_increase(AreaIncreaseSettings{ AvoidanceType::FAST, 0, false, false, false, false }), missing_roof_layers(force_tips_to_roof ? dont_move_until : 0), skip_ovalisation(skip_ovalisation)
        {
        }


        explicit SupportElement(const SupportElement& elem, Polygons* newArea = nullptr)
            : // copy constructor with possibility to set a new area
              target_height(elem.target_height),
              target_position(elem.target_position),
              next_position(elem.next_position),
              next_height(elem.next_height),
              effective_radius_height(elem.effective_radius_height),
              to_buildplate(elem.to_buildplate),
              distance_to_top(elem.distance_to_top),
              area(newArea != nullptr ? newArea : elem.area),
              result_on_layer(elem.result_on_layer),
              increased_to_model_radius(elem.increased_to_model_radius),
              to_model_gracious(elem.to_model_gracious),
              elephant_foot_increases(elem.elephant_foot_increases),
              use_min_xy_dist(elem.use_min_xy_dist),
              supports_roof(elem.supports_roof),
              dont_move_until(elem.dont_move_until),
              can_use_safe_radius(elem.can_use_safe_radius),
              last_area_increase(elem.last_area_increase),
              missing_roof_layers(elem.missing_roof_layers),
              skip_ovalisation(elem.skip_ovalisation)

        {
            parents.insert(parents.begin(), elem.parents.begin(), elem.parents.end());
        }

        /*!
         * \brief Create a new Element for one layer below the element of the pointer supplied.
         */

        explicit SupportElement(SupportElement* element_above)
            : target_height(element_above->target_height),
              target_position(element_above->target_position),
              next_position(element_above->next_position),
              next_height(element_above->next_height),
              effective_radius_height(element_above->effective_radius_height),
              to_buildplate(element_above->to_buildplate),
              distance_to_top(element_above->distance_to_top + 1),
              area(element_above->area),
              result_on_layer(Point(-1, -1)), // set to invalid as we are a new node on a new layer
              increased_to_model_radius(element_above->increased_to_model_radius),
              to_model_gracious(element_above->to_model_gracious),
              elephant_foot_increases(element_above->elephant_foot_increases),
              use_min_xy_dist(element_above->use_min_xy_dist),
              supports_roof(element_above->supports_roof),
              dont_move_until(element_above->dont_move_until),
              can_use_safe_radius(element_above->can_use_safe_radius),
              last_area_increase(element_above->last_area_increase),
              missing_roof_layers(element_above->missing_roof_layers),
              skip_ovalisation(false)
        {
            parents = { element_above };
        }

        // ONLY to be called in merge as it assumes a few assurances made by it.
        explicit SupportElement(const SupportElement& first, const SupportElement& second, size_t next_height, Point next_position, coord_t increased_to_model_radius, const TreeSupportSettings& config) : next_position(next_position), next_height(next_height), area(nullptr), increased_to_model_radius(increased_to_model_radius), use_min_xy_dist(first.use_min_xy_dist || second.use_min_xy_dist), supports_roof(first.supports_roof || second.supports_roof), dont_move_until(std::max(first.dont_move_until, second.dont_move_until)), can_use_safe_radius(first.can_use_safe_radius || second.can_use_safe_radius), missing_roof_layers(std::min(first.missing_roof_layers, second.missing_roof_layers)), skip_ovalisation(false)

        {
            if (first.target_height > second.target_height)
            {
                target_height = first.target_height;
                target_position = first.target_position;
            }
            else
            {
                target_height = second.target_height;
                target_position = second.target_position;
            }
            effective_radius_height = std::max(first.effective_radius_height, second.effective_radius_height);
            distance_to_top = std::max(first.distance_to_top, second.distance_to_top);

            to_buildplate = first.to_buildplate && second.to_buildplate;
            to_model_gracious = first.to_model_gracious && second.to_model_gracious; // valid as we do not merge non-gracious with gracious

            AddParents(first.parents);
            AddParents(second.parents);

            elephant_foot_increases = 0;
            if (config.diameter_scale_bp_radius > 0)
            {
                coord_t foot_increase_radius = std::abs(std::max(config.getCollisionRadius(second), config.getCollisionRadius(first)) - config.getCollisionRadius(*this));
                // elephant_foot_increases has to be recalculated, as when a smaller tree with a larger elephant_foot_increases merge with a larger branch 
                // the elephant_foot_increases may have to be lower as otherwise the radius suddenly increases. This results often in a non integer value.
                elephant_foot_increases = foot_increase_radius / (config.branch_radius * (config.diameter_scale_bp_radius - config.diameter_angle_scale_factor));
            }

            // set last settings to the best out of both parents. If this is wrong, it will only cause a small performance penalty instead of weird behavior.
            last_area_increase = {
                std::min(first.last_area_increase.type, second.last_area_increase.type),
                std::min(first.last_area_increase.increase_speed, second.last_area_increase.increase_speed),
                first.last_area_increase.increase_radius || second.last_area_increase.increase_radius,
                first.last_area_increase.no_error || second.last_area_increase.no_error,
                first.last_area_increase.use_min_distance && second.last_area_increase.use_min_distance,
                first.last_area_increase.move || second.last_area_increase.move };
        }

        /*!
         * \brief The layer this support elements wants reach
         */
        LayerIndex target_height;

        /*!
         * \brief The position this support elements wants to support on layer=target_height
         */
        Point target_position;

        /*!
         * \brief The next position this support elements wants to reach. NOTE: This is mainly a suggestion regarding direction inside the influence area.
         */
        Point next_position;


        /*!
         * \brief The next height this support elements wants to reach
         */
        LayerIndex next_height;

        /*!
         * \brief The Effective distance to top of this element regarding radius increases and collision calculations.
         */

        size_t effective_radius_height;

        /*!
         * \brief The element trys to reach the buildplate
         */

        bool to_buildplate;

        /*!
         * \brief All elements in the layer above the current one that are supported by this element
         */
        std::vector<SupportElement*> parents;

        /*!
         * \brief The amount of layers this element is below the topmost layer of this branch.
         */
        size_t distance_to_top;

        /*!
         * \brief The resulting influence area.
         * Will only be set in the results of createLayerPathing, and will be nullptr inside!
         */
        Polygons* area;

        /*!
         * \brief The resulting center point around which a circle will be drawn later.
         * Will be set by setPointsOnAreas
         */
        Point result_on_layer = Point(-1, -1);
        /*!
         * \brief The amount of extra radius we got from merging branches that could have reached the buildplate, but merged with ones that can not.
         */
        coord_t increased_to_model_radius; // how much to model we increased only relevant for merging
        /*!
         * \brief Will the branch be able to rest completely on a flat surface, be it buildplate or model ?
         */
        bool to_model_gracious;

        /*!
         * \brief Counter about the times the elephant foot was increased. Can be fractions for merge reasons.
         */
        double elephant_foot_increases;

        /*!
         * \brief Whether the min_xy_distance can be used to get avoidance or similar. Will only be true if support_xy_overrides_z=Z overrides X/Y.
         */
        bool use_min_xy_dist;

        /*!
         * \brief True if this Element or any parent provides support to a support roof.
         */
        bool supports_roof;

        /*!
         * \brief The element trys not to move until this dtt is reached, is set to 0 if the element had to move.
         */
        size_t dont_move_until;

        /*!
         * \brief An influence area is considered safe when it can use the holefree avoidance <=> It will not have to encounter holes on its way downward.
         */
        bool can_use_safe_radius;

        /*!
         * \brief Settings used to increase the influence area to its current state.
         */
        AreaIncreaseSettings last_area_increase;

        /*!
         * \brief Amount of roof layers that were not yet added, because the branch needed to move.
         */
        size_t missing_roof_layers;

        /*!
         * \brief Skip the ovalisation to parent and children when generating the final circles.
         */
        bool skip_ovalisation;

        bool operator==(const SupportElement& other) const
        {
            return target_position == other.target_position && target_height == other.target_height;
        }

        bool operator<(const SupportElement& other) const // true if me < other
        {
            return !(*this == other) && !(*this > other);
        }
        bool operator>(const SupportElement& other) const
        {
            // Doesn't really have to make sense, only required for ordering in maps to ensure deterministic behavior.
            if (*this == other)
                return false;
            if (other.target_height != target_height)
                return other.target_height < target_height;
            return other.target_position.x() == target_position.x() ? other.target_position.y() < target_position.y() : other.target_position.x() < target_position.x();
        }

        void AddParents(const std::vector<SupportElement*>& adding)
        {
            for (SupportElement* ptr : adding)
            {
                parents.emplace_back(ptr);
            }
        }
    };

    /*!
     * \brief This struct contains settings used in the tree support. Thanks to this most functions do not need to know of meshes etc. Also makes the code shorter.
     */
    struct TreeSupportSettings
    {
        TreeSupportSettings() = default; // required for the definition of the config variable in the TreeSupport class.

        explicit TreeSupportSettings(const TreeSupportMeshGroupSettings& mesh_group_settings)
            : angle(mesh_group_settings.support_tree_angle),
              angle_slow(mesh_group_settings.support_tree_angle_slow),
              support_line_width(mesh_group_settings.support_line_width),
              layer_height(mesh_group_settings.layer_height),
              branch_radius(mesh_group_settings.support_tree_branch_diameter / 2),
              min_radius(mesh_group_settings.support_tree_tip_diameter / 2), // The actual radius is 50 microns larger as the resulting branches will be increased by 50 microns to avoid rounding errors effectively increasing the xydistance
              maximum_move_distance((angle < M_PI / 2.) ? (coord_t)(tan(angle) * layer_height) : std::numeric_limits<coord_t>::max()),
              maximum_move_distance_slow((angle_slow < M_PI / 2.) ? (coord_t)(tan(angle_slow) * layer_height) : std::numeric_limits<coord_t>::max()),
              support_bottom_layers(mesh_group_settings.support_bottom_enable ? (mesh_group_settings.support_bottom_height + layer_height / 2) / layer_height : 0),
              tip_layers(std::max((branch_radius - min_radius) / (support_line_width / 3), branch_radius / layer_height)), // Ensure lines always stack nicely even if layer height is large
              diameter_angle_scale_factor(sin(mesh_group_settings.support_tree_branch_diameter_angle) * layer_height / branch_radius),
              max_to_model_radius_increase(mesh_group_settings.support_tree_max_diameter_increase_by_merges_when_support_to_model / 2),
              min_dtt_to_model(round_up_divide(mesh_group_settings.support_tree_min_height_to_model, layer_height)),
              increase_radius_until_radius(mesh_group_settings.support_tree_branch_diameter / 2),
              increase_radius_until_layer(increase_radius_until_radius <= branch_radius ? tip_layers * (increase_radius_until_radius / branch_radius) : (increase_radius_until_radius - branch_radius) / (branch_radius * diameter_angle_scale_factor)),
              support_rests_on_model(! mesh_group_settings.support_material_buildplate_only),
              xy_distance(mesh_group_settings.support_xy_distance),
              bp_radius(mesh_group_settings.support_tree_bp_diameter / 2),
              diameter_scale_bp_radius(std::min(sin(0.7) * layer_height / branch_radius, 1.0 / (branch_radius / (support_line_width / 2.0)))), // Either 40? or as much as possible so that 2 lines will overlap by at least 50%, whichever is smaller.
              support_xy_overrides_z(mesh_group_settings.support_xy_overrides_z),
              xy_min_distance(support_xy_overrides_z ? xy_distance : mesh_group_settings.support_xy_distance_overhang),
              z_distance_top_layers(round_up_divide(mesh_group_settings.support_top_distance, layer_height)),
              z_distance_bottom_layers(round_up_divide(mesh_group_settings.support_bottom_distance, layer_height)),
              performance_interface_skip_layers(round_up_divide(mesh_group_settings.support_interface_skip_height, layer_height)),
              support_infill_angles(mesh_group_settings.support_infill_angles),
              support_roof_angles(mesh_group_settings.support_roof_angles),
              roof_pattern(mesh_group_settings.support_roof_pattern),
              support_pattern(mesh_group_settings.support_pattern),
              support_roof_line_width(mesh_group_settings.support_roof_line_width),
              support_line_spacing(mesh_group_settings.support_line_spacing),
              support_bottom_offset(mesh_group_settings.support_bottom_offset),
              support_wall_count(mesh_group_settings.support_wall_count),
              resolution(mesh_group_settings.resolution),
              support_roof_line_distance(mesh_group_settings.support_roof_line_distance), // in the end the actual infill has to be calculated to subtract interface from support areas according to interface_preference.
              settings(mesh_group_settings),
              min_feature_size(mesh_group_settings.min_feature_size)


        {
            layer_start_bp_radius = (bp_radius - branch_radius) / (branch_radius * diameter_scale_bp_radius);

            // safeOffsetInc can only work in steps of the size xy_min_distance in the worst case => xy_min_distance has to be a bit larger than 0 in this worst case and should be large enough for performance to not suffer extremely
            // When for all meshes the z bottom and top distance is more than one layer though the worst case is xy_min_distance + min_feature_size
            // This is not the best solution, but the only one to ensure areas can not lag though walls at high maximum_move_distance.
            if (has_to_rely_on_min_xy_dist_only)
                xy_min_distance = std::max(coord_t(100), xy_min_distance); // If set to low rounding errors WILL cause errors. Best to keep it above 25.

            xy_distance = std::max(xy_distance, xy_min_distance);

            // (logic) from getInterfaceAngles in FFFGcodeWriter.
            auto getInterfaceAngles = [&](std::vector<AngleRadians>& angles, SupportMaterialInterfacePattern pattern) {
                if (angles.empty())
                {
                    if (pattern == SupportMaterialInterfacePattern::smipConcentric)
                        angles.push_back(0); // Concentric has no rotation.
                    /*
                    else if (pattern == SupportMaterialInterfacePattern::TRIANGLES)
                        angles.push_back(90); // Triangular support interface shouldn't alternate every layer.
                    */
                    else {
                        if (TreeSupportSettings::some_model_contains_thick_roof) {
                            // Some roofs are quite thick.
                            // Alternate between the two kinds of diagonal: / and \ .
                            angles.push_back(M_PI / 4.);
                            angles.push_back(3. * M_PI / 4.);
                        }
                        if (angles.empty())
                            angles.push_back(M_PI / 2.); // Perpendicular to support lines.
                    }
                }
            };

            //getInterfaceAngles(support_infill_angles, support_pattern);
            support_infill_angles = { M_PI / 2. };
            getInterfaceAngles(support_roof_angles, roof_pattern);
//            const std::unordered_map<std::string, InterfacePreference> interface_map = { { "support_area_overwrite_interface_area", InterfacePreference::SUPPORT_AREA_OVERWRITES_INTERFACE }, { "interface_area_overwrite_support_area", InterfacePreference::INTERFACE_AREA_OVERWRITES_SUPPORT }, { "support_lines_overwrite_interface_area", InterfacePreference::SUPPORT_LINES_OVERWRITE_INTERFACE }, { "interface_lines_overwrite_support_area", InterfacePreference::INTERFACE_LINES_OVERWRITE_SUPPORT }, { "nothing", InterfacePreference::NOTHING } };
//            interface_preference = interface_map.at(mesh_group_settings.get<std::string>("support_interface_priority"));
//FIXME this was the default
//            interface_preference = InterfacePreference::SUPPORT_LINES_OVERWRITE_INTERFACE;
            interface_preference = InterfacePreference::SUPPORT_AREA_OVERWRITES_INTERFACE;
        }

    private:
        double angle;
        double angle_slow;
        std::vector<coord_t> known_z;

    public:
        // some static variables dependent on other meshes that are not currently processed.
        // Has to be static because TreeSupportConfig will be used in TreeModelVolumes as this reduces redundancy.
        inline static bool some_model_contains_thick_roof = false;
        inline static bool has_to_rely_on_min_xy_dist_only = false;
        /*!
         * \brief Width of a single line of support.
         */
        coord_t support_line_width;
        /*!
         * \brief Height of a single layer
         */
        coord_t layer_height;
        /*!
         * \brief Radius of a branch when it has left the tip.
         */
        coord_t branch_radius;
        /*!
         * \brief smallest allowed radius, required to ensure that even at DTT 0 every circle will still be printed
         */
        coord_t min_radius;
        /*!
         * \brief How far an influence area may move outward every layer at most.
         */
        coord_t maximum_move_distance;
        /*!
         * \brief How far every influence area will move outward every layer if possible.
         */
        coord_t maximum_move_distance_slow;
        /*!
         * \brief Amount of bottom layers. 0 if disabled.
         */
        size_t support_bottom_layers;
        /*!
         * \brief Amount of effectiveDTT increases are required to reach branch radius.
         */
        size_t tip_layers;
        /*!
         * \brief Factor by which to increase the branch radius.
         */
        double diameter_angle_scale_factor;
        /*!
         * \brief How much a branch resting on the model may grow in radius by merging with branches that can reach the buildplate.
         */
        coord_t max_to_model_radius_increase;
        /*!
         * \brief If smaller (in layers) than that, all branches to model will be deleted
         */
        size_t min_dtt_to_model;
        /*!
         * \brief Increase radius in the resulting drawn branches, even if the avoidance does not allow it. Will be cut later to still fit.
         */
        coord_t increase_radius_until_radius;
        /*!
         * \brief Same as increase_radius_until_radius, but contains the DTT at which the radius will be reached.
         */
        size_t increase_radius_until_layer;
        /*!
         * \brief True if the branches may connect to the model.
         */
        bool support_rests_on_model;
        /*!
         * \brief How far should support be from the model.
         */
        coord_t xy_distance;
        /*!
         * \brief Radius a branch should have when reaching the buildplate.
         */
        coord_t bp_radius;
        /*!
         * \brief The layer index at which an increase in radius may be required to reach the bp_radius.
         */
        coord_t layer_start_bp_radius;
        /*!
         * \brief Factor by which to increase the branch radius to reach the required bp_radius at layer 0. Note that this radius increase will not happen in the tip, to ensure the tip is structurally sound.
         */
        double diameter_scale_bp_radius;
        /*!
         * \brief Should Z distance override X/Y distance, or the other way around.
         */
        bool support_xy_overrides_z;
        /*!
         * \brief minimum xy_distance. Only relevant when Z overrides XY, otherwise equal to xy_distance-
         */
        coord_t xy_min_distance;
        /*!
         * \brief Amount of layers distance required the top of the support to the model
         */
        size_t z_distance_top_layers;
        /*!
         * \brief Amount of layers distance required from the top of the model to the bottom of a support structure.
         */
        size_t z_distance_bottom_layers;
        /*!
         * \brief used for performance optimization at the support floor. Should have no impact on the resulting tree.
         */
        size_t performance_interface_skip_layers;
        /*!
         * \brief User specified angles for the support infill.
         */
        std::vector<AngleRadians> support_infill_angles;
        /*!
         * \brief User specified angles for the support roof infill.
         */
        std::vector<AngleRadians> support_roof_angles;
        /*!
         * \brief Pattern used in the support roof. May contain non relevant data if support roof is disabled.
         */
        SupportMaterialInterfacePattern roof_pattern;
        /*!
         * \brief Pattern used in the support infill.
         */
        SupportMaterialPattern support_pattern;
        /*!
         * \brief Line width of the support roof.
         */
        coord_t support_roof_line_width;
        /*!
         * \brief Distance between support infill lines.
         */
        coord_t support_line_spacing;
        /*!
         * \brief Offset applied to the support floor area.
         */
        coord_t support_bottom_offset;
        /*
         * \brief Amount of walls the support area will have.
         */
        int support_wall_count;
        /*
         * \brief Maximum allowed deviation when simplifying.
         */
        coord_t resolution;
        /*
         * \brief Distance between the lines of the roof.
         */
        coord_t support_roof_line_distance;
        /*
         * \brief How overlaps of an interface area with a support area should be handled.
         */
        InterfacePreference interface_preference;

        /*
         * \brief The infill class wants a settings object. This one will be the correct one for all settings it uses.
         */
        TreeSupportMeshGroupSettings settings;

        /*
         * \brief Minimum thickness of any model features.
         */
        coord_t min_feature_size;

      public:
        bool operator==(const TreeSupportSettings& other) const
        {
            return branch_radius == other.branch_radius && tip_layers == other.tip_layers && diameter_angle_scale_factor == other.diameter_angle_scale_factor && layer_start_bp_radius == other.layer_start_bp_radius && bp_radius == other.bp_radius && diameter_scale_bp_radius == other.diameter_scale_bp_radius && min_radius == other.min_radius && xy_min_distance == other.xy_min_distance && // as a recalculation of the collision areas is required to set a new min_radius.
                   xy_distance - xy_min_distance == other.xy_distance - other.xy_min_distance && // if the delta of xy_min_distance and xy_distance is different the collision areas have to be recalculated.
                   support_rests_on_model == other.support_rests_on_model && increase_radius_until_layer == other.increase_radius_until_layer && min_dtt_to_model == other.min_dtt_to_model && max_to_model_radius_increase == other.max_to_model_radius_increase && maximum_move_distance == other.maximum_move_distance && maximum_move_distance_slow == other.maximum_move_distance_slow && z_distance_bottom_layers == other.z_distance_bottom_layers && support_line_width == other.support_line_width && 
                   support_xy_overrides_z == other.support_xy_overrides_z && support_line_spacing == other.support_line_spacing && support_roof_line_width == other.support_roof_line_width && // can not be set on a per-mesh basis currently, so code to enable processing different roof line width in the same iteration seems useless.
                   support_bottom_offset == other.support_bottom_offset && support_wall_count == other.support_wall_count && support_pattern == other.support_pattern && roof_pattern == other.roof_pattern && // can not be set on a per-mesh basis currently, so code to enable processing different roof patterns in the same iteration seems useless.
                   support_roof_angles == other.support_roof_angles && support_infill_angles == other.support_infill_angles && increase_radius_until_radius == other.increase_radius_until_radius && support_bottom_layers == other.support_bottom_layers && layer_height == other.layer_height && z_distance_top_layers == other.z_distance_top_layers && resolution == other.resolution && // Infill generation depends on deviation and resolution.
                   support_roof_line_distance == other.support_roof_line_distance && interface_preference == other.interface_preference
                   && min_feature_size == other.min_feature_size // interface_preference should be identical to ensure the tree will correctly interact with the roof.
                   // The infill class now wants the settings object and reads a lot of settings, and as the infill class is used to calculate support roof lines for interface-preference. Not all of these may be required to be identical, but as I am not sure, better safe than sorry
#if 0
                    && (interface_preference == InterfacePreference::INTERFACE_AREA_OVERWRITES_SUPPORT || interface_preference == InterfacePreference::SUPPORT_AREA_OVERWRITES_INTERFACE
                    // Perimeter generator parameters
                       || 
                            (settings.get<bool>("fill_outline_gaps") == other.settings.get<bool>("fill_outline_gaps") && 
                             settings.get<coord_t>("min_bead_width") == other.settings.get<coord_t>("min_bead_width") && 
                             settings.get<AngleRadians>("wall_transition_angle") == other.settings.get<AngleRadians>("wall_transition_angle") && 
                             settings.get<coord_t>("wall_transition_length") == other.settings.get<coord_t>("wall_transition_length") && 
                             settings.get<Ratio>("wall_split_middle_threshold") == other.settings.get<Ratio>("wall_split_middle_threshold") && 
                             settings.get<Ratio>("wall_add_middle_threshold") == other.settings.get<Ratio>("wall_add_middle_threshold") && 
                             settings.get<int>("wall_distribution_count") == other.settings.get<int>("wall_distribution_count") && 
                             settings.get<coord_t>("wall_transition_filter_distance") == other.settings.get<coord_t>("wall_transition_filter_distance") && 
                             settings.get<coord_t>("wall_transition_filter_deviation") == other.settings.get<coord_t>("wall_transition_filter_deviation") && 
                             settings.get<coord_t>("wall_line_width_x") == other.settings.get<coord_t>("wall_line_width_x") && 
                             settings.get<int>("meshfix_maximum_extrusion_area_deviation") == other.settings.get<int>("meshfix_maximum_extrusion_area_deviation"))
                        )
#endif
                ;
        }


        /*!
         * \brief Get the Distance to top regarding the real radius this part will have. This is different from distance_to_top, which is can be used to calculate the top most layer of the branch.
         * \param elem[in] The SupportElement one wants to know the effectiveDTT
         * \return The Effective DTT.
         */
        [[nodiscard]] inline size_t getEffectiveDTT(const TreeSupport::SupportElement& elem) const
        {
            return elem.effective_radius_height < increase_radius_until_layer ? (elem.distance_to_top < increase_radius_until_layer ? elem.distance_to_top : increase_radius_until_layer) : elem.effective_radius_height;
        }

        /*!
         * \brief Get the Radius part will have based on numeric values.
         * \param distance_to_top[in] The effective distance_to_top of the element
         * \param elephant_foot_increases[in] The elephant_foot_increases of the element.
         * \return The radius an element with these attributes would have.
         */
        [[nodiscard]] inline coord_t getRadius(size_t distance_to_top, const double elephant_foot_increases = 0) const
        {
            return (distance_to_top <= tip_layers ? min_radius + (branch_radius - min_radius) * distance_to_top / tip_layers : // tip
                           branch_radius + // base
                               branch_radius * (distance_to_top - tip_layers) * diameter_angle_scale_factor)
                   + // gradual increase
                   branch_radius * elephant_foot_increases * (std::max(diameter_scale_bp_radius - diameter_angle_scale_factor, 0.0));
        }

        /*!
         * \brief Get the Radius, that this element will have.
         * \param elem[in] The Element.
         * \return The radius the element has.
         */
        [[nodiscard]] inline coord_t getRadius(const TreeSupport::SupportElement& elem) const
        {
            return getRadius(getEffectiveDTT(elem), elem.elephant_foot_increases);
        }

        /*!
         * \brief Get the collision Radius of this Element. This can be smaller then the actual radius, as the drawAreas will cut off areas that may collide with the model.
         * \param elem[in] The Element.
         * \return The collision radius the element has.
         */
        [[nodiscard]] inline coord_t getCollisionRadius(const TreeSupport::SupportElement& elem) const
        {
            return getRadius(elem.effective_radius_height, elem.elephant_foot_increases);
        }

        /*!
         * \brief Get the Radius an element should at least have at a given layer.
         * \param layer_idx[in] The layer.
         * \return The radius every element should aim to achieve.
         */
        [[nodiscard]] inline coord_t recommendedMinRadius(LayerIndex layer_idx) const
        {
            double scale = (layer_start_bp_radius - int(layer_idx)) * diameter_scale_bp_radius;
            return scale > 0 ? branch_radius + branch_radius * scale : 0;
        }

        /*!
         * \brief Return on which z in microns the layer will be printed. Used only for support infill line generation.
         * \param layer_idx[in] The layer.
         * \return The radius every element should aim to achieve.
         */
        [[nodiscard]] inline coord_t getActualZ(LayerIndex layer_idx)
        {
            return layer_idx < coord_t(known_z.size()) ? known_z[layer_idx] : (layer_idx - known_z.size()) * layer_height + known_z.size() ? known_z.back() : 0;
        }

        /*!
         * \brief Set the z every Layer is printed at. Required for getActualZ to work
         * \param z[in] The z every LayerIndex is printed. Vector is used as a map<LayerIndex,coord_t> with the index of each element being the corresponding LayerIndex
         * \return The radius every element should aim to achieve.
         */
        void setActualZ(std::vector<coord_t>& z)
        {
            known_z = z;
        }
    };

private:
    /*!
     * \brief Creates the initial influence areas (that can later be propagated down) by placing them below the overhang.
     *
     * Generates Points where the Model should be supported and creates the areas where these points have to be placed.
     *
     * \param mesh[in] The mesh that is currently processed.
     * \param move_bounds[out] Storage for the influence areas.
     * \param storage[in] Background storage, required for adding roofs.
     */
    void generateInitialAreas(const PrintObject &print_object, 
        std::vector<std::set<SupportElement*>> &move_bounds,
        SupportGeneratorLayersPtr               &top_contacts,
        SupportGeneratorLayersPtr               &top_interface_layers,
        SupportGeneratorLayerStorage            &layer_storage);

    /*!
     * \brief Checks if an influence area contains a valid subsection and returns the corresponding metadata and the new Influence area.
     *
     * Calculates an influence areas of the layer below, based on the influence area of one element on the current layer.
     * Increases every influence area by maximum_move_distance_slow. If this is not enough, as in we would change our gracious or to_buildplate status the influence areas are instead increased by maximum_move_distance_slow.
     * Also ensures that increasing the radius of a branch, does not cause it to change its status (like to_buildplate ). If this were the case, the radius is not increased instead.
     *
     * Warning: The used format inside this is different as the SupportElement does not have a valid area member. Instead this area is saved as value of the dictionary. This was done to avoid not needed heap allocations.
     *
     * \param settings[in] Which settings have to be used to check validity.
     * \param layer_idx[in] Number of the current layer.
     * \param parent[in] The metadata of the parents influence area.
     * \param relevant_offset[in] The maximal possible influence area. No guarantee regarding validity with current layer collision required, as it is ensured in-function!
     * \param to_bp_data[out] The part of the Influence area that can reach the buildplate.
     * \param to_model_data[out] The part of the Influence area that do not have to reach the buildplate. This has overlap with new_layer_data.
     * \param increased[out]  Area than can reach all further up support points. No assurance is made that the buildplate or the model can be reached in accordance to the user-supplied settings.
     * \param overspeed[in] How much should the already offset area be offset again. Usually this is 0.
     * \param mergelayer[in] Will the merge method be called on this layer. This information is required as some calculation can be avoided if they are not required for merging.
     * \return A valid support element for the next layer regarding the calculated influence areas. Empty if no influence are can be created using the supplied influence area and settings.
     */
    std::optional<TreeSupport::SupportElement> increaseSingleArea(AreaIncreaseSettings settings, LayerIndex layer_idx, SupportElement* parent, const Polygons& relevant_offset, Polygons& to_bp_data, Polygons& to_model_data, Polygons& increased, const coord_t overspeed, const bool mergelayer);
    /*!
     * \brief Increases influence areas as far as required.
     *
     * Calculates influence areas of the layer below, based on the influence areas of the current layer.
     * Increases every influence area by maximum_move_distance_slow. If this is not enough, as in it would change the gracious or to_buildplate status, the influence areas are instead increased by maximum_move_distance.
     * Also ensures that increasing the radius of a branch, does not cause it to change its status (like to_buildplate ). If this were the case, the radius is not increased instead.
     *
     * Warning: The used format inside this is different as the SupportElement does not have a valid area member. Instead this area is saved as value of the dictionary. This was done to avoid not needed heap allocations.
     *
     * \param to_bp_areas[out] Influence areas that can reach the buildplate
     * \param to_model_areas[out] Influence areas that do not have to reach the buildplate. This has overlap with new_layer_data, as areas that can reach the buildplate are also considered valid areas to the model.
     * This redundancy is required if a to_buildplate influence area is allowed to merge with a to model influence area.
     * \param influence_areas[out] Area than can reach all further up support points. No assurance is made that the buildplate or the model can be reached in accordance to the user-supplied settings.
     * \param bypass_merge_areas[out] Influence areas ready to be added to the layer below that do not need merging.
     * \param last_layer[in] Influence areas of the current layer.
     * \param layer_idx[in] Number of the current layer.
     * \param mergelayer[in] Will the merge method be called on this layer. This information is required as some calculation can be avoided if they are not required for merging.
     */
    void increaseAreas(std::unordered_map<SupportElement, Polygons>& to_bp_areas, std::unordered_map<SupportElement, Polygons>& to_model_areas, std::map<SupportElement, Polygons>& influence_areas, std::vector<SupportElement*>& bypass_merge_areas, const std::vector<SupportElement*>& last_layer, const LayerIndex layer_idx, const bool mergelayer);

    /*!
     * \brief Propagates influence downwards, and merges overlapping ones.
     *
     * \param move_bounds[in,out] All currently existing influence areas
     */
    void createLayerPathing(std::vector<std::set<SupportElement*>>& move_bounds);


    /*!
     * \brief Sets the result_on_layer for all parents based on the SupportElement supplied.
     *
     * \param elem[in] The SupportElements, which parent's position should be determined.
     */
    void setPointsOnAreas(const SupportElement* elem);
    /*!
     * \brief Get the best point to connect to the model and set the result_on_layer of the relevant SupportElement accordingly.
     *
     * \param move_bounds[in,out] All currently existing influence areas
     * \param first_elem[in,out] SupportElement that did not have its result_on_layer set meaning that it does not have a child element.
     * \param layer_idx[in] The current layer.
     * \return Should elem be deleted.
     */
    bool setToModelContact(std::vector<std::set<SupportElement*>>& move_bounds, SupportElement* first_elem, const LayerIndex layer_idx);

    /*!
     * \brief Set the result_on_layer point for all influence areas
     *
     * \param move_bounds[in,out] All currently existing influence areas
     */
    void createNodesFromArea(std::vector<std::set<SupportElement*>>& move_bounds);

    /*!
     * \brief Draws circles around result_on_layer points of the influence areas
     *
     * \param linear_data[in] All currently existing influence areas with the layer they are on
     * \param layer_tree_polygons[out] Resulting branch areas with the layerindex they appear on. layer_tree_polygons.size() has to be at least linear_data.size() as each Influence area in linear_data will save have at least one (that's why it's a vector<vector>) corresponding branch area in layer_tree_polygons.
     * \param inverse_tree_order[in] A mapping that returns the child of every influence area.
     */
    void generateBranchAreas(std::vector<std::pair<LayerIndex, SupportElement*>>& linear_data, std::vector<std::unordered_map<SupportElement*, Polygons>>& layer_tree_polygons, const std::map<SupportElement*, SupportElement*>& inverse_tree_order);

    /*!
     * \brief Applies some smoothing to the outer wall, intended to smooth out sudden jumps as they can happen when a branch moves though a hole.
     *
     * \param layer_tree_polygons[in,out] Resulting branch areas with the layerindex they appear on.
     */
    void smoothBranchAreas(std::vector<std::unordered_map<SupportElement*, Polygons>>& layer_tree_polygons);

    /*!
     * \brief Drop down areas that do rest non-gracefully on the model to ensure the branch actually rests on something.
     *
     * \param layer_tree_polygons[in] Resulting branch areas with the layerindex they appear on.
     * \param linear_data[in] All currently existing influence areas with the layer they are on
     * \param dropped_down_areas[out] Areas that have to be added to support all non-graceful areas.
     * \param inverse_tree_order[in] A mapping that returns the child of every influence area.
     */
    void dropNonGraciousAreas(std::vector<std::unordered_map<SupportElement*, Polygons>>& layer_tree_polygons, const std::vector<std::pair<LayerIndex, SupportElement*>>& linear_data, std::vector<std::vector<std::pair<LayerIndex, Polygons>>>& dropped_down_areas, const std::map<SupportElement*, SupportElement*>& inverse_tree_order);

    /*!
     * \brief Generates Support Floor, ensures Support Roof can not cut of branches, and saves the branches as support to storage
     *
     * \param support_layer_storage[in] Areas where support should be generated.
     * \param support_roof_storage[in] Areas where support was replaced with roof.
     * \param storage[in,out] The storage where the support should be stored.
     */
    void finalizeInterfaceAndSupportAreas(
        const PrintObject               &print_object,
        std::vector<Polygons>           &support_layer_storage,
        std::vector<Polygons>           &support_roof_storage,

        SupportGeneratorLayersPtr   	&bottom_contacts,
        SupportGeneratorLayersPtr   	&top_contacts,
        SupportGeneratorLayersPtr       &intermediate_layers,
        SupportGeneratorLayerStorage    &layer_storage);

    /*!
     * \brief Draws circles around result_on_layer points of the influence areas and applies some post processing.
     *
     * \param move_bounds[in] All currently existing influence areas
     * \param storage[in,out] The storage where the support should be stored.
     */
    void drawAreas(
        PrintObject                             &print_object, 
        std::vector<std::set<SupportElement*>>  &move_bounds,

        SupportGeneratorLayersPtr            	&bottom_contacts,
        SupportGeneratorLayersPtr   	        &top_contacts,
        SupportGeneratorLayersPtr               &intermediate_layers,
        SupportGeneratorLayerStorage            &layer_storage);

    /*!
     * \brief Settings with the indexes of meshes that use these settings.
     *
     */
    std::vector<std::pair<TreeSupportSettings, std::vector<size_t>>> m_grouped_meshes;

    /*!
     * \brief Generator for model collision, avoidance and internal guide volumes.
     *
     */
    TreeModelVolumes m_volumes;

    /*!
     * \brief Contains config settings to avoid loading them in every function. This was done to improve readability of the code.
     */
    TreeSupportSettings m_config;

    /*!
     * \brief The progress multiplier of all values added progress bar.
     * Required for the progress bar the behave as expected when areas have to be calculated multiple times
     */
    double m_progress_multiplier = 1;

    /*!
     * \brief The progress offset added to all values communicated to the progress bar.
     * Required for the progress bar the behave as expected when areas have to be calculated multiple times
     */
    double m_progress_offset = 0;
};


} // namespace Slic3r

namespace std
{
template <>
struct hash<Slic3r::TreeSupport::SupportElement>
{
    size_t operator()(const Slic3r::TreeSupport::SupportElement& node) const
    {
        size_t hash_node = Slic3r::PointHash{}(node.target_position);
        boost::hash_combine(hash_node, size_t(node.target_height));
        return hash_node;
    }
};
} // namespace std

#endif /* slic3r_TreeSupport_hpp */
