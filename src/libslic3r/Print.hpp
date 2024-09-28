///|/ Copyright (c) Prusa Research 2016 - 2023 Lukáš Matěna @lukasmatena, Enrico Turri @enricoturri1966, Vojtěch Bubník @bubnikv, Tomáš Mészáros @tamasmeszaros, Pavel Mikuš @Godrak, Lukáš Hejl @hejllukas, Filip Sykala @Jony01, Oleksandra Iushchenko @YuSanka, Vojtěch Král @vojtechkral
///|/ Copyright (c) BambuStudio 2023 manch1n @manch1n
///|/ Copyright (c) SuperSlicer 2022 Remi Durand @supermerill
///|/ Copyright (c) 2019 Bryan Smith
///|/ Copyright (c) 2017 Eyal Soha @eyal0
///|/ Copyright (c) Slic3r 2013 - 2016 Alessandro Ranellucci @alranel
///|/ Copyright (c) 2017 Joseph Lenox @lordofhyphens
///|/
///|/ ported from lib/Slic3r/Print.pm:
///|/ Copyright (c) Prusa Research 2016 - 2018 Vojtěch Bubník @bubnikv, Tomáš Mészáros @tamasmeszaros
///|/ Copyright (c) Slic3r 2011 - 2016 Alessandro Ranellucci @alranel
///|/ Copyright (c) 2012 - 2013 Mark Hindess
///|/ Copyright (c) 2013 Devin Grady
///|/ Copyright (c) 2012 - 2013 Mike Sheldrake @mesheldrake
///|/ Copyright (c) 2012 Henrik Brix Andersen @henrikbrixandersen
///|/ Copyright (c) 2012 Michael Moon
///|/ Copyright (c) 2011 Richard Goodwin
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_Print_hpp_
#define slic3r_Print_hpp_

#include "Fill/FillAdaptive.hpp"
#include "Fill/FillLightning.hpp"
#include "PrintBase.hpp"

#include "BoundingBox.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Flow.hpp"
#include "Point.hpp"
#include "Slicing.hpp"
#include "SupportSpotsGenerator.hpp"
#include "TriangleMeshSlicer.hpp"
#include "Surface.hpp"
#include "GCode/ToolOrdering.hpp"
#include "GCode/WipeTower.hpp"
#include "GCode/ThumbnailData.hpp"
#include "MultiMaterialSegmentation.hpp"

#include "libslic3r.h"

#include <Eigen/Geometry>

#include <atomic>
#include <ctime>
#include <functional>
#include <optional>
#include <set>
#include <tcbspan/span.hpp>

namespace Slic3r {

class GCodeGenerator;
struct GCodeProcessorResult;
class Layer;
class ModelObject;
class Print;
class PrintObject;
class SupportLayer;

namespace FillAdaptive {
    struct Octree;
    struct OctreeDeleter;
    using OctreePtr = std::unique_ptr<Octree, OctreeDeleter>;
}; // namespace FillAdaptive

namespace FillLightning {
    class Generator;
    struct GeneratorDeleter;
    using GeneratorPtr = std::unique_ptr<Generator, GeneratorDeleter>;
}; // namespace FillLightning


// Print step IDs for keeping track of the print state.
// The Print steps are applied in this order.
enum PrintStep : uint8_t {
    psWipeTower,
    // Ordering of the tools on PrintObjects for a multi-material print.
    // psToolOrdering is a synonym to psWipeTower, as the Wipe Tower calculates and modifies the ToolOrdering,
    // while if printing without the Wipe Tower, the ToolOrdering is calculated as well.
    psToolOrdering = psWipeTower,
    psAlertWhenSupportsNeeded,
    psSkirtBrim,
    // Last step before G-code export, after this step is finished, the initial extrusion path preview
    // should be refreshed.
    psSlicingFinished = psSkirtBrim,
    psGCodeExport,
   //TODO: psGCodeLoader (for params that are only used for time display and such)
    psCount,
};

enum PrintObjectStep : uint8_t {
    posSlice,
    posPerimeters,
    posPrepareInfill,
    posInfill,
    posIroning,
    posSupportSpotsSearch,
    posSupportMaterial, 
    posEstimateCurledExtrusions,
    posCalculateOverhangingPerimeters,
    posSimplifyPath, // simplify &  arc fitting from BBS
    posCount,
};

/**
* order:
*            m_objects[idx]->make_perimeters();
*                   -> slice()
*                   -> make_perimeters()
*            m_objects[idx]->infill();
*            m_objects[idx]->ironing();
*            obj->generate_support_spots();
*            psAlertWhenSupportsNeeded
*            obj.generate_support_material();
*            obj.estimate_curled_extrusions();
*            obj.calculate_overhanging_perimeters();
*            _make_wipe_tower();
*            _make_skirt();
*            make_brim();
*            simplify_extrusion_path();
* 
*           then export_gcode();
* */

// step % for starting this step
inline std::map<PrintObjectStep, int> objectstep_2_percent = {{PrintObjectStep::posSlice, 0},
                                                       {PrintObjectStep::posPerimeters, 10},
                                                       {PrintObjectStep::posPrepareInfill, 20},
                                                       {PrintObjectStep::posInfill, 30},
                                                       {PrintObjectStep::posIroning, 40},
                                                       {PrintObjectStep::posSupportSpotsSearch, 45},
                                                       {PrintObjectStep::posSupportMaterial, 50},
                                                       {PrintObjectStep::posEstimateCurledExtrusions, 60},
                                                       {PrintObjectStep::posCalculateOverhangingPerimeters, 65},
                                                       {PrintObjectStep::posSimplifyPath, 80},
                                                       {PrintObjectStep::posCount, 85}};

inline std::map<PrintStep, int> printstep_2_percent = {
    {PrintStep::psAlertWhenSupportsNeeded, 45},
    {PrintStep::psSkirtBrim, 70},
    {PrintStep::psWipeTower, 75},
    {PrintStep::psGCodeExport, 85},
    {PrintStep::psCount, 100},
};
// A PrintRegion object represents a group of volumes to print
// sharing the same config (including the same assigned extruder(s))
class PrintRegion
{
public:
    PrintRegion() = default;
    PrintRegion(const PrintRegionConfig &config);
    PrintRegion(const PrintRegionConfig &config, const size_t config_hash, int print_object_region_id = -1) : m_config(config), m_config_hash(config_hash), m_print_object_region_id(print_object_region_id) {}
    PrintRegion(PrintRegionConfig &&config);
    PrintRegion(PrintRegionConfig &&config, const size_t config_hash, int print_object_region_id = -1) : m_config(std::move(config)), m_config_hash(config_hash), m_print_object_region_id(print_object_region_id) {}
    ~PrintRegion() = default;

// Methods NOT modifying the PrintRegion's state:
public:
    const PrintRegionConfig&    config() const throw() { return m_config; }
    size_t                      config_hash() const throw() { return m_config_hash; }
    // Identifier of this PrintRegion in the list of Print::m_print_regions.
    int                         print_region_id() const throw() { return m_print_region_id; }
    int                         print_object_region_id() const throw() { return m_print_object_region_id; }
	// 1-based extruder identifier for this region and role.
    uint16_t 				    extruder(FlowRole role, const PrintObject& object) const;
    Flow                        flow(const PrintObject &object, FlowRole role, double layer_height, size_t layer_id) const;
    float                       width(FlowRole role, bool first_layer, const PrintObject& object) const;
    // Average diameter of nozzles participating on extruding this region.
    coordf_t                    nozzle_dmr_avg(const PrintConfig &print_config) const;

    // Collect 0-based extruder indices used to print this region's object.
	void                        collect_object_printing_extruders(const Print& print, std::set<uint16_t> &object_extruders) const;
	static void                 collect_object_printing_extruders(const PrintConfig &print_config, const PrintObjectConfig &object_config, const PrintRegionConfig &region_config, std::set<uint16_t> &object_extruders);

// Methods modifying the PrintRegion's state:
public:
    void                        set_config(const PrintRegionConfig &config) { m_config = config; m_config_hash = m_config.hash(); }
    void                        set_config(PrintRegionConfig &&config) { m_config = std::move(config); m_config_hash = m_config.hash(); }
    void                        config_apply_only(const ConfigBase &other, const t_config_option_keys &keys, bool ignore_nonexistent = false) 
                                        { m_config.apply_only(other, keys, ignore_nonexistent); m_config_hash = m_config.hash(); }
private:
    friend Print;
    friend void print_region_ref_inc(PrintRegion&);
    friend void print_region_ref_reset(PrintRegion&);
    friend int  print_region_ref_cnt(const PrintRegion&);

    PrintRegionConfig  m_config;
    size_t             m_config_hash;
    int                m_print_region_id { -1 };
    int                m_print_object_region_id { -1 };
    int                m_ref_cnt { 0 };
};

inline bool operator==(const PrintRegion &lhs, const PrintRegion &rhs) { return lhs.config_hash() == rhs.config_hash() && lhs.config() == rhs.config(); }
inline bool operator!=(const PrintRegion &lhs, const PrintRegion &rhs) { return ! (lhs == rhs); }

// For const correctness: Wrapping a vector of non-const pointers as a span of const pointers.
template<class T>
using SpanOfConstPtrs           = tcb::span<const T* const>;

using LayerPtrs                 = std::vector<Layer*>;
using SupportLayerPtrs          = std::vector<SupportLayer*>;

class BoundingBoxf3;        // TODO: for temporary constructor parameter

// Single instance of a PrintObject.
// As multiple PrintObjects may be generated for a single ModelObject (their instances differ in rotation around Z),
// ModelObject's instancess will be distributed among these multiple PrintObjects.
struct PrintInstance
{
    // Parent PrintObject
    PrintObject 		*print_object;
    // Source ModelInstance of a ModelObject, for which this print_object was created.
	const ModelInstance *model_instance;
	// Shift of this instance's center into the world coordinates.
	Point 				 shift;
};

typedef std::vector<PrintInstance> PrintInstances;

class PrintObjectRegions
{
public:
    // Bounding box of a ModelVolume transformed into the working space of a PrintObject, possibly
    // clipped by a layer range modifier.
    // Only Eigen types of Nx16 size are vectorized. This bounding box will not be vectorized.
    static_assert(sizeof(Eigen::AlignedBox<float, 3>) == 24, "Eigen::AlignedBox<float, 3> is not being vectorized, thus it does not need to be aligned");
    using BoundingBox = Eigen::AlignedBox<float, 3>;
    struct VolumeExtents {
        ObjectID             volume_id;
        BoundingBox          bbox;
    };

    struct VolumeRegion
    {
        // ID of the associated ModelVolume.
        const ModelVolume   *model_volume { nullptr };
        // Index of a parent VolumeRegion.
        int                  parent { -1 };
        // Pointer to PrintObjectRegions::all_regions, null for a negative volume.
        PrintRegion         *region { nullptr };
        // Pointer to VolumeExtents::bbox.
        const BoundingBox   *bbox { nullptr };
        // To speed up merging of same regions.
        const VolumeRegion  *prev_same_region { nullptr };
    };

    struct PaintedRegion
    {
        // 1-based extruder identifier.
        unsigned int     extruder_id;
        // Index of a parent VolumeRegion.
        int              parent { -1 };
        // Pointer to PrintObjectRegions::all_regions.
        PrintRegion     *region { nullptr };
    };

    // One slice over the PrintObject (possibly the whole PrintObject) and a list of ModelVolumes and their bounding boxes
    // possibly clipped by the layer_height_range.
    struct LayerRangeRegions
    {
        t_layer_height_range        layer_height_range;
        // Config of the layer range, null if there is just a single range with no config override.
        // Config is owned by the associated ModelObject.
        const DynamicPrintConfig*   config { nullptr };
        // Volumes sorted by ModelVolume::id().
        std::vector<VolumeExtents>  volumes;

        // Sorted in the order of their source ModelVolumes, thus reflecting the order of region clipping, modifier overrides etc.
        std::vector<VolumeRegion>   volume_regions;
        std::vector<PaintedRegion>  painted_regions;

        bool has_volume(const ObjectID id) const {
            auto it = lower_bound_by_predicate(this->volumes.begin(), this->volumes.end(), [id](const VolumeExtents &l) { return l.volume_id < id; });
            return it != this->volumes.end() && it->volume_id == id;
        }
    };

    struct GeneratedSupportPoints{
        Transform3d object_transform; // for frontend object mapping
        SupportSpotsGenerator::SupportPoints support_points;
        SupportSpotsGenerator::PartialObjects partial_objects;
    };

    std::vector<std::unique_ptr<PrintRegion>>   all_regions;
    std::vector<LayerRangeRegions>              layer_ranges;
    // Transformation of this ModelObject into one of the associated PrintObjects (all PrintObjects derived from a single modelObject differ by a Z rotation only).
    // This transformation is used to calculate VolumeExtents.
    Transform3d                                 trafo_bboxes;
    std::vector<ObjectID>                       cached_volume_ids;

    std::optional<GeneratedSupportPoints> generated_support_points;

    void ref_cnt_inc() { ++ m_ref_cnt; }
    void ref_cnt_dec() { if (-- m_ref_cnt == 0) delete this; }
    void clear() {
        all_regions.clear();
        layer_ranges.clear();
        cached_volume_ids.clear();
    }

private:
    friend class PrintObject;
    // Number of PrintObjects generated from the same ModelObject and sharing the regions.
    // ref_cnt could only be modified by the main thread, thus it does not need to be atomic.
    size_t                                      m_ref_cnt{ 0 };
};

class PrintObject : public PrintObjectBaseWithState<Print, PrintObjectStep, posCount>
{
private: // Prevents erroneous use by other classes.
    typedef PrintObjectBaseWithState<Print, PrintObjectStep, posCount> Inherited;

public:
    // Size of an object: XYZ in scaled coordinates. The size might not be quite snug in XY plane.
    const Vec3crd&               size() const			{ return m_size; }
    const PrintObjectConfig&     config() const         { return m_config; }
    const PrintRegionConfig&     default_region_config(const PrintRegionConfig &from_print) const;
    auto                         layers() const         { return SpanOfConstPtrs<Layer>(const_cast<const Layer* const* const>(m_layers.data()), m_layers.size()); }
    auto                         support_layers() const { return SpanOfConstPtrs<SupportLayer>(const_cast<const SupportLayer* const* const>(m_support_layers.data()), m_support_layers.size()); }
    const Transform3d&           trafo() const          { return m_trafo; }
    // Trafo with the center_offset() applied after the transformation, to center the object in XY before slicing.
    Transform3d                  trafo_centered() const 
        { Transform3d t = this->trafo(); t.pretranslate(Vec3d(- unscale<double>(m_center_offset.x()), - unscale<double>(m_center_offset.y()), 0)); return t; }
    const PrintInstances&        instances() const      { return m_instances; }

    // Whoever will get a non-const pointer to PrintObject will be able to modify its layers.
    LayerPtrs&                   layers()               { return m_layers; }
    SupportLayerPtrs&            edit_support_layers()       { return m_support_layers; }

    // Bounding box is used to align the object infill patterns, and to calculate attractor for the rear seam.
    // The bounding box may not be quite snug.
    BoundingBox                  bounding_box() const   { return BoundingBox(Point(- m_size.x() / 2, - m_size.y() / 2), Point(m_size.x() / 2, m_size.y() / 2)); }
    // Height is used for slicing, for sorting the objects by height for sequential printing and for checking vertical clearence in sequential print mode.
    // The height is snug.
    coord_t 				     height() const         { return m_size.z(); }
    // Centering offset of the sliced mesh from the scaled and rotated mesh of the model.
    const Point& 			     center_offset() const  { return m_center_offset; }

    bool                         has_brim() const       {
        return (this->config().brim_width.value > 0 && this->config().brim_width_interior.value > 0)
            && ! this->has_raft();
    }

    // This is the *total* layer count (including support layers)
    // this value is not supposed to be compared with Layer::id
    // since they have different semantics.
    size_t 			total_layer_count() const { return this->layer_count() + this->support_layer_count(); }
    size_t 			layer_count() const { return m_layers.size(); }
    void 			clear_layers();
    const Layer* 	get_layer(int idx) const { return m_layers[idx]; }
    Layer* 			get_layer(int idx) 		 { return m_layers[idx]; }
    // Get a layer exactly at print_z.
    const Layer*	get_layer_at_printz(coordf_t print_z) const;
    Layer*			get_layer_at_printz(coordf_t print_z);
    // Get a layer approximately at print_z.
    const Layer*	get_layer_at_printz(coordf_t print_z, coordf_t epsilon) const;
    Layer*			get_layer_at_printz(coordf_t print_z, coordf_t epsilon);
    // Get the first layer approximately bellow print_z.
    const Layer*	get_first_layer_bellow_printz(coordf_t print_z, coordf_t epsilon) const;
    // For sparse infill, get the max spasing avaialable in this object (avaialable after prepare_infill)
    coord_t         get_sparse_max_spacing() const { return m_max_sparse_spacing; }

    // print_z: top of the layer; slice_z: center of the layer.
    Layer*          add_layer(int id, coordf_t height, coordf_t print_z, coordf_t slice_z);

    size_t          support_layer_count() const { return m_support_layers.size(); }
    void            clear_support_layers();
    const SupportLayer*   get_support_layer(int idx) { return m_support_layers[idx]; }
    void            add_support_layer(int id, int interface_id, coordf_t height, coordf_t print_z);
    SupportLayerPtrs::iterator insert_support_layer(SupportLayerPtrs::const_iterator pos, size_t id, size_t interface_id, coordf_t height, coordf_t print_z, coordf_t slice_z);
    //void            delete_support_layer(int idx);
    
    // Initialize the layer_height_profile from the model_object's layer_height_profile, from model_object's layer height table, or from slicing parameters.
    // Returns true, if the layer_height_profile was changed.
    static bool     update_layer_height_profile(const ModelObject &model_object, const SlicingParameters &slicing_parameters, std::vector<coordf_t> &layer_height_profile);

    // Collect the slicing parameters, to be used by variable layer thickness algorithm,
    // by the interactive layer height editor and by the printing process itself.
    // The slicing parameters are dependent on various configuration values
    // (layer height, first layer height, raft settings, print nozzle diameter etc).
    const SlicingParameters&                    slicing_parameters() const { return *m_slicing_params; }
    static std::shared_ptr<SlicingParameters>   slicing_parameters(const DynamicPrintConfig &full_config, const ModelObject &model_object, float object_max_z);

    size_t                      num_printing_regions() const throw() { assert(m_shared_regions); return m_shared_regions->all_regions.size(); }
    const PrintRegion&          printing_region(size_t idx) const throw() { assert(m_shared_regions); return *m_shared_regions->all_regions[idx].get(); }
    //FIXME returing all possible regions before slicing, thus some of the regions may not be slicing at the end.
    std::vector<std::reference_wrapper<const PrintRegion>> all_regions() const;
    const PrintObjectRegions*   shared_regions() const throw() { return m_shared_regions; }

    bool                        has_support()           const { return m_config.support_material || m_config.support_material_enforce_layers > 0; }
    bool                        has_raft()              const { return m_config.raft_layers > 0; }
    bool                        has_support_material()  const { return this->has_support() || this->has_raft(); }
    // Checks if the model object is painted using the multi-material painting gizmo.
    bool                        is_mm_painted()         const { return this->model_object()->is_mm_painted(); }

    // returns 0-based indices of extruders used to print the object (without brim, support and other helper extrusions)
    std::set<uint16_t>   object_extruders() const;
    double               get_first_layer_height() const;

    // Called by make_perimeters()
    void slice();

    // Helpers to slice support enforcer / blocker meshes by the support generator.
    std::vector<Polygons>       slice_support_volumes(const ModelVolumeType model_volume_type) const;
    std::vector<Polygons>       slice_support_blockers() const { return this->slice_support_volumes(ModelVolumeType::SUPPORT_BLOCKER); }
    std::vector<Polygons>       slice_support_enforcers() const { return this->slice_support_volumes(ModelVolumeType::SUPPORT_ENFORCER); }

    // Helpers to project custom facets on slices
    void project_and_append_custom_facets(bool seam, EnforcerBlockerType type, std::vector<Polygons>& expolys) const;

    /// skirts if done per copy and not per platter
    const std::optional<ExtrusionEntityCollection>& skirt_first_layer() const { return m_skirt_first_layer; }
    const ExtrusionEntityCollection& skirt() const { return m_skirt; }
    const ExtrusionEntityCollection& brim() const { return m_brim; }

protected:
    // to be called from Print only.
    friend class Print;
    friend class PrintBaseWithState<PrintStep, psCount>;

	PrintObject(Print* print, ModelObject* model_object, const Transform3d& trafo, PrintInstances&& instances);
    ~PrintObject() override {
        if (m_shared_regions && --m_shared_regions->m_ref_cnt == 0)
            delete m_shared_regions;
        clear_layers();
        clear_support_layers();
    }

    void                    config_apply(const ConfigBase &other, bool ignore_nonexistent = false) { m_config.apply(other, ignore_nonexistent); }
    void                    config_apply_only(const ConfigBase &other, const t_config_option_keys &keys, bool ignore_nonexistent = false) { m_config.apply_only(other, keys, ignore_nonexistent); }
    PrintBase::ApplyStatus  set_instances(PrintInstances &&instances);
    // Invalidates the step, and its depending steps in PrintObject and Print.
    bool                    invalidate_step(PrintObjectStep step);
    // Invalidates all PrintObject and Print steps.
    bool                    invalidate_all_steps();
    // Invalidate steps based on a set of parameters changed.
    // It may be called for both the PrintObjectConfig and PrintRegionConfig.
    bool                    invalidate_state_by_config_options(
        const ConfigOptionResolver &old_config, const ConfigOptionResolver &new_config, const std::vector<t_config_option_key> &opt_keys);
    // If ! m_slicing_params.valid, recalculate.
    void                    update_slicing_parameters();

    // Called on main thread with stopped or paused background processing to let PrintObject release data for its milestones that were invalidated or canceled.
    void                    cleanup();

    static PrintObjectConfig object_config_from_model_object(const PrintObjectConfig &default_object_config, const ModelObject &object, size_t num_extruders);

private:
    void make_perimeters();
    void prepare_infill();
    void clear_fills();
    void infill();
    void ironing();
    void generate_support_spots();
    void generate_support_material();
    void estimate_curled_extrusions();
    void calculate_overhanging_perimeters();
    void simplify_extrusion_path();

    void slice_volumes();
    // Has any support (not counting the raft).
    ExPolygons _shrink_contour_holes(double contour_delta, double default_delta, double convex_delta, const ExPolygons& input) const;
    void _transform_hole_to_polyholes();
    void _min_overhang_threshold();
    ExPolygons _smooth_curves(const ExPolygons &input, const PrintRegionConfig &conf) const;
    void detect_surfaces_type();
    void apply_solid_infill_below_layer_area();
    void process_external_surfaces();
    void discover_vertical_shells();
    void bridge_over_infill();
    void replaceSurfaceType(SurfaceType st_to_replace, SurfaceType st_replacement, SurfaceType st_under_it);
    // void clip_fill_surfaces(); //infill_only_where_needed
    void tag_under_bridge();
    void discover_horizontal_shells();
    void clean_surfaces();
    void combine_infill();
    void _generate_support_material();
    void _compute_max_sparse_spacing();
    std::pair<FillAdaptive::OctreePtr, FillAdaptive::OctreePtr> prepare_adaptive_infill_data(
        const std::vector<std::pair<const Surface*, float>>& surfaces_w_bottom_z) const;
    FillLightning::GeneratorPtr prepare_lightning_infill_data();

    // XYZ in scaled coordinates
    Vec3crd									m_size;
    PrintObjectConfig                       m_config;
    // Translation in Z + Rotation + Scaling / Mirroring.
    Transform3d                             m_trafo = Transform3d::Identity();
    // Slic3r::Point objects in scaled G-code coordinates
    std::vector<PrintInstance>              m_instances;
    // The mesh is being centered before thrown to Clipper, so that the Clipper's fixed coordinates require less bits.
    // This is the adjustment of the  the Object's coordinate system towards PrintObject's coordinate system.
    Point                                   m_center_offset;

    // Object split into layer ranges and regions with their associated configurations.
    // Shared among PrintObjects created for the same ModelObject.
    PrintObjectRegions                     *m_shared_regions { nullptr };

    std::shared_ptr<SlicingParameters>      m_slicing_params;
    LayerPtrs                               m_layers;
    SupportLayerPtrs                        m_support_layers;

    // Ordered collections of extrusion paths to build skirt loops and brim.
    // have to be duplicated per copy
    std::optional<ExtrusionEntityCollection> m_skirt_first_layer;
    ExtrusionEntityCollection               m_skirt;
    ExtrusionEntityCollection               m_brim;

    // this is set to true when LayerRegion->slices is split in top/internal/bottom
    // so that next call to make_perimeters() performs a union() before computing loops
    bool                    				m_typed_slices = false;

    //this setting allow fill_aligned_z to get the max sparse spacing spacing.
    coord_t                                 m_max_sparse_spacing = 0;

    // pair < adaptive , support>, filled by prepare_adaptive_infill_data() (in bridge_over_infill() in prepare_infill()) and used in infill()
    std::pair<FillAdaptive::OctreePtr, FillAdaptive::OctreePtr> m_adaptive_fill_octrees;
    // filled by prepare_lightning_infill_data() (in bridge_over_infill() in prepare_infill()) and used in infill()
    FillLightning::GeneratorPtr m_lightning_generator;
};



struct WipeTowerData
{
    // Following section will be consumed by the GCodeGenerator.
    // Tool ordering of a non-sequential print has to be known to calculate the wipe tower.
    // Cache it here, so it does not need to be recalculated during the G-code generation.
    ToolOrdering                                         &tool_ordering;
    // Cache of tool changes per print layer.
    std::unique_ptr<std::vector<WipeTower::ToolChangeResult>> priming;
    std::vector<std::vector<WipeTower::ToolChangeResult>> tool_changes;
    std::unique_ptr<WipeTower::ToolChangeResult>          final_purge;
    std::vector<std::pair<float, std::vector<float>>>     used_filament_until_layer;
    int                                                   number_of_toolchanges;

    // Depth of the wipe tower to pass to GLCanvas3D for exact bounding box:
    float                                                 depth;
    std::vector<std::pair<float, float>>                  z_and_depth_pairs;
    float                                                 brim_width;
    float                                                 height;

    // Data needed to generate fake extrusions for conflict checking.
    float                                                 width;
    float                                                 first_layer_height;
    float                                                 cone_angle;
    Vec2d                                                 position;
    float                                                 rotation_angle;

    void clear() {
        priming.reset(nullptr);
        tool_changes.clear();
        final_purge.reset(nullptr);
        used_filament_until_layer.clear();
        number_of_toolchanges = -1;
        depth = 0.f;
        z_and_depth_pairs.clear();
        brim_width = 0.f;
        height = 0.f;
        width = 0.f;
        first_layer_height = 0.f;
        cone_angle = 0.f;
        position = Vec2d::Zero();
        rotation_angle = 0.f;
    }

private:
	// Only allow the WipeTowerData to be instantiated internally by Print, 
	// as this WipeTowerData shares reference to Print::m_tool_ordering.
	friend class Print;
	WipeTowerData(ToolOrdering &tool_ordering) : tool_ordering(tool_ordering) { clear(); }
	WipeTowerData(const WipeTowerData & /* rhs */) = delete;
	WipeTowerData &operator=(const WipeTowerData & /* rhs */) = delete;
};

struct PrintStatistics
{
    PrintStatistics() { clear(); }
    // PrintEstimatedStatistics::ETimeMode::Normal -> time
    std::map<uint8_t, double>       estimated_print_time;
    std::map<uint8_t, std::string>  estimated_print_time_str;
    double                          total_used_filament;
    std::vector<std::pair<size_t, double>> color_extruderid_to_used_filament; // id -> mm (length)
    double                          total_extruded_volume;
    double                          total_cost;
    int                             total_toolchanges;
    double                          total_weight;
    std::vector<std::pair<size_t, double>> color_extruderid_to_used_weight;
    double                          total_wipe_tower_cost;
    double                          total_wipe_tower_filament;
    double                          total_wipe_tower_filament_weight;
    std::vector<unsigned int>       printing_extruders;
    unsigned int                    initial_extruder_id;
    std::string                     initial_filament_type;
    std::string                     printing_filament_types;
    std::map<size_t, double>        filament_stats; // extruder id -> volume in mm3
    std::vector<std::pair<double, float>> layer_area_stats; // print_z to area

    std::atomic_bool is_computing_gcode;

    // Config with the filled in print statistics.
    DynamicConfig           config() const;
    // Config with the statistics keys populated with placeholder strings.
    static DynamicConfig    placeholders();
    // Replace the print statistics placeholders in the path.
    std::string             finalize_output_path(const std::string &path_in) const;

    void clear() {
        total_used_filament    = 0.;
        total_extruded_volume  = 0.;
        total_cost             = 0.;
        total_toolchanges      = 0;
        total_weight           = 0.;
        total_wipe_tower_cost  = 0.;
        total_wipe_tower_filament = 0.;
        total_wipe_tower_filament_weight = 0.;
        initial_extruder_id    = 0;
        initial_filament_type.clear();
        printing_filament_types.clear();
        filament_stats.clear();
        printing_extruders.clear();
        is_computing_gcode = false;
    }

    static const std::string FilamentUsedG;
    static const std::string FilamentUsedGMask;
    static const std::string TotalFilamentUsedG;
    static const std::string TotalFilamentUsedGMask;
    static const std::string TotalFilamentUsedGValueMask;
    static const std::string FilamentUsedCm3;
    static const std::string FilamentUsedCm3Mask;
    static const std::string FilamentUsedMm;
    static const std::string FilamentUsedMmMask;
    static const std::string FilamentCost;
    static const std::string FilamentCostMask;
    static const std::string TotalFilamentCost;
    static const std::string TotalFilamentCostMask;
    static const std::string TotalFilamentCostValueMask;
    static const std::string TotalFilamentUsedWipeTower;
    static const std::string TotalFilamentUsedWipeTowerValueMask;
};

struct ConflictResult
{
    std::string _objName1;
    std::string _objName2;
    double      _height;
    const void* _obj1; // nullptr means wipe tower
    const void* _obj2;
    int         layer = -1;
    ConflictResult(const std::string& objName1, const std::string& objName2, double height, const void* obj1, const void* obj2)
        : _objName1(objName1), _objName2(objName2), _height(height), _obj1(obj1), _obj2(obj2)
    {}
    ConflictResult() = default;
};

using ConflictResultOpt = std::optional<ConflictResult>;

class BrimLoop {
public:
    BrimLoop(const Polygon& p) : lines(Polylines{ p.split_at_first_point() }), is_loop(true) {}
    BrimLoop(const Polyline& l) : lines(Polylines{l}), is_loop(false) {}
    Polylines lines;
    std::vector<BrimLoop> children;
    bool is_loop; // has only one polyline stored and front == back
    Polygon polygon() const{
        assert(is_loop);
        Polygon poly = Polygon(lines.front().points);
        if (poly.points.front() == poly.points.back())
            poly.points.resize(poly.points.size() - 1);
        return poly;
    }
};

using PrintObjectPtrs          = std::vector<PrintObject*>;
using ConstPrintObjectPtrs     = std::vector<const PrintObject*>;

using PrintRegionPtrs          = std::vector<PrintRegion*>;

// The complete print tray with possibly multiple objects.
class Print : public PrintBaseWithState<PrintStep, psCount>
{
private: // Prevents erroneous use by other classes.
    typedef PrintBaseWithState<PrintStep, psCount> Inherited;
    // Bool indicates if supports of PrintObject are top-level contour.
    typedef std::pair<PrintObject *, bool>         PrintObjectInfo;

public:
    Print() {
        //create config hierachy
        m_default_object_config.parent = &m_config;
        m_default_region_config.parent = &m_default_object_config;
    };
	virtual ~Print() { this->clear(); }

	PrinterTechnology	technology() const noexcept override { return ptFFF; }

    // Methods, which change the state of Print / PrintObject / PrintRegion.
    // The following methods are synchronized with process() and export_gcode(),
    // so that process() and export_gcode() may be called from a background thread.
    // In case the following methods need to modify data processed by process() or export_gcode(),
    // a cancellation callback is executed to stop the background processing before the operation.
    void                clear() override;
    bool                empty() const override { return m_objects.empty(); }
    // List of existing PrintObject IDs, to remove notifications for non-existent IDs.
    std::vector<ObjectID> print_object_ids() const override;

    ApplyStatus         apply(const Model &model, DynamicPrintConfig config) override;
    void                set_task(const TaskParams &params) override { PrintBaseWithState<PrintStep, psCount>::set_task_impl(params, m_objects); }
    void                process() override;
    void                finalize() override { PrintBaseWithState<PrintStep, psCount>::finalize_impl(m_objects); }
    void                cleanup() override;

    // Exports G-code into a file name based on the path_template, returns the file path of the generated G-code file.
    // If preview_data is not null, the preview_data is filled in for the G-code visualization (not used by the command line Slic3r).
    std::string         export_gcode(const std::string& path_template, GCodeProcessorResult* result, ThumbnailsGeneratorCallback thumbnail_cb = nullptr);

    // methods for handling state
    bool                is_step_done(PrintStep step) const { return Inherited::is_step_done(step); }
    // Returns true if an object step is done on all objects and there's at least one object.    
    bool                is_step_done(PrintObjectStep step) const;
    // Returns true if the last step was finished with success.
    bool                finished() const override { return this->is_step_done(psGCodeExport); }

    bool                has_infinite_skirt() const;
    bool                has_skirt() const;
    bool                has_brim() const;

    // Returns an empty string if valid, otherwise returns an error message.
    std::pair<PrintValidationError, std::string> validate(std::vector<std::string>* warnings = nullptr) const override;
    Flow                brim_flow(size_t extruder_id, const PrintObjectConfig &brim_config) const;
    Flow                skirt_flow(size_t extruder_id, bool first_layer=false) const;
    double              get_min_first_layer_height() const;
    double              get_object_first_layer_height(const PrintObject& object) const;

    // get the extruders of these obejcts
    std::set<uint16_t>  object_extruders(const PrintObjectPtrs &objects) const;
    // get all extruders from the list of objects in this print ( same as print.object_extruders(print.objects()) )
    std::set<uint16_t>  object_extruders() const;
    std::set<uint16_t>  support_material_extruders() const;
    std::set<uint16_t>  extruders() const;
    double              max_allowed_layer_height() const;
    bool                has_support_material() const;
    // Make sure the background processing has no access to this model_object during this call!
    void                auto_assign_extruders(ModelObject* model_object) const;

    const PrintConfig&          config() const { return m_config; }
    const PrintObjectConfig&    default_object_config() const { return m_default_object_config; }
    const PrintRegionConfig&    default_region_config() const { return m_default_region_config; }
    SpanOfConstPtrs<PrintObject> objects() const { return SpanOfConstPtrs<PrintObject>(const_cast<const PrintObject* const* const>(m_objects.data()), m_objects.size()); }
    PrintObject*                get_object(size_t idx) { return const_cast<PrintObject*>(m_objects[idx]); }
    const PrintObject*          get_object(size_t idx) const { return m_objects[idx]; }
    const PrintObject*          get_print_object_by_model_object_id(ObjectID object_id) const {
        auto it = std::find_if(m_objects.begin(), m_objects.end(),
                               [object_id](const PrintObject* obj) { return obj->model_object()->id() == object_id; });
        return (it == m_objects.end()) ? nullptr : *it;
    }
    // PrintObject by its ObjectID, to be used to uniquely bind slicing warnings to their source PrintObjects
    // in the notification center.
    const PrintObject*          get_object(ObjectID object_id) const { 
        auto it = std::find_if(m_objects.begin(), m_objects.end(), 
            [object_id](const PrintObject *obj) { return obj->id() == object_id; });
        return (it == m_objects.end()) ? nullptr : *it;
    }
    // How many of PrintObject::copies() over all print objects are there?
    // If zero, then the print is empty and the print shall not be executed.
    uint16_t                    num_object_instances() const;

    const std::optional<ExtrusionEntityCollection>& skirt_first_layer() const { return m_skirt_first_layer; }

    const ExtrusionEntityCollection& skirt() const { return m_skirt; }
    const ExtrusionEntityCollection& brim() const { return m_brim; }
    // Convex hull of the 1st layer extrusions, for bed leveling and placing the initial purge line.
    // It encompasses the object extrusions, support extrusions, skirt, brim, wipe tower.
    // It does NOT encompass user extrusions generated by custom G-code,
    // therefore it does NOT encompass the initial purge line.
    // It does NOT encompass MMU/MMU2 starting (wipe) areas.
    const Polygon&                   first_layer_convex_hull() const { return m_first_layer_convex_hull; }

    const PrintStatistics&      print_statistics() const { return m_print_statistics; }
    PrintStatistics&            print_statistics() { return m_print_statistics; }
    std::time_t                 timestamp_last_change() const { return m_timestamp_last_change; }

    // Wipe tower support.
    bool                        has_wipe_tower() const;
    const WipeTowerData&        wipe_tower_data(size_t extruders_cnt, double nozzle_diameter) const;
    const WipeTowerData&        wipe_tower_data() const { return wipe_tower_data(0,0); }
    const ToolOrdering& 		tool_ordering() const { return m_tool_ordering; }

	std::string                 output_filename(const std::string &filename_base = std::string()) const override;

    size_t                      num_print_regions() const throw() { return m_print_regions.size(); }
    const PrintRegion&          get_print_region(size_t idx) const  { return *m_print_regions[idx]; }
    const ToolOrdering&         get_tool_ordering() const { return m_wipe_tower_data.tool_ordering; }

    const Polygons& get_sequential_print_clearance_contours() const { return m_sequential_print_clearance_contours; }
//TODO: decide to use this one or the printconfig one.
    static bool sequential_print_horizontal_clearance_valid(const Print& print, Polygons* polygons = nullptr);

    //put this in public to be accessible for tests, it was in private before.
    bool                invalidate_state_by_config_options(const ConfigOptionResolver& new_config, const std::vector<t_config_option_key> &opt_keys);

    // Invalidates the step, and its depending steps in Print.
    //in public to invalidate gcode when the physical printer change. It's needed if we allow the gcode macro to read these values.
    bool                invalidate_step(PrintStep step);

    // just a little wrapper to let the user know that this print can only be modified to emit warnings & update advancement status, change stats.
    // TODO: have the status out of the printbase class and into another one, so we can have a const print & a mutable statusmonitor
    class StatusMonitor
    {
    private:
        Print& print;

    public:
        StatusMonitor(Print &print_mutable) : print(print_mutable) {}

        // need this extra method because active_step_add_warning is protected and so need the friend status, and Gcode has it.
        void active_step_add_warning(PrintStateBase::WarningLevel warning_level, const std::string &message, int message_id = 0)
        {
            print.active_step_add_warning(warning_level, message, message_id);
        }
        PrintStatistics &stats() { return print.m_print_statistics; }
        bool             set_started(PrintStep step) { return print.set_started(step); }
        PrintStateBase::TimeStamp set_done(PrintStep step) { return print.set_done(step); }
        
    };

protected:
private:

    void                _make_skirt_brim();
    void                _make_skirt(const PrintObjectPtrs &objects, ExtrusionEntityCollection &out, std::optional<ExtrusionEntityCollection> &out_first_layer);
    void                _make_wipe_tower();
    void                finalize_first_layer_convex_hull();
    void                alert_when_supports_needed();

    // Islands of objects and their supports extruded at the 1st layer.
    Polygons            first_layer_islands() const;
    // Return 4 wipe tower corners in the world coordinates (shifted and rotated), including the wipe tower brim.
    Points              first_layer_wipe_tower_corners() const;

    // Returns true if any of the print_objects has print_object_step valid.
    // That means data shared by all print objects of the print_objects span may still use the shared data.
    // Otherwise the shared data shall be released.
    // Unguarded variant, thus it shall only be called from main thread with background processing stopped.
    static bool         is_shared_print_object_step_valid_unguarded(SpanOfConstPtrs<PrintObject> print_objects, PrintObjectStep print_object_step);

    PrintConfig                             m_config;
    PrintObjectConfig                       m_default_object_config;
    PrintRegionConfig                       m_default_region_config;
    PrintObjectPtrs                         m_objects;
    PrintRegionPtrs                         m_print_regions;

    // Ordered collections of extrusion paths to build skirt loops and brim.
    std::optional<ExtrusionEntityCollection> m_skirt_first_layer;
    ExtrusionEntityCollection               m_skirt;
    ExtrusionEntityCollection               m_brim;
    // Convex hull of the 1st layer extrusions.
    // It encompasses the object extrusions, support extrusions, skirt, brim, wipe tower.
    // It does NOT encompass user extrusions generated by custom G-code,
    // therefore it does NOT encompass the initial purge line.
    // It does NOT encompass MMU/MMU2 starting (wipe) areas.
    Polygon                                 m_first_layer_convex_hull;
    Points                                  m_skirt_convex_hull;

    // Following section will be consumed by the GCodeGenerator.
    ToolOrdering 							m_tool_ordering;
    WipeTowerData                           m_wipe_tower_data {m_tool_ordering};

    // Estimated print time, filament consumed.
    PrintStatistics                         m_print_statistics;
    // time of last change, used by the gui to see if it needs to be updated
    std::time_t                             m_timestamp_last_change;

    // Cache to store sequential print clearance contours
    Polygons m_sequential_print_clearance_contours;

    // To allow GCode to set the Print's GCodeExport step status.
    //friend class GCodeGenerator;
    // To allow GCodeProcessor to emit warnings.
    //friend class GCodeProcessor;
    // Allow PrintObject to access m_mutex and m_cancel_callback.
    friend class PrintObject;

    std::optional<ConflictResult> m_conflict_result;
};

//for testing purpose (in printobject)
ExPolygons dense_fill_fit_to_size(const ExPolygon &polygon_to_cover,
    const ExPolygon& growing_area, const coord_t offset, float coverage);

} /* slic3r_Print_hpp_ */

#endif
