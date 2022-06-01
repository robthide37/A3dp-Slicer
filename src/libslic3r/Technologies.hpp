#ifndef _prusaslicer_technologies_h_
#define _prusaslicer_technologies_h_

//=============
// debug techs
//=============
// Shows camera target in the 3D scene
#define ENABLE_SHOW_CAMERA_TARGET 0
// Log debug messages to console when changing selection
#define ENABLE_SELECTION_DEBUG_OUTPUT 0
// Renders a small sphere in the center of the bounding box of the current selection when no gizmo is active
#define ENABLE_RENDER_SELECTION_CENTER 0
// Shows an imgui dialog with camera related data
#define ENABLE_CAMERA_STATISTICS 0
// Render the picking pass instead of the main scene (use [T] key to toggle between regular rendering and picking pass only rendering)
#define ENABLE_RENDER_PICKING_PASS 0
// Enable extracting thumbnails from selected gcode and save them as png files
#define ENABLE_THUMBNAIL_GENERATOR_DEBUG 0
// Disable synchronization of unselected instances
#define DISABLE_INSTANCES_SYNCH 0
// Use wxDataViewRender instead of wxDataViewCustomRenderer
#define ENABLE_NONCUSTOM_DATA_VIEW_RENDERING 0
// Enable G-Code viewer statistics imgui dialog
#define ENABLE_GCODE_VIEWER_STATISTICS 0
// Enable G-Code viewer comparison between toolpaths height and width detected from gcode and calculated at gcode generation 
#define ENABLE_GCODE_VIEWER_DATA_CHECKING 0
// Enable project dirty state manager debug window
#define ENABLE_PROJECT_DIRTY_STATE_DEBUG_WINDOW 0
// Disable using instanced models to render options in gcode preview
#define DISABLE_GCODEVIEWER_INSTANCED_MODELS 1


// Enable rendering of objects using environment map
#define ENABLE_ENVIRONMENT_MAP 0
// Enable smoothing of objects normals
#define ENABLE_SMOOTH_NORMALS 0


//====================
// 2.5.0.alpha1 techs
//====================
#define ENABLE_2_5_0_ALPHA1 1

// Enable changes in preview layout
#define ENABLE_PREVIEW_LAYOUT (1 && ENABLE_2_5_0_ALPHA1)
// Enable drawing the items in legend toolbar using icons
#define ENABLE_LEGEND_TOOLBAR_ICONS (1 && ENABLE_PREVIEW_LAYOUT)
// Enable coloring of toolpaths in preview by layer time
#define ENABLE_PREVIEW_LAYER_TIME (1 && ENABLE_2_5_0_ALPHA1)
// Enable showing time estimate for travel moves in legend
#define ENABLE_TRAVEL_TIME (1 && ENABLE_2_5_0_ALPHA1)
// Enable not killing focus in object manipulator fields when hovering over 3D scene
#define ENABLE_OBJECT_MANIPULATOR_FOCUS (0 && ENABLE_2_5_0_ALPHA1)
// Enable removal of wipe tower magic object_id equal to 1000
#define ENABLE_WIPETOWER_OBJECTID_1000_REMOVAL (1 && ENABLE_2_5_0_ALPHA1)
// Enable removal of legacy OpenGL calls
#define ENABLE_LEGACY_OPENGL_REMOVAL (1 && ENABLE_2_5_0_ALPHA1)
// Enable using vertex attributes and matrices in shaders
#define ENABLE_GL_SHADERS_ATTRIBUTES (1 && ENABLE_LEGACY_OPENGL_REMOVAL)
// Enable rendering imgui using shaders
#define ENABLE_GL_IMGUI_SHADERS (1 && ENABLE_GL_SHADERS_ATTRIBUTES)
// Shows an imgui dialog with GLModel statistics data
#define ENABLE_GLMODEL_STATISTICS (0 && ENABLE_LEGACY_OPENGL_REMOVAL)
// Enable show non-manifold edges
#define ENABLE_SHOW_NON_MANIFOLD_EDGES (1 && ENABLE_2_5_0_ALPHA1)
// Enable rework of Reload from disk command
#define ENABLE_RELOAD_FROM_DISK_REWORK (1 && ENABLE_2_5_0_ALPHA1)
// Enable showing toolpaths center of gravity
#define ENABLE_SHOW_TOOLPATHS_COG (1 && ENABLE_2_5_0_ALPHA1)
// Enable recalculating toolpaths when switching to/from volumetric rate visualization
#define ENABLE_VOLUMETRIC_RATE_TOOLPATHS_RECALC (1 && ENABLE_2_5_0_ALPHA1)
// Enable editing volumes transformation in world coordinates and instances in local coordinates
#define ENABLE_WORLD_COORDINATE (1 && ENABLE_2_5_0_ALPHA1)
// Enable modified camera control using mouse
#define ENABLE_NEW_CAMERA_MOVEMENTS (1 && ENABLE_2_5_0_ALPHA1)
// Enable modified rectangle selection
#define ENABLE_NEW_RECTANGLE_SELECTION (1 && ENABLE_2_5_0_ALPHA1)
// Enable alternative version of file_wildcards()
#define ENABLE_ALTERNATIVE_FILE_WILDCARDS_GENERATOR (1 && ENABLE_2_5_0_ALPHA1)
// Enable processing of gcode G2 and G3 lines
#define ENABLE_PROCESS_G2_G3_LINES (1 && ENABLE_2_5_0_ALPHA1)
// Enable fix of used filament data exported to gcode file
#define ENABLE_USED_FILAMENT_POST_PROCESS (1 && ENABLE_2_5_0_ALPHA1)
// Enable gizmo grabbers to share common models
#define ENABLE_GIZMO_GRABBER_REFACTOR (1 && ENABLE_2_5_0_ALPHA1)
// Enable copy of custom bed model and texture
#define ENABLE_COPY_CUSTOM_BED_MODEL_AND_TEXTURE (1 && ENABLE_2_5_0_ALPHA1)

#endif // _prusaslicer_technologies_h_
