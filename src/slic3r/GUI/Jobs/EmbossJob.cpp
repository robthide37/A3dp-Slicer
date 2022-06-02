#include "EmbossJob.hpp"

#include <libslic3r/Model.hpp>
#include <libslic3r/Format/OBJ.hpp> // load_obj for default mesh
#include <libslic3r/CutSurface.hpp> // use surface cuts

#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/NotificationManager.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/MainFrame.hpp"
#include "slic3r/GUI/GUI.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"
#include "slic3r/GUI/Gizmos/GLGizmoEmboss.hpp"
#include "slic3r/GUI/CameraUtils.hpp"
#include "slic3r/GUI/format.hpp"
#include "slic3r/Utils/UndoRedo.hpp"

using namespace Slic3r;
using namespace GUI;

// private namespace
namespace priv{
/// <summary>
/// Assert check of inputs data
/// </summary>
/// <param name="input"></param>
/// <returns></returns>
bool check(const EmbossDataBase &input, bool check_fontfile = true);
bool check(const EmbossDataCreateVolume &input, bool is_main_thread = false);
bool check(const EmbossDataCreateObject &input);
bool check(const EmbossDataUpdate &input, bool is_main_thread = false);
bool check(const UseSurfaceData &input, bool is_main_thread = false);

// <summary>
/// Try to create mesh from text
/// </summary>
/// <param name="input">Text to convert on mesh
/// + Shape of characters + Property of font</param>
/// <param name="font">Font file with cache
/// NOTE: Cache glyphs is changed</param>
/// <param name="was_canceled">To check if process was canceled</param>
/// <returns>Triangle mesh model</returns>
template<typename Fnc>
static TriangleMesh try_create_mesh(const EmbossDataBase      &input,
                                    Emboss::FontFileWithCache &font,
                                    Fnc                        was_canceled);
template<typename Fnc>
static TriangleMesh create_mesh(EmbossDataBase &input,
                                Fnc             was_canceled,
                                Job::Ctl       &ctl);

/// <summary>
/// Create default mesh for embossed text
/// </summary>
/// <returns>Not empty model(index trinagle set - its)</returns>
static TriangleMesh create_default_mesh();

/// <summary>
/// Must be called on main thread
/// </summary>
/// <param name="mesh">New mesh data</param>
/// <param name="data">Text configuration, ...</param>
static void update_volume(TriangleMesh &&mesh, const EmbossDataUpdate &data);

/// <summary>
/// Select Volume from objects
/// </summary>
/// <param name="objects">All objects in scene</param>
/// <param name="volume_id">Identifier of volume in object</param>
/// <returns>Pointer to volume when exist otherwise nullptr</returns>
static ModelVolume *get_volume(ModelObjectPtrs &objects,
                               const ObjectID  &volume_id);
/// <summary>
/// extract scale in 2d
/// </summary>
/// <param name="fp">Property of font style</param>
/// <param name="ff">Font file for size --> unit per em</param>
/// <returns>scaling factor</returns>
static double get_shape_scale(const FontProp &fp, const Emboss::FontFile &ff);

/// <summary>
/// Create projection for cut surface from mesh
/// </summary>
/// <param name="tr">Volume transformation in object</param>
/// <param name="shape_scale">Convert shape to milimeters</param>
/// <param name="shape_bb">Bounding box 2d of shape to center result</param>
/// <param name="z_range">Bounding box 3d of model volume for projection ranges</param> 
/// <returns>Orthogonal cut_projection</returns>
static Emboss::OrthoProject create_projection_for_cut(
    Transform3d                    tr,
    double                         shape_scale,
    const BoundingBox             &shape_bb,
    const std::pair<float, float> &z_range);

/// <summary>
/// Create tranformation for emboss Cutted surface
/// </summary>
/// <param name="is_outside">True .. raise, False .. engrave</param>
/// <param name="emboss">Depth of embossing</param>
/// <param name="tr">Text voliume transformation inside object</param>
/// <param name="cut">Cutted surface from model</param>
/// <returns>Projection</returns>
static Emboss::OrthoProject3f create_emboss_projection(
    bool is_outside, float emboss, Transform3d tr, SurfaceCut &cut);

static void create_message(const std::string &message); // only in finalize
static bool process(std::exception_ptr &eptr);

class EmbossJobException: public std::exception {
public: EmbossJobException(char const *const message)
        : std::exception(message)
    {}
};
}

/////////////////
/// Create Volume
EmbossCreateVolumeJob::EmbossCreateVolumeJob(EmbossDataCreateVolume &&input)
    : m_input(std::move(input))
{
    assert(priv::check(m_input, true));
}

void EmbossCreateVolumeJob::process(Ctl &ctl) {
    if (!priv::check(m_input)) throw std::runtime_error("Bad input data for EmbossCreateVolumeJob.");
    auto was_canceled = [&ctl]()->bool { return ctl.was_canceled(); };
    m_result = priv::create_mesh(m_input, was_canceled, ctl);
    if (was_canceled()) return;
            
    // Create new volume inside of object
    const FontProp &font_prop = m_input.text_configuration.font_item.prop;
    Transform3d surface_trmat = Emboss::create_transformation_onto_surface(
            m_input.hit.position, m_input.hit.normal);
    Emboss::apply_transformation(font_prop, surface_trmat);
    m_transformation = m_input.hit_instance_tr.inverse() *
                       m_input.hit_object_tr * surface_trmat;    
}

void EmbossCreateVolumeJob::finalize(bool canceled, std::exception_ptr &eptr) {
    // doesn't care about exception when process was canceled by user
    if (canceled) {
        eptr = nullptr;
        return;
    }
    if (priv::process(eptr)) return;
    if (m_result.its.empty()) 
        return priv::create_message(_u8L("Can't create empty volume."));

    GUI_App         &app      = wxGetApp();
    Plater          *plater   = app.plater();
    ObjectList      *obj_list = app.obj_list();
    GLCanvas3D      *canvas   = plater->canvas3D();
    ModelObjectPtrs &objects  = plater->model().objects;

     // create volume in object
    size_t object_idx = m_input.object_idx;

    // Parent object for text volume was propably removed.
    // Assumption: User know what he does, so text volume is no more needed.
    if (objects.size() <= object_idx) 
        return;

    Plater::TakeSnapshot snapshot(plater, _L("Add Emboss text Volume"));

    ModelObject    *obj    = objects[object_idx];
    ModelVolumeType type   = m_input.volume_type;
    
    ModelVolume    *volume = obj->add_volume(std::move(m_result), type);

    // set a default extruder value, since user can't add it manually
    volume->config.set_key_value("extruder", new ConfigOptionInt(0));

    // do not allow model reload from disk
    volume->source.is_from_builtin_objects = true;

    volume->name               = m_input.volume_name;
    volume->text_configuration = std::move(m_input.text_configuration);
    volume->set_transformation(m_transformation);

    // update volume name in object list
    // updata selection after new volume added
    // change name of volume in right panel
    // select only actual volume
    // when new volume is created change selection to this volume
    auto add_to_selection = [volume](const ModelVolume *vol) {
        return vol == volume;
    };
    wxDataViewItemArray sel = obj_list->reorder_volumes_and_get_selection(
        m_input.object_idx, add_to_selection);
    if (!sel.IsEmpty()) obj_list->select_item(sel.front());

    // update printable state on canvas
    if (type == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object(object_idx);

    obj_list->selection_changed();

    // WHY selection_changed set manipulation to world ???
    // so I set it back to local --> RotationGizmo need it
    ObjectManipulation *manipul = wxGetApp().obj_manipul();
    manipul->set_coordinates_type(ECoordinatesType::Local);

    // redraw scene
    canvas->reload_scene(true);
}


/////////////////
/// Create Object
EmbossCreateObjectJob::EmbossCreateObjectJob(EmbossDataCreateObject &&input)
    : m_input(std::move(input))
{
    assert(priv::check(m_input));
}

void EmbossCreateObjectJob::process(Ctl &ctl) 
{
    if (!priv::check(m_input))
        throw std::runtime_error("Bad input data for EmbossCreateObjectJob.");

    auto was_canceled = [&ctl]()->bool { return ctl.was_canceled(); };
    m_result = priv::create_mesh(m_input, was_canceled, ctl);
    if (was_canceled()) return;

    // Create new object
    // calculate X,Y offset position for lay on platter in place of
    // mouse click
    Vec2d bed_coor = CameraUtils::get_z0_position(
        m_input.camera, m_input.screen_coor);

    // check point is on build plate:
    Points bed_shape_;
    bed_shape_.reserve(m_input.bed_shape.size());
    for (const Vec2d &p : m_input.bed_shape)
        bed_shape_.emplace_back(p.cast<int>());
    Polygon bed(bed_shape_);
    if (!bed.contains(bed_coor.cast<int>()))
        // mouse pose is out of build plate so create object in center of plate
        bed_coor = bed.centroid().cast<double>();

    double z = m_input.text_configuration.font_item.prop.emboss / 2;
    Vec3d  offset(bed_coor.x(), bed_coor.y(), z);
    offset -= m_result.center();
    Transform3d::TranslationType tt(offset.x(), offset.y(), offset.z());
    m_transformation = Transform3d(tt);
}

void EmbossCreateObjectJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    // doesn't care about exception when process was canceled by user
    if (canceled) {
        eptr = nullptr;
        return;
    }
    if (priv::process(eptr)) return;

    // only for sure
    if (m_result.empty()) 
        return priv::create_message(_u8L("Can't create empty object."));

    GUI_App    &app      = wxGetApp();
    Plater     *plater   = app.plater();
    ObjectList *obj_list = app.obj_list();
    GLCanvas3D *canvas   = plater->canvas3D();

    Plater::TakeSnapshot snapshot(plater, _L("Add Emboss text object"));

    // Create new object and change selection
    bool center = false;
    obj_list->load_mesh_object(std::move(m_result), m_input.volume_name,
                                center, &m_input.text_configuration,
                                &m_transformation);

    // When add new object selection is empty.
    // When cursor move and no one object is selected than
    // Manager::reset_all() So Gizmo could be closed before end of creation object
    GLGizmosManager &manager = canvas->get_gizmos_manager();
    if (manager.get_current_type() != GLGizmosManager::Emboss)
        manager.open_gizmo(GLGizmosManager::Emboss);   

    // redraw scene
    canvas->reload_scene(true);
}

/////////////////
/// Update Volume
EmbossUpdateJob::EmbossUpdateJob(EmbossDataUpdate&& input)
    : m_input(std::move(input))
{
    assert(priv::check(m_input, true));
}

void EmbossUpdateJob::process(Ctl &ctl)
{
    if (!priv::check(m_input))
        throw std::runtime_error("Bad input data for EmbossUpdateJob.");

    auto was_canceled = [&ctl, &cancel = m_input.cancel]()->bool {
        if (cancel->load()) return true;
        return ctl.was_canceled();
    };
    m_result = priv::try_create_mesh(m_input, m_input.font_file, was_canceled);
    if (was_canceled()) return;
    if (m_result.its.empty())
        throw priv::EmbossJobException(
            _u8L("Created text volume is empty. Change text or "
                 "font.").c_str());

    // center triangle mesh
    Vec3d shift = m_result.bounding_box().center();
    m_result.translate(-shift.cast<float>());    
}

void EmbossUpdateJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    // doesn't care about exception when process was canceled by user
    if (canceled || m_input.cancel->load()) {
        eptr = nullptr;
        return;
    }
    if (priv::process(eptr)) return;
    priv::update_volume(std::move(m_result), m_input);    
}

/////////////////
/// Cut Surface
UseSurfaceJob::UseSurfaceJob(UseSurfaceData &&input)
    : m_input(std::move(input))
{
    assert(priv::check(m_input, true));
}

void UseSurfaceJob::process(Ctl &ctl) {
    if (!priv::check(m_input)) 
        throw std::runtime_error("Bad input data for UseSurfaceJob.");
    
    // check cancelation of process
    auto was_canceled = [&ctl, &cancel = m_input.cancel]()->bool {
        if (cancel->load()) return true;
        return ctl.was_canceled();
    };

    const TextConfiguration &tc   = m_input.text_configuration;
    const char              *text = tc.text.c_str();
    const FontProp          &fp   = tc.font_item.prop;
    ExPolygons shapes = Emboss::text2shapes(m_input.font_file, text, fp);
    if (shapes.empty() || shapes.front().contour.empty())
        throw priv::EmbossJobException(
            _u8L("Font doesn't have any shape for given text.").c_str());
    
    if (was_canceled()) return;

    BoundingBox bb = get_extents(shapes);
    // TODO: merge input sources somehow
    const UseSurfaceData::ModelSource &source = m_input.sources[0];

    Transform3d   mesh_tr_inv       = source.tr.inverse();
    Transform3d   cut_projection_tr = mesh_tr_inv * m_input.text_tr;
    Transform3d   emboss_tr         = cut_projection_tr.inverse();
    BoundingBoxf3 mesh_bb_tr        = source.bb.transformed(emboss_tr);
    std::pair<float, float> z_range{mesh_bb_tr.min.z(), mesh_bb_tr.max.z()};

    const Emboss::FontFile &ff = *m_input.font_file.font_file;
    double shape_scale =  priv::get_shape_scale(fp, ff);
    Emboss::OrthoProject cut_projection = priv::create_projection_for_cut(
        cut_projection_tr, shape_scale, bb, z_range);

    // Use CGAL to cut surface from triangle mesh
    SurfaceCut cut = cut_surface(source.its, shapes, cut_projection);
    if (cut.empty())
        throw priv::EmbossJobException(
            _u8L("There is no valid surface for text projection.").c_str());
    if (was_canceled()) return;

    // !! Projection needs to transform cut
    Emboss::OrthoProject3f projection = priv::create_emboss_projection(
        m_input.is_outside, fp.emboss, emboss_tr, cut);
    
    indexed_triangle_set new_its = cut2model(cut, projection);
    assert(!new_its.empty());

    if (was_canceled()) return;
    //its_write_obj(new_its, "C:/data/temp/projected.obj"); // only debug
    m_result = TriangleMesh(std::move(new_its));
}

void UseSurfaceJob::finalize(bool canceled, std::exception_ptr &eptr)
{
    // doesn't care about exception when process was canceled by user
    if (canceled || m_input.cancel->load()) { 
        eptr = nullptr;
        return;
    }
    if (priv::process(eptr)) return;
    priv::update_volume(std::move(m_result), m_input);
}

////////////////////////////
/// private namespace implementation
bool priv::check(const EmbossDataBase &input, bool check_fontfile){
    bool res = true;
    if (check_fontfile) {
        assert(input.font_file.has_value());
        res &= input.font_file.has_value();
    }
    assert(!input.text_configuration.fix_3mf_tr.has_value());
    res &= !input.text_configuration.fix_3mf_tr.has_value();
    assert(!input.text_configuration.text.empty());
    res &= !input.text_configuration.text.empty();
    assert(!input.volume_name.empty());
    res &= !input.volume_name.empty();
    return res; 
}
bool priv::check(const EmbossDataCreateVolume &input, bool is_main_thread) {
    bool check_fontfile = false;
    bool res = check((EmbossDataBase) input, check_fontfile);
    assert(input.volume_type != ModelVolumeType::INVALID);
    res &= input.volume_type != ModelVolumeType::INVALID;
    assert(input.object_idx >= 0);
    res &= input.object_idx >= 0;
    if (is_main_thread)
        assert(input.object_idx < wxGetApp().model().objects.size());
    assert(input.screen_coor.x() >= 0.);
    res &= input.screen_coor.x() >= 0.;
    assert(input.screen_coor.y() >= 0.);
    res &= input.screen_coor.y() >= 0.;
    return res; 
}
bool priv::check(const EmbossDataCreateObject &input) {
    bool check_fontfile = false;
    bool res = check((EmbossDataBase) input, check_fontfile);
    assert(input.screen_coor.x() >= 0.);
    res &= input.screen_coor.x() >= 0.;
    assert(input.screen_coor.y() >= 0.);
    res &= input.screen_coor.y() >= 0.;
    assert(input.bed_shape.size() >= 3); // at least triangle
    res &= input.bed_shape.size() >= 3;
    return res;
}
bool priv::check(const EmbossDataUpdate &input, bool is_main_thread){
    bool res = check((EmbossDataBase) input);
    assert(input.volume_id.id >= 0);
    res &= input.volume_id.id >= 0;
    if (is_main_thread)
        assert(get_volume(wxGetApp().model().objects, input.volume_id) != nullptr);
    assert(input.cancel != nullptr);
    res &= input.cancel != nullptr;
    if (is_main_thread)
        assert(!input.cancel->load());
    return res;
}
bool priv::check(const UseSurfaceData &input, bool is_main_thread){
    bool res = check((EmbossDataUpdate) input, is_main_thread);
    assert(!input.sources.empty());
    res &= !input.sources.empty();
    return res;
}

template<typename Fnc>
TriangleMesh priv::try_create_mesh(const EmbossDataBase &input, Emboss::FontFileWithCache &font, Fnc was_canceled)
{
    const TextConfiguration &tc = input.text_configuration;
    const char *text = tc.text.c_str();
    const FontProp &prop = tc.font_item.prop;

    assert(font.has_value());
    if (!font.has_value()) return {};

    ExPolygons shapes = Emboss::text2shapes(font, text, prop);
    if (shapes.empty()) return {};
    if (was_canceled()) return {};

    const auto  &cn = prop.collection_number;
    unsigned int font_index = (cn.has_value()) ? *cn : 0;
    assert(font_index < font.font_file->infos.size());
    int unit_per_em = font.font_file->infos[font_index].unit_per_em;
    float scale    = prop.size_in_mm / unit_per_em;
    float depth    = prop.emboss / scale;
    auto  projectZ = std::make_unique<Emboss::ProjectZ>(depth);
    Emboss::ProjectScale project(std::move(projectZ), scale);
    if (was_canceled()) return {};
    return TriangleMesh(Emboss::polygons2model(shapes, project));
}

template<typename Fnc>
TriangleMesh priv::create_mesh(EmbossDataBase &input, Fnc was_canceled, Job::Ctl& ctl)
{
    // It is neccessary to create some shape
    // Emboss text window is opened by creation new emboss text object
    TriangleMesh result;
    if (input.font_file.has_value()) {
        result = try_create_mesh(input, input.font_file, was_canceled);
        if (was_canceled()) return {};
    }

    if (result.its.empty()) {
        result = priv::create_default_mesh();
        if (was_canceled()) return {};
        // only info
        ctl.call_on_main_thread([]() {
            create_message(_u8L("It is used default volume for embossed "
                                "text, try to change text or font for fix it."));
        });
    }

    assert(!result.its.empty());
    return result;
}

TriangleMesh priv::create_default_mesh()
{
    // When cant load any font use default object loaded from file
    std::string  path = Slic3r::resources_dir() + "/data/embossed_text.stl";
    TriangleMesh triangle_mesh;
    if (!load_obj(path.c_str(), &triangle_mesh)) {
        // when can't load mesh use cube
        return TriangleMesh(its_make_cube(36., 4., 2.5));
    }
    return triangle_mesh;
}

void priv::update_volume(TriangleMesh &&mesh,
                         const EmbossDataUpdate &data)
{
    // for sure that some object will be created
    if (mesh.its.empty())
        return priv::create_message("Empty mesh can't be created.");

    GUI_App &        app      = wxGetApp(); // may be move to input
    Plater *         plater   = app.plater();
    GLCanvas3D *     canvas   = plater->canvas3D();

    // Check emboss gizmo is still open
    GLGizmosManager &manager  = canvas->get_gizmos_manager();
    if (manager.get_current_type() != GLGizmosManager::Emboss) return;

    std::string snap_name = GUI::format(_L("Text: %1%"), data.text_configuration.text);
    Plater::TakeSnapshot snapshot(plater, snap_name, UndoRedo::SnapshotType::GizmoAction);
    ModelVolume *volume = get_volume(plater->model().objects, data.volume_id);
    // could appear when user delete edited volume
    if (volume == nullptr)
        return;

    // apply fix matrix made by store to .3mf
    const auto &tc = volume->text_configuration;
    assert(tc.has_value());
    if (tc.has_value() && tc->fix_3mf_tr.has_value())
        volume->set_transformation(volume->get_matrix() * tc->fix_3mf_tr->inverse());

    // update volume
    volume->set_mesh(std::move(mesh));
    volume->set_new_unique_id();
    volume->calculate_convex_hull();
    volume->get_object()->invalidate_bounding_box();
    volume->name               = data.volume_name;
    volume->text_configuration = data.text_configuration;

    // update volume in right panel( volume / object name)
    const Selection &selection = canvas->get_selection();
    const GLVolume * gl_volume = selection.get_volume(
        *selection.get_volume_idxs().begin());
    int object_idx = gl_volume->object_idx();
    int volume_idx = gl_volume->volume_idx();
    ObjectList *obj_list = app.obj_list();
    obj_list->update_name_in_list(object_idx, volume_idx);

    // update printable state on canvas
    if (volume->type() == ModelVolumeType::MODEL_PART)
        canvas->update_instance_printable_state_for_object((size_t) object_idx);
    
    // Move object on bed
    if (GLGizmoEmboss::is_text_object(volume)) {
        // Must be called after because GLVolume is invalidated by new_unique_id
        plater->CallAfter([object_idx]() {
            wxGetApp().plater()->canvas3D()->ensure_on_bed(object_idx, false);            
        });
    }

    // redraw scene
    bool refresh_immediately = false;
    canvas->reload_scene(refresh_immediately);
}

ModelVolume *priv::get_volume(ModelObjectPtrs &objects,
                               const ObjectID  &volume_id)
{
    for (ModelObject *obj : objects)
        for (ModelVolume *vol : obj->volumes)
            if (vol->id() == volume_id) return vol;
    return nullptr;
};


double priv::get_shape_scale(const FontProp &fp, const Emboss::FontFile &ff)
{
    const auto  &cn          = fp.collection_number;
    unsigned int font_index  = (cn.has_value()) ? *cn : 0;
    int          unit_per_em = ff.infos[font_index].unit_per_em;
    double       scale       = fp.size_in_mm / unit_per_em;
    // Shape is scaled for store point coordinate as integer
    return scale * Emboss::SHAPE_SCALE;
}

Emboss::OrthoProject priv::create_projection_for_cut(
    Transform3d                    tr,
    double                         shape_scale,
    const BoundingBox             &shape_bb,
    const std::pair<float, float> &z_range)
{
    // create sure that emboss object is bigger than source object
    const float safe_extension = 1.0f;
    float min_z = z_range.first - safe_extension;
    float max_z = z_range.second + safe_extension;
    assert(min_z < max_z);
    // range between min and max value
    double projection_size = max_z - min_z; 
    Matrix3d transformation_for_vector = tr.linear();
    // Projection must be negative value.
    // System of text coordinate
    // X .. from left to right
    // Y .. from bottom to top
    // Z .. from text to eye
    Vec3d untransformed_direction(0., 0., projection_size);
    Vec3f project_direction =
        (transformation_for_vector * untransformed_direction).cast<float>();

    // Projection is in direction from far plane
    tr.translate(Vec3d(0., 0., min_z));

    tr.scale(shape_scale);
    // Text alignemnt to center 2D
    Vec2d move = -(shape_bb.max + shape_bb.min).cast<double>() / 2.;
    //Vec2d move = -shape_bb.center().cast<double>(); // not precisse
    tr.translate(Vec3d(move.x(), move.y(), 0.));
    return Emboss::OrthoProject(tr, project_direction);
}

Emboss::OrthoProject3f priv::create_emboss_projection(
    bool is_outside, float emboss, Transform3d tr, SurfaceCut &cut)
{
    // Offset of clossed side to model
    const float surface_offset = 1e-3f; // [in mm]
    float 
        front_move = (is_outside) ? emboss : surface_offset,
        back_move  = -((is_outside) ? surface_offset : emboss);    
    its_transform(cut, tr.pretranslate(Vec3d(0., 0., front_move)));    
    Vec3f from_front_to_back(0.f, 0.f, back_move - front_move);
    return Emboss::OrthoProject3f(from_front_to_back);
}

bool priv::process(std::exception_ptr &eptr) { 
    if (!eptr) return false;
    try {
        std::rethrow_exception(eptr);
    } catch (priv::EmbossJobException &e) {
        create_message(e.what());
        eptr = nullptr;
    }
    return true;
}

#include <wx/msgdlg.h>

void priv::create_message(const std::string &message) {
    wxMessageBox(wxString(message), _L("Issue during embossing the text."),
                 wxOK | wxICON_WARNING);
}
