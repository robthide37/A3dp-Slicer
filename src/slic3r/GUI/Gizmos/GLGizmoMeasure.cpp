// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/PresetBundle.hpp"

#include <imgui/imgui_internal.h>

#include <numeric>

#include <GL/glew.h>

#include <wx/clipbrd.h>

#if ENABLE_MEASURE_GIZMO

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA SELECTED_1ST_COLOR   = { 0.25f, 0.75f, 0.75f, 1.0f };
static const Slic3r::ColorRGBA SELECTED_2ND_COLOR   = { 0.75f, 0.25f, 0.75f, 1.0f };

static const int POINT_ID         = 100;
static const int EDGE_ID          = 200;
static const int CIRCLE_ID        = 300;
static const int PLANE_ID         = 400;
static const int SELECTION_1_ID   = 501;
static const int SELECTION_2_ID   = 502;

static const float TRIANGLE_BASE = 10.0f;
static const float TRIANGLE_HEIGHT = TRIANGLE_BASE * 1.618033f;

static const std::string CTRL_STR =
#ifdef __APPLE__
"⌘"
#else
"Ctrl"
#endif //__APPLE__
;

static std::string surface_feature_type_as_string(Measure::SurfaceFeatureType type)
{
    switch (type)
    {
    default:
    case Measure::SurfaceFeatureType::Undef:  { return _u8L("Undefined"); }
    case Measure::SurfaceFeatureType::Point:  { return _u8L("Vertex"); }
    case Measure::SurfaceFeatureType::Edge:   { return _u8L("Edge"); }
    case Measure::SurfaceFeatureType::Circle: { return _u8L("Circle"); }
    case Measure::SurfaceFeatureType::Plane:  { return _u8L("Plane"); }
    }
}

static std::string point_on_feature_type_as_string(Measure::SurfaceFeatureType type, int hover_id)
{
    std::string ret;
    switch (type) {
    case Measure::SurfaceFeatureType::Point:  { ret = _u8L("Vertex"); break; }
    case Measure::SurfaceFeatureType::Edge:   { ret = (hover_id == POINT_ID) ? _u8L("Center of edge") : _u8L("Point on edge"); break; }
    case Measure::SurfaceFeatureType::Circle: { ret = (hover_id == POINT_ID) ? _u8L("Center of circle") : _u8L("Point on circle"); break; }
    case Measure::SurfaceFeatureType::Plane:  { ret = _u8L("Point on plane"); break; }
    default: { assert(false); break; }
    }
    return ret;
}

static GLModel::Geometry init_plane_data(const indexed_triangle_set& its, const std::vector<std::vector<int>>& planes_triangles, int idx)
{
    assert(0 <= idx && idx < (int)planes_triangles.size());
    const std::vector<int>& triangle_indices = planes_triangles[idx];

    GLModel::Geometry init_data;
    init_data.format = { GUI::GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3 };
    unsigned int i = 0;
    for (int idx : triangle_indices) {
        const Vec3f& v0 = its.vertices[its.indices[idx][0]];
        const Vec3f& v1 = its.vertices[its.indices[idx][1]];
        const Vec3f& v2 = its.vertices[its.indices[idx][2]];

        const Vec3f n = (v1 - v0).cross(v2 - v0).normalized();
        init_data.add_vertex(v0, n);
        init_data.add_vertex(v1, n);
        init_data.add_vertex(v2, n);
        init_data.add_triangle(i, i + 1, i + 2);
        i += 3;
    }

    return init_data;
}

GLGizmoMeasure::GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{
    GLModel::Geometry sphere_geometry = smooth_sphere(16, 7.5f);
    m_sphere.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(sphere_geometry.get_as_indexed_triangle_set()));
    m_sphere.model.init_from(std::move(sphere_geometry));

    GLModel::Geometry cylinder_geometry = smooth_cylinder(16, 5.0f, 1.0f);
    m_cylinder.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(cylinder_geometry.get_as_indexed_triangle_set()));
    m_cylinder.model.init_from(std::move(cylinder_geometry));
}

bool GLGizmoMeasure::on_mouse(const wxMouseEvent &mouse_event)
{
    m_mouse_pos = { double(mouse_event.GetX()), double(mouse_event.GetY()) };

    if (mouse_event.Moving()) {
        // only for sure 
        m_mouse_left_down = false;
        return false;
    }
    else if (mouse_event.LeftDown()) {
        if (m_hover_id != -1) {
            SelectedFeatures selected_features_old = m_selected_features;
            m_mouse_left_down = true;
            
            auto item_from_feature = [this]() {
                SelectedFeatures::Item item;
                if (m_hover_id == SELECTION_1_ID && m_selected_features.first.feature.has_value())
                    item = m_selected_features.first;
                else if (m_hover_id == SELECTION_2_ID && m_selected_features.second.feature.has_value())
                    item = m_selected_features.second;
                else {
                    item = {
                        (m_mode == EMode::ExtendedSelection) ? point_on_feature_type_as_string(m_curr_feature->get_type(), m_hover_id) : surface_feature_type_as_string(m_curr_feature->get_type()),
                        (m_mode == EMode::ExtendedSelection) ? Measure::SurfaceFeature(*m_curr_point_on_feature_position) : m_curr_feature
                    };
                }
                return item;
            };

            if (m_selected_features.first.feature.has_value()) {
                auto it = std::find_if(m_selection_raycasters.begin(), m_selection_raycasters.end(),
                    [](std::shared_ptr<SceneRaycasterItem> item) { return SceneRaycaster::decode_id(SceneRaycaster::EType::Gizmo, item->get_id()) == SELECTION_2_ID; });
                if (it != m_selection_raycasters.end())
                    m_selection_raycasters.erase(it);
                m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, SELECTION_2_ID);

                const SelectedFeatures::Item item = item_from_feature();
                if (m_selected_features.first != item) {
                    if (m_selected_features.second == item)
                        m_selected_features.second.reset();
                    else {
                        m_selected_features.second = item;
                        if (m_mode == EMode::ExtendedSelection)
                            m_selection_raycasters.push_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, SELECTION_2_ID, *m_sphere.mesh_raycaster));
                    }
                }
            }
            else {
                const SelectedFeatures::Item item = item_from_feature();
                m_selected_features.first = item;
                if (m_mode == EMode::ExtendedSelection)
                    m_selection_raycasters.push_back(m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, SELECTION_1_ID, *m_sphere.mesh_raycaster));
            }

            if (m_selected_features != selected_features_old && m_selected_features.second.feature.has_value())
                m_measurement_result = Measure::get_measurement(*m_selected_features.first.feature, *m_selected_features.second.feature);

            return true;
        }

        // fix: prevent restart gizmo when reselect object
        // take responsibility for left up
        if (m_parent.get_first_hover_volume_idx() >= 0)
            m_mouse_left_down = true;
    }
    else if (mouse_event.LeftUp()) {
        if (m_mouse_left_down) {
            // responsible for mouse left up after selecting plane
            m_mouse_left_down = false;
            return true;
        }
    }
    else if (mouse_event.RightDown() && mouse_event.CmdDown()) {
        m_selected_features.reset();
        m_selection_raycasters.clear();
        m_imgui->set_requires_extra_frame();
    }
    else if (mouse_event.Leaving())
        m_mouse_left_down = false;

    return false;
}

void GLGizmoMeasure::data_changed()
{
    const Selection&  selection     = m_parent.get_selection();
    const ModelObject* model_object = nullptr;
    const ModelVolume* model_volume = nullptr;
    if (selection.is_single_full_instance() ||
        selection.is_from_single_object() ) {        
        model_object = selection.get_model()->objects[selection.get_object_idx()];
        model_volume = model_object->volumes[selection.get_first_volume()->volume_idx()];
    }
    if (model_object != m_old_model_object || model_volume != m_old_model_volume)
        update_if_needed();

    m_last_inv_zoom = 0.0f;
    m_last_plane_idx = -1;
    m_selected_features.reset();
    m_selection_raycasters.clear();
}

bool GLGizmoMeasure::gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down)
{
    if (action == SLAGizmoEventType::CtrlDown) {
        if (m_ctrl_kar_filter.is_first()) {
            if (m_curr_feature.has_value()) {
                m_mode = EMode::ExtendedSelection;
                disable_scene_raycasters();
            }
        }

        m_ctrl_kar_filter.increase_count();
    }
    else if (action == SLAGizmoEventType::CtrlUp) {
        m_ctrl_kar_filter.reset_count();
        m_mode = EMode::BasicSelection;
        restore_scene_raycasters_state();
    }

    return true;
}

bool GLGizmoMeasure::on_init()
{
    m_shortcut_key = WXK_CONTROL_U;
    return true;
}

void GLGizmoMeasure::on_set_state()
{
    if (m_state == Off) {
        m_ctrl_kar_filter.reset_count();
        m_curr_feature.reset();
        m_curr_point_on_feature_position.reset();
        restore_scene_raycasters_state();
    }
    else {
        m_mode = EMode::BasicSelection;
        // store current state of scene raycaster for later use
        m_scene_raycasters.clear();
        auto scene_raycasters = m_parent.get_raycasters_for_picking(SceneRaycaster::EType::Volume);
        if (scene_raycasters != nullptr) {
            m_scene_raycasters.reserve(scene_raycasters->size());
            for (auto r : *scene_raycasters) {
                SceneRaycasterState state = { r, r->is_active() };
                m_scene_raycasters.emplace_back(state);
            }
        }
    }
}

CommonGizmosDataID GLGizmoMeasure::on_get_requirements() const
{
    return CommonGizmosDataID(int(CommonGizmosDataID::SelectionInfo) | int(CommonGizmosDataID::Raycaster));
}

std::string GLGizmoMeasure::on_get_name() const
{
    return _u8L("Measure");
}

bool GLGizmoMeasure::on_is_activable() const
{
    const Selection& selection = m_parent.get_selection();
    bool res = (wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() == ptSLA) ?
        selection.is_single_full_instance() :
        selection.is_single_volume() || selection.is_single_volume_instance();
    if (res)
        res &= !selection.get_first_volume()->is_sinking();
    return res;
}

void GLGizmoMeasure::on_render()
{
#if ENABLE_MEASURE_GIZMO_DEBUG
    render_debug_dialog();
#endif // ENABLE_MEASURE_GIZMO_DEBUG

    // do not render if the user is panning/rotating the 3d scene
    if (m_parent.is_mouse_dragging())
        return;

    const Selection& selection = m_parent.get_selection();

    if ((wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() == ptSLA && selection.is_single_full_instance()) ||
        (selection.is_single_volume() || selection.is_single_volume_instance())) {
        update_if_needed();

        m_volume_matrix = selection.get_first_volume()->world_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const float inv_zoom = (float)camera.get_inv_zoom();

        Vec3f position_on_model;
        Vec3f normal_on_model;
        size_t model_facet_idx;
        const bool mouse_on_object = m_c->raycaster()->raycasters().front()->unproject_on_mesh(m_mouse_pos, m_volume_matrix, camera, position_on_model, normal_on_model, nullptr, &model_facet_idx);
        const bool is_hovering_on_locked_feature = m_mode == EMode::ExtendedSelection && m_hover_id != -1;

        if (m_mode == EMode::BasicSelection) {
            std::optional<Measure::SurfaceFeature> curr_feature = mouse_on_object ? m_measuring->get_feature(model_facet_idx, position_on_model.cast<double>()) : std::nullopt;
            m_curr_point_on_feature_position.reset();
            if (m_curr_feature != curr_feature) {
                m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, POINT_ID);
                m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, EDGE_ID);
                m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, PLANE_ID);
                m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo, CIRCLE_ID);
                m_raycasters.clear();
                m_curr_feature = curr_feature;
                if (!m_curr_feature.has_value())
                    return;

                switch (m_curr_feature->get_type()) {
                default: { assert(false); break; }
                case Measure::SurfaceFeatureType::Point:
                {
                    m_raycasters.insert({ POINT_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, POINT_ID, *m_sphere.mesh_raycaster) });
                    break;
                }
                case Measure::SurfaceFeatureType::Edge:
                {
                    m_raycasters.insert({ EDGE_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, EDGE_ID, *m_cylinder.mesh_raycaster) });
                    if (m_curr_feature->get_extra_point().has_value())
                        m_raycasters.insert({ POINT_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, POINT_ID, *m_sphere.mesh_raycaster) });
                    break;
                }
                case Measure::SurfaceFeatureType::Circle:
                {
                    const auto [center, radius, normal] = m_curr_feature->get_circle();
                    if (m_last_inv_zoom != inv_zoom || m_last_circle != m_curr_feature) {
                        m_last_inv_zoom = inv_zoom;
                        m_last_circle   = m_curr_feature;
                        m_circle.reset();
                        GLModel::Geometry circle_geometry = smooth_torus(64, 16, float(radius), 5.0f * inv_zoom);
                        m_circle.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(circle_geometry.get_as_indexed_triangle_set()));
                        m_circle.model.init_from(std::move(circle_geometry));
                    }

                    m_raycasters.insert({ CIRCLE_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, CIRCLE_ID, *m_circle.mesh_raycaster) });
                    m_raycasters.insert({ POINT_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, POINT_ID, *m_sphere.mesh_raycaster) });
                    break;
                }
                case Measure::SurfaceFeatureType::Plane:
                {
                    const auto [idx, normal, point] = m_curr_feature->get_plane();
                    if (m_last_plane_idx != idx) {
                        m_last_plane_idx = idx;
                        const indexed_triangle_set its = (m_old_model_volume != nullptr) ? m_old_model_volume->mesh().its : m_old_model_object->volumes.front()->mesh().its;
                        const std::vector<std::vector<int>> planes_triangles = m_measuring->get_planes_triangle_indices();
                        GLModel::Geometry init_data = init_plane_data(its, planes_triangles, idx);
                        m_plane.reset();
                        m_plane.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(init_data.get_as_indexed_triangle_set()));
                        m_plane.model.init_from(std::move(init_data));
                    }

                    m_raycasters.insert({ PLANE_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, PLANE_ID, *m_plane.mesh_raycaster) });
                    break;
                }
                }
            }
        }
        else if (is_hovering_on_locked_feature) {
            auto position_on_feature = [this](int feature_type_id, const Camera& camera, std::function<Vec3f(const Vec3f&)> callback = nullptr) -> Vec3d {
                auto it = m_raycasters.find(feature_type_id);
                if (it != m_raycasters.end() && it->second != nullptr) {
                    Vec3f p;
                    Vec3f n;
                    const Transform3d& trafo = it->second->get_transform();
                    bool res = it->second->get_raycaster()->closest_hit(m_mouse_pos, trafo, camera, p, n);
                    assert(res);
                    if (res) {
                        if (callback)
                            p = callback(p);
                        return trafo * p.cast<double>();
                    }
                }
                return Vec3d::Zero();
            };

            switch (m_curr_feature->get_type())
            {
            default: { assert(false); break; }
            case Measure::SurfaceFeatureType::Point:
            {
                m_curr_point_on_feature_position = m_curr_feature->get_point();
                break;
            }
            case Measure::SurfaceFeatureType::Edge:
            {
                const std::optional<Vec3d> extra = m_curr_feature->get_extra_point();
                if (extra.has_value() && m_hover_id == POINT_ID)
                    m_curr_point_on_feature_position = *extra;
                else
                    m_curr_point_on_feature_position = m_volume_matrix.inverse() * position_on_feature(EDGE_ID, camera, [](const Vec3f& v) { return Vec3f(0.0f, 0.0f, v.z()); });
                break;
            }
            case Measure::SurfaceFeatureType::Plane:
            {
                m_curr_point_on_feature_position = m_volume_matrix.inverse() * position_on_feature(PLANE_ID, camera);
                break;
            }
            case Measure::SurfaceFeatureType::Circle:
            {
                const auto [center, radius, normal] = m_curr_feature->get_circle();
                if (m_hover_id == POINT_ID)
                    m_curr_point_on_feature_position = center;
                else {
                    const float r = radius; // needed for the following lambda
                    m_curr_point_on_feature_position = m_volume_matrix.inverse() * position_on_feature(CIRCLE_ID, camera, [r](const Vec3f& v) {
                        float angle = std::atan2(v.y(), v.x());
                        if (angle < 0.0f)
                            angle += 2.0f * float(M_PI);
                        return Vec3f(float(r) * std::cos(angle), float(r) * std::sin(angle), 0.0f);
                        });
                }
                break;
            }
            }
        }

        if (!m_curr_feature.has_value() && !m_selected_features.first.feature.has_value())
            return;

        GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
        if (shader == nullptr)
            return;

        shader->start_using();
        shader->set_uniform("emission_factor", 0.25f);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());

        glsafe(::glClear(GL_DEPTH_BUFFER_BIT));
        glsafe(::glEnable(GL_DEPTH_TEST));

        const Transform3d& view_matrix = camera.get_view_matrix();

        auto set_matrix_uniforms = [shader, &view_matrix](const Transform3d& model_matrix) {
            const Transform3d view_model_matrix = view_matrix * model_matrix;
            shader->set_uniform("view_model_matrix", view_model_matrix);
            const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
            shader->set_uniform("view_normal_matrix", view_normal_matrix);
        };

        auto render_feature = [this, set_matrix_uniforms](const Measure::SurfaceFeature& feature, const std::vector<ColorRGBA>& colors,
            const Transform3d& model_matrix, float inv_zoom, bool update_raycasters) {
            switch (feature.get_type())
            {
            default: { assert(false); break; }
            case Measure::SurfaceFeatureType::Point:
            {
                const Vec3d& position = feature.get_point();
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(position) * Geometry::scale_transform(inv_zoom);
                set_matrix_uniforms(feature_matrix);
                m_sphere.model.set_color(colors.front());
                m_sphere.model.render();
                if (update_raycasters) {
                    auto it = m_raycasters.find(POINT_ID);
                    if (it != m_raycasters.end() && it->second != nullptr)
                        it->second->set_transform(feature_matrix);
                }
                break;
            }
            case Measure::SurfaceFeatureType::Circle:
            {
                const auto& [center, radius, normal] = feature.get_circle();
                // render center
                const Transform3d center_matrix = model_matrix * Geometry::translation_transform(center) * Geometry::scale_transform(inv_zoom);
                set_matrix_uniforms(center_matrix);
                m_sphere.model.set_color(colors.front());
                m_sphere.model.render();
                if (update_raycasters) {
                    auto it = m_raycasters.find(POINT_ID);
                    if (it != m_raycasters.end() && it->second != nullptr)
                        it->second->set_transform(center_matrix);
                }
                // render circle
                const Transform3d circle_matrix = model_matrix * Geometry::translation_transform(center) * Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), normal);
                set_matrix_uniforms(circle_matrix);
                m_circle.model.set_color(colors.back());
                m_circle.model.render();
                if (update_raycasters) {
                    auto it = m_raycasters.find(CIRCLE_ID);
                    if (it != m_raycasters.end() && it->second != nullptr)
                        it->second->set_transform(circle_matrix);
                }
                break;
            }
            case Measure::SurfaceFeatureType::Edge:
            {
                const auto& [start, end] = feature.get_edge();
                // render extra point
                const std::optional<Vec3d> extra = feature.get_extra_point();
                if (extra.has_value()) {
                    const Transform3d point_matrix = model_matrix * Geometry::translation_transform(*extra) * Geometry::scale_transform(inv_zoom);
                    set_matrix_uniforms(point_matrix);
                    m_sphere.model.set_color(colors.front());
                    m_sphere.model.render();
                    if (update_raycasters) {
                        auto it = m_raycasters.find(POINT_ID);
                        if (it != m_raycasters.end() && it->second != nullptr)
                            it->second->set_transform(point_matrix);
                    }
                }
                // render edge
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(start) *
                    Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), end - start) *
                    Geometry::scale_transform({ (double)inv_zoom, (double)inv_zoom, (end - start).norm() });
                set_matrix_uniforms(feature_matrix);
                m_cylinder.model.set_color(colors.back());
                m_cylinder.model.render();
                if (update_raycasters) {
                    auto it = m_raycasters.find(EDGE_ID);
                    if (it != m_raycasters.end() && it->second != nullptr)
                        it->second->set_transform(feature_matrix);
                }
                break;
            }
            case Measure::SurfaceFeatureType::Plane:
            {
                const auto& [idx, normal, pt] = feature.get_plane();
                assert(idx < m_plane_models_cache.size());
                set_matrix_uniforms(model_matrix);
                m_plane_models_cache[idx].set_color(colors.front());
                m_plane_models_cache[idx].render();
                if (update_raycasters) {
                    auto it = m_raycasters.find(PLANE_ID);
                    if (it != m_raycasters.end() && it->second != nullptr)
                        it->second->set_transform(model_matrix);
                }
                break;
            }
            }
        };

        auto hover_selection_color = [this]() {
            return saturate(!m_selected_features.first.feature.has_value() ? SELECTED_1ST_COLOR : SELECTED_2ND_COLOR, 1.5f);
        };

        auto hovering_color = [this, hover_selection_color, &selection]() {
            return (m_mode == EMode::ExtendedSelection) ? selection.get_first_volume()->render_color : hover_selection_color();
        };

        if (m_curr_feature.has_value()) {
            std::vector<ColorRGBA> colors;
            if (m_selected_features.first.feature.has_value() && *m_curr_feature == *m_selected_features.first.feature)
                colors.emplace_back(SELECTED_1ST_COLOR);
            else if (m_selected_features.second.feature.has_value() && *m_curr_feature == *m_selected_features.second.feature)
                colors.emplace_back(SELECTED_2ND_COLOR);
            else {
                switch (m_curr_feature->get_type())
                {
                default: { assert(false); break; }
                case Measure::SurfaceFeatureType::Point:
                {
                    colors.emplace_back(hover_selection_color());
                    break;
                }
                case Measure::SurfaceFeatureType::Edge:
                case Measure::SurfaceFeatureType::Circle:
                {
                    colors.emplace_back((m_hover_id == POINT_ID) ? hover_selection_color() : hovering_color());
                    colors.emplace_back(hovering_color());
                    break;
                }
                case Measure::SurfaceFeatureType::Plane:
                {
                    colors.emplace_back(hovering_color());
                    break;
                }
                }
            }

            render_feature(*m_curr_feature, colors, m_volume_matrix, inv_zoom, true);
        }

        if (m_selected_features.first.feature.has_value() && (!m_curr_feature.has_value() || *m_curr_feature != *m_selected_features.first.feature)) {
            std::vector<ColorRGBA> colors;
            colors.emplace_back(SELECTED_1ST_COLOR);
            render_feature(*m_selected_features.first.feature, colors, m_volume_matrix, inv_zoom, false);
            if (m_selected_features.first.feature->get_type() == Measure::SurfaceFeatureType::Point) {
                auto it = std::find_if(m_selection_raycasters.begin(), m_selection_raycasters.end(),
                    [](std::shared_ptr<SceneRaycasterItem> item) { return SceneRaycaster::decode_id(SceneRaycaster::EType::Gizmo, item->get_id()) == SELECTION_1_ID; });
                if (it != m_selection_raycasters.end())
                    (*it)->set_transform(m_volume_matrix * Geometry::translation_transform(m_selected_features.first.feature->get_point()) * Geometry::scale_transform(inv_zoom));
            }
        }
        if (m_selected_features.second.feature.has_value() && (!m_curr_feature.has_value() || *m_curr_feature != *m_selected_features.second.feature)) {
            std::vector<ColorRGBA> colors;
            colors.emplace_back(SELECTED_2ND_COLOR);
            render_feature(*m_selected_features.second.feature, colors, m_volume_matrix, inv_zoom, false);
            if (m_selected_features.second.feature->get_type() == Measure::SurfaceFeatureType::Point) {
                auto it = std::find_if(m_selection_raycasters.begin(), m_selection_raycasters.end(),
                    [](std::shared_ptr<SceneRaycasterItem> item) { return SceneRaycaster::decode_id(SceneRaycaster::EType::Gizmo, item->get_id()) == SELECTION_2_ID; });
                if (it != m_selection_raycasters.end())
                    (*it)->set_transform(m_volume_matrix * Geometry::translation_transform(m_selected_features.second.feature->get_point()) * Geometry::scale_transform(inv_zoom));
            }
        }

        if (is_hovering_on_locked_feature && m_curr_point_on_feature_position.has_value()) {
            if (m_hover_id != POINT_ID) {
                const Transform3d matrix = m_volume_matrix * Geometry::translation_transform(*m_curr_point_on_feature_position) * Geometry::scale_transform(inv_zoom);
                set_matrix_uniforms(matrix);
                m_sphere.model.set_color(hover_selection_color());
                m_sphere.model.render();
            }
        }

        shader->stop_using();
    }

    render_dimensioning();
}

void GLGizmoMeasure::update_if_needed()
{
    auto update_plane_models_cache = [this](const indexed_triangle_set& its) {
        m_plane_models_cache.clear();
        const std::vector<std::vector<int>> planes_triangles = m_measuring->get_planes_triangle_indices();
        for (int idx = 0; idx < (int)planes_triangles.size(); ++idx) {
            m_plane_models_cache.emplace_back(GLModel());
            GLModel::Geometry init_data = init_plane_data(its, planes_triangles, idx);
            m_plane_models_cache.back().init_from(std::move(init_data));
        }
    };

    auto do_update = [this, update_plane_models_cache](const ModelObject* object, const ModelVolume* volume) {
        const indexed_triangle_set& its = (volume != nullptr) ? volume->mesh().its : object->volumes.front()->mesh().its;
        m_measuring.reset(new Measure::Measuring(its));

        update_plane_models_cache(its);

        // Let's save what we calculated it from:
        m_volumes_matrices.clear();
        m_volumes_types.clear();
        m_first_instance_scale  = Vec3d::Ones();
        m_first_instance_mirror = Vec3d::Ones();
        if (object != nullptr) {
            for (const ModelVolume* vol : object->volumes) {
                m_volumes_matrices.push_back(vol->get_matrix());
                m_volumes_types.push_back(vol->type());
            }
            m_first_instance_scale  = object->instances.front()->get_scaling_factor();
            m_first_instance_mirror = object->instances.front()->get_mirror();
        }
        m_old_model_object = object;
        m_old_model_volume = volume;
    };

    const ModelObject* mo = m_c->selection_info()->model_object();
    const ModelVolume* mv = m_c->selection_info()->model_volume();
    if (m_state != On || (mo == nullptr && mv == nullptr))
        return;

    if (mo == nullptr)
        mo = mv->get_object();

    if (mo->instances.empty())
        return;

    if (!m_measuring || mo != m_old_model_object || mv != m_old_model_volume || mo->volumes.size() != m_volumes_matrices.size())
        do_update(mo, mv);

    // We want to recalculate when the scale changes - some planes could (dis)appear.
    if (!mo->instances.front()->get_scaling_factor().isApprox(m_first_instance_scale) ||
        !mo->instances.front()->get_mirror().isApprox(m_first_instance_mirror))
        do_update(mo, mv);

    for (unsigned int i = 0; i < mo->volumes.size(); ++i) {
        if (!mo->volumes[i]->get_matrix().isApprox(m_volumes_matrices[i]) ||
            mo->volumes[i]->type() != m_volumes_types[i]) {
            do_update(mo, mv);
            break;
        }
    }
}

void GLGizmoMeasure::disable_scene_raycasters()
{
    for (auto r : m_scene_raycasters) {
        r.raycaster->set_active(false);
    }
}

void GLGizmoMeasure::restore_scene_raycasters_state()
{
    for (auto r : m_scene_raycasters) {
        r.raycaster->set_active(r.state);
    }
}

class DimensioningHelper
{
    struct Cache
    {
        std::array<int, 4> viewport;
        Matrix4d ndc_to_ss_matrix;
        Transform3d ndc_to_ss_matrix_inverse;
    };

    static Cache s_cache;

public:
    static Vec3d model_to_world(const Vec3d& model, const Transform3d& world_matrix) {
        return world_matrix * model;
    }

    static Vec4d world_to_clip(const Vec3d& world, const Matrix4d& projection_view_matrix) {
        return projection_view_matrix * Vec4d(world.x(), world.y(), world.z(), 1.0);
    }

    static Vec3d clip_to_ndc(const Vec4d& clip) {
        return Vec3d(clip.x(), clip.y(), clip.z()) / clip.w();
    }

    static Vec2d ndc_to_ss(const Vec3d& ndc, const std::array<int, 4>& viewport) {
        const double half_w = 0.5 * double(viewport[2]);
        const double half_h = 0.5 * double(viewport[3]);
        return { half_w * ndc.x() + double(viewport[0]) + half_w, half_h * ndc.y() + double(viewport[1]) + half_h };
    };

    static Vec4d model_to_clip(const Vec3d& model, const Transform3d& world_matrix, const Matrix4d& projection_view_matrix) {
        return world_to_clip(model_to_world(model, world_matrix), projection_view_matrix);
    }

    static Vec3d model_to_ndc(const Vec3d& model, const Transform3d& world_matrix, const Matrix4d& projection_view_matrix) {
        return clip_to_ndc(world_to_clip(model_to_world(model, world_matrix), projection_view_matrix));
    }

    static Vec2d model_to_ss(const Vec3d& model, const Transform3d& world_matrix, const Matrix4d& projection_view_matrix, const std::array<int, 4>& viewport) {
        return ndc_to_ss(clip_to_ndc(world_to_clip(model_to_world(model, world_matrix), projection_view_matrix)), viewport);
    }

    static const Matrix4d& ndc_to_ss_matrix(const std::array<int, 4>& viewport) {
        update(viewport);
        return s_cache.ndc_to_ss_matrix;
    }

    static const Transform3d ndc_to_ss_matrix_inverse(const std::array<int, 4>& viewport) {
        update(viewport);
        return s_cache.ndc_to_ss_matrix_inverse;
    }

private:
    static void update(const std::array<int, 4>& viewport) {
        if (s_cache.viewport == viewport)
            return;

        const double half_w = 0.5 * double(viewport[2]);
        const double half_h = 0.5 * double(viewport[3]);
        s_cache.ndc_to_ss_matrix << half_w, 0.0, 0.0, double(viewport[0]) + half_w,
                                    0.0, half_h, 0.0, double(viewport[1]) + half_h,
                                    0.0, 0.0, 1.0, 0.0,
                                    0.0, 0.0, 0.0, 1.0;

        s_cache.ndc_to_ss_matrix_inverse = s_cache.ndc_to_ss_matrix.inverse();
        s_cache.viewport = viewport;
    }
};

DimensioningHelper::Cache DimensioningHelper::s_cache = { { 0, 0, 0, 0 }, Matrix4d::Identity(), Transform3d::Identity() };

void GLGizmoMeasure::render_dimensioning()
{
    static SelectedFeatures last_selected_features;

    if (!m_selected_features.first.feature.has_value() || !m_selected_features.second.feature.has_value())
        return;

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    auto point_point = [this, shader](const Vec3d& v1, const Vec3d& v2) {
        if (v1.isApprox(v2))
            return;

        const Camera& camera = wxGetApp().plater()->get_camera();
        const Matrix4d projection_view_matrix = camera.get_projection_matrix().matrix() * camera.get_view_matrix().matrix();
        const std::array<int, 4>& viewport = camera.get_viewport();

        // screen coordinates
        const Vec2d v1ss = DimensioningHelper::model_to_ss(v1, m_volume_matrix, projection_view_matrix, viewport);
        const Vec2d v2ss = DimensioningHelper::model_to_ss(v2, m_volume_matrix, projection_view_matrix, viewport);

        if (v1ss.isApprox(v2ss))
            return;

        const Vec2d v12ss = v2ss - v1ss;
        const double v12ss_len = v12ss.norm();

        const bool overlap = v12ss_len - 2.0 * TRIANGLE_HEIGHT < 0.0;

        const auto q12ss = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitX(), Vec3d(v12ss.x(), v12ss.y(), 0.0));
        const auto q21ss = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitX(), Vec3d(-v12ss.x(), -v12ss.y(), 0.0));

        shader->set_uniform("projection_matrix", Transform3d::Identity());

        const Vec3d v1ss_3 = { v1ss.x(), v1ss.y(), 0.0 };
        const Vec3d v2ss_3 = { v2ss.x(), v2ss.y(), 0.0 };

        const Transform3d ss_to_ndc_matrix = DimensioningHelper::ndc_to_ss_matrix_inverse(viewport);

        // stem
        shader->set_uniform("view_model_matrix", overlap ?
            ss_to_ndc_matrix * Geometry::translation_transform(v1ss_3) * q12ss * Geometry::translation_transform(-2.0 * TRIANGLE_HEIGHT * Vec3d::UnitX()) * Geometry::scale_transform({ v12ss_len + 4.0 * TRIANGLE_HEIGHT, 1.0f, 1.0f }) :
            ss_to_ndc_matrix * Geometry::translation_transform(v1ss_3) * q12ss * Geometry::scale_transform({ v12ss_len, 1.0f, 1.0f }));
            m_dimensioning.line.render();

        // arrow 1
        shader->set_uniform("view_model_matrix", overlap ?
            ss_to_ndc_matrix * Geometry::translation_transform(v1ss_3) * q12ss :
            ss_to_ndc_matrix * Geometry::translation_transform(v1ss_3) * q21ss);
        m_dimensioning.triangle.render();

        // arrow 2
        shader->set_uniform("view_model_matrix", overlap ?
            ss_to_ndc_matrix * Geometry::translation_transform(v2ss_3) * q21ss :
            ss_to_ndc_matrix * Geometry::translation_transform(v2ss_3) * q12ss);
        m_dimensioning.triangle.render();
    };

    auto point_edge = [this, shader](const Measure::SurfaceFeature& f1, const Measure::SurfaceFeature& f2) {
        assert(f1.get_type() == Measure::SurfaceFeatureType::Point && f2.get_type() == Measure::SurfaceFeatureType::Edge);
        const std::pair<Vec3d, Vec3d> e = f2.get_edge();
        const Vec3d v_proj    = m_measurement_result.distance_infinite->to;

        const Vec3d e1e2 = e.second - e.first;
        const Vec3d v_proje1 = v_proj - e.first;
        const bool on_e1_side = v_proje1.dot(e1e2) < -EPSILON;
        const bool on_e2_side = !on_e1_side && v_proje1.norm() > e1e2.norm();
        if (on_e1_side || on_e2_side) {
            const Camera& camera = wxGetApp().plater()->get_camera();
            const Matrix4d projection_view_matrix = camera.get_projection_matrix().matrix() * camera.get_view_matrix().matrix();
            const std::array<int, 4>& viewport = camera.get_viewport();
            const Transform3d ss_to_ndc_matrix = DimensioningHelper::ndc_to_ss_matrix_inverse(viewport);

            const Vec2d v_projss = DimensioningHelper::model_to_ss(v_proj, m_volume_matrix, projection_view_matrix, viewport);
            auto render_extension = [this, &v_projss, &projection_view_matrix, &viewport, &ss_to_ndc_matrix, shader](const Vec3d& p) {
                const Vec2d pss = DimensioningHelper::model_to_ss(p, m_volume_matrix, projection_view_matrix, viewport);
                if (!pss.isApprox(v_projss)) {
                    const Vec2d pv_projss = v_projss - pss;
                    const double pv_projss_len = pv_projss.norm();

                    const auto q = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitX(), Vec3d(pv_projss.x(), pv_projss.y(), 0.0));

                    shader->set_uniform("projection_matrix", Transform3d::Identity());
                    shader->set_uniform("view_model_matrix", ss_to_ndc_matrix * Geometry::translation_transform({ pss.x(), pss.y(), 0.0 }) * q *
                        Geometry::scale_transform({ pv_projss_len, 1.0f, 1.0f }));
                    m_dimensioning.line.render();
                }
            };

            render_extension(on_e1_side ? e.first : e.second);
        }
    };

    auto arc_edge_edge = [this, shader](const Measure::SurfaceFeature& f1, const Measure::SurfaceFeature& f2, double radius = 0.0) {
        assert(f1.get_type() == Measure::SurfaceFeatureType::Edge && f2.get_type() == Measure::SurfaceFeatureType::Edge);
        const Measure::MeasurementResult res = Measure::get_measurement(f1, f2);
        const double angle    = res.angle->angle;
        const Vec3d  center   = res.angle->center;
        const std::pair<Vec3d, Vec3d> e1 = res.angle->e1;
        const std::pair<Vec3d, Vec3d> e2 = res.angle->e2;
        const double calc_radius = res.angle->radius;
        const bool   coplanar = res.angle->coplanar;

        if (calc_radius == 0.0)
            return;

        const double draw_radius = (radius > 0.0) ? radius : calc_radius;

        const Vec3d e1_unit = Measure::edge_direction(e1);
        const Vec3d e2_unit = Measure::edge_direction(e2);

        if (!m_dimensioning.arc.is_initialized()) {
            const unsigned int resolution = std::max<unsigned int>(2, 64 * angle / double(PI));
            GLModel::Geometry init_data;
            init_data.format = { GLModel::Geometry::EPrimitiveType::LineStrip, GLModel::Geometry::EVertexLayout::P3 };
            init_data.color = ColorRGBA::WHITE();
            init_data.reserve_vertices(resolution + 1);
            init_data.reserve_indices(resolution + 1);

            // vertices + indices
            const Vec3d normal = e1_unit.cross(e2_unit).normalized();
            const double step = angle / double(resolution);
            for (unsigned int i = 0; i <= resolution; ++i) {
                const double a = step * double(i);
                const Vec3d v = draw_radius * (Eigen::Quaternion<double>(Eigen::AngleAxisd(a, normal)) * e1_unit);
                init_data.add_vertex((Vec3f)v.cast<float>());
                init_data.add_index(i);
            }

            m_dimensioning.arc.init_from(std::move(init_data));
        }

        // render arc
        const Camera& camera = wxGetApp().plater()->get_camera();
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        shader->set_uniform("view_model_matrix", camera.get_view_matrix() * m_volume_matrix * Geometry::translation_transform(center));
        m_dimensioning.arc.render();

        // render edge 1 extension
        const Vec3d e11e12 = e1.second - e1.first;
        const Vec3d e11center = center - e1.first;
        const double e11center_len = e11center.norm();
        if (e11center_len > EPSILON && e11center.dot(e11e12) < 0.0) {
            const Camera& camera = wxGetApp().plater()->get_camera();
            shader->set_uniform("projection_matrix", camera.get_projection_matrix());
            shader->set_uniform("view_model_matrix", camera.get_view_matrix() * m_volume_matrix * Geometry::translation_transform(center) *
                Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitX(), Measure::edge_direction(e1.first, e1.second)) *
                Geometry::scale_transform({ e11center_len, 1.0f, 1.0f }));
            m_dimensioning.line.render();
        }

        // render edge 2 extension
        const Vec3d e21center = center - e2.first;
        const double e21center_len = e21center.norm();
        if (e21center_len > EPSILON) {
            const Camera& camera = wxGetApp().plater()->get_camera();
            shader->set_uniform("projection_matrix", camera.get_projection_matrix());
            shader->set_uniform("view_model_matrix", camera.get_view_matrix() * m_volume_matrix * Geometry::translation_transform(center) *
                Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitX(), Measure::edge_direction(e2.first, e2.second)) *
                Geometry::scale_transform({ (coplanar && radius > 0.0) ? e21center_len : draw_radius, 1.0f, 1.0f }));
            m_dimensioning.line.render();
        }
    };

    auto arc_edge_plane = [this, arc_edge_edge](const Measure::SurfaceFeature& f1, const Measure::SurfaceFeature& f2) {
        assert(f1.get_type() == Measure::SurfaceFeatureType::Edge && f2.get_type() == Measure::SurfaceFeatureType::Plane);
        const Measure::MeasurementResult res = Measure::get_measurement(f1, f2);
        const std::pair<Vec3d, Vec3d> e1 = res.angle->e1;
        const std::pair<Vec3d, Vec3d> e2 = res.angle->e2;
        const double calc_radius = res.angle->radius;
        if (calc_radius == 0.0)
            return;

        arc_edge_edge(Measure::SurfaceFeature(Measure::SurfaceFeatureType::Edge, e1.first, e1.second),
            Measure::SurfaceFeature(Measure::SurfaceFeatureType::Edge, e2.first, e2.second), calc_radius);
    };

    shader->start_using();

    if (!m_dimensioning.line.is_initialized()) {
        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3 };
        init_data.color = ColorRGBA::WHITE();
        init_data.reserve_vertices(2);
        init_data.reserve_indices(2);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex(Vec3f(1.0f, 0.0f, 0.0f));

        // indices
        init_data.add_line(0, 1);

        m_dimensioning.line.init_from(std::move(init_data));
    }

    if (!m_dimensioning.triangle.is_initialized()) {
        GLModel::Geometry init_data;
        init_data.format = { GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3 };
        init_data.color = ColorRGBA::WHITE();
        init_data.reserve_vertices(3);
        init_data.reserve_indices(3);

        // vertices
        init_data.add_vertex(Vec3f(0.0f, 0.0f, 0.0f));
        init_data.add_vertex(Vec3f(-TRIANGLE_HEIGHT, 0.5f * TRIANGLE_BASE, 0.0f));
        init_data.add_vertex(Vec3f(-TRIANGLE_HEIGHT, -0.5f * TRIANGLE_BASE, 0.0f));

        // indices
        init_data.add_triangle(0, 1, 2);

        m_dimensioning.triangle.init_from(std::move(init_data));
    }

    if (last_selected_features != m_selected_features)
        m_dimensioning.arc.reset();

    glsafe(::glDisable(GL_DEPTH_TEST));

    if (m_selected_features.second.feature.has_value()) {
        // Render the arrow between the points that the backend passed:
        const Measure::DistAndPoints& dap = m_measurement_result.distance_infinite.has_value()
                                          ? *m_measurement_result.distance_infinite
                                          : *m_measurement_result.distance_strict;
        point_point(dap.from, dap.to);

        const Measure::SurfaceFeature* f1 = &(*m_selected_features.first.feature);
        const Measure::SurfaceFeature* f2 = &(*m_selected_features.second.feature);
        Measure::SurfaceFeatureType ft1 = f1->get_type();
        Measure::SurfaceFeatureType ft2 = f2->get_type();

        // Order features by type so following conditions are simple.
        if (ft1 > ft2) {
            std::swap(ft1, ft2);
            std::swap(f1, f2);
        }

        // Where needed, draw also the extension of the edge to where the dist is measured:
        if (ft1 == Measure::SurfaceFeatureType::Point && ft2 == Measure::SurfaceFeatureType::Edge)
            point_edge(*f1, *f2);

        // Now if there is an angle to show, draw the arc:
        if (ft1 == Measure::SurfaceFeatureType::Edge && ft2 == Measure::SurfaceFeatureType::Edge)
            arc_edge_edge(*f1, *f2);
        else if (ft1 == Measure::SurfaceFeatureType::Edge && ft2 == Measure::SurfaceFeatureType::Plane)
            arc_edge_plane(*f1, *f2);
    }
    
    glsafe(::glEnable(GL_DEPTH_TEST));

    shader->stop_using();
}

static void add_row_to_table(std::function<void(void)> col_1 = nullptr, std::function<void(void)> col_2 = nullptr)
{
    assert(col_1 != nullptr && col_2 != nullptr);
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    col_1();
    ImGui::TableSetColumnIndex(1);
    col_2();
}

static void add_strings_row_to_table(ImGuiWrapper& imgui, const std::string& col_1, const ImVec4& col_1_color, const std::string& col_2, const ImVec4& col_2_color)
{
    add_row_to_table([&]() { imgui.text_colored(col_1_color, col_1); }, [&]() { imgui.text_colored(col_2_color, col_2); });
};

static std::string format_double(double value)
{
    char buf[1024];
    sprintf(buf, "%.3f", value);
    return std::string(buf);
}

static std::string format_vec3(const Vec3d& v)
{
    char buf[1024];
    sprintf(buf, "X: %.3f, Y: %.3f, Z: %.3f", v.x(), v.y(), v.z());
    return std::string(buf);
}

#if ENABLE_MEASURE_GIZMO_DEBUG
void GLGizmoMeasure::render_debug_dialog()
{
    auto add_feature_data = [this](const SelectedFeatures::Item& item) {
        add_strings_row_to_table(*m_imgui, "Type", ImGuiWrapper::COL_ORANGE_LIGHT, item.source, ImGui::GetStyleColorVec4(ImGuiCol_Text));
        switch (item.feature->get_type())
        {
        case Measure::SurfaceFeatureType::Point:
        {
            add_strings_row_to_table(*m_imgui, "m_pt1", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(item.feature->get_point()), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            break;
        }
        case Measure::SurfaceFeatureType::Edge:
        {
            auto [from, to] = item.feature->get_edge();
            add_strings_row_to_table(*m_imgui, "m_pt1", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(from), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            add_strings_row_to_table(*m_imgui, "m_pt2", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(to), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            break;
        }
        case Measure::SurfaceFeatureType::Plane:
        {
            auto [idx, normal, origin] = item.feature->get_plane();
            add_strings_row_to_table(*m_imgui, "m_pt1", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            add_strings_row_to_table(*m_imgui, "m_pt2", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(origin), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            add_strings_row_to_table(*m_imgui, "m_value", ImGuiWrapper::COL_ORANGE_LIGHT, format_double(idx), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            break;
        }
        case Measure::SurfaceFeatureType::Circle:
        {
            auto [center, radius, normal] = item.feature->get_circle();
            add_strings_row_to_table(*m_imgui, "m_pt1", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(center), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            add_strings_row_to_table(*m_imgui, "m_pt2", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            add_strings_row_to_table(*m_imgui, "m_value", ImGuiWrapper::COL_ORANGE_LIGHT, format_double(radius), ImGui::GetStyleColorVec4(ImGuiCol_Text));
            break;
        }
        }
        std::optional<Vec3d> extra_point = item.feature->get_extra_point();
        if (extra_point.has_value())
            add_strings_row_to_table(*m_imgui, "m_pt3", ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(*extra_point), ImGui::GetStyleColorVec4(ImGuiCol_Text));
    };

    m_imgui->begin(_L("Measure tool debug"), ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
    if (!m_selected_features.first.feature.has_value() && !m_selected_features.second.feature.has_value())
        m_imgui->text("Empty selection");
    else {
        const ImGuiTableFlags flags = ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersH;
        if (m_selected_features.first.feature.has_value()) {
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, "Selection 1");
            if (ImGui::BeginTable("Selection 1", 2, flags)) {
                add_feature_data(m_selected_features.first);
                ImGui::EndTable();
            }
        }
        if (m_selected_features.second.feature.has_value()) {
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, "Selection 2");
            if (ImGui::BeginTable("Selection 2", 2, flags)) {
                add_feature_data(m_selected_features.second);
                ImGui::EndTable();
            }
        }
    }
    m_imgui->end();
}
#endif // ENABLE_MEASURE_GIZMO_DEBUG

void GLGizmoMeasure::on_render_input_window(float x, float y, float bottom_limit)
{
    static std::optional<Measure::SurfaceFeature> last_feature;
    static EMode last_mode = EMode::BasicSelection;
    static SelectedFeatures last_selected_features;

    static float last_y = 0.0f;
    static float last_h = 0.0f;

    m_imgui->begin(_L("Measure tool"), ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

    // adjust window position to avoid overlap the view toolbar
    const float win_h = ImGui::GetWindowHeight();
    y = std::min(y, bottom_limit - win_h);
    ImGui::SetWindowPos(ImVec2(x, y), ImGuiCond_Always);
    if (last_h != win_h || last_y != y) {
        // ask canvas for another frame to render the window in the correct position
        m_imgui->set_requires_extra_frame();
        if (last_h != win_h)
            last_h = win_h;
        if (last_y != y)
            last_y = y;
    }

    if (ImGui::BeginTable("Commands", 2)) {
        add_row_to_table(
            [this]() {
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Left mouse button"));
            },
            [this]() {
                m_imgui->text_colored(ImGui::GetStyleColorVec4(ImGuiCol_Text),
                    m_selected_features.second.feature.has_value() ?
                    ((m_mode == EMode::BasicSelection) ? _u8L("Select/Unselect feature") : _u8L("Select/Unselect point")) :
                    ((m_mode == EMode::BasicSelection) ? _u8L("Select feature") : _u8L("Select point")));
                ImGui::SameLine();
                const ImVec2 pos = ImGui::GetCursorScreenPos();
                const float rect_size = ImGui::GetTextLineHeight();
                ImGui::GetWindowDrawList()->AddRectFilled({ pos.x + 1.0f, pos.y + 1.0f }, { pos.x + rect_size, pos.y + rect_size },
                    ImGuiWrapper::to_ImU32(m_selected_features.first.feature.has_value() ? SELECTED_2ND_COLOR : SELECTED_1ST_COLOR));
                ImGui::Dummy(ImVec2(rect_size, rect_size));
            }
            );

        if (m_selected_features.first.feature.has_value())
            add_strings_row_to_table(*m_imgui, CTRL_STR + "+" + _u8L("Right mouse button"), ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Restart selection"), ImGui::GetStyleColorVec4(ImGuiCol_Text));

        if (m_mode == EMode::BasicSelection && m_hover_id != -1)
            add_strings_row_to_table(*m_imgui, CTRL_STR, ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Enable point selection"), ImGui::GetStyleColorVec4(ImGuiCol_Text));
        ImGui::EndTable();
    }

    const bool use_inches = wxGetApp().app_config->get("use_inches") == "1";
    const std::string units = use_inches ? _u8L(" (in)") : _u8L(" (mm)");

    if (m_curr_feature.has_value()) {
        const Measure::SurfaceFeatureType feature_type = m_curr_feature->get_type();
        if (m_mode == EMode::BasicSelection) {
            if (feature_type != Measure::SurfaceFeatureType::Undef) {
                ImGui::Separator();
                m_imgui->text(surface_feature_type_as_string(feature_type) + ":");
                if (ImGui::BeginTable("Data", 2)) {
                    switch (feature_type)
                    {
                    default: { assert(false); break; }
                    case Measure::SurfaceFeatureType::Point:
                    {
                        Vec3d position = m_volume_matrix * m_curr_feature->get_point();
                        if (use_inches)
                            position = ObjectManipulation::mm_to_in * position;
                        add_strings_row_to_table(*m_imgui, _u8L("Position"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(position), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Edge:
                    {
                        auto [from, to] = m_curr_feature->get_edge();
                        from = m_volume_matrix * from;
                        to   = m_volume_matrix * to;
                        if (use_inches) {
                            from = ObjectManipulation::mm_to_in * from;
                            to   = ObjectManipulation::mm_to_in * to;
                        }
                        add_strings_row_to_table(*m_imgui, _u8L("From"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(from), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(*m_imgui, _u8L("To"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(to), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(*m_imgui, _u8L("Length") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double((to - from).norm()), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Circle:
                    {
                        auto [center, radius, normal] = m_curr_feature->get_circle();
                        center = m_volume_matrix * center;
                        normal = m_volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                        if (use_inches) {
                            center = ObjectManipulation::mm_to_in * center;
                            radius = ObjectManipulation::mm_to_in * radius;
                        }
                        add_strings_row_to_table(*m_imgui, _u8L("Center"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(center), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(*m_imgui, _u8L("Radius") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(radius), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(*m_imgui, _u8L("Normal"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Plane:
                    {
                        auto [idx, normal, origin] = m_curr_feature->get_plane();
                        origin = m_volume_matrix * origin;
                        normal = m_volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                        if (use_inches)
                            origin = ObjectManipulation::mm_to_in * origin;
                        add_strings_row_to_table(*m_imgui, _u8L("Origin"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(origin), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(*m_imgui, _u8L("Normal"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    }
                    ImGui::EndTable();
                }
            }
        }
        else if (m_mode == EMode::ExtendedSelection) {
            if (m_hover_id != -1 && m_curr_point_on_feature_position.has_value()) {
                ImGui::Separator();
                m_imgui->text(point_on_feature_type_as_string(feature_type, m_hover_id) + ":");
                if (ImGui::BeginTable("Data", 2)) {
                    Vec3d position = m_volume_matrix * *m_curr_point_on_feature_position;
                    if (use_inches)
                        position = ObjectManipulation::mm_to_in * position;
                    add_strings_row_to_table(*m_imgui, _u8L("Position"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(position), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::EndTable();
                }
            }
        }
    }

    ImGui::Separator();
    const ImGuiTableFlags flags = ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersH;
    if (ImGui::BeginTable("Selection", 2, flags)) {
        add_strings_row_to_table(*m_imgui, _u8L("Selection") + " 1:", ImGuiWrapper::to_ImVec4(SELECTED_1ST_COLOR), m_selected_features.first.feature.has_value() ?
            m_selected_features.first.source : _u8L("None"), ImGuiWrapper::to_ImVec4(SELECTED_1ST_COLOR));
        add_strings_row_to_table(*m_imgui, _u8L("Selection") + " 2:", ImGuiWrapper::to_ImVec4(SELECTED_2ND_COLOR), m_selected_features.second.feature.has_value() ?
            m_selected_features.second.source : _u8L("None"), ImGuiWrapper::to_ImVec4(SELECTED_2ND_COLOR));
        ImGui::EndTable();
    }

    //if (m_selected_features.first.feature.has_value()) {
    //    if (m_imgui->button(_u8L("Restart"))) {
    //        m_selected_features.reset();
    //        m_selection_raycasters.clear();
    //        m_imgui->set_requires_extra_frame();
    //    }
    //}

    auto add_measure_row_to_table = [this](const std::string& col_1, const ImVec4& col_1_color, const std::string& col_2, const ImVec4& col_2_color) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        m_imgui->text_colored(col_1_color, col_1);
        ImGui::TableSetColumnIndex(1);
        m_imgui->text_colored(col_2_color, col_2);
        ImGui::TableSetColumnIndex(2);
        ImGuiIO& io = ImGui::GetIO();
        const ImTextureID tex_id = io.Fonts->TexID;
        assert(io.Fonts->TexWidth > 0 && io.Fonts->TexHeight > 0);
        float inv_tex_w = 1.0f / float(io.Fonts->TexWidth);
        float inv_tex_h = 1.0f / float(io.Fonts->TexHeight);
        const ImFontAtlasCustomRect* const rect = m_imgui->GetTextureCustomRect(ImGui::ClipboardBtnIcon);
        const ImVec2 size = { float(rect->Width), float(rect->Height) };
        const ImVec2 uv0 = ImVec2(float(rect->X) * inv_tex_w, float(rect->Y) * inv_tex_h);
        const ImVec2 uv1 = ImVec2(float(rect->X + rect->Width) * inv_tex_w, float(rect->Y + rect->Height) * inv_tex_h);
        ImGui::PushStyleColor(ImGuiCol_Button,        { 0.25f, 0.25f, 0.25f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, { 0.4f, 0.4f, 0.4f, 1.0f });
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  { 0.25f, 0.25f, 0.25f, 1.0f });
        if (m_imgui->image_button(tex_id, size, uv0, uv1)) {
            wxTheClipboard->Open();
            wxTheClipboard->SetData(new wxTextDataObject(col_1 + ": " + col_2));
            wxTheClipboard->Close();
        }
        ImGui::PopStyleColor(3);
        if (ImGui::IsItemHovered()) {
            const float max_tooltip_width = ImGui::GetFontSize() * 20.0f;
            m_imgui->tooltip(into_u8(_L("Copy to clipboard")).c_str(), max_tooltip_width);
        }
    };

    if (m_selected_features.second.feature.has_value()) {
        const Measure::MeasurementResult& measure = m_measurement_result;

        ImGui::Separator();
        if (measure.has_any_data()) {
            m_imgui->text(_u8L("Measure") + ":");
            if (ImGui::BeginTable("Measure", 3)) {
                if (measure.angle.has_value()) {
                    ImGui::PushID((void*)(intptr_t)1);
                    add_measure_row_to_table(_u8L("Angle") + _u8L(" (°)"), ImGuiWrapper::COL_ORANGE_LIGHT, format_double(Geometry::rad2deg(measure.angle->angle)),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::PopID();
                }
                if (measure.distance_infinite.has_value()) {
                    double distance = measure.distance_infinite->dist;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    ImGui::PushID((void*)(intptr_t)2);
                    add_measure_row_to_table(_u8L("Distance Infinite") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::PopID();
                }
                if (measure.distance_strict.has_value()) {
                    double distance = measure.distance_strict->dist;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    ImGui::PushID((void*)(intptr_t)3);
                    add_measure_row_to_table(_u8L("Distance Strict") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::PopID();
                }
                if (measure.distance_xyz.has_value()) {
                    Vec3d distance = *measure.distance_xyz;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    ImGui::PushID((void*)(intptr_t)4);
                    add_measure_row_to_table(_u8L("Distance XYZ") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
        }
        else
            m_imgui->text(_u8L("No measure available"));
    }

    if (last_feature != m_curr_feature || last_mode != m_mode || last_selected_features != m_selected_features) {
        // the dialog may have changed its size, ask for an extra frame to render it properly
        last_feature = m_curr_feature;
        last_mode = m_mode;
        last_selected_features = m_selected_features;
        m_imgui->set_requires_extra_frame();
    }

    m_imgui->end();
}

void GLGizmoMeasure::on_register_raycasters_for_picking()
{
    // the features are rendered on top of the scene, so the raytraced picker should take it into account
    m_parent.set_raycaster_gizmos_on_top(true);
}

void GLGizmoMeasure::on_unregister_raycasters_for_picking()
{
    m_parent.remove_raycasters_for_picking(SceneRaycaster::EType::Gizmo);
    m_parent.set_raycaster_gizmos_on_top(false);
    m_raycasters.clear();
    m_selection_raycasters.clear();
}

} // namespace GUI
} // namespace Slic3r

#endif // ENABLE_MEASURE_GIZMO
