// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/GUI_ObjectManipulation.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/PresetBundle.hpp"

#include <numeric>

#include <GL/glew.h>

#if ENABLE_MEASURE_GIZMO

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA SELECTED_1ST_COLOR   = { 0.25f, 0.75f, 0.75f, 1.0f };
static const Slic3r::ColorRGBA SELECTED_2ND_COLOR   = { 0.75f, 0.25f, 0.75f, 1.0f };

static const int POINT_ID         = 100;
static const int EDGE_ID          = 200;
static const int CIRCLE_ID        = 300;
static const int PLANE_ID         = 400;

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
    case Measure::SurfaceFeatureType::Undef:  { return L("Undefined"); }
    case Measure::SurfaceFeatureType::Point:  { return L("Vertex"); }
    case Measure::SurfaceFeatureType::Edge:   { return L("Edge"); }
    case Measure::SurfaceFeatureType::Circle: { return L("Circle"); }
    case Measure::SurfaceFeatureType::Plane:  { return L("Plane"); }
    }
}

static std::string point_on_feature_type_as_string(Measure::SurfaceFeatureType type, int hover_id)
{
    std::string ret;
    switch (type) {
    case Measure::SurfaceFeatureType::Point:  { ret = _u8L("Vertex"); break; }
    case Measure::SurfaceFeatureType::Edge:   { ret = _u8L("Point on edge"); break; }
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
            m_mouse_left_down = true;
            
            auto item_from_feature = [this]() {
                const SelectedFeatures::Item item = { m_mode,
                    (m_mode == EMode::ExtendedSelection) ? point_on_feature_type_as_string(m_curr_feature->get_type(), m_hover_id) : surface_feature_type_as_string(m_curr_feature->get_type()),
                    (m_mode == EMode::ExtendedSelection) ? Measure::SurfaceFeature(*m_curr_point_on_feature_position) : m_curr_feature };
                return item;
            };

            if (m_selected_features.first.feature.has_value()) {
                const SelectedFeatures::Item item = item_from_feature();
                if (m_selected_features.first != item)
                    m_selected_features.second = item;
            }
            else
                m_selected_features.first = item_from_feature();

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
        m_scene_raycaster_state.clear();
        m_scene_raycasters = m_parent.get_raycasters_for_picking(SceneRaycaster::EType::Volume);
        if (m_scene_raycasters != nullptr) {
            for (const auto& r : *m_scene_raycasters) {
                m_scene_raycaster_state.emplace_back(r->is_active());
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
    return (wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() == ptSLA) ?
        selection.is_single_full_instance() :
        selection.is_single_volume() || selection.is_single_volume_instance();
}

void GLGizmoMeasure::on_render()
{
    // do not render if the user is panning/rotating the 3d scene
    if (m_parent.is_mouse_dragging())
        return;

    const Selection& selection = m_parent.get_selection();

    if ((wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() == ptSLA && selection.is_single_full_instance()) ||
        (selection.is_single_volume() || selection.is_single_volume_instance())) {
        update_if_needed();

        const Transform3d& model_matrix = selection.get_first_volume()->world_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const float inv_zoom = (float)camera.get_inv_zoom();

        Vec3f position_on_model;
        Vec3f normal_on_model;
        size_t model_facet_idx;
        const bool mouse_on_object = m_c->raycaster()->raycasters().front()->unproject_on_mesh(m_mouse_pos, model_matrix, camera, position_on_model, normal_on_model, nullptr, &model_facet_idx);
        const bool is_hovering_on_locked_feature = m_mode == EMode::ExtendedSelection && m_hover_id != -1;

        if (m_mode == EMode::BasicSelection) {
            std::optional<Measure::SurfaceFeature> curr_feature = mouse_on_object ? m_measuring->get_feature(model_facet_idx, position_on_model.cast<double>()) : std::nullopt;
            m_curr_point_on_feature_position.reset();
            if (m_curr_feature != curr_feature) {
                GLGizmoMeasure::on_unregister_raycasters_for_picking();
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
                m_curr_point_on_feature_position = position_on_feature(POINT_ID, camera);
                break;
            }
            case Measure::SurfaceFeatureType::Edge:
            {
                m_curr_point_on_feature_position = position_on_feature(EDGE_ID, camera, [](const Vec3f& v) { return Vec3f(0.0f, 0.0f, v.z()); });
                break;
            }
            case Measure::SurfaceFeatureType::Plane:
            {
                m_curr_point_on_feature_position = position_on_feature(PLANE_ID, camera);
                break;
            }
            case Measure::SurfaceFeatureType::Circle:
            {
                const auto [center, radius, normal] = m_curr_feature->get_circle();
                if (m_hover_id == POINT_ID)
                    m_curr_point_on_feature_position = model_matrix * center;
                else {
                    const float r = radius; // needed for the following lambda
                    m_curr_point_on_feature_position = position_on_feature(CIRCLE_ID, camera, [r](const Vec3f& v) {
                        float angle = std::atan2(v.y(), v.x());
                        if (angle < 0.0f)
                            angle += 2.0f * float(M_PI);
                        return float(r) * Vec3f(std::cos(angle), std::sin(angle), 0.0f);
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
                const auto& [center, radius, n] = feature.get_circle();
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
                const Transform3d circle_matrix = model_matrix * Geometry::translation_transform(center);
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
                auto q = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), end - start);
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(start) * q *
                    Geometry::scale_transform({ (double)inv_zoom, (double)inv_zoom, (end - start).norm() });
                set_matrix_uniforms(feature_matrix);
                m_cylinder.model.set_color(colors.front());
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
                case Measure::SurfaceFeatureType::Circle:
                {
                    colors.emplace_back((m_hover_id == POINT_ID) ? hover_selection_color() : hovering_color());
                    colors.emplace_back(hovering_color());
                    break;
                }
                case Measure::SurfaceFeatureType::Edge:
                case Measure::SurfaceFeatureType::Plane:
                {
                    colors.emplace_back(hovering_color());
                    break;
                }
                }
            }

            render_feature(*m_curr_feature, colors, model_matrix, inv_zoom, true);
        }

        if (m_selected_features.first.feature.has_value() && (!m_curr_feature.has_value() || *m_curr_feature != *m_selected_features.first.feature)) {
            std::vector<ColorRGBA> colors;
            colors.emplace_back(SELECTED_1ST_COLOR);
            render_feature(*m_selected_features.first.feature, colors,
                (m_selected_features.first.mode == EMode::BasicSelection) ? model_matrix : Transform3d::Identity(), inv_zoom, false);
        }
        if (m_selected_features.second.feature.has_value() && (!m_curr_feature.has_value() || *m_curr_feature != *m_selected_features.second.feature)) {
            std::vector<ColorRGBA> colors;
            colors.emplace_back(SELECTED_2ND_COLOR);
            render_feature(*m_selected_features.second.feature, colors,
                (m_selected_features.second.mode == EMode::BasicSelection) ? model_matrix : Transform3d::Identity(), inv_zoom, false);
        }

        if (is_hovering_on_locked_feature && m_curr_point_on_feature_position.has_value()) {
            if (m_hover_id != POINT_ID) {
                const Transform3d matrix = Geometry::translation_transform(*m_curr_point_on_feature_position) * Geometry::scale_transform(inv_zoom);
                set_matrix_uniforms(matrix);
                m_sphere.model.set_color(hover_selection_color());
                m_sphere.model.render();
            }
        }

        shader->stop_using();
    }
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
    if (m_scene_raycasters != nullptr) {
        for (auto r : *m_scene_raycasters) {
            r->set_active(false);
        }
    }
}

void GLGizmoMeasure::restore_scene_raycasters_state()
{
    if (m_scene_raycasters != nullptr) {
        assert(m_scene_raycasters->size() == m_scene_raycaster_state.size());
        for (size_t i = 0; i < m_scene_raycasters->size(); ++i) {
            (*m_scene_raycasters)[i]->set_active(m_scene_raycaster_state[i]);
        }
    }
}

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

    auto add_row_to_table = [this](std::function<void(void)> col_1 = nullptr, std::function<void(void)> col_2 = nullptr) {
        assert(col_1 != nullptr && col_2 != nullptr);
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        col_1();
        ImGui::TableSetColumnIndex(1);
        col_2();
    };

    auto add_strings_row_to_table = [this, add_row_to_table](const std::string& col_1, const ImVec4& col_1_color, const std::string& col_2, const ImVec4& col_2_color) {
        add_row_to_table(
            [this, &col_1, &col_1_color]() { m_imgui->text_colored(col_1_color, col_1); },
            [this, &col_2, &col_2_color]() { m_imgui->text_colored(col_2_color, col_2); }
            );
    };

    auto format_double = [](double value) {
        char buf[1024];
        sprintf(buf, "%.3f", value);
        return std::string(buf);
    };

    auto format_vec3 = [](const Vec3d& v) {
        char buf[1024];
        sprintf(buf, "X: %.3f, Y: %.3f, Z: %.3f", v.x(), v.y(), v.z());
        return std::string(buf);
    };

    if (ImGui::BeginTable("Commands", 2)) {
        add_row_to_table(
            [this]() {
            m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Left mouse button"));
            },
            [this]() {
                m_imgui->text_colored(ImGui::GetStyleColorVec4(ImGuiCol_Text), (m_mode == EMode::BasicSelection) ? _u8L("Select feature") : _u8L("Select point"));
                ImGui::SameLine();
                const ImVec2 pos = ImGui::GetCursorScreenPos();
                const float rect_size = ImGui::GetTextLineHeight();
                ImGui::GetWindowDrawList()->AddRectFilled({ pos.x + 1.0f, pos.y + 1.0f }, { pos.x + rect_size, pos.y + rect_size },
                    ImGuiWrapper::to_ImU32(m_selected_features.first.feature.has_value() ? SELECTED_2ND_COLOR : SELECTED_1ST_COLOR));
                ImGui::Dummy(ImVec2(rect_size, rect_size));
            }
            );

        if (m_selected_features.first.feature.has_value())
            add_strings_row_to_table(CTRL_STR + "+" + _u8L("Right mouse button"), ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Restart selection"), ImGui::GetStyleColorVec4(ImGuiCol_Text));

        if (m_mode == EMode::BasicSelection && m_hover_id != -1)
            add_strings_row_to_table(CTRL_STR, ImGuiWrapper::COL_ORANGE_LIGHT, _u8L("Enable point selection"), ImGui::GetStyleColorVec4(ImGuiCol_Text));
        ImGui::EndTable();
    }

    const bool use_inches = wxGetApp().app_config->get("use_inches") == "1";
    const std::string units = use_inches ? _u8L(" (in)") : _u8L(" (mm)");

    if (m_curr_feature.has_value()) {
        const Transform3d volume_matrix = m_parent.get_selection().get_first_volume()->world_matrix();
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
                        Vec3d position = volume_matrix * m_curr_feature->get_point();
                        if (use_inches)
                            position = ObjectManipulation::mm_to_in * position;
                        add_strings_row_to_table(_u8L("Position"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(position), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Edge:
                    {
                        auto [from, to] = m_curr_feature->get_edge();
                        from = volume_matrix * from;
                        to = volume_matrix * to;
                        if (use_inches) {
                            from = ObjectManipulation::mm_to_in * from;
                            to   = ObjectManipulation::mm_to_in * to;
                        }
                        add_strings_row_to_table(_u8L("From"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(from), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(_u8L("To"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(to), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(_u8L("Length") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double((to - from).norm()), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Circle:
                    {
                        auto [center, radius, normal] = m_curr_feature->get_circle();
                        center = volume_matrix * center;
                        normal = volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                        if (use_inches) {
                            center = ObjectManipulation::mm_to_in * center;
                            radius = ObjectManipulation::mm_to_in * radius;
                        }
                        add_strings_row_to_table(_u8L("Center"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(center), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(_u8L("Radius") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(radius), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(_u8L("Normal"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        break;
                    }
                    case Measure::SurfaceFeatureType::Plane:
                    {
                        auto [idx, normal, origin] = m_curr_feature->get_plane();
                        origin = volume_matrix * origin;
                        normal = volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                        if (use_inches)
                            origin = ObjectManipulation::mm_to_in * origin;
                        add_strings_row_to_table(_u8L("Origin"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(origin), ImGui::GetStyleColorVec4(ImGuiCol_Text));
                        add_strings_row_to_table(_u8L("Normal"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(normal), ImGui::GetStyleColorVec4(ImGuiCol_Text));
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
                    Vec3d position = ObjectManipulation::mm_to_in * *m_curr_point_on_feature_position;
                    if (use_inches)
                        position = ObjectManipulation::mm_to_in * position;
                    add_strings_row_to_table(_u8L("Position"), ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(position),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                    ImGui::EndTable();
                }
            }
        }
    }

    ImGui::Separator();
    const ImGuiTableFlags flags = ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersH;
    if (ImGui::BeginTable("Selection", 2, flags)) {
        add_strings_row_to_table(_u8L("Selection") + " 1:", ImGuiWrapper::to_ImVec4(SELECTED_1ST_COLOR), m_selected_features.first.feature.has_value() ?
            m_selected_features.first.source : _u8L("None"), ImGuiWrapper::to_ImVec4(SELECTED_1ST_COLOR));
        add_strings_row_to_table(_u8L("Selection") + " 2:", ImGuiWrapper::to_ImVec4(SELECTED_2ND_COLOR), m_selected_features.second.feature.has_value() ?
            m_selected_features.second.source : _u8L("None"), ImGuiWrapper::to_ImVec4(SELECTED_2ND_COLOR));
        ImGui::EndTable();
    }

    //if (m_selected_features.first.feature.has_value()) {
    //    if (m_imgui->button(_u8L("Restart"))) {
    //        m_selected_features.reset();
    //        m_imgui->set_requires_extra_frame();
    //    }
    //}

    if (m_selected_features.second.feature.has_value()) {
        const Measure::MeasurementResult measure = Measure::get_measurement(*m_selected_features.first.feature, *m_selected_features.second.feature);
        ImGui::Separator();
        if (measure.has_any_data()) {
            m_imgui->text(_u8L("Measure") + ":");
            if (ImGui::BeginTable("Measure", 2)) {
                if (measure.angle.has_value()) {
                    add_strings_row_to_table(_u8L("Angle") + _u8L(" (°)"), ImGuiWrapper::COL_ORANGE_LIGHT, format_double(Geometry::rad2deg(*measure.angle)),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                }
                if (measure.distance_infinite.has_value()) {
                    double distance = ObjectManipulation::mm_to_in * *measure.distance_infinite;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    add_strings_row_to_table(_u8L("Distance Infinite") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                }
                if (measure.distance_strict.has_value()) {
                    double distance = ObjectManipulation::mm_to_in * *measure.distance_strict;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    add_strings_row_to_table(_u8L("Distance Strict") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_double(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
                }
                if (measure.distance_xyz.has_value()) {
                    Vec3d distance = ObjectManipulation::mm_to_in * *measure.distance_xyz;
                    if (use_inches)
                        distance = ObjectManipulation::mm_to_in * distance;
                    add_strings_row_to_table(_u8L("Distance XYZ") + units, ImGuiWrapper::COL_ORANGE_LIGHT, format_vec3(distance),
                        ImGui::GetStyleColorVec4(ImGuiCol_Text));
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
}

} // namespace GUI
} // namespace Slic3r

#endif // ENABLE_MEASURE_GIZMO
