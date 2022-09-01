// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/PresetBundle.hpp"

#include <numeric>

#include <GL/glew.h>

#if ENABLE_MEASURE_GIZMO

std::string surface_feature_type_as_string(Slic3r::Measure::SurfaceFeatureType type)
{
    switch (type)
    {
    default:
    case Slic3r::Measure::SurfaceFeatureType::Undef:  { return L("Undefined"); }
    case Slic3r::Measure::SurfaceFeatureType::Point:  { return L("Vertex"); }
    case Slic3r::Measure::SurfaceFeatureType::Edge:   { return L("Edge"); }
    case Slic3r::Measure::SurfaceFeatureType::Circle: { return L("Circle"); }
    case Slic3r::Measure::SurfaceFeatureType::Plane:  { return L("Plane"); }
    }
}

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA BASIC_HOVER_COLOR    = { 0.8f, 0.2f, 0.2f, 1.0f };
static const Slic3r::ColorRGBA EXTENDED_HOVER_COLOR = { 0.8f, 0.8f, 0.2f, 1.0f };
static const Slic3r::ColorRGBA LOCK_COLOR           = { 0.5f, 0.5f, 0.5f, 1.0f };

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

GLGizmoMeasure::GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{
    GLModel::Geometry sphere_geometry = smooth_sphere(16, 7.5f);
    m_sphere.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(sphere_geometry.get_as_indexed_triangle_set())));
    m_sphere.model.init_from(std::move(sphere_geometry));

    GLModel::Geometry cylinder_geometry = smooth_cylinder(16, 5.0f, 1.0f);
    m_cylinder.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(cylinder_geometry.get_as_indexed_triangle_set())));
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
}

bool GLGizmoMeasure::gizmo_event(SLAGizmoEventType action, const Vec2d& mouse_position, bool shift_down, bool alt_down, bool control_down)
{
    if (action == SLAGizmoEventType::CtrlDown) {
        if (m_ctrl_kar_filter.is_first()) {
            if (m_curr_feature.has_value())
                m_mode = EMode::ExtendedSelection;
        }

        m_ctrl_kar_filter.increase_count();
    }
    else if (action == SLAGizmoEventType::CtrlUp) {
        m_ctrl_kar_filter.reset_count();
        m_mode = EMode::BasicSelection;
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
        m_curr_ex_feature_position.reset();
    }
    else
        m_mode = EMode::BasicSelection;
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
        size_t facet_idx;
        m_c->raycaster()->raycasters().front()->unproject_on_mesh(m_mouse_pos, model_matrix, camera, position_on_model, normal_on_model, nullptr, &facet_idx);

        const bool is_hovering_on_extended_selection = m_mode == EMode::ExtendedSelection && m_hover_id != -1;

        std::optional<Measure::SurfaceFeature> curr_feature;
        if (m_mode == EMode::BasicSelection)
            curr_feature = m_measuring->get_feature(facet_idx, position_on_model.cast<double>());

        if (m_mode == EMode::BasicSelection) {
            m_curr_ex_feature_position.reset();
            if (m_curr_feature != curr_feature) {
                GLGizmoMeasure::on_unregister_raycasters_for_picking();
                m_curr_feature = curr_feature;
                if (!m_curr_feature.has_value())
                    return;

                switch (m_curr_feature->get_type()) {
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
                    const auto& [center, radius, normal] = m_curr_feature->get_circle();

                    if (m_last_inv_zoom != inv_zoom) {
                        m_last_inv_zoom = inv_zoom;
                        m_circle.reset();
                        GLModel::Geometry circle_geometry = smooth_torus(64, 16, float(radius), 5.0f * inv_zoom);
                        m_circle.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(circle_geometry.get_as_indexed_triangle_set())));
                        m_circle.model.init_from(std::move(circle_geometry));
                    }

                    m_raycasters.insert({ CIRCLE_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, CIRCLE_ID, *m_circle.mesh_raycaster) });
                    m_raycasters.insert({ POINT_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, POINT_ID, *m_sphere.mesh_raycaster) });
                    break;
                }
                case Measure::SurfaceFeatureType::Plane:
                {
                    const auto& [idx, normal, point] = m_curr_feature->get_plane();
                    if (m_last_plane_idx != idx) {
                        m_last_plane_idx = idx;
                        const std::vector<std::vector<int>> planes_triangles = m_measuring->get_planes_triangle_indices();
                        const std::vector<int>& triangle_indices = planes_triangles[idx];

                        const indexed_triangle_set its = (m_old_model_volume != nullptr) ? m_old_model_volume->mesh().its : m_old_model_object->volumes.front()->mesh().its;
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

                        m_plane.reset();
                        m_plane.mesh_raycaster = std::make_unique<MeshRaycaster>(std::make_shared<const TriangleMesh>(std::move(init_data.get_as_indexed_triangle_set())));
                        m_plane.model.init_from(std::move(init_data));
                    }

                    m_raycasters.insert({ PLANE_ID, m_parent.add_raycaster_for_picking(SceneRaycaster::EType::Gizmo, PLANE_ID, *m_plane.mesh_raycaster) });
                    break;
                }
                }
            }
        }
        else if (is_hovering_on_extended_selection) {
            switch (m_curr_feature->get_type())
            {
            case Measure::SurfaceFeatureType::Point:
            {
                m_curr_ex_feature_position = model_matrix * m_curr_feature->get_point();
                break;
            }
            case Measure::SurfaceFeatureType::Edge:
            {
                auto it = m_raycasters.find(EDGE_ID);
                if (it != m_raycasters.end() && it->second != nullptr) {
                    Vec3f p;
                    Vec3f n;
                    const Transform3d& trafo = it->second->get_transform();
                    it->second->get_raycaster()->unproject_on_mesh(m_mouse_pos, trafo, camera, p, n);
                    p = { 0.0f, 0.0f, p.z() };
                    m_curr_ex_feature_position = trafo * p.cast<double>();
                }
                break;
            }
            case Measure::SurfaceFeatureType::Plane:
            {
                m_curr_ex_feature_position = model_matrix * position_on_model.cast<double>();
                break;
            }
            case Measure::SurfaceFeatureType::Circle:
            {
                const auto& [center, radius, normal] = m_curr_feature->get_circle();
                if (m_hover_id == POINT_ID)
                    m_curr_ex_feature_position = model_matrix * center;
                else {
                    auto it = m_raycasters.find(CIRCLE_ID);
                    if (it != m_raycasters.end() && it->second != nullptr) {
                        Vec3f p;
                        Vec3f n;
                        const Transform3d& trafo = it->second->get_transform();
                        it->second->get_raycaster()->unproject_on_mesh(m_mouse_pos, trafo, camera, p, n);
                        float angle = std::atan2(p.y(), p.x());
                        if (angle < 0.0f)
                            angle += 2.0f * float(M_PI);
                        p = float(radius) * Vec3f(std::cos(angle), std::sin(angle), 0.0f);
                        m_curr_ex_feature_position = trafo * p.cast<double>();
                    }
                }
                break;
            }
            }
        }

        if (!m_curr_feature.has_value())
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

        ColorRGBA feature_color;
        switch (m_mode)
        {
        case EMode::BasicSelection:    { feature_color = BASIC_HOVER_COLOR; break; }
        case EMode::ExtendedSelection: { feature_color = LOCK_COLOR; break; }
        }

        auto set_matrix_uniforms = [shader, &view_matrix](const Transform3d& model_matrix) {
            const Transform3d view_model_matrix = view_matrix * model_matrix;
            shader->set_uniform("view_model_matrix", view_model_matrix);
            const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
            shader->set_uniform("view_normal_matrix", view_normal_matrix);
        };

        switch (m_curr_feature->get_type()) {
        case Measure::SurfaceFeatureType::Point:
        {
            const Vec3d& position = m_curr_feature->get_point();
            const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(position) * Geometry::scale_transform(inv_zoom);
            set_matrix_uniforms(feature_matrix);
            m_sphere.model.set_color(is_hovering_on_extended_selection ? EXTENDED_HOVER_COLOR : feature_color);
            m_sphere.model.render();
            auto it = m_raycasters.find(POINT_ID);
            if (it != m_raycasters.end() && it->second != nullptr)
                it->second->set_transform(feature_matrix);
            break;
        }
        case Measure::SurfaceFeatureType::Circle:
        {
            const auto& [center, radius, n] = m_curr_feature->get_circle();
            // render center
            const Transform3d center_matrix = model_matrix * Geometry::translation_transform(center) * Geometry::scale_transform(inv_zoom);
            set_matrix_uniforms(center_matrix);
            m_sphere.model.set_color((is_hovering_on_extended_selection && m_hover_id == POINT_ID) ? EXTENDED_HOVER_COLOR : feature_color);
            m_sphere.model.render();
            auto it = m_raycasters.find(POINT_ID);
            if (it != m_raycasters.end() && it->second != nullptr)
                it->second->set_transform(center_matrix);

            // render circle
            const Transform3d circle_matrix = model_matrix * Geometry::translation_transform(center);
            set_matrix_uniforms(circle_matrix);
            m_circle.model.set_color(feature_color);
            m_circle.model.render();
            it = m_raycasters.find(CIRCLE_ID);
            if (it != m_raycasters.end() && it->second != nullptr)
                it->second->set_transform(circle_matrix);
            break;
        }
        case Measure::SurfaceFeatureType::Edge:
        {
            const auto& [start, end] = m_curr_feature->get_edge();
            auto q = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), end - start);
            const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(start) * q *
                Geometry::scale_transform({ (double)inv_zoom, (double)inv_zoom, (end - start).norm() });
            set_matrix_uniforms(feature_matrix);
            m_cylinder.model.set_color(feature_color);
            m_cylinder.model.render();
            auto it = m_raycasters.find(EDGE_ID);
            if (it != m_raycasters.end() && it->second != nullptr)
                it->second->set_transform(feature_matrix);
            break;
        }
        case Measure::SurfaceFeatureType::Plane:
        {
            const auto& [idx, normal, pt] = m_curr_feature->get_plane();
            assert(idx < m_plane_models_cache.size());
            set_matrix_uniforms(model_matrix);
            m_plane_models_cache[idx].set_color(feature_color);
            m_plane_models_cache[idx].render();
            auto it = m_raycasters.find(PLANE_ID);
            if (it != m_raycasters.end() && it->second != nullptr)
                it->second->set_transform(model_matrix);
            break;
        }
        }

        if (is_hovering_on_extended_selection && m_curr_ex_feature_position.has_value()) {
            if (m_hover_id != POINT_ID) {
                const Transform3d matrix = Geometry::translation_transform(*m_curr_ex_feature_position) * Geometry::scale_transform(inv_zoom);
                set_matrix_uniforms(matrix);
                m_sphere.model.set_color(EXTENDED_HOVER_COLOR);
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
        for (const std::vector<int>& triangle_indices : planes_triangles) {
            m_plane_models_cache.emplace_back(GLModel());
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
        m_first_instance_scale = Vec3d::Ones();
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

void GLGizmoMeasure::on_render_input_window(float x, float y, float bottom_limit)
{
    static std::optional<Measure::SurfaceFeature> last_feature;
    static EMode last_mode = EMode::BasicSelection;

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

    auto add_row_to_table = [this](const wxString& label, const std::string& value) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        m_imgui->text_colored(ImGuiWrapper::COL_ORANGE_LIGHT, label);
        ImGui::TableSetColumnIndex(1);
        m_imgui->text(value);
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
        add_row_to_table(_u8L("Left mouse button"), (m_mode == EMode::BasicSelection) ? _u8L("Select feature") : _u8L("Select point"));
        if (m_mode == EMode::BasicSelection && m_hover_id != -1)
            add_row_to_table(CTRL_STR, _u8L("Enable point selection"));
        ImGui::EndTable();
    }

    if (m_curr_feature.has_value()) {
        const Transform3d volume_matrix = m_parent.get_selection().get_first_volume()->world_matrix();
        const Measure::SurfaceFeatureType feature_type = m_curr_feature->get_type();
        if (m_mode == EMode::BasicSelection) {
            if (feature_type != Measure::SurfaceFeatureType::Undef) {
                if (ImGui::CollapsingHeader(surface_feature_type_as_string(feature_type).c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
                    if (ImGui::BeginTable("Data", 2)) {
                        switch (feature_type)
                        {
                        case Measure::SurfaceFeatureType::Point:
                        {
                            const Vec3d position = volume_matrix * m_curr_feature->get_point();
                            add_row_to_table(_u8L("Position") + ":", format_vec3(position));
                            break;
                        }
                        case Measure::SurfaceFeatureType::Edge:
                        {
                            auto [from, to] = m_curr_feature->get_edge();
                            from = volume_matrix * from;
                            to = volume_matrix * to;
                            add_row_to_table(_u8L("From") + ":", format_vec3(from));
                            add_row_to_table(_u8L("To") + ":", format_vec3(to));
                            break;
                        }
                        case Measure::SurfaceFeatureType::Circle:
                        {
                            auto [center, radius, normal] = m_curr_feature->get_circle();
                            center = volume_matrix * center;
                            normal = volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                            add_row_to_table(_u8L("Center") + ":", format_vec3(center));
                            add_row_to_table(_u8L("Radius") + ":", format_double(radius));
                            add_row_to_table(_u8L("Normal") + ":", format_vec3(normal));
                            break;
                        }
                        case Measure::SurfaceFeatureType::Plane:
                        {
                            auto [idx, normal, origin] = m_curr_feature->get_plane();
                            origin = volume_matrix * origin;
                            normal = volume_matrix.matrix().block(0, 0, 3, 3).inverse().transpose() * normal;
                            add_row_to_table(_u8L("Origin") + ":", format_vec3(origin));
                            add_row_to_table(_u8L("Normal") + ":", format_vec3(normal));
                            break;
                        }
                        }
                        ImGui::EndTable();
                    }
                }
            }
        }
        else if (m_mode == EMode::ExtendedSelection) {
            if (m_hover_id != -1 && m_curr_ex_feature_position.has_value()) {
                std::string header;
                switch (feature_type) {
                case Measure::SurfaceFeatureType::Point:  { header = _u8L("Vertex"); break; }
                case Measure::SurfaceFeatureType::Edge:   { header = _u8L("Point on edge"); break; }
                case Measure::SurfaceFeatureType::Circle: { header = (m_hover_id == POINT_ID) ? _u8L("Center of circle") : _u8L("Point on circle"); break; }
                case Measure::SurfaceFeatureType::Plane:  { header = _u8L("Point on plane"); break; }
                default: { assert(false); break; }
                }

                if (ImGui::CollapsingHeader(header.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
                    if (ImGui::BeginTable("Data", 2)) {
                        add_row_to_table(_u8L("Position") + ":", format_vec3(*m_curr_ex_feature_position));
                        ImGui::EndTable();
                    }
                }
            }
        }
    }

    if (last_feature != m_curr_feature || last_mode != m_mode) {
        // the dialog may have changed its size, ask for an extra frame to render it properly
        last_feature = m_curr_feature;
        last_mode = m_mode;
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
