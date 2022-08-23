// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/Measure.hpp"
#include "libslic3r/PresetBundle.hpp"

#include <numeric>

#include <GL/glew.h>

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA HOVER_COLOR = { 0.8f, 0.2f, 0.2f, 1.0f };

GLGizmoMeasure::GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{
    m_vbo_sphere.init_from(smooth_sphere(16, 7.5f));
    m_vbo_cylinder.init_from(smooth_cylinder(16, 5.0f, 1.0f));
}

bool GLGizmoMeasure::on_mouse(const wxMouseEvent &mouse_event)
{
    m_mouse_pos_x = mouse_event.GetX();
    m_mouse_pos_y = mouse_event.GetY();


    if (mouse_event.Moving()) {
        // only for sure 
        m_mouse_left_down = false;
        return false;
    }
    if (mouse_event.LeftDown()) {
        if (m_hover_id != -1) {
            m_mouse_left_down = true;
            
            return true;
        }

        // fix: prevent restart gizmo when reselect object
        // take responsibility for left up
        if (m_parent.get_first_hover_volume_idx() >= 0) m_mouse_left_down = true;
        
    } else if (mouse_event.LeftUp()) {
        if (m_mouse_left_down) {
            // responsible for mouse left up after selecting plane
            m_mouse_left_down = false;
            return true;
        }
    } else if (mouse_event.Leaving()) {
        m_mouse_left_down = false;
    }
    return false;
}



void GLGizmoMeasure::data_changed()
{
    const Selection &  selection    = m_parent.get_selection();
    const ModelObject* model_object = nullptr;
    const ModelVolume* model_volume = nullptr;
    if (selection.is_single_full_instance() ||
        selection.is_from_single_object() ) {        
        model_object = selection.get_model()->objects[selection.get_object_idx()];
        model_volume = model_object->volumes[selection.get_first_volume()->volume_idx()];
    }
    if (model_object != m_old_model_object || model_volume != m_old_model_volume)
        update_if_needed();
}



bool GLGizmoMeasure::on_init()
{
    // FIXME m_shortcut_key = WXK_CONTROL_F;
    return true;
}



void GLGizmoMeasure::on_set_state()
{
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

    GLShaderProgram* shader = wxGetApp().get_shader("gouraud_light");
    if (shader == nullptr)
        return;
    
    shader->start_using();
    shader->set_uniform("emission_factor", 0.25f);

    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glEnable(GL_BLEND));

    if ((wxGetApp().preset_bundle->printers.get_edited_preset().printer_technology() == ptSLA && selection.is_single_full_instance()) ||
        (selection.is_single_volume() || selection.is_single_volume_instance())) {
        const Transform3d& model_matrix = selection.get_first_volume()->world_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d& view_matrix = camera.get_view_matrix();
        const float inv_zoom = camera.get_inv_zoom();

        shader->set_uniform("projection_matrix", camera.get_projection_matrix());

        update_if_needed();

        Vec3f pos;
        Vec3f normal;
        size_t facet_idx;
        m_c->raycaster()->raycasters().front()->unproject_on_mesh(Vec2d(m_mouse_pos_x, m_mouse_pos_y), model_matrix, camera, pos, normal, nullptr, &facet_idx);

#if ENABLE_DEBUG_DIALOG
        m_imgui->begin(std::string("DEBUG"));
        m_imgui->checkbox(wxString("Show all features"), m_show_all);
        m_imgui->checkbox(wxString("Show all planes"), m_show_planes);
        ImGui::Separator();
        m_imgui->text(std::string("face_idx: ") + std::to_string(facet_idx));
        m_imgui->text(std::string("pos_x: ") + std::to_string(pos.x()));
        m_imgui->text(std::string("pos_y: ") + std::to_string(pos.y()));
        m_imgui->text(std::string("pos_z: ") + std::to_string(pos.z()));
        m_imgui->end();
#endif // ENABLE_DEBUG_DIALOG

        std::vector<Measure::SurfaceFeature> features;
#if ENABLE_DEBUG_DIALOG
        if (m_show_all) {
            features = m_measuring->get_all_features(); // EXPENSIVE - debugging only.
            features.erase(std::remove_if(features.begin(), features.end(),
                           [](const Measure::SurfaceFeature& f) {
                            return f.get_type() == Measure::SurfaceFeatureType::Plane;
                            }), features.end());
        }
        else {
#endif // ENABLE_DEBUG_DIALOG
            std::optional<Measure::SurfaceFeature> feat = m_measuring->get_feature(facet_idx, pos.cast<double>());
            if (feat.has_value())
                features.emplace_back(*feat);
#if ENABLE_DEBUG_DIALOG
        }
#endif // ENABLE_DEBUG_DIALOG
            
        for (const Measure::SurfaceFeature& feature : features) {
            switch (feature.get_type()) {
            case Measure::SurfaceFeatureType::Point:
            {
                const Vec3d& position = feature.get_point();
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(position) * Geometry::scale_transform(inv_zoom);
                const Transform3d view_model_matrix = view_matrix * feature_matrix;
                shader->set_uniform("view_model_matrix", view_model_matrix);
                const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * feature_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                shader->set_uniform("view_normal_matrix", view_normal_matrix);
                m_vbo_sphere.set_color(HOVER_COLOR);
                m_vbo_sphere.render();
                break;
            }
            case Measure::SurfaceFeatureType::Circle:
            {
                const auto& [center, radius, n] = feature.get_circle();
                // render center
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(center) * Geometry::scale_transform(inv_zoom);
                const Transform3d view_model_matrix = view_matrix * feature_matrix;
                shader->set_uniform("view_model_matrix", view_model_matrix);
                const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * feature_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                shader->set_uniform("view_normal_matrix", view_normal_matrix);
                m_vbo_sphere.set_color(HOVER_COLOR);
                m_vbo_sphere.render();

                // Now draw the circle itself - let's take a funny shortcut:
                Vec3d rad = n.cross(Vec3d::UnitX());
                if (rad.squaredNorm() < 0.1)
                    rad = n.cross(Vec3d::UnitY());
                rad *= radius * rad.norm();
                const int N = 32;
                m_vbo_sphere.set_color(HOVER_COLOR);
                for (int i = 0; i < N; ++i) {
                    rad = Eigen::AngleAxisd(2.0 * M_PI / N, n) * rad;
                    const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(center + rad) * Geometry::scale_transform(N / 100.0 * inv_zoom);
                    const Transform3d view_model_matrix = view_matrix * feature_matrix;
                    shader->set_uniform("view_model_matrix", view_model_matrix);
                    const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * feature_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                    shader->set_uniform("view_normal_matrix", view_normal_matrix);
                    m_vbo_sphere.render();
                }
                break;
            }
            case Measure::SurfaceFeatureType::Edge:
            {
                const auto& [start, end] = feature.get_edge();
                auto q = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), end - start);
                const Transform3d feature_matrix = model_matrix * Geometry::translation_transform(start) * q *
                    Geometry::scale_transform({ (double)inv_zoom, (double)inv_zoom, (end - start).norm() });
                const Transform3d view_model_matrix = view_matrix * feature_matrix;
                shader->set_uniform("view_model_matrix", view_model_matrix);
                const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * feature_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                shader->set_uniform("view_normal_matrix", view_normal_matrix);
                m_vbo_cylinder.set_color(HOVER_COLOR);
                m_vbo_cylinder.render();

/*
                std::optional<Vec3d> extra_point = feature.get_extra_point();
                if (extra_point.has_value()) {
                    const Vec3d pin = *extra_point;
                    const Transform3d feature_matrix = Geometry::translation_transform(pin + selection.get_first_volume()->get_sla_shift_z() * Vec3d::UnitZ()) *
                        Geometry::scale_transform(inv_zoom) * model_matrix;
                    const Transform3d view_model_matrix = view_matrix * feature_matrix;
                    shader->set_uniform("view_model_matrix", view_model_matrix);
                    const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * feature_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                    shader->set_uniform("view_normal_matrix", view_normal_matrix);
                    m_vbo_sphere.set_color(HOVER_COLOR);
                    m_vbo_sphere.render();
                }
*/
                break;
            }
            case Measure::SurfaceFeatureType::Plane:
            {
                const auto& [idx, normal, pt] = feature.get_plane();
                assert(idx < m_plane_models.size());
                const Transform3d view_model_matrix = view_matrix * model_matrix;
                shader->set_uniform("view_model_matrix", view_model_matrix);
                const Matrix3d view_normal_matrix = view_matrix.matrix().block(0, 0, 3, 3) * model_matrix.matrix().block(0, 0, 3, 3).inverse().transpose();
                shader->set_uniform("view_normal_matrix", view_normal_matrix);
                m_plane_models[idx]->set_color(HOVER_COLOR);
                m_plane_models[idx]->render();
                break;
            }
            }
        }
#if ENABLE_DEBUG_DIALOG
        if (m_show_planes)
            for (const auto& glmodel : m_plane_models) {
                glmodel->set_color(HOVER_COLOR);
                glmodel->render();
            }
#endif // ENABLE_DEBUG_DIALOG
    }

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}





#if ! ENABLE_LEGACY_OPENGL_REMOVAL
    #error NOT IMPLEMENTED
#endif


void GLGizmoMeasure::update_if_needed()
{
    auto do_update = [this](const ModelObject* object, const ModelVolume* volume) {
        const indexed_triangle_set& its = (volume != nullptr) ? volume->mesh().its : object->volumes.front()->mesh().its;
        m_measuring.reset(new Measure::Measuring(its));
        m_plane_models.clear();
        const std::vector<std::vector<int>> planes_triangles = m_measuring->get_planes_triangle_indices();
        for (const std::vector<int>& triangle_indices : planes_triangles) {
            m_plane_models.emplace_back(std::unique_ptr<GLModel>(new GLModel()));
            GLModel::Geometry init_data;
            init_data.format = { GUI::GLModel::Geometry::EPrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3 };
            init_data.color = ColorRGBA(0.9f, 0.9f, 0.9f, 0.5f);
            int i = 0;
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
            m_plane_models.back()->init_from(std::move(init_data));
        }

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
            m_first_instance_scale = object->instances.front()->get_scaling_factor();
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

} // namespace GUI
} // namespace Slic3r
