// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Model.hpp"
#include "libslic3r/Measure.hpp"

#include <numeric>

#include <GL/glew.h>

namespace Slic3r {
namespace GUI {

static const Slic3r::ColorRGBA DEFAULT_PLANE_COLOR       = { 0.9f, 0.9f, 0.9f, 0.9f };
static const Slic3r::ColorRGBA DEFAULT_HOVER_PLANE_COLOR = { 0.9f, 0.2f, 0.2f, 1.f };

GLGizmoMeasure::GLGizmoMeasure(GLCanvas3D& parent, const std::string& icon_filename, unsigned int sprite_id)
    : GLGizmoBase(parent, icon_filename, sprite_id)
{
    m_vbo_sphere.init_from(its_make_sphere(1., M_PI/32.));
    m_vbo_cylinder.init_from(its_make_cylinder(1., 1.));
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
    const ModelObject *model_object = nullptr;
    if (selection.is_single_full_instance() ||
        selection.is_from_single_object() ) {        
        model_object = selection.get_model()->objects[selection.get_object_idx()];
    }    
    set_flattening_data(model_object);
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
    // This is assumed in GLCanvas3D::do_rotate, do not change this
    // without updating that function too.
    return m_parent.get_selection().is_single_full_instance();
}



void GLGizmoMeasure::on_render()
{
    const Selection& selection = m_parent.get_selection();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;
    
    shader->start_using();

    glsafe(::glClear(GL_DEPTH_BUFFER_BIT));

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glEnable(GL_BLEND));
    glsafe(::glLineWidth(2.f));

    if (selection.is_single_full_instance()) {
        const Transform3d& m = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        
        
        update_if_needed();
            

        m_imgui->begin(std::string("DEBUG"));
        
        m_imgui->checkbox(wxString("Show all features"), m_show_all);
        m_imgui->checkbox(wxString("Show all planes"), m_show_planes);

        Vec3f pos;
        Vec3f normal;
        size_t facet_idx;
        m_c->raycaster()->raycasters().front()->unproject_on_mesh(Vec2d(m_mouse_pos_x, m_mouse_pos_y), m, camera, pos, normal, nullptr, &facet_idx);
        ImGui::Separator();
        m_imgui->text(std::string("face_idx: ") + std::to_string(facet_idx));
        m_imgui->text(std::string("pos_x: ") + std::to_string(pos.x()));
        m_imgui->text(std::string("pos_y: ") + std::to_string(pos.y()));
        m_imgui->text(std::string("pos_z: ") + std::to_string(pos.z()));



        std::vector<Measure::SurfaceFeature> features;
         if (m_show_all) {
            features = m_measuring->get_all_features(); // EXPENSIVE - debugging only.
            features.erase(std::remove_if(features.begin(), features.end(),
                           [](const Measure::SurfaceFeature& f) {
                            return f.get_type() == Measure::SurfaceFeatureType::Plane;
                            }), features.end());
         } else {
            std::optional<Measure::SurfaceFeature> feat = m_measuring->get_feature(facet_idx, pos.cast<double>());
            if (feat)
                features.emplace_back(*feat);
         }

            
            
        for (const Measure::SurfaceFeature& feature : features) {

            if (feature.get_type() == Measure::SurfaceFeatureType::Point) {
                Transform3d view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(feature.get_point()));
                view_feature_matrix.scale(0.5);
                shader->set_uniform("view_model_matrix", view_feature_matrix);
                m_vbo_sphere.set_color(ColorRGBA(0.8f, 0.2f, 0.2f, 1.f));
                m_vbo_sphere.render();
            }
            else if (feature.get_type() == Measure::SurfaceFeatureType::Circle) {
                const auto& [c, radius, n] = feature.get_circle();
                Transform3d view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(c));
                view_feature_matrix.scale(0.5);
                shader->set_uniform("view_model_matrix", view_feature_matrix);
                m_vbo_sphere.set_color(ColorRGBA(0.8f, 0.2f, 0.2f, 1.f));
                m_vbo_sphere.render();

                // Now draw the circle itself - let's take a funny shortcut:
                Vec3d rad = n.cross(Vec3d::UnitX());
                if (rad.squaredNorm() < 0.1)
                    rad = n.cross(Vec3d::UnitY());
                rad *= radius * rad.norm();
                const int N = 20;
                for (int i=0; i<N; ++i) {
                    rad = Eigen::AngleAxisd(6.28/N, n) * rad;
                    view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(c));
                    view_feature_matrix = view_feature_matrix * Transform3d(Eigen::Translation3d(rad));
                    view_feature_matrix.scale(N/100.);
                    shader->set_uniform("view_model_matrix", view_feature_matrix);
                    m_vbo_sphere.render();
                }
            }
            else if (feature.get_type() == Measure::SurfaceFeatureType::Edge) {
                const auto& [start, end] = feature.get_edge();
                Transform3d view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(start));
                auto q  = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), end - start);
                view_feature_matrix *= q;
                view_feature_matrix.scale(Vec3d(0.075, 0.075, (end - start).norm()));
                shader->set_uniform("view_model_matrix", view_feature_matrix);
                m_vbo_cylinder.set_color(ColorRGBA(0.8f, 0.2f, 0.2f, 1.f));
                m_vbo_cylinder.render();
                if (feature.get_extra_point()) {
                    Vec3d pin = *feature.get_extra_point();
                    view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(pin));
                    view_feature_matrix.scale(0.5);
                    shader->set_uniform("view_model_matrix", view_feature_matrix);
                    m_vbo_sphere.set_color(ColorRGBA(0.8f, 0.2f, 0.2f, 1.f));
                    m_vbo_sphere.render();
                }
            }
            else if (feature.get_type() == Measure::SurfaceFeatureType::Plane) {
                const auto& [idx, normal, pt] = feature.get_plane();
                assert(idx < m_plane_models.size());
                m_plane_models[idx]->render();
            }   
        }
        shader->set_uniform("view_model_matrix", view_model_matrix);
        if (m_show_planes)
            for (const auto& glmodel : m_plane_models)
                glmodel->render();
        
        m_imgui->end();
    }

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}





#if ! ENABLE_LEGACY_OPENGL_REMOVAL
    #error NOT IMPLEMENTED
#endif

void GLGizmoMeasure::set_flattening_data(const ModelObject* model_object)
{
    if (model_object != m_old_model_object)
        update_if_needed();
}


void GLGizmoMeasure::update_if_needed()
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    if (m_state != On || ! mo || mo->instances.empty())
        return;

    if (! m_measuring || mo != m_old_model_object
     || mo->volumes.size() != m_volumes_matrices.size())
        goto UPDATE;

    // We want to recalculate when the scale changes - some planes could (dis)appear.
    if (! mo->instances.front()->get_scaling_factor().isApprox(m_first_instance_scale)
     || ! mo->instances.front()->get_mirror().isApprox(m_first_instance_mirror))
        goto UPDATE;

    for (unsigned int i=0; i < mo->volumes.size(); ++i)
        if (! mo->volumes[i]->get_matrix().isApprox(m_volumes_matrices[i])
         || mo->volumes[i]->type() != m_volumes_types[i])
            goto UPDATE;

    return;

UPDATE:
    const indexed_triangle_set& its = mo->volumes.front()->mesh().its;
    m_measuring.reset(new Measure::Measuring(its));
    m_plane_models.clear();
    const std::vector<std::vector<int>> planes_triangles = m_measuring->get_planes_triangle_indices();
    for (const std::vector<int>& triangle_indices : planes_triangles) {
        m_plane_models.emplace_back(std::unique_ptr<GLModel>(new GLModel()));
        GUI::GLModel::Geometry init_data;
        init_data.format = { GUI::GLModel::Geometry::EPrimitiveType::Triangles, GUI::GLModel::Geometry::EVertexLayout::P3 };
        init_data.color = ColorRGBA(0.9f, 0.9f, 0.9f, 0.5f);
        int i = 0;
        for (int idx : triangle_indices) {  
            init_data.add_vertex(its.vertices[its.indices[idx][0]]);
            init_data.add_vertex(its.vertices[its.indices[idx][1]]);
            init_data.add_vertex(its.vertices[its.indices[idx][2]]);
            init_data.add_triangle(i, i+1, i+2);
            i+=3;
        }
        m_plane_models.back()->init_from(std::move(init_data));
    }

    // Let's save what we calculated it from:
    m_volumes_matrices.clear();
    m_volumes_types.clear();
    for (const ModelVolume* vol : mo->volumes) {
        m_volumes_matrices.push_back(vol->get_matrix());
        m_volumes_types.push_back(vol->type());
    }
    m_first_instance_scale = mo->instances.front()->get_scaling_factor();
    m_first_instance_mirror = mo->instances.front()->get_mirror();
    m_old_model_object = mo;
}

} // namespace GUI
} // namespace Slic3r
