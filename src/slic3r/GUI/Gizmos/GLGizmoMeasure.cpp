// Include GLGizmoBase.hpp before I18N.hpp as it includes some libigl code, which overrides our localization "L" macro.
#include "GLGizmoMeasure.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"

#include "slic3r/GUI/Gizmos/GLGizmosCommon.hpp"

#include "libslic3r/Geometry/ConvexHull.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/SurfaceMesh.hpp"
#include "libslic3r/Geometry/Circle.hpp"

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
    if (mouse_event.Moving()) {
        // only for sure 
        m_mouse_left_down = false;
        return false;
    }
    if (mouse_event.LeftDown()) {
        if (m_hover_id != -1) {
            m_mouse_left_down = true;
            Selection &selection = m_parent.get_selection();
            if (selection.is_single_full_instance()) {
                // Rotate the object so the normal points downward:
                selection.flattening_rotate(m_planes[m_hover_id].normal);
                m_parent.do_rotate(L("Gizmo-Place on Face"));
            }
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
    return CommonGizmosDataID::SelectionInfo;
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
        if (this->is_plane_update_necessary())
            update_planes();

        m_imgui->begin(std::string("DEBUG"));
        if (m_imgui->button("<-"))
            --m_currently_shown_plane;
        ImGui::SameLine();
        if (m_imgui->button("->"))
            ++m_currently_shown_plane;
        m_currently_shown_plane = std::clamp(m_currently_shown_plane, 0, std::max(0, int(m_planes.size())-1));
        m_imgui->text(std::to_string(m_currently_shown_plane));
        m_imgui->checkbox(wxString("Show all"), m_show_all_planes);
        m_imgui->checkbox(wxString("Show points"), m_show_points);
        m_imgui->checkbox(wxString("Show edges"), m_show_edges);
        m_imgui->checkbox(wxString("Show circles"), m_show_circles);
        m_imgui->end();


        int i = m_show_all_planes ? 0 : m_currently_shown_plane;
        for (; i < (int)m_planes.size(); ++i) {
            // Render all the borders.
            for (int j=0; j<(int)m_planes[i].vbos.size(); ++j) {
                m_planes[i].vbos[j].set_color(j == m_hover_id ? DEFAULT_HOVER_PLANE_COLOR : DEFAULT_PLANE_COLOR);
                m_planes[i].vbos[j].render();
            }


            // Render features:
            for (const SurfaceFeature& feature : m_planes[i].surface_features) {
                Transform3d view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(feature.pos));
                if (m_show_edges && feature.type == SurfaceFeature::Line) {
                    auto q  = Eigen::Quaternion<double>::FromTwoVectors(Vec3d::UnitZ(), feature.endpoint - feature.pos);
                    view_feature_matrix *= q;
                    view_feature_matrix.scale(Vec3d(0.3, 0.3, (feature.endpoint - feature.pos).norm()));
                    shader->set_uniform("view_model_matrix", view_feature_matrix);
                    m_vbo_cylinder.set_color(ColorRGBA(0.7f, 0.7f, 0.f, 1.f));
                    m_vbo_cylinder.render();
                }

                if ((m_show_points && feature.type == SurfaceFeature::Line) || m_show_circles && feature.type == SurfaceFeature::Circle) {
                    view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(feature.pos));
                    view_feature_matrix.scale(0.5);
                    shader->set_uniform("view_model_matrix", view_feature_matrix);
                    m_vbo_sphere.set_color(feature.type == SurfaceFeature::Line
                                            ? ColorRGBA(1.f, 0.f, 0.f, 1.f)
                                            : ColorRGBA(0.f, 1.f, 0.f, 1.f));
                    m_vbo_sphere.render();

                    /*view_feature_matrix = view_model_matrix * Transform3d(Eigen::Translation3d(feature.endpoint));
                    view_feature_matrix.scale(0.5);
                    shader->set_uniform("view_model_matrix", view_feature_matrix);
                    m_vbo_sphere.set_color(feature.type == SurfaceFeature::Line
                                            ? ColorRGBA(1.f, 0.f, 0.f, 1.f)
                                            : ColorRGBA(1.f, 1.f, 0.f, 1.f));
                    m_vbo_sphere.render();*/
                }

                
                shader->set_uniform("view_model_matrix", view_model_matrix);
            }

            if (! m_show_all_planes)
                break;
        }
    }

    glsafe(::glEnable(GL_CULL_FACE));
    glsafe(::glDisable(GL_BLEND));

    shader->stop_using();
}





#if ! ENABLE_LEGACY_OPENGL_REMOVAL
    #error NOT IMPLEMENTED
#endif
#if ! ENABLE_GL_SHADERS_ATTRIBUTES
    #error NOT IMPLEMENTED
#endif





void GLGizmoMeasure::on_render_for_picking()
{
    const Selection& selection = m_parent.get_selection();

    GLShaderProgram* shader = wxGetApp().get_shader("flat");
    if (shader == nullptr)
        return;

    shader->start_using();

    glsafe(::glDisable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_BLEND));
    glsafe(::glLineWidth(2.f));

    if (selection.is_single_full_instance() && !wxGetKeyState(WXK_CONTROL)) {
        const Transform3d& m = selection.get_volume(*selection.get_volume_idxs().begin())->get_instance_transformation().get_matrix();
        const Camera& camera = wxGetApp().plater()->get_camera();
        const Transform3d view_model_matrix = camera.get_view_matrix() *
            Geometry::assemble_transform(selection.get_volume(*selection.get_volume_idxs().begin())->get_sla_shift_z() * Vec3d::UnitZ()) * m;

        shader->set_uniform("view_model_matrix", view_model_matrix);
        shader->set_uniform("projection_matrix", camera.get_projection_matrix());
        if (this->is_plane_update_necessary())
            update_planes();
        //for (int i = 0; i < (int)m_planes.size(); ++i) {
        int i = m_currently_shown_plane;
        if (i < int(m_planes.size())) {
            for (int j=0; j<(int)m_planes[i].vbos.size(); ++j) {
                m_planes[i].vbos[j].set_color(picking_color_component(j));
                m_planes[i].vbos[j].render();
            }
        }
    }

    glsafe(::glEnable(GL_CULL_FACE));

    shader->stop_using();
}



void GLGizmoMeasure::set_flattening_data(const ModelObject* model_object)
{
    if (model_object != m_old_model_object) {
        m_planes.clear();
        m_planes_valid = false;
    }
}



static std::pair<Vec3d, double> get_center_and_radius(const std::vector<Vec3d>& border, int start_idx, int end_idx, const Transform3d& trafo)
{
    Vec2ds pts;
    double z = 0.;
    for (int i=start_idx; i<=end_idx; ++i) {
        Vec3d pt_transformed = trafo * border[i];
        z = pt_transformed.z();
        pts.emplace_back(pt_transformed.x(), pt_transformed.y());        
    }

    auto circle = Geometry::circle_ransac(pts, 20); // FIXME: iterations?

    return std::make_pair(trafo.inverse() * Vec3d(circle.center.x(), circle.center.y(), z), circle.radius);
}



void GLGizmoMeasure::extract_features(GLGizmoMeasure::PlaneData& plane)
{
    plane.surface_features.clear();
    const Vec3d& normal = plane.normal;

    const double edge_threshold = 25. * (M_PI/180.);
    std::vector<double> angles;

    Eigen::Quaterniond q;
    q.setFromTwoVectors(plane.normal, Vec3d::UnitZ());
    Transform3d trafo = Transform3d::Identity();
    trafo.rotate(q);

    
    
    for (const std::vector<Vec3d>& border : plane.borders) {
        assert(border.size() > 1);
        assert(! border.front().isApprox(border.back()));
        int start_idx = -1;


        // First calculate angles at all the vertices.
        angles.clear();
        for (int i=0; i<int(border.size()); ++i) {
            const Vec3d& v2 = (i == 0 ? border[0] - border[border.size()-1]
                                    : border[i] - border[i-1]);
            const Vec3d& v1 = i == border.size()-1 ? border[0] - border.back()
                                                : border[i+1] - border[i];
            double angle = -atan2(normal.dot(v1.cross(v2)), v1.dot(v2));
            if (angle < -M_PI/2.)
                angle += M_PI;
            angles.push_back(angle);
        }
        assert(border.size() == angles.size());


        bool circle = false;
        std::vector<std::pair<size_t, size_t>> circles;
        for (int i=1; i<angles.size(); ++i) {
            if (angles[i] < edge_threshold && Slic3r::is_approx(angles[i], angles[i-1]) && i != angles.size()-1 ) {
                // circle
                if (! circle) {
                    circle = true;
                    start_idx = std::max(0, i-2);
                }
            } else {
                if (circle) {
                    circles.emplace_back(start_idx, i);
                    circle = false;
                }
            }
        }

        for (const auto& [start_idx, end_idx] : circles) {
            std::pair<Vec3d, double> center_and_radius = get_center_and_radius(border, start_idx, end_idx, trafo);
            plane.surface_features.emplace_back(SurfaceFeature{
                SurfaceFeature::Circle,
                // border[start_idx], border[end_idx],
                center_and_radius.first, center_and_radius.first, center_and_radius.second                
            });
        }


        std::cout << "==================== " << std::endl;
    }


     for (const SurfaceFeature& f : plane.surface_features) {
            std::cout << "- detected " << (f.type == SurfaceFeature::Line ? "Line" : "Circle") << std::endl;
            std::cout<< f.pos << std::endl << std::endl << f.endpoint << std::endl;
            std::cout << "----------------" << std::endl;
        }



}



void GLGizmoMeasure::update_planes()
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    TriangleMesh ch;
    for (const ModelVolume* vol : mo->volumes) {
        if (vol->type() != ModelVolumeType::MODEL_PART)
            continue;
        TriangleMesh vol_ch = vol->mesh();
        vol_ch.transform(vol->get_matrix());
        ch.merge(vol_ch);
    }
    m_planes.clear();
    

    // Now we'll go through all the facets and append Points of facets sharing the same normal.
    // This part is still performed in mesh coordinate system.
    const size_t             num_of_facets  = ch.facets_count();
    std::vector<size_t>      face_to_plane(num_of_facets, size_t(-1));
    const std::vector<Vec3f> face_normals   = its_face_normals(ch.its);
    const std::vector<Vec3i> face_neighbors = its_face_neighbors(ch.its);
    std::vector<int>         facet_queue(num_of_facets, 0);
    int                      facet_queue_cnt = 0;
    const stl_normal*        normal_ptr      = nullptr;
    size_t seed_facet_idx = 0;

    auto is_same_normal = [](const stl_normal& a, const stl_normal& b) -> bool {
        return (std::abs(a(0) - b(0)) < 0.001 && std::abs(a(1) - b(1)) < 0.001 && std::abs(a(2) - b(2)) < 0.001);
    };

    while (1) {
        // Find next unvisited triangle:
        for (; seed_facet_idx < num_of_facets; ++ seed_facet_idx)
            if (face_to_plane[seed_facet_idx] == size_t(-1)) {
                facet_queue[facet_queue_cnt ++] = seed_facet_idx;
                normal_ptr = &face_normals[seed_facet_idx];
                face_to_plane[seed_facet_idx] = m_planes.size();
                m_planes.emplace_back();                
                break;
            }
        if (seed_facet_idx == num_of_facets)
            break; // Everything was visited already

        while (facet_queue_cnt > 0) {
            int facet_idx = facet_queue[-- facet_queue_cnt];
            const stl_normal& this_normal = face_normals[facet_idx];
            if (is_same_normal(this_normal, *normal_ptr)) {
                const Vec3i& face = ch.its.indices[facet_idx];

                face_to_plane[facet_idx] = m_planes.size() - 1;
                m_planes.back().facets.emplace_back(facet_idx);
                for (int j = 0; j < 3; ++ j)
                    if (int neighbor_idx = face_neighbors[facet_idx][j]; neighbor_idx >= 0 && face_to_plane[neighbor_idx] == size_t(-1))
                        facet_queue[facet_queue_cnt ++] = neighbor_idx;
            }
        }

        m_planes.back().normal = normal_ptr->cast<double>();
        std::sort(m_planes.back().facets.begin(), m_planes.back().facets.end());
    }

    assert(std::none_of(face_to_plane.begin(), face_to_plane.end(), [](size_t val) { return val == size_t(-1); }));

    SurfaceMesh sm(ch.its);
    for (int plane_id=0; plane_id < int(m_planes.size()); ++plane_id) {
    //int plane_id = 5; {
        const auto& facets = m_planes[plane_id].facets;
        m_planes[plane_id].borders.clear();
        std::vector<std::array<bool, 3>> visited(facets.size(), {false, false, false});
        
        for (int face_id=0; face_id<int(facets.size()); ++face_id) {
            assert(face_to_plane[facets[face_id]] == plane_id);
            for (int edge_id=0; edge_id<3; ++edge_id) {
                if (visited[face_id][edge_id] || face_to_plane[face_neighbors[facets[face_id]][edge_id]] == plane_id) {
                    visited[face_id][edge_id] = true;
                    continue;
                }

                Halfedge_index he = sm.halfedge(Face_index(facets[face_id]));
                while (he.side() != edge_id)
                    he = sm.next(he);
            
                // he is the first halfedge on the border. Now walk around and append the points.
                //const Halfedge_index he_orig = he;
                m_planes[plane_id].borders.emplace_back();
                std::vector<Vec3d>& last_border = m_planes[plane_id].borders.back();
                last_border.emplace_back(sm.point(sm.source(he)).cast<double>());
                //Vertex_index target = sm.target(he);
                const Halfedge_index he_start = he;
                
                Face_index fi = he.face();
                auto face_it = std::lower_bound(facets.begin(), facets.end(), int(fi));
                assert(face_it != facets.end());
                assert(*face_it == int(fi));
                visited[face_it - facets.begin()][he.side()] = true;

                do {
                    const Halfedge_index he_orig = he;
                    he = sm.next_around_target(he);
                    while ( face_to_plane[sm.face(he)] == plane_id && he != he_orig)
                        he = sm.next_around_target(he);
                    he = sm.opposite(he);
                    
                    Face_index fi = he.face();
                    auto face_it = std::lower_bound(facets.begin(), facets.end(), int(fi));
                    assert(face_it != facets.end());
                    assert(*face_it == int(fi));
                    if (visited[face_it - facets.begin()][he.side()] && he != he_start) {
                        last_border.resize(1);
                        break;
                    }
                    visited[face_it - facets.begin()][he.side()] = true;

                    last_border.emplace_back(sm.point(sm.source(he)).cast<double>());
                } while (he != he_start);

                if (last_border.size() == 1)
                    m_planes[plane_id].borders.pop_back();
            }
        }
    }


    // DEBUGGING:
    m_planes.erase(std::remove_if(m_planes.begin(), m_planes.end(), [](const PlaneData& p) { return p.borders.empty(); }), m_planes.end());





 
    

    // Planes are finished - let's save what we calculated it from:
    m_volumes_matrices.clear();
    m_volumes_types.clear();
    for (const ModelVolume* vol : mo->volumes) {
        m_volumes_matrices.push_back(vol->get_matrix());
        m_volumes_types.push_back(vol->type());
    }
    m_first_instance_scale = mo->instances.front()->get_scaling_factor();
    m_first_instance_mirror = mo->instances.front()->get_mirror();
    m_old_model_object = mo;

    // And finally create respective VBOs. The polygon is convex with
    // the vertices in order, so triangulation is trivial.
    for (PlaneData& plane : m_planes) {
        for (std::vector<Vec3d>& vertices : plane.borders) {
            GLModel::Geometry init_data;
            init_data.format = { GLModel::Geometry::EPrimitiveType::LineStrip, GLModel::Geometry::EVertexLayout::P3N3 };
            init_data.reserve_vertices(vertices.size());
            init_data.reserve_indices(vertices.size());
            // vertices + indices
            for (size_t i = 0; i < vertices.size(); ++i) {
                init_data.add_vertex((Vec3f)vertices[i].cast<float>(), (Vec3f)plane.normal.cast<float>());
                init_data.add_index((unsigned int)i);
            }
            plane.vbos.emplace_back();
            plane.vbos.back().init_from(std::move(init_data));
            vertices.pop_back(); // first and last are the same
        }

        static int n=0;
        std::cout << "==================== " << std::endl;
        std::cout << "==================== " << std::endl;
        std::cout << "==================== " << std::endl;
        std::cout << "Plane num. " << n++ << std::endl;
        extract_features(plane);


        // FIXME: vertices should really be local, they need not
        // persist now when we use VBOs
        plane.borders.clear();
        plane.borders.shrink_to_fit();
    }

    m_planes_valid = true;
}



bool GLGizmoMeasure::is_plane_update_necessary() const
{
    const ModelObject* mo = m_c->selection_info()->model_object();
    if (m_state != On || ! mo || mo->instances.empty())
        return false;

    if (! m_planes_valid || mo != m_old_model_object
     || mo->volumes.size() != m_volumes_matrices.size())
        return true;

    // We want to recalculate when the scale changes - some planes could (dis)appear.
    if (! mo->instances.front()->get_scaling_factor().isApprox(m_first_instance_scale)
     || ! mo->instances.front()->get_mirror().isApprox(m_first_instance_mirror))
        return true;

    for (unsigned int i=0; i < mo->volumes.size(); ++i)
        if (! mo->volumes[i]->get_matrix().isApprox(m_volumes_matrices[i])
         || mo->volumes[i]->type() != m_volumes_types[i])
            return true;

    return false;
}

} // namespace GUI
} // namespace Slic3r
