#include "Measure.hpp"

#include "libslic3r/Geometry/Circle.hpp"
#include "libslic3r/SurfaceMesh.hpp"



namespace Slic3r {
namespace Measure {


constexpr double feature_hover_limit = 0.5; // how close to a feature the mouse must be to highlight it
constexpr double edge_endpoint_limit = 0.5; // how close to an edge endpoint the mouse ...



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




class MeasuringImpl {
public:
    explicit MeasuringImpl(const indexed_triangle_set& its);
    struct PlaneData {
        std::vector<int> facets;
        std::vector<std::vector<Vec3d>> borders; // FIXME: should be in fact local in update_planes()
        std::vector<SurfaceFeature> surface_features;
        Vec3d normal;
        float area;
    };

    std::vector<SurfaceFeature> get_all_features() const;
    std::optional<SurfaceFeature> get_feature(size_t face_idx, const Vec3d& point) const;
    std::vector<std::vector<int>> get_planes_triangle_indices() const;

private:
    void update_planes();
    void extract_features();
    
    std::vector<PlaneData> m_planes;
    std::vector<size_t>    m_face_to_plane;
    const indexed_triangle_set& m_its;
};






MeasuringImpl::MeasuringImpl(const indexed_triangle_set& its)
: m_its{its}
{
    update_planes();
    extract_features();
}


void MeasuringImpl::update_planes()
{
    m_planes.clear();

    // Now we'll go through all the facets and append Points of facets sharing the same normal.
    // This part is still performed in mesh coordinate system.
    const size_t             num_of_facets  = m_its.indices.size();
    m_face_to_plane.resize(num_of_facets, size_t(-1));
    const std::vector<Vec3f> face_normals   = its_face_normals(m_its);
    const std::vector<Vec3i> face_neighbors = its_face_neighbors(m_its);
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
            if (m_face_to_plane[seed_facet_idx] == size_t(-1)) {
                facet_queue[facet_queue_cnt ++] = seed_facet_idx;
                normal_ptr = &face_normals[seed_facet_idx];
                m_face_to_plane[seed_facet_idx] = m_planes.size();
                m_planes.emplace_back();                
                break;
            }
        if (seed_facet_idx == num_of_facets)
            break; // Everything was visited already

        while (facet_queue_cnt > 0) {
            int facet_idx = facet_queue[-- facet_queue_cnt];
            const stl_normal& this_normal = face_normals[facet_idx];
            if (is_same_normal(this_normal, *normal_ptr)) {
                const Vec3i& face = m_its.indices[facet_idx];

                m_face_to_plane[facet_idx] = m_planes.size() - 1;
                m_planes.back().facets.emplace_back(facet_idx);
                for (int j = 0; j < 3; ++ j)
                    if (int neighbor_idx = face_neighbors[facet_idx][j]; neighbor_idx >= 0 && m_face_to_plane[neighbor_idx] == size_t(-1))
                        facet_queue[facet_queue_cnt ++] = neighbor_idx;
            }
        }

        m_planes.back().normal = normal_ptr->cast<double>();
        std::sort(m_planes.back().facets.begin(), m_planes.back().facets.end());
    }

    assert(std::none_of(m_face_to_plane.begin(), m_face_to_plane.end(), [](size_t val) { return val == size_t(-1); }));

    SurfaceMesh sm(m_its);
    for (int plane_id=0; plane_id < int(m_planes.size()); ++plane_id) {
    //int plane_id = 5; {
        const auto& facets = m_planes[plane_id].facets;
        m_planes[plane_id].borders.clear();
        std::vector<std::array<bool, 3>> visited(facets.size(), {false, false, false});
        
        for (int face_id=0; face_id<int(facets.size()); ++face_id) {
            assert(m_face_to_plane[facets[face_id]] == plane_id);
            for (int edge_id=0; edge_id<3; ++edge_id) {
                if (visited[face_id][edge_id] || m_face_to_plane[face_neighbors[facets[face_id]][edge_id]] == plane_id) {
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
                    while ( m_face_to_plane[sm.face(he)] == plane_id && he != he_orig)
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

    m_planes.erase(std::remove_if(m_planes.begin(), m_planes.end(),
                       [](const PlaneData& p) { return p.borders.empty(); }),
                       m_planes.end());
}






void MeasuringImpl::extract_features()
{
    auto N_to_angle = [](double N) -> double { return 2.*M_PI / N; };
    constexpr double polygon_upper_threshold = N_to_angle(4.5);
    constexpr double polygon_lower_threshold = N_to_angle(8.5);
    std::vector<double> angles;
    std::vector<double> lengths;


    for (int i=0; i<m_planes.size(); ++i) {
        PlaneData& plane = m_planes[i];
        plane.surface_features.clear();
        const Vec3d& normal = plane.normal;

        Eigen::Quaterniond q;
        q.setFromTwoVectors(plane.normal, Vec3d::UnitZ());
        Transform3d trafo = Transform3d::Identity();
        trafo.rotate(q);    
        
        for (const std::vector<Vec3d>& border : plane.borders) {
            assert(border.size() > 1);
            int start_idx = -1;

            // First calculate angles at all the vertices.
            angles.clear();
            lengths.clear();
            for (int i=0; i<int(border.size()); ++i) {
                const Vec3d& v2 = (i == 0 ? border[0] - border[border.size()-1]
                                        : border[i] - border[i-1]);
                const Vec3d& v1 = i == border.size()-1 ? border[0] - border.back()
                                                    : border[i+1] - border[i];
                double angle = atan2(-normal.dot(v1.cross(v2)), -v1.dot(v2)) + M_PI;
                if (angle > M_PI)
                    angle = 2*M_PI - angle;

                angles.push_back(angle);
                lengths.push_back(v2.squaredNorm());
            }
            assert(border.size() == angles.size());
            assert(border.size() == lengths.size());


            bool circle = false;
            std::vector<SurfaceFeature> circles;
            std::vector<std::pair<size_t, size_t>> circles_idxs;
            for (int i=1; i<angles.size(); ++i) {
                if (Slic3r::is_approx(lengths[i], lengths[i-1])
                    && Slic3r::is_approx(angles[i], angles[i-1])
                    && i != angles.size()-1 ) {
                    // circle
                    if (! circle) {
                        circle = true;
                        start_idx = std::max(0, i-2);
                    }
                } else {
                    if (circle) {
                        // Add the circle and remember indices into borders.
                        const auto& [center, radius] = get_center_and_radius(border, start_idx, i, trafo);
                        circles_idxs.emplace_back(start_idx, i);
                        circles.emplace_back(SurfaceFeature(SurfaceFeatureType::Circle, center, plane.normal, std::optional<Vec3d>(), radius));
                        circle = false;
                    }
                }
            }

            // Some of the "circles" may actually be polygons. We want them detected as
            // edges, but also to remember the center and save it into those edges.
            // We will add all such edges manually and delete the detected circles,
            // leaving it in circles_idxs so they are not picked again:
            assert(circles.size() == circles_idxs.size());
            for (int i=circles.size()-1; i>=0; --i) {
                assert(circles_idxs[i].first + 1 < angles.size() - 1); // Check that this is internal point of the circle, not the first, not the last.
                double angle = angles[circles_idxs[i].first + 1];
                if (angle > polygon_lower_threshold) {
                    if (angle < polygon_upper_threshold) {
                        const Vec3d center = std::get<0>(circles[i].get_circle());
                        for (int j=circles_idxs[i].first + 1; j<=circles_idxs[i].second; ++j)
                            plane.surface_features.emplace_back(SurfaceFeature(SurfaceFeatureType::Edge,
                                border[j-1], border[j], std::make_optional(center), 0.));
                    } else {
                        // This will be handled just like a regular edge.
                        circles_idxs.erase(circles_idxs.begin() + i);
                    }
                    circles.erase(circles.begin() + i);
                }
            }






            // We have the circles. Now go around again and pick edges.
            int cidx = 0; // index of next circle in the way
            for (int i=1; i<int(border.size()); ++i) {
                if (cidx < circles_idxs.size() && i > circles_idxs[cidx].first)
                    i = circles_idxs[cidx++].second;
                else plane.surface_features.emplace_back(SurfaceFeature(
                        SurfaceFeatureType::Edge, border[i-1], border[i], std::optional<Vec3d>(), 0.));
            }

            // FIXME Throw away / do not create edges which are parts of circles or
            // which lead to circle points (unless they belong to the same plane.)

            // FIXME Check and merge first and last circle if needed.

            // Now move the circles into the feature list.
            assert(std::all_of(circles.begin(), circles.end(), [](const SurfaceFeature& f) {
                return f.get_type() == SurfaceFeatureType::Circle;
            }));
            plane.surface_features.insert(plane.surface_features.end(), std::make_move_iterator(circles.begin()),
            std::make_move_iterator(circles.end()));
        }

        // The last surface feature is the plane itself.
        plane.surface_features.emplace_back(SurfaceFeature(SurfaceFeatureType::Plane,
            plane.normal, plane.borders.front().front(), std::optional<Vec3d>(), i + 0.0001));

        plane.borders.clear();
        plane.borders.shrink_to_fit();
    }
}



std::vector<SurfaceFeature> MeasuringImpl::get_all_features() const
{
    std::vector<SurfaceFeature> features;
    //PlaneData& plane = m_planes[0];
    for (const PlaneData& plane : m_planes)    
        for (const SurfaceFeature& feature : plane.surface_features)
            features.emplace_back(feature);
    return features;        
}






std::optional<SurfaceFeature> MeasuringImpl::get_feature(size_t face_idx, const Vec3d& point) const
{
    if (face_idx >= m_face_to_plane.size())
        return std::optional<SurfaceFeature>();

    const PlaneData& plane = m_planes[m_face_to_plane[face_idx]];
    
    size_t closest_feature_idx = size_t(-1);
    double min_dist = std::numeric_limits<double>::max();

    MeasurementResult res;
    SurfaceFeature point_sf(point);

    for (size_t i=0; i<plane.surface_features.size() - 1; ++i) {
        // The -1 is there to prevent measuring distance to the plane itself,
        // which is needless and relatively expensive.
        res = get_measurement(plane.surface_features[i], point_sf);
        if (res.distance_strict) { // TODO: this should become an assert after all combinations are implemented.
            double dist = *res.distance_strict;
            if (dist < feature_hover_limit && dist < min_dist) {
                min_dist = std::min(dist, min_dist);
                closest_feature_idx = i;
            }
        }
    }

    if (closest_feature_idx != size_t(-1)) {
        const SurfaceFeature& f = plane.surface_features[closest_feature_idx];
        if (f.get_type() == SurfaceFeatureType::Edge) {
            // If this is an edge, check if we are not close to the endpoint. If so,
            // we will include the endpoint as well.
            constexpr double limit_sq = edge_endpoint_limit * edge_endpoint_limit;
            const auto& [sp, ep] = f.get_edge();
            if ((point-sp).squaredNorm() < limit_sq)
                return std::make_optional(SurfaceFeature(sp));
            if ((point-ep).squaredNorm() < limit_sq)
                return std::make_optional(SurfaceFeature(ep));
        }
        return std::make_optional(f);
    }

    // Nothing detected, return the plane as a whole.
    assert(plane.surface_features.back().get_type() == SurfaceFeatureType::Plane);
    return std::make_optional(plane.surface_features.back());
}





std::vector<std::vector<int>> MeasuringImpl::get_planes_triangle_indices() const
{
    std::vector<std::vector<int>> out;
    for (const PlaneData& plane : m_planes)
        out.emplace_back(plane.facets);        
    return out;
}













Measuring::Measuring(const indexed_triangle_set& its)
: priv{std::make_unique<MeasuringImpl>(its)}
{}

Measuring::~Measuring() {}


std::vector<SurfaceFeature> Measuring::get_all_features() const
{
    return priv->get_all_features();
}


std::optional<SurfaceFeature> Measuring::get_feature(size_t face_idx, const Vec3d& point) const
{
    return priv->get_feature(face_idx, point);
}



std::vector<std::vector<int>> Measuring::get_planes_triangle_indices() const
{
    return priv->get_planes_triangle_indices();
}














MeasurementResult get_measurement(const SurfaceFeature& a, const SurfaceFeature& b)
{
    assert(a.get_type() != SurfaceFeatureType::Undef && b.get_type() != SurfaceFeatureType::Undef);

    const bool swap = int(a.get_type()) > int(b.get_type());
    const SurfaceFeature& f1 = swap ? b : a;
    const SurfaceFeature& f2 = swap ? a : b;

    MeasurementResult result;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    if (f1.get_type() == SurfaceFeatureType::Point) {
        if (f2.get_type() == SurfaceFeatureType::Point) {
            Vec3d diff = (f2.get_point() - f1.get_point());
            result.distance_strict = diff.norm();
            result.distance_xyz = diff;
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Edge) {
            const auto& [s,e] = f2.get_edge();
            Eigen::ParametrizedLine<double, 3> line(s, (e-s).normalized());
            double dist_inf = line.distance(f1.get_point());
            Vec3d proj = line.projection(f1.get_point());
            double len_sq = (e-s).squaredNorm();
            double dist_start_sq = (proj-s).squaredNorm();
            double dist_end_sq = (proj-e).squaredNorm();
            if (dist_start_sq < len_sq && dist_end_sq < len_sq) {
                // projection falls on the line - the strict distance is the same as infinite
                result.distance_strict = std::make_optional(dist_inf);
            } else { // the result is the closer of the endpoints
                result.distance_strict = std::make_optional(std::sqrt(std::min(dist_start_sq, dist_end_sq) + dist_inf));
            }
            result.distance_infinite = std::make_optional(dist_inf);
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Circle) {
            // Find a plane containing normal, center and the point.
            const auto& [c, radius, n] = f2.get_circle();
            Eigen::Hyperplane<double, 3> circle_plane(n, c);
            Vec3d proj = circle_plane.projection(f1.get_point());
            double dist = std::sqrt(std::pow((proj - c).norm() - radius, 2.) +
                (f1.get_point() - proj).squaredNorm());

            result.distance_strict = std::make_optional(dist);
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Plane) {
            const auto& [idx, normal, pt] = f2.get_plane();
            Eigen::Hyperplane<double, 3> plane(normal, pt);
            result.distance_infinite = plane.absDistance(f1.get_point());
            // TODO: result.distance_strict =
        }
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    } else if (f1.get_type() == SurfaceFeatureType::Edge) {
        if (f2.get_type() == SurfaceFeatureType::Edge) {
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Circle) {
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Plane) {
        }
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    } else if (f1.get_type() == SurfaceFeatureType::Circle) {
        if (f2.get_type() == SurfaceFeatureType::Circle) {
    ///////////////////////////////////////////////////////////////////////////
        } else if (f2.get_type() == SurfaceFeatureType::Plane) {
        }
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    } else if (f1.get_type() == SurfaceFeatureType::Plane) {
        assert(f2.get_type() == SurfaceFeatureType::Plane);

        const auto& [idx1, normal1, pt1] = f1.get_plane();
        const auto& [idx2, normal2, pt2] = f2.get_plane();
        double angle = 0.;

        if (! normal1.isApprox(normal2)) {
            // The planes are parallel, calculate distance.
            Eigen::Hyperplane<double, 3> plane(normal1, pt1);
            result.distance_infinite = plane.absDistance(pt2);
        } else {
            // Planes are not parallel, calculate angle.
            angle = std::acos(std::abs(normal1.dot(normal2)));
        }
        result.angle = angle;
    }


    
    return result;
}









} // namespace Measure
} // namespace Slic3r
