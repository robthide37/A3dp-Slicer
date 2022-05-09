#include "PointCloud.hpp"

#include "libslic3r/Geometry.hpp"
#include "libslic3r/Tesselate.hpp"

#include <igl/random_points_on_mesh.h>

namespace Slic3r { namespace branchingtree {

std::optional<Vec3f> find_merge_pt(const Vec3f &A,
                                   const Vec3f &B,
                                   float        max_slope)
{
    Vec3f Da = (B - A).normalized(), Db = -Da;
    auto [polar_da, azim_da] = Geometry::dir_to_spheric(Da);
    auto [polar_db, azim_db] = Geometry::dir_to_spheric(Db);
    polar_da = std::max(polar_da, float(PI) - max_slope);
    polar_db = std::max(polar_db, float(PI) - max_slope);

    Da = Geometry::spheric_to_dir<float>(polar_da, azim_da);
    Db = Geometry::spheric_to_dir<float>(polar_db, azim_db);

    double t1 =
        (A.z() * Db.x() + Db.z() * B.x() - B.z() * Db.x() - Db.z() * A.x()) /
        (Da.x() * Db.z() - Da.z() * Db.x());

    double t2 = (A.x() + Da.x() * t1 - B.x()) / Db.x();

    return t1 > 0. && t2 > 0. ? A + t1 * Da : std::optional<Vec3f>{};
}

void to_eigen_mesh(const indexed_triangle_set &its,
                   Eigen::MatrixXd            &V,
                   Eigen::MatrixXi            &F)
{
    V.resize(its.vertices.size(), 3);
    F.resize(its.indices.size(), 3);
    for (unsigned int i = 0; i < its.indices.size(); ++i)
        F.row(i) = its.indices[i];

    for (unsigned int i = 0; i < its.vertices.size(); ++i)
        V.row(i) = its.vertices[i].cast<double>();
}

std::vector<Node> sample_mesh(const indexed_triangle_set &its,
                                  double                      radius)
{
    std::vector<Node> ret;

    double surface_area = 0.;
    for (const Vec3i &face : its.indices) {
        std::array<Vec3f, 3> tri = {its.vertices[face(0)],
                                    its.vertices[face(1)],
                                    its.vertices[face(2)]};

        auto U = tri[1] - tri[0], V = tri[2] - tri[0];
        surface_area += 0.5 * U.cross(V).norm();
    }

    int N = surface_area / (PI * radius * radius);

    Eigen::MatrixXd B;
    Eigen::MatrixXi FI;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    to_eigen_mesh(its, V, F);
    igl::random_points_on_mesh(N, V, F, B, FI);

    ret.reserve(size_t(N));
    for (int i = 0; i < FI.size(); i++) {
        Vec3i face = its.indices[FI(i)];

        Vec3f c = B.row(i)(0) * its.vertices[face(0)] +
                  B.row(i)(1) * its.vertices[face(1)] +
                  B.row(i)(2) * its.vertices[face(2)];

        ret.emplace_back(c);
    }

    return ret;
}

std::vector<Node> sample_bed(const ExPolygons &bed, float z, double radius)
{
    std::vector<Vec3f> ret;

    auto triangles = triangulate_expolygons_3d(bed, z);
    indexed_triangle_set its;
    its.vertices.reserve(triangles.size());

    for (size_t i = 0; i < triangles.size(); i += 3) {
        its.vertices.emplace_back(triangles[i].cast<float>());
        its.vertices.emplace_back(triangles[i + 1].cast<float>());
        its.vertices.emplace_back(triangles[i + 2].cast<float>());

        its.indices.emplace_back(i, i + 1, i + 2);
    }

    return sample_mesh(its, radius);
}

PointCloud::PointCloud(const indexed_triangle_set &M,
                       std::vector<Node>       support_leafs,
                       const Properties           &props)
    : m_leafs{std::move(support_leafs)}
    , m_meshpoints{sample_mesh(M, props.sampling_radius())}
    , m_bedpoints{sample_bed(props.bed_shape(),
                             props.ground_level(),
                             props.sampling_radius())}
    , m_props{props}
    , cos2bridge_slope{std::cos(props.max_slope()) *
                       std::abs(std::cos(props.max_slope()))}
    , MESHPTS_BEGIN{m_bedpoints.size()}
    , SUPP_BEGIN{MESHPTS_BEGIN + m_meshpoints.size()}
    , JUNCTIONS_BEGIN{SUPP_BEGIN + m_leafs.size()}
    , m_searchable_indices(JUNCTIONS_BEGIN, true)
    , m_queue_indices(JUNCTIONS_BEGIN, UNQUEUED)
    , m_reachable_cnt{JUNCTIONS_BEGIN}
    , m_ktree{CoordFn{this}, SUPP_BEGIN} // Only for bed and mesh points
{
    for (size_t i = 0; i < m_bedpoints.size(); ++i)
        m_bedpoints[i].id = int(i);

    for (size_t i = 0; i < m_meshpoints.size(); ++i)
        m_meshpoints[i].id = int(MESHPTS_BEGIN + i);

    for (size_t i = 0; i < m_leafs.size(); ++i)
        m_leafs[i].id = int(SUPP_BEGIN + i);
}

float PointCloud::get_distance(const Vec3f &p, size_t node)
{
    auto t = get_type(node);
    auto ret = std::numeric_limits<float>::infinity();

    switch (t) {
    case MESH:
    case BED: {
        // Points of mesh or bed which are outside of the support cone of
        // 'pos' must be discarded.
        if (is_outside_support_cone(p, get(node).pos))
            ret = std::numeric_limits<float>::infinity();
        else
            ret  = (get(node).pos - p).norm();

        break;
    }
    case SUPP:
    case JUNCTION:{
        auto mergept = find_merge_pt(p, get(node).pos, m_props.max_slope());
        if (!mergept || mergept->z() < (m_props.ground_level() + 2 * get(node).Rmin))
            ret = std::numeric_limits<float>::infinity();
        else
            ret = (p - *mergept).norm();

        break;
    }
    case NONE:
        ;
    }

       // Setting the ret val to infinity will effectively discard this
       // connection of nodes. max_branch_length property if used here
       // to discard node=>node and node=>mesh connections longer than this
       // property.
    if (t != BED && ret > m_props.max_branch_length())
        ret = std::numeric_limits<float>::infinity();

    return ret;
}

}} // namespace Slic3r::branchingtree
