#include "PointCloud.hpp"

#include "libslic3r/Geometry.hpp"
#include "libslic3r/Tesselate.hpp"

#include <igl/random_points_on_mesh.h>

namespace Slic3r { namespace branchingtree {

std::optional<Vec3f> find_merge_pt(const Vec3f &A,
                                   const Vec3f &B,
                                   float        critical_angle)
{
    // The idea is that A and B both have their support cones. But searching
    // for the intersection of these support cones is difficult and its enough
    // to reduce this problem to 2D and search for the intersection of two
    // rays that merge somewhere between A and B. The 2D plane is a vertical
    // slice of the 3D scene where the 2D Y axis is equal to the 3D Z axis and
    // the 2D X axis is determined by the XY direction of the AB vector.
    //
    // Z^
    //  |    A *
    //  |     . .   B *
    //  |    .   .   . .
    //  |   .     . .   .
    //  |  .       x     .
    //  -------------------> XY

    // Determine the transformation matrix for the 2D projection:
    Vec3f diff = {B.x() - A.x(), B.y() - A.y(), 0.f};
    Vec3f dir  = diff.normalized(); // TODO: avoid normalization

    Eigen::Matrix<float, 2, 3> tr2D;
    tr2D.row(0) = Vec3f{dir.x(), dir.y(), dir.z()};
    tr2D.row(1) = Vec3f{0.f, 0.f, 1.f};

    // Transform the 2 vectors A and B into 2D vector 'a' and 'b'. Here we can
    // omit 'a', pretend that its the origin and use BA as the vector b.
    Vec2f b = tr2D * (B - A);

    // Get the square sine of the ray emanating from 'a' towards 'b'. This ray might
    // exceed the allowed angle but that is corrected subsequently.
    // The sign of the original sine is also needed, hence b.y is multiplied by
    // abs(b.y)
    float b_sqn = b.squaredNorm();
    float sin2sig_a = b_sqn > EPSILON ? (b.y() * std::abs(b.y())) / b_sqn : 0.f;

    // sine2 from 'b' to 'a' is the opposite of sine2 from a to b
    float sin2sig_b = -sin2sig_a;

    // Derive the allowed angles from the given critical angle.
    // critical_angle is measured from the horizontal X axis.
    // The rays need to go downwards which corresponds to negative angles

    float sincrit = std::sin(critical_angle);  // sine of the critical angle
    float sin2crit = -sincrit * sincrit;       // signed sine squared
    sin2sig_a = std::min(sin2sig_a, sin2crit); // Do the angle saturation of both rays
    sin2sig_b = std::min(sin2sig_b, sin2crit); //
    float sin2_a = std::abs(sin2sig_a);        // Get cosine squared values
    float sin2_b = std::abs(sin2sig_b);
    float cos2_a = 1.f - sin2_a;
    float cos2_b = 1.f - sin2_b;

    // Derive the new direction vectors. This is by square rooting the sin2
    // and cos2 values and restoring the original signs
    Vec2f Da = {std::copysign(std::sqrt(cos2_a), b.x()), std::copysign(std::sqrt(sin2_a), sin2sig_a)};
    Vec2f Db = {-std::copysign(std::sqrt(cos2_b), b.x()), std::copysign(std::sqrt(sin2_b), sin2sig_b)};

    // Determine where two rays ([0, 0], Da), (b, Db) intersect.
    // Based on
    // https://stackoverflow.com/questions/27459080/given-two-points-and-two-direction-vectors-find-the-point-where-they-intersect
    // One ray is emanating from (0, 0) so the formula is simplified
    double t1 = (Db.y() * b.x() - b.y() * Db.x()) /
                (Da.x() * Db.y() - Da.y() * Db.x());

    Vec2f mp = t1 * Da;
    Vec3f Mp = A + tr2D.transpose() * mp;

    return t1 >= 0.f ? Mp : Vec3f{};
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

std::vector<Node> sample_mesh(const indexed_triangle_set &its, double radius)
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
        int face_id = FI(i);

        if (face_id < 0 || face_id >= int(its.indices.size()))
            continue;

        Vec3i face = its.indices[face_id];

        if (face(0) >= int(its.vertices.size()) ||
            face(1) >= int(its.vertices.size()) ||
            face(2) >= int(its.vertices.size()))
            continue;

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
                       std::vector<Node>           support_leafs,
                       const Properties           &props)
    : PointCloud{sample_mesh(M, props.sampling_radius()),
                 sample_bed(props.bed_shape(),
                            props.ground_level(),
                            props.sampling_radius()),
                 std::move(support_leafs), props}
{}

PointCloud::PointCloud(std::vector<Node> meshpts,
                       std::vector<Node> bedpts,
                       std::vector<Node> support_leafs,
                       const Properties &props)
    : m_leafs{std::move(support_leafs)}
    , m_meshpoints{std::move(meshpts)}
    , m_bedpoints{std::move(bedpts)}
    , m_props{props}
    , cos2bridge_slope{std::cos(props.max_slope()) *
                       std::abs(std::cos(props.max_slope()))}
    , MESHPTS_BEGIN{m_bedpoints.size()}
    , LEAFS_BEGIN{MESHPTS_BEGIN + m_meshpoints.size()}
    , JUNCTIONS_BEGIN{LEAFS_BEGIN + m_leafs.size()}
    , m_searchable_indices(JUNCTIONS_BEGIN + m_junctions.size(), true)
    , m_queue_indices(JUNCTIONS_BEGIN + m_junctions.size(), Unqueued)
    , m_reachable_cnt{JUNCTIONS_BEGIN + m_junctions.size()}
{
    for (size_t i = 0; i < m_bedpoints.size(); ++i) {
        m_bedpoints[i].id = int(i);
        m_ktree.insert({m_bedpoints[i].pos, i});
    }
    
    for (size_t i = 0; i < m_meshpoints.size(); ++i) {
        Node &n = m_meshpoints[i];
        n.id = int(MESHPTS_BEGIN + i);
        m_ktree.insert({n.pos, n.id});
    }
    
    for (size_t i = 0; i < m_leafs.size(); ++i) {
        Node &n = m_leafs[i];
        n.id = int(LEAFS_BEGIN + i);
        m_ktree.insert({n.pos, n.id});
    }
}

float PointCloud::get_distance(const Vec3f &p, size_t node_id) const
{
    auto t = get_type(node_id);
    auto ret = std::numeric_limits<float>::infinity();
    const auto &node = get(node_id);
    
    switch (t) {
    case MESH:
    case BED: {
        // Points of mesh or bed which are outside of the support cone of
        // 'pos' must be discarded.
        if (is_outside_support_cone(p, node.pos))
            ret = std::numeric_limits<float>::infinity();
        else
            ret  = (node.pos - p).norm();
        
        break;
    }
    case LEAF:
    case JUNCTION:{
        auto mergept = find_merge_pt(p, node.pos, m_props.max_slope());
        double maxL2 = m_props.max_branch_length() * m_props.max_branch_length();

        if (!mergept || mergept->z() < (m_props.ground_level() + 2 * node.Rmin))
            ret = std::numeric_limits<float>::infinity();
        else if (double a = (node.pos - *mergept).squaredNorm(),
                 b        = (p - *mergept).squaredNorm();
                 a < maxL2 && b < maxL2)
            ret = std::sqrt(b);

        break;
    }
    case NONE:
        ;
    }
    
    // Setting the ret val to infinity will effectively discard this
    // connection of nodes. max_branch_length property is used here
    // to discard node=>node and node=>mesh connections longer than this
    // property.
    if (t != BED && ret > m_props.max_branch_length())
        ret = std::numeric_limits<float>::infinity();
    
    return ret;
}

}} // namespace Slic3r::branchingtree
