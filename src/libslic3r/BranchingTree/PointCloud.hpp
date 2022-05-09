#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include "BranchingTree.hpp"

#include "libslic3r/KDTreeIndirect.hpp"
#include "libslic3r/MutablePriorityQueue.hpp"

namespace Slic3r { namespace branchingtree {

std::optional<Vec3f> find_merge_pt(const Vec3f &A,
                                   const Vec3f &B,
                                   float        max_slope);

void to_eigen_mesh(const indexed_triangle_set &its,
                   Eigen::MatrixXd            &V,
                   Eigen::MatrixXi            &F);

std::vector<Node> sample_mesh(const indexed_triangle_set &its,
                                  double                      radius);

std::vector<Node> sample_bed(const ExPolygons &bed,
                                 float             z,
                                 double            radius = 1.);

enum PtType { SUPP, MESH, BED, JUNCTION, NONE };

// A cloud of points including support points, mesh points, junction points
// and anchor points on the bed. Junction points can be added or removed, all
// the other point types are established on creation and remain unchangeable.
class PointCloud {
    std::vector<Node> m_leafs, m_junctions, m_meshpoints, m_bedpoints;

    const branchingtree::Properties &m_props;

    const double cos2bridge_slope;
    const size_t MESHPTS_BEGIN, SUPP_BEGIN, JUNCTIONS_BEGIN;

private:

    // These vectors have the same size as there are indices for nodes to keep
    // access complexity constant. WARN: there might be cache non-locality costs
    std::vector<bool> m_searchable_indices; // searchable flag value of a node
    std::vector<size_t> m_queue_indices;    // queue id of a node if queued

    size_t m_reachable_cnt;

    struct CoordFn
    {
        const PointCloud *self;
        CoordFn(const PointCloud *s) : self{s} {}
        float operator()(size_t nodeid, size_t dim) const
        {
            return self->get(nodeid).pos(int(dim));
        }
    };

    KDTreeIndirect<3, float, CoordFn> m_ktree;

    bool is_outside_support_cone(const Vec3f &supp, const Vec3f &pt)
    {
        Vec3d D = (pt - supp).cast<double>();
        double dot_sq = -D.z() * std::abs(-D.z());

        return dot_sq < D.squaredNorm() * cos2bridge_slope;
    }

    static constexpr auto UNQUEUED = size_t(-1);

    template<class PC>
    static auto *get_node(PC &&pc, size_t id)
    {
        auto *ret = decltype(pc.m_bedpoints.data())(nullptr);

        switch(pc.get_type(id)) {
        case BED: ret = &pc.m_bedpoints[id]; break;
        case MESH: ret = &pc.m_meshpoints[id - pc.MESHPTS_BEGIN]; break;
        case SUPP: ret = &pc.m_leafs [id - pc.SUPP_BEGIN]; break;
        case JUNCTION: ret = &pc.m_junctions[id - pc.JUNCTIONS_BEGIN]; break;
        case NONE: assert(false);
        }

        return ret;
    }

public:

    struct ZCompareFn
    {
        const PointCloud *self;
        ZCompareFn(const PointCloud *s) : self{s} {}
        bool operator()(size_t node_a, size_t node_b) const
        {
            return self->get(node_a).pos.z() > self->get(node_b).pos.z();
        }
    };

    PointCloud(const indexed_triangle_set & M,
               std::vector<Node>        support_leafs,
               const Properties &           props);

    PtType get_type(size_t node_id) const
    {
        PtType ret = NONE;

        if (node_id < MESHPTS_BEGIN && !m_bedpoints.empty()) ret = BED;
        else if (node_id < SUPP_BEGIN && !m_meshpoints.empty()) ret = MESH;
        else if (node_id < JUNCTIONS_BEGIN && !m_leafs.empty()) ret = SUPP;
        else if (node_id >= JUNCTIONS_BEGIN && !m_junctions.empty()) ret = JUNCTION;

        return ret;
    }

    const Node &get(size_t node_id) const
    {
        return *get_node(*this, node_id);
    }

    Node &get(size_t node_id)
    {
        return *get_node(*this, node_id);
    }

    const Node *find(size_t node_id) const { return get_node(*this, node_id); }
    Node *find(size_t node_id) { return get_node(*this, node_id); }

    // Return the original index of a leaf in the input array, if the given
    // node id is indeed of type SUPP
    int get_leaf_id(size_t node_id) const
    {
        return node_id >= SUPP_BEGIN && node_id < JUNCTIONS_BEGIN ?
                   node_id - SUPP_BEGIN :
                   Node::ID_NONE;
    }

    size_t get_queue_idx(size_t node_id) const { return m_queue_indices[node_id]; }

    float get_distance(const Vec3f &p, size_t node);

    size_t next_junction_id() const
    {
        return JUNCTIONS_BEGIN + m_junctions.size();
    }

    size_t insert_junction(const Node &p)
    {
        size_t new_id = next_junction_id();
        m_junctions.emplace_back(p);
        m_junctions.back().id = int(new_id);
        m_searchable_indices.emplace_back(true);
        m_queue_indices.emplace_back(UNQUEUED);
        ++m_reachable_cnt;

        return new_id;
    }

    void remove_node(size_t node_id)
    {
        m_searchable_indices[node_id] = false;
        m_queue_indices[node_id] = UNQUEUED;
        --m_reachable_cnt;
    }

    size_t reachable_count() const { return m_reachable_cnt; }

    template<class Fn> void foreach_reachable(const Vec3f &pos, Fn &&visitor)
    {
        auto closest_anchors =
            find_closest_points<5>(m_ktree, pos, [this, &pos](size_t id) {
                return m_searchable_indices[id] &&
                       !is_outside_support_cone(pos, get(id).pos);
            });

        for (size_t anchor : closest_anchors)
            if (anchor != m_ktree.npos)
                visitor(anchor, get_distance(pos, anchor));

        for (size_t i = SUPP_BEGIN; i < m_searchable_indices.size(); ++i)
            if (m_searchable_indices[i])
                visitor(i, get_distance(pos, i));
    }

    auto start_queue()
    {
        auto ptsqueue = make_mutable_priority_queue<size_t, false>(
            [this](size_t el, size_t idx) { m_queue_indices[el] = idx; },
            ZCompareFn{this});

        ptsqueue.reserve(m_leafs.size());
        size_t iend = SUPP_BEGIN + m_leafs.size();
        for (size_t i = SUPP_BEGIN; i < iend; ++i)
            ptsqueue.push(i);

        return ptsqueue;
    }

    const Properties & properties() const { return m_props; }
};

template<class PC, class Fn> void traverse(PC &&pc, size_t root, Fn &&fn)
{
    if (auto nodeptr = pc.find(root); nodeptr != nullptr) {
        auto &nroot = *nodeptr;
        fn(nroot);
        if (nroot.left  >= 0) traverse(pc, nroot.left, fn);
        if (nroot.right >= 0) traverse(pc, nroot.right, fn);
    }
}

bool build_tree(PointCloud &pcloud, Builder &builder);

}} // namespace Slic3r::branchingtree

#endif // POINTCLOUD_HPP
