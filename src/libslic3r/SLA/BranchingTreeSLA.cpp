#include "BranchingTreeSLA.hpp"

#include "libslic3r/Execution/ExecutionTBB.hpp"

#include "libslic3r/KDTreeIndirect.hpp"

#include "SupportTreeUtils.hpp"
#include "BranchingTree/PointCloud.hpp"

#include "Pad.hpp"

#include <map>

namespace Slic3r { namespace sla {

class BranchingTreeBuilder: public branchingtree::Builder {
    SupportTreeBuilder &m_builder;
    const SupportableMesh  &m_sm;
    const branchingtree::PointCloud &m_cloud;

    std::set<int /*ground node_id that was already processed*/> m_ground_mem;

    // Establish an index of
    using PointIndexEl = std::pair<Vec3f, unsigned>;
    boost::geometry::index::
        rtree<PointIndexEl, boost::geometry::index::rstar<16, 4> /* ? */>
            m_pillar_index;

    std::vector<branchingtree::Node> m_pillars; // to put an index over them

    // cache succesfull ground connections
    std::map<size_t, GroundConnection>    m_gnd_connections;

    // Scaling of the input value 'widening_factor:<0, 1>' to produce resonable
    // widening behaviour
    static constexpr double WIDENING_SCALE = 0.02;

    double get_radius(const branchingtree::Node &j) const
    {
        double w = WIDENING_SCALE * m_sm.cfg.pillar_widening_factor * j.weight;

        return std::min(m_sm.cfg.base_radius_mm, double(j.Rmin) + w);
    }

    std::vector<size_t>  m_unroutable_pinheads;

    void build_subtree(size_t root)
    {
        traverse(m_cloud, root, [this](const branchingtree::Node &node) {
            if (node.left >= 0 && node.right >= 0) {
                auto nparent = m_cloud.get(node.id);
                auto nleft = m_cloud.get(node.left);
                auto nright = m_cloud.get(node.right);
                Vec3d from1d = nleft.pos.cast<double>();
                Vec3d from2d = nright.pos.cast<double>();
                Vec3d tod    = nparent.pos.cast<double>();
                double mergeR = get_radius(nparent);
                double leftR  = get_radius(nleft);
                double rightR = get_radius(nright);

                m_builder.add_diffbridge(from1d, tod, leftR, mergeR);
                m_builder.add_diffbridge(from2d, tod, rightR, mergeR);
                m_builder.add_junction(tod, mergeR);
            } else if (int child = node.left + node.right + 1; child >= 0) {
                auto from = m_cloud.get(child);
                auto to   = m_cloud.get(node.id);
                auto tod  = to.pos.cast<double>();
                double toR = get_radius(to);
                m_builder.add_diffbridge(from.pos.cast<double>(),
                                         tod,
                                         get_radius(from),
                                         toR);
                m_builder.add_junction(tod, toR);
            }
        });
    }

    void discard_subtree(size_t root)
    {
        // Discard all the support points connecting to this branch.
        // As a last resort, try to route child nodes to ground and stop
        // traversing if any child branch succeeds.
        traverse(m_cloud, root, [this](const branchingtree::Node &node) {
            branchingtree::TraverseReturnT ret{true, true};

            int suppid_parent = m_cloud.get_leaf_id(node.id);
            int suppid_left   = branchingtree::Node::ID_NONE;
            int suppid_right  = branchingtree::Node::ID_NONE;

            branchingtree::Node dst = node;
            dst.weight += node.pos.z();
            dst.Rmin = std::max(node.Rmin, dst.Rmin);

            if (node.left >= 0 && add_ground_bridge(m_cloud.get(node.left), dst))
                ret.to_left = false;
            else
                suppid_left = m_cloud.get_leaf_id(node.left);

            if (node.right >= 0 && add_ground_bridge(m_cloud.get(node.right), dst))
                ret.to_right = false;
            else
                suppid_right = m_cloud.get_leaf_id(node.right);

            if (suppid_parent >= 0)
                m_unroutable_pinheads.emplace_back(suppid_parent);
            if (suppid_left >= 0)
                m_unroutable_pinheads.emplace_back(suppid_left);
            if (suppid_right >= 0)
                m_unroutable_pinheads.emplace_back(suppid_right);

            return ret;
        });
    }

    std::optional<PointIndexEl>
    search_for_existing_pillar(const branchingtree::Node &from) const
    {
        namespace bgi = boost::geometry::index;

        struct Output
        {
            std::optional<PointIndexEl> &res;

            Output &operator*() { return *this; }
            Output &operator=(const PointIndexEl &el) { res = el; return *this; }
            Output &operator++() { return *this; }
        };

        std::optional<PointIndexEl> result;

        auto filter = bgi::satisfies([this, &from](const PointIndexEl &e) {
            assert(e.second < m_pillars.size());

            auto len = (from.pos - e.first).norm();
            const branchingtree::Node &to = m_pillars[e.second];
            double sd = m_sm.cfg.safety_distance_mm;

            Beam beam{Ball{from.pos.cast<double>(), get_radius(from)},
                      Ball{e.first.cast<double>(), get_radius(to)}};

            return !branchingtree::is_occupied(to) &&
                   len < m_sm.cfg.max_bridge_length_mm &&
                   !m_cloud.is_outside_support_cone(from.pos, e.first) &&
                   beam_mesh_hit(ex_tbb, m_sm.emesh, beam, sd).distance() > len;
        });

        m_pillar_index.query(filter && bgi::nearest(from.pos, 1),
                             Output{result});

        return result;
    }

public:
    BranchingTreeBuilder(SupportTreeBuilder          &builder,
                     const SupportableMesh       &sm,
                     const branchingtree::PointCloud &cloud)
        : m_builder{builder}, m_sm{sm}, m_cloud{cloud}
    {}

    bool add_bridge(const branchingtree::Node &from,
                    const branchingtree::Node &to) override;

    bool add_merger(const branchingtree::Node &node,
                    const branchingtree::Node &closest,
                    const branchingtree::Node &merge_node) override;

    bool add_ground_bridge(const branchingtree::Node &from,
                           const branchingtree::Node &/*to*/) override;

    bool add_mesh_bridge(const branchingtree::Node &from,
                         const branchingtree::Node &to) override;

    void report_unroutable(const branchingtree::Node &j) override
    {
        BOOST_LOG_TRIVIAL(error) << "Cannot route junction at " << j.pos.x()
                                 << " " << j.pos.y() << " " << j.pos.z();

        // Discard all the support points connecting to this branch.
        discard_subtree(j.id);
    }

    const std::vector<size_t>& unroutable_pinheads() const
    {
        return m_unroutable_pinheads;
    }

    bool is_valid() const override { return !m_builder.ctl().stopcondition(); }

    const std::vector<branchingtree::Node> & pillars() const { return m_pillars; }

    const GroundConnection *ground_conn(size_t pillar) const
    {
        const GroundConnection *ret = nullptr;

        auto it = m_gnd_connections.find(pillar);
        if (it != m_gnd_connections.end())
            ret = &it->second;

        return ret;
    }
};

bool BranchingTreeBuilder::add_bridge(const branchingtree::Node &from,
                                      const branchingtree::Node &to)
{
    Vec3d fromd = from.pos.cast<double>(), tod = to.pos.cast<double>();
    double fromR = get_radius(from), toR = get_radius(to);
    Beam beam{Ball{fromd, fromR}, Ball{tod, toR}};
    auto   hit = beam_mesh_hit(ex_tbb, m_sm.emesh, beam,
                               m_sm.cfg.safety_distance_mm);

    bool ret = hit.distance() > (tod - fromd).norm();

    return ret;
}

bool BranchingTreeBuilder::add_merger(const branchingtree::Node &node,
                                      const branchingtree::Node &closest,
                                      const branchingtree::Node &merge_node)
{
    Vec3d from1d = node.pos.cast<double>(),
          from2d = closest.pos.cast<double>(),
          tod    = merge_node.pos.cast<double>();

    double mergeR   = get_radius(merge_node);
    double nodeR    = get_radius(node);
    double closestR = get_radius(closest);
    Beam beam1{Ball{from1d, nodeR}, Ball{tod, mergeR}};
    Beam beam2{Ball{from2d, closestR}, Ball{tod, mergeR}};

    auto sd = m_sm.cfg.safety_distance_mm ;
    auto hit1 = beam_mesh_hit(ex_tbb, m_sm.emesh, beam1, sd);
    auto hit2 = beam_mesh_hit(ex_tbb, m_sm.emesh, beam2, sd);

    bool ret = hit1.distance() > (tod - from1d).norm() &&
               hit2.distance() > (tod - from2d).norm();

    return ret;
}

bool BranchingTreeBuilder::add_ground_bridge(const branchingtree::Node &from,
                                             const branchingtree::Node &to)
{
    bool ret = false;

    namespace bgi = boost::geometry::index;

    auto it = m_ground_mem.find(from.id);
    if (it == m_ground_mem.end()) {
//        std::optional<PointIndexEl> result = search_for_existing_pillar(from);

        sla::Junction j{from.pos.cast<double>(), get_radius(from)};
//        if (!result) {
            auto conn = optimize_ground_connection(
                                           ex_tbb,
                                           m_sm,
                                           j,
                                           get_radius(to));

            if (conn) {
//                Junction connlast = conn.path.back();
//                branchingtree::Node n{connlast.pos.cast<float>(), float(connlast.r)};
//                n.left = from.id;
                m_pillars.emplace_back(from);
//                m_pillar_index.insert({n.pos, m_pillars.size() - 1});
                m_gnd_connections[m_pillars.size() - 1] = conn;

                ret = true;
            }
//        } else {
//            const auto &resnode = m_pillars[result->second];
//            m_builder.add_diffbridge(j.pos, resnode.pos.cast<double>(), j.r, get_radius(resnode));
//            m_pillars[result->second].right = from.id;
//            ret = true;
//        }

        // Remember that this node was tested if can go to ground, don't
        // test it with any other destination ground point because
        // it is unlikely that search_ground_route would find a better solution
        m_ground_mem.insert(from.id);
    }

    if (ret) {
        build_subtree(from.id);
    }

    return ret;
}

bool BranchingTreeBuilder::add_mesh_bridge(const branchingtree::Node &from,
                                           const branchingtree::Node &to)
{
    if (from.weight > m_sm.cfg.max_weight_on_model_support)
        return false;

    sla::Junction fromj = {from.pos.cast<double>(), get_radius(from)};

    auto anchor = m_sm.cfg.ground_facing_only ?
                      std::optional<Anchor>{} : // If no mesh connections are allowed
                      calculate_anchor_placement(ex_tbb, m_sm, fromj,
                                                 to.pos.cast<double>());

    if (anchor) {
        sla::Junction toj = {anchor->junction_point(), anchor->r_back_mm};

        auto hit = beam_mesh_hit(ex_tbb, m_sm.emesh,
                                 Beam{{fromj.pos, fromj.r}, {toj.pos, toj.r}}, 0.);

        if (hit.distance() > distance(fromj.pos, toj.pos)) {
            m_builder.add_diffbridge(fromj.pos, toj.pos, fromj.r, toj.r);
            m_builder.add_anchor(*anchor);

            build_subtree(from.id);
        } else {
            anchor.reset();
        }
    }

    return bool(anchor);
}

inline void build_pillars(SupportTreeBuilder &builder,
                          BranchingTreeBuilder &vbuilder,
                          const SupportableMesh &sm)
{
    for (size_t pill_id = 0; pill_id < vbuilder.pillars().size(); ++pill_id) {
        auto * conn = vbuilder.ground_conn(pill_id);
        if (conn)
            build_ground_connection(builder, sm, *conn);
    }
}

void create_branching_tree(SupportTreeBuilder &builder, const SupportableMesh &sm)
{
    auto coordfn = [&sm](size_t id, size_t dim) { return sm.pts[id].pos(dim); };
    KDTreeIndirect<3, float, decltype (coordfn)> tree{coordfn, sm.pts.size()};

    auto nondup_idx = non_duplicate_suppt_indices(tree, sm.pts, 0.1);
    std::vector<std::optional<Head>> heads(nondup_idx.size());
    auto leafs = reserve_vector<branchingtree::Node>(nondup_idx.size());

    execution::for_each(
        ex_tbb, size_t(0), nondup_idx.size(),
        [&sm, &heads, &nondup_idx, &builder](size_t i) {
            if (!builder.ctl().stopcondition())
                heads[i] = calculate_pinhead_placement(ex_seq, sm, nondup_idx[i]);
        },
        execution::max_concurrency(ex_tbb)
    );

    if (builder.ctl().stopcondition())
        return;

    for (auto &h : heads)
        if (h && h->is_valid()) {
            leafs.emplace_back(h->junction_point().cast<float>(), h->r_back_mm);
            h->id = long(leafs.size() - 1);
            builder.add_head(h->id, *h);
        }

    auto &its = *sm.emesh.get_triangle_mesh();
    ExPolygons bedpolys = {branchingtree::make_bed_poly(its)};

    auto props = branchingtree::Properties{}
                     .bed_shape(bedpolys)
                     .ground_level(sla::ground_level(sm))
                     .max_slope(sm.cfg.bridge_slope)
                     .max_branch_length(sm.cfg.max_bridge_length_mm);

    auto meshpts = sm.cfg.ground_facing_only ?
                       std::vector<branchingtree::Node>{} :
                       branchingtree::sample_mesh(its,
                                                  props.sampling_radius());

    auto bedpts  = branchingtree::sample_bed(props.bed_shape(),
                                             float(props.ground_level()),
                                             props.sampling_radius());

    branchingtree::PointCloud nodes{std::move(meshpts), std::move(bedpts),
                                    std::move(leafs), props};

    BranchingTreeBuilder vbuilder{builder, sm, nodes};
    branchingtree::build_tree(nodes, vbuilder);

    if constexpr (props.group_pillars()) {

        std::vector<branchingtree::Node> bedleafs;
        for (auto n : vbuilder.pillars()) {
            n.left =  branchingtree::Node::ID_NONE;
            n.right = branchingtree::Node::ID_NONE;
            bedleafs.emplace_back(n);
        }

        branchingtree::PointCloud gndnodes{{}, nodes.get_bedpoints(), bedleafs, props};
        BranchingTreeBuilder gndbuilder{builder, sm, gndnodes};
        branchingtree::build_tree(gndnodes, gndbuilder);

        build_pillars(builder, gndbuilder, sm);
    } else {
        build_pillars(builder, vbuilder, sm);
    }

    for (size_t id : vbuilder.unroutable_pinheads())
        builder.head(id).invalidate();

}

}} // namespace Slic3r::sla
