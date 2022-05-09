#include "BranchingTree.hpp"
#include "PointCloud.hpp"

#include <numeric>
#include <optional>
#include <algorithm>

#include "libslic3r/SLA/SupportTreeUtils.hpp"

namespace Slic3r { namespace branchingtree {

bool build_tree(PointCloud &nodes, Builder &builder)
{
    auto ptsqueue = nodes.start_queue();
    auto &properties = nodes.properties();

    struct NodeDistance { size_t node_id; float distance; };
    auto distances = reserve_vector<NodeDistance>(nodes.reachable_count());

    while (!ptsqueue.empty()) {
        size_t node_id = ptsqueue.top();
        ptsqueue.pop();

        Node node = nodes.get(node_id);
        nodes.remove_node(node_id);

        distances.clear();
        distances.reserve(nodes.reachable_count());

        nodes.foreach_reachable(node.pos, [&distances](size_t id, float distance) {
            if (!std::isinf(distance))
                distances.emplace_back(NodeDistance{id, distance});
        });

        std::sort(distances.begin(), distances.end(),
                  [](auto &a, auto &b) { return a.distance < b.distance; });

        if (distances.empty()) {
            builder.report_unroutable(node);
            continue;
        }

        auto closest_it = distances.begin();
        bool routed = false;
        while (closest_it != distances.end() && !routed) {
            size_t closest_node_id = closest_it->node_id;
            Node closest_node = nodes.get(closest_node_id);

            auto type = nodes.get_type(closest_node_id);
            float w = nodes.get(node_id).weight + closest_it->distance;
            closest_node.Rmin = std::max(node.Rmin, closest_node.Rmin);

            switch (type) {
            case BED: {
                closest_node.weight = w;
                if ((routed = builder.add_ground_bridge(node, closest_node))) {
                    closest_node.left = closest_node.right = node_id;
                    nodes.get(closest_node_id) = closest_node;
                    nodes.remove_node(closest_node_id);
                }

                break;
            }
            case MESH: {
                closest_node.weight = w;
                if ((routed = builder.add_mesh_bridge(node, closest_node))) {
                    closest_node.left = closest_node.right = node_id;
                    nodes.get(closest_node_id) = closest_node;
                    nodes.remove_node(closest_node_id);
                }

                break;
            }
            case SUPP:
            case JUNCTION: {
                auto max_slope = float(properties.max_slope());

                if (auto mergept = find_merge_pt(node.pos, closest_node.pos, max_slope)) {

                    float mergedist_closest = (*mergept - closest_node.pos).norm();
                    float mergedist_node = (*mergept - node.pos).norm();
                    float Wnode = nodes.get(node_id).weight;
                    float Wclosest = nodes.get(closest_node_id).weight;
                    float Wsum = std::max(Wnode, Wclosest);
                    float distsum = std::max(mergedist_closest, mergedist_node);
                    w = Wsum + distsum;

                    if (mergedist_closest > EPSILON) {
                        Node mergenode{*mergept, closest_node.Rmin};
                        mergenode.weight = w;
                        mergenode.id = int(nodes.next_junction_id());

                        if ((routed = builder.add_merger(node, closest_node, mergenode))) {
                            mergenode.left = node_id;
                            mergenode.right = closest_node_id;
                            size_t new_idx = nodes.insert_junction(mergenode);
                            ptsqueue.push(new_idx);
                            ptsqueue.remove(nodes.get_queue_idx(closest_node_id));
                            nodes.remove_node(closest_node_id);
                        }
                    } else if (closest_node.left == Node::ID_NONE ||
                               closest_node.right == Node::ID_NONE)
                    {
                        closest_node.weight = w;
                        if ((routed = builder.add_bridge(node, closest_node))) {
                            if (closest_node.left == Node::ID_NONE)
                                closest_node.left = node_id;
                            else if (closest_node.right == Node::ID_NONE)
                                closest_node.right = node_id;

                            nodes.get(closest_node_id) = closest_node;
                        }
                    }
                }

                break;
            }
            case NONE:;
            }

            ++closest_it;
        }

        if (!routed)
            builder.report_unroutable(node);
    }

    return true;
}

bool build_tree(const indexed_triangle_set & its,
                const std::vector<Node> &support_roots,
                Builder &                    builder,
                const Properties &           properties)
{
    PointCloud nodes(its, support_roots, properties);

    return build_tree(nodes, builder);
}

}} // namespace Slic3r::branchingtree
