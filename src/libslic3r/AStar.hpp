#ifndef ASTAR_HPP
#define ASTAR_HPP

#include <unordered_map>

#include "libslic3r/Point.hpp"
#include "libslic3r/MutablePriorityQueue.hpp"

namespace Slic3r { namespace astar {

// Input interface for the Astar algorithm. Specialize this struct for a
// particular type and implement all the 4 methods and specify the Node type
// to register the new type for the astar implementation.
template<class T> struct TracerTraits_
{
    // The type of a node used by this tracer. Usually a point in space.
    using Node = typename T::Node;

    // Call fn for every new node reachable from node 'src'. fn should have the
    // candidate node as its only argument.
    template<class Fn>
    static void foreach_reachable(const T &tracer, const Node &src, Fn &&fn)
    {
        tracer.foreach_reachable(src, fn);
    }

    // Get the distance from node 'a' to node 'b'. This is sometimes referred
    // to as the g value of a node in AStar context.
    static float distance(const T &tracer, const Node &a, const Node &b)
    {
        return tracer.distance(a, b);
    }

    // Get the estimated distance heuristic from node 'n' to the destination.
    // This is referred to as the h value in AStar context.
    // If node 'n' is the goal, this function should return a negative value.
    // Note that this heuristic should be admissible (never bigger than the real
    // cost) in order for Astar to work.
    static float goal_heuristic(const T &tracer, const Node &n)
    {
        return tracer.goal_heuristic(n);
    }

    // Return a unique identifier (hash) for node 'n'.
    static size_t unique_id(const T &tracer, const Node &n)
    {
        return tracer.unique_id(n);
    }
};

// Helper definition to get the node type of a tracer
template<class T>
using TracerNodeT = typename TracerTraits_<remove_cvref_t<T>>::Node;

namespace detail {
// Helper functions dispatching calls through the TracerTraits_ interface

template<class T> using TracerTraits = TracerTraits_<remove_cvref_t<T>>;

template<class T, class Fn>
void foreach_reachable(const T &tracer, const TracerNodeT<T> &from, Fn &&fn)
{
    TracerTraits<T>::foreach_reachable(tracer, from, fn);
}

template<class T>
float trace_distance(const T &tracer, const TracerNodeT<T> &a, const TracerNodeT<T> &b)
{
    return TracerTraits<T>::distance(tracer, a, b);
}

template<class T>
float goal_heuristic(const T &tracer, const TracerNodeT<T> &n)
{
    return TracerTraits<T>::goal_heuristic(tracer, n);
}

template<class T>
size_t unique_id(const T &tracer, const TracerNodeT<T> &n)
{
    return TracerTraits<T>::unique_id(tracer, n);
}

} // namespace astar_detail

constexpr size_t UNASSIGNED = size_t(-1);

template<class Tracer>
struct QNode // Queue node. Keeps track of scores g, and h
{
    TracerNodeT<Tracer>  node;      // The actual node itself
    size_t queue_id;  // Position in the open queue or UNASSIGNED if closed
    size_t parent;    // unique id of the parent or UNASSIGNED

    float  g, h;
    float f() const { return g + h; }

    QNode(TracerNodeT<Tracer> n = {},
          size_t              p    = UNASSIGNED,
          float               gval = std::numeric_limits<float>::infinity(),
          float               hval = 0.f)
        : node{std::move(n)}, parent{p}, queue_id{UNASSIGNED}, g{gval}, h{hval}
    {}
};

// Run the AStar algorithm on a tracer implementation.
// The 'tracer' argument encapsulates the domain (grid, point cloud, etc...)
// The 'source' argument is the starting node.
// The 'out' argument is the output iterator into which the output nodes are
// written. For performance reasons, the order is reverse, from the destination
// to the source -- (destination included, source is not).
// The 'cached_nodes' argument is an optional associative container to hold a
// QNode entry for each visited node. Any compatible container can be used
// (like std::map or maps with different allocators, even a sufficiently large
// std::vector).
//
// Note that no destination node is given in the signature. The tracer's
// goal_heuristic() method should return a negative value if a node is a
// destination node.
template<class Tracer,
         class It,
         class NodeMap = std::unordered_map<size_t, QNode<Tracer>>>
bool search_route(const Tracer              &tracer,
                  const TracerNodeT<Tracer> &source,
                  It                         out,
                  NodeMap                  &&cached_nodes = {})
{
    using namespace detail;

    using Node  = TracerNodeT<Tracer>;
    using QNode = QNode<Tracer>;

    struct LessPred { // Comparison functor needed by the priority queue
        NodeMap &m;
        bool operator ()(size_t node_a, size_t node_b) {
            return m[node_a].f() < m[node_b].f();
        }
    };

    auto qopen = make_mutable_priority_queue<size_t, true>(
        [&cached_nodes](size_t el, size_t qidx) {
            cached_nodes[el].queue_id = qidx;
        },
        LessPred{cached_nodes});

    QNode initial{source, /*parent = */ UNASSIGNED, /*g = */0.f};
    size_t source_id = unique_id(tracer, source);
    cached_nodes[source_id] = initial;
    qopen.push(source_id);

    size_t goal_id = goal_heuristic(tracer, source) < 0.f ? source_id :
                                                            UNASSIGNED;

    while (goal_id == UNASSIGNED && !qopen.empty()) {
        size_t q_id = qopen.top();
        qopen.pop();
        QNode &q = cached_nodes[q_id];

        // This should absolutely be initialized in the cache already
        assert(!std::isinf(q.g));

        foreach_reachable(tracer, q.node, [&](const Node &succ_nd) {
            if (goal_id != UNASSIGNED)
                return true;

            float  h       = goal_heuristic(tracer, succ_nd);
            float  dst     = trace_distance(tracer, q.node, succ_nd);
            size_t succ_id = unique_id(tracer, succ_nd);
            QNode  qsucc_nd{succ_nd, q_id, q.g + dst, h};

            if (h < 0.f) {
                goal_id = succ_id;
                cached_nodes[succ_id] = qsucc_nd;
            } else {
                // If succ_id is not in cache, it gets created with g = infinity
                QNode &prev_nd = cached_nodes[succ_id];

                if (qsucc_nd.g < prev_nd.g) {
                    // new route is better, apply it:

                    // Save the old queue id, it would be lost after the next line
                    size_t queue_id = prev_nd.queue_id;

                    // The cache needs to be updated either way
                    prev_nd = qsucc_nd;

                    if (queue_id == UNASSIGNED)
                        // was in closed or unqueued, rescheduling
                        qopen.push(succ_id);
                    else // was in open, updating
                        qopen.update(queue_id);
                }
            }

            return goal_id != UNASSIGNED;
        });
    }

    // Write the output, do not reverse. Clients can do so if they need to.
    if (goal_id != UNASSIGNED) {
        const QNode *q = &cached_nodes[goal_id];
        while (!std::isinf(q->g) && q->parent != UNASSIGNED) {
            *out = q->node;
            ++out;
            q = &cached_nodes[q->parent];
        }

        if (std::isinf(q->g)) // Something went wrong
            goal_id = UNASSIGNED;
    }

    return goal_id != UNASSIGNED;
}

}} // namespace Slic3r::astar

#endif // ASTAR_HPP
