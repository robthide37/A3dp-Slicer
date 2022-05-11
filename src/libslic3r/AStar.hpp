#ifndef ASTAR_HPP
#define ASTAR_HPP

#include "libslic3r/Point.hpp"
#include "libslic3r/MutablePriorityQueue.hpp"

#include <unordered_map>

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

// Run the AStar algorithm on a tracer implementation.
// The 'tracer' argument encapsulates the domain (grid, point cloud, etc...)
// The 'source' argument is the starting node.
// The 'out' argument is the output iterator into which the output nodes are
// written.
// Note that no destination node is given. The tracer's goal_heuristic() method
// should return a negative value if a node is a destination node.
template<class Tracer, class It>
bool search_route(const Tracer &tracer, const TracerNodeT<Tracer> &source, It out)
{
    using namespace detail;

    using Node = TracerNodeT<Tracer>;
    enum  class QueueType { Open, Closed, None };

    struct QNode // Queue node. Keeps track of scores g, and h
    {
        Node      node;                    // The actual node itself
        QueueType qtype = QueueType::None; // Which queue holds this node

        float  g = 0.f, h = 0.f;
        float f() const { return g + h; }
    };

       // TODO: apply a linear memory allocator
    using QMap = std::unordered_map<size_t, QNode>;

       // The traversed nodes are stored here encapsulated in QNodes
    QMap cached_nodes;

    struct LessPred { // Comparison functor needed by MutablePriorityQueue
        QMap &m;
        bool operator ()(size_t node_a, size_t node_b) {
            auto ait = m.find(node_a);
            auto bit = m.find(node_b);
            assert (ait != m.end() && bit != m.end());

            return ait->second.f() < bit->second.f();
        }
    };

    auto qopen =
        make_mutable_priority_queue<size_t, false>([](size_t, size_t){},
                                                   LessPred{cached_nodes});

    auto qclosed =
        make_mutable_priority_queue<size_t, false>([](size_t, size_t){},
                                                   LessPred{cached_nodes});

    QNode initial{source, QueueType::Open};
    cached_nodes.insert({unique_id(tracer, source), initial});
    qopen.push(unique_id(tracer, source));

    bool goal_reached = false;

    while (!goal_reached && !qopen.empty()) {
        size_t q_id = qopen.top();
        qopen.pop();
        QNode q = cached_nodes.at(q_id);

        foreach_reachable(tracer, q.node, [&](const Node &nd) {
            if (goal_reached) return goal_reached;

            float h = goal_heuristic(tracer, nd);
            if (h < 0.f) {
                goal_reached = true;
            } else {
                float  dst = trace_distance(tracer, q.node, nd);
                QNode  qnd{nd, QueueType::None, q.g + dst, h};
                size_t qnd_id = unique_id(tracer, nd);

                auto it = cached_nodes.find(qnd_id);

                if (it == cached_nodes.end() ||
                    (it->second.qtype != QueueType::None && qnd.f() < it->second.f())) {
                    qnd.qtype = QueueType::Open;
                    cached_nodes.insert_or_assign(qnd_id, qnd);
                    qopen.push(qnd_id);
                }
            }

            return goal_reached;
        });

        q.qtype = QueueType::Closed;
        cached_nodes.insert_or_assign(q_id, q);
        qclosed.push(q_id);

           // write the output
        *out = q.node;
        ++out;
    }

    return goal_reached;
}

}} // namespace Slic3r::astar

#endif // ASTAR_HPP
