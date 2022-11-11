#ifndef ORGANICTREE_HPP
#define ORGANICTREE_HPP

#include <type_traits>
#include <utility>
#include <optional>

namespace Slic3r { namespace organictree {

enum class NodeType { Bed, Mesh, Junction };

template <class T> struct DomainTraits_ {
    using Node = typename T::Node;

    static void push(const T &dom, const Node &n)
    {
        dom.push_junction(n);
    }

    static Node pop(T &dom) { return dom.pop(); }

    static bool empty(const T &dom) { return dom.empty(); }

    static std::optional<std::pair<Node, NodeType>>
    closest(const T &dom, const Node &n)
    {
        return dom.closest(n);
    }

    static Node merge_node(const T &dom, const Node &a, const Node &b)
    {
        return dom.merge_node(a, b);
    }

    static void bridge(T &dom, const Node &from, const Node &to)
    {
        dom.bridge(from, to);
    }

    static void anchor(T &dom, const Node &from, const Node &to)
    {
        dom.anchor(from, to);
    }

    static void pillar(T &dom, const Node &from, const Node &to)
    {
        dom.pillar(from, to);
    }

    static void merge (T &dom, const Node &n1, const Node &n2, const Node &mrg)
    {
        dom.merge(n1, n2, mrg);
    }

    static void report_fail(T &dom, const Node &n) { dom.report_fail(n); }
};

template<class Domain>
void build_tree(Domain &&D)
{
    using Dom = DomainTraits_<std::remove_cv_t<std::remove_reference_t<Domain>>>;
    using Node = typename Dom::Node;

    while (! Dom::empty(D)) {
        Node n = Dom::pop(D);

        std::optional<std::pair<Node, NodeType>> C = Dom::closest(D, n);

        if (!C) {
            Dom::report_fail(D, n);
        } else switch (C->second) {
            case NodeType::Bed:
                Dom::pillar(D, n, C->first);
                break;
            case NodeType::Mesh:
                Dom::anchor(D, n, C->first);
                break;
            case NodeType::Junction: {
                Node M = Dom::merge_node(D, n, C->first);

                if (M == C->first) {
                    Dom::bridge(D, n, C->first);
                } else {
                    Dom::push(D, M);
                    Dom::merge(D, n, M, C->first);
                }
                break;
            }
        }
    }
}

}} // namespace Slic3r::organictree

#endif // ORGANICTREE_HPP
