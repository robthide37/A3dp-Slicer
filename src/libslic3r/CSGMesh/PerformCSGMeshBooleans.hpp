#ifndef PERFORMCSGMESHBOOLEANS_HPP
#define PERFORMCSGMESHBOOLEANS_HPP

#include <stack>
#include <vector>

#include "CSGMesh.hpp"

#include "libslic3r/Execution/ExecutionTBB.hpp"
#include "libslic3r/MeshBoolean.hpp"

namespace Slic3r { namespace csg {

// This method can be overriden when a specific CSGPart type supports caching
// of the voxel grid
template<class CSGPartT>
MeshBoolean::cgal::CGALMeshPtr get_cgalmesh(const CSGPartT &csgpart)
{
    const indexed_triangle_set *its = csg::get_mesh(csgpart);
    MeshBoolean::cgal::CGALMeshPtr ret;

    if (its) {
        indexed_triangle_set m = *its;
        its_transform(m, get_transform(csgpart));
        ret = MeshBoolean::cgal::triangle_mesh_to_cgal(m);
    }

    return ret;
}

namespace detail_cgal {

using MeshBoolean::cgal::CGALMeshPtr;

inline void perform_csg(CSGType op, CGALMeshPtr &dst, CGALMeshPtr &src)
{
    switch (op) {
    case CSGType::Union:
        MeshBoolean::cgal::plus(*dst, *src);
        break;
    case CSGType::Difference:
        MeshBoolean::cgal::minus(*dst, *src);
        break;
    case CSGType::Intersection:
        MeshBoolean::cgal::intersect(*dst, *src);
        break;
    }
}

template<class Ex, class It>
std::vector<CGALMeshPtr> get_cgalptrs(Ex policy, const Range<It> &csgrange)
{
    std::vector<CGALMeshPtr> ret(csgrange.size());
    execution::for_each(policy, size_t(0), csgrange.size(),
                        [&csgrange, &ret](size_t i) {
        auto it = csgrange.begin();
        std::advance(it, i);
        auto &csgpart = *it;
        ret[i]        = get_cgalmesh(csgpart);
    });

    return ret;
}

} // namespace detail

template<class It>
void perform_csgmesh_booleans(MeshBoolean::cgal::CGALMeshPtr &cgalm,
                              const Range<It>                &csgrange)
{
    using MeshBoolean::cgal::CGALMesh;
    using MeshBoolean::cgal::CGALMeshPtr;
    using namespace detail_cgal;

    struct Frame { CSGType op = CSGType::Union; CGALMeshPtr cgalptr; };

    std::stack opstack{std::vector<Frame>{}};

    if (!csgrange.empty() &&
        csg::get_stack_operation(*csgrange.begin()) != CSGStackOp::Push)
        opstack.push({});

    std::vector<CGALMeshPtr> cgalmeshes = get_cgalptrs(ex_tbb, csgrange);

    size_t csgidx = 0;
    for (auto &csgpart : csgrange) {

        auto op = get_operation(csgpart);
        CGALMeshPtr &cgalptr = cgalmeshes[csgidx++];

        if (get_stack_operation(csgpart) == CSGStackOp::Push) {
            opstack.push({op, nullptr});
            op = CSGType::Union;
        }

        Frame *top = &opstack.top();

        if (!top->cgalptr && op == CSGType::Union) {
            top->cgalptr = std::move(cgalptr);
        } else if (top->cgalptr && cgalptr) {
            perform_csg(get_operation(csgpart), top->cgalptr, cgalptr);
        }

        if (get_stack_operation(csgpart) == CSGStackOp::Pop) {
            CGALMeshPtr src = std::move(top->cgalptr);
            auto popop = opstack.top().op;
            opstack.pop();
            CGALMeshPtr &dst = opstack.top().cgalptr;
            perform_csg(popop, dst, src);
        }
    }

    cgalm = std::move(opstack.top().cgalptr);
}

template<class It, class Visitor>
It check_csgmesh_booleans(const Range<It> &csgrange, Visitor &&vfn)
{
    using namespace detail_cgal;

    std::vector<CGALMeshPtr> cgalmeshes(csgrange.size());
    auto check_part = [&csgrange, &cgalmeshes](size_t i)
    {
        auto it = csgrange.begin();
        std::advance(it, i);
        auto &csgpart = *it;
        auto m = get_cgalmesh(csgpart);

        if (!m || MeshBoolean::cgal::empty(*m))
            return;

        if (!MeshBoolean::cgal::does_bound_a_volume(*m) ||
            MeshBoolean::cgal::does_self_intersect(*m))
            return;

        cgalmeshes[i] = std::move(m);
    };
    execution::for_each(ex_tbb, size_t(0), csgrange.size(), check_part);

    It ret = csgrange.end();
    for (size_t i = 0; i < csgrange.size(); ++i) {
        if (!cgalmeshes[i]) {
            auto it = csgrange.begin();
            std::advance(it, i);
            vfn(it);

            if (it == csgrange.end())
                ret = it;
        }
    }

    return csgrange.end();
}

template<class It>
It check_csgmesh_booleans(const Range<It> &csgrange)
{
    return check_csgmesh_booleans(csgrange, [](auto &) {});
}

template<class It>
MeshBoolean::cgal::CGALMeshPtr perform_csgmesh_booleans(const Range<It> &csgparts)
{
    auto ret = MeshBoolean::cgal::triangle_mesh_to_cgal(indexed_triangle_set{});
    if (ret)
        perform_csgmesh_booleans(ret, csgparts);

    return ret;
}

} // namespace csg
} // namespace Slic3r

#endif // PERFORMCSGMESHBOOLEANS_HPP
