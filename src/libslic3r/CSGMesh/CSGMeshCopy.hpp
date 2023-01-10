#ifndef CSGMESHCOPY_HPP
#define CSGMESHCOPY_HPP

#include "CSGMesh.hpp"

namespace Slic3r { namespace csg {

// Copy a csg range but for the meshes, only copy the pointers.
template<class It, class OutIt>
void copy_csgrange_shallow(const Range<It> &csgrange, OutIt out)
{
    for (const auto &part : csgrange) {
        CSGPart cpy{AnyPtr<const indexed_triangle_set>{get_mesh(part)},
                    get_operation(part),
                    get_transform(part)};

        cpy.stack_operation = get_stack_operation(part);

        *out = std::move(cpy);
        ++out;
    }
}

// Copy the csg range, allocating new meshes
template<class It, class OutIt>
void copy_csgrange_deep(const Range<It> &csgrange, OutIt out)
{
    for (const auto &part : csgrange) {

        CSGPart cpy{{}, get_operation(part), get_transform(part)};

        if (auto meshptr = get_mesh(part)) {
            cpy.its_ptr = std::make_unique<const indexed_triangle_set>(*meshptr);
        }

        cpy.stack_operation = get_stack_operation(part);

        *out = std::move(cpy);
        ++out;
    }
}

template<class ItA, class ItB>
bool is_same(const Range<ItA> &A, const Range<ItB> &B)
{
    bool ret = true;

    size_t s = A.size();

    if (B.size() != s)
        ret = false;

    size_t i = 0;
    auto itA = A.begin();
    auto itB = B.begin();
    for (; ret && i < s; ++itA, ++itB, ++i) {
        ret = ret &&
              get_mesh(*itA) == get_mesh(*itB) &&
              get_operation(*itA) == get_operation(*itB) &&
              get_stack_operation(*itA) == get_stack_operation(*itB) &&
              get_transform(*itA).isApprox(get_transform(*itB));
    }

    return ret;
}

}} // namespace Slic3r::csg

#endif // CSGCOPY_HPP
