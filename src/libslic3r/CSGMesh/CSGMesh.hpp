#ifndef CSGMESH_HPP
#define CSGMESH_HPP

#include <libslic3r/MTUtils.hpp> // for AnyPtr
#include <admesh/stl.h>

namespace Slic3r { namespace csg {

// A CSGPartT should be an object that can provide at least a mesh + trafo and an
// associated csg operation. A collection of CSGPartT objects can then
// be interpreted as one model and used in various contexts. It can be assembled
// with CGAL or OpenVDB, rendered with OpenCSG or provided to a ray-tracer to
// deal with various parts of it according to the supported CSG types...
//
// A few simple templated interface functions are provided here and a default
// CSGPart class that implements the necessary means to be usable as a
// CSGPartT object.

// Supported CSG operation types
enum class CSGType { Union, Difference, Intersection };
enum class CSGStackOp { Push, Continue, Pop };

// Get the CSG operation of the part. Can be overriden for any type
template<class CSGPartT> CSGType get_operation(const CSGPartT &part)
{
    return part.operation;
}

template<class CSGPartT> CSGStackOp get_stack_operation(const CSGPartT &part)
{
    return part.stack_operation;
}

// Get the mesh for the part. Can be overriden for any type
template<class CSGPartT>
const indexed_triangle_set *get_mesh(const CSGPartT &part)
{
    return part.its_ptr.get();
}

// Get the transformation associated with the mesh inside a CSGPartT object.
// Can be overriden for any type.
template<class CSGPartT>
Transform3f get_transform(const CSGPartT &part)
{
    return part.trafo;
}

// Provide default overloads for indexed_triangle_set to be usable as a plain
// CSGPart with an implicit union operation

inline CSGType get_operation(const indexed_triangle_set &part)
{
    return CSGType::Union;
}

inline CSGStackOp get_stack_operation(const indexed_triangle_set &part)
{
    return CSGStackOp::Continue;
}

inline const indexed_triangle_set * get_mesh(const indexed_triangle_set &part)
{
    return &part;
}

inline Transform3f get_transform(const indexed_triangle_set &part)
{
    return Transform3f::Identity();
}

// Default implementation
struct CSGPart {
    AnyPtr<const indexed_triangle_set> its_ptr;
    Transform3f trafo;
    CSGType operation;
    CSGStackOp stack_operation;

    CSGPart(AnyPtr<const indexed_triangle_set> ptr = {},
            CSGType                            op  = CSGType::Union,
            const Transform3f                 &tr  = Transform3f::Identity())
        : its_ptr{std::move(ptr)}
        , operation{op}
        , stack_operation{CSGStackOp::Continue}
        , trafo{tr}
    {}
};

}} // namespace Slic3r::csg

#endif // CSGMESH_HPP
