#ifndef PERFORMCSGMESHBOOLEANS_HPP
#define PERFORMCSGMESHBOOLEANS_HPP

#include "CSGMesh.hpp"

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

template<class It>
void perform_csgmesh_booleans(MeshBoolean::cgal::CGALMesh &cgalm,
                              const Range<It>             &csgparts)
{
    for (auto &csgpart : csgparts) {
        auto m = get_cgalmesh(csgpart);
        if (m) {
            switch (get_operation(csgpart)) {
            case CSGType::Union:
                MeshBoolean::cgal::plus(cgalm, *m);
                break;
            case CSGType::Difference:
                MeshBoolean::cgal::minus(cgalm, *m);
                break;
            case CSGType::Intersection:
                MeshBoolean::cgal::intersect(cgalm, *m);
                break;
            }
        }
    }
}

template<class It>
MeshBoolean::cgal::CGALMeshPtr perform_csgmesh_booleans(const Range<It> &csgparts)
{
    auto ret = MeshBoolean::cgal::triangle_mesh_to_cgal(indexed_triangle_set{});
    if (ret)
        perform_csgmesh_booleans(*ret, csgparts);

    return ret;
}

} // namespace csg
} // namespace Slic3r

#endif // PERFORMCSGMESHBOOLEANS_HPP
