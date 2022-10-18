#ifndef PERFORMCSGMESHBOOLEANS_HPP
#define PERFORMCSGMESHBOOLEANS_HPP

#include "CSGMesh.hpp"

#include "libslic3r/MeshBoolean.hpp"

namespace Slic3r {

template<class It>
void perform_csgmesh_booleans(MeshBoolean::cgal::CGALMesh &cgalm,
                              const Range<It>             &csg)
{

}

} // namespace Slic3r

#endif // PERFORMCSGMESHBOOLEANS_HPP
