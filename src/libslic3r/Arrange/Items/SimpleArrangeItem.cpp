#include "SimpleArrangeItem.hpp"
#include "libslic3r/Arrange/ArrangeImpl.hpp"

namespace Slic3r { namespace arr2 {

Polygon SimpleArrangeItem::outline() const
{
    Polygon ret = shape();
    ret.rotate(m_rotation);
    ret.translate(m_translation);

    return ret;
}

}} // namespace Slic3r::arr2
