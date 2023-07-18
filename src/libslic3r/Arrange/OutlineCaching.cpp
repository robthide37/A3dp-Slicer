#include "OutlineCaching.hpp"

namespace Slic3r { namespace arr2 {

const CacheEntryFull *OutlineCache::full_outline(const ObjectID &id)
{
    const CacheEntryFull *ret = nullptr;

    auto it = m_cache_full.find(id.id);
    if (it != m_cache_full.end())
        ret = &(it->second);

    return ret;
}

void OutlineCache::set_full_outline(const ObjectID &id, ExPolygons outline, std::any ctx)
{
    m_cache_full[id.id] = CacheEntryFull{std::move(outline), std::move(ctx)};
}

}} // namespace Slic3r::arr2
