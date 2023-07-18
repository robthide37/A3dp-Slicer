#ifndef OUTLINECACHING_HPP
#define OUTLINECACHING_HPP

#include <any>
#include <unordered_map>

#include "libslic3r/ObjectID.hpp"
#include "libslic3r/ExPolygon.hpp"

namespace Slic3r { namespace arr2 {

struct CacheEntryConvex { Polygon outline; std::any context; };
struct CacheEntryFull
{
    ExPolygons outline;
    std::any   context;

    CacheEntryFull() = default;
    CacheEntryFull(ExPolygons outl, std::any ctx)
        : outline{std::move(outl)}, context{std::move(ctx)}
    {}
};

class OutlineCache
{
    std::unordered_map<size_t, CacheEntryFull>   m_cache_full;
    std::unordered_map<size_t, CacheEntryConvex> m_cache_convex;

public:
    const CacheEntryFull * full_outline(const ObjectID &id);

    void set_full_outline(const ObjectID &id, ExPolygons outline, std::any ctx);
};

}} // namespace Slic3r::arr2

#endif // OUTLINECACHING_HPP
