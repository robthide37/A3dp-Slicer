#ifndef WRITABLEITEMTRAITS_HPP
#define WRITABLEITEMTRAITS_HPP

#include "libslic3r/Arrange/Core/ArrangeItemTraits.hpp"
#include "libslic3r/Arrange/Core/DataStoreTraits.hpp"

#include "libslic3r/ExPolygon.hpp"

namespace Slic3r { namespace arr2 {

template<class Itm> struct IsWritableItem_ : public std::false_type
{};

// Using this interface to set up any arrange item. Provides default
// implementation but it needs to be explicitly switched on with
// IsWritableItem_ or completely reimplement a specialization.
template<class Itm, class En = void> struct WritableItemTraits_
{
    static_assert(IsWritableItem_<Itm>::value, "Not a Writable item type!");

    static void set_priority(Itm &itm, int p) { itm.set_priority(p); }

    static void set_convex_shape(Itm &itm, const Polygon &shape)
    {
        itm.set_convex_shape(shape);
    }

    static void set_shape(Itm &itm, const ExPolygons &shape)
    {
        itm.set_shape(shape);
    }

    static void set_convex_envelope(Itm &itm, const Polygon &envelope)
    {
        itm.set_convex_envelope(envelope);
    }

    static void set_envelope(Itm &itm, const ExPolygons &envelope)
    {
        itm.set_envelope(envelope);
    }

    template<class T>
    static void set_arbitrary_data(Itm &itm, const std::string &key, T &&data)
    {
        if constexpr (IsWritableDataStore<Itm>)
            set_data(itm, key, std::forward<T>(data));
    }

    static void set_allowed_rotations(Itm                       &itm,
                                      const std::vector<double> &rotations)
    {
        itm.set_allowed_rotations(rotations);
    }
};

template<class T>
using WritableItemTraits = WritableItemTraits_<StripCVRef<T>>;

template<class T> constexpr bool IsWritableItem = IsWritableItem_<T>::value;
template<class T, class TT = T>
using WritableItemOnly = std::enable_if_t<IsWritableItem<T>, TT>;

template<class Itm> void set_priority(Itm &itm, int p)
{
    WritableItemTraits<Itm>::set_priority(itm, p);
}

template<class Itm> void set_convex_shape(Itm &itm, const Polygon &shape)
{
    WritableItemTraits<Itm>::set_convex_shape(itm, shape);
}

template<class Itm> void set_shape(Itm &itm, const ExPolygons &shape)
{
    WritableItemTraits<Itm>::set_shape(itm, shape);
}

template<class Itm>
void set_convex_envelope(Itm &itm, const Polygon &envelope)
{
    WritableItemTraits<Itm>::set_convex_envelope(itm, envelope);
}

template<class Itm> void set_envelope(Itm &itm, const ExPolygons &envelope)
{
    WritableItemTraits<Itm>::set_envelope(itm, envelope);
}

template<class T, class Itm>
void set_arbitrary_data(Itm &itm, const std::string &key, T &&data)
{
    WritableItemTraits<Itm>::set_arbitrary_data(itm, key, std::forward<T>(data));
}

template<class Itm>
void set_allowed_rotations(Itm &itm, const std::vector<double> &rotations)
{
    WritableItemTraits<Itm>::set_allowed_rotations(itm, rotations);
}

template<class ArrItem> int raise_priority(ArrItem &itm)
{
    int ret = get_priority(itm) + 1;
    set_priority(itm, ret);

    return ret;
}

template<class ArrItem> int reduce_priority(ArrItem &itm)
{
    int ret = get_priority(itm) - 1;
    set_priority(itm, ret);

    return ret;
}

template<class It> int lowest_priority(const Range<It> &item_range)
{
    auto minp_it = std::min_element(item_range.begin(),
                                    item_range.end(),
                                    [](auto &itm1, auto &itm2) {
                                        return get_priority(itm1) <
                                               get_priority(itm2);
                                    });

    int min_priority = 0;
    if (minp_it != item_range.end())
        min_priority = get_priority(*minp_it);

    return min_priority;
}

}} // namespace Slic3r::arr2

#endif // WRITABLEITEMTRAITS_HPP
