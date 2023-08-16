#ifndef ARRANGEFIRSTFIT_HPP
#define ARRANGEFIRSTFIT_HPP

#include <iterator>

#include <libslic3r/Arrange/Core/ArrangeBase.hpp>

namespace Slic3r { namespace arr2 { namespace firstfit {

struct SelectionTag {};

// Can be specialized by Items
template<class ArrItem, class En = void>
struct ItemArrangedVisitor {
    template<class Bed, class PIt, class RIt>
    static void on_arranged(ArrItem &itm,
                            const Bed &bed,
                            const Range<PIt> &packed_items,
                            const Range<RIt> &remaining_items)
    {}
};

// Use the the visitor baked into the ArrItem type by default
struct DefaultOnArrangedFn {
    template<class ArrItem, class Bed, class PIt, class RIt>
    void operator()(ArrItem &itm,
                    const Bed &bed,
                    const Range<PIt> &packed,
                    const Range<RIt> &remaining)
    {
        ItemArrangedVisitor<StripCVRef<ArrItem>>::on_arranged(itm, bed, packed,
                                                              remaining);
    }
};

struct DefaultItemCompareFn {
    template<class ArrItem>
    bool operator() (const ArrItem &ia, const ArrItem &ib)
    {
        return get_priority(ia) > get_priority(ib);
    }
};

template<class CompareFn  = DefaultItemCompareFn,
         class OnArrangedFn = DefaultOnArrangedFn,
         class StopCondition = DefaultStopCondition>
struct SelectionStrategy
{
    CompareFn  cmpfn;
    OnArrangedFn on_arranged_fn;
    StopCondition cancel_fn;

    SelectionStrategy(CompareFn cmp = {},
             OnArrangedFn on_arranged = {},
             StopCondition stopcond = {})
        : cmpfn{cmp},
          on_arranged_fn{std::move(on_arranged)},
          cancel_fn{std::move(stopcond)}
    {}
};

} // namespace firstfit

template<class... Args> struct SelStrategyTag_<firstfit::SelectionStrategy<Args...>> {
    using Tag = firstfit::SelectionTag;
};

template<class It,
         class ConstIt,
         class TBed,
         class SelStrategy,
         class PackStrategy>
void arrange(
    SelStrategy &&sel,
    PackStrategy &&ps,
    const Range<It> &items,
    const Range<ConstIt> &fixed,
    const TBed &bed,
    const firstfit::SelectionTag &)
{
    using ArrItem = typename std::iterator_traits<It>::value_type;
    using ArrItemRef = std::reference_wrapper<ArrItem>;

    auto sorted_items = reserve_vector<ArrItemRef>(items.size());

    for (auto &itm : items) {
        set_bed_index(itm, Unarranged);
        sorted_items.emplace_back(itm);
    }

    int max_bed_idx = get_bed_count(fixed);

    using Context = PackStrategyContext<PackStrategy, ArrItem>;

    auto bed_contexts = reserve_vector<Context>(max_bed_idx + 1);

    for (auto &itm : fixed) {
        if (get_bed_index(itm) >= 0) {
            auto bedidx = static_cast<size_t>(get_bed_index(itm));

            while (bed_contexts.size() <= bedidx)
                bed_contexts.emplace_back(
                    create_context<ArrItem>(ps, bed, bedidx));

            add_fixed_item(bed_contexts[bedidx], itm);
        }
    }

    if constexpr (!std::is_null_pointer_v<decltype(sel.cmpfn)>) {
        std::stable_sort(sorted_items.begin(), sorted_items.end(), sel.cmpfn);
    }

    auto is_cancelled = [&sel]() {
        return sel.cancel_fn();
    };

    remove_unpackable_items(ps, sorted_items, bed, [&is_cancelled]() {
        return is_cancelled();
    });

    auto it = sorted_items.begin();

    using SConstIt = typename std::vector<ArrItemRef>::const_iterator;

    while (it != sorted_items.end() && !is_cancelled()) {
        bool was_packed = false;
        size_t j = 0;
        while (!was_packed && !is_cancelled()) {
            for (; j < bed_contexts.size() && !was_packed && !is_cancelled(); j++) {
                set_bed_index(*it, int(j));

                auto remaining = Range{std::next(static_cast<SConstIt>(it)),
                                       sorted_items.cend()};

                was_packed = pack(ps, bed, *it, bed_contexts[j], remaining);

                if(was_packed) {
                    add_packed_item(bed_contexts[j], *it);

                    auto packed_range = Range{sorted_items.cbegin(),
                                              static_cast<SConstIt>(it)};

                    sel.on_arranged_fn(*it, bed, packed_range, remaining);
                } else {
                    set_bed_index(*it, Unarranged);
                }
            }

            if (!was_packed) {
                bed_contexts.emplace_back(
                    create_context<ArrItem>(ps, bed, bed_contexts.size()));
                j = bed_contexts.size() - 1;
            }
        }
        ++it;
    }
}

}} // namespace Slic3r::arr2

#endif // ARRANGEFIRSTFIT_HPP
