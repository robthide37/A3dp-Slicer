#ifndef ARR2_SCENE_HPP
#define ARR2_SCENE_HPP

#include <any>
#include <string_view>

#include "libslic3r/ObjectID.hpp"
#include "libslic3r/AnyPtr.hpp"
#include "libslic3r/Arrange/ArrangeSettingsView.hpp"
#include "libslic3r/Arrange/SegmentedRectangleBed.hpp"

namespace Slic3r { namespace arr2 {

class AnyWritable
{
public:
    virtual ~AnyWritable() = default;

    virtual void write(std::string_view key, std::any d) = 0;
};

class Arrangeable
{
public:
    virtual ~Arrangeable() = default;

    virtual ObjectID   id() const             = 0;
    virtual ObjectID   geometry_id() const    = 0;
    virtual ExPolygons full_outline() const   = 0;
    virtual Polygon    convex_outline() const = 0;

    virtual ExPolygons full_envelope() const { return {}; }
    virtual Polygon    convex_envelope() const { return {}; }

    virtual void transform(const Vec2d &transl, double rot) = 0;

    virtual bool is_printable() const { return true; }
    virtual bool is_selected() const { return true; }
    virtual int  priority() const { return 0; }

    virtual void imbue_data(AnyWritable &datastore) const {}
    void imbue_data(AnyWritable &&datastore) const { imbue_data(datastore); }

    // Returns the bed index on which the given ModelInstance is sitting.
    virtual int get_bed_index() const = 0;

    // Assign the ModelInstance to the given bed index. Note that this
    // method can return false, indicating that the given bed is not available
    // to be occupied (e.g. the handler has a limited amount of logical bed)
    virtual bool assign_bed(int bed_idx) = 0;
};

class ArrangeableModel
{
public:
    virtual ~ArrangeableModel() = default;

    virtual void for_each_arrangeable(std::function<void(Arrangeable &)>) = 0;
    virtual void for_each_arrangeable(std::function<void(const Arrangeable&)>) const = 0;

    virtual void visit_arrangeable(const ObjectID &id, std::function<void(const Arrangeable &)>) const = 0;
    virtual void visit_arrangeable(const ObjectID &id, std::function<void(Arrangeable &)>) = 0;

    virtual ObjectID add_arrangeable(const ObjectID &prototype_id) = 0;

    size_t arrangeable_count() const
    {
        size_t cnt = 0;
        for_each_arrangeable([&cnt](auto &) { ++cnt; });

        return cnt;
    }
};

using XLBed = SegmentedRectangleBed<std::integral_constant<size_t, 4>,
                                    std::integral_constant<size_t, 4>>;

template<class... Args> struct ExtendedBed_
{
    using Type =
        boost::variant<XLBed, /* insert other types if needed*/ Args...>;
};

template<class... Args> struct ExtendedBed_<boost::variant<Args...>>
{
    using Type = boost::variant<XLBed, Args...>;
};

using ExtendedBed = typename ExtendedBed_<ArrangeBed>::Type;

template<class BedFn> void visit_bed(BedFn &&fn, const ExtendedBed &bed)
{
    boost::apply_visitor(fn, bed);
}

template<class BedFn> void visit_bed(BedFn &&fn, ExtendedBed &bed)
{
    boost::apply_visitor(fn, bed);
}

inline ExtendedBed to_extended_bed(const ArrangeBed &bed)
{
    ExtendedBed ret;
    boost::apply_visitor([&ret](auto &rawbed) { ret = rawbed; }, bed);

    return ret;
}

class Scene;

// A little CRTP to implement fluent interface returning Subclass references
template<class Subclass>
class SceneBuilderBase
{
protected:
    AnyPtr<ArrangeableModel> m_arrangeable_model;

    AnyPtr<const ArrangeSettingsView> m_settings;

    ExtendedBed m_bed = arr2::InfiniteBed{};

    coord_t m_brims_offs = 0;
    coord_t m_skirt_offs = 0;

public:

    virtual ~SceneBuilderBase() = default;

    SceneBuilderBase() = default;
    SceneBuilderBase(const SceneBuilderBase &) = delete;
    SceneBuilderBase& operator=(const SceneBuilderBase &) = delete;
    SceneBuilderBase(SceneBuilderBase &&) = default;
    SceneBuilderBase& operator=(SceneBuilderBase &&) = default;

    // All setters return an rvalue reference so that at the end, the
    // build_scene method can be called fluently

    Subclass &&set_arrange_settings(AnyPtr<const ArrangeSettingsView> settings)
    {
        m_settings = std::move(settings);
        return std::move(static_cast<Subclass&>(*this));
    }

    Subclass &&set_arrange_settings(const ArrangeSettingsView &settings)
    {
        m_settings = std::make_unique<ArrangeSettings>(settings);
        return std::move(static_cast<Subclass&>(*this));
    }

    Subclass &&set_bed(const Points &pts)
    {
        m_bed = arr2::to_arrange_bed(pts);
        return std::move(static_cast<Subclass&>(*this));
    }

    Subclass && set_bed(const arr2::ArrangeBed &bed)
    {
        m_bed = bed;
        return std::move(static_cast<Subclass&>(*this));
    }

    Subclass &&set_bed(const XLBed &bed)
    {
        m_bed = bed;
        return std::move(static_cast<Subclass&>(*this));
    }

    Subclass &&set_arrangeable_model(AnyPtr<ArrangeableModel> model)
    {
        m_arrangeable_model = std::move(model);
        return std::move(static_cast<Subclass&>(*this));
    }

    // Can only be called on an rvalue instance (hence the && at the end),
    // the method will potentially move its content into sc
    virtual void build_scene(Scene &sc) &&;
};

class BasicSceneBuilder: public SceneBuilderBase<BasicSceneBuilder> {};

class Scene
{
    template <class Sub> friend class SceneBuilderBase;

    AnyPtr<ArrangeableModel>                m_amodel;
    AnyPtr<const ArrangeSettingsView>       m_settings;
    ExtendedBed m_bed;

public:
    // Can only be built from an rvalue SceneBuilder, as it's content will
    // potentially be moved to the constructed ArrangeScene object
    template<class Sub>
    explicit Scene(SceneBuilderBase<Sub> &&bld)
    {
        std::move(bld).build_scene(*this);
    }

    const ArrangeableModel &model() const noexcept { return *m_amodel; }
    ArrangeableModel       &model() noexcept { return *m_amodel; }

    const ArrangeSettingsView &settings() const noexcept { return *m_settings; }

    template<class BedFn> void visit_bed(BedFn &&fn) const
    {
        arr2::visit_bed(fn, m_bed);
    }

    const ExtendedBed & bed() const { return m_bed; }

    std::vector<ObjectID> selected_ids() const;
};

std::set<ObjectID> selected_geometry_ids(const Scene &sc);

class EmptyArrangeableModel: public ArrangeableModel
{
public:
    void for_each_arrangeable(std::function<void(Arrangeable &)>) override {}
    void for_each_arrangeable(std::function<void(const Arrangeable&)>) const override {}
    void visit_arrangeable(const ObjectID &id, std::function<void(const Arrangeable &)>) const override {}
    void visit_arrangeable(const ObjectID &id, std::function<void(Arrangeable &)>) override {}
    ObjectID add_arrangeable(const ObjectID &prototype_id) override { return {}; }
};

template<class Subclass>
void SceneBuilderBase<Subclass>::build_scene(Scene &sc) &&
{
    if (!m_arrangeable_model)
        m_arrangeable_model = std::make_unique<EmptyArrangeableModel>();

    if (!m_settings)
        m_settings = std::make_unique<arr2::ArrangeSettings>();

    coord_t inset = std::max(scaled(m_settings->get_distance_from_bed()),
                             m_skirt_offs + m_brims_offs);

    coord_t md = scaled(m_settings->get_distance_from_objects());
    md = md / 2 - inset;

    visit_bed([md](auto &rawbed) { rawbed = offset(rawbed, md); }, m_bed);

    sc.m_settings = std::move(m_settings);
    sc.m_amodel = std::move(m_arrangeable_model);
    sc.m_bed = std::move(m_bed);
}

class ArrangeResult
{
public:
    virtual ~ArrangeResult() = default;

    virtual bool apply_on(ArrangeableModel &mdlwt) = 0;
};

enum class Tasks { Arrange, FillBed };

class ArrangeTaskCtl
{
public:
    virtual ~ArrangeTaskCtl() = default;

    virtual void update_status(int st) = 0;

    virtual bool was_canceled() const = 0;
};

class DummyCtl : public ArrangeTaskCtl
{
public:
    void update_status(int) override {}
    bool was_canceled() const override { return false; }
};

class ArrangeTaskBase
{
public:
    using Ctl = ArrangeTaskCtl;

    virtual ~ArrangeTaskBase() = default;

    [[nodiscard]] virtual std::unique_ptr<ArrangeResult> process(Ctl &ctl) = 0;

    [[nodiscard]] virtual int item_count_to_process() const = 0;

    [[nodiscard]] static std::unique_ptr<ArrangeTaskBase> create(
        Tasks task_type, const Scene &sc);

    [[nodiscard]] std::unique_ptr<ArrangeResult> process(Ctl &&ctl)
    {
        return process(ctl);
    }

    [[nodiscard]] std::unique_ptr<ArrangeResult> process()
    {
        return process(DummyCtl{});
    }
};

void arrange(Scene &scene, ArrangeTaskCtl &ctl);
inline void arrange(Scene &scene, ArrangeTaskCtl &&ctl = DummyCtl{})
{
    arrange(scene, ctl);
}

} // namespace arr2
} // namespace Slic3r

#endif // ARR2_SCENE_HPP
