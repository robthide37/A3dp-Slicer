#ifndef ARRANGESETTINGSDB_APPCFG_HPP
#define ARRANGESETTINGSDB_APPCFG_HPP

#include "ArrangeSettingsView.hpp"
#include "libslic3r/AppConfig.hpp"
#include "libslic3r/PrintConfig.hpp"

namespace Slic3r {

class ArrangeSettingsDb_AppCfg: public arr2::ArrangeSettingsDb
{
    AppConfig *m_appcfg;
    std::function<const DynamicPrintConfig*(void)> m_config_getter;
    std::function<PrinterTechnology(void)> m_printtech_getter;

    struct Slot { Values vals; Values defaults; std::string postfix; };

    // Settings and their defaults are stored separately for fff,
    // sla and fff sequential mode
    Slot m_settings_fff, m_settings_fff_seq, m_settings_sla;

    PrinterTechnology current_printer_technology() const;
    const DynamicPrintConfig *config() const;

    template<class Self>
    static auto & get_slot(Self *self) {
        PrinterTechnology ptech = self->current_printer_technology();

        auto *ptr = &self->m_settings_fff;

        if (ptech == ptSLA) {
            ptr = &self->m_settings_sla;
        } else if (ptech == ptFFF && self->config()) {
            auto co_opt = self->config()->template option<ConfigOptionBool>(
                "complete_objects");
            if (co_opt && co_opt->value)
                ptr = &self->m_settings_fff_seq;
            else
                ptr = &self->m_settings_fff;
        }

        return *ptr;
    }

    template<class Self>
    static auto& get_ref(Self *self) { return get_slot(self).vals; }

public:
    explicit ArrangeSettingsDb_AppCfg(
        AppConfig *appcfg,
        std::function<const DynamicPrintConfig *(void)> cfgfn,
        std::function<PrinterTechnology(void)> printtech_getter);

    float get_distance_from_objects() const override { return get_ref(this).d_obj; }
    float get_distance_from_bed() const  override { return get_ref(this).d_bed; }
    bool  is_rotation_enabled() const override { return get_ref(this).rotations; }

    XLPivots get_xl_alignment() const override { return m_settings_fff.vals.xl_align; }
    GeometryHandling get_geometry_handling() const override { return m_settings_fff.vals.geom_handling; }
    ArrangeStrategy get_arrange_strategy() const override { return m_settings_fff.vals.arr_strategy; }

    void distance_from_obj_range(float &min, float &max) const override;
    void distance_from_bed_range(float &min, float &max) const override;

    ArrangeSettingsDb& set_distance_from_objects(float v) override;
    ArrangeSettingsDb& set_distance_from_bed(float v) override;
    ArrangeSettingsDb& set_rotation_enabled(bool v) override;

    ArrangeSettingsDb& set_xl_alignment(XLPivots v) override;
    ArrangeSettingsDb& set_geometry_handling(GeometryHandling v) override;
    ArrangeSettingsDb& set_arrange_strategy(ArrangeStrategy v) override;

    Values get_defaults() const override { return get_slot(this).defaults; }
};

} // namespace Slic3r

#endif // ARRANGESETTINGSDB_APPCFG_HPP
