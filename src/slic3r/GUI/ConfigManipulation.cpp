// #include "libslic3r/GCodeSender.hpp"
#include "ConfigManipulation.hpp"
#include "I18N.hpp"
#include "GUI_App.hpp"
#include "format.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "MsgDialog.hpp"

#include <wx/msgdlg.h>

namespace Slic3r {
namespace GUI {

void ConfigManipulation::apply(DynamicPrintConfig* config, DynamicPrintConfig* new_config)
{
    bool modified = false;
    for (auto opt_key : config->diff(*new_config)) {
        config->set_key_value(opt_key, new_config->option(opt_key)->clone());
        modified = true;
    }

    if (modified && load_config != nullptr)
        load_config();
}

void ConfigManipulation::toggle_field(const std::string& opt_key, const bool toggle, int opt_index/* = -1*/)
{
    if (local_config) {
        if (local_config->option(opt_key) == nullptr)
            return;
    }
    cb_toggle_field(opt_key, toggle, opt_index);
}

// !! if using cb_value_change(X) or somthgin like that, you need a special code in Field.cpp (search for 'update_print_fff_config')
void ConfigManipulation::update_print_fff_config(DynamicPrintConfig* config, const bool is_global_config)
{
    // #ys_FIXME_to_delete
    //! Temporary workaround for the correct updates of the TextCtrl (like "layer_height"):
    // KillFocus() for the wxSpinCtrl use CallAfter function. So,
    // to except the duplicate call of the update() after dialog->ShowModal(),
    // let check if this process is already started.
    if (is_msg_dlg_already_exist)
        return;

    // layer_height shouldn't be equal to zero
    if (config->opt_float("layer_height") < EPSILON)
    {
        const wxString msg_text = _(L("Layer height is not valid.\n\nThe layer height will be reset to 0.01."));
        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("Layer height")), wxICON_WARNING | wxOK);
        DynamicPrintConfig new_conf = *config;
        is_msg_dlg_already_exist = true;
        dialog.ShowModal();
        new_conf.set_key_value("layer_height", new ConfigOptionFloat(0.01));
        apply(config, &new_conf);
        is_msg_dlg_already_exist = false;
    }

    if (config->option<ConfigOptionFloatOrPercent>("first_layer_height")->value < EPSILON)
    {
        const wxString msg_text = _(L("First layer height is not valid.\n\nThe first layer height will be reset to 0.01."));
        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("First layer height")), wxICON_WARNING | wxOK);
        DynamicPrintConfig new_conf = *config;
        is_msg_dlg_already_exist = true;
        dialog.ShowModal();
        new_conf.set_key_value("first_layer_height", new ConfigOptionFloatOrPercent(0.01, false));
        apply(config, &new_conf);
        is_msg_dlg_already_exist = false;
    }

    double fill_density = config->option<ConfigOptionPercent>("fill_density")->value;

    if (config->opt_bool("spiral_vase") && !(
        config->opt_int("top_solid_layers") == 0
        && fill_density == 0
        && config->opt_bool("support_material") == false
        && config->opt_int("support_material_enforce_layers") == 0
        && config->opt_enum<PerimeterGeneratorType>("perimeter_generator") == PerimeterGeneratorType::Classic
        && config->opt_bool("exact_last_layer_height") == false
        && config->opt_bool("ensure_vertical_shell_thickness") == true
        && config->opt_bool("infill_dense") == false
        && config->opt_bool("extra_perimeters") == false
        && config->opt_bool("extra_perimeters_overhangs") == false
        && config->opt_bool("extra_perimeters_odd_layers") == false
        && config->opt_bool("overhangs_reverse") == false
        && config->opt_bool("gap_fill_last") == false
        && config->opt_int("solid_over_perimeters") == 0
        && config->option("seam_notch_all")->get_float() == 0
        && config->option("seam_notch_inner")->get_float() == 0
        && config->option("seam_notch_outer")->get_float() == 0
        )) {
        wxString msg_text = _(L("The Spiral Vase mode requires:\n"
            "- no top solid layers\n"
            "- 0% fill density\n"
            "- classic perimeter slicing\n"
            "- no support material\n"
            "- Ensure vertical shell thickness enabled\n"
            "- disabled 'no solid infill over perimeters'\n"
            "- unchecked 'exact last layer height'\n"
            "- unchecked 'dense infill'\n"
            "- unchecked 'extra perimeters'"
            "- unchecked 'gap fill after last perimeter'"
            "- disabled  'no solid fill over X perimeters'"
            "- disabled 'seam notch'"));
        if (is_global_config)
            msg_text += "\n\n" + _(L("Shall I adjust those settings in order to enable Spiral Vase?"));
        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("Spiral Vase")),
            wxICON_WARNING | (is_global_config ? wxYES | wxNO : wxOK));
        DynamicPrintConfig new_conf = *config;
        auto answer = dialog.ShowModal();
        bool support = true;
        if (!is_global_config) {
            if (this->local_config->get().optptr("spiral_vase"))
                new_conf.set_key_value("spiral_vase", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("top_solid_layers"))
                new_conf.set_key_value("top_solid_layers", new ConfigOptionInt(0));
            else if (this->local_config->get().optptr("fill_density"))
                new_conf.set_key_value("fill_density", new ConfigOptionPercent(0));
            else if (this->local_config->get().optptr("perimeter_generator"))
                new_conf.set_key_value("perimeter_generator", new ConfigOptionEnum<PerimeterGeneratorType>(PerimeterGeneratorType::Classic));
            else if (this->local_config->get().optptr("support_material"))
                new_conf.set_key_value("support_material", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("support_material_enforce_layers"))
                new_conf.set_key_value("support_material_enforce_layers", new ConfigOptionInt(0));
            else if (this->local_config->get().optptr("exact_last_layer_height"))
                new_conf.set_key_value("exact_last_layer_height", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("ensure_vertical_shell_thickness"))
                new_conf.set_key_value("ensure_vertical_shell_thickness", new ConfigOptionBool(true));
            else if (this->local_config->get().optptr("infill_dense"))
                new_conf.set_key_value("infill_dense", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("extra_perimeters"))
                new_conf.set_key_value("extra_perimeters", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("extra_perimeters_overhangs"))
                new_conf.set_key_value("extra_perimeters_overhangs", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("extra_perimeters_odd_layers"))
                new_conf.set_key_value("extra_perimeters_odd_layers", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("overhangs_reverse"))
                new_conf.set_key_value("overhangs_reverse", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("gap_fill_last"))
                new_conf.set_key_value("gap_fill_last", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("solid_over_perimeters"))
                new_conf.set_key_value("solid_over_perimeters", new ConfigOptionInt(0));
            else if (this->local_config->get().optptr("seam_notch_all"))
                new_conf.set_key_value("seam_notch_all", new ConfigOptionFloatOrPercent(0, false));
            else if (this->local_config->get().optptr("seam_notch_inner"))
                new_conf.set_key_value("seam_notch_all", new ConfigOptionFloatOrPercent(0, false));
            else if (this->local_config->get().optptr("seam_notch_outer"))
                new_conf.set_key_value("seam_notch_all", new ConfigOptionFloatOrPercent(0, false));
            this->local_config->apply_only(new_conf, this->local_config->keys(), true);
        } else if (answer == wxID_YES) {
            new_conf.set_key_value("top_solid_layers", new ConfigOptionInt(0));
            new_conf.set_key_value("fill_density", new ConfigOptionPercent(0));
            new_conf.set_key_value("perimeter_generator", new ConfigOptionEnum<PerimeterGeneratorType>(PerimeterGeneratorType::Classic));
            new_conf.set_key_value("support_material", new ConfigOptionBool(false));
            new_conf.set_key_value("support_material_enforce_layers", new ConfigOptionInt(0));
            new_conf.set_key_value("exact_last_layer_height", new ConfigOptionBool(false));
            new_conf.set_key_value("ensure_vertical_shell_thickness", new ConfigOptionBool(true));
            new_conf.set_key_value("infill_dense", new ConfigOptionBool(false));
            new_conf.set_key_value("extra_perimeters", new ConfigOptionBool(false));
            new_conf.set_key_value("extra_perimeters_overhangs", new ConfigOptionBool(false));
            new_conf.set_key_value("extra_perimeters_odd_layers", new ConfigOptionBool(false));
            new_conf.set_key_value("overhangs_reverse", new ConfigOptionBool(false));
            new_conf.set_key_value("gap_fill_last", new ConfigOptionBool(false));
            new_conf.set_key_value("solid_over_perimeters", new ConfigOptionInt(0));
            new_conf.set_key_value("seam_notch_all", new ConfigOptionFloatOrPercent(0, false));
            new_conf.set_key_value("seam_notch_inner", new ConfigOptionFloatOrPercent(0, false));
            new_conf.set_key_value("seam_notch_outer", new ConfigOptionFloatOrPercent(0, false));
            fill_density = 0;
            support = false;
        } else {
            new_conf.set_key_value("spiral_vase", new ConfigOptionBool(false));
        }
        apply(config, &new_conf);
        if (cb_value_change) {
            cb_value_change("fill_density", fill_density);
            if (!support)
                cb_value_change("support_material", false);
    }
    }

    if (config->opt_bool("wipe_tower") && config->opt_bool("support_material") &&
        ((ConfigOptionEnumGeneric*)config->option("support_material_contact_distance_type"))->value != zdNone &&
        (config->opt_int("support_material_extruder") != 0 || config->opt_int("support_material_interface_extruder") != 0)) {
        wxString msg_text = _(L("The Wipe Tower currently supports the non-soluble supports only\n"
            "if they are printed with the current extruder without triggering a tool change.\n"
            "(both support_material_extruder and support_material_interface_extruder need to be set to 0)."));
        if (is_global_config)
            msg_text += "\n\n" + _(L("Shall I adjust those settings in order to enable the Wipe Tower?"));
        MessageDialog dialog (m_msg_dlg_parent, msg_text, _(L("Wipe Tower")),
            wxICON_WARNING | (is_global_config ? wxYES | wxNO : wxOK));
        DynamicPrintConfig new_conf = *config;
        auto answer = dialog.ShowModal();
        if (!is_global_config) {
            if (this->local_config->get().optptr("wipe_tower"))
                new_conf.set_key_value("wipe_tower", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("support_material_extruder"))
                new_conf.set_key_value("support_material_extruder", new ConfigOptionInt(0));
            else if (this->local_config->get().optptr("support_material_interface_extruder"))
                new_conf.set_key_value("support_material_interface_extruder", new ConfigOptionInt(0));
            else if (this->local_config->get().optptr("support_material_contact_distance_type"))
                new_conf.set_key_value("support_material_contact_distance_type", new ConfigOptionEnum<SupportZDistanceType>(zdNone));
            else if (this->local_config->get().optptr("support_material"))
                new_conf.set_key_value("support_material", new ConfigOptionBool(false));
            this->local_config->apply_only(new_conf, this->local_config->keys(), true);
        } else if (answer == wxID_YES) {
            new_conf.set_key_value("support_material_extruder", new ConfigOptionInt(0));
            new_conf.set_key_value("support_material_interface_extruder", new ConfigOptionInt(0));
        } else
            new_conf.set_key_value("wipe_tower", new ConfigOptionBool(false));
        apply(config, &new_conf);
    }

    if (config->opt_bool("wipe_tower") && config->opt_bool("support_material") &&
        ((ConfigOptionEnumGeneric*)config->option("support_material_contact_distance_type"))->value == zdNone &&
        !config->opt_bool("support_material_synchronize_layers")) {
        wxString msg_text = _(L("For the Wipe Tower to work with the soluble supports, the support layers\n"
            "need to be synchronized with the object layers."));
        if (is_global_config)
            msg_text += "\n\n" + _(L("Shall I synchronize support layers in order to enable the Wipe Tower?"));
        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("Wipe Tower")),
            wxICON_WARNING | (is_global_config ? wxYES | wxNO : wxOK));
        DynamicPrintConfig new_conf = *config;
        auto answer = dialog.ShowModal();
        if (!is_global_config) {
            if (this->local_config->get().optptr("wipe_tower"))
                new_conf.set_key_value("wipe_tower", new ConfigOptionBool(false));
            else if (this->local_config->get().optptr("support_material_synchronize_layers"))
                new_conf.set_key_value("support_material_synchronize_layers", new ConfigOptionBool(true));
            else if (this->local_config->get().optptr("support_material_contact_distance_type"))
                new_conf.set_key_value("support_material_contact_distance_type", new ConfigOptionEnum<SupportZDistanceType>(zdFilament));
            else if (this->local_config->get().optptr("support_material"))
                new_conf.set_key_value("support_material", new ConfigOptionBool(false));
            this->local_config->apply_only(new_conf, this->local_config->keys(), true);
        } else if (answer == wxID_YES) {
            new_conf.set_key_value("support_material_synchronize_layers", new ConfigOptionBool(true));
        } else {
            new_conf.set_key_value("wipe_tower", new ConfigOptionBool(false));
        }
        apply(config, &new_conf);
    }

    if (config->opt_float("brim_width") > 0 && config->opt_float("brim_separation") >= config->opt_float("brim_width")) {
        wxString msg_text = _(L("It's not possible to use a bigger value for the brim offset than the brim width, as it won't extrude anything."
            " Brim offset have to be lower than the brim width."));
        if (is_global_config) {
            msg_text += "\n\n" + _(L("Shall I switch the brim offset to 0?"));
            wxMessageDialog dialog(nullptr, msg_text, _(L("Brim configuration")),
                wxICON_WARNING | (is_global_config ? wxYES | wxNO : wxOK));
            auto answer = dialog.ShowModal();
            if (!is_global_config || answer == wxID_YES) {
                DynamicPrintConfig new_conf = *config;
                new_conf.set_key_value("brim_separation", new ConfigOptionFloat(0));
                apply(config, &new_conf);
            }
        }
    }

    static bool support_material_overhangs_queried = false;

    // Check "support_material" and "overhangs" relations only on global settings level
    if (is_global_config && config->opt_bool("support_material")) {
        // Ask only once.
        if (!m_support_material_overhangs_queried) {
            m_support_material_overhangs_queried = true;
            if (config->option<ConfigOptionFloatOrPercent>("overhangs_width_speed") == 0) {
                wxString msg_text = _(L("Supports work better, if the following feature is enabled:\n"
                    "- overhangs with bridge speed & fan"));
                if (is_global_config) {
                    msg_text += "\n\n" + _(L("Shall I adjust those settings for supports?"));
                MessageDialog dialog(m_msg_dlg_parent, msg_text, _L("Support Generator"), wxICON_WARNING | wxYES | wxNO);
                    DynamicPrintConfig new_conf = *config;
                    auto answer = dialog.ShowModal();
                if (answer == wxID_YES) {
                        // Enable "detect bridging perimeters".
                        new_conf.set_key_value("overhangs_width_speed", new ConfigOptionFloatOrPercent(50, true));
                    } else if (answer == wxID_NO) {
                        // Do nothing, leave supports on and "detect bridging perimeters" off.
                    } else if (answer == wxID_CANCEL) {
                        // Disable supports.
                        new_conf.set_key_value("support_material", new ConfigOptionBool(false));
                        support_material_overhangs_queried = false;
                    }
                    apply(config, &new_conf);
                }
            }
        }
    } else {
        m_support_material_overhangs_queried = false;
    }

    if (config->option<ConfigOptionPercent>("fill_density")->value == 100) {
        std::string  fill_pattern = config->option<ConfigOptionEnum<InfillPattern>>("fill_pattern")->serialize();
        const std::vector<std::string>& top_fill_pattern_values = config->def()->get("top_fill_pattern")->enum_values;
        bool correct_100p_fill = std::find(top_fill_pattern_values.begin(), top_fill_pattern_values.end(), fill_pattern) != top_fill_pattern_values.end();
        if (!correct_100p_fill) {
            const std::vector<std::string>& bottom_fill_pattern_values = config->def()->get("bottom_fill_pattern")->enum_values;
            correct_100p_fill = std::find(bottom_fill_pattern_values.begin(), bottom_fill_pattern_values.end(), fill_pattern) != bottom_fill_pattern_values.end();
        }
        if (!correct_100p_fill) {
            const std::vector<std::string>& bottom_fill_pattern_values = config->def()->get("solid_fill_pattern")->enum_values;
            correct_100p_fill = std::find(bottom_fill_pattern_values.begin(), bottom_fill_pattern_values.end(), fill_pattern) != bottom_fill_pattern_values.end();
        }
        if (!correct_100p_fill) {
            // get fill_pattern name from enum_labels for using this one at dialog_msg
            const ConfigOptionDef* fill_pattern_def = config->def()->get("fill_pattern");
            assert(fill_pattern_def != nullptr);
            auto it_pattern = std::find(fill_pattern_def->enum_values.begin(), fill_pattern_def->enum_values.end(), fill_pattern);
            assert(it_pattern != fill_pattern_def->enum_values.end());
            if (it_pattern != fill_pattern_def->enum_values.end()) {
                wxString msg_text = GUI::format_wxstr(_L("The %1% infill pattern is not supposed to work at 100%% density."),
                    _(fill_pattern_def->enum_labels[it_pattern - fill_pattern_def->enum_values.begin()]));
                if (is_global_config) {
                    msg_text += "\n\n" + _L("Shall I switch to rectilinear fill pattern?");
                MessageDialog dialog(m_msg_dlg_parent, msg_text, _L("Infill"),
                        wxICON_WARNING | (is_global_config ? wxYES | wxNO : wxOK));
                    DynamicPrintConfig new_conf = *config;
                    auto answer = dialog.ShowModal();
                    if (!is_global_config || answer == wxID_YES) {
                        new_conf.set_key_value("fill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinear));
                        fill_density = 100;
                    } else
                        fill_density = wxGetApp().preset_bundle->fff_prints.get_selected_preset().config.option<ConfigOptionPercent>("fill_density")->value;
                    new_conf.set_key_value("fill_density", new ConfigOptionPercent(fill_density));
                    apply(config, &new_conf);
                    if (cb_value_change)
                        cb_value_change("fill_density", fill_density);
                }
            }
        }
    }
}

//note: printer options are toogled in TabPrinter::toggle_options()

void ConfigManipulation::toggle_print_fff_options(DynamicPrintConfig* config)
{
    bool have_perimeters = config->opt_int("perimeters") > 0;
    for (auto el : { "ensure_vertical_shell_thickness", "external_perimeter_speed", "extra_perimeters", "extra_perimeters_odd_layers",
        "external_perimeters_first", "external_perimeter_extrusion_width", "external_perimeter_extrusion_spacing","external_perimeter_extrusion_change_odd_layers",
        "overhangs","perimeter_speed",
        "seam_position", "small_perimeter_speed", "small_perimeter_min_length", " small_perimeter_max_length", "spiral_vase",
        "perimeter_generator", "seam_notch_all", "seam_notch_inner", "seam_notch_outer"})
        toggle_field(el, have_perimeters);

    bool has_spiral_vase = have_perimeters && config->opt_bool("spiral_vase");

    for (auto el : { "external_perimeters_vase", "external_perimeters_nothole", "external_perimeters_hole", "perimeter_bonding"})
        toggle_field(el, config->opt_bool("external_perimeters_first"));

    bool have_arachne = have_perimeters && config->opt_enum<PerimeterGeneratorType>("perimeter_generator") == PerimeterGeneratorType::Arachne;
    for (auto el : { "wall_transition_length", "wall_transition_filter_deviation", "wall_transition_angle", "wall_distribution_count", "min_feature_size", "min_bead_width", "aaa" })
       toggle_field(el, have_arachne);


    for (auto el : {"perimeter_loop", "extra_perimeters_overhangs", "no_perimeter_unsupported_algo",
        "thin_perimeters", "overhangs_reverse", "perimeter_round_corners"})
        toggle_field(el, have_perimeters && !have_arachne);

    toggle_field("only_one_perimeter_top", have_perimeters);
    toggle_field("only_one_perimeter_first_layer", config->opt_int("perimeters") > 1);
    toggle_field("overhangs_width", config->option<ConfigOptionFloatOrPercent>("overhangs_width_speed")->value > 0);
    toggle_field("overhangs_reverse_threshold", have_perimeters && config->opt_bool("overhangs_reverse"));
    toggle_field("overhangs_speed_enforce", have_perimeters && !config->opt_bool("perimeter_loop"));
    toggle_field("min_width_top_surface", have_perimeters && config->opt_bool("only_one_perimeter_top"));
    toggle_field("thin_perimeters_all", have_perimeters && config->option("thin_perimeters")->get_float() != 0 && !have_arachne);
    bool have_thin_wall = !have_arachne && have_perimeters;
    toggle_field("thin_walls", have_thin_wall);
    for (auto el : { "thin_walls_min_width", "thin_walls_overlap", "thin_walls_merge" })
        toggle_field(el, have_thin_wall && config->opt_bool("thin_walls"));

    for (auto el : { "seam_angle_cost", "seam_travel_cost", "seam_visibility" })
        toggle_field(el, have_perimeters && config->option<ConfigOptionEnum<SeamPosition>>("seam_position")->value == SeamPosition::spCost);

    toggle_field("perimeter_loop_seam", config->opt_bool("perimeter_loop"));

    bool have_notch = have_perimeters && (config->option("seam_notch_all")->get_float() != 0 ||
                                          config->option("seam_notch_inner")->get_float() != 0 ||
                                          config->option("seam_notch_outer")->get_float() != 0);
    toggle_field("seam_notch_angle", have_notch);

    bool have_gap_fill = !have_arachne;
    toggle_field("gap_fill_enabled", have_gap_fill);
    for (auto el : { "gap_fill_extension", "gap_fill_last", "gap_fill_max_width", "gap_fill_min_area", "gap_fill_min_length", "gap_fill_min_width" })
        toggle_field(el, config->opt_bool("gap_fill_enabled") && have_gap_fill);
    // gap fill  can appear in infill
    //toggle_field("gap_fill_speed", have_perimeters && config->opt_bool("gap_fill_enabled"));

    for (auto el : { "fuzzy_skin_thickness", "fuzzy_skin_point_dist" })
        toggle_field(el, config->option<ConfigOptionEnum<FuzzySkinType>>("fuzzy_skin")->value != FuzzySkinType::None);

    toggle_field("avoid_crossing_not_first_layer", config->opt_bool("avoid_crossing_perimeters"));

    bool have_infill = config->option<ConfigOptionPercent>("fill_density")->value > 0;
    // infill_extruder uses the same logic as in Print::extruders()
    for (auto el : { "fill_pattern", "infill_connection", "infill_every_layers", "infill_only_where_needed",
                    "solid_infill_every_layers", "solid_infill_below_area", "infill_extruder", "infill_anchor_max" })
        toggle_field(el, have_infill);
    // Only allow configuration of open anchors if the anchoring is enabled.
    bool has_infill_anchors = have_infill && config->option<ConfigOptionEnum<InfillConnection>>("infill_connection")->value != InfillConnection::icNotConnected;
    toggle_field("infill_anchor_max", has_infill_anchors);
    has_infill_anchors = has_infill_anchors && config->option<ConfigOptionFloatOrPercent>("infill_anchor_max")->value > 0;
    toggle_field("infill_anchor", has_infill_anchors);

    bool can_have_infill_dense = config->option<ConfigOptionPercent>("fill_density")->value < 50;
    for (auto el : { "infill_dense" })
        toggle_field(el, can_have_infill_dense);

    bool have_infill_dense = config->opt_bool("infill_dense") && can_have_infill_dense;
    for (auto el : { "infill_dense_algo" })
        toggle_field(el, have_infill_dense);
    if(have_infill)
        for (auto el : { "infill_every_layers", "infill_only_where_needed" })
            toggle_field(el, !have_infill_dense);

    bool has_top_solid_infill 	 = config->opt_int("top_solid_layers") > 0 || has_spiral_vase;
    bool has_bottom_solid_infill = config->opt_int("bottom_solid_layers") > 0;
    bool has_solid_infill 		 = has_top_solid_infill || has_bottom_solid_infill || (have_infill && (config->opt_int("solid_infill_every_layers") > 0 || config->opt_float("solid_infill_below_area") > 0));
    // solid_infill_extruder uses the same logic as in Print::extruders()
    for (auto el : { "top_fill_pattern", "bottom_fill_pattern", "solid_fill_pattern", "enforce_full_fill_volume", "external_infill_margin", "bridged_infill_margin",
        "solid_infill_extruder", "solid_infill_extrusion_width", "solid_nfill_extrusion_change_odd_layers", "solid_infill_extrusion_spacing", "solid_infill_speed" })
        toggle_field(el, has_solid_infill);

    toggle_field("infill_first", (has_solid_infill || have_infill));

    for (auto el : {"fill_angle_cross","fill_angle_increment", "fill_angle_template", "fill_angle_cross", "bridge_angle", "infill_extrusion_width",
                    "infill_extrusion_spacing", "infill_extrusion_change_odd_layers", "infill_speed" })
        toggle_field(el, have_infill || has_solid_infill);
        
    toggle_field("fill_angle", (have_infill || has_solid_infill) && ((ConfigOptionVectorBase*)config->option("fill_angle_template"))->size() == 0);

    toggle_field("top_solid_min_thickness", ! has_spiral_vase && has_top_solid_infill);
    toggle_field("bottom_solid_min_thickness", ! has_spiral_vase && has_bottom_solid_infill);

    //speed
    for (auto el : { "small_perimeter_min_length", "small_perimeter_max_length" })
        toggle_field(el, config->option("small_perimeter_speed")->get_float() > 0);

    bool has_ironing_pattern = config->opt_enum<InfillPattern>("top_fill_pattern") == InfillPattern::ipSmooth
        || config->opt_enum<InfillPattern>("bottom_fill_pattern") == InfillPattern::ipSmooth
        || config->opt_enum<InfillPattern>("solid_fill_pattern") == InfillPattern::ipSmooth;
    for (auto el : {"fill_smooth_width, fill_smooth_distribution" })
        toggle_field(el, has_ironing_pattern);

    for (auto el : { "ironing", "top_fill_pattern", "infill_connection_top",  "top_infill_extrusion_width",  "top_infill_extrusion_spacing", "top_solid_infill_speed" })
        toggle_field(el, has_top_solid_infill || (has_spiral_vase && has_bottom_solid_infill));

    for (auto el : { "bottom_fill_pattern", "infill_connection_bottom" })
        toggle_field(el, has_bottom_solid_infill);

    for (auto el : { "solid_fill_pattern", "infill_connection_solid", "bridge_fill_pattern", "infill_connection_bridge" })
        toggle_field(el, has_solid_infill); // should be top_solid_layers") > 1 || bottom_solid_layers") > 1

    for (auto el : { "hole_to_polyhole_threshold", "hole_to_polyhole_twisted" })
        toggle_field(el, config->opt_bool("hole_to_polyhole"));

    bool have_skirt = config->opt_int("skirts") > 0;
    toggle_field("skirt_height", have_skirt && config->opt_enum<DraftShield>("draft_shield") != dsEnabled);
    toggle_field("skirt_width", have_skirt);
    for (auto el : { "skirt_brim", "skirt_distance", "skirt_distance_from_brim", "draft_shield", "min_skirt_length" })
        toggle_field(el, have_skirt);

    bool have_brim = config->opt_float("brim_width") > 0 || config->opt_float("brim_width_interior") > 0;
    toggle_field("brim_separation", have_brim);
    // perimeter_extruder uses the same logic as in Print::extruders()
    toggle_field("perimeter_extruder", have_perimeters || have_brim);

    toggle_field("brim_ears", config->opt_float("brim_width") > 0);
    toggle_field("brim_inside_holes", config->opt_float("brim_width") > 0 && config->opt_float("brim_width_interior") == 0);
    toggle_field("brim_ears_max_angle", have_brim && config->opt_bool("brim_ears"));
    toggle_field("brim_ears_pattern", have_brim && config->opt_bool("brim_ears"));

    bool have_raft = config->opt_int("raft_layers") > 0;
    bool have_support_material = config->opt_bool("support_material") || have_raft;
    bool have_support_material_auto = have_support_material && config->opt_bool("support_material_auto");
    bool have_support_interface = have_support_material && (config->opt_int("support_material_interface_layers") > 0 || config->opt_int("support_material_bottom_interface_layers") > 0);
    bool have_support_soluble = have_support_material && ((ConfigOptionEnumGeneric*)config->option("support_material_contact_distance_type"))->value == zdNone;
    auto support_material_style = config->opt_enum<SupportMaterialStyle>("support_material_style");
    for (auto el : { "support_material_style", "support_material_pattern", "support_material_with_sheath",
                    "support_material_spacing", "support_material_angle", "support_material_angle_height", 
                    "support_material_bottom_interface_layers", "support_material_interface_layers",
                    "dont_support_bridges", "support_material_extrusion_width",
                    "support_material_contact_distance_type",
                    "support_material_xy_spacing",
                    "support_material_layer_height"})
        toggle_field(el, have_support_material);
    toggle_field("support_material_threshold", have_support_material_auto);
    toggle_field("support_material_bottom_contact_distance", have_support_material && ! have_support_soluble);
    toggle_field("support_material_closing_radius", have_support_material && support_material_style == smsSnug);

    for (auto el : { "support_material_contact_distance", "support_material_bottom_contact_distance" })
        toggle_field(el, have_support_material && !have_support_soluble);

    for (auto el : { "support_material_interface_pattern", "support_material_interface_spacing", "support_material_interface_extruder",
                    "support_material_interface_speed", "support_material_interface_contact_loops", "support_material_interface_layer_height"
                    "support_material_interface_angle", "support_material_interface_angle_increment"})
        toggle_field(el, have_support_material && have_support_interface);
    toggle_field("support_material_synchronize_layers", have_support_soluble);

    toggle_field("perimeter_extrusion_width", have_perimeters || have_brim);
    toggle_field("perimeter_extrusion_spacing", have_perimeters || have_brim);
    toggle_field("perimeter_extrusion_change_odd_layers", have_perimeters || have_brim);
    toggle_field("skirt_extrusion_width", have_skirt);
    toggle_field("support_material_extruder", have_support_material || have_skirt);
    toggle_field("support_material_speed", have_support_material || have_skirt || have_brim);
    toggle_field("brim_speed", have_brim || have_skirt);

    toggle_field("raft_contact_distance", have_raft && !have_support_soluble);
    for (auto el : { "raft_expansion", "first_layer_acceleration_over_raft", "first_layer_speed_over_raft",
        "raft_layer_height", "raft_interface_layer_height"})
        toggle_field(el, have_raft);

    //for default_extrusion_width/spacing, you need to ahve at least an extrusion_width with 0
    bool have_default_width = config->option("first_layer_extrusion_width")->get_float() == 0 ||
        (config->option("perimeter_extrusion_width")->get_float() == 0 && (have_perimeters || have_brim)) ||
        (config->option("external_perimeter_extrusion_width")->get_float() == 0 && have_perimeters) ||
        (config->option("infill_extrusion_width")->get_float() == 0 && (have_infill || has_solid_infill)) ||
        (config->option("solid_infill_extrusion_width")->get_float() == 0 && has_solid_infill) ||
        (config->option("top_infill_extrusion_width")->get_float() == 0 && has_top_solid_infill) ||
        (config->option("support_material_extrusion_width")->get_float() == 0 && have_support_material) ||
        (config->option("skirt_extrusion_width")->get_float() == 0 && have_skirt);
    toggle_field("extrusion_width", have_default_width);
    toggle_field("extrusion_spacing", have_default_width);

    bool has_PP_ironing = has_top_solid_infill && config->opt_bool("ironing");
    for (auto el : { "ironing_type", "ironing_flowrate", "ironing_spacing", "ironing_angle" })
        toggle_field(el, has_PP_ironing);

    bool has_ironing = has_PP_ironing || has_ironing_pattern;
    for (auto el : { "ironing_speed" })
        toggle_field(el, has_ironing);
    

    bool have_sequential_printing = config->opt_bool("complete_objects");
    for (auto el : { /*"extruder_clearance_radius", "extruder_clearance_height",*/ "complete_objects_one_skirt",
        "complete_objects_sort"})
        toggle_field(el, have_sequential_printing);
    toggle_field("parallel_objects_step", !have_sequential_printing);

    bool have_ooze_prevention = config->opt_bool("ooze_prevention");
    toggle_field("standby_temperature_delta", have_ooze_prevention);

    bool have_wipe_tower = config->opt_bool("wipe_tower");
    for (auto el : { "wipe_tower_x", "wipe_tower_y", "wipe_tower_width", "wipe_tower_rotation_angle", "wipe_tower_brim_width",
                     "wipe_tower_bridging", "wipe_tower_brim", "wipe_tower_no_sparse_layers", "single_extruder_multi_material_priming",
                     "wipe_tower_speed", "wipe_tower_wipe_starting_speed" })
        toggle_field(el, have_wipe_tower);

    bool have_avoid_crossing_perimeters = config->opt_bool("avoid_crossing_perimeters");
    toggle_field("avoid_crossing_perimeters_max_detour", have_avoid_crossing_perimeters);

    for (auto el : { "fill_smooth_width", "fill_smooth_distribution" })
        toggle_field(el, (has_top_solid_infill && config->option<ConfigOptionEnum<InfillPattern>>("top_fill_pattern")->value == InfillPattern::ipSmooth)
            || (has_bottom_solid_infill && config->option<ConfigOptionEnum<InfillPattern>>("bottom_fill_pattern")->value == InfillPattern::ipSmooth)
            || (has_solid_infill && config->option<ConfigOptionEnum<InfillPattern>>("solid_fill_pattern")->value == InfillPattern::ipSmooth)
            || (have_support_material && config->option<ConfigOptionEnum<InfillPattern>>("support_material_interface_pattern")->value == InfillPattern::ipSmooth));

    //TODO: can the milling_diameter or the milling_cutter be check to enable/disable this?
    for (auto el : { "milling_after_z", "milling_extra_size", "milling_speed" })
        toggle_field(el, config->opt_bool("milling_post_process"));
	
    //TODO: can the laser_power or the laser_tool be check to enable/disable this?
    for (auto el : { "laser_energy" })
        toggle_field(el, config->opt_bool("laser_support_interface_pp"));

    bool have_default_acceleration = config->option<ConfigOptionFloatOrPercent>("default_acceleration")->value > 0;
    for (auto el : { "perimeter_acceleration", "external_perimeter_acceleration", "thin_walls_acceleration" })
        toggle_field(el, have_default_acceleration && have_perimeters);
    toggle_field("infill_acceleration", have_default_acceleration && have_infill);
    toggle_field("solid_infill_acceleration", have_default_acceleration && has_solid_infill);
    toggle_field("top_solid_infill_acceleration", have_default_acceleration && has_top_solid_infill);
    toggle_field("ironing_acceleration", have_default_acceleration && has_ironing);
    toggle_field("support_material_acceleration", have_default_acceleration && (have_support_material || have_brim || have_skirt));
    toggle_field("support_material_interface_acceleration", have_default_acceleration && have_support_material && have_support_interface);
    toggle_field("brim_acceleration", have_default_acceleration && (have_brim || have_skirt));
    for (auto el : { "bridge_acceleration", "bridge_internal_acceleration", "overhangs_acceleration", "gap_fill_acceleration", "travel_acceleration", "travel_deceleration_use_target", "first_layer_acceleration" })
        toggle_field(el, have_default_acceleration);

    // for default speed, it needs at least a dependent field with a %
    toggle_field("default_speed", config->option<ConfigOptionFloatOrPercent>("perimeter_speed")->percent || 
        config->option<ConfigOptionFloatOrPercent>("solid_infill_speed")->percent || 
        config->option<ConfigOptionFloatOrPercent>("bridge_speed")->percent || 
        config->option<ConfigOptionFloatOrPercent>("support_material_speed")->percent);
    toggle_field("max_print_speed", config->opt_float("max_volumetric_speed") != 0);
}

void ConfigManipulation::update_print_sla_config(DynamicPrintConfig* config, const bool is_global_config/* = false*/)
{
    double head_penetration = config->opt_float("support_head_penetration");
    double head_width = config->opt_float("support_head_width");
    if (head_penetration > head_width) {
        wxString msg_text = _(L("Head penetration should not be greater than the head width."));

        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("Invalid Head penetration")), wxICON_WARNING | wxOK);
        DynamicPrintConfig new_conf = *config;
        if (dialog.ShowModal() == wxID_OK) {
            new_conf.set_key_value("support_head_penetration", new ConfigOptionFloat(head_width));
            apply(config, &new_conf);
        }
    }

    double pinhead_d = config->opt_float("support_head_front_diameter");
    double pillar_d = config->opt_float("support_pillar_diameter");
    if (pinhead_d > pillar_d) {
        wxString msg_text = _(L("Pinhead diameter should be smaller than the pillar diameter."));

        MessageDialog dialog(m_msg_dlg_parent, msg_text, _(L("Invalid pinhead diameter")), wxICON_WARNING | wxOK);

        DynamicPrintConfig new_conf = *config;
        if (dialog.ShowModal() == wxID_OK) {
            new_conf.set_key_value("support_head_front_diameter", new ConfigOptionFloat(pillar_d / 2.0));
            apply(config, &new_conf);
        }
    }
}

void ConfigManipulation::toggle_print_sla_options(DynamicPrintConfig* config)
{
    bool supports_en = config->opt_bool("supports_enable");

    toggle_field("support_head_front_diameter", supports_en);
    toggle_field("support_head_penetration", supports_en);
    toggle_field("support_head_width", supports_en);
    toggle_field("support_pillar_diameter", supports_en);
    toggle_field("support_small_pillar_diameter_percent", supports_en);
    toggle_field("support_max_bridges_on_pillar", supports_en);
    toggle_field("support_pillar_connection_mode", supports_en);
    toggle_field("support_buildplate_only", supports_en);
    toggle_field("support_base_diameter", supports_en);
    toggle_field("support_base_height", supports_en);
    toggle_field("support_base_safety_distance", supports_en);
    toggle_field("support_critical_angle", supports_en);
    toggle_field("support_max_bridge_length", supports_en);
    toggle_field("support_max_pillar_link_distance", supports_en);
    toggle_field("support_points_density_relative", supports_en);
    toggle_field("support_points_minimal_distance", supports_en);

    bool pad_en = config->opt_bool("pad_enable");

    toggle_field("pad_wall_thickness", pad_en);
    toggle_field("pad_wall_height", pad_en);
    toggle_field("pad_brim_size", pad_en);
    toggle_field("pad_max_merge_distance", pad_en);
 // toggle_field("pad_edge_radius", supports_en);
    toggle_field("pad_wall_slope", pad_en);
    toggle_field("pad_around_object", pad_en);
    toggle_field("pad_around_object_everywhere", pad_en);

    bool zero_elev = config->opt_bool("pad_around_object") && pad_en;

    toggle_field("support_object_elevation", supports_en && !zero_elev);
    toggle_field("pad_object_gap", zero_elev);
    toggle_field("pad_around_object_everywhere", zero_elev);
    toggle_field("pad_object_connector_stride", zero_elev);
    toggle_field("pad_object_connector_width", zero_elev);
    toggle_field("pad_object_connector_penetration", zero_elev);
}


} // GUI
} // Slic3r
