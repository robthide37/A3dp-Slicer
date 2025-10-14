///|/ Copyright (c) SuperSlicer 2025 Durand Rémi @supermerill
///|/
///|/ SuperSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "RegionSettings.hpp"

#include "Print.hpp"

namespace Slic3r {

bool RegionSettings::SettingsValue::operator==(const SettingsValue &rhs) const {
    if (val_options.size() != rhs.val_options.size()) {
        assert(false);
        return false;
    }
    for (size_t i = 0; i < val_options.size(); i++) {
        if (*val_options[i] != *rhs.val_options[i]) {
            return false;
        }
    }
    return true;
}
bool RegionSettings::SettingsValue::operator!=(const SettingsValue &rhs) const { return !(*this == rhs); }
bool RegionSettings::SettingsValue::operator<(const SettingsValue &rhs) const {
    if (val_options.size() != rhs.val_options.size()) {
        assert(false);
        return val_options.size() < rhs.val_options.size();
    }
    for (size_t i = 0; i < val_options.size(); i++) {
        assert(val_options[i]->type() == rhs.val_options[i]->type());
        if (*val_options[i] != *rhs.val_options[i]) {
            return *val_options[i] < *rhs.val_options[i];
        }
    }
    return false;
}

void RegionSettings::segregate_regions(const ExPolygon &my_srf, const std::set<LayerRegion *> lregions) {
    this->key_areas.clear();
    BoundingBox my_srf_bb(my_srf.contour.points);
    my_srf_bb.offset(SCALED_EPSILON * 3);
    std::map<const Surface *, ExPolygons> srf_to_optimized_overlap;
    for (const t_config_option_keys &opt_keys : opt_key_set) {
        std::map<SettingsValue, ClipExpoly> &opt_2_areas =
            this->key_areas[this->config.option(opt_keys.front())];
        // check if it needs to be added
        bool many_values = false;
        RegionSettings::SettingsValue default_value = RegionSettings::SettingsValue::create(this->config, this->config, opt_keys);
        if (!lregions.empty()) {
            for (const LayerRegion *region : lregions) {
                if (RegionSettings::SettingsValue::create(this->config, region->region().config(), opt_keys) != default_value) {
                    many_values = true;
                    break;
                }
            }
        }
        if (many_values) {
            for (const LayerRegion *lregion : lregions) {
                SettingsValue settings_group_key = RegionSettings::SettingsValue::create(this->config, lregion->region().config(), opt_keys);
                ClipExpoly &areas = opt_2_areas[settings_group_key];
                // only get surfaces that overlap with my bb
                for (const Surface &surface : lregion->slices()) {
                    auto it_srf = srf_to_optimized_overlap.find(&surface);
                    if (it_srf == srf_to_optimized_overlap.end()) {
                        BoundingBox test_bb(surface.expolygon.contour.points);
                        if (test_bb.overlap(my_srf_bb) /*&& surface.expolygon.overlaps(my_srf)done by clip*/) {
                            srf_to_optimized_overlap[&surface] =
                                offset_ex(ClipperUtils::clip_clipper_polygons_with_subject_bbox(surface.expolygon,
                                                                                                my_srf_bb),
                                          SCALED_EPSILON * 10);
                        } else {
                            //it_srf = srf_to_optimized_overlap.emplace(&surface, ExPolygons{});
                            srf_to_optimized_overlap[&surface] = ExPolygons{};
                        }
                        it_srf = srf_to_optimized_overlap.find(&surface);
                        assert(it_srf != srf_to_optimized_overlap.end());
                    }
                    if (!it_srf->second.empty()) {
                        append(areas.expolys, it_srf->second);
                    }
                }
            }
            // delete empty ones
            //std::erase_if(opt_2_areas, [](auto &kv) { return kv.second.expolys.empty(); });
            for (auto it = opt_2_areas.begin(); it != opt_2_areas.end();) {
                if (it->second.expolys.empty()) {
                    it = opt_2_areas.erase(it);
                } else {
                    ++it;
                }
            }
            // check if there is really an overlap with many regions
            if (opt_2_areas.size() == 1) {
                //no overlap. don't remove it to use as the default, just remove the area
                opt_2_areas.begin()->second.clear();
            }
            // union the surfaces and un-offset them.
            for (auto &[opt, clip_expolys] : opt_2_areas) {
                clip_expolys.expolys = offset_ex(union_ex(clip_expolys.expolys), -SCALED_EPSILON * 10);
                clip_expolys.compute_bb();
            }
        }
        if (opt_2_areas.empty()) {
            // only one value for evrything
            opt_2_areas[default_value] = {};
        }
    }
}
void RegionSettings::ClipExpoly::compute_bb() {
    for (ExPolygon &expoly : this->expolys) {
        bboxes.emplace_back(expoly.contour.points);
    }
    assert(bboxes.size() == expolys.size());
}
ExPolygons RegionSettings::ClipExpoly::intersections(const ExPolygons &to_clip) const {
    if (this->expolys.empty()) {
        return to_clip;
    }
    ExPolygons intersections;
    for (const ExPolygon &expoly : to_clip) {
        BoundingBox bb_contour(expoly.contour.points);
        for (size_t i = 0; i < this->expolys.size(); ++i) {
            if (bb_contour.overlap(this->bboxes[i])) {
                append(intersections, intersection_ex(expoly, this->expolys[i]));
            }
        }
    }
    return intersections;
}

ExPolygons RegionSettings::ClipExpoly::intersections(coord_t offset, const ExPolygons &to_clip) const {
    if (this->expolys.empty()) {
        return to_clip;
    }
    ExPolygons intersections;
    for (const ExPolygon &expoly : to_clip) {
        BoundingBox bb_contour(expoly.contour.points);
        for (size_t i = 0; i < this->expolys.size(); ++i) {
            BoundingBox bb_offseted(this->bboxes[i]);
            bb_offseted.offset(offset);
            if (bb_offseted.overlap(bb_contour)) {
                append(intersections, intersection_ex({expoly}, offset_ex(this->expolys[i], offset)));
            }
        }
    }
    return intersections;
}
void RegionSettings::ClipExpoly::clear() {
    expolys.clear();
    bboxes.clear();
}

} // namespace Slic3r
