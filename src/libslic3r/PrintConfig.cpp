///|/ Copyright (c) Prusa Research 2016 - 2023 Vojtěch Bubník @bubnikv, Lukáš Matěna @lukasmatena, Lukáš Hejl @hejllukas, Tomáš Mészáros @tamasmeszaros, Oleksandra Iushchenko @YuSanka, Pavel Mikuš @Godrak, David Kocík @kocikdav, Enrico Turri @enricoturri1966, Filip Sykala @Jony01, Vojtěch Král @vojtechkral
///|/ Copyright (c) 2023 Pedro Lamas @PedroLamas
///|/ Copyright (c) 2023 Mimoja @Mimoja
///|/ Copyright (c) 2020 - 2021 Sergey Kovalev @RandoMan70
///|/ Copyright (c) 2021 Niall Sheridan @nsheridan
///|/ Copyright (c) 2021 Martin Budden
///|/ Copyright (c) 2021 Ilya @xorza
///|/ Copyright (c) 2020 Paul Arden @ardenpm
///|/ Copyright (c) 2020 rongith
///|/ Copyright (c) 2019 Spencer Owen @spuder
///|/ Copyright (c) 2019 Stephan Reichhelm @stephanr
///|/ Copyright (c) 2018 Martin Loidl @LoidlM
///|/ Copyright (c) SuperSlicer 2018 Remi Durand @supermerill
///|/ Copyright (c) 2016 - 2017 Joseph Lenox @lordofhyphens
///|/ Copyright (c) Slic3r 2013 - 2016 Alessandro Ranellucci @alranel
///|/ Copyright (c) 2016 Vanessa Ezekowitz @VanessaE
///|/ Copyright (c) 2015 Alexander Rössler @machinekoder
///|/ Copyright (c) 2014 Petr Ledvina @ledvinap
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "PrintConfig.hpp"
#include "Config.hpp"
#include "Flow.hpp"
#include "format.hpp"
#include "I18N.hpp"
#include "Semver.hpp"
#include "Utils.hpp"

#include "SLA/SupportTree.hpp"
#include "GCode/Thumbnails.hpp"

#include <set>
#include <unordered_set>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>
#include <boost/nowide/iostream.hpp>

#include <algorithm>
#include <float.h>

namespace Slic3r {

static t_config_enum_names enum_names_from_keys_map(const t_config_enum_values &enum_keys_map)
{
    t_config_enum_names names;
    int cnt = 0;
    for (const auto& kvp : enum_keys_map)
        cnt = std::max(cnt, kvp.second);
    cnt += 1;
    names.assign(cnt, "");
    for (const auto& kvp : enum_keys_map)
        names[kvp.second] = kvp.first;
    return names;
}

#define CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(NAME) \
    static t_config_enum_names s_keys_names_##NAME = enum_names_from_keys_map(s_keys_map_##NAME); \
    template<> const t_config_enum_values& ConfigOptionEnum<NAME>::get_enum_values() { return s_keys_map_##NAME; } \
    template<> const t_config_enum_names& ConfigOptionEnum<NAME>::get_enum_names() { return s_keys_names_##NAME; }

static const t_config_enum_values s_keys_map_ArcFittingType {
    { "disabled",       int(ArcFittingType::Disabled) },
    { "bambu",          int(ArcFittingType::Bambu) },
    { "emit_center",    int(ArcFittingType::EmitCenter) } // arwelder
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(ArcFittingType)

static const t_config_enum_values s_keys_map_PrinterTechnology{
    {"FFF",             ptFFF },
    {"SLA",             ptSLA },
    {"SLS",             ptSLS},
    {"CNC",             ptMill},
    {"LSR",             ptLaser},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(PrinterTechnology)


static const t_config_enum_values s_keys_map_CompleteObjectSort {
    {"object", cosObject},
    {"lowy", cosY},
    {"lowz", cosZ},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(CompleteObjectSort)

static const t_config_enum_values s_keys_map_OutputFormat {
    {"SL1", ofSL1},
    {"SL1_SVG", ofSL1_SVG},
    {"mCWS", ofMaskedCWS},
    {"AnyMono", ofAnycubicMono},
    {"AnyMonoX", ofAnycubicMonoX},
    {"AnyMonoSE", ofAnycubicMonoSE},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(OutputFormat)

static const t_config_enum_values s_keys_map_WipeAlgo {
    {"linear", waLinear},
    {"quadra", waQuadra},
    {"expo", waHyper},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(WipeAlgo)

static const t_config_enum_values s_keys_map_GCodeFlavor{
    {"reprapfirmware",  gcfRepRap},
    {"repetier",        gcfRepetier},
    {"teacup",          gcfTeacup},
    {"makerware",       gcfMakerWare},
    {"marlin",          gcfMarlinLegacy },
    {"marlin2",         gcfMarlinFirmware },
    {"klipper",         gcfKlipper},
    {"sailfish",        gcfSailfish},
    {"smoothie",        gcfSmoothie},
    {"sprinter",        gcfSprinter},
    {"mach3",           gcfMach3},
    {"machinekit",      gcfMachinekit},
    {"no-extrusion",    gcfNoExtrusion},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(GCodeFlavor)

static const t_config_enum_values s_keys_map_MachineLimitsUsage {
    {"emit_to_gcode",       int(MachineLimitsUsage::EmitToGCode)},
    {"time_estimate_only",  int(MachineLimitsUsage::TimeEstimateOnly)},
    {"limits",              int(MachineLimitsUsage::Limits)},
    {"ignore",              int(MachineLimitsUsage::Ignore)},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(MachineLimitsUsage)

static const t_config_enum_values s_keys_map_PrintHostType {
    { "prusalink",      htPrusaLink },
    { "prusaconnect",   htPrusaConnect },
    { "octoprint",      htOctoPrint },
    { "moonraker",      htMoonraker },
    { "duet",           htDuet },
    { "flashair",       htFlashAir },
    { "astrobox",       htAstroBox },
    { "repetier",       htRepetier },
    {"klipper",         htKlipper},
    {"mpmdv2",          htMPMDv2},
    { "mks",            htMKS },
    {"monoprice",       htMiniDeltaLCD },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(PrintHostType)

static const t_config_enum_values s_keys_map_AuthorizationType {
    {"key", atKeyPassword},
    {"user", atUserPassword},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(AuthorizationType)

static const t_config_enum_values s_keys_map_BridgeType {
    {"nozzle",  uint8_t(BridgeType::btFromNozzle)},
    {"height",  uint8_t(BridgeType::btFromHeight)},
    {"flow",    uint8_t(BridgeType::btFromFlow)},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(BridgeType)

static const t_config_enum_values s_keys_map_FuzzySkinType {
    { "none",           int(FuzzySkinType::None) },
    { "external",       int(FuzzySkinType::External) },
    { "shell",          int(FuzzySkinType::Shell) },
    { "all",            int(FuzzySkinType::All) }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(FuzzySkinType)

static const t_config_enum_values s_keys_map_InfillPattern {
    {"rectilinear",         ipRectilinear},
    {"rectilineargapfill",  ipRectilinearWGapFill},
    {"alignedrectilinear",  ipAlignedRectilinear},
    {"monotonic",           ipMonotonic},
    {"monotonicgapfill",    ipMonotonicWGapFill},
    {"grid",                ipGrid},
    {"triangles",           ipTriangles},
    {"stars",               ipStars},
    {"cubic",               ipCubic},
    {"line",                ipLine},
    {"monotoniclines",      ipMonotonicLines },
    {"concentric",          ipConcentric},
    {"concentricgapfill",   ipConcentricGapFill},
    {"honeycomb",           ipHoneycomb},
    {"3dhoneycomb",         ip3DHoneycomb},
    {"gyroid",              ipGyroid},
    {"hilbertcurve",        ipHilbertCurve},
    {"archimedeanchords",   ipArchimedeanChords},
    {"octagramspiral",      ipOctagramSpiral},
    {"smooth",              ipSmooth},
    {"smoothtriple",        ipSmoothTriple},
    {"smoothhilbert",       ipSmoothHilbert},
    {"rectiwithperimeter",  ipRectiWithPerimeter},
    {"scatteredrectilinear", ipScatteredRectilinear},
    {"sawtooth",            ipSawtooth},
    {"adaptivecubic",       ipAdaptiveCubic},
    {"supportcubic",        ipSupportCubic},
    {"lightning",           ipLightning},
    {"auto",                ipAuto}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(InfillPattern)

static const t_config_enum_values s_keys_map_IroningType {
    { "top",            int(IroningType::TopSurfaces) },
    { "topmost",        int(IroningType::TopmostOnly) },
    { "solid",          int(IroningType::AllSolid) }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(IroningType)

static const t_config_enum_values s_keys_map_PerimeterDirection{
    {"ccw_cw",  pdCCW_CW},
    {"ccw_ccw", pdCCW_CCW},
    {"cw_ccw",  pdCW_CCW},
    {"cw_cw",   pdCW_CW},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(PerimeterDirection);

static const t_config_enum_values s_keys_map_SlicingMode {
    { "regular",        int(SlicingMode::Regular) },
    { "even_odd",       int(SlicingMode::EvenOdd) },
    { "close_holes",    int(SlicingMode::CloseHoles) }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SlicingMode)

static const t_config_enum_values s_keys_map_SupportMaterialPattern {
    { "rectilinear",        smpRectilinear },
    { "rectilinear-grid",   smpRectilinearGrid },
    { "honeycomb",          smpHoneycomb }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SupportMaterialPattern)

static const t_config_enum_values s_keys_map_SupportMaterialStyle {
    { "grid",           smsGrid },
    { "snug",           smsSnug },
    { "tree",           smsTree },
    { "organic",        smsOrganic }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SupportMaterialStyle)

//unused
//static const t_config_enum_values s_keys_map_SupportMaterialInterfacePattern {
//    { "auto",           smipAuto },
//    { "rectilinear",    smipRectilinear },
//    { "concentric",     smipConcentric }
//};
//CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SupportMaterialInterfacePattern)

static const t_config_enum_values s_keys_map_SeamPosition {
        {"random",    spRandom},
        {"allrandom", spAllRandom},
        {"nearest",   spNearest}, // unused, replaced by cost
        {"cost",      spCost},
        {"aligned", spAligned},
        {"contiguous", spExtremlyAligned},
        {"rear", spRear},
        {"custom", spCustom}, // for seam object
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SeamPosition);

static const t_config_enum_values s_keys_map_DenseInfillAlgo{
        { "automatic", dfaAutomatic },
        { "autonotfull", dfaAutoNotFull },
        { "autoenlarged", dfaAutoOrEnlarged },
        { "autosmall",  dfaAutoOrNothing},
        { "enlarged", dfaEnlarged },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(DenseInfillAlgo)

static const t_config_enum_values s_keys_map_NoPerimeterUnsupportedAlgo{
        { "none", npuaNone },
        { "noperi", npuaNoPeri },
        { "bridges", npuaBridges },
        { "bridgesoverhangs", npuaBridgesOverhangs },
        { "filled", npuaFilled },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(NoPerimeterUnsupportedAlgo)

static const t_config_enum_values s_keys_map_InfillConnection{
        { "connected", icConnected },
        { "holes", icHoles },
        { "outershell", icOuterShell },
        { "notconnected", icNotConnected },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(InfillConnection)

static const t_config_enum_values s_keys_map_RemainingTimeType{
    { "m117", rtM117 },
    { "m73", rtM73 },
    { "m73q", rtM73_Quiet },
    { "m73m117", rtM73_M117 },
    { "none", rtNone },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(RemainingTimeType)

static const t_config_enum_values s_keys_map_SupportZDistanceType{
    { "filament", zdFilament },
    { "plane", zdPlane },
    { "none", zdNone },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SupportZDistanceType)

static const t_config_enum_values s_keys_map_SLADisplayOrientation{
    { "landscape",      sladoLandscape},
    { "portrait",       sladoPortrait}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SLADisplayOrientation)

static const t_config_enum_values s_keys_map_SLAPillarConnectionMode{
    {"zigzag",          int(SLAPillarConnectionMode::zigzag)},
    {"cross",           int(SLAPillarConnectionMode::cross)},
    {"dynamic",         int(SLAPillarConnectionMode::dynamic)}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SLAPillarConnectionMode)

static const t_config_enum_values s_keys_map_SLAMaterialSpeed = {
    {"slow",            slamsSlow},
    {"fast",            slamsFast},
    {"high_viscosity",  slamsHighViscosity}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SLAMaterialSpeed);

static inline const t_config_enum_values s_keys_map_SLASupportTreeType = {
    {"default", int(sla::SupportTreeType::Default)},
    {"branching",   int(sla::SupportTreeType::Branching)},
    //TODO: {"organic", int(sla::SupportTreeType::Organic)}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(SLASupportTreeType);

// unused
static const t_config_enum_values s_keys_map_BrimType = {
    {"no_brim",         btNoBrim},
    {"outer_only",      btOuterOnly},
    {"inner_only",      btInnerOnly},
    {"outer_and_inner", btOuterAndInner}
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(BrimType)

static const t_config_enum_values s_keys_map_DraftShield = {
    { "disabled", dsDisabled },
    { "limited",  dsLimited  },
    { "enabled",  dsEnabled  }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(DraftShield)

static const t_config_enum_values s_keys_map_LabelObjectsStyle = {
    { "disabled",  int(LabelObjectsStyle::Disabled)  },
    { "octoprint", int(LabelObjectsStyle::Octoprint) },
    { "firmware",  int(LabelObjectsStyle::Firmware)  },
    { "both",      int(LabelObjectsStyle::Both)},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(LabelObjectsStyle)

static const t_config_enum_values s_keys_map_GCodeThumbnailsFormat = {
    { "PNG", int(GCodeThumbnailsFormat::PNG) },
    { "JPG", int(GCodeThumbnailsFormat::JPG) },
    { "QOI", int(GCodeThumbnailsFormat::QOI) },
    { "BIQU",int(GCodeThumbnailsFormat::BIQU) },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(GCodeThumbnailsFormat)

static const t_config_enum_values s_keys_map_ZLiftTop{
    {"everywhere", zltAll},
    {"onlytop", zltTop},
    {"nottop", zltNotTop},
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(ZLiftTop);

static const t_config_enum_values s_keys_map_ForwardCompatibilitySubstitutionRule{
    { "disable",        ForwardCompatibilitySubstitutionRule::Disable },
    { "enable",         ForwardCompatibilitySubstitutionRule::Enable },
    { "enable_silent",  ForwardCompatibilitySubstitutionRule::EnableSilent },
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(ForwardCompatibilitySubstitutionRule);

static t_config_enum_values s_keys_map_PerimeterGeneratorType {
    { "classic", int(PerimeterGeneratorType::Classic) },
    { "arachne", int(PerimeterGeneratorType::Arachne) }
};
CONFIG_OPTION_ENUM_DEFINE_STATIC_MAPS(PerimeterGeneratorType)

static void assign_printer_technology_to_unknown(t_optiondef_map &options, PrinterTechnology printer_technology)
{
    for (std::pair<const t_config_option_key, ConfigOptionDef> &kvp : options)
        if (kvp.second.printer_technology == ptUnknown)
            kvp.second.printer_technology = printer_technology;
}

// Maximum extruder temperature, bumped to 1500 to support printing of glass.
namespace {
    const int max_temp = 1500;
};

ConfigOption *disable_defaultoption(ConfigOption *option, bool default_is_disabled = true) {
    return option->set_can_be_disabled(default_is_disabled);
}

ConfigOptionVectorBase *disable_defaultoption(ConfigOptionVectorBase *option, bool default_is_disabled = true) {
    return (ConfigOptionVectorBase *)option->set_can_be_disabled(default_is_disabled);
}

PrintConfigDef::PrintConfigDef()
{
    this->init_common_params();
    //assign params that are not already allocated to FFF+SLA (default from slic3rPE)
    assign_printer_technology_to_unknown(this->options, ptFFF | ptSLA);
    this->init_extruder_option_keys();
    this->init_fff_params();
    assign_printer_technology_to_unknown(this->options, ptFFF);
    this->init_sla_params();
    assign_printer_technology_to_unknown(this->options, ptSLA);
    this->init_milling_params();
    assign_printer_technology_to_unknown(this->options, ptMill);
    this->init_laser_params();
    assign_printer_technology_to_unknown(this->options, ptLaser);
    this->finalize();
}

void PrintConfigDef::init_common_params()
{
    ConfigOptionDef* def;

    // version of the settings:
    // 4 letters for the software (SUSI, PRSA, BMBU, ORCA)
    // a '_'
    // version in X.X.X.X
    def = this->add("print_version", coString);
    // defautl to none : only set if loaded. only write our version
    def->set_default_value(new ConfigOptionStringVersion());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("printer_technology", coEnum);
    def->label = L("Printer technology");
    def->tooltip = L("Printer technology");
    def->category = OptionCategory::general;
    def->set_enum<PrinterTechnology>({ "FFF", "SLA" });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<PrinterTechnology>(ptFFF));

    def = this->add("bed_shape", coPoints);
    def->label = L("Bed shape");
    def->category = OptionCategory::general;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionPoints{ Vec2d(0, 0), Vec2d(200, 0), Vec2d(200, 200), Vec2d(0, 200) });

    def = this->add("bed_custom_texture", coString);
    def->label = L("Bed custom texture");
    def->category = OptionCategory::general;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("bed_custom_model", coString);
    def->label = L("Bed custom model");
    def->category = OptionCategory::general;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("thumbnails", coPoints);
    def->label = L("Thumbnails size");
    def->tooltip = L("Picture sizes to be stored into a .gcode / .bgcode and .sl1 / .sl1s files, in the following format: \"XxY/EXT, XxY/EXT, ...\"\n"
                     "Currently supported extensions are PNG, QOI and JPG.");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->max = 2048;
    def->set_default_value(new ConfigOptionPoints{ std::initializer_list<Vec2d>{ Vec2d(0,0), Vec2d(0,0) } });

    def = this->add("thumbnails_color", coString);
    def->label = L("Color");
    def->full_label = L("Thumbnail color");
    def->category = OptionCategory::filament;
    def->tooltip = L("This is the color that will be enforced on objects in the thumbnails.");
    def->gui_type = ConfigOptionDef::GUIType::color;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionString("#018aff"));

    def = this->add("thumbnails_custom_color", coBool);
    def->label = L("Enforce thumbnail color");
    def->tooltip = L("Enforce a specific color on thumbnails."
        " If not enforced, their color will be the one defined by the filament.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("thumbnails_end_file", coBool);
    def->label = L("Print at the end");
    def->tooltip = L("Print the thumbnail code at the end of the gcode file instead of the front."
        "\nBe careful! Most firmwares expect it at the front, so be sure that your firmware support it.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false)); 

    def = this->add("thumbnails_with_bed", coBool);
    def->label = L("Bed on thumbnail");
    def->tooltip = L("Show the bed texture on the thumbnail picture.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("thumbnails_format", coEnum);
    def->label = L("Format of G-code thumbnails");
    def->tooltip = L("Format of G-code thumbnails: PNG for best quality, JPG for smallest size, QOI for low memory firmware");
    def->mode = comExpert | comPrusa;
    def->set_enum<GCodeThumbnailsFormat>({ "PNG", "JPG", "QOI", "BIQU" });
    def->set_default_value(new ConfigOptionEnum<GCodeThumbnailsFormat>(GCodeThumbnailsFormat::PNG));

    def          = this->add("thumbnails_tag_format", coBool);
    def->label   = L("Write the thumbnail type in gcode.");
    def->tooltip = L("instead of writing 'thumbnails' as tag in the gcode, it will write 'thumbnails_PNG', thumbnails_JPG', 'thumbnail_QOI', etc.."
        "\n Some firmware need it to know how to decode the thumbnail, some others don't support it.");
    def->mode    = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("thumbnails_with_support", coBool);
    def->label = L("Support on thumbnail");
    def->tooltip = L("Show the supports (and pads) on the thumbnail picture.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("layer_height", coFloat);
    def->label = L("Base Layer height");
    def->category = OptionCategory::slicing;
    def->tooltip = L("This setting controls the height (and thus the total number) of the slices/layers. "
        "Thinner layers give better accuracy but take more time to print.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.2));

    def = this->add("max_print_height", coFloat);
    def->label = L("Max print height");
    def->category = OptionCategory::general;
    def->tooltip = L("Set this to the maximum height that can be reached by your extruder while printing.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(200.0));

    def = this->add("slice_closing_radius", coFloat);
    def->label = L("Slice gap closing radius");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Cracks smaller than 2x gap closing radius are being filled during the triangle mesh slicing. "
        "The gap closing operation may reduce the final print resolution, therefore it is advisable to keep the value reasonably low.");
    def->sidetext = L("mm");
    def->min = 0;
    def->precision = 8;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.049));

    def = this->add("print_host", coString);
    def->label = L("Hostname, IP or URL");
    def->category = OptionCategory::general;
    def->tooltip = L("Slic3r can upload G-code files to a printer host. This field should contain "
                   "the hostname, IP address or URL of the printer host instance. "
                   "Print host behind HAProxy with basic auth enabled can be accessed by putting the user name and password into the URL "
                   "in the following format: https://username:password@your-octopi-address/");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("printhost_apikey", coString);
    def->label = L("API Key / Password");
    def->category = OptionCategory::general;
    def->tooltip = L("Slic3r can upload G-code files to a printer host. This field should contain "
                   "the API Key or the password required for authentication.");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));
    
    // for repetier
    def = this->add("printhost_port", coString);
    def->label = L("Printer");
    def->tooltip = L("Name of the printer");
    def->gui_type = ConfigOptionDef::GUIType::select_close;
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));
    
    // only if there isn't a native SSL support
    def = this->add("printhost_cafile", coString);
    def->label = L("HTTPS CA File");
    def->category = OptionCategory::general;
    def->tooltip = L("Custom CA certificate file can be specified for HTTPS OctoPrint connections, in crt/pem format. "
                   "If left blank, the default OS CA certificate repository is used.");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("printhost_client_cert", coString);
    def->label = L("Client Certificate File");
    def->category = OptionCategory::general;
    def->tooltip = L("A Client certificate file for use with 2-way ssl authentication, in p12/pfx format. "
                   "If left blank, no client certificate is used.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("printhost_client_cert_password", coString);
    def->label = L("Client Certificate Password");
    def->category = OptionCategory::general;
    def->tooltip = L("Password for client certificate for 2-way ssl authentication. "
                   "Leave blank if no password is needed");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionString(""));

    // For PrusaLink
    def = this->add("printhost_user", coString);
    def->label = L("User");
//    def->tooltip = L("");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));

    // For PrusaLink
    def = this->add("printhost_password", coString);
    def->label = L("Password");
//    def->tooltip = L("");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionString(""));

    // Only available on Windows.
    def = this->add("printhost_ssl_ignore_revoke", coBool);
    def->label = L("Ignore HTTPS certificate revocation checks");
    def->tooltip = L("Ignore HTTPS certificate revocation checks in case of missing or offline distribution points. "
        "One may want to enable this option for self signed certificates if connection fails.");
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionBool(false));
    
    def = this->add("preset_names", coStrings);
    def->label = L("Printer preset names");
    def->tooltip = L("Names of presets related to the physical printer");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionStrings{});

    // For PrusaLink
    def = this->add("printhost_authorization_type", coEnum);
    def->label = L("Authorization Type");
//    def->tooltip = L("");
    def->set_enum<AuthorizationType>({
        { "key", L("API key") },
        { "user", L("HTTP digest") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionEnum<AuthorizationType>(atKeyPassword));

    // temporary workaround for compatibility with older Slicer
    {
        def = this->add("preset_name", coString);
        def->set_default_value(new ConfigOptionString());
    }
}

void PrintConfigDef::init_fff_params()
{
    // Maximum extruder temperature, bumped to 1500 to support printing of glass.
    const int max_temp = 1500;

    ConfigOptionDef* def;

    def = this->add("allow_empty_layers", coBool);
    def->label = L("Allow empty layers");
    def->full_label = L("Allow empty layers");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Prevent the gcode builder from triggering an exception if a full layer is empty, and allow the print to start from thin air afterward.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("arc_fitting", coEnum);
    def->label = L("Arc fitting");
    def->category = OptionCategory::firmware;
    def->tooltip = L("Enable to get a G-code file which has G2 and G3 moves. "
                     "G-code resolution will be used as the fitting tolerance.");
    def->set_enum<ArcFittingType>({
        { "disabled",       "Disabled" },
        { "emit_center",    "Enabled: G2/3 I J (ArcWelder)" },
        { "bambu",       "Enabled: G2/3 I J (Bambu)" },
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<ArcFittingType>(ArcFittingType::Disabled));


    def = this->add("arc_fitting_tolerance", coFloatOrPercent);
    def->label = L("Arc fitting tolerance");
    def->sidetext = L("mm or %");
    def->category = OptionCategory::firmware;
    def->tooltip = L("When using the arc_fitting option, allow the curve to deviate a cetain % from the collection of strait paths."
        "\nCan be a mm value or a percentage of the current extrusion width.");
    def->mode = comAdvancedE | comSuSi;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloatOrPercent(5, true));

    def = this->add("avoid_crossing_curled_overhangs", coBool);
    def->label = L("Avoid crossing curled overhangs (Experimental)");
    def->category = OptionCategory::perimeter;
    // TRN PrintSettings: "Avoid crossing curled overhangs (Experimental)"
    def->tooltip = L("Plan travel moves such that the extruder avoids areas where the filament may be curled up. "
                   "This is mostly happening on steeper rounded overhangs and may cause a crash with the nozzle. "
                   "This feature slows down both the print and the G-code generation.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("avoid_crossing_perimeters", coBool);
    def->label = L("Avoid crossing perimeters");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Optimize travel moves in order to minimize the crossing of perimeters. "
        "This is mostly useful with Bowden extruders which suffer from oozing. "
        "This feature slows down both the print and the G-code generation.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("avoid_crossing_not_first_layer", coBool);
    def->label = L("Don't avoid crossing on 1st layer");
    def->full_label = L("Don't avoid crossing on 1st layer");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Disable 'Avoid crossing perimeters' for the first layer.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("avoid_crossing_perimeters_max_detour", coFloatOrPercent);
    def->label = L("Avoid crossing perimeters - Max detour length");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("The maximum detour length for avoid crossing perimeters. "
                     "If the detour is longer than this value, avoid crossing perimeters is not applied for this travel path. "
                     "Detour length can be specified either as an absolute value or as percentage (for example 50%) of a direct travel path.");
    def->sidetext = L("mm or % (zero to disable)");
    def->min = 0;
    def->max_literal = { 1000, false };
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0., false));

    def = this->add("avoid_crossing_top", coBool);
    def->label = L("Avoid top surface for travels");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("When using 'Avoid crossing perimeters', consider the top surfaces as a void, to avoid travelling over them if possible.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("bed_temperature", coInts);
    def->label = L("Other layers");
    def->category = OptionCategory::filament;
    def->tooltip = L("Bed temperature for layers after the first one. "
                   "Set zero to disable bed temperature control commands in the output.");
    def->sidetext = L("°C");
    def->full_label = L("Bed temperature");
    def->sidetext = L("°C");
    def->min = 0;
    def->max = 300;
    def->is_vector_extruder = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInts { 0 });

    def = this->add("before_layer_gcode", coString);
    def->label = L("Before layer change G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This custom code is inserted at every layer change, right before the Z move. "
                   "Note that you can use placeholder variables for all Slic3r settings as well "
                   "as {layer_num} and {layer_z}.");
    def->multiline = true;
    def->full_width = true;
    def->height = 5;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("between_objects_gcode", coString);
    def->label = L("Between objects G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This code is inserted between objects when using sequential printing. By default extruder and bed temperature are reset using non-wait command; however if M104, M109, M140 or M190 are detected in this custom code, Slic3r will not add temperature commands. Note that you can use placeholder variables for all Slic3r settings, so you can put a \"M109 S{first_layer_temperature}\" command wherever you want.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("bottom_solid_layers", coInt);
    //TRN Print Settings: "Bottom solid layers"
    def->label = L_CONTEXT("Bottom", "Layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Number of solid layers to generate on bottom surfaces.");
    def->full_label = L("Bottom solid layers");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(3));

    def = this->add("bottom_solid_min_thickness", coFloat);
    def->label = L_CONTEXT("Bottom", "Layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("The number of bottom solid layers is increased above bottom_solid_layers if necessary to satisfy "
    				 "minimum thickness of bottom shell.");
    def->full_label = L("Minimum bottom shell thickness");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("bridge_acceleration", coFloatOrPercent);
    def->label = L("Bridges");
    def->full_label = L("Bridge acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for bridges."
                "\nCan be a % of the default acceleration"
                "\nSet zero to use default acceleration for bridges.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "default_acceleration";
    def->min = 0;
    def->max_literal = { -220, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("bridge_angle", coFloat);
    def->label = L("Bridging");
    def->full_label = L("Bridging angle");
    def->category = OptionCategory::infill;
    def->tooltip = L("Bridging angle override."
                   "\nIf disabled, the bridging angle will be calculated automatically."
                   " Otherwise the provided angle will be used for all bridges."
                   "Note: 180° is the same as zero angle.");
    def->sidetext = L("°");
    def->min = 0;
    def->can_be_disabled = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloat(0.), false));

    def = this->add("bridge_fan_speed", coInts);
    def->label = L("Bridges fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during bridges and overhangs. It won't slow down the fan if it's currently running at a higher speed."
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts{ 100 }, false));

    def = this->add("bridge_type", coEnum);
    def->label = L("Bridge flow baseline");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allow you to choose the base for the bridge flow compute, the result will be multiplied by the bridge flow to have the final result."
        "\nA bridge is an extrusion with nothing under it to flatten it, and so it can't have a 'rectangle' shape but a circle one."
        "\n * The default way to compute a bridge flow is to use the nozzle diameter as the diameter of the extrusion cross-section. It shouldn't be higher than that to prevent sagging."
        "\n * A second way to compute a bridge flow is to use the current layer height, so it shouldn't protrude below it. Note that may create too thin extrusions and so a bad bridge quality."
        "\n * A Third way to compute a bridge flow is to continue to use the current flow/section (mm3 per mm). If there is no current flow, it will use the solid infill one."
        " To use if you have some difficulties with the big flow changes from perimeter and infill flow to bridge flow and vice-versa, the bridge flow ratio let you compensate for the change in speed."
        " \nThe preview will display the expected shape of the bridge extrusion (cylinder), don't expect a magical thick and solid air to flatten the extrusion magically.");
    def->sidetext = L("%");
    def->set_enum<BridgeType>({
        { "nozzle", L("Nozzle diameter") },
        { "height", L("Layer height") },
        { "flow", L("Keep current flow") },
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<BridgeType>{ BridgeType::btFromNozzle });

    def = this->add("bridge_flow_ratio", coPercent);
    def->label = L("Bridge");
    def->full_label = L("Bridge flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("This factor affects the amount of plastic for bridging. "
                   "You can decrease it slightly to pull the extrudates and prevent sagging, "
                   "although default settings are usually good and you should experiment "
                   "with cooling (use a fan) before tweaking this."
                   "\nFor reference, the default bridge flow is (in mm3/mm): (nozzle diameter) * (nozzle diameter) * PI/4");
    def->min = 2;
    def->max = 1000;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("over_bridge_flow_ratio", coPercent);
    def->label = L("Above the bridges");
    def->full_label = L("Above bridge flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("Flow ratio to compensate for the gaps in a bridged top surface. Used for ironing infill"
        "pattern to prevent regions where the low-flow pass does not provide a smooth surface due to a lack of plastic."
        " You can increase it slightly to pull the top layer at the correct height. Recommended maximum: 120%.");
    def->min = 2;
    def->max = 1000;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("bridge_precision", coFloatOrPercent);
    def->label = L("Bridge precision");
    def->category = OptionCategory::slicing;
    def->tooltip = L("This is the precision of the bridge detection. If you put it too low, the bridge detection will be very inneficient."
                    "\nCan be a % of the bridge spacing.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(25, true));

    def = this->add("bridge_overlap_min", coPercent);
    def->label = L("Min");
    def->full_label = L("Min bridge density");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("Minimum density for bridge lines. If Lower than bridge_overlap, then the overlap value can be lowered automatically down to this value."
        " If the value is higher, this parameter has no effect."
        "\nDefault to 87.5% to allow a little void between the lines.");
    def->min = 2;
    def->max = 2000;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(80));

    def = this->add("bridge_overlap", coPercent);
    def->label = L("Max");
    def->full_label = L("Max bridge density");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("Maximum density for bridge lines. If you want more space between line (or less), you can modify it."
        " A value of 50% will create two times less lines, and a value of 200% will create two time more lines that overlap each other.");
    def->min = 2;
    def->max = 2000;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(90));

    def = this->add("bridge_speed", coFloatOrPercent);
    def->label = L("Bridges");
    def->full_label = L("Bridge speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing bridges."
        "\nThis can be expressed as a percentage (for example: 60%) over the Default speed."
        "\nSet zero to use the autospeed for this feature");
    def->sidetext = L("mm/s or %");
    def->aliases = { "bridge_feed_rate" };
    def->ratio_over = "default_speed";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(60, true));

    def = this->add("brim_inside_holes", coBool);
    def->label = L("Brim inside holes");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Allow to create a brim over an island when it's inside a hole (or surrounded by an object)."
        "\nIncompatible with brim_width_interior, as it enables it with brim_width width.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("brim_per_object", coBool);
    def->label = L("Brim per object");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Create a brim per object instead of a brim for the plater."
        " Useful for complete_object or if you have your brim detaching before printing the object."
        "\nBe aware that the brim may be truncated if objects are too close together..");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("brim_width", coFloat);
    def->label = L("Brim width");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Horizontal width of the brim that will be printed around each object on the first layer."
        "\nWhen raft is used, no brim is generated (use raft_first_layer_expansion).");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("brim_width_interior", coFloat);
    def->label = L("Interior Brim width");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Horizontal width of the brim that will be printed inside each object on the first layer.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("brim_ears", coBool);
    def->label = L("Brim ears");
    def->full_label = L("Brim ears");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Only draw brim over the sharp edges of the model.");
    def->mode = comSimpleAE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("brim_ears_max_angle", coFloat);
    def->label = L("Max angle");
    def->full_label = L("Brim ear max angle");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Maximum angle to let a brim ear appear. \nIf set to 0, no brim will be created. \nIf set to ~178, brim will be created on everything but straight sections.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 180;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(125)); 
    
    def = this->add("brim_acceleration", coFloatOrPercent);
    def->label = L("Brim & Skirt");
    def->full_label = L("Brim & Skirt acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for brim and skirt. "
        "\nCan be a % of the support acceleration"
        "\nSet zero to use support acceleration.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "support_material_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("brim_ears_detection_length", coFloat);
    def->label = L("Detection radius");
    def->full_label = L("Brim ear detection length");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("The geometry will be decimated before dectecting sharp angles. This parameter indicates the minimum length of the deviation for the decimation."
                    "\n0 to deactivate");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(1));

    def = this->add("brim_ears_pattern", coEnum);
    def->label = L("Pattern");
    def->full_label = L("Ear pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Pattern for the ear. The concentric is the default one."
                    " The rectilinear has a perimeter around it, you can try it if the concentric has too many problems to stick to the build plate.");
    def->set_enum<InfillPattern>({
        { "concentric", L("Concentric") },
        { "rectilinear", L("Rectilinear") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipConcentric));

    def = this->add("brim_separation", coFloat);
    def->label = L("Brim separation gap");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Offset of brim from the printed object. Should be kept at 0 unless you encounter great difficulties to separate them."
        "\nIt's subtracted to brim_width and brim_width_interior, so it has to be lower than them. The offset is applied after the first layer XY compensation (elephant foot).");
    def->sidetext = L("mm");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));
    def->aliases = { "brim_offset" }; // from superslicer 2.3

    def = this->add("brim_speed", coFloatOrPercent);
    def->label = L("Brim & Skirt");
    def->full_label = L("Brim & Skirt speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("This separate setting will affect the speed of brim and skirt. "
        "\nIf expressed as percentage (for example: 80%) it will be calculated over the Support speed setting."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "support_material_speed";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

#if 0
    def = this->add("brim_type", coEnum);
    def->label = L("Brim type");
    def->category = L("Skirt and brim");
    def->tooltip = L("The places where the brim will be printed around each object on the first layer.");
    def->set_enum<BrimType>({
        { "no_brim",         L("No brim") },
        { "outer_only",      L("Outer brim only") },
        { "inner_only",      L("Inner brim only") },
        { "outer_and_inner", L("Outer and inner brim") } 
    });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<BrimType>(btOuterOnly));
#endif

    def = this->add("chamber_temperature", coInts);
    def->label = L("Chamber");
    def->full_label = L("Chamber temperature");
    def->category = OptionCategory::cooling;
    def->tooltip = L("Chamber temperature. Note that this setting doesn't do anything, but you can access it in Start G-code, Tool change G-code and the other ones, like for other temperature settings.");
    def->sidetext = L("°C");
    def->min = 0;
    def->max = 300;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts{ 0 });

    def = this->add("colorprint_heights", coFloats);
    def->label = L("Colorprint height");
    def->category = OptionCategory::slicing;
    def->mode = comExpert | comPrusa; // note: hidden setting
    def->tooltip = L("Heights at which a filament change is to occur. ");
    def->set_default_value(new ConfigOptionFloats { });

    def = this->add("compatible_printers", coStrings);
    def->label = L("Compatible printers");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionStrings());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("compatible_printers_condition", coString);
    def->label = L("Compatible printers condition");
    def->tooltip = L("A boolean expression using the configuration values of an active printer profile. "
                   "If this expression evaluates to true, this profile is considered compatible "
                   "with the active printer profile.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("compatible_prints", coStrings);
    def->label = L("Compatible print profiles");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionStrings());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("compatible_prints_condition", coString);
    def->label = L("Compatible print profiles condition");
    def->tooltip = L("A boolean expression using the configuration values of an active print profile. "
                   "If this expression evaluates to true, this profile is considered compatible "
                   "with the active print profile.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    // The following value is to be stored into the project file (AMF, 3MF, Config ...)
    // and it contains a sum of "compatible_printers_condition" values over the print and filament profiles.
    def = this->add("compatible_printers_condition_cummulative", coStrings);
    def->set_default_value(new ConfigOptionStrings());
    def->cli = ConfigOptionDef::nocli;
    def->mode = comNone | comPrusa; // note: hidden setting
    def = this->add("compatible_prints_condition_cummulative", coStrings);
    def->set_default_value(new ConfigOptionStrings());
    def->cli = ConfigOptionDef::nocli;
    def->mode = comNone | comPrusa; // note: hidden setting

    def = this->add("complete_objects", coBool);
    def->label = L("Complete individual objects");
    def->category = OptionCategory::output;
    def->tooltip = L("When printing multiple objects or copies, this feature will complete "
        "each object before moving onto next one (and starting it from its bottom layer). "
        "This feature is useful to avoid the risk of ruined prints. "
        "Slic3r should warn and prevent you from extruder collisions, but beware.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("parallel_objects_step", coFloat);
    def->label = L("Parallel printing step");
    def->category = OptionCategory::output;
    def->tooltip = L("When multiple objects are present, instead of jumping form one to another at each layer"
        " the printer will continue to print the current object layers up to this height before moving to the next object."
        " (first layers will be still printed one by one)."
        "\nThis feature also use the same extruder clearance radius field as 'complete individual objects' (complete_objects)"
        ", but you can modify them to instead reflect the clerance of the nozzle, if this field reflect the z-clearance of it."
        "\nThis field is exclusive with 'complete individual objects' (complete_objects). Set to 0 to deactivate.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("complete_objects_one_skirt", coBool);
    def->label = L("Allow only one skirt loop");
    def->category = OptionCategory::output;
    def->tooltip = L("When using 'Complete individual objects', the default behavior is to draw the skirt around each object."
        " if you prefer to have only one skirt for the whole platter, use this option.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("complete_objects_sort", coEnum);
    def->label = L("Object sort");
    def->category = OptionCategory::output;
    def->tooltip = L("When printing multiple objects or copies on after another, this will help you to choose how it's ordered."
        "\nObject will sort them by the order of the right panel."
        "\nLowest Y will sort them by their lowest Y point. Useful for printers with a X-bar."
        "\nLowest Z will sort them by their height, useful for delta printers.");
    def->mode = comAdvancedE | comSuSi;
    def->set_enum<CompleteObjectSort>({
        { "object", L("Right panel") },
        { "lowy", L("lowest Y") },
        { "lowz", L("lowest Z") },
    });
    def->set_default_value(new ConfigOptionEnum<CompleteObjectSort>(cosObject));

#if 0
    //not used anymore, to remove !! @DEPRECATED
    def = this->add("cooling", coBools);
    def->label = L("Enable auto cooling");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This flag enables the automatic cooling logic that adjusts print speed "
                   "and fan speed according to layer printing time.");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { true });
#endif

    def = this->add("cooling_tube_retraction", coFloat);
    def->label = L("Cooling tube position");
    def->category = OptionCategory::mmsetup;
    def->tooltip = L("Distance of the center-point of the cooling tube from the extruder tip.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(91.5));

    def = this->add("cooling_tube_length", coFloat);
    def->label = L("Cooling tube length");
    def->category = OptionCategory::mmsetup;
    def->tooltip = L("Length of the cooling tube to limit space for cooling moves inside it.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(5.));

    def = this->add("curve_smoothing_angle_convex", coFloat);
    def->label = L("Min convex angle");
    def->full_label = L("Curve smoothing minimum angle (convex)");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Minimum (convex) angle at a vertex to enable smoothing"
        " (trying to create a curve around the vertex). "
        "180 : nothing will be smooth, 0 : all angles will be smoothened.");
    def->sidetext = L("°");
    def->aliases = { "curve_smoothing_angle" };
    def->cli = "curve-smoothing-angle-convex=f";
    def->min = 0;
    def->max = 180;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("curve_smoothing_angle_concave", coFloat);
    def->label = L("Min concave angle");
    def->full_label = L("Curve smoothing minimum angle (concave)");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Minimum (concave) angle at a vertex to enable smoothing"
        " (trying to create a curve around the vertex). "
        "180 : nothing will be smooth, 0 : all angles will be smoothened.");
    def->sidetext = L("°");
    def->cli = "curve-smoothing-angle-concave=f";
    def->min = 0;
    def->max = 180;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("curve_smoothing_precision", coFloat);
    def->label = L("Precision");
    def->full_label = L("Curve smoothing precision");
    def->category = OptionCategory::slicing;
    def->tooltip = L("These parameters allow the slicer to smooth the angles in each layer. "
        "The precision will be at least the new precision of the curve. Set to 0 to deactivate."
        "\nNote: as it uses the polygon's edges and only works in the 2D planes, "
        "you must have a very clean or hand-made 3D model."
        "\nIt's really only useful to smoothen functional models or very wide angles.");
    def->sidetext = L("mm");
    def->min = 0;
    def->precision = 8;
    def->cli = "curve-smoothing-precision=f";
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("curve_smoothing_cutoff_dist", coFloat);
    def->label = L("cutoff");
    def->full_label = L("Curve smoothing cutoff dist");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Maximum distance between two points to allow adding new ones. Allow to avoid distorting long strait areas.\nSet zero to disable.");
    def->sidetext = L("mm");
    def->min = 0;
    def->cli = "curve-smoothing-cutoff-dist=f";
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(2));

    def = this->add("default_acceleration", coFloatOrPercent);
    def->label = L("Default");
    def->category = OptionCategory::speed;
    def->full_label = L("Default acceleration");
    def->tooltip = L("This is the acceleration your printer will be reset to after "
        "the role-specific acceleration values are used (perimeter/infill). "
        "\nAccelerations from the left column can also be expressed as a percentage of this value."
        "\nThis can be expressed as a percentage (for example: 80%) over the machine Max Acceleration for X axis."
        "\nSet zero to prevent resetting acceleration at all.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "machine_max_acceleration_x";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("default_filament_profile", coStrings);
    def->label = L("Default filament profile");
    def->tooltip = L("Default filament profile associated with the current printer profile. "
                   "On selection of the current printer profile, this filament profile will be activated.");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionStrings());
    def->cli = ConfigOptionDef::nocli;

    def           = this->add("default_fan_speed", coInts);
    def->label    = L("Default fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip  = L(
        "Default speed for the fan, to set the speed for features where there is no fan control. Useful for PLA and other low-temp filament."
        "\nSet 0 to disable the fan by default. Useful for ABS and other high-temp filaments."
        "\nIf disabled, no fan speed command will be emmited when possible (if a feature set a speed, it won't be reverted).");
    def->mode               = comSimpleAE | comSuSi;
    def->min                = 0;
    def->max                = 100;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));
    def->aliases = { "min_fan_speed" }; // only if "fan_always_on"

    def = this->add("default_print_profile", coString);
    def->label = L("Default print profile");
    def->tooltip = L("Default print profile associated with the current printer profile. "
                   "On selection of the current printer profile, this print profile will be activated.");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("default_speed", coFloatOrPercent);
    def->label = L("Default");
    def->category = OptionCategory::speed;
    def->full_label = L("Default speed");
    def->tooltip = L("This is the reference speed that other 'main' speed can reference to by a %."
        "\nThis setting doesn't do anything by itself, and so is deactivated unless a speed depends on it (a % from the left column)."
        "\nThis can be expressed as a percentage (for example: 80%) over the machine Max Feedrate for X axis."
        "\nSet zero to use autospeed for speed fields using a % of this setting.");
    def->sidetext = L("mm/s for %-based speed");
    def->sidetext_width = 40;
    def->ratio_over = "machine_max_feedrate_x";
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, false));

    def = this->add("disable_fan_first_layers", coInts);
    def->label = L("Disable fan for the first");
    def->category = OptionCategory::cooling;
    def->tooltip = L("You can set this to a positive value to disable fan at all "
                   "during the first layers, so that it does not make adhesion worse.");
    def->sidetext = L("layers");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 1 });

    def = this->add("dont_support_bridges", coBool);
    def->label = L("Don't support bridges");
    def->category = OptionCategory::support;
    def->tooltip = L("Experimental option for preventing support material from being generated "
                   "under bridged areas.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("draft_shield", coEnum);
    def->label = L("Draft shield");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("With draft shield active, the skirt will be printed skirt_distance from the object, possibly intersecting brim.\n"
        "Enabled = skirt is as tall as the highest printed object.\n"
        "Limited = skirt is as tall as specified by skirt_height.\n"
        "This is useful to protect an ABS or ASA print from warping and detaching from print bed due to wind draft.");
    def->set_enum<DraftShield>({
        { "disabled", L("Disabled") },
        { "limited", L("Limited") },
        { "enabled", L("Enabled") },
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<DraftShield>(dsDisabled));

    def = this->add("duplicate_distance", coFloat);
    def->label = L("Default distance between objects");
    def->category = OptionCategory::output;
    def->tooltip = L("Default distance used for the auto-arrange feature of the platter.\nSet to 0 to use the last value instead.");
    def->sidetext = L("mm");
    def->aliases = { "multiply_distance" };
    def->min = 0;
    def->mode = comExpert | comPrusa | comSuSi;
    def->set_default_value(new ConfigOptionFloat(6));

    def = this->add("end_gcode", coString);
    def->label = L("End G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This end procedure is inserted at the end of the output file. "
                   "Note that you can use placeholder variables for all Slic3r settings.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString("M104 S0 ; turn off temperature\nG28 X0  ; home X axis\nM84     ; disable motors\n"));

    def = this->add("end_filament_gcode", coStrings);
    def->label = L("End G-code");
    def->full_label = L("Filament end G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This end procedure is inserted at the end of the output file, before the printer end gcode (and "
                   "before any toolchange from this filament in case of multimaterial printers). "
                   "Note that you can use placeholder variables for all Slic3r settings. "
                   "If you have multiple extruders, the gcode is processed in extruder order.");
    def->multiline = true;
    def->full_width = true;
    def->height = 120;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings { "; Filament-specific end gcode \n;END gcode for filament\n" });

    def = this->add("top_fill_pattern", coEnum);
    def->label = L("Top");
    def->full_label = L("Top fill Pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Fill pattern for top infill. This only affects the top visible layer, and not its adjacent solid shells."
        "\nIf you want an 'aligned' pattern, set 90° to the fill angle increment setting.");
    def->cli = "top-fill-pattern|external-fill-pattern=s";
    def->aliases = { "external_fill_pattern" };
    def->set_enum<InfillPattern>({
        { "rectilinear",        L("Rectilinear") },
        { "monotonic",          L("Monotonic") },
        { "monotonicgapfill",   L("Monotonic (filled)") },
        { "monotoniclines",     L("Monotonic Lines") },
        { "alignedrectilinear", L("Aligned Rectilinear") },
        { "concentric",         L("Concentric") },
        { "concentricgapfill",  L("Concentric (filled)") },
        { "hilbertcurve",       L("Hilbert Curve") },
        { "archimedeanchords",  L("Archimedean Chords") },
        { "octagramspiral",     L("Octagram Spiral") },
        { "sawtooth",     L("Sawtooth") },
        { "smooth",     L("Ironing") },
    });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipMonotonic));

    def = this->add("bottom_fill_pattern", coEnum);
    def->label = L("Bottom");
    def->full_label = L("Bottom fill pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Fill pattern for bottom infill. This only affects the bottom visible layer, and not its adjacent solid shells."
        "\nIf you want an 'aligned' pattern, set 90° to the fill angle increment setting.");
    def->cli = "bottom-fill-pattern|external-fill-pattern=s";
    def->aliases = { "external_fill_pattern" };
    def->set_enum<InfillPattern>({
        { "rectilinear",        L("Rectilinear") },
        { "monotonic",          L("Monotonic") },
        { "monotonicgapfill",   L("Monotonic (filled)") },
        { "monotoniclines",     L("Monotonic Lines") },
        { "alignedrectilinear", L("Aligned Rectilinear") },
        { "concentric",         L("Concentric") },
        { "concentricgapfill",  L("Concentric (filled)") },
        { "hilbertcurve",       L("Hilbert Curve") },
        { "archimedeanchords",  L("Archimedean Chords") },
        { "octagramspiral",     L("Octagram Spiral") },
        { "smooth",     L("Ironing") },
    });

    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipMonotonic));

    def = this->add("solid_fill_pattern", coEnum);
    def->label = L("Solid fill pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Fill pattern for solid (internal) infill. This only affects the solid not-visible layers. You should use rectilinear in most cases. You can try ironing for translucent material."
        " Rectilinear (filled) replaces zig-zag patterns by a single big line & is more efficient for filling little spaces."
        "\nIf you want an 'aligned' pattern, set 90° to the fill angle increment setting.");
    def->set_enum<InfillPattern>({
        { "rectilinear",        L("Rectilinear") },
        { "rectilineargapfill", L("Rectilinear (filled)") },
        { "monotonic",          L("Monotonic") },
        { "monotonicgapfill",   L("Monotonic (filled)") },
        { "monotoniclines",     L("Monotonic Lines") },
        { "alignedrectilinear", L("Aligned Rectilinear") },
        { "concentric",         L("Concentric") },
        { "concentricgapfill",  L("Concentric (filled)") },
        { "hilbertcurve",       L("Hilbert Curve") },
        { "archimedeanchords",  L("Archimedean Chords") },
        { "octagramspiral",     L("Octagram Spiral") },
        { "smooth",             L("Ironing") },
    });

    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipRectilinearWGapFill));

    def = this->add("bridge_fill_pattern", coEnum);
    def->label = L("Bridging fill pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Fill pattern for bridges and internal bridge infill.");
    def->set_enum<InfillPattern>({
        { "rectilinear", L("Rectilinear") },
        { "monotonic", L("Monotonic") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipRectilinear));

    def = this->add("enforce_full_fill_volume", coBool);
    def->label = L("Enforce 100% fill volume");
    def->category = OptionCategory::infill;
    def->tooltip = L("Experimental option which modifies (in solid infill) fill flow to have the exact amount of plastic inside the volume to fill "
        "(it generally changes the flow from -7% to +4%, depending on the size of the surface to fill and the overlap parameters, "
        "but it can go as high as +50% for infill in very small areas where rectilinear doesn't have good coverage). It has the advantage "
        "to remove the over-extrusion seen in thin infill areas, from the overlap ratio");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("enforce_retract_first_layer", coBool);
    def->label = L("But on first layer");
    def->full_label = L("Don't check crossings for retraction on first layer");
    def->category = OptionCategory::extruders;
    def->tooltip = L("let the retraction happens on the first layer even if the travel path does not exceed the upper layer's perimeters.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("external_infill_margin", coFloatOrPercent);
    def->label = L("Default");
    def->full_label = L("Default infill margin");
    def->category = OptionCategory::infill;
    def->tooltip = L("This parameter grows the top/bottom/solid layers by the specified mm to anchor them into the sparse infill and support the perimeters above."
        " Put 0 to deactivate it. Can be a % of the width of the perimeters.");
    def->sidetext = L("mm or %");
    def->ratio_over = "perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 50, true };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(150, true));

    def = this->add("bridged_infill_margin", coFloatOrPercent);
    def->label = L("Bridged");
    def->full_label = L("Bridge margin");
    def->category = OptionCategory::infill;
    def->tooltip = L("This parameter grows the bridged solid infill layers by the specified mm to anchor them into the sparse infill and over the perimeters below. Put 0 to deactivate it. Can be a % of the width of the external perimeter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "external_perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 50, true };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(200, true));

    def = this->add("external_perimeter_extrusion_width", coFloatOrPercent);
    def->label = L("External perimeters");
    def->full_label = L("External perimeters width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for external perimeters. "
        "If left zero, default extrusion width will be used if set, otherwise 1.05 x nozzle diameter will be used. "
        "If expressed as percentage (for example 112.5%), it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(105, true));

    def = this->add("external_perimeter_extrusion_spacing", coFloatOrPercent);
    def->label = L("External perimeters");
    def->full_label = L("External perimeters spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like the External perimeters width, but this value is the distance between the edge and the 'frontier' to the next perimeter."
                "\nSetting the spacing will deactivate the width setting, and vice versa."
                "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("external_perimeter_extrusion_change_odd_layers", coFloatOrPercent);
    def->label = L("External perimeters");
    def->full_label = L("External perimeters spacing change on even layers");
    def->category = OptionCategory::width;
    def->tooltip = L("Change width on every even layer (and not on odd layers like the first one) for better overlap with adjacent layers and getting stringer shells. "
                     "Try values about +/- 0.1 with different sign for external and internal perimeters."
                     "\nThis could be combined with extra permeters on even layers."
                     "\nWorks as absolute spacing or a % of the spacing."
                     "\nset 0 to disable");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(false, 0));

    def = this->add("external_perimeter_cut_corners", coPercent);
    def->label = L("Cutting corners");
    def->full_label = L("Ext. peri. cut corners");
    def->category = OptionCategory::width;
    def->tooltip = L("Activate this option to modify the flow to acknowledge that the nozzle is round and the corners will have a round shape, and so change the flow to realize that and avoid over-extrusion."
        " 100% is activated, 0% is deactivated and 50% is half-activated."
        "\nNote: At 100% this changes the flow by ~5% over a very small distance (~nozzle diameter), so it shouldn't be noticeable unless you have a very big nozzle and a very precise printer.");
    def->sidetext = L("%");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(0));

    def = this->add("external_perimeter_fan_speed", coInts);
    def->label = L("External perimeter fan speed");
    def->tooltip = L("When set to a non-zero value this fan speed is used only for external perimeters (visible ones) and thin walls."
                    "\nSet to 0 to stop the fan."
                    "\nIf disabled, the default fan speed will be used."
                    "\nExternal perimeters can benefit from higher fan speed to improve surface finish, "
                    "while internal perimeters, infill, etc. benefit from lower fan speed to improve layer adhesion."
                    "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("external_perimeter_overlap", coPercent);
    def->label = L("external perimeter overlap");
    def->full_label = L("Ext. peri. overlap");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allows you to reduce the overlap between the perimeters and the external one, to reduce the impact of the perimeters' artifacts."
        " 100% means that no gap is left, and 0% means that the external perimeter isn't contributing to the overlap with the 'inner' one.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("external_perimeter_acceleration", coFloatOrPercent);
    def->label = L("External");
    def->full_label = L("External Perimeter acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for external perimeters. "
                "\nCan be a % of the internal perimeter acceleration"
                "\nSet zero to use internal perimeter acceleration for external perimeters.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "perimeter_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("external_perimeter_speed", coFloatOrPercent);
    def->label = L("External");
    def->full_label = L("External perimeters speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("This separate setting will affect the speed of external perimeters (the visible ones). "
                   "\nIf expressed as percentage (for example: 80%) it will be calculated over the Internal Perimeters speed setting."
                   "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "perimeter_speed";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("external_perimeters_first", coBool);
    def->label = L("first");
    def->full_label = L("External perimeters first");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Print contour perimeters from the outermost one to the innermost one "
        "instead of the default inverse order.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("external_perimeters_vase", coBool);
    def->label = L("In vase mode (no seam)");
    def->full_label = L("External perimeters in vase mode");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Print contour perimeters in two circles, in a continuous way, like for a vase mode. It needs the external_perimeters_first parameter to work."
        " \nDoesn't work for the first layer, as it may damage the bed overwise."
        " \nNote that it will use min_layer_height from your hardware setting as the base height (it doesn't start at 0)"
        ", so be sure to put here the lowest value your printer can handle."
        " if it's not lower than two times the current layer height, it falls back to the normal algorithm, as there is not enough room to do two loops.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("external_perimeters_nothole", coBool);
    def->label = L("Only for contours");
    def->full_label = L("Ext peri first for outer side");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Only do the vase trick on the external side. Useful when the thickness is too low.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("external_perimeters_hole", coBool);
    def->label = L("Only for holes");
    def->full_label = L("ext peri first for inner side");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Only do the vase trick on the external side. Useful when you only want to remove seam from screw hole.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("extra_perimeters", coBool);
    def->label = L("filling horizontal gaps on slopes");
    def->full_label = L("Extra perimeters (do nothing)");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Add more perimeters when needed for avoiding gaps in sloping walls. "
        "Slic3r keeps adding perimeters, until more than 70% of the loop immediately above "
        "is supported."
        "\nIf you succeed in triggering the algorithm behind this setting, please send me a message."
        " Personally, I think it's useless.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("extra_perimeters_on_overhangs", coBool);
    def->label = L("Extra perimeters on overhangs (Experimental)");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Detect overhang areas where bridges cannot be anchored, and fill them with "
                    "extra perimeter paths. These paths are anchored to the nearby non-overhang area when possible."
                    "\nIf you use this setting, strongly consider also using overhangs_reverse.");
    def->aliases = {"extra_perimeters_overhangs"};
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("extra_perimeters_odd_layers", coBool);
    def->label = L("On even layers");
    def->full_label = L("Extra perimeter on even layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Add one perimeter every even layer (and not on odd layers like the first one). With this, infill is taken into the sandwich"
        " and you may be able to reduce drastically the infill/perimeter overlap setting. ");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("extruder", coInt);
    def->label = L("Extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use (unless more specific extruder settings are specified). "
        "This value overrides perimeter and infill extruders, but not the support extruders.");
    def->mode = comAdvancedE | comPrusa; // note: hidden setting
    def->min = 0;  // 0 = inherit defaults
    def->set_enum_labels(ConfigOptionDef::GUIType::i_enum_open, 
        { L("default"), "1", "2", "3", "4", "5", "6", "7", "8", "9" }); // override label for item 0


    def = this->add("first_layer_extruder", coInt);
    def->gui_type = ConfigOptionDef::GUIType::i_enum_open;
    def->label = L("First layer extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use (unless more specific extruder settings are specified) for the first layer.");
    def->min = 0;  // 0 = inherit defaults
    def->set_enum_labels(ConfigOptionDef::GUIType::i_enum_open, 
        { L("default"), "1", "2", "3", "4", "5", "6", "7", "8", "9" }); // override label for item 0

    def = this->add("extruder_clearance_height", coFloat);
    def->label = L("Height");
    def->full_label = L("Extruder clearance height");
    def->category = OptionCategory::output;
    def->tooltip = L("Set this to the vertical distance between your nozzle tip and (usually) the X carriage rods. "
                   "In other words, this is the height of the clearance cylinder around your extruder, "
                   "and it represents the maximum depth the extruder can peek before colliding with "
                   "other printed objects."); // TODO: "peek?" is this the correct word?
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(20));

    def = this->add("extruder_clearance_radius", coFloat);
    def->label = L("Radius");
    def->category = OptionCategory::output;
    def->full_label = L("Extruder clearance radius");
    def->tooltip = L("Set this to the clearance radius around your extruder. "
                   "If the extruder is not centered, choose the largest value for safety. "
                   "This setting is used to check for collisions and to display the graphical preview "
                   "in the platter."
                   "\nSet zero to disable clearance checking.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(20));

    def = this->add("extruder_colour", coStrings);
    def->label = L("Extruder Color");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This is only used in Slic3r interface as a visual help.");
    def->gui_type = ConfigOptionDef::GUIType::color;
    // Empty string means no color assigned yet.
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings{ "" });

    def = this->add("extruder_extrusion_multiplier_speed", coGraphs);
    def->label = L("Extrusion multipler");
    def->tooltip = L("This string is edited by a Dialog and contains extusion multiplier for different speeds.");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionGraphs( GraphData(0,10, GraphData::GraphType::LINEAR,
        {{10,1.},{20,1.},{30,1.},{40,1.},{60,1.},{80,1.},{120,1.},{160,1.},{240,1.},{320,1.},{480,1.},{640,1.},{960,1.},{1280,1.}}
    )));
    def->graph_settings = std::make_shared<GraphSettings>();
    def->graph_settings->title       = L("Extrusion multiplier per extrusion speed");
    def->graph_settings->description = L("Choose the extrusion multipler value for multiple speeds.\nYou can add/remove points with a right clic.");
    def->graph_settings->x_label     = L("Print speed (mm/s)");
    def->graph_settings->y_label     = L("Extrusion multiplier");
    def->graph_settings->null_label  = L("No compensation");
    def->graph_settings->label_min_x = L("Graph min speed");
    def->graph_settings->label_max_x = L("Graph max speed");
    def->graph_settings->label_min_y = L("Minimum flow");
    def->graph_settings->label_max_y = L("Maximum flow");
    def->graph_settings->min_x       = 10;
    def->graph_settings->max_x       = 2000;
    def->graph_settings->step_x      = 1.;
    def->graph_settings->min_y       = 0.1;
    def->graph_settings->max_y       = 2;
    def->graph_settings->step_y      = 0.1;
    def->graph_settings->allowed_types = {GraphData::GraphType::LINEAR, GraphData::GraphType::SQUARE};

    def = this->add("extruder_offset", coPoints);
    def->label = L("Extruder offset");
    def->category = OptionCategory::extruders;
    def->tooltip = L("If your firmware doesn't handle the extruder displacement you need the G-code "
        "to take it into account. This option lets you specify the displacement of each extruder "
        "with respect to the first one. It expects positive coordinates (they will be subtracted "
        "from the XY coordinate).");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPoints{ Vec2d(0,0) });

    def = this->add("extruder_temperature_offset", coFloats);
    def->label = L("Extruder temp offset");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This offset will be added to all extruder temperatures set in the filament settings."
        "\nNote that you should set 'M104 S{first_layer_temperature{initial_extruder} + extruder_temperature_offset{initial_extruder}}'"
        "\ninstead of 'M104 S{first_layer_temperature}' in the start_gcode");
    def->sidetext = L("°C");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0 });

    def = this->add("extruder_fan_offset", coPercents);
    def->label = L("Extruder fan offset");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This offset wil be added to all fan values set in the filament properties. It won't make them go higher than 100% nor lower than 0%.");
    def->sidetext = L("%");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents{ 0 });


    def = this->add("extrusion_axis", coString);
    def->label = L("Extrusion axis");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Use this option to set the axis letter associated with your printer's extruder "
                   "(usually E but some printers use A).");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString("E"));

    def = this->add("extrusion_multiplier", coFloats);
    def->label = L("Extrusion multiplier");
    def->category = OptionCategory::filament;
    def->tooltip = L("This factor changes the amount of flow proportionally. You may need to tweak "
        "this setting to get nice surface finish and correct single wall widths. "
        "Usual values are between 0.9 and 1.1. If you think you need to change this more, "
        "check filament diameter and your firmware E steps.");
    def->mode = comSimpleAE | comPrusa;
    def->min = 0;
    def->max = 2;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 1. });

    def = this->add("print_extrusion_multiplier", coPercent);
    def->label = L("Extrusion multiplier");
    def->category = OptionCategory::filament;
    def->tooltip = L("This factor changes the amount of flow proportionally. You may need to tweak "
        "this setting to get nice surface finish and correct single wall widths. "
        "Usual values are between 90% and 110%. If you think you need to change this more, "
        "check filament diameter and your firmware E steps."
        " This print setting is multiplied against the extrusion_multiplier from the filament tab."
        " Its only purpose is to offer the same functionality but on a per-object basis."); // TODO: replace "against" with "with"?
    def->sidetext = L("%");
    def->mode = comSimpleAE | comSuSi;
    def->min = 0;
    def->max = 200;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("extrusion_width", coFloatOrPercent);
    def->label = L("Default extrusion width");
    def->category = OptionCategory::width;
    def->tooltip = L("This is the DEFAULT extrusion width. It's ONLY used to REPLACE 0-width fields. It's useless when all other width fields have a value."
        "\nSet this to a non-zero value to allow a manual extrusion width. "
        "If left to zero, Slic3r derives extrusion widths from the nozzle diameter "
        "(see the tooltips for perimeter extrusion width, infill extrusion width etc). "
        "If expressed as percentage (for example: 105%), it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comExpert | comPrusa;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("extrusion_spacing", coFloatOrPercent);
    def->label = L("Default extrusion spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("This is the DEFAULT extrusion spacing. It's convert to a width and this width can be used to REPLACE 0-width fields. It's useless when all  width fields have a value."
        "Like Default extrusion width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
                "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

#if 0
    //not used anymore, to remove !! @DEPRECATED (replaces by default_fan_speed)
    def = this->add("fan_always_on", coBools);
    def->label = L("Keep fan always on");
    def->category = OptionCategory::cooling;
    def->tooltip = L("If this is enabled, fan will continuously run at base speed if no other setting overrides that speed."
                " Useful for PLA, harmful for ABS.");
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ true });
#endif

    def = this->add("fan_below_layer_time", coFloats);
    def->label = L("Enable fan if layer print time is below");
    def->category = OptionCategory::cooling;
    def->tooltip = L("If layer print time is estimated below this number of seconds, fan will be enabled "
        "and its speed will be calculated by interpolating the default and maximum speeds."
        "\nSet zero to disable.");
    def->sidetext = L("approximate seconds");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 60 });

    def = this->add("fan_name", coStrings);
    def->label = L("Fan configuration name");
    def->category = OptionCategory::cooling;
    def->tooltip = L("If this field is not empty, the gcode will use this name for this extruder fan, instead of the generic one (only klipper firmware support it right now).");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings{ "" });

    def = this->add("filament_colour", coStrings);
    def->label = L("Color");
    def->full_label = L("Filament color");
    def->category = OptionCategory::filament;
    def->tooltip = L("This is only used in the Slic3r interface as a visual help.");
    def->gui_type = ConfigOptionDef::GUIType::color;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings{ "#29B2B2" });

    def = this->add("filament_custom_variables", coStrings);
    def->label = L("Custom variables");
    def->full_label = L("Custom Filament variables");
    def->category = OptionCategory::filament;
    def->tooltip = L("You can add data accessible to custom-gcode macros."
        "\nEach line can define one variable."
        "\nThe format is 'variable_name=value'. The variable name should only have [a-zA-Z0-9] characters or '_'."
        "\nA value that can be parsed as a int or float will be avaible as a numeric value."
        "\nA value that is enclosed by double-quotes will be available as a string (without the quotes)"
        "\nA value that only takes values as 'true' or 'false' will be a boolean)"
        "\nEvery other value will be parsed as a string as-is."
        "\nThese variables will be available as an array in the custom gcode (one item per extruder), don't forget to use them with the {current_extruder} index to get the current value."
        " If a filament has a typo on the variable that change its type, then the parser will convert everything to strings."
        "\nAdvice: before using a variable, it's safer to use the function 'default_XXX(variable_name, default_value)'"
        " (enclosed in bracket as it's a script) in case it's not set. You can replace XXX by 'int' 'bool' 'double' 'string'.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings{ "" });

    def = this->add("filament_fill_top_flow_ratio", coPercents);
    def->label = L("Top fill");
    def->full_label = L("Top fill flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("You can increase this to over-extrude on the top layer if there is not enough plastic to make a good fill."
                    "\nThis setting multiply the percentage available in the print setting."
                    " You should only add the little percentage difference that this filament has versus your main one.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercents{100});

    def = this->add("filament_first_layer_flow_ratio", coPercents);
    def->label = L("First layer");
    def->full_label = L("First layer flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("You can increase this to over/under-extrude on the first layer if there is not enough / too many plastic because your bed isn't levelled / flat."
                    "\nThis setting multiply the percentage available in the print setting."
                    " You should only add the little percentage difference that this filament has versus your main one.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercents{100});

    def = this->add("filament_notes", coStrings);
    def->label = L("Filament notes");
    def->category = OptionCategory::notes;
    def->tooltip = L("You can put your notes regarding the filament here.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings { "" });

    def = this->add("filament_max_speed", coFloats);
    def->label = L("Max speed");
    def->category = OptionCategory::filament;
    def->tooltip = L("Maximum speed allowed for this filament. Limits the maximum "
        "speed of a print to the minimum of the print speed and the filament speed. "
        "Set zero for no limit.");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0. });

    def = this->add("filament_max_volumetric_speed", coFloats);
    def->label = L("Max volumetric speed");
    def->category = OptionCategory::filament;
    def->tooltip = L("Maximum volumetric speed allowed for this filament. Limits the maximum volumetric "
        "speed of a print to the minimum of print and filament volumetric speed. "
        "Set zero for no limit.");
    def->sidetext = L("mm³/s");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0. });

    def = this->add("filament_max_wipe_tower_speed", coFloats);
    def->label = L("Max speed on the wipe tower");
    def->tooltip = L("This setting is used to set the maximum speed when extruding inside the wipe tower (use M220)."
        " In %, set 0 to disable and use the Filament type instead."
        "\nIf disabled, these filament types will have a defaut value of:"
        "\n - PVA: 80% to 60%"
        "\n - SCAFF: 35%"
        "\n - FLEX: 35%"
        "\n - OTHERS: 100%"
        "\nNote that the wipe tower reset the speed at 100% for the unretract in any case." // TODO: "reset" -> "resets"?
        "\nIf using marlin, M220 B/R is used to save the speed override before the wipe tower print.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 400;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0 });

    def = this->add("filament_loading_speed", coFloats);
    def->label = L("Loading speed");
    def->tooltip = L("Speed used for loading the filament on the wipe tower. ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 28. });

    //skinnydip section starts
    def = this->add("filament_enable_toolchange_temp", coBools);
    def->label = L("Toolchange temperature enabled");
    def->tooltip = L("Determines whether toolchange temperatures will be applied");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_use_fast_skinnydip", coBools);
    def->label = L("Fast mode");
    def->tooltip = L("Experimental: drops nozzle temperature during cooling moves instead of prior to extraction to reduce wait time.");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_enable_toolchange_part_fan", coBools);
    def->label = L("Use part fan to cool hotend");
    def->tooltip = L("Experimental setting.  May enable the hotend to cool down faster during toolchanges");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_toolchange_part_fan_speed", coInts);
    def->label = L("Toolchange part fan speed");
    def->tooltip = L("Experimental setting.  Fan speeds that are too high can clash with the hotend's PID routine.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 50 });

    def = this->add("filament_use_skinnydip", coBools);
    def->label = L("Enable Skinnydip string reduction");
    def->tooltip = L("Skinnydip performs a secondary dip into the meltzone to burn off fine strings of filament");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_melt_zone_pause", coInts);
    def->label = L("Pause in melt zone");
    def->tooltip = L("Stay in melt zone for this amount of time before extracting the filament.  Not usually necessary.");
    def->sidetext = L("milliseconds");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 0 });

    def = this->add("filament_cooling_zone_pause", coInts);
    def->label = L("Pause before extraction ");
    def->tooltip = L("Can be useful to avoid bondtech gears deforming hot tips, but not ordinarily needed");
    def->sidetext = L("milliseconds");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 0 });

    def = this->add("filament_dip_insertion_speed", coFloats);
    def->label = L("Speed to move into melt zone");
    def->tooltip = L("usually not necessary to change this");
    def->sidetext = L("mm/sec");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 33. });

    def = this->add("filament_dip_extraction_speed", coFloats);
    def->label = L("Speed to extract from melt zone");
    def->tooltip = L("usually not necessary to change this");
    def->sidetext = L("mm/sec");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 70. });

    def = this->add("filament_toolchange_temp", coInts);
    def->label = L("Toolchange temperature");
    def->tooltip = L("To further reduce stringing, it can be helpful to set a lower temperature just prior to extracting filament from the hotend.");
    def->sidetext = L("°C");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 200 });

    def = this->add("filament_skinnydip_distance", coFloats);
    def->label = L("Insertion distance");
    def->tooltip = L("For stock extruders, usually 40-42mm.  For bondtech extruder upgrade, usually 30-32mm.  Start with a low value and gradually increase it until strings are gone.  If there are blobs on your wipe tower, your value is too high.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 31. });
    //skinnydip section ends

    def = this->add("filament_loading_speed_start", coFloats);
    def->label = L("Loading speed at the start");
    def->tooltip = L("Speed used at the very beginning of loading phase. ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 3. });

    def = this->add("filament_unloading_speed", coFloats);
    def->label = L("Unloading speed");
    def->tooltip = L("Speed used for unloading the filament on the wipe tower (does not affect "
                      " initial part of unloading just after ramming). ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 90. });

    def = this->add("filament_unloading_speed_start", coFloats);
    def->label = L("Unloading speed at the start");
    def->tooltip = L("Speed used for unloading the tip of the filament immediately after ramming. ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 100. });

    def = this->add("filament_toolchange_delay", coFloats);
    def->label = L("Delay after unloading");
    def->tooltip = L("Time to wait after the filament is unloaded. "
                   "May help to get reliable toolchanges with flexible materials "
                   "that may need more time to shrink to original dimensions. ");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("filament_cooling_moves", coInts);
    def->label = L("Number of cooling moves");
    def->tooltip = L("Filament is cooled by being moved back and forth in the "
                   "cooling tubes. Specify desired number of these moves.");
    def->max = 0;
    def->max = 20;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 4 });

    def = this->add("filament_cooling_initial_speed", coFloats);
    def->label = L("Speed of the first cooling move");
    def->tooltip = L("Cooling moves are gradually accelerated, starting at this speed. ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 2.2 });

    def = this->add("filament_minimal_purge_on_wipe_tower", coFloats);
    def->label = L("Minimal purge on wipe tower");
    def->tooltip = L("After a tool change, the exact position of the newly loaded filament inside "
                     "the nozzle may not be known, and the filament pressure is likely not yet stable. "
                     "Before purging the print head into an infill or a sacrificial object, Slic3r will always prime "
                     "this amount of material into the wipe tower to produce successive infill or sacrificial object extrusions reliably.");
    def->sidetext = L("mm³");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 15. });

    def = this->add("filament_cooling_final_speed", coFloats);
    def->label = L("Speed of the last cooling move");
    def->tooltip = L("Cooling moves are gradually accelerated towards this speed. ");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 3.4 });

    def = this->add("filament_load_time", coFloats);
    def->label = L("Filament load time");
    def->tooltip = L("Time for the printer firmware (or the Multi Material Unit 2.0) to load a new filament during a tool change (when executing the T code). This time is added to the total print time by the G-code time estimator.");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0.0 });

    def = this->add("filament_ramming_parameters", coStrings);
    def->label = L("Ramming parameters");
    def->tooltip = L("This string is edited by RammingDialog and contains ramming specific parameters.");
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings { "120 100 6.6 6.8 7.2 7.6 7.9 8.2 8.7 9.4 9.9 10.0|"
       " 0.05 6.6 0.45 6.8 0.95 7.8 1.45 8.3 1.95 9.7 2.45 10 2.95 7.6 3.45 7.6 3.95 7.6 4.45 7.6 4.95 7.6" });

    def = this->add("filament_unload_time", coFloats);
    def->label = L("Filament unload time");
    def->tooltip = L("Time for the printer firmware (or the Multi Material Unit 2.0) to unload a filament during a tool change (when executing the T code). This time is added to the total print time by the G-code time estimator.");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0.0 });

    def = this->add("filament_multitool_ramming", coBools);
    def->label = L("Enable ramming for multitool setups");
    def->tooltip = L("Perform ramming when using multitool printer (i.e. when the 'Single Extruder Multimaterial' in Printer Settings is unchecked). "
                     "When checked, a small amount of filament is rapidly extruded on the wipe tower just before the toolchange. "
                     "This option is only used when the wipe tower is enabled.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_multitool_ramming_volume", coFloats);
    def->label = L("Multitool ramming volume");
    def->tooltip = L("The volume to be rammed before the toolchange.");
    def->sidetext = L("mm³");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloats { 10. });

    def = this->add("filament_multitool_ramming_flow", coFloats);
    def->label = L("Multitool ramming flow");
    def->tooltip = L("Flow used for ramming the filament before the toolchange.");
    def->sidetext = L("mm³/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloats { 10. });

    def = this->add("filament_diameter", coFloats);
    def->label = L("Diameter");
    def->tooltip = L("Enter your filament diameter here. Good precision is required, so use a caliper "
                   "and do multiple measurements along the filament, then compute the average.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 1.75 });

    def = this->add("filament_shrink", coPercents);
    def->label = L("Shrinkage");
    def->tooltip = L("Enter the shrinkage percentage that the filament will get after cooling (94% if you measure 94mm instead of 100mm)."
        " The part will be scaled in xy to compensate."
        " Only the filament used for the perimeter is taken into account."
        "\nBe sure to allow enough space between objects, as this compensation is done after the checks.");
    def->sidetext = L("%");
    def->ratio_over = "";
    def->min = 10;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents{ 100 });

    def = this->add("filament_max_overlap", coPercents);
    def->label = L("Max line overlap");
    def->tooltip = L("This setting will ensure that all 'overlap' are not higher than this value."
        " This is useful for filaments that are too viscous, as the line can't flow under the previous one."
        "\nNote: top solid infill lines are excluded, to prevent visual defects.");
    def->sidetext = L("%");
    def->ratio_over = "";
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents{ 100 });

    def = this->add("filament_density", coFloats);
    def->label = L("Density");
    def->category = OptionCategory::filament;
    def->tooltip = L("Enter your filament density here. This is only for statistical information. "
                   "A decent way is to weigh a known length of filament and compute the ratio "
                   "of the length to volume. Better is to calculate the volume directly through displacement.");
    def->sidetext = L("g/cm³");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0. });

    def = this->add("filament_type", coStrings);
    def->label = L("Filament type");
    def->category = OptionCategory::filament;
    def->tooltip = L("The filament material type for use in custom G-codes.");
    def->gui_flags = "show_value";
    def->set_enum_values(ConfigOptionDef::GUIType::select_open, {
        "PLA", 
        "PET",
        "ABS",
        "ASA",
        "FLEX", 
        "HIPS",
        "EDGE",
        "NGEN",
        "PA",
        "NYLON",
        "PVA",
        "PC",
        "PP",
        "PEI",
        "PEEK",
        "PEKK",
        "POM",
        "PSU",
        "PVDF",
        "SCAFF"
    });
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings { "PLA" });

    def = this->add("filament_soluble", coBools);
    def->label = L("Soluble material");
    def->category = OptionCategory::filament;
    def->tooltip = L("Soluble material is most likely used for a soluble support.");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("filament_cost", coFloats);
    def->label = L("Cost");
    def->full_label = L("Filament cost");
    def->category = OptionCategory::filament;
    def->tooltip = L("Enter your filament cost per kg here. This is only for statistical information.");
    def->sidetext = L("money/kg");
    def->min = 0;
    def->is_vector_extruder = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("filament_spool_weight", coFloats);
    def->label = L("Spool weight");
    def->category = OptionCategory::filament;
    def->tooltip = L("Enter weight of the empty filament spool. "
                     "One may weigh a partially consumed filament spool before printing and one may compare the measured weight "
                     "with the calculated weight of the filament with the spool to find out whether the amount "
                     "of filament on the spool is sufficient to finish the print.");
    def->sidetext = L("g");
    def->min = 0;
    def->is_vector_extruder = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("filament_settings_id", coStrings);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionStrings { "" });
    def->cli = ConfigOptionDef::nocli;

    def = this->add("filament_settings_modified", coBools);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionBools({false}));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("filament_vendor", coString);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString(L("(Unknown)")));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("fill_aligned_z", coBool);
    def->label = L("Align sparse infill in z");
    def->category = OptionCategory::infill;
    def->tooltip = L("The voids of the pattern of sparse infill grows with the extrusion width, to keep the same percentage of fill."
                    " With this setting set to true, the algorithms will now only use the highest sparse infill width avaialble to create the pattern."
                    " This way, the pattern can still be aligned even if the width is changing (from first layer width, from a modifier, from aanother extruder with different diameter)."
                    "\nExperimental: works only for infill that won't depends on the fill area. So for infill where it's not useful (Hilbert, Archimedean, Octagram, Scattered, Lightning), this setting is disabled."
                    "\n This setting is useful for rectilinear, monotonic, grid, trianlge, star, cubic, gyroid, honeycomb, 3D honeycomb, adaptative cubic, support cubic patterns.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("fill_angle", coFloat);
    def->label = L("Fill");
    def->full_label = L("Fill angle");
    def->category = OptionCategory::infill;
    def->tooltip = L("Default base angle for infill orientation. Cross-hatching will be applied to this. "
                   "Bridges will be infilled using the best direction Slic3r can detect, so this setting "
                   "does not affect them.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 360;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(45));

    def             = this->add("fill_angle_cross", coBool);
    def->label      = L("Alternate Fill Angle");
    def->category   = OptionCategory::infill;
    def->tooltip    = L("It's better for some infill like rectilinear to rotate 90° each layer. If this setting is deactivated, they won't do that anymore.");
    def->mode       = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def             = this->add("fill_angle_follow_model", coBool);
    def->label      = L("Rotate with object");
    def->category   = OptionCategory::infill;
    def->tooltip    = L("If your object has a z-rotation, then the infill will also be rotated by this value.");
    def->mode       = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("fill_angle_increment", coFloat);
    def->label = L("Fill");
    def->full_label = L("Fill angle increment");
    def->category = OptionCategory::infill;
    def->tooltip = L("Add this angle each layer to the base angle for infill. "
                    "May be useful for art, or to be sure to hit every object's feature even with very low infill. "
                    "Still experimental, tell me what makes it useful, or the problems that arise using it.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 360;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def             = this->add("fill_angle_template", coFloats);   
    def->label      = L("Fill angle template");
    def->full_label = L("Fill angle template");
    def->category   = OptionCategory::infill;
    def->tooltip    = L("This define the succetion of infill angle. When defined, it replaces the fill_angle"
        ", and there won't be any extra 90° for each layer added, but the fill_angle_increment will still be used."
        " The first layer start with the first angle. If a new pattern is used in a modifier"
        ", it will choose the layer angle from the pattern as if it has started from the first layer."
        "Empty this settings to disable and recover the old behavior.");
    def->sidetext   = L("°");
    def->min        = 0;
    def->max        = 360;
    def->full_width = true;
    def->mode       = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloats(0.));

    def = this->add("fill_density", coPercent);
    def->gui_flags = "show_value";
    def->label = L("Fill density");
    def->category = OptionCategory::infill;
    def->tooltip = L("Density of internal infill, expressed in the range 0% - 100%."
        "\nSet 0 to remove any sparse infill."
        "\nNote that using a value of 100% won't change the type of infill from sparse to solid."
        " If you want only solid infill, you can set the 'Solid infill every X layers' (solid_infill_every_layers) to 1 instead.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->set_enum_values(ConfigOptionDef::GUIType::f_enum_open, {
        { "0", "0%" },
        { "4", "4%" },
        { "5.5", "5.5%" },
        { "7.5", "7.5%" },
        { "10", "10%" },
        { "13", "13%" },
        { "18", "18%" },
        { "23", "23%" },
        { "31", "31%" },
        { "42", "42%" },
        { "55", "55%" },
        { "75", "75%" },
        //{ "100", "100%" } // can still be entered, but not showing it may make people increase solid layer count instead (which is the proper way).
    });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionPercent(18));

    def = this->add("fill_pattern", coEnum);
    def->label = L("Pattern");
    def->full_label = L("Sparse fill pattern");
    def->category = OptionCategory::infill;
    def->tooltip = L("Fill pattern for general low-density infill."
        "\nIf you want an 'aligned' pattern, set 90° to the fill angle increment setting.");
    def->set_enum<InfillPattern>({
        { "rectilinear",        L("Rectilinear") },
        { "alignedrectilinear", L("Aligned Rectilinear") },
        { "monotonic",          L("Monotonic") },
        { "grid",               L("Grid") }, 
        { "triangles",          L("Triangles")},
        { "stars",              L("Stars")},
        { "cubic",              L("Cubic")},
        { "line",               L("Line")},
        { "concentric",         L("Concentric")},
        { "honeycomb",          L("Honeycomb")},
        { "3dhoneycomb",        L("3D Honeycomb")},
        { "gyroid",             L("Gyroid")},
        { "hilbertcurve",       L("Hilbert Curve")},
        { "archimedeanchords",  L("Archimedean Chords")},
        { "octagramspiral",     L("Octagram Spiral")},
        {"scatteredrectilinear",L("Scattered Rectilinear")},
        { "adaptivecubic",      L("Adaptive Cubic")},
        { "supportcubic",       L("Support Cubic")},
        { "lightning",          L("Lightning")}
    });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value( new ConfigOptionEnum<InfillPattern>(ipStars));

    def = this->add("fill_top_flow_ratio", coPercent);
    def->label = L("Top fill");
    def->full_label = L("Top fill flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("You can increase this to over-extrude on the top layer if there is not enough plastic to make a good fill.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("first_layer_flow_ratio", coPercent);
    def->label = L("First layer");
    def->full_label = L("First layer flow ratio");
    def->sidetext = L("%");
    def->category = OptionCategory::width;
    def->tooltip = L("You can increase this to over-extrude on the first layer if there is not enough plastic because your bed isn't levelled."
                    "\nNote: DON'T USE THIS if your only problem is bed leveling, LEVEL YOUR BED!"
                    " Use this setting only as last resort after all calibrations failed.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("first_layer_size_compensation", coFloat);
    def->label = L("First layer");
    def->full_label = L("XY First layer compensation");
    def->category = OptionCategory::slicing;
    def->tooltip = L("The first layer will be grown / shrunk in the XY plane by the configured value "
        "to compensate for the 1st layer squish aka an Elephant Foot effect. (should be negative = inwards = remove area)");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi | comPrusa; // just a rename & inverted of prusa 's elefant_foot
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("first_layer_size_compensation_layers", coInt);
    def->label = L("height in layers");
    def->full_label = L("XY First layer compensation height in layers");
    def->category = OptionCategory::slicing;
    def->tooltip = L("The number of layers on which the elephant foot compensation will be active. "
        "The first layer will be shrunk by the elephant foot compensation value, then "
        "the next layers will be gradually shrunk less, up to the layer indicated by this value.");
    def->sidetext = L("layers");
    def->min = 1;
    def->max = 30;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("fill_smooth_width", coFloatOrPercent);
    def->label = L("Width");
    def->full_label = L("Ironing width");
    def->category = OptionCategory::infill;
    def->tooltip = L("This is the width of the ironing pass, in a % of the top infill extrusion width, should not be more than 50%"
        " (two times more lines, 50% overlap). It's not necessary to go below 25% (four times more lines, 75% overlap). \nIf you have problems with your ironing process,"
        " don't forget to look at the flow->above bridge flow, as this setting should be set to min 110% to let you have enough plastic in the top layer."
        " A value too low will make your extruder eat the filament.");
    def->ratio_over = "top_infill_extrusion_width";
    def->min = 0;
    def->max_literal = { 1, true };
    def->mode = comExpert | comSuSi;
    def->sidetext = L("mm/%");
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("fill_smooth_distribution", coPercent);
    def->label = L("Distribution");
    def->full_label = L("Ironing flow distribution");
    def->category = OptionCategory::infill;
    def->tooltip = L("This is the percentage of the flow that is used for the second ironing pass. Typical 10-20%. "
        "Should not be higher than 20%, unless you have your top extrusion width greatly superior to your nozzle width. "
        "A value too low and your extruder will eat the filament. A value too high and the first pass won't print well.");
    //def->min = 0;
    //def->max = 0.9;
    def->mode = comExpert | comSuSi;
    def->sidetext = L("%");
    def->set_default_value(new ConfigOptionPercent(10));

    def = this->add("small_area_infill_flow_compensation", coBool);
    def->label = L("Enable small area flow compensation");
    def->category = OptionCategory::infill;
    def->tooltip = L("Enable flow compensation for small infill areas."
                    "\nFirst layer is always disabled, to not compromise adhesion.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("small_area_infill_flow_compensation_model", coGraph);
    def->label = L("Flow Compensation Model");
    def->category = OptionCategory::infill;
    def->tooltip = L("Flow Compensation Model, used to adjust the flow for small solid infill "
                     "lines. The model is a graph of flow correction factors (between 0 and 1) per extrusion length (in mm)."
                     "\nThe first point length has to be 0mm. the last point need to have a flow correction of 1.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionGraph(GraphData(0,10, GraphData::GraphType::SPLINE,
        {{0,0},{0.2,0.44},{0.4,0.61},{0.6,0.7},{0.8,0.76},{1.5,0.86},{2,0.89},{3,0.92},{5,0.95},{10,1}}
    )));
    def->graph_settings = std::make_shared<GraphSettings>();
    def->graph_settings->title       = L("Flow Compensation Model");
    def->graph_settings->description = def->tooltip;
    def->graph_settings->x_label     = L("Length of an extrusion (mm)");
    def->graph_settings->y_label     = L("Flow correction (ratio between 0 and 1)");
    def->graph_settings->null_label  = L("No values");
    def->graph_settings->label_min_x = "";
    def->graph_settings->label_max_x = L("Maximum length");
    def->graph_settings->label_min_y = L("Minimum ratio");
    def->graph_settings->label_max_y = L("Maximum ratio");
    def->graph_settings->min_x       = 0;
    def->graph_settings->max_x       = 100;
    def->graph_settings->step_x      = 0.1;
    def->graph_settings->min_y       = 0;
    def->graph_settings->max_y       = 1;
    def->graph_settings->step_y      = 0.01;
    def->graph_settings->allowed_types = {GraphData::GraphType::LINEAR, GraphData::GraphType::SPLINE, GraphData::GraphType::SQUARE};

    def = this->add("first_layer_acceleration", coFloatOrPercent);
    def->label = L("Max");
    def->full_label = L("First layer acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the maximum acceleration your printer will use for first layer."
                "\nIf set to %, all accelerations will be reduced by that ratio."
                "\nSet zero to disable acceleration control for first layer.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "depends";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("first_layer_acceleration_over_raft", coFloatOrPercent);
    def->label = L("First object layer over raft interface");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for first layer of object above raft interface."
        "\nIf set to %, all accelerations will be reduced by that ratio."
        "\nSet zero to disable acceleration control for first layer of object above raft interface.");
    def->sidetext = L("mm/s²");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("first_layer_bed_temperature", coInts);
    def->label = L("First layer");
    def->full_label = L("First layer bed temperature");
    def->category = OptionCategory::filament;
    def->tooltip = L("Heated build plate temperature for the first layer. Set zero to disable "
                   "bed temperature control commands in the output.");
    def->sidetext = L("°C");
    def->max = 0;
    def->max = 300;
    def->is_vector_extruder = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInts { 0 });

    def = this->add("first_layer_extrusion_width", coFloatOrPercent);
    def->label = L("First layer");
    def->full_label = L("First layer width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for first layer. "
        "You can use this to force fatter extrudates for better adhesion. If expressed "
        "as percentage (for example 140%) it will be computed over the nozzle diameter "
        "of the nozzle used for the type of extrusion. "
        "If set to zero, it will use the default extrusion width."
        "If disabled, nothing is changed compared to a normal layer."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->can_be_disabled = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloatOrPercent(140, true), false));

    def = this->add("first_layer_extrusion_spacing", coFloatOrPercent);
    def->label = L("First layer");
    def->full_label = L("First layer spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like First layer width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));
    

    def = this->add("first_layer_infill_extrusion_width", coFloatOrPercent);
    def->label = L("First layer");
    def->full_label = L("First layer infill width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for first layer infill (sparse and solid). "
        "You can use this to force fatter extrudates for better adhesion. If expressed "
        "as percentage (for example 140%) it will be computed over the nozzle diameter "
        "of the nozzle used for the type of extrusion. "
        "If set to zero, it will use the default extrusion width."
        "If disabled, the first layer width is also used for first layer infills (if enabled)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->can_be_disabled = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloatOrPercent(140, true)));

    def = this->add("first_layer_infill_extrusion_spacing", coFloatOrPercent);
    def->label = L("First layer");
    def->full_label = L("First layer infill spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like First layer infill width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("first_layer_height", coFloatOrPercent);
    def->label = L("First layer height");
    def->category = OptionCategory::slicing;
    def->tooltip = L("When printing with very low layer heights, you might still want to print a thicker "
                   "bottom layer to improve adhesion and tolerance for non perfect build plates. "
                   "This can be expressed as an absolute value or as a percentage (for example: 75%) "
                   "over the lowest nozzle diameter used in by the object.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 20, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(75, true));

    def = this->add("first_layer_speed", coFloatOrPercent);
    def->label = L("Max");
    def->full_label = L("Default first layer speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("If expressed as absolute value in mm/s, this speed will be applied as a maximum to all the print moves (but infill) of the first layer."
        "\nIf expressed as a percentage it will scale the current speed."
        "\nSet it at 100% to remove any first layer speed modification (but for infill).");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "depends";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(30, false));

    def = this->add("first_layer_speed_over_raft", coFloatOrPercent);
    def->label = L("Max speed of object first layer over raft interface");
    def->category = OptionCategory::speed;
    def->tooltip = L("If expressed as absolute value in mm/s, this speed will be usedf as a max over all the print moves "
        "of the first object layer above raft interface, regardless of their type."
        "\nIf expressed as a percentage it will scale the current speed (max 100%)."
        "\nSet it at 100% to remove this speed modification.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "depends";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(30, false));
    
    def = this->add("first_layer_infill_speed", coFloatOrPercent);
    def->label = L("Max infill");
    def->full_label = L("Infill max first layer speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("If expressed as absolute value in mm/s, this speed will be applied as a maximum for all infill print moves of the first layer."
                   "\nIf expressed as a percentage it will scale the current infill speed."
                   "\nSet it at 100% to remove any infill first layer speed modification."
                   "\nSet zero to disable (using first_layer_speed instead).");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "depends";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("first_layer_min_speed", coFloat);
    def->label = L("Min");
    def->full_label = L("Min first layer speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Minimum speed when printing the first layer."
        "\nSet zero to disable.");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));
    
    def = this->add("first_layer_temperature", coInts);
    def->label = L("First layer");
    def->full_label = L("First layer nozzle temperature");
    def->category = OptionCategory::filament;
    def->tooltip = L("Extruder nozzle temperature for first layer. If you want to control temperature manually "
                   "during print, set zero to disable temperature control commands in the output file.");
    def->sidetext = L("°C");
    def->min = 0;
    def->max = max_temp;
    def->is_vector_extruder = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInts { 200 });

    def = this->add("full_fan_speed_layer", coInts);
    def->label = L("Full fan speed at layer");
    def->category = OptionCategory::filament;
    def->tooltip = L("Fan speed will be ramped up linearly from zero at layer \"disable_fan_first_layers\" "
                   "to maximum at layer \"full_fan_speed_layer\". "
                   "\"full_fan_speed_layer\" will be ignored if equal or lower than \"disable_fan_first_layers\", in which case "
                   "the fan will be running at maximum allowed speed at layer \"disable_fan_first_layers\" + 1."
                   "\nset 0 to disable");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 4 });

    def = this->add("fuzzy_skin", coEnum);
    def->label = L("Fuzzy Skin");
    def->category = OptionCategory::fuzzy_skin;
    def->tooltip = L("Fuzzy skin type.");
    def->set_enum<FuzzySkinType>({
        { "none",       L("None") },
        { "external",   L("Outside walls") },
        { "shell",      L("External walls") },
        { "all",        L("All walls") }
    });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<FuzzySkinType>(FuzzySkinType::None));

    def = this->add("fuzzy_skin_thickness", coFloatOrPercent);
    def->label = L("Fuzzy skin thickness");
    def->category = OptionCategory::fuzzy_skin;
    def->tooltip = L("The maximum distance that each skin point can be offset (both ways), "
        "measured perpendicular to the perimeter wall."
        "\nCan be a % of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(150, true));

    def = this->add("fuzzy_skin_point_dist", coFloatOrPercent);
    def->label = L("Fuzzy skin point distance");
    def->category = OptionCategory::fuzzy_skin;
    def->tooltip = L("Perimeters will be split into multiple segments by inserting Fuzzy skin points. "
        "Lowering the Fuzzy skin point distance will increase the number of randomly offset points on the perimeter wall."
        "\nCan be a % of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(200, true));

    def = this->add("gap_fill_enabled", coBool);
    def->label = L("Gap fill");
    def->full_label = L("Enable Gap fill");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Enable gap fill algorithm. It will extrude small lines between perimeters "
        "when there is not enough space for another perimeter or an infill.");
    def->mode = comAdvancedE | comPrusa;
    def->aliases = { "gap_fill" }; //superslicer 2.3 or older
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("gap_fill_acceleration", coFloatOrPercent);
    def->label = L("Gap fill");
    def->full_label = L("Gap fill acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for gap fills. "
                "\nThis can be expressed as a percentage over the perimeter acceleration."
                "\nSet zero to use perimeter acceleration for gap fills.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "perimeter_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("gap_fill_extension", coFloatOrPercent);
    def->label = L("Extension");
    def->full_label = L("Gap fill: extra extension");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Increase the length of all gapfills by this amount (may overextrude a little bit)\nCan be a % of the perimeter width");
    def->ratio_over = "perimeter_width";
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max_literal = { 50, true };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent{ 0, false });

    def = this->add("gap_fill_fan_speed", coInts);
    def->label = L("Gap fill fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all gap fill Perimeter moves"
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("gap_fill_flow_match_perimeter", coPercent);
    def->label = L("Cap with perimeter flow");
    def->full_label = L("Gapfill: cap speed with perimeter flow");
    def->category = OptionCategory::speed;
    def->tooltip = L("A percentage of the perimeter flow (mm3/s) is used as a limit for the gap fill flow, and so the gapfill may reduce its speed when the gap fill extrusions became too thick."
                " This allow you to use a high gapfill speed, to print the thin gapfill quickly and reduce the difference in flow rate for the gapfill."
                "\nSet zero to deactivate.");
    def->sidetext = L("%");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(0));

    def = this->add("gap_fill_last", coBool);
    def->label = L("after last perimeter");
    def->full_label = L("Gapfill: after last perimeter");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("All gaps, between the last perimeter and the infill, which are thinner than a perimeter will be filled by gapfill.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("gap_fill_max_width", coFloatOrPercent);
    def->label = L("Max width");
    def->full_label = L("Gapfill: Max width");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This setting represents the maximum width of a gapfill. Points wider than this threshold won't be created.\nCan be a % of the perimeter width\n0 to auto");
    def->ratio_over = "perimeter_width";
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent{ 0, false });

    def = this->add("gap_fill_min_area", coFloatOrPercent);
    def->label = L("Min surface");
    def->full_label = L("Gapfill: Min surface");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This setting represents the minimum mm² for a gapfill extrusion to be created.\nCan be a % of (perimeter width)²");
    def->ratio_over = "perimeter_width_square";
    def->sidetext = L("mm² or %");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent{ 100, true });

    def = this->add("gap_fill_min_length", coFloatOrPercent);
    def->label = L("Min length");
    def->full_label = L("Gapfill: Min length");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This setting represents the minimum mm for a gapfill extrusion to be extruded.\nCan be a % of the perimeter width\n0 to auto");
    def->ratio_over = "perimeter_width";
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent{ 0, false });

    def = this->add("gap_fill_min_width", coFloatOrPercent);
    def->label = L("Min width");
    def->full_label = L("Gapfill: Min width");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This setting represents the minimum width of a gapfill. Points thinner than this threshold won't be created.\nCan be a % of the perimeter width\n0 to auto");
    def->ratio_over = "perimeter_width";
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent{ 0, false });

    def = this->add("gap_fill_overlap", coPercent);
    def->label = L("Gap fill overlap");
    def->full_label = L("Gap fill overlap");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allows you to reduce the overlap between the perimeters and the gap fill."
        " 100% means that no gaps are left, and 0% means that the gap fill won't touch the perimeters."
        "\nMay be useful if you can see the gapfill on the exterrnal surface, to reduce that artifact.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("gap_fill_speed", coFloatOrPercent);
    def->label = L("Gap fill");
    def->full_label = L("Gap fill speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for filling small gaps using short zigzag moves. Keep this reasonably low "
        "to avoid too much shaking and resonance issues."
        "\nGap fill extrusions are ignored from the automatic volumetric speed computation, unless you set it to 0."
        "\nThis can be expressed as a percentage (for example: 80%) over the Internal Perimeter speed.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "perimeter_speed";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50,true));

    def = this->add("gcode_ascii", coBool);
    def->label = L("Only ascii characters in gcode");
    def->category = OptionCategory::firmware;
    def->tooltip = L("When printing the gcode file, replace any non-ascii character by a '_'."
        " Can be useful if the firmware or a software in a workflow doesn't support uft-8.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("gcode_comments", coBool);
    def->label = L("Verbose G-code");
    def->category = OptionCategory::output;
    def->tooltip = L("Enable this to get a commented G-code file, with each line explained by descriptive text. "
        "If you print from an SD card, the additional weight of the file could make your firmware "
        "slow down.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(0));

    def = this->add("gcode_filename_illegal_char", coString);
    def->label = L("Illegal characters");
    def->full_label = L("Illegal characters for filename");
    def->category = OptionCategory::output;
    def->tooltip = L("All characters that are written here will be replaced by '_' when writing the gcode file name."
        "\nIf the first charater is '[' or '(', then this field will be considered as a regexp (enter '[^a-zA-Z0-9]' to only use ascii char).");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionString("[<>:\"/\\\\|?*]"));

    def = this->add("gcode_flavor", coEnum);
    def->label = L("G-code flavor");
    def->category = OptionCategory::general;
    def->tooltip = L("Some G/M-code commands, including temperature control and others, are not universal. "
                   "Set this option to your printer's firmware to get a compatible output. "
                   "The \"No extrusion\" flavor prevents Slic3r from exporting any extrusion value at all.");
    def->set_enum<GCodeFlavor>({
        { "sprinter",       "RepRap/Sprinter" },
        { "reprapfirmware", "RepRapFirmware" },
        { "repetier",       "Repetier" },
        { "teacup",         "Teacup" },
        { "makerware",      "MakerWare (MakerBot)" },
        { "marlin",         "Marlin (legacy)" },
        { "marlin2",        "Marlin 2" },
        { "klipper",        "Klipper" },
        { "sailfish",       "Sailfish (MakerBot)" },
        { "mach3",          "Mach3/LinuxCNC" },
        { "machinekit",     "Machinekit" },
        { "smoothie",       "Smoothie" },
        { "no-extrusion",   L("No extrusion") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<GCodeFlavor>(gcfMarlinLegacy));

    def = this->add("gcode_label_objects", coEnum);
    def->label = L("Label objects");
    def->tooltip = L("Selects whether labels should be exported at object boundaries and in what format.\n"
                     "OctoPrint = comments to be consumed by OctoPrint CancelObject plugin.\n"
                     "Firmware = firmware specific G-code (it will be chosen based on firmware flavor and it can end up to be empty).\n\n"
                     "This settings is NOT compatible with Single Extruder Multi Material setup and Wipe into Object / Wipe into Infill.");

    def->set_enum<LabelObjectsStyle>({
        { "disabled",   L("Disabled") },
        { "octoprint",  L("OctoPrint comments") },
        { "firmware",   L("Firmware-specific") },
        { "both",       L("Print both") }
        });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<LabelObjectsStyle>(LabelObjectsStyle::Octoprint));

    def = this->add("gcode_precision_xyz", coInt);
    def->label = L("xyz decimals");
    def->category = OptionCategory::output;
    def->tooltip = L("Choose how many digits after the dot for xyz coordinates.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionInt(3));

    def = this->add("gcode_precision_e", coInt);
    def->label = L("Extruder decimals");
    def->category = OptionCategory::output;
    def->tooltip = L("Choose how many digits after the dot for extruder moves.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionInt(5));

    def = this->add("gcode_substitutions", coStrings);
    def->label = L("G-code substitutions");
    def->tooltip = L("Find / replace patterns in G-code lines and substitute them.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionStrings());

    def = this->add("high_current_on_filament_swap", coBool);
    def->label = L("High extruder current on filament swap");
    def->category = OptionCategory::general;
    def->tooltip = L("It may be beneficial to increase the extruder motor current during the filament exchange"
                   " sequence to allow for rapid ramming feed rates and to overcome resistance when loading"
                   " a filament with an ugly shaped tip.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(0));

    def = this->add("infill_acceleration", coFloatOrPercent);
    def->label = L("Sparse");
    def->full_label = L("Infill acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for Sparse infill."
                "\nCan be a % of the solid infill acceleration"
                "\nSet zero to use solid infill acceleration for infill.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "solid_infill_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("travel_acceleration", coFloat);
    def->label = L("Travel");
    def->tooltip = L("This is the acceleration your printer will use for travel moves. Set zero to disable "
                     "acceleration control for travel.");
    def->sidetext = L("mm/s²");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("infill_every_layers", coInt);
    def->label = L("Combine infill every");
    def->category = OptionCategory::infill;
    def->tooltip = L("This feature allows you to combine infill and speed up your print by extruding thicker "
                   "infill layers while preserving thin perimeters, thus accuracy.");
    def->sidetext = L("layers");
    def->full_label = L("Combine infill every n layers");
    def->min = 1;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    auto def_infill_anchor_min = def = this->add("infill_anchor", coFloatOrPercent);
    def->label = L("Length of the infill anchor");
    def->category = OptionCategory::infill;
    def->tooltip = L("Connect an infill line to an internal perimeter with a short segment of an additional perimeter. "
                     "If expressed as percentage (example: 15%) it is calculated over infill extrusion width. Slic3r tries to connect two close infill lines to a short perimeter segment. If no such perimeter segment "
                     "shorter than infill_anchor_max is found, the infill line is connected to a perimeter segment at just one side "
                     "and the length of the perimeter segment taken is limited to this parameter, but no longer than anchor_length_max. "
                     "\nSet this parameter to zero to disable anchoring perimeters connected to a single infill line.");
    def->sidetext = L("mm or %");
    def->ratio_over = "infill_extrusion_width";
    def->max_literal = { 1000, false };
    def->set_enum_values(ConfigOptionDef::GUIType::f_enum_open, {
        { "0",      L("0 (no open anchors)") },
        { "1",      L("1 mm") },
        { "2",      L("2 mm") },
        { "5",      L("5 mm") },
        { "10",     L("10 mm") },
        { "1000",   L("1000 (unlimited)") }
    });
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(600, true));

    def = this->add("infill_anchor_max", coFloatOrPercent);
    def->label = L("Maximum length of the infill anchor");
    def->category    = def_infill_anchor_min->category;
    def->tooltip = L("Connect an infill line to an internal perimeter with a short segment of an additional perimeter. "
                     "If expressed as percentage (example: 15%) it is calculated over infill extrusion width. Slic3r tries to connect two close infill lines to a short perimeter segment. If no such perimeter segment "
                     "shorter than this parameter is found, the infill line is connected to a perimeter segment at just one side "
                     "and the length of the perimeter segment taken is limited to infill_anchor, but no longer than this parameter. "
                     "\nIf set to 0, the old algorithm for infill connection will be used, it should create the same result as with 1000 & 0.");
    def->sidetext    = def_infill_anchor_min->sidetext;
    def->ratio_over  = def_infill_anchor_min->ratio_over;
    def->set_enum_values(ConfigOptionDef::GUIType::f_enum_open, {
        { "0",      L("0 (Simple connect)") },
        { "1",      L("1 mm") },
        { "2",      L("2 mm") },
        { "5",      L("5 mm") },
        { "10",     L("10 mm") },
        { "1000",   L("1000 (unlimited)") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("infill_connection", coEnum);
    def->label = L("Connection of sparse infill lines");
    def->category = OptionCategory::infill;
    def->tooltip = L("Give to the infill algorithm if the infill needs to be connected, and on which perimeters"
        " Can be useful for art or with high infill/perimeter overlap."
        " The result may vary between infill types.");
    def->set_enum<InfillConnection>({
        { "connected", L("Connected") },
        { "holes", L("Connected to hole perimeters") },
        { "outershell", L("Connected to outer perimeters") },
        { "notconnected", L("Not connected") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillConnection>(icConnected));

    def = this->add("infill_connection_top", coEnum);
    def->label = L("Connection of top infill lines");
    def->category = OptionCategory::infill;
    def->tooltip = L("Give to the infill algorithm if the infill needs to be connected, and on which perimeters"
        " Can be useful for art or with high infill/perimeter overlap."
        " The result may vary between infill types.");
    def->set_enum<InfillConnection>({
        { "connected", L("Connected") },
        { "holes", L("Connected to hole perimeters") },
        { "outershell", L("Connected to outer perimeters") },
        { "notconnected", L("Not connected") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillConnection>(icConnected));

    def = this->add("infill_connection_bottom", coEnum);
    def->label = L("Connection of bottom infill lines");
    def->category = OptionCategory::infill;
    def->tooltip = L("Give to the infill algorithm if the infill needs to be connected, and on which perimeters"
        " Can be useful for art or with high infill/perimeter overlap."
        " The result may vary between infill types.");
    def->set_enum<InfillConnection>({
        { "connected", L("Connected") },
        { "holes", L("Connected to hole perimeters") },
        { "outershell", L("Connected to outer perimeters") },
        { "notconnected", L("Not connected") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillConnection>(icConnected));

    def = this->add("infill_connection_bridge", coEnum);
    def->label = L("Connection of bridged infill lines");
    def->category = OptionCategory::infill;
    def->tooltip = L("Give to the bridge infill algorithm if the infill needs to be connected, and on which perimeters."
        " Can be useful to disconnect to reduce a little bit the pressure buildup when going over the bridge's anchors.");
    def->set_enum<InfillConnection>({
        { "connected", L("Connected") },
        { "holes", L("Connected to hole perimeters") },
        { "outershell", L("Connected to outer perimeters") },
        { "notconnected", L("Not connected") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillConnection>(icNotConnected));

    def = this->add("infill_connection_solid", coEnum);
    def->label = L("Connection of solid infill lines");
    def->category = OptionCategory::infill;
    def->tooltip = L("Give to the infill algorithm if the infill needs to be connected, and on which perimeters"
        " Can be useful for art or with high infill/perimeter overlap."
        " The result may vary between infill types.");
    def->set_enum<InfillConnection>({
        { "connected", L("Connected") },
        { "holes", L("Connected to hole perimeters") },
        { "outershell", L("Connected to outer perimeters") },
        { "notconnected", L("Not connected") },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillConnection>(icConnected));

    def = this->add("infill_dense", coBool);
    def->label = L("Dense infill layer");
    def->full_label = L("Dense infill layer");
    def->category = OptionCategory::infill;
    def->tooltip = L("Enables the creation of a support layer under the first solid layer. This allows you to use a lower infill ratio without compromising the top quality."
        " The dense infill is laid out with a 50% infill density.");
    def->mode = comSimpleAE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));
    
    def = this->add("infill_dense_algo", coEnum);
    def->label = L("Algorithm");
    def->full_label = L("Dense infill algorithm");
    def->category = OptionCategory::infill;
    def->tooltip = L("Choose the way the dense layer is laid out."
        " The automatic option lets it try to draw the smallest surface with only strait lines inside the sparse infill."
        " The Anchored option just slightly enlarges (by 'Default infill margin') the surfaces that need a better support.");
    def->set_enum<DenseInfillAlgo>({
        { "automatic", L("Automatic") },
        { "autonotfull", L("Automatic, unless full") },
        { "autosmall", L("Automatic, only for small areas") },
        { "autoenlarged", L("Automatic, or anchored if too big") },
        { "enlarged", L("Anchored") },
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<DenseInfillAlgo>(dfaAutoOrEnlarged));

    def = this->add("infill_extruder", coInt);
    def->label = L("Infill extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing infill.");
    def->min = 1;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("infill_extrusion_width", coFloatOrPercent);
    def->label = L("Infill");
    def->full_label = L("Infill width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for infill. "
        "If left as zero, default extrusion width will be used if set, otherwise 1.125 x nozzle diameter will be used. "
        "You may want to use fatter extrudates to speed up the infill and make your parts stronger. "
        "If expressed as percentage (for example 110%) it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("infill_extrusion_change_odd_layers", coFloatOrPercent);
    def->label = L("Infill");
    def->full_label = L("Infill spacing change on even layers");
    def->category = OptionCategory::width;
    def->tooltip = L("Change width on every even layer (and not on odd layers like the first one) for better overlap with adjacent layers and getting stringer shells. "
                     "Try values about +/- 0.1 with different sign."
                     "\nThis could be combined with extra permeters on even layers."
                     "\nWorks as absolute spacing or a % of the spacing."
                     "\nset 0 to disable");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(false, 0));

    def = this->add("infill_extrusion_spacing", coFloatOrPercent);
    def->label = L("Infill");
    def->full_label = L("Infill spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like First layer width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
         "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("infill_fan_speed", coInts);
    def->label = L("Internal Infill fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all Internal Infill moves"
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("infill_first", coBool);
    def->label = L("Infill before perimeters");
    def->category = OptionCategory::infill;
    def->tooltip = L("This option will switch the print order of perimeters and infill, making the latter first.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    //def = this->add("infill_only_where_needed", coBool);
    //def->label = L("Only infill where needed");
    //def->category = OptionCategory::infill;
    //def->tooltip = L("This option will limit infill to the areas actually needed for supporting ceilings "
    //               "(it will act as internal support material). If enabled, this slows down the G-code generation "
    //               "due to the multiple checks involved.");
    //def->mode = comAdvancedE | comPrusa;
    //def->set_default_value(new ConfigOptionBool(false));

    def = this->add("infill_overlap", coFloatOrPercent);
    def->label = L("Infill/perimeters encroachment");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting applies an additional overlap between infill and perimeters for better bonding. "
                   "Theoretically this shouldn't be needed, but backlash might cause gaps. If expressed "
                   "as percentage (example: 15%) it is calculated over perimeter extrusion width."
                    "\nDon't put a value higher than 50% (of the perimeter width), as it will fuse with it and follow the perimeter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 0.5, true };
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(25, true));

    def = this->add("infill_speed", coFloatOrPercent);
    def->label = L("Sparse");
    def->full_label = L("Sparse infill speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing the internal fill."
        "\nThis can be expressed as a percentage (for example: 80%) over the Solid Infill speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "solid_infill_speed";
    def->aliases = { "print_feed_rate", "infill_feed_rate" };
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(300, true));

    def = this->add("inherits", coString);
    def->label = L("Inherits profile");
    def->tooltip = L("Name of the profile, from which this profile inherits.");
    def->full_width = true;
    def->height = 5;
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    // The following value is to be stored into the project file (AMF, 3MF, Config ...)
    // and it contains a sum of "inherits" values over the print and filament profiles.
    def = this->add("inherits_cummulative", coStrings);
    def->set_default_value(new ConfigOptionStrings());
    def->mode = comPrusa;
    def->cli = ConfigOptionDef::nocli;

    def = this->add("interface_shells", coBool);
    def->label = L("Interface shells");
    def->tooltip = L("Force the generation of solid shells between adjacent materials/volumes. "
                   "Useful for multi-extruder prints with translucent materials or manual soluble "
                   "support material.");
    def->category = OptionCategory::perimeter;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("internal_bridge_acceleration", coFloatOrPercent);
    def->label = L("Internal bridges ");
    def->full_label = L("Internal bridges acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for internal bridges. "
                "\nCan be a % of the default acceleration"
                "\nSet zero to use bridge acceleration for internal bridges.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "bridge_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));
    def->aliases = { "bridge_internal_acceleration" };

    def = this->add("internal_bridge_fan_speed", coInts);
    def->label = L("Infill bridges fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all infill bridges. It won't slow down the fan if it's currently running at a higher speed."
        "\nSet to 0 to stop the fan."
        "\nIf disabled, Bridge fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));
    def->aliases = { "bridge_internal_fan_speed" };

    def = this->add("internal_bridge_speed", coFloatOrPercent);
    def->label = L("Internal bridges");
    def->full_label = L("Internal bridge speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing the bridges that support the top layer.\nCan be a % of the bridge speed.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "bridge_speed";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(150,true));
    def->aliases = { "bridge_speed_internal" };

    def = this->add("mmu_segmented_region_max_width", coFloat);
    def->label = L("Maximum width of a segmented region");
    def->tooltip = L("Maximum width of a segmented region. Zero disables this feature.");
    def->sidetext = L("mm (zero to disable)");
    def->min = 0;
    def->category = OptionCategory::mmsetup;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("mmu_segmented_region_interlocking_depth", coFloat);
    def->label = L("Interlocking depth of a segmented region");
    def->tooltip = L("Interlocking depth of a segmented region. It will be ignored if "
                       "\"mmu_segmented_region_max_width\" is zero or if \"mmu_segmented_region_interlocking_depth\""
                       "is bigger then \"mmu_segmented_region_max_width\". Zero disables this feature.");
    def->sidetext = L("mm (zero to disable)");
    def->min = 0;
    def->category = OptionCategory::mmsetup;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("ironing", coBool);
    def->label = L("Enable ironing");
    def->tooltip = L("Enable ironing of the top layers with the hot print head for smooth surface");
    def->category = OptionCategory::ironing;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("ironing_acceleration", coFloatOrPercent);
    def->label = L("Ironing");
    def->full_label = L("Ironing acceleration");
    def->category = OptionCategory::ironing;
    def->tooltip = L("This is the acceleration your printer will use for ironing. "
                "\nCan be a % of the top solid infill acceleration"
                "\nSet zero or 100% to use top solid infill acceleration for ironing.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "top_solid_infill_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("ironing_angle", coFloat);
    def->label = L("Ironing angle");
    def->category = OptionCategory::ironing;
    def->tooltip = L("Ironing post-process angle."
        "\nIf positive, the ironing will use this angle."
        "\nIf -1, it will use the fill angle."
        "\nIf lower than -1, it will use the fill angle minus this angle.");
    def->sidetext = L("°");
    def->min = -360;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(-45));

    def = this->add("ironing_type", coEnum);
    def->label = L("Ironing Type");
    def->category = OptionCategory::ironing;
    def->tooltip = L("Ironing Type");
    def->set_enum<IroningType>({
        { "top",        L("All top surfaces") },
        { "topmost",    L("Topmost surface only") },
        { "solid",      L("All solid surfaces") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<IroningType>(IroningType::TopSurfaces));

    def = this->add("ironing_flowrate", coPercent);
    def->label = L("Flow rate");
    def->category = OptionCategory::ironing;
    def->tooltip = L("Percent of a flow rate relative to object's normal layer height."
                " It's the percentage of the layer that will be over-extruded on top to do the ironing.");
    def->sidetext = L("%");
    def->ratio_over = "layer_height";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionPercent(15));

    def = this->add("ironing_spacing", coFloatOrPercent);
    def->label = L("Spacing between ironing lines");
    def->category = OptionCategory::ironing;
    def->tooltip = L("Distance between ironing lines."
                    "\nCan be a % of the nozzle diameter used for ironing.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(25, true));

    def = this->add("ironing_speed", coFloatOrPercent);
    def->label = L("Ironing");
    def->full_label = L("Ironing speed");
    def->category = OptionCategory::ironing;
    def->tooltip = L("Ironing speed. Used for the ironing pass of the ironing infill pattern, and the post-process infill."
        "\nThis can be expressed as a percentage (for example: 80%) over the Top Solid Infill speed."
        "\nIroning extrusions are ignored from the automatic volumetric speed computation.");
    def->sidetext = L("mm/s");
    def->ratio_over = "top_solid_infill_speed";
    def->min = 0.1;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("layer_gcode", coString);
    def->label = L("After layer change G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This custom code is inserted at every layer change, right after the Z move "
        "and before the extruder moves to the first layer point. Note that you can use "
        "placeholder variables for all Slic3r settings as well as {layer_num} and {layer_z}.");
    def->cli = "after-layer-gcode|layer-gcode";
    def->multiline = true;
    def->full_width = true;
    def->height = 5;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("feature_gcode", coString);
    def->label = L("After layer change G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This custom code is inserted at every extrusion type change."
        "Note that you can use placeholder variables for all Slic3r settings as well as {last_extrusion_role}, {extrusion_role}, {layer_num} and {layer_z}."
        " The 'extrusion_role' strings can take these string values:"
        " { Perimeter, ExternalPerimeter, OverhangPerimeter, InternalInfill, SolidInfill, TopSolidInfill, BridgeInfill, GapFill, Skirt, SupportMaterial, SupportMaterialInterface, WipeTower, Mixed }."
        " Mixed is only used when the role of the extrusion is not unique, not exactly inside another category or not known.");
    def->cli = "feature-gcode";
    def->multiline = true;
    def->full_width = true;
    def->height = 5;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionString(""));

    // def = this->add("exact_last_layer_height", coBool);
    // def->label = L("Exact last layer height");
    // def->category = OptionCategory::perimeter;
    // def->tooltip = L("This setting controls the height of last object layers to put the last layer at the exact highest height possible. Experimental.");
    // def->mode = comHidden;
    // def->set_default_value(new ConfigOptionBool(false));

    def = this->add("remaining_times", coBool);
    def->label = L("Supports remaining times");
    def->category = OptionCategory::firmware;
    def->tooltip = L("Emit something at 1 minute intervals into the G-code to let the firmware show accurate remaining time.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("remaining_times_type", coEnum);
    def->label = L("Method");
    def->full_label = L("Supports remaining times method");
    def->category = OptionCategory::firmware;
    def->tooltip = L("M73: Emit M73 P{percent printed} R{remaining time in minutes} at 1 minute"
        " intervals into the G-code to let the firmware show accurate remaining time."
        " As of now only the Prusa i3 MK3 firmware recognizes M73."
        " Also the i3 MK3 firmware supports M73 Qxx Sxx for the silent mode."
        "\nM117: Send a command to display a message to the printer, this is 'Time Left .h..m..s'." );
    def->mode = comExpert | comSuSi;
    def->set_enum<RemainingTimeType>({
        { "m117", L("M117") },
        { "m73", L("M73") },
        { "m73m117", L("M73 & M117") },
    });
    def->set_default_value(new ConfigOptionEnum<RemainingTimeType>(RemainingTimeType::rtM73));

    def = this->add("silent_mode", coBool);
    def->label = L("Supports stealth mode");
    def->category = OptionCategory::firmware;
    def->tooltip = L("The firmware supports stealth mode");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("fan_speedup_time", coFloat);
    def->label = L("Fan startup delay");
    def->category = OptionCategory::firmware;
    def->tooltip = L("Move the fan start in the past by at least this delay (in seconds, you can use decimals)."
        " It assumes infinite acceleration for this time estimation, and will only take into account G1 and G0 moves."
        "\nIt won't move fan comands from custom gcodes (they act as a sort of 'barrier')."
        "\nIt won't move fan comands into the start gcode if the 'only custom start gcode' is activated."
        "\nUse 0 to deactivate.");
    def->sidetext = L("s");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("fan_speedup_overhangs", coBool);
    def->label = L("Fan delay only for overhangs");
    def->category = OptionCategory::firmware;
    def->tooltip = L("Will only take into account the delay for the cooling of overhangs.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("binary_gcode", coBool);
    def->label = L("Supports binary G-code");
    def->tooltip = L("Enable, if the firmware supports binary G-code format (bgcode). "
                     "To generate .bgcode files, make sure you have binary G-code enabled in Configuration->Preferences->Other.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("fan_kickstart", coFloat);
    def->label = L("Fan KickStart time");
    def->category = OptionCategory::firmware;
    def->tooltip = L("Add a M106 S255 (max speed for fan) for this amount of seconds before going down to the desired speed to kick-start the cooling fan."
                    "\nThis value is used for a 0->100% speedup, it will go down if the delta is lower."
                    "\nSet to 0 to deactivate.");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("lift_min", coFloat);
    def->label = L("Min height for travel");
    def->category = OptionCategory::extruders;
    def->tooltip = L("When an extruder travels to an object (from the start position or from an object to another), the nozzle height is guaranteed to be at least at this value."
        "\nIt's made to ensure the nozzle won't hit clips or things you have on your bed. But be careful to not put a clip in the 'convex shape' of an object."
        "\nSet to 0 to disable.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("machine_limits_usage", coEnum);
    def->label = L("How to apply limits");
    def->full_label = L("Purpose of Machine Limits");
    def->category = OptionCategory::limits;
    def->tooltip = L("How to apply the Machine Limits."
                    "\n* In every case, they will be used as safeguards: Even if you use a print profile that sets an acceleration of 5000,"
                    " if in your machine limits the acceleration is 4000, the outputted gcode will use the 4000 limit."
                    "\n* You can also use it as a safeguard and to have a better printing time estimate."
                    "\n* You can also use it as a safeguard, to have a better printing time estimate and emit the limits at the begining of the gcode file, with M201 M202 M203 M204 and M205 commands."
                    " If you want only to write a sub-set, choose the 'for time estimate' option and write your own gcodes in the custom gcode section.");
    def->set_enum<MachineLimitsUsage>({
        { "emit_to_gcode",      L("Also emit limits to G-code") },
        { "time_estimate_only", L("Use also for time estimate") },
        { "limits",             L("Use only as safeguards") },
        { "ignore",             L("Disable") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<MachineLimitsUsage>(MachineLimitsUsage::TimeEstimateOnly));

    {
        struct AxisDefault {
            std::string         name;
            std::vector<double> max_feedrate;
            std::vector<double> max_acceleration;
            std::vector<double> max_jerk;
        };
        std::vector<AxisDefault> axes {
            // name, max_feedrate,  max_acceleration, max_jerk
            { "x", { 500., 200. }, {  9000., 1000. }, { 10. , 10.  } },
            { "y", { 500., 200. }, {  9000., 1000. }, { 10. , 10.  } },
            { "z", {  12.,  12. }, {   500.,  200. }, {  0.2,  0.4 } },
            { "e", { 120., 120. }, { 10000., 5000. }, {  2.5,  2.5 } }
        };
        for (const AxisDefault &axis : axes) {
            std::string axis_upper = boost::to_upper_copy<std::string>(axis.name);
            // Add the machine feedrate limits for XYZE axes. (M203)
            def = this->add("machine_max_feedrate_" + axis.name, coFloats);
            def->full_label = (boost::format("Maximum feedrate %1%") % axis_upper).str();
            (void)L("Maximum feedrate X");
            (void)L("Maximum feedrate Y");
            (void)L("Maximum feedrate Z");
            (void)L("Maximum feedrate E");
            def->category = OptionCategory::limits;
            def->tooltip  = (boost::format("Maximum feedrate of the %1% axis") % axis_upper).str();
            (void)L("Maximum feedrate of the X axis");
            (void)L("Maximum feedrate of the Y axis");
            (void)L("Maximum feedrate of the Z axis");
            (void)L("Maximum feedrate of the E axis");
            def->sidetext = L("mm/s");
            def->min = 0;
            def->mode = comAdvancedE | comPrusa;
            def->set_default_value(new ConfigOptionFloats(axis.max_feedrate));
            // Add the machine acceleration limits for XYZE axes (M201)
            def = this->add("machine_max_acceleration_" + axis.name, coFloats);
            def->full_label = (boost::format("Maximum acceleration %1%") % axis_upper).str();
            (void)L("Maximum acceleration X");
            (void)L("Maximum acceleration Y");
            (void)L("Maximum acceleration Z");
            (void)L("Maximum acceleration E");
            def->category = OptionCategory::limits;
            def->tooltip  = (boost::format("Maximum acceleration of the %1% axis") % axis_upper).str();
            (void)L("Maximum acceleration of the X axis");
            (void)L("Maximum acceleration of the Y axis");
            (void)L("Maximum acceleration of the Z axis");
            (void)L("Maximum acceleration of the E axis");
            def->sidetext = L("mm/s²");
            def->min = 0;
            def->mode = comAdvancedE | comPrusa;
            def->set_default_value(new ConfigOptionFloats(axis.max_acceleration));
            // Add the machine jerk limits for XYZE axes (M205)
            def = this->add("machine_max_jerk_" + axis.name, coFloats);
            def->full_label = (boost::format("Maximum jerk %1%") % axis_upper).str();
            (void)L("Maximum jerk X");
            (void)L("Maximum jerk Y");
            (void)L("Maximum jerk Z");
            (void)L("Maximum jerk E");
            def->category = OptionCategory::limits;
            def->tooltip  = (boost::format("Maximum jerk of the %1% axis") % axis_upper).str();
            (void)L("Maximum jerk of the X axis");
            (void)L("Maximum jerk of the Y axis");
            (void)L("Maximum jerk of the Z axis");
            (void)L("Maximum jerk of the E axis");
            def->sidetext = L("mm/s");
            def->min = 0;
            def->mode = comAdvancedE | comPrusa;
            def->set_default_value(new ConfigOptionFloats(axis.max_jerk));
        }
    }

    // M205 S... [mm/sec]
    def = this->add("machine_min_extruding_rate", coFloats);
    def->full_label = L("Minimum feedrate when extruding");
    def->category = OptionCategory::limits;
    def->tooltip = L("Minimum feedrate when extruding (M205 S)");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloats{ 0., 0. });

    // M205 T... [mm/sec]
    def = this->add("machine_min_travel_rate", coFloats);
    def->full_label = L("Minimum travel feedrate");
    def->category = OptionCategory::limits;
    def->tooltip = L("Minimum travel feedrate (M205 T)");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloats{ 0., 0. });

    // M204 P... [mm/sec^2]
    def = this->add("machine_max_acceleration_extruding", coFloats);
    def->full_label = L("Maximum acceleration when extruding");
    def->category = OptionCategory::limits;
    def->tooltip = L("Maximum acceleration when extruding");
    def->sidetext = L("mm/s²");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloats{ 1500., 1250. });

    // M204 R... [mm/sec^2]
    def = this->add("machine_max_acceleration_retracting", coFloats);
    def->full_label = L("Maximum acceleration when retracting");
    def->category = OptionCategory::limits;
    def->tooltip = L("Maximum acceleration when retracting.\n\n"
                     "Not used for RepRapFirmware, which does not support it.");
    def->sidetext = L("mm/s²");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloats{ 1500., 1250. });

    // M204 T... [mm/sec^2]
    def = this->add("machine_max_acceleration_travel", coFloats);
    def->full_label = L("Maximum acceleration for travel moves");
    def->category = OptionCategory::limits;
    def->tooltip = L("Maximum acceleration for travel moves.");
    def->sidetext = L("mm/s²");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloats{ 1500., 1250. });

    def = this->add("max_gcode_per_second", coFloat);
    def->label = L("Maximum G1 per second");
    def->category = OptionCategory::speed;
    def->tooltip = L("If your firmware stops while printing, it may have its gcode queue full."
        " Set this parameter to merge extrusions into bigger ones to reduce the number of gcode commands the printer has to process each second."
        "\nOn 8bit controlers, a value of 150 is typical."
        "\nNote that reducing your printing speed (at least for the external extrusions) will reduce the number of time this will triggger and so increase quality."
        "\nSet zero to disable.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0/*1500*/));

    def = this->add("max_fan_speed", coInts);
    def->label = L("Max");
    def->full_label = L("Max fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This setting represents the maximum speed of your fan, used when the layer print time is Very short.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 100 });

    def = this->add("max_layer_height", coFloatsOrPercents);
    def->label = L("Max");
    def->full_label = L("Max layer height");
    def->category = OptionCategory::general;
    def->tooltip = L("This is the highest printable layer height for this extruder, used to cap "
                   "the variable layer height and support layer height. Maximum recommended layer height "
                   "is 75% of the extrusion width to achieve reasonable inter-layer adhesion. "
                   "\nCan be a % of the nozzle diameter."
                   "\nIf disabled, layer height is limited to 75% of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 1, true };
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloatsOrPercents{ FloatOrPercent{ 75, true} }, false));

    def = this->add("max_print_speed", coFloatOrPercent);
    def->label = L("Max auto-speed");
    def->full_label = L("Max print speed for Autospeed");
    def->category = OptionCategory::speed;
    def->tooltip = L("When setting other speed settings to 0, Slic3r will autocalculate the optimal speed "
        "in order to keep constant extruder pressure. This experimental setting is used "
        "to set the highest print speed you want to allow."
        "\nThis can be expressed as a percentage (for example: 100%) over the machine Max Feedrate for X axis.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "machine_max_feedrate_x";
    def->min = 1;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(80, false));

    def = this->add("max_speed_reduction", coPercents);
    def->label = L("Max speed reduction");
    def->category = OptionCategory::speed;
    def->tooltip = L("This setting control by how much the speed can be reduced to increase the layer time."
        " It's a maximum reduction, so a lower value makes the minimum speed higher."
        " Set to 90% if you don't want the speed to go below 10% of the current speed."
        "\nSet zero to disable");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents{ 90 });

    def = this->add("max_volumetric_speed", coFloat);
    def->label = L("Volumetric auto-speed");
    def->full_label = L("Maximum Volumetric print speed for Autospeed");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This setting allows you to set the maximum flowrate for your print, and so cap the desired flow rate for the autospeed algorithm."
        " The autospeed tries to keep a constant feedrate for the entire object, and so can lower the volumetric speed for some features."
        "\nThe autospeed is only enable on speed fields that have a value of 0. If a speed field is a % of a 0 field, then it will be a % of the value it should have got from the autospeed."
        "\nIf this field is set to 0, then there is no autospeed nor maximum flowrate. If a speed value i still set to 0, it will get the max speed allwoed by the printer.");
    def->sidetext = L("mm³/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("max_volumetric_extrusion_rate_slope_positive", coFloat);
    def->label = L("Max volumetric slope positive");
    def->tooltip = L("This experimental setting is used to limit the speed of change in extrusion rate "
                       "for a transition from lower speed to higher speed. "
                   "A value of 1.8 mm³/s² ensures, that a change from the extrusion rate "
                   "of 1.8 mm³/s (0.45mm extrusion width, 0.2mm extrusion height, feedrate 20 mm/s) "
                   "to 5.4 mm³/s (feedrate 60 mm/s) will take at least 2 seconds.");
    def->sidetext = L("mm³/s²");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("max_volumetric_extrusion_rate_slope_negative", coFloat);
    def->label = L("Max volumetric slope negative");
    def->tooltip = L("This experimental setting is used to limit the speed of change in extrusion rate "
                       "for a transition from higher speed to lower speed. "
                   "A value of 1.8 mm³/s² ensures, that a change from the extrusion rate "
                   "of 5.4 mm³/s (0.45 mm extrusion width, 0.2 mm extrusion height, feedrate 60 mm/s) "
                   "to 1.8 mm³/s (feedrate 20 mm/s) will take at least 2 seconds.");
    def->sidetext = L("mm³/s²");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

#if 0
    // replaced by fan_printer_min_speed, to remove !! @DEPRECATED
    def = this->add("min_fan_speed", coInts);
    def->label = L("Default fan speed");
    def->full_label = L("Default fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This setting represents the base fan speed this filament needs, or at least the minimum PWM your fan needs to work.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts{ 35 });
#endif

    def = this->add("fan_percentage", coBool);
    def->label = L("Fan PWM from 0-100");
    def->category = OptionCategory::output;
    def->tooltip = L("Set this if your printer uses control values from 0-100 instead of 0-255.");
    def->cli = "fan-percentage";
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("min_layer_height", coFloatsOrPercents);
    def->label = L("Min");
    def->full_label = L("Min layer height");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This is the lowest printable layer height for this extruder and limits "
        "the resolution for variable layer height. Typical values are between 0.05 mm and 0.1 mm."
        "\nCan be a % of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 5, false };
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{ 5, true} });

    def = this->add("min_width_top_surface", coFloatOrPercent);
    def->label = L("Minimum top width for infill");
    def->category = OptionCategory::speed;
    def->tooltip = L("If a top surface has to be printed and it's partially covered by another layer, it won't be considered at a top layer where its width is below this value."
        " This can be useful to not let the 'one perimeter on top' trigger on surface that should be covered only by perimeters."
        " This value can be a mm or a % of the perimeter extrusion width."
        "\nWarning: If enabled, artifacts can be created is you have some thin features on the next layer, like letters. Set this setting to 0 to remove these artifacts.");
    def->sidetext = L("mm or %");
    def->ratio_over = "perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 15, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("min_print_speed", coFloats);
    def->label = L("Min print speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Slic3r will never scale the speed below this one.");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 10. });

    def = this->add("min_skirt_length", coFloat);
    def->label = L("Minimal filament extrusion length");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Generate no less than the number of skirt loops required to consume "
                   "the specified amount of filament on the bottom layer. For multi-extruder machines, "
                   "this minimum applies to each extruder.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("model_precision", coFloat);
    def->label = L("Model rounding precision");
    def->full_label = L("Model rounding precision");
    def->category = OptionCategory::slicing;
    def->tooltip = L("This is the rounding error of the input object."
        " It's used to align points that should be in the same line."
        "\nSet zero to disable.");
    def->sidetext = L("mm");
    def->min = 0;
    def->precision = 8;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0.0001));

    def = this->add("notes", coString);
    def->label = L("Configuration notes");
    def->category = OptionCategory::notes;
    def->tooltip = L("Here you can put your personal notes. This text will be added to the G-code "
                   "header comments.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("nozzle_diameter", coFloats);
    def->label = L("Nozzle diameter");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This is the diameter of your extruder nozzle (for example: 0.5, 0.35 etc.)");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0.4 });

    def = this->add("host_type", coEnum);
    def->label = L("Host Type");
    def->category = OptionCategory::general;
    def->tooltip = L("Slic3r can upload G-code files to a printer host. This field must contain "
                   "the kind of the host."
                   "\nPrusaLink is only available for prusa printer.");
    def->set_enum<PrintHostType>({
        { "prusalink",      "PrusaLink" },
        { "prusaconnect",   "PrusaConnect" },
        { "octoprint",      "OctoPrint" },
        { "moonraker",      "Klipper (via Moonraker)" },
        { "duet",           "Duet" },
        { "flashair",       "FlashAir" },
        { "astrobox",       "AstroBox" },
        { "repetier",       "Repetier" },
        { "klipper",        "Klipper" },
        { "mpmdv2",         "MPMDv2" },
        { "mks",            "MKS" },
        { "monoprice",      "Monoprice lcd" },
    });
    def->mode = comAdvancedE | comPrusa;
    def->cli = ConfigOptionDef::nocli;
    def->set_default_value(new ConfigOptionEnum<PrintHostType>(htPrusaLink));

    def = this->add("print_custom_variables", coString);
    def->label = L("Custom variables");
    def->full_label = L("Custom Print variables");
    def->category = OptionCategory::filament;
    def->tooltip = L("You can add data accessible to custom-gcode macros."
        "\nEach line can define one variable."
        "\nThe format is 'variable_name=value'. the variable name should only have [a-zA-Z0-9] characters or '_'."
        "\nA value that can be parsed as a int or float will be avaible as a numeric value."
        "\nA value that is enclosed by double-quotes will be available as a string (without the quotes)"
        "\nA value that only takes values as 'true' or 'false' will be a boolean)"
        "\nEvery other value will be parsed as a string as-is."
        "\nAdvice: before using a variable, it's safer to use the function 'default_XXX(variable_name, default_value)'"
        " (enclosed in bracket as it's a script) in case it's not set. You can replace XXX by 'int' 'bool' 'double' 'string'.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionString{ "" });

    def = this->add("object_gcode", coString);
    def->label = L("Per object G-code");
    def->category = OptionCategory::advanced;
    def->tooltip = L("This code is inserted each layer, when the object began to print (just after the label if any)."
                     " It's main advantage is when you use it as a object modifer (right click on a model)."
                     "\nSpecial variables: 'layer_num','layer_z'");
    def->multiline = true;
    def->full_width = true;
    def->height = 10;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("only_one_perimeter_first_layer", coBool);
    def->label = L("Only one perimeter on First layer");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Use only one perimeter on first layer, to give more space to the top infill pattern.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("only_one_perimeter_top", coBool);
    def->label = L("Only one perimeter on Top surfaces");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Use only one perimeter on flat top surface, to give more space to the top infill pattern.");
    def->mode = comSimpleAE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("only_one_perimeter_top_other_algo", coBool);
    def->label = L("Only one peri - other algo");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("If you have some problem with the 'Only one perimeter on Top surfaces' option, you can try to activate this on the problematic layer.");
    def->mode = comHidden;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("only_retract_when_crossing_perimeters", coBool);
    def->label = L("Only retract when crossing perimeters");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Disables retraction when the travel path does not exceed the upper layer's perimeters "
        "(and thus any ooze will probably be invisible).");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("ooze_prevention", coBool);
    def->label = L("Enable");
    // TRN PrintSettings: Enable ooze prevention
    def->tooltip = L("This option will drop the temperature of the inactive extruders to prevent oozing.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("output_filename_format", coString);
    def->label = L("Output filename format");
    def->category = OptionCategory::output;
    def->tooltip = L("You can use all configuration options as variables inside this template. "
                   "For example: {layer_height}, {fill_density} etc. You can also use {timestamp}, "
                   "{year}, {month}, {day}, {hour}, {minute}, {second}, {version}, {input_filename}, "
                   "{input_filename_base}, {default_output_extension}.");
    def->full_width = true;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString("{input_filename_base}.gcode"));

    def = this->add("overhangs_acceleration", coFloatOrPercent);
    def->label = L("Overhangs");
    def->full_label = L("Overhang acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for overhangs."
                "\nCan be a % of the bridge acceleration"
                "\nSet zero to to use bridge acceleration for overhangs.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "bridge_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("overhangs_bridge_threshold", coFloat);
    def->label = L("Bridge max length");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Maximum distance for bridges. If the distance is over that, it will be considered as overhangs for 'overhangs_max_slope'."
                    "\nIf disabled, accept all distances."
                    "\nSet to 0 to ignore bridges.");
    def->sidetext = L("mm");
    def->min = 0;
    def->can_be_disabled = true;
    def->mode = comExpert | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloat(0)));

    def = this->add("overhangs_bridge_upper_layers", coInt);
    def->label = L("Consider upper bridges");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Don't put overhangs if the area will filled in next layer by bridges."
                    "\nIf disabled, accept all upper layers."
                    "\nSet to 0 to only consider our layer bridges.");
    def->sidetext = L("layers");
    def->min = 0;
    def->can_be_disabled = true;
    def->mode = comExpert | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionInt(2), false));

    def             = this->add("overhangs_dynamic_fan_speed", coGraphs);
    def->label      = L("Dynamic overhang speeds");
    def->category   = OptionCategory::speed;
    def->tooltip    = L("This setting can only works correctly if dynamic speed is also enabled (overhangs_dynamic_fan_speed)."
        "\nOverhang size is expressed as a percentage of overlap of the extrusion with the previous layer: "
        "100% would be full overlap (no overhang), while 0% represents full overhang (floating extrusion, bridge)."
        "\nFan speeds for overhang sizes in between are calculated via linear interpolation."
        "\nIf enabled, overhangs_fan_speed is disabled, as the fan speed for full overhang is used.");
    def->sidetext   = L("%");
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->mode       = comExpert | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionGraphs({GraphData(0,4, GraphData::GraphType::LINEAR,
        {{0,100},{25,80},{50,60},{75,40},{100,20}}
    )})));
    def->graph_settings = std::make_shared<GraphSettings>();
    def->graph_settings->title       = L("Overhangs fan speed by % of overlap");
    def->graph_settings->description = L("Choose the Overhangs maximu fan speed for each percentage of overlap with the layer below."
        "If the current fan speed (from perimeter, external, of default) is higher, then this setting won't slow the fan."
        "\n100% overlap is when the extrusion is fully on top of the previous layer's extrusion."
        "\n0% overlap is when the extrusion centerline is at a distance of 'overhangs threshold for speed'(overhangs_bridge_threshold)"
        "\nfrom the nearest extrusion of the previous layer.");
    def->graph_settings->x_label     = L("overlap % with previous layer");
    def->graph_settings->y_label     = L("Fan speed (%)");
    def->graph_settings->null_label  = L("No fan speed");
    def->graph_settings->label_min_x = L("");
    def->graph_settings->label_max_x = L("");
    def->graph_settings->label_min_y = L("");
    def->graph_settings->label_max_y = L("");
    def->graph_settings->min_x       = 0;
    def->graph_settings->max_x       = 100;
    def->graph_settings->step_x      = 1.;
    def->graph_settings->min_y       = 0;
    def->graph_settings->max_y       = 100;
    def->graph_settings->step_y      = 1.;
    def->graph_settings->allowed_types = {GraphData::GraphType::LINEAR, GraphData::GraphType::SQUARE, GraphData::GraphType::SPLINE};

    def             = this->add("overhangs_dynamic_speed", coGraph);
    def->label      = L("Dynamic overhang speeds");
    def->category   = OptionCategory::speed;
    def->tooltip    = L("Overhang size is expressed as a percentage of overlap of the extrusion with the previous layer:"
                        " 100% would be full overlap (no overhang), while 0% represents full overhang (floating extrusion, bridge)."
                        " Speeds for overhang sizes in between are calculated via linear interpolation,"
                        " as a percentage between the (external) perimeter speed and the overhang speed."
                        "\nNote that the speeds generated to gcode will never exceed the max volumetric speed value.");
    def->sidetext   = L("mm/s");
    def->can_be_disabled = true;
    def->mode       = comExpert | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionGraph(GraphData(0,4, GraphData::GraphType::LINEAR,
        {{0,0},{25,10},{50,40},{75,70},{100,100}}
    ))));
    def->graph_settings = std::make_shared<GraphSettings>();
    def->graph_settings->title       = L("Overhangs speed ratio by % of overlap");
    def->graph_settings->description = L("Choose the Overhangs speed for each percentage of overlap with the layer below."
        "\nThe speed is a percentage ratio between overhangs speed (for 0% overlap) and perimeter / external perimeter speed (for 100% overlap)."
        "\n100% overlap is when the extrusion is fully on top of the previous layer's extrusion."
        "\n0% overlap is when the extrusion centerline is at a distance of 'overhangs threshold for speed'(overhangs_bridge_threshold)"
        "\nfrom the nearest extrusion of the previous layer.");
    def->graph_settings->x_label     = L("overlap % with previous layer");
    def->graph_settings->y_label     = L("Speed ratio (%)");
    def->graph_settings->null_label  = L("Uses overhangs speed");
    def->graph_settings->label_min_x = L("");
    def->graph_settings->label_max_x = L("");
    def->graph_settings->label_min_y = L("");
    def->graph_settings->label_max_y = L("");
    def->graph_settings->min_x       = 0;
    def->graph_settings->max_x       = 100;
    def->graph_settings->step_x      = 1.;
    def->graph_settings->min_y       = 0;
    def->graph_settings->max_y       = 100;
    def->graph_settings->step_y      = 1.;
    def->graph_settings->allowed_types = {GraphData::GraphType::LINEAR, GraphData::GraphType::SQUARE, GraphData::GraphType::SPLINE};

    def = this->add("overhangs_fan_speed", coInts);
    def->label = L("Overhangs Perimeter fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all Overhang Perimeter moves"
        "\nIf disabled, the previous (perimeter) fan speed will be used."
        "\nCan be overriden by disable_fan_first_layers and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("overhangs_max_slope", coFloatOrPercent);
    def->label = L("Overhangs max slope");
    def->full_label = L("Overhangs max slope");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Maximum slope for overhangs. if at each layer, the overhangs hangs by more than this value, then the geometry will be cut."
                    " It doesn't cut into detected bridgeable areas."
                    "\nCan be a % of the highest nozzle diameter."
                    "\nSet to 0 to disable.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("overhangs_speed", coFloatOrPercent);
    def->label = L("Overhangs");
    def->full_label = L("Overhangs speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing overhangs."
        "\nCan be a % of the bridge speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s");
    def->ratio_over = "bridge_speed";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("overhangs_speed_enforce", coInt);
    def->label = L("Enforce overhangs speed");
    def->full_label = L("Enforce overhangs speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Set the speed of the full perimeters to the overhang speed, and also the next one(s) if any."
                "\nSet to 0 to disable."
                "\nSet to 1 to set the overhang speed to the full perimeter if there is any overhang detected inside it."
                "\nSet to more than 1 to also set the overhang speed to the next perimeter(s) (only in classic mode)."
                );
    def->sidetext = L("perimeters");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("overhangs_width_speed", coFloatOrPercent);
    def->label = L("'As bridge' speed threshold");
    def->full_label = L("Overhang bridge speed threshold");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Minimum unsupported width for an extrusion to apply the bridge fan & overhang speed to it."
        "\nCan be in mm or in a % of the nozzle diameter."
        "\nCan be overriden by the overhang flow threshold if its value lower than this threshold."
        "\nIf dynamic speed is used, then the dynamic speed will be computed between 0% and this threshold.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->can_be_disabled = true;
    def->mode = comExpert | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloatOrPercent(55,true), false));

    def = this->add("overhangs_width", coFloatOrPercent);
    def->label = L("'As bridge' flow threshold");
    def->full_label = L("Overhang bridge flow threshold");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Minimum unsupported width for an extrusion to apply the overhang bridge flow to it."
        "\nCan be in mm or in a % of the nozzle diameter."
        "\nIf lower than the threshold for overhangs speed, then this threshold is used for both."
        "\nIf dynamic speed is used, and the overhangs speed threshold isn't enabled or is higher than this one,"
        " then the dynamic speed will be computed between 0% and this threshold.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 10, true };
    def->can_be_disabled = true;
    def->mode = comExpert | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionFloatOrPercent(75, true), false));

    def = this->add("overhangs_reverse", coBool);
    def->label = L("Reverse on even");
    def->full_label = L("Overhang reversal on even layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Extrude perimeters that have an overhanging part in the reverse direction on even layers (not on odd layers like the first one)."
        " This alternating pattern can significantly improve steep overhangs."
        "\n!! this is a very slow algorithm !!");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("overhangs_reverse_threshold", coFloatOrPercent);
    def->label = L("Reverse threshold");
    def->full_label = L("Overhang reversal threshold");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Number of mm the overhang need to be for the reversal to be considered useful. Can be a % of the perimeter width.");
    def->ratio_over = "perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 20, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(250, true));

    def = this->add("no_perimeter_unsupported_algo", coEnum);
    def->label = L("No perimeters on bridge areas");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Experimental option to remove perimeters where there is nothing under them and where a bridged infill should be better. "
        "\n * Remove perimeters: remove the unsupported perimeters, leave the bridge area as-is."
        "\n * Keep only bridges: remove the perimeters in the bridge areas, keep only bridges that end in solid area."
        "\n * Keep bridges and overhangs: remove the unsupported perimeters, keep only bridges that end in solid area, fill the rest with overhang perimeters+bridges."
        "\n * Fill the voids with bridges: remove the unsupported perimeters, draw bridges over the whole hole.*"
        " !! this one can escalate to problems with overhangs shaped like  /\\, so you should use it only on one layer at a time via the height-range modifier!"
        "\n!!Computationally intensive!!. ");
    def->set_enum<NoPerimeterUnsupportedAlgo>({
        { "none",             "Disabled" },
        { "noperi",           "Remove perimeters" },
        { "bridges",          "Keep only bridges" },
        { "bridgesoverhangs", "Keep bridges and overhangs" },
        { "filled",           "Fill the voids with bridges" },
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<NoPerimeterUnsupportedAlgo>(npuaNone));

    def = this->add("parking_pos_retraction", coFloat);
    def->label = L("Filament parking position");
    def->tooltip = L("Distance of the extruder tip from the position where the filament is parked "
                      "when unloaded. This should match the value in printer firmware. ");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(92.));

    def = this->add("extra_loading_move", coFloat);
    def->label = L("Extra loading distance");
    def->tooltip = L("When set to zero, the distance the filament is moved from parking position during load "
                      "is exactly the same as it was moved back during unload. When positive, it is loaded further, "
                      " if negative, the loading move is shorter than unloading. ");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(-2.));

    def = this->add("perimeter_acceleration", coFloatOrPercent);
    def->label = L("Internal");
    def->full_label = L("Internal Perimeter acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for internal perimeters. "
                "\nCan be a % of the default acceleration"
                "\nSet zero to use default acceleration for internal perimeters.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "default_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("perimeter_bonding", coPercent);
    def->label = L("Better bonding");
    def->full_label = L("Perimeter bonding");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This setting may slightly degrade the quality of your external perimeter, in exchange for a better bonding between perimeters."
        "Use it if you have great difficulties with perimeter bonding, for example with high temperature filaments."
        "\nThis percentage is the % of overlap between perimeters, a bit like perimeter_overlap and external_perimeter_overlap, but in reverse."
        " You have to set perimeter_overlap and external_perimeter_overlap to 100%, or this setting has no effect."
        " 0: no effect, 50%: half of the nozzle will be over an already extruded perimeter while extruding a new one"
        ", unless it's an external one)."
        "\nNote: it needs the external and perimeter overlap to be at 100% and to print the external perimeter first.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 50;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(0));

    def = this->add("perimeter_extruder", coInt);
    def->label = L("Perimeter extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing perimeters and brim. First extruder is 1.");
    def->aliases = { "perimeters_extruder" };
    def->min = 1;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("perimeter_extrusion_width", coFloatOrPercent);
    def->label = L("Perimeters");
    def->full_label = L("Perimeter width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for perimeters. "
        "You may want to use thinner extrudates to get more accurate surfaces. "
        "If left zero, default extrusion width will be used if set, otherwise 1.125 x nozzle diameter will be used. "
        "If expressed as percentage (for example 105%) it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->aliases = { "perimeters_extrusion_width" };
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("perimeter_extrusion_spacing", coFloatOrPercent);
    def->label = L("Perimeters");
    def->full_label = L("Perimeter spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like Perimeter width but spacing is the distance between two perimeter lines (as they overlap a bit, it's not the same)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using the perimeter 'Overlap' percentages and default layer height.");
    def->sidetext = L("mm or %");
    def->aliases = { "perimeters_extrusion_width" };
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("perimeter_extrusion_change_odd_layers", coFloatOrPercent);
    def->label = L("Perimeters");
    def->full_label = L("Perimeters spacing change on even layers");
    def->category = OptionCategory::width;
    def->tooltip = L("Change width on every even layer (and not on odd layers like the first one) for better overlap with adjacent layers and getting stringer shells. "
                     "Try values about +/- 0.1 with different sign for external and internal perimeters."
                     "\nThis could be combined with extra permeters on even layers."
                     "\nWorks as absolute spacing or a % of the spacing."
                     "\nset 0 to disable");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("perimeter_direction", coEnum);
    def->label = L("Perimeter direction");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Default direction to print the perimeters (contours and holes): clockwise (CW) or counter-clockwise (CCW).");
    def->set_enum<PerimeterDirection>({
        { "ccw_cw", "Contour: CCW, Holes: CW" },
        { "ccw_ccw","Contour & holes: CCW" },
        { "cw_ccw", "Contour: CW, Holes: CCW" },
        { "cw_cw","Contour & holes: CW" },
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<PerimeterDirection>(pdCCW_CW));

    def = this->add("perimeter_fan_speed", coInts);
    def->label = L("Internal Perimeter fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all Perimeter moves"
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("perimeter_loop", coBool);
    def->label = L("Perimeters loop");
    def->full_label = L("Perimeters loop");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Join the perimeters to create only one continuous extrusion without any z-hop."
        " Long inside travel (from external to holes) are not extruded to give some space to the infill.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("perimeter_loop_seam", coEnum);
    def->label = L("Seam position");
    def->full_label = L("Perimeter loop seam");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Position of perimeters starting points.");
    def->set_enum<SeamPosition>({
        { "nearest", "Nearest" },
        { "rear",    "Rear" },
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<SeamPosition>(spRear));

    def = this->add("perimeter_overlap", coPercent);
    def->label = L("perimeter overlap");
    def->full_label = L("Perimeter overlap");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allows you to reduce the overlap between the perimeters, to reduce the impact of the perimeters' artifacts."
        " 100% means that no gap is left, and 0% means that perimeters are not touching each other anymore."
        "\nIt's very experimental, please report about the usefulness. It may be removed if there is no use for it.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("perimeter_reverse", coBool);
    def->label = L("Perimeter reversal on even layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("On even layers, all perimeter loops are reversed (it disables the overhang reversal, so it doesn't double-reverse)."
                    "That setting will likely create defects on the perimeters, so it's only useful is for materials that have some direction-dependent properties (stress lines).");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("perimeter_round_corners", coBool);
    def->label = L("Round corners");
    def->full_label = L("Round corners for perimeters");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Internal perimeters will go around sharp corners by turning around instead of making the same sharp corner."
                        " This can help when there are visible holes in sharp corners on internal perimeters."
                        "\nCan incur some more processing time, and corners are a bit less sharp.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("perimeter_speed", coFloatOrPercent);
    def->label = L("Internal");
    def->full_label = L("Internal perimeters speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for perimeters (contours, aka vertical shells)."
        "\nThis can be expressed as a percentage (for example: 80%) over the Default speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->aliases = { "perimeter_feed_rate" };
    def->ratio_over = "default_speed";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(60, true));

    def = this->add("perimeters", coInt);
    def->label = L("Perimeters");
    def->full_label = L("Perimeters count");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This option sets the number of perimeters to generate for each layer."
                   "\nIf perimeters_hole is activated, then this number is only for contour perimeters."
                   "Note that if a contour perimeter encounter a hole, it will go around like a hole perimeter."
                   "\nNote that Slic3r may increase this number automatically when it detects "
                   "sloping surfaces which benefit from a higher number of perimeters "
                   "if the Extra Perimeters option is enabled.");
    def->sidetext = L("(minimum).");
    def->aliases = { "perimeter_offsets" };
    def->min = 0;
    def->max = 10000;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(3));

    def = this->add("perimeters_hole", coInt);
    def->label = L("Max perimeter count for holes");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This option sets the number of perimeters to have over holes."
                   " Note that if a hole-perimeter fuse with the contour, then it will go around like a contour perimeter.."
                   "\nIf disabled, holes will have the same number of perimeters as contour."
                   "\nNote that Slic3r may increase this number automatically when it detects "
                   "sloping surfaces which benefit from a higher number of perimeters "
                   "if the Extra Perimeters option is enabled.");
    def->sidetext = L("(minimum).");
    def->min = 0;
    def->max = 10000;
    def->can_be_disabled = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(disable_defaultoption(new ConfigOptionInt(0)));

    def = this->add("post_process", coStrings);
    def->label = L("Post-processing scripts");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("If you want to process the output G-code through custom scripts, "
                   "just list their absolute paths here."
                   "\nSeparate multiple scripts with a semicolon or a line return.\n!! please use '\\;' here if you want a not-line-separation ';'!!"
                   "\nScripts will be passed the absolute path to the G-code file as the first argument, "
                   "and they can access the Slic3r config settings by reading environment variables."
                   "\nThe script, if passed as a relative path, will also be searched from the slic3r directory, "
                   "the slic3r configuration directory and the user directory.");
    def->multiline = true;
    def->full_width = true;
    def->height = 6;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionStrings());

    def = this->add("priming_position", coPoint);
    def->label = L("Priming position");
    def->full_label = L("Priming position");
    def->tooltip = L("Coordinates of the left front corner of the priming patch."
                     "\nIf set to 0,0 then the position is computed automatically.");
    //TODO: enable/disable
    def->category = OptionCategory::customgcode;
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionPoint(Vec2d(0,0)));

    def = this->add("printer_custom_variables", coString);
    def->label = L("Custom variables");
    def->full_label = L("Custom Printer variables");
    def->category = OptionCategory::filament;
    def->tooltip = L("You can add data accessible to custom-gcode macros."
        "\nEach line can define one variable."
        "\nThe format is 'variable_name=value'. the variable name should only have [a-zA-Z0-9] characters or '_'."
        "\nA value that can be parsed as a int or float will be avaible as a numeric value."
        "\nA value that is enclosed by double-quotes will be available as a string (without the quotes)"
        "\nA value that only takes values as 'true' or 'false' will be a boolean)"
        "\nEvery other value will be parsed as a string as-is."
        "\nAdvice: before using a variable, it's safer to use the function 'default_XXX(variable_name, default_value)'"
        " (enclosed in bracket as it's a script) in case it's not set. You can replace XXX by 'int' 'bool' 'double' 'string'.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionString{ "" });

    def = this->add("fan_printer_min_speed", coInt);
    def->label = L("Minimum fan speed");
    def->full_label = L("Minimum fan speed");
    def->category = OptionCategory::general;
    def->tooltip = L("This setting represents the minimum fan speed (like minimum PWM) your fan needs to work.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("printer_model", coString);
    def->label = L("Printer type");
    def->tooltip = L("Type of the printer.");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("printer_notes", coString);
    def->label = L("Printer notes");
    def->category = OptionCategory::notes;
    def->tooltip = L("You can put your notes regarding the printer here.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("printer_vendor", coString);
    def->label = L("Printer vendor");
    def->tooltip = L("Name of the printer vendor.");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("printer_variant", coString);
    def->label = L("Printer variant");
    def->tooltip = L("Name of the printer variant. For example, the printer variants may be differentiated by a nozzle diameter.");
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("print_settings_id", coString);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString(""));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("print_settings_modified", coBool);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionBool(false));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("printer_settings_id", coString);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString(""));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("printer_settings_modified", coBool);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionBool(false));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("physical_printer_settings_id", coString);
    def->mode = comNone | comPrusa; // note: hidden setting
    def->set_default_value(new ConfigOptionString(""));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("filament_default_pa", coFloats);
    def->label = L("default");
    def->category = OptionCategory::filament;
    def->tooltip = L("Default linear/pressure advance. This is only activated for marlin, klipper and reprap firmware. For the other ones, you have to add it yourself in the gcode."
        "Note that the meaning of this value differ for each firmware, so if you set it, it means that this profile may be incompatible with some."
        "\nSet 0 to deactivate.");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0. });

    def = this->add("filament_bridge_pa", coFloatsOrPercents);
    def->label = L("bridge");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for bridge sections. Can be a % over default pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_bridge_internal_pa", coFloatsOrPercents);
    def->label = L("internal bridge");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for internal bridge sections. Can be a % over default pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_brim_pa", coFloatsOrPercents);
    def->label = L("brim");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for brim. Can be a % over support pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_support_material_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_external_perimeter_pa", coFloatsOrPercents);
    def->label = L("external perimeter");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for external perimeter. Can be a % over support pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_support_material_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_first_layer_pa", coFloatsOrPercents);
    def->label = L("first layer");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for first layer sections. If %, it's a % over the current feature");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "depends";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_first_layer_pa_over_raft", coFloatsOrPercents);
    def->label = L("over raft");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for first layer sections over raft . If %, it's a % over the current feature");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "depends";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_gap_fill_pa", coFloatsOrPercents);
    def->label = L("gap fill");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for gap fill sections. Can be a % over perimeter pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_perimeter_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_infill_pa", coFloatsOrPercents);
    def->label = L("infill");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for infill sections. Can be a % over solid infill pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_solid_infill_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_ironing_pa", coFloatsOrPercents);
    def->label = L("ironing");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for ironing sections. Can be a % over top solid infill pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_top_solid_infill_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_overhangs_pa", coFloatsOrPercents);
    def->label = L("overhangs");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for overhang sections. Can be a % over bridge pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_bridge_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_perimeter_pa", coFloatsOrPercents);
    def->label = L("perimeters");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for perimeter sections. Can be a % over default pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_solid_infill_pa", coFloatsOrPercents);
    def->label = L("solid infill");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for solid infill sections. Can be a % over default pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_support_material_pa", coFloatsOrPercents);
    def->label = L("support");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for support sections. Can be a % over default pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_support_material_interface_pa", coFloatsOrPercents);
    def->label = L("support interface");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for support interface sections. Can be a % over support pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_support_material_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_thin_walls_pa", coFloatsOrPercents);
    def->label = L("thin walls");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for thin wall sections. Can be a % over external perimeter pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_external_perimeter_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_top_solid_infill_pa", coFloatsOrPercents);
    def->label = L("top solid infill");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for top solid infill sections. Can be a % over solid infill pa");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_solid_infill_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{100,true} });

    def = this->add("filament_travel_pa", coFloatsOrPercents);
    def->label = L("travel");
    def->category = OptionCategory::filament;
    def->tooltip = L("Pressure advance for travel sections, may help retraction and unretraction."
            " Can be a % over default pa."
            "\nSet -1 to let the previous pa continue in the travel.");
    def->mode = comExpert | comSuSi;
    def->ratio_over = "filament_default_pa";
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{-1,false} });

    def = this->add("raft_contact_distance", coFloat);
    def->label = L("Raft contact Z distance");
    def->category = OptionCategory::support;
    def->tooltip = L("The vertical distance between object and raft. Ignored for soluble interface. It uses the same type as the support z-offset type.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.1));

    def = this->add("raft_expansion", coFloat);
    def->label = L("Raft expansion");
    def->category = OptionCategory::support;
    def->tooltip = L("Expansion of the raft in XY plane for better stability.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.5));

    def = this->add("raft_first_layer_density", coPercent);
    def->label = L("First layer density");
    def->full_label = L("Raft first layer density");
    def->category = OptionCategory::support;
    def->tooltip = L("Density of the first raft or support layer.");
    def->sidetext = L("%");
    def->min = 10;
    def->max = 100;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionPercent(90));

    def = this->add("raft_first_layer_expansion", coFloat);
    def->label = L("First layer expansion");
    def->full_label = L("Raft first layer expansion");
    def->category = OptionCategory::support;
    def->tooltip = L("Expansion of the first raft or support layer to improve adhesion to print bed.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(3.));

    def = this->add("raft_layers", coInt);
    def->label = L("Raft layers");
    def->category = OptionCategory::support;
    def->tooltip = L("The object will be raised by this number of layers, and support material "
        "will be generated under it.");
    def->sidetext = L("layers");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("raft_layer_height", coFloatOrPercent);
    def->label = L("Raft layer height");
    def->category = OptionCategory::support;
    def->tooltip = L("Maximum layer height for the raft, after the first layer that uses the first layer height, and before the interface layers."
        "\nCan be a % of the nozzle diameter"
        "\nIf set to 0, the support layer height will be used.");
    def->sidetext = L("mm");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("raft_interface_layer_height", coFloatOrPercent);
    def->label = L("Raft interface layer height");
    def->category = OptionCategory::support;
    def->tooltip = L("Maximum layer height for the raft interface."
        "\nCan be a % of the nozzle diameter"
        "\nIf set to 0, the support layer height will be used.");
    def->sidetext = L("mm");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("region_gcode", coString);
    def->label = L("Per region G-code");
    def->category = OptionCategory::output;
    def->tooltip = L("This code is inserted when a region is starting to print something (infill, perimeter, ironing)."
                     " It's main advantage is when you use it as a object modifer(right click on a model to add it there)"
                     "\nSpecial variables: 'layer_num','layer_z'");
    def->multiline = true;
    def->full_width = true;
    def->height = 10;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("resolution", coFloat);
    def->label = L("Slice resolution");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Minimum detail resolution, used to simplify the input file for speeding up "
        "the slicing job and reducing memory usage. High-resolution models often carry "
        "more details than printers can render. Set zero to disable any simplification "
        "and use full resolution from input. "
        "\nNote: Slic3r has an internal working resolution of 0.0001mm."
        "\nInfill & Thin areas are simplified up to 0.0125mm.");
    def->sidetext = L("mm");
    def->min = 0;
    def->precision = 8;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.0125));

    //note: replaced by gcode_min_resolution in printer's profile
    //def = this->add("gcode_resolution", coFloat);
    //def->label = L("G-code resolution");
    //def->tooltip = L("Maximum deviation of exported G-code paths from their full resolution counterparts. "
    //    "Very high resolution G-code requires huge amount of RAM to slice and preview, "
    //    "also a 3D printer may stutter not being able to process a high resolution G-code in a timely manner. "
    //    "On the other hand, a low resolution G-code will produce a low poly effect and because "
    //    "the G-code reduction is performed at each layer independently, visible artifacts may be produced.");
    //def->sidetext = L("mm");
    //def->min = 0;
    //def->mode = comExpert | comPrusa;
    //def->set_default_value(new ConfigOptionFloat(0.0));

    def = this->add("gcode_command_buffer", coInt);
    def->label = L("Command buffer");
    def->category = OptionCategory::speed;
    def->tooltip = L("Buffer the firmware has for gcode comamnds. Allow to have some burst of command with a rate over 'max_gcode_per_second' for a few instant.");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionInt(10));

    def = this->add("gcode_min_length", coFloatOrPercent);
    def->label = L("Minimum extrusion length");
    def->category = OptionCategory::speed;
    def->tooltip = L("When outputting gcode, this setting ensure that there is almost no commands more than this value apart."
        " Be sure to also use max_gcode_per_second instead, as it's much better when you have very different speeds for features"
        " (Too many too small commands may overload the firmware / connection)."
        "\nSet zero to disable.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->precision = 6;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0/*0.02*/, false));
    def->aliases = {"min_length"};

    def = this->add("gcode_min_resolution", coFloatOrPercent);
    def->label = L("minimum resolution");
    def->category = OptionCategory::speed;
    def->tooltip = L("Maximum deviation of exported G-code paths from their full resolution counterparts"
        " when some commands are culled by 'gcode_min_length' or 'max_gcode_per_second' to not overload the firmware."
        "\nCan be a % of perimeter width.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->precision = 6;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("resolution_internal", coFloat);
    def->label = L("Internal resolution");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Minimum detail resolution, used for internal structures (gapfill and some infill patterns)."
            "\nDon't put a too-small value (0.05mm is way too low for many printers A3dp machines .1 rec), as it may create too many very small segments that may be difficult to display and print.");
    def->sidetext = L("mm");
    def->min = 0.001;
    def->precision = 8;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0.1));

    def = this->add("retract_before_travel", coFloats);
    def->label = L("Minimum travel after retraction");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Retraction is not triggered when travel moves are shorter than this length.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->min = 0;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 2. });

    def = this->add("retract_lift_before_travel", coFloats);
    def->label = L("Minimum travel after z lift");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Z lift is not triggered when travel moves are shorter than this length.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->min = 0;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 2. });

    def = this->add("retract_before_wipe", coPercents);
    def->label = L("Retract amount before wipe");
    def->category = OptionCategory::extruders;
    def->tooltip = L("With bowden extruders, it may be wise to do some amount of quick retract "
                   "before doing the wipe movement.");
    def->sidetext = L("%");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents { 0. });

    def = this->add("retract_layer_change", coBools);
    def->label = L("Retract on layer change");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This flag enforces a retraction whenever a Z move is done (before it).");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools { false });

    def = this->add("retract_length", coFloats);
    def->label = L("Retraction length");
    def->full_label = L("Retraction Length");
    def->category = OptionCategory::extruders;
    def->tooltip = L("When retraction is triggered, filament is pulled back by the specified amount "
                   "(the length is measured on raw filament, before it enters the extruder).");
    def->sidetext = L("mm (zero to disable)");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 2. });

    def = this->add("print_retract_length", coFloat);
    def->label = L("Retraction length");
    def->category = OptionCategory::filament;
    def->tooltip = L("Override the retract_length setting from the printer config. Used for calibration. Set negative to disable");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat( -1.f));

    def = this->add("retract_length_toolchange", coFloats);
    def->label = L("Length");
    def->full_label = L("Retraction Length (Toolchange)");
    def->tooltip = L("When retraction is triggered before changing tool, filament is pulled back "
                   "by the specified amount (the length is measured on raw filament, before it enters "
                   "the extruder)."
                    "\nNote: This value will be unretracted when this extruder will load the next time.");
    def->sidetext = L("mm (zero to disable)");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 10. });

    def = this->add("travel_slope", coFloats);
    def->label = L("Ramping slope angle");
    def->tooltip = L("Minimum slope of the ramp in the initial phase of the travel."
                    " If the travel isn't long enough, the angle will be increased."
                    "\n90° means a direct lift, like if there was no ramp."
                    "\n0° means that the lift will always be hit at the end of the travel.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 90;
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{0.0});

    def = this->add("travel_ramping_lift", coBools);
    def->label = L("Use ramping lift");
    def->tooltip = L("Generates a ramping lift instead of lifting the extruder directly upwards. "
                     "The travel is split into two phases: the ramp and the standard horizontal travel. "
                     "This option helps reduce stringing."
                     "\nAlso works for the z move when a layer change occurs.");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ false });

    // why not reuse rretract_lift ? because it's a max? My current impl enforced the lift, so it's okay for me to remove it.
    // def = this->add("travel_max_lift", coFloats);
    // def->label = L("Maximum ramping lift");
    // def->tooltip = L("Maximum lift height of the ramping lift. It may not be reached if the next position "
                     // "is close to the old one.");
    // def->sidetext = L("mm");
    // def->min = 0;
    // def->max_literal = {1000, false};
    // def->mode = comAdvancedE | comPrusa;
    // def->is_vector_extruder = true;
    // def->set_default_value(new ConfigOptionFloats{0.0});

    def = this->add("travel_lift_before_obstacle", coBools);
    def->label = L("Steeper ramp before obstacles");
    def->tooltip = L("If enabled, PrusaSlicer detects obstacles along the travel path and makes the slope steeper "
                     "in case an obstacle might be hit during the initial phase of the travel.");
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{false});

    def = this->add("retract_lift", coFloats);
    def->label = L("Lift height");
    def->category = OptionCategory::extruders;
    def->tooltip = L("If you set this to a positive value, the extruder is quickly raised every time a retraction "
                   "is triggered. When using multiple extruders, only the setting for the first extruder "
                   "will be considered.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max_literal = {1000, false};
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{0.});

    def = this->add("retract_lift_above", coFloats);
    def->label = L("Above Z");
    def->full_label = L("Only lift Z above");
    def->category = OptionCategory::extruders;
    def->tooltip = L("If you set this to a positive value, Z lift will only take place above the specified "
                   "absolute Z. You can tune this setting for skipping lift on the first layers.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });



    def = this->add("retract_lift_below", coFloats);
    def->label = L("Below Z");
    def->full_label = L("Only lift Z below");
    def->category = OptionCategory::extruders;
    def->tooltip = L("If you set this to a positive value, Z lift will only take place below "
                   "the specified absolute Z. You can tune this setting for limiting lift "
                   "to the first layers.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("retract_lift_first_layer", coBools);
    def->label = L("Enforce on first layer");
    def->full_label = L("Enforce lift on first layer");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Select this option to enforce z-lift on the first layer."
        "\nUseful to still use the lift on the first layer even if the 'Only lift Z below' (retract_lift_above) is higher than 0.");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ false });

    def = this->add("retract_lift_top", coStrings);
    def->label = L("On surfaces");
    def->full_label = L("Lift only on");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Select this option to not use/enforce the z-lift on a top surface.");
    def->gui_type = ConfigOptionDef::GUIType::f_enum_open;
    def->gui_flags = "show_value";
    def->set_enum_values(ConfigOptionDef::GUIType::select_open,
        { "All surfaces", "Not on top", "Only on top"});
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings{ "All surfaces" });


    def = this->add("retract_restart_extra", coFloats);
    def->label = L("Deretraction extra length");
    def->tooltip = L("When the retraction is compensated after the travel move, the extruder will push "
                   "this additional amount of filament. This setting is rarely needed.");
    def->sidetext = L("mm");
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("retract_restart_extra_toolchange", coFloats);
    def->label = L("Extra length on restart");
    def->full_label = L("Extrat length on toolchange restart");
    def->tooltip = L("When the retraction is compensated after changing tool, the extruder will push "
                    "this additional amount of filament"
                    " (but not on the first extruder after start, as it should already be loaded).");
    def->sidetext = L("mm");
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("retract_speed", coFloats);
    def->label = L("Retraction Speed");
    def->full_label = L("Retraction Speed");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The speed for retractions (this only applies to the extruder motor).");
    def->sidetext = L("mm/s");
    def->mode = comAdvancedE | comPrusa;
    def->min = 0.001;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 40. });

    def = this->add("deretract_speed", coFloats);
    def->label = L("Deretraction Speed");
    def->full_label = L("Deretraction Speed");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The speed for loading of a filament into extruder after retraction "
                   "(this only applies to the extruder motor). If left as zero, the retraction speed is used.");
    def->sidetext = L("mm/s");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats { 0. });

    def = this->add("seam_position", coEnum);
    def->label = L("Seam position");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Position of perimeters' starting points."
                    "\nCost-based option let you choose the angle and travel cost. A high angle cost will place the seam where it can be hidden by a corner"
                    ", the travel cost place the seam near the last position (often at the end of the previous infill). Default is 60 % and 100 %."
                    " There is also the visibility and the overhang cost, but they are static."
                    "\n Scattered: seam is placed at a random position on external perimeters"
                    "\n Random: seam is placed at a random position for all perimeters"
                    "\n Aligned: seams are grouped in the best place possible (minimum 6 layers per group)"
                    "\n Contiguous: seam is placed over a seam from the previous layer (useful with enforcers)"
                    "\n Rear: seam is placed at the far side (highest Y coordinates)");
    def->set_enum<SeamPosition>({
        { "cost",       L("Cost-based") },
        { "random",     L("Scattered") },
        { "allrandom",  L("Random") },
        { "aligned",    L("Aligned") },
        { "contiguous", L("Contiguous") },

        { "rear",       L("Rear") }
    });
    def->mode = comSimpleAE | comPrusa | comSuSi;
    def->set_default_value(new ConfigOptionEnum<SeamPosition>(spCost));

    def = this->add("seam_angle_cost", coPercent);
    def->label = L("Angle cost");
    def->full_label = L("Seam angle cost");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Cost of placing the seam at a bad angle. The worst angle (max penalty) is when it's flat."
        "\n100% is the default penalty");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(60));

    def = this->add("seam_gap", coFloatsOrPercents);
    def->label = L("Seam gap");
    def->category = OptionCategory::extruders;
    def->tooltip = L("To avoid visible seam, the extrusion can be stoppped a bit before the end of the loop."
        "\nCan be a mm or a % of the current extruder diameter.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max_literal = { 5, false };
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{15,true} });

    def = this->add("seam_gap_external", coFloatsOrPercents);
    def->label = L("Seam gap for external perimeters");
    def->category = OptionCategory::extruders;
    def->tooltip = L("To avoid visible seam, the extrusion can be stoppped a bit before the end of the loop."
        "\n this setting is enforced only for external perimeter. It overrides 'seam_gap' if different than 0"
        "\nCan be a mm or a % of the current seam gap.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max_literal = { 5, false };
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloatsOrPercents{ FloatOrPercent{0,false} });

    def = this->add("seam_notch_all", coFloatOrPercent);
    def->label = L("for everything");
    def->full_label = L("Seam notch");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("It's sometimes very problematic to have a little buldge from the seam."
        " This setting move the seam inside the part, in a little cavity (for every seams in external perimeters, unless it's in an overhang)."
        "\nThe size of the cavity is in mm or a % of the external perimeter width. It's overriden by the two other 'seam notch' setting when applicable."
        "\nSet zero to disable.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max_literal = { 5, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("seam_notch_angle", coFloat);
    def->label = L("max angle");
    def->full_label = L("Seam notch maximum angle");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("If the (external) angle at the seam is higher than this value, then no notch will be set. If the angle is too high, there isn't enough room for the notch."
                    "\nCan't be lower than 180° or it filters everything. At 360, it allows everything.");
    def->sidetext = L("°");
    def->min = 180;
    def->max = 360;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(250));

    def = this->add("seam_notch_inner", coFloatOrPercent);
    def->label = L("for round holes");
    def->full_label = L("Seam notch for round holes");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("In convex holes (circular/oval), it's sometimes very problematic to have a little buldge from the seam."
        " This setting move the seam inside the part, in a little cavity (for all external perimeters in convex holes, unless it's in an overhang)."
        "\nThe size of the cavity is in mm or a % of the external perimeter width"
        "\nSet zero to disable.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max = 50;
    def->max_literal = { 5, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("seam_notch_outer", coFloatOrPercent);
    def->label = L("for round perimeters");
    def->full_label = L("Seam notch for round perimeters");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("In convex perimeters (circular/oval), it's sometimes very problematic to have a little buldge from the seam."
        " This setting move the seam inside the part, in a little cavity (for all external perimeters if the path is convex, unless it's in an overhang)."
        "\nThe size of the cavity is in mm or a % of the external perimeter width"
        "\nSet zero to disable.");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->max = 50;
    def->max_literal = { 5, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("seam_travel_cost", coPercent);
    def->label = L("Travel cost");
    def->full_label = L("Seam travel cost");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Cost of moving the extruder. The highest penalty is when the point is the furthest from the position of the extruder before extruding the external perimeter");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("seam_visibility", coBool);
    def->label = L("use visibility check");
    def->full_label = L("Seam visibility check");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Check and penalize seams that are the most visible. launch rays to check from how many direction a point is visible."
        "\nThis is a compute-intensive option.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("staggered_inner_seams", coBool);
    def->label = L("Staggered inner seams");
    // TRN PrintSettings: "Staggered inner seams"
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This option causes the inner seams to be shifted backwards based on their depth, forming a zigzag pattern.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

#if 0
    def = this->add("seam_preferred_direction", coFloat);
//    def->gui_type = ConfigOptionDef::GUIType::slider;
    def->label = L("Direction");
    def->sidetext = L("°");
    def->full_label = L("Preferred direction of the seam");
    def->tooltip = L("Seam preferred direction");
    def->min = 0;
    def->max = 360;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("seam_preferred_direction_jitter", coFloat);
//    def->gui_type = ConfigOptionDef::GUIType::slider;
    def->label = L("Jitter");
    def->sidetext = L("°");
    def->full_label = L("Seam preferred direction jitter");
    def->tooltip = L("Preferred direction of the seam - jitter");
    def->min = 0;
    def->max = 360;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(30));
#endif
    def = this->add("skirt_brim", coInt);
    def->label = L("Brim");
    def->full_label = L("Skirt brim");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Extra skirt lines on the first layer.");
    def->sidetext = L("lines");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("skirt_distance", coFloat);
    def->label = L("Distance from object");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Distance between skirt and object(s) ; or from the brim if using draft shield or you set 'skirt_distance_from_brim'.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(6));

    def = this->add("skirt_distance_from_brim", coBool);
    def->label = L("from brim");
    def->full_label = L("Skirt distance from brim");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("The distance is computed from the brim and not from the objects");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("skirt_height", coInt);
    def->label = L("Skirt height");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Height of skirt expressed in layers. Set this to a tall value to use skirt "
        "as a shield against drafts.");
    def->sidetext = L("layers");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("skirt_extrusion_width", coFloatOrPercent);
    def->label = L("Skirt");
    def->full_label = L("Skirt width");
    def->category = OptionCategory::width;
    def->tooltip = L("Horizontal width of the skirt that will be printed around each object."
        " If left as zero, first layer extrusion width will be used if set and the skirt is only 1 layer height"
        ", or perimeter extrusion width will be used (using the computed value if not set).");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(130, true));

    def = this->add("skirts", coInt);
    def->label = L("Loops (minimum)");
    def->full_label = L("Skirt Loops");
    def->category = OptionCategory::skirtBrim;
    def->tooltip = L("Number of loops for the skirt. If the Minimum Extrusion Length option is set, "
                   "the number of loops might be greater than the one configured here. Set zero "
                   "to disable skirt completely.");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("slicing_mode", coEnum);
    def->label = L("Slicing Mode");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Use \"Even-odd\" for 3DLabPrint airplane models. Use \"Close holes\" to close all holes in the model.");
    def->set_enum<SlicingMode>({
        { "regular",        L("Regular") },
        { "even_odd",       L("Even-odd") },
        { "close_holes",    L("Close holes") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<SlicingMode>(SlicingMode::Regular));

    def = this->add("slowdown_below_layer_time", coFloats);
    def->label = L("Slow down if layer print time is below");
    def->category = OptionCategory::cooling;
    def->tooltip = L("If layer print time is estimated below this number of seconds, print moves "
        "speed will be scaled down to extend duration to this value, if possible."
        "\nSet zero to disable.");
    def->sidetext = L("approximate seconds");
    def->min = 0;
    def->max = 1000;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 5 });

    def = this->add("small_perimeter_speed", coFloatOrPercent);
    def->label = L("Speed");
    def->full_label = L("Small perimeters speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("This separate setting will affect the speed of perimeters having radius <= 6.5mm (usually holes)."
                   "\nIf expressed as percentage (for example: 80%) it will be calculated on the Internal Perimeters speed setting above."
                   "\nSet zero to disable.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "perimeter_speed";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("small_perimeter_min_length", coFloatOrPercent);
    def->label = L("Min length");
    def->full_label = L("Min small perimeters length");
    def->category = OptionCategory::speed;
    def->tooltip = L("This sets the threshold for small perimeter length. Every loop with a length lower than this will be printed at small perimeter speed"
        "\nCan be a mm value or a % of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 100, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(6, false));


    def = this->add("small_perimeter_max_length", coFloatOrPercent);
    def->label = L("Max length");
    def->full_label = L("Max small perimeters length");
    def->category = OptionCategory::speed;
    def->tooltip = L("This sets the end of the threshold for small perimeter length."
        " Every perimeter loop lower than this will see their speed reduced a bit, from their normal speed at this length down to small perimeter speed."
        "\nCan be a mm or a % of the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max_literal = { 500, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(20, false));

    def = this->add("solid_infill_below_area", coFloat);
    def->label = L("Solid infill threshold area");
    def->category = OptionCategory::infill;
    def->tooltip = L("Force solid infill for regions having a smaller area than the specified threshold.");
    def->sidetext = L("mm²");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(4));

    def = this->add("solid_infill_below_layer_area", coFloat);
    def->label = L("Solid infill layer threshold area");
    def->category = OptionCategory::infill;
    def->tooltip = L("Force solid infill for the whole layer when the combined area of all objects that are printed at the same layer is smaller than this value.");
    def->sidetext = L("mm²");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("solid_infill_below_width", coFloatOrPercent);
    def->label = L("Solid infill threshold width");
    def->category = OptionCategory::infill;
    def->tooltip = L("Force solid infill for parts of regions having a smaller width than the specified threshold."
                    "\nCan be a % of the current solid infill spacing."
                    "\nSet 0 to disable");
    def->sidetext = L("mm or %");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("solid_infill_overlap", coPercent);
    def->label = L("Solid infill overlap");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allows you to reduce the overlap between the lines of the solid fill, to reduce the % filled if you see overextrusion signs on solid areas."
        " Note that you should be sure that your flow (filament extrusion multiplier) is well calibrated and your filament max overlap is set before thinking to modify this."
        "\nNote: top surfaces are still extruded with 100% overlap to prevent gaps.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("solid_infill_extruder", coInt);
    def->label = L("Solid infill extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing solid infill.");
    def->min = 1;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("solid_infill_every_layers", coInt);
    def->label = L("Solid infill every");
    def->category = OptionCategory::infill;
    def->tooltip = L("This feature allows you to force a solid layer every given number of layers. "
                   "Zero to disable. You can set this to any value (for example 9999); "
                   "Slic3r will automatically choose the maximum possible number of layers "
                   "to combine according to nozzle diameter and layer height.");
    def->sidetext = L("layers");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("solid_infill_extrusion_width", coFloatOrPercent);
    def->label = L("Solid infill");
    def->full_label = L("Solid infill width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for infill for solid surfaces. "
        "If left as zero, default extrusion width will be used if set, otherwise 1.125 x nozzle diameter will be used. "
        "If expressed as percentage (for example 110%) it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("solid_infill_extrusion_change_odd_layers", coFloatOrPercent);
    def->label = L("Infill");
    def->full_label = L("Solid infill spacing change on even layers");
    def->category = OptionCategory::width;
    def->tooltip = L("Change width on every even layer (and not on odd layers like the first one) for better overlap with adjacent layers and getting stringer shells. "
        "Try values about +/- 0.1 with different sign."
        "\nThis could be combined with extra permeters on even layers."
        "\nWorks as absolute spacing or a % of the spacing."
        "\nset 0 to disable");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(false, 0));

    def = this->add("solid_infill_extrusion_spacing", coFloatOrPercent);
    def->label = L("Solid spacing");
    def->full_label = L("Solid infill spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like Solid infill width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("solid_infill_fan_speed", coInts);
    def->label = L("Solid Infill fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all Solid Infill moves"
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer and increased by low layer time.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("solid_infill_speed", coFloatOrPercent);
    def->label = L("Solid");
    def->full_label = L("Solid infill speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing solid regions (top/bottom/internal horizontal shells). "
        "\nThis can be expressed as a percentage (for example: 80%) over the Default speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "default_speed";
    def->aliases = { "solid_infill_feed_rate" };
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(30, true));

#if 0
    def = this->add("solid_layers", coInt);
    def->label = L("Solid layers");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Number of solid layers to generate on top and bottom surfaces.");
    def->shortcut.push_back("top_solid_layers");
    def->shortcut.push_back("bottom_solid_layers");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;

    def = this->add("solid_min_thickness", coFloat);
    def->label = L("Minimum thickness of a top / bottom shell");
    def->tooltip = L("Minimum thickness of a top / bottom shell");
    def->shortcut.push_back("top_solid_min_thickness");
    def->shortcut.push_back("bottom_solid_min_thickness");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
#endif

    def = this->add("spiral_vase", coBool);
    def->label = L("Spiral vase");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("This feature will raise Z gradually while printing a single-walled object "
                   "in order to remove any visible seam. This option requires "
                   "no infill, no top solid layers and no support material. You can still set "
                   "any number of bottom solid layers as well as skirt/brim loops."
                   " After the bottom solid layers, the number of perimeters is enforce to 1."
                   "It won't work when printing more than one single object.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("standby_temperature_delta", coInt);
    def->label = L("Temperature variation");
    // TRN PrintSettings : "Ooze prevention" > "Temperature variation"
    def->tooltip = L("Temperature difference to be applied when an extruder is not active. "
                     "The value is not used when 'idle_temperature' in filament settings "
                     "is defined.");
    def->sidetext = "∆°C";
    def->min = -max_temp;
    def->max = max_temp;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionInt(-5));

    def = this->add("autoemit_temperature_commands", coBool);
    def->label = L("Emit temperature commands automatically");
    def->tooltip = L("When enabled, Slic3r will check whether your Custom Start G-Code contains M104 or M190. "
                     "If so, the temperatures will not be emitted automatically so you're free to customize "
                     "the order of heating commands and other custom actions. Note that you can use "
                     "placeholder variables for all Slic3r settings, so you can put "
                     "a \"M109 S[first_layer_temperature]\" command wherever you want.\n"
                     "If your Custom Start G-Code does NOT contain M104 or M190, "
                     "Slic3r will execute the Start G-Code after bed reached its target temperature "
                     "and extruder just started heating.\n\n"
                     "When disabled, Slic3r will NOT emit commands to heat up extruder and bed, "
                     "leaving both to Custom Start G-Code.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("start_gcode", coString);
    def->label = L("Start G-code");
    def->tooltip = L("This start procedure is inserted at the beginning, possibly prepended by "
                     "temperature-changing commands. See 'autoemit_temperature_commands'.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString("G28 ; home all axes\nG1 Z5 F5000 ; lift nozzle\n"));

    def = this->add("start_gcode_manual", coBool);
    def->label = L("Only custom Start G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("Ensure that the slicer won't add heating, fan, extruder... commands before or just after your start-gcode."
                    "\nIf set to true, you have to write a good and complete start_gcode, as no checks are made anymore."
                    "\nExemple:\nG21 ; set units to millimeters\nG90 ; use absolute coordinates\n{if use_relative_e_distances}M83{else}M82{endif}\nG92 E0 ; reset extrusion distance");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("start_filament_gcode", coStrings);
    def->label = L("Start G-code");
    def->full_label = L("Filament start G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This start procedure is inserted at the beginning, after any printer start gcode (and "
                   "after any toolchange to this filament in case of multi-material printers). "
                   "This is used to override settings for a specific filament. If Slic3r detects "
                   "M104, M109, M140 or M190 in your custom codes, such commands will "
                   "not be prepended automatically so you're free to customize the order "
                   "of heating commands and other custom actions. Note that you can use placeholder variables "
                   "for all Slic3r settings, so you can put a \"M109 S{first_layer_temperature}\" command "
                   "wherever you want. If you have multiple extruders, the gcode is processed "
                   "in extruder order.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings { "; Filament gcode\n" });

    def = this->add("color_change_gcode", coString);
    def->label = L("Color change G-code");
    def->tooltip = L("This G-code will be used as a code for the color change"
                     " If empty, the default color change print command for the selected G-code flavor will be used (if any).");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("pause_print_gcode", coString);
    def->label = L("Pause Print G-code");
    def->tooltip = L("This G-code will be used as a code for the pause print."
                    " If empty, the default pause print command for the selected G-code flavor will be used (if any).");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("template_custom_gcode", coString);
    def->label = L("Custom G-code");
    def->tooltip = L("This G-code will be used as a custom code");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("single_extruder_multi_material", coBool);
    def->label = L("Single Extruder Multi Material");
    def->category = OptionCategory::mmsetup;
    def->tooltip = L("The printer multiplexes filaments into a single hot end.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("single_extruder_multi_material_priming", coBool);
    def->label = L("Prime all printing extruders");
    def->category = OptionCategory::mmsetup;
    def->tooltip = L("If enabled, all printing extruders will be primed at the front edge of the print bed at the start of the print.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("solid_infill_acceleration", coFloatOrPercent);
    def->label = L("Solid ");
    def->full_label = L("Solid acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for solid infill. "
                "\nCan be a % of the default acceleration"
                "\nSet zero or 100% to use default acceleration for solid infill.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "default_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("solid_over_perimeters", coInt);
    def->label = L("No solid infill over");
    def->full_label = L("No solid infill over perimeters");
    def->sidetext = L("perimeters");
    def->sidetext_width = 20;
    def->category = OptionCategory::perimeter;
    def->tooltip = L("In sloping areas, when you have a number of top / bottom solid layers and few perimeters, "
        " it may be necessary to put some solid infill above/below the perimeters to fulfill the top/bottom layers criteria."
        "\nBy setting this to something higher than 0, you can control this behaviour, which might be desirable if "
        "\nundesirable solid infill is being generated on slopes."
        "\nThe number set here indicates the number of layers between the inside of the part and the air"
        " at and beyond which solid infill should no longer be added above/below. If this setting is equal or higher than "
        " the top/bottom solid layer count, it won't do anything. If this setting is set to 1, it will evict "
        " all solid fill above/below perimeters. "
        "\nSet zero to disable.");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionInt(2));

    def = this->add("support_material", coBool);
    def->label = L("Generate support material");
    def->category = OptionCategory::support;
    def->tooltip = L("Enable support material generation.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("support_material_acceleration", coFloatOrPercent);
    def->label = L("Default");
    def->full_label = L("Support acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for support material. "
                "\nCan be a % of the default acceleration"
                "\nSet zero to use default acceleration for support material.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "default_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("support_material_auto", coBool);
    def->label = L("Auto generated supports");
    def->category = OptionCategory::support;
    def->tooltip = L("If checked, supports will be generated automatically based on the overhang threshold value."\
                     " If unchecked, supports will be generated inside the \"Support Enforcer\" volumes only.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("support_material_interface_acceleration", coFloatOrPercent);
    def->label = L("Interface");
    def->full_label = L("Support interface acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for support material interfaces. "
                "\nCan be a % of the support material acceleration"
                "\nSet zero to use support acceleration for support material interfaces.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "support_material_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("support_material_xy_spacing", coFloatOrPercent);
    def->label = L("XY separation between an object and its support");
    def->category = OptionCategory::support;
    def->tooltip = L("XY separation between an object and its support. If expressed as percentage "
                   "(for example 50%), it will be calculated over external perimeter width.");
    def->sidetext = L("mm or %");
    def->ratio_over = "external_perimeter_extrusion_width";
    def->min = 0;
    def->max_literal = { 10, false};
    def->mode = comAdvancedE | comPrusa;
    // Default is half the external perimeter width.
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("support_material_angle", coFloat);
    def->label = L("Pattern angle");
    def->full_label = L("Support pattern angle");
    def->category = OptionCategory::support;
    def->tooltip = L("Use this setting to rotate the support material pattern on the horizontal plane.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 359;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("support_material_angle_height", coFloat);
    def->label = L("Pattern angle swap height");
    def->category = OptionCategory::support;
    def->tooltip = L("Use this setting to rotate the support material pattern by 90° at this height (in mm). Set 0 to disable.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("support_material_buildplate_only", coBool);
    def->label = L("Support on build plate only");
    def->category = OptionCategory::support;
    def->tooltip = L("Only create support if it lies on a build plate. Don't create support on a print.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("support_material_contact_distance_type", coEnum);
    def->label = L("Type");
    def->full_label = L("Support contact distance type");
    def->category = OptionCategory::support;
    def->tooltip = L("How to compute the vertical z-distance.\n"
        "From filament: it uses the nearest bit of the filament. When a bridge is extruded, it goes below the current plane.\n"
        "From plane: it uses the plane-z. Same as 'from filament' if no 'bridge' is extruded.\n"
        "None: No z-offset. Useful for Soluble supports.\n");
    def->set_enum<SupportZDistanceType>({
        { "filament", L("From filament") },
        { "plane",    L("From plane") },
        { "none",     L("None (soluble)") }
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<SupportZDistanceType>(zdFilament));

    def = this->add("support_material_contact_distance", coFloatOrPercent);
    def->label = L("Top");
    def->full_label = L("Contact distance on top of supports");
    def->category = OptionCategory::support;
    def->tooltip = L("The vertical distance between support material interface and the object"
        "(when the object is printed on top of the support). "
        "Setting this to 0 will also prevent Slic3r from using bridge flow and speed "
        "for the first object layer. Can be a % of the nozzle diameter.");
    def->ratio_over = "nozzle_diameter";
    def->sidetext = L("mm");
    def->min = 0;
    def->max_literal = { 20, true };
    def->mode = comAdvancedE | comPrusa;
    def->aliases = { "support_material_contact_distance_top" }; // Sli3r, PS
    def->set_default_value(new ConfigOptionFloatOrPercent(0.2, false));

    def = this->add("support_material_bottom_contact_distance", coFloatOrPercent);
    def->label = L("Bottom");
    def->full_label = L("Contact distance under the bottom of supports");
    def->category = OptionCategory::support;
    def->tooltip = L("The vertical distance between object and support material interface"
        "(when the support is printed on top of the object). Can be a % of the nozzle diameter."
        "\nIf set to zero, support_material_contact_distance will be used for both top and bottom contact Z distances.");
    def->ratio_over = "nozzle_diameter";
    def->sidetext = L("mm");
    def->set_enum_values(ConfigOptionDef::GUIType::f_enum_open, {
    //TRN Print Settings: "Bottom contact Z distance". Have to be as short as possible
        { "0",      L("Same as top") },
        { "0.1",    "0.1" },
        { "0.2",    "0.2" },
        { "50%",    "50%" },
    });
    def->min = 0;
    def->max_literal = { 20, true };
    def->mode = comAdvancedE | comPrusa;
    def->aliases = { "support_material_contact_distance_bottom" }; //since PS 2.4
    def->set_default_value(new ConfigOptionFloatOrPercent(0.2,false));

    def = this->add("support_material_enforce_layers", coInt);
    def->label = L("Enforce support for the first");
    def->category = OptionCategory::support;
    def->tooltip = L("Generate support material for the specified number of layers counting from bottom, "
                   "regardless of whether normal support material is enabled or not and regardless "
                   "of any angle threshold. This is useful for getting more adhesion of objects "
                   "having a very thin or poor footprint on the build plate.");
    def->sidetext = L("layers");
    def->full_label = L("Enforce support for the first n layers");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("support_material_extruder", coInt);
    def->label = L("Support material extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing support material "
                   "(1+, 0 to use the current extruder to minimize tool changes).");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("support_material_extrusion_width", coFloatOrPercent);
    def->label = L("Support material");
    def->full_label = L("Support material width");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for support material. "
        "If left as zero, default extrusion width will be used if set, otherwise nozzle diameter will be used. "
        "If expressed as percentage (for example 110%) it will be computed over nozzle diameter.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("support_material_fan_speed", coInts);
    def->label = L("Support Material fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all support moves"
        "\nSet to 0 to stop the fan."
        "\nIf disabled, default fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("support_material_interface_angle", coFloat);
    def->label = L("Pattern angle");
    def->full_label = L("Support interface pattern angle");
    def->category = OptionCategory::support;
    def->tooltip = L("Use this setting to rotate the support material pattern on the horizontal plane.\n0 to use the support_material_angle.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 360;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(90));

    def = this->add("support_material_interface_angle_increment", coFloat);
    def->label = L("Support interface angle increment");
    def->category = OptionCategory::support;
    def->tooltip = L("Each layer, add this angle to the interface pattern angle. 0 to keep the same angle, 90 to cross.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 360;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("support_material_interface_fan_speed", coInts);
    def->label = L("Support interface fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all support interfaces, to be able to weaken their bonding with a high fan speed."
        "\nSet to 0 to stop the fan."
        "\nIf disabled, Support Material fan speed will be used."
        "\nCan only be overriden by disable_fan_first_layers.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));


    def = this->add("support_material_interface_contact_loops", coBool);
    def->label = L("Interface loops");
    def->category = OptionCategory::support;
    def->tooltip = L("Cover the top contact layer of the supports with loops. Disabled by default.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("support_material_interface_extruder", coInt);
    def->label = L("Support material/raft interface extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing support material interface "
        "(1+, 0 to use the current extruder to minimize tool changes). This affects raft too.");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("support_material_interface_layers", coInt);
    def->label = L("Top interface layers");
    def->category = OptionCategory::support;
    def->tooltip = L("Number of interface layers to insert between the object(s) and support material.");
    def->sidetext = L("layers");
    def->min = 0;
    def->set_enum_values(ConfigOptionDef::GUIType::i_enum_open, {
        { "0", L("0 (off)") },
        { "1", L("1 (light)") },
        { "2", L("2 (default)") },
        { "3", L("3 (heavy)") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(3));

    def = this->add("support_material_bottom_interface_layers", coInt);
    def->label = L("Bottom interface layers");
    def->category = OptionCategory::support;
    def->tooltip = L("Number of interface layers to insert between the object(s) and support material."
        "\nIf disabled, support_material_interface_layers value is used");
    def->sidetext = L("layers");
    def->min = 0;
    def->can_be_disabled = true;
    def->set_enum_values(ConfigOptionDef::GUIType::i_enum_open, {
    //TRN Print Settings: "Bottom interface layers". Have to be as short as possible
        { "0", L("0 (off)") },
        { "1", L("1 (light)") },
        { "2", L("2 (default)") },
        { "3", L("3 (heavy)") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionInt(0)));

    def = this->add("support_material_closing_radius", coFloat);
    def->label = L("Closing radius");
    def->category = OptionCategory::support;
    def->tooltip = L("For snug supports, the support regions will be merged using morphological closing operation."
        " Gaps smaller than the closing radius will be filled in.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(2));

    def = this->add("support_material_interface_layer_height", coFloatOrPercent);
    def->label = L("Support interface layer height");
    def->category = OptionCategory::support;
    def->tooltip = L("Maximum layer height for the support interface."
        "\nCan be a % of the nozzle diameter"
        "\nIf set to 0, the extruder maximum height will be used.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));

    def = this->add("support_material_interface_spacing", coFloat);
    def->label = L("Interface pattern spacing");
    def->category = OptionCategory::support;
    def->tooltip = L("Spacing between interface lines. Set zero to get a solid interface.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("support_material_interface_speed", coFloatOrPercent);
    def->label = L("Interface");
    def->full_label = L("Support interface speed");
    def->category = OptionCategory::support;
    def->tooltip = L("Speed for printing support material interface layers."
        "\nIf expressed as percentage (for example 50%) it will be calculated over support material speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "support_material_speed";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("support_material_pattern", coEnum);
    def->label = L("Pattern");
    def->full_label = L("Support pattern");
    def->category = OptionCategory::support;
    def->tooltip = L("Pattern used to generate support material.");
    def->set_enum<SupportMaterialPattern>({
        { "rectilinear",        L("Rectilinear") },
        { "rectilinear-grid",   L("Rectilinear grid") },
        { "honeycomb",          L("Honeycomb") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<SupportMaterialPattern>(smpRectilinear));

    def = this->add("support_material_bottom_interface_pattern", coEnum);
    def->label = L("Bottom Pattern");
    def->full_label = L("Support bottom interface pattern");
    def->category = OptionCategory::support;
    def->tooltip = L("Pattern for the bottom interface layers (the ones that start on the object)."
        "\nDefault pattern is the same as the top interface, unless it's Hilbert Curve or Ironing."
        "\nNote that 'Hilbert', 'Ironing' , '(filled)' patterns are really discouraged, and meant to be used with soluble supports and 100% fill interface layer.");
    def->set_enum<InfillPattern>({
        { "auto",              L("Default") },
        { "rectilinear",       L("Rectilinear") },
        { "monotonic",         L("Monotonic") },
        { "concentric",        L("Concentric") },
        { "hilbertcurve",      L("Hilbert Curve") },
        { "smooth",            L("Ironing") },
    });
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipAuto));

    def = this->add("support_material_top_interface_pattern", coEnum);
    def->label = L("Top Pattern");
    def->full_label = L("Support top interface pattern");
    def->category = OptionCategory::support;
    def->tooltip = L("Pattern for the top interface layers."
        "\nNote that 'Hilbert', 'Ironing' and '(filled)' patterns are meant to be used with soluble supports and 100% fill interface layer.");
    def->set_enum<InfillPattern>({
        { "auto",              L("Default") },
        { "rectilinear",       L("Rectilinear") },
        { "monotonic",         L("Monotonic") },
        { "concentric",        L("Concentric") },
        { "sawtooth",          L("Sawtooth") },
        { "hilbertcurve",      L("Hilbert Curve") },
        { "concentricgapfill", L("Concentric (filled)") },
        { "smooth",            L("Ironing") },
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<InfillPattern>(ipRectilinear));
    def->aliases = {"support_material_interface_pattern"};

    def = this->add("support_material_layer_height", coFloatOrPercent);
    def->label = L("Support layer height");
    def->category = OptionCategory::support;
    def->tooltip = L("Maximum layer height for the support, after the first layer that uses the first layer height, and before the interface layers."
        "\nCan be a % of the nozzle diameter"
        "\nIf set to 0, the extruder maximum height will be used.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0, false));
    
    def = this->add("support_material_spacing", coFloat);
    def->label = L("Pattern spacing");
    def->category = OptionCategory::support;
    def->tooltip = L("Spacing between support material lines.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(2.5));

    def = this->add("support_material_speed", coFloatOrPercent);
    def->label = L("Default");
    def->full_label = L("Support speed");
    def->category = OptionCategory::support;
    def->tooltip = L("Speed for printing support material."
        "\nThis can be expressed as a percentage (for example: 80%) over the Default speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "default_speed";
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(60, true));

    def = this->add("support_material_style", coEnum);
    def->label = L("Style");
    def->full_label = L("Support tower style");
    def->category = OptionCategory::support;
    def->tooltip = L("Style and shape of the support towers. Projecting the supports into a regular grid "
        "will create more stable supports, while snug support towers will save material and reduce "
        "object scarring.");
    def->set_enum<SupportMaterialStyle>({
        { "grid", L("Grid") }, 
        { "snug", L("Snug") },
        { "organic", L("Organic") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<SupportMaterialStyle>(smsGrid));

    def = this->add("support_material_synchronize_layers", coBool);
    def->label = L("Synchronize with object layers");
    def->category = OptionCategory::support;
    // TRN PrintSettings : "Synchronize with object layers"
    def->tooltip = L("Synchronize support layers with the object print layers. This is useful "
                   "with multi-material printers, where the extruder switch is expensive. "
                   "This option is only available when top contact Z distance is set to zero.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("support_material_threshold", coInt);
    def->label = L("Overhang threshold");
    def->category = OptionCategory::support;
    def->tooltip = L("Support material will not be generated for overhangs whose slope angle "
                   "(90° = vertical) is above the given threshold. In other words, this value "
                   "represent the most horizontal slope (measured from the horizontal plane) "
                   "that you can print without support material. Set zero for automatic detection "
                   "(recommended).");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 90;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("support_material_with_sheath", coBool);
    def->label = L("With sheath around the support");
    def->category = OptionCategory::support;
    def->tooltip = L("Add a sheath (a single perimeter line) around the base support. This makes "
                   "the support more reliable, but also more difficult to remove.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("support_tree_angle", coFloat);
    def->label = L("Maximum Branch Angle");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Maximum Branch Angle"
    def->tooltip = L("The maximum angle of the branches, when the branches have to avoid the model. "
                     "Use a lower angle to make them more vertical and more stable. Use a higher angle to be able to have more reach.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 85;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(40));

    def = this->add("support_tree_angle_slow", coFloat);
    def->label = L("Preferred Branch Angle");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Preferred Branch Angle"
    def->tooltip = L("The preferred angle of the branches, when they do not have to avoid the model. "
                     "Use a lower angle to make them more vertical and more stable. Use a higher angle for branches to merge faster.");
    def->sidetext = L("°");
    def->min = 10;
    def->max = 85;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(25));

    def = this->add("support_tree_tip_diameter", coFloat);
    def->label = L("Tip Diameter");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Tip Diameter"
    def->tooltip = L("Branch tip diameter for organic supports.");
    def->sidetext = L("mm");
    def->min = 0.1f;
    def->max = 100.f;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.8));

    def = this->add("support_tree_branch_diameter", coFloat);
    def->label = L("Branch Diameter");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Branch Diameter"
    def->tooltip = L("The diameter of the thinnest branches of organic support. Thicker branches are more sturdy. "
                     "Branches towards the base will be thicker than this.");
    def->sidetext = L("mm");
    def->min = 0.1f;
    def->max = 100.f;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(2));

    def = this->add("support_tree_branch_diameter_angle", coFloat);
    // TRN PrintSettings: #lmFIXME 
    def->label = L("Branch Diameter Angle");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Branch Diameter Angle"
    def->tooltip = L("The angle of the branches' diameter as they gradually become thicker towards the bottom. "
                     "An angle of 0 will cause the branches to have uniform thickness over their length. "
                     "A bit of an angle can increase stability of the organic support.");
    def->sidetext = L("°");
    def->min = 0;
    def->max = 15;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(5));

    def = this->add("support_tree_branch_diameter_double_wall", coFloat);
    def->label = L("Branch Diameter with double walls");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Branch Diameter"
    def->tooltip = L("Branches with area larger than the area of a circle of this diameter will be printed with double walls for stability. "
                     "Set this value to zero for no double walls.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 100.f;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(3));

    // Tree Support Branch Distance
    // How far apart the branches need to be when they touch the model. Making this distance small will cause 
    // the tree support to touch the model at more points, causing better overhang but making support harder to remove.
    def = this->add("support_tree_branch_distance", coFloat);
    // TRN PrintSettings: #lmFIXME 
    def->label = L("Branch Distance");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Branch Distance"
    def->tooltip = L("How far apart the branches need to be when they touch the model. "
                     "Making this distance small will cause the tree support to touch the model at more points, "
                     "causing better overhang but making support harder to remove.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("support_tree_top_rate", coPercent);
    def->label = L("Branch Density");
    def->category = OptionCategory::support;
    // TRN PrintSettings: "Organic supports" > "Branch Density"
    def->tooltip = L("Adjusts the density of the support structure used to generate the tips of the branches. "
                     "A higher value results in better overhangs but the supports are harder to remove, "
                     "thus it is recommended to enable top support interfaces instead of a high branch density value "
                     "if dense interfaces are needed.");
    def->sidetext = L("%");
    def->min = 5;
    def->max_literal = {35, false};
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionPercent(15));

    def = this->add("temperature", coInts);
    def->label = L("Other layers");
    def->full_label = L("Temperature");
    def->category = OptionCategory::filament;
    def->tooltip = L("Extruder nozzle temperature for layers after the first one. Set zero to disable "
                   "temperature control commands in the output G-code.");
    def->sidetext = L("°C");
    def->full_label = L("Nozzle temperature");
    def->min = 0;
    def->max = max_temp;
    def->mode = comSimpleAE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionInts { 200 });

    def = this->add("print_temperature", coInt);
    def->label = L("Temperature");
    def->category = OptionCategory::filament;
    def->tooltip = L("Override the temperature of the extruder. Avoid making too many changes, it won't stop for cooling/heating. 0 to disable. May only work on Height range modifiers.");
    def->mode = comExpert | comSuSi;
    def->min = 0;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("print_first_layer_temperature", coInt);
    def->label = L("First Layer Temperature");
    def->category = OptionCategory::filament;
    def->tooltip = L("Override the temperature of the extruder (for the first layer). Avoid making too many changes, it won't stop for cooling/heating. 0 to disable (using print_temperature if defined). May only work on Height range modifiers.");
    def->mode = comExpert | comSuSi;
    def->min = 0;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("print_retract_lift", coFloat);
    def->label = L("Z-lift override");
    def->category = OptionCategory::filament;
    def->tooltip = L("Set the new lift-z value for this override. 0 will disable the z-lift. -& to disable. May only work on Height range modifiers.");
    def->sidetext = L("mm");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(-1));

    def = this->add("thin_perimeters", coPercent);
    def->label = L("Overlapping external perimeter");
    def->full_label = L("Overlapping external perimeter");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Allow outermost perimeter to overlap itself to avoid the use of thin walls. Note that flow isn't adjusted and so this will result in over-extruding and undefined behavior."
                "\n100% means that perimeters can overlap completly on top of each other."
                "\n0% will deactivate this setting."
                "\nValues below 2% don't have any effect."
                "\n-1% will also deactivate the anti-hysteris checks for external perimeters.");
    def->sidetext = "%";
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(80));

    def = this->add("thin_perimeters_all", coPercent);
    def->label = L("Overlapping all perimeters");
    def->full_label = L("Overlapping all perimeters");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Allow all perimeters to overlap, instead of just external ones."
                "\n100% means that perimeters can overlap completly on top of each other."
                "\n0% will deactivate this setting."
                "\nValues below 2% don't have any effect."
                "\n-1% will also deactivate the anti-hysteris checks for internal perimeters.");
    def->sidetext = "%";
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(20));

    def = this->add("thin_walls", coBool);
    def->label = L("Thin walls");
    def->full_label = L("Thin walls");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Detect single-width walls (parts where two extrusions don't fit and we need "
        "to collapse them into a single trace). If unchecked, Slic3r may try to fit perimeters "
        "where it's not possible, creating some overlap leading to over-extrusion.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("thin_walls_min_width", coFloatOrPercent);
    def->label = L("Min width");
    def->full_label = L("Thin walls min width");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Minimum width for the extrusion to be extruded (widths lower than the nozzle diameter will be over-extruded at the nozzle diameter)."
        " If expressed as percentage (for example 110%) it will be computed over nozzle diameter."
        " The default behavior of PrusaSlicer is with a 33% value. Put 100% to avoid any sort of over-extrusion.");
    def->ratio_over = "nozzle_diameter";
    def->mode = comExpert | comSuSi;
    def->min = 0;
    def->max_literal = { 20, true };
    def->set_default_value(new ConfigOptionFloatOrPercent(33, true));

    def = this->add("thin_walls_overlap", coFloatOrPercent);
    def->label = L("Overlap");
    def->full_label = L("Thin wall overlap");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Overlap between the thin wall and the perimeters. Can be a % of the external perimeter width (default 50%)");
    def->ratio_over = "external_perimeter_extrusion_width";
    def->mode = comExpert | comSuSi;
    def->min = 0;
    def->max_literal = { 10, true };
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("thin_walls_merge", coBool);
    def->label = L("Merging with perimeters");
    def->full_label = L("Thin wall merge");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Allow the external perimeter to merge the thin walls in the path."
                    " You can deactivate this if you are using thin walls as a custom support, to reduce adhesion a little.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("thin_walls_acceleration", coFloatOrPercent);
    def->label = L("Thin Walls");
    def->full_label = L("Thin walls acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for thin walls. "
                "\nCan be a % of the external perimeter acceleration"
                "\nSet zero to use external perimeter acceleration for thin walls.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "external_perimeter_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("thin_walls_speed", coFloatOrPercent);
    def->label = L("Thin walls");
    def->full_label = L("Thin walls speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for thin walls (external extrusions that are alone because the obect is too thin at these places)."
        "\nThis can be expressed as a percentage (for example: 80%) over the External Perimeter speed."
        "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "external_perimeter_speed";
    def->min = 0;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("threads", coInt);
    def->label = L("Threads");
    def->tooltip = L("Threads are used to parallelize long-running tasks. Optimal threads number "
                   "is slightly above the number of available cores/processors.");
    def->readonly = true;
    def->min = 1;
    def->mode = comExpert | comPrusa; // note: hidden setting (and should be a preference)
    {
        int threads = (unsigned int)boost::thread::hardware_concurrency();
        def->set_default_value(new ConfigOptionInt(threads > 0 ? threads : 2));
        def->cli = ConfigOptionDef::nocli;
    }

    def = this->add("time_cost", coFloat);
    def->label = L("Time cost");
    def->category = OptionCategory::firmware;
    def->tooltip = L("This setting allows you to set how much an hour of printing time is costing you in printer maintenance, loan, human albor, etc.");
    def->mode = comExpert | comSuSi;
    def->sidetext = L("$ per hour");
    def->min = 0;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("time_estimation_compensation", coPercent);
    def->label = L("Time estimation compensation");
    def->category = OptionCategory::firmware;
    def->tooltip = L("This setting allows you to modify the time estimation by a % amount. As Slic3r only uses the Marlin algorithm, it's not precise enough if another firmware is used.");
    def->mode = comExpert | comSuSi;
    def->sidetext = L("%");
    def->min = 0;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("time_start_gcode", coFloat);
    def->label = L("Time for start custom gcode");
    def->category = OptionCategory::firmware;
    def->tooltip = L("This setting allows you to modify the time estimation by a flat amount to compensate for start script, the homing routine, and other things.");
    def->mode = comExpert | comSuSi;
    def->sidetext = L("s");
    def->min = 0;
    def->set_default_value(new ConfigOptionFloat(20));

    def = this->add("time_toolchange", coFloat);
    def->label = L("Time for toolchange");
    def->category = OptionCategory::firmware;
    def->tooltip = L("This setting allows you to modify the time estimation by a flat amount for each toolchange.");
    def->mode = comExpert | comSuSi;
    def->sidetext = L("s");
    def->min = 0;
    def->set_default_value(new ConfigOptionFloat(30));

    def = this->add("toolchange_gcode", coString);
    def->label = L("Tool change G-code");
    def->category = OptionCategory::customgcode;
    def->tooltip = L("This custom code is inserted at every extruder change. If you don't leave this empty, you are "
        "expected to take care of the toolchange yourself - Slic3r will not output any other G-code to "
        "change the filament. You can use placeholder variables for all Slic3r settings as well as {toolchange_z}, {layer_z}, {layer_num}, {max_layer_z}, {previous_extruder} "
        "and {next_extruder}, so e.g. the standard toolchange command can be scripted as T{next_extruder}."
        "!! Warning !!: if any character is written here, Slic3r won't output any toolchange command by itself.");
    def->multiline = true;
    def->full_width = true;
    def->height = 5;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("tool_name", coStrings);
    def->label = L("Tool name");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Only used for Klipper, where you can name the extruder. If not set, will be 'extruderX' with 'X' replaced by the extruder number.");
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("top_fan_speed", coInts);
    def->label = L("Top Solid fan speed");
    def->category = OptionCategory::cooling;
    def->tooltip = L("This fan speed is enforced during all top fills (including ironing)."
        "\nSet to 0 to stop the fan."
        "\nIf disabled, Solid Infill fan speed will be used."
        "\nCan be disabled by disable_fan_first_layers, slowed down by full_fan_speed_layer.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->can_be_disabled = true;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts({ 100 })));

    def = this->add("top_infill_extrusion_width", coFloatOrPercent);
    def->label = L("Top solid infill");
    def->category = OptionCategory::width;
    def->tooltip = L("Set this to a non-zero value to set a manual extrusion width for infill for top surfaces. "
        "You may want to use thinner extrudates to fill all narrow regions and get a smoother finish. "
        "If left as zero, default extrusion width will be used if set, otherwise nozzle diameter will be used. "
        "If expressed as percentage (for example 110%) it will be computed over nozzle diameter."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(105, true));

    def = this->add("top_infill_extrusion_spacing", coFloatOrPercent);
    def->label = L("Top solid spacing");
    def->category = OptionCategory::width;
    def->tooltip = L("Like Top solid infill width but spacing is the distance between two lines (as they overlap a bit, it's not the same)."
        "\nYou can set either 'Spacing', or 'Width'; the other will be calculated, using default layer height.");
    def->sidetext = L("mm or %");
    def->ratio_over = "nozzle_diameter";
    def->min = 0;
    def->max = 1000;
    def->max_literal = { 10, true };
    def->precision = 6;
    def->can_phony = true;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value((new ConfigOptionFloatOrPercent(0, false))->set_phony(true));

    def = this->add("top_solid_infill_acceleration", coFloatOrPercent);
    def->label = L("Top solid ");
    def->full_label = L("Top solid acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("This is the acceleration your printer will use for top solid infill. "
                "\nCan be a % of the solid infill acceleration"
                "\nSet zero or 100% to use solid infill acceleration for top solid infill.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "solid_infill_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(0,false));

    def = this->add("top_solid_infill_overlap", coPercent);
    def->label = L("Top solid infill overlap");
    def->category = OptionCategory::width;
    def->tooltip = L("This setting allows you to reduce the overlap between the lines of the top solid fill, to reduce the % filled if you see overextrusion signs on solid areas."
        "\nNote that you should be sure that your flow (filament extrusion multiplier) is well calibrated and your filament max overlap is set before thinking to modify this."
        "\nAlso, lowering it below 100% may create visible gaps in the top surfaces"
        "\nSet overlap setting is the only one that can't be reduced by the filament's max overlap.");
    def->sidetext = L("%");
    def->min = 0;
    def->max = 100;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionPercent(100));

    def = this->add("top_solid_infill_speed", coFloatOrPercent);
    def->label = L("Top solid");
    def->full_label = L("Top solid speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for printing top solid layers (it only applies to the uppermost "
                   "external layers and not to their internal solid layers). You may want "
                   "to slow down this to get a nicer surface finish."
                   "\nThis can be expressed as a percentage (for example: 80%) over the Solid Infill speed."
                   "\nSet zero to use autospeed for this feature.");
    def->sidetext = L("mm/s or %");
    def->ratio_over = "solid_infill_speed";
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloatOrPercent(50, true));

    def = this->add("top_solid_layers", coInt);
    //TRN Print Settings: "Top solid layers"
    def->label = L_CONTEXT("Top", "Layers");
    def->full_label = L("Top layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Number of solid layers to generate on top surfaces.");
    def->full_label = L("Top solid layers");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(3));

    def = this->add("top_solid_min_thickness", coFloat);
    def->label = L_CONTEXT("Top", "Layers");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("The number of top solid layers is increased above top_solid_layers if necessary to satisfy "
                     "minimum thickness of top shell."
                     " This is useful to prevent pillowing effect when printing with variable layer height.");
    def->full_label = L("Minimum top shell thickness");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("travel_acceleration", coFloatOrPercent);
    def->label = L("Travel");
    def->full_label = L("Travel acceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("Acceleration for travel moves (jumps between distant extrusion points)."
                     "\nCan be a % of the default acceleration"
                     "\nSet zero to use default acceleration for travel moves.");
    def->sidetext = L("mm/s² or %");
    def->ratio_over = "default_acceleration";
    def->min = 0;
    def->max_literal = { -200, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(1500, false));

    def = this->add("travel_deceleration_use_target", coBool);
    def->label = L("Decelerate with target acceleration");
    def->full_label = L("Use target acceleration for travel deceleration");
    def->category = OptionCategory::speed;
    def->tooltip = L("If selected, the deceleration of a travel will use the acceleration value of the extrusion that will be printed after it (if any) ");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("travel_speed", coFloat);
    def->label = L("Travel");
    def->full_label = L("Travel speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for travel moves (jumps between distant extrusion points).");
    def->sidetext = L("mm/s");
    def->aliases = { "travel_feed_rate" };
    def->min = 1;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(130));

    def = this->add("travel_speed_z", coFloat);
    def->label = L("Z Travel");
    def->full_label = L("Z travel speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Speed for movements along the Z axis.\nWhen set to zero, the value "
                     "is ignored and regular travel speed is used instead.");
    def->sidetext = L("mm/s");
    def->aliases = { "travel_feed_rate_z" };
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("use_firmware_retraction", coBool);
    def->label = L("Use firmware retraction");
    def->category = OptionCategory::general;
    def->tooltip = L("This setting uses G10 and G11 commands to have the firmware "
                   "handle the retraction. Note that this has to be supported by firmware.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("use_relative_e_distances", coBool);
    def->label = L("Use relative E distances");
    def->category = OptionCategory::general;
    def->tooltip = L("If your firmware requires relative E values, check this, "
                   "otherwise leave it unchecked. Most firmwares use absolute values.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("use_volumetric_e", coBool);
    def->label = L("Use volumetric E");
    def->category = OptionCategory::general;
    def->tooltip = L("This experimental setting uses outputs the E values in cubic millimeters "
                   "instead of linear millimeters. If your firmware doesn't already know "
                   "filament diameter(s), you can put commands like 'M200 D{filament_diameter_0} T0' "
                   "in your start G-code in order to turn volumetric mode on and use the filament "
                   "diameter associated to the filament selected in Slic3r. This is only supported "
                   "in recent Marlin.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("variable_layer_height", coBool);
    def->label = L("Enable variable layer height feature");
    def->category = OptionCategory::general;
    def->tooltip = L("Some printers or printer setups may have difficulties printing "
                   "with a variable layer height. Enabled by default.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("wipe", coBools);
    def->label = L("Wipe while retracting");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This flag will move the nozzle while retracting to minimize the possible blob on leaky extruders."
        "\nNote that as a wipe only happens when there is a retraction, the 'only retract when crossing perimeters' print setting can greatly reduce the number of wipes.");
    def->mode = comAdvancedE | comPrusa;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ false });

    def = this->add("wipe_inside_start", coBools);
    def->label = L("Wipe inside at start");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Before extruding an external perimeter, this flag will place the nozzle a bit inward and in advance of the seam position before unretracting."
        " It will then move to the seam position before extruding.");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ false });

    def = this->add("wipe_inside_end", coBools);
    def->label = L("Wipe inside at end");
    def->category = OptionCategory::extruders;
    def->tooltip = L("This flag will wipe the nozzle a bit inward after extruding an external perimeter."
        " The wipe_extra_perimeter is executed first, then this move inward before the retraction wipe."
        " Note that the retraction wipe will follow the exact external perimeter (center) line if this parameter is disabled, and will follow the inner side of the external perimeter line if enabled");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ true });

    def = this->add("wipe_inside_depth", coPercents);
    def->label = L("Max Wipe deviation");
    def->full_label = L("Maximum Wipe deviation to the inside");
    def->category = OptionCategory::extruders;
    def->tooltip = L("By how much the 'wipe inside' can dive inside the object (if possible)?"
        "\nIn % of the perimeter width."
        "\nNote: don't put a value higher than 50% if you have only one perimeter, or 150% for two perimeter, etc... or it will ooze instead of wipe.");
    def->sidetext = L("%");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPercents{ 50 });

    def = this->add("wipe_only_crossing", coBools);
    def->label = L("Wipe only when crossing perimeters");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Don't wipe when you don't cross a perimeter. Need 'avoid_crossing_perimeters' and 'wipe' enabled.");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionBools{ true });

    def = this->add("wipe_speed", coFloats);
    def->label = L("Wipe speed");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Speed in mm/s of the wipe. If it's faster, it will try to go further away, as the wipe time is set by ( 100% - 'retract before wipe') * 'retaction length' / 'retraction speed'."
        "\nIf set to zero, the travel speed is used.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloats{ 0 });

    def = this->add("wipe_tower", coBool);
    def->label = L("Enable");
    def->full_label = L("Enable wipe tower");
    def->category = OptionCategory::general;
    def->tooltip = L("Multi material printers may need to prime or purge extruders on tool changes. "
                   "Extrude the excess material into the wipe tower.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("wipe_tower_speed", coFloat);
    def->label = L("Wipe Tower Speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Printing speed of the wipe tower. Capped by filament_max_volumetric_speed (if set)."
        "\nIf set to zero, a value of 80mm/s is used.");
    def->sidetext = L("mm/s");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(80.));

    def = this->add("wipe_tower_wipe_starting_speed", coFloatOrPercent);
    def->label = L("Wipe tower wipe starting speed");
    def->category = OptionCategory::speed;
    def->tooltip = L("Start of the wiping speed ramp up (for wipe tower)."
        "\nCan be a % of the 'Wipe tower main speed'."
        "\nSet to 0 to disable.");
    def->sidetext = L("mm/s or %");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(33, true));


    def = this->add("wiping_volumes_extruders", coFloats);
    def->label = L("Purging volumes - load/unload volumes");
    def->tooltip = L("This vector saves required volumes to change from/to each tool used on the "
                     "wipe tower. These values are used to simplify creation of the full purging "
                     "volumes below. ");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloats { 70., 70., 70., 70., 70., 70., 70., 70., 70., 70.  });

    def = this->add("wiping_volumes_matrix", coFloats);
    def->label = L("Purging volumes - matrix");
    def->tooltip = L("This matrix describes volumes (in cubic milimetres) required to purge the"
                     " new filament on the wipe tower for any given pair of tools. ");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloats {   0., 140., 140., 140., 140.,
                                                    140.,   0., 140., 140., 140.,
                                                    140., 140.,   0., 140., 140.,
                                                    140., 140., 140.,   0., 140.,
                                                    140., 140., 140., 140.,   0. });


    def = this->add("wipe_advanced", coBool);
    def->label = L("Enable advanced wiping volume");
    def->tooltip = L("Allow Slic3r to compute the purge volume via smart computations. Use the pigment% of each filament and following parameters");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("wipe_advanced_nozzle_melted_volume", coFloat);
    def->label = L("Nozzle volume");
    def->tooltip = L("The volume of melted plastic inside your nozzle. Used by 'advanced wiping'.");
    def->sidetext = L("mm3");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(120));

    def = this->add("filament_wipe_advanced_pigment", coFloats);
    def->label = L("Pigment percentage");
    def->tooltip = L("The pigment % for this filament (bewteen 0 and 1, 1=100%). 0 for translucent/natural, 0.2-0.5 for white and 1 for black.");
    def->min = 0;
    def->max = 1;
    def->mode = comExpert | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0.5 });

    def = this->add("wipe_advanced_multiplier", coFloat);
    def->label = L("Multiplier");
    def->full_label = L("Auto-wipe multiplier");
    def->tooltip = L("The volume multiplier used to compute the final volume to extrude by the algorithm.");
    def->sidetext = L("mm3");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(60));


    def = this->add("wipe_advanced_algo", coEnum);
    def->label = L("Algorithm");
    def->full_label = L("Auto-wipe algorithm");
    def->tooltip = L("Algorithm for the advanced wipe.\n"
        "Linear : volume = nozzle + volume_mult * (pigmentBefore-pigmentAfter)\n"
        "Quadratic: volume = nozzle + volume_mult * (pigmentBefore-pigmentAfter)+ volume_mult * (pigmentBefore-pigmentAfter)^3\n"
        "Hyperbola: volume = nozzle + volume_mult * (0.5+pigmentBefore) / (0.5+pigmentAfter)");
    def->set_enum<WipeAlgo>({
        { "linear", L("Linear") }, 
        { "quadra", L("Quadratric") },
        { "expo",   L("Hyperbola") }
    });
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionEnum<WipeAlgo>(waLinear));

    def = this->add("wipe_tower_brim_width", coFloatOrPercent);
    def->label = L("Wipe tower brim width");
    def->tooltip = L("Width of the brim for the wipe tower. Can be in mm or in % of the (assumed) only one nozzle diameter.");
    def->ratio_over = "nozzle_diameter";
    def->mode = comAdvancedE | comPrusa;
    def->min = 0;
    def->max_literal = { 100, true };
    def->aliases = { "wipe_tower_brim" }; // SuperSlicer 2.3 and before
    def->set_default_value(new ConfigOptionFloatOrPercent(2,false));

    def = this->add("wipe_tower_no_sparse_layers", coBool);
    def->label = L("No sparse layers (EXPERIMENTAL)");
    def->category = OptionCategory::mmsetup;
    def->tooltip = L("If enabled, the wipe tower will not be printed on layers with no toolchanges. "
                     "On layers with a toolchange, extruder will travel downward to print the wipe tower. "
                     "User is responsible for ensuring there is no collision with the print.");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));
    def = this->add("wipe_tower_x", coFloat);
    def->label = L("X");
    def->full_label = L("Wipe tower X");
    def->tooltip = L("X coordinate of the left front corner of a wipe tower");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(180.));


    def = this->add("wipe_tower_y", coFloat);
    def->label = L("Y");
    def->full_label = L("Wipe tower Y");
    def->tooltip = L("Y coordinate of the left front corner of a wipe tower");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(140.));

    def = this->add("wipe_tower_width", coFloat);
    def->label = L("Width");
    def->full_label = L("Wipe tower Width");
    def->tooltip = L("Width of a wipe tower");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(60.));

    def = this->add("wipe_tower_rotation_angle", coFloat);
    def->label = L("Wipe tower rotation angle");
    def->tooltip = L("Wipe tower rotation angle with respect to x-axis.");
    def->sidetext = L("°");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("wipe_tower_cone_angle", coFloat);
    def->label = L("Stabilization cone apex angle");
    def->tooltip = L("Angle at the apex of the cone that is used to stabilize the wipe tower. "
                     "Larger angle means wider base.");
    def->sidetext = L("°");
    def->mode = comAdvancedE |comPrusa;
    def->min = 0.;
    def->max = 90.;
    def->set_default_value(new ConfigOptionFloat(0.));

    def = this->add("wipe_tower_extra_spacing", coPercent);
    def->label = L("Wipe tower purge lines spacing");
    def->tooltip = L("Spacing of purge lines on the wipe tower.");
    def->sidetext = L("%");
    def->mode = comExpert |comPrusa;
    def->min = 100.;
    def->max = 300.;
    def->set_default_value(new ConfigOptionPercent(100.));

    def = this->add("wipe_into_infill", coBool);
    def->category = OptionCategory::wipe;
    def->label = L("Wipe into this object's infill");
    def->tooltip = L("Purging after toolchange will be done inside this object's infills. "
                     "This lowers the amount of waste but may result in longer print time "
                     " due to additional travel moves.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("wipe_into_objects", coBool);
    def->category = OptionCategory::wipe;
    def->label = L("Wipe into this object");
    def->tooltip = L("Object will be used to purge the nozzle after a toolchange to save material "
        "that would otherwise end up in the wipe tower and decrease print time. "
        "Colours of the objects will be mixed as a result.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("wipe_extra_perimeter", coFloats);
    def->category = OptionCategory::extruders;
    def->label = L("Extra Wipe for external perimeters");
    def->tooltip = L("When the external perimeter loop extrusion ends, a wipe is done, going slightly inside the print."
        " The number in this settting increases the wipe by moving the nozzle along the loop again before the final wipe.");
    def->min = 0;
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats{ 0.f });

    def = this->add("wipe_tower_bridging", coFloat);
    def->label = L("Maximal bridging distance");
    def->tooltip = L("Maximal distance between supports on sparse infill sections. ");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(10.));

    def = this->add("wipe_tower_extruder", coInt);
    def->label = L("Wipe tower extruder");
    def->category = OptionCategory::extruders;
    def->tooltip = L("The extruder to use when printing perimeter of the wipe tower. "
                     "Set to 0 to use the one that is available (non-soluble would be preferred).");
    def->min = 0;
    def->mode = comAdvancedE |comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("solid_infill_every_layers", coInt);
    def->label = L("Solid infill every");
    def->category = OptionCategory::infill;
    def->tooltip = L("This feature allows to force a solid layer every given number of layers. "
                   "Zero to disable. You can set this to any value (for example 9999); "
                   "Slic3r will automatically choose the maximum possible number of layers "
                   "to combine according to nozzle diameter and layer height.");
    def->sidetext = L("layers");
    def->min = 0;
    def->mode = comExpert |comPrusa;
    def->set_default_value(new ConfigOptionInt(0));

    def = this->add("xy_size_compensation", coFloat);
    def->label = L("Outer");
    def->full_label = L("Outer XY size compensation");
    def->category = OptionCategory::slicing;
    def->tooltip = L("The object will be grown/shrunk in the XY plane by the configured value "
        "(negative = inwards = remove area, positive = outwards = add area). This might be useful for fine-tuning sizes."
        "\nThis one only applies to the 'exterior' shell of the object."
        "\n !!! it's recommended you put the same value into the 'Inner XY size compensation', unless you are sure you don't have horizontal holes. !!! ");
    def->sidetext = L("mm");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("xy_inner_size_compensation", coFloat);
    def->label = L("Inner");
    def->full_label = L("Inner XY size compensation");
    def->category = OptionCategory::slicing;
    def->tooltip = L("The object will be grown/shrunk in the XY plane by the configured value "
        "(negative = inwards = remove area, positive = outwards = add area). This might be useful for fine-tuning sizes."
        "\nThis one only applies to the 'inner' shell of the object (!!! horizontal holes break the shell !!!)");
    def->sidetext = L("mm");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("hole_size_compensation", coFloat);
    def->label = L("XY compensation");
    def->full_label = L("XY holes compensation");
    def->category = OptionCategory::slicing;
    def->tooltip = L("The convex holes will be grown / shrunk in the XY plane by the configured value"
        " (negative = inwards = remove area, positive = outwards = add area, should be negative as the holes are always a bit smaller irl)."
        " This might be useful for fine-tuning hole sizes."
        "\nThis setting behaves the same as 'Inner XY size compensation' but only for convex shapes. It's added to 'Inner XY size compensation', it does not replace it. ");
    def->sidetext = L("mm");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("hole_size_threshold", coFloat);
    def->label = L("Threshold");
    def->full_label = L("XY holes threshold");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Maximum area for the hole where the hole_size_compensation will apply fully."
            " After that, it will decrease down to 0 for four times this area."
            " Set to 0 to let the hole_size_compensation apply fully for all detected holes");
    def->sidetext = L("mm²");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(100));

    def = this->add("hole_to_polyhole", coBool);
    def->label = L("Convert round holes to polyholes");
    def->full_label = L("Convert round holes to polyholes");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Search for almost-circular holes that span more than one layer and convert the geometry to polyholes."
        " Use the nozzle size and the (biggest) diameter to compute the polyhole."
        "\nSee http://hydraraptor.blogspot.com/2011/02/polyholes.html");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("hole_to_polyhole_threshold", coFloatOrPercent);
    def->label = L("Roundness margin");
    def->full_label = L("Polyhole detection margin");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Maximum deflection of a point to the estimated radius of the circle."
        "\nAs cylinders are often exported as triangles of varying size, points may not be on the circle circumference."
        " This setting allows you some leeway to broaden the detection."
        "\nIn mm or in % of the radius.");
    def->sidetext = L("mm or %");
    def->max_literal = { 10, false};
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(0.01, false));

    def = this->add("hole_to_polyhole_twisted", coBool);
    def->label = L("Twisting");
    def->full_label = L("Polyhole twist");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Rotate the polyhole every layer.");
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("z_offset", coFloat);
    def->label = L("Z offset");
    def->category = OptionCategory::general;
    def->tooltip = L("This value will be added (or subtracted) from all the Z coordinates "
                   "in the output G-code. It is used to compensate for bad Z endstop position: "
                   "for example, if your endstop zero actually leaves the nozzle 0.3mm far "
                   "from the print bed, set this to -0.3 (or fix your endstop).");
    def->sidetext = L("mm");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("z_step", coFloat);
    def->label = L("Z full step");
    def->tooltip = L("Set this to the height moved when your Z motor (or equivalent) turns one step."
                    "If your motor needs 200 steps to move your head/platter by 1mm, this field should be 1/200 = 0.005."
                    "\nNote that the gcode will write the z values with 6 digits after the dot if z_step is set (it's 3 digits if it's disabled)."
                    "\nSet zero to disable.");
    def->cli = "z-step=f";
    def->sidetext = L("mm");
    def->min = 0;
    def->precision = 8;
    def->mode = comExpert | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0.005));

    def = this->add("init_z_rotate", coFloat);
    def->label = L("Preferred orientation");
    def->category = OptionCategory::general;
    def->tooltip = L("Rotate stl around z axes while adding them to the bed.");
    def->sidetext = L("°");
    def->min = -360;
    def->max = 360;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(0.0));

    //// Arachne slicing, put it in alpha order when out of experimental

    def = this->add("perimeter_generator", coEnum);
    def->label = L("Perimeter generator");
    def->category = OptionCategory::perimeter;
    def->tooltip = L("Classic perimeter generator produces perimeters with constant extrusion width and for "
                      "very thin areas is used gap-fill. "
                      "Arachne engine produces perimeters with variable extrusion width. "
                      "This setting also affects the Concentric infill.");
    def->set_enum<PerimeterGeneratorType>({
        { "classic", L("Classic") },
        { "arachne", L("Arachne") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<PerimeterGeneratorType>(PerimeterGeneratorType::Classic));

    def = this->add("wall_transition_length", coFloatOrPercent);
    def->label = L("Perimeter transition length");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("When transitioning between different numbers of perimeters as the part becomes"
                       "thinner, a certain amount of space is allotted to split or join the perimeter segments. "
                       "If expressed as a percentage (for example 100%), it will be computed based on the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloatOrPercent(100, true));

    def = this->add("wall_transition_filter_deviation", coFloatOrPercent);
    def->label = L("Perimeter transitioning filter margin");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("Prevent transitioning back and forth between one extra perimeter and one less. This "
                       "margin extends the range of extrusion widths which follow to [Minimum perimeter width "
                       "- margin, 2 * Minimum perimeter width + margin]. Increasing this margin "
                       "reduces the number of transitions, which reduces the number of extrusion "
                       "starts/stops and travel time. However, large extrusion width variation can lead to "
                       "under- or overextrusion problems."
                       "If expressed as percentage (for example 25%), it will be computed over nozzle diameter.");
    def->sidetext = L("mm");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloatOrPercent(25, true));

    def = this->add("wall_transition_angle", coFloat);
    def->label = L("Perimeter transitioning threshold angle");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("When to create transitions between even and odd numbers of perimeters. A wedge shape with"
                       " an angle greater than this setting will not have transitions and no perimeters will be "
                       "printed in the center to fill the remaining space. Reducing this setting reduces "
                       "the number and length of these center perimeters, but may leave gaps or overextrude.");
    def->sidetext = L("°");
    def->mode = comExpert | comPrusa;
    def->min = 1.;
    def->max = 59.;
    def->set_default_value(new ConfigOptionFloat(10.));

    def = this->add("wall_distribution_count", coInt);
    def->label = L("Perimeter distribution count");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("The number of perimeters, counted from the center, over which the variation needs to be "
                       "spread. Lower values mean that the outer perimeters don't change in width.");
    def->mode = comExpert | comPrusa;
    def->min = 1;
    def->set_default_value(new ConfigOptionInt(1));

    def = this->add("min_feature_size", coFloatOrPercent);
    def->label = L("Minimum feature size");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("Minimum thickness of thin features. Model features that are thinner than this value will "
                       "not be printed, while features thicker than the Minimum feature size will be widened to "
                       "the Minimum perimeter width. "
                       "If expressed as a percentage (for example 25%), it will be computed based on the nozzle diameter.");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloatOrPercent(25, true));

    def = this->add("min_bead_width", coFloatOrPercent);
    def->label = L("Minimum perimeter width");
    def->category = OptionCategory::advanced;
    def->tooltip  = L("Width of the perimeter that will replace thin features (according to the Minimum feature size) "
                       "of the model. If the Minimum perimeter width is thinner than the thickness of the feature,"
                       " the perimeter will become as thick as the feature itself. "
                       "If expressed as percentage (for example 85%), it will be computed over nozzle diameter.");
    def->sidetext = L("mm or %");
    def->mode = comExpert | comPrusa;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloatOrPercent(85, true));

    /////// End of Arachne settings /////

    // Declare retract values for filament profile, overriding the printer's extruder profile.
    for (const std::string &opt_key : m_filament_override_option_keys
        //{
        //// floats
        //"retract_length", "retract_lift", "retract_lift_above", "retract_lift_below", "retract_speed", 
        //"travel_max_lift",
        //"deretract_speed", "retract_restart_extra", "retract_before_travel", "retract_lift_before_travel",
        //"retract_length_toolchange", "retract_restart_extra_toolchange",
        //"wipe_extra_perimeter", "wipe_speed",
        //"wipe_inside_depth", "wipe_inside_end", "wipe_inside_start",
        //// bools
        //"retract_layer_change", "wipe", "wipe_only_crossing",
        //"travel_lift_before_obstacle", "travel_ramping_lift", "travel_slope",
        //// percents
        //"retract_before_wipe", "travel_slope",
        //// floatsOrPercents
        //"seam_gap"
        //}

        ) {
        auto it_opt = options.find(opt_key);
        assert(it_opt != options.end());
        def = this->add(std::string("filament_") + opt_key, it_opt->second.type);
        def->can_be_disabled = true;
        def->is_optional = true;
        def->label      = it_opt->second.label;
        def->full_label = it_opt->second.full_label;
        def->tooltip    = it_opt->second.tooltip;
        def->sidetext   = it_opt->second.sidetext;
        def->mode       = it_opt->second.mode;
        // create default value with the default value is taken from the default value of the config.
        // put a disbaled value as first entry.
        switch (def->type) {
        case coBools: {
            ConfigOptionBools *opt = new ConfigOptionBools({it_opt->second.default_value.get()->get_bool()});
            opt->set_can_be_disabled(true);
            def->set_default_value(opt);
            break;
        }
        case coFloats: {
            ConfigOptionFloats *opt = new ConfigOptionFloats({it_opt->second.default_value.get()->get_float()});
            opt->set_can_be_disabled(true);
            def->set_default_value(opt);
            break;
        }
        case coPercents: {
            ConfigOptionPercents *opt = new ConfigOptionPercents({it_opt->second.default_value.get()->get_float()});
            opt->set_can_be_disabled(true);
            def->set_default_value(opt);
            break;
        }
        case coFloatsOrPercents: {
            ConfigOptionFloatsOrPercents*opt = new ConfigOptionFloatsOrPercents(
                {static_cast<const ConfigOptionFloatsOrPercents*>(it_opt->second.default_value.get())->get_at(0)});
            opt->set_can_be_disabled(true);
            def->set_default_value(opt);
            break;
        }
        default: assert(false);
        }
        assert(!def->default_value->is_enabled());
    }
}

void PrintConfigDef::init_extruder_option_keys()
{
    // ConfigOptionFloats, ConfigOptionPercents, ConfigOptionBools, ConfigOptionStrings
    m_extruder_option_keys = {
        "default_filament_profile",
        "deretract_speed",
        "extruder_colour",
        "extruder_extrusion_multiplier_speed",
        "extruder_fan_offset",
        "extruder_offset",
        "extruder_temperature_offset",
        "fan_name",
        "max_layer_height",
        "min_layer_height",
        "nozzle_diameter",
        "retract_before_travel",
        "retract_before_wipe",
        "retract_layer_change",
        "retract_length",
        "retract_length_toolchange",
        "retract_lift",
        "retract_lift_above",
        "retract_lift_before_travel",
        "retract_lift_below",
        "retract_lift_first_layer",
        "retract_lift_top",
        "retract_restart_extra",
        "retract_restart_extra_toolchange",
        "retract_speed",
        "seam_gap",
        "seam_gap_external",
        "tool_name",
        "travel_lift_before_obstacle",
        // "travel_max_lift",
        "travel_ramping_lift",
        "travel_slope",
        "wipe",
        "wipe_extra_perimeter",
        "wipe_inside_depth",
        "wipe_inside_end",
        "wipe_inside_start",
        "wipe_only_crossing",
        "wipe_speed",
    };
    assert(std::is_sorted(m_extruder_option_keys.begin(), m_extruder_option_keys.end()));

    m_extruder_retract_keys = {
        "deretract_speed",
        "retract_before_travel",
        "retract_before_wipe",
        "retract_layer_change",
        "retract_length",
        "retract_length_toolchange",
        "retract_lift",
        "retract_lift_above",
        "retract_lift_before_travel",
        "retract_lift_below",
        "retract_lift_first_layer",
        "retract_lift_top",
        "retract_restart_extra",
        "retract_restart_extra_toolchange",
        "retract_speed",
        "seam_gap",
        "seam_gap_external",
        "travel_lift_before_obstacle",
        // "travel_max_lift",
        "travel_ramping_lift",
        "travel_slope",
        "wipe",
        "wipe_extra_perimeter",
        "wipe_inside_depth",
        "wipe_inside_end",
        "wipe_inside_start",
        "wipe_only_crossing",
        "wipe_speed",
    };
    assert(std::is_sorted(m_extruder_retract_keys.begin(), m_extruder_retract_keys.end()));
    m_filament_override_option_keys = {
        "deretract_speed",
        "retract_before_travel",
        "retract_before_wipe",
        "retract_layer_change",
        "retract_length",
        "retract_length_toolchange",
        "retract_lift",
        "retract_lift_above",
        "retract_lift_before_travel",
        "retract_lift_below",
        "retract_restart_extra",
        "retract_restart_extra_toolchange",
        "retract_speed",
        "seam_gap",
        "travel_lift_before_obstacle",
        // "travel_max_lift",
        "travel_ramping_lift",
        "travel_slope",
        "wipe",
        "wipe_extra_perimeter",
        "wipe_inside_depth",
        "wipe_inside_end",
        "wipe_inside_start",
        "wipe_only_crossing",
        "wipe_speed",
    };
}

void PrintConfigDef::init_milling_params()
{
    // ConfigOptionFloats, ConfigOptionPercents, ConfigOptionBools, ConfigOptionStrings
    m_milling_option_keys = {
        "milling_diameter",
        "milling_toolchange_end_gcode",
        "milling_toolchange_start_gcode",
        //"milling_offset",
        //"milling_z_offset",
        "milling_z_lift",

    };

    ConfigOptionDef* def;

    // Milling Printer settings

    def = this->add("milling_cutter", coInt);
    def->label = L("Milling cutter");
    def->category = OptionCategory::general;
    def->tooltip = L("The milling cutter to use (unless more specific extruder settings are specified). ");
    def->min = 0;  // 0 = inherit defaults
    def->set_enum_values(ConfigOptionDef::GUIType::i_enum_open, {
    //TRN Print Settings: "Bottom contact Z distance". Have to be as short as possible
        { "default",      L("Default") },
        { "1",    "1" },
        { "2",    "2" },
        { "3",    "3" },
        { "4",    "4" },
        { "5",    "5" },
        { "6",    "6" },
        { "7",    "7" },
        { "8",    "8" },
        { "9",    "9" },
    });

    def = this->add("milling_diameter", coFloats);
    def->label = L("Milling diameter");
    def->category = OptionCategory::milling_extruders;
    def->tooltip = L("This is the diameter of your cutting tool.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats(3.14));

    def = this->add("milling_offset", coPoints);
    def->label = L("Tool offset");
    def->category = OptionCategory::extruders;
    def->tooltip = L("If your firmware doesn't handle the extruder displacement you need the G-code "
        "to take it into account. This option lets you specify the displacement of each extruder "
        "with respect to the first one. It expects positive coordinates (they will be subtracted "
        "from the XY coordinate).");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPoints{ Vec2d(0, 0) });

    def = this->add("milling_z_offset", coFloats);
    def->label = L("Tool z offset");
    def->category = OptionCategory::extruders;
    def->tooltip = L(".");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats(0));

    def = this->add("milling_z_lift", coFloats);
    def->label = L("Tool z lift");
    def->category = OptionCategory::extruders;
    def->tooltip = L("Amount of lift for travel.");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloats(2));

    def = this->add("milling_toolchange_start_gcode", coStrings);
    def->label = L("G-Code to switch to this toolhead");
    def->category = OptionCategory::milling_extruders;
    def->tooltip = L("Put here the gcode to change the toolhead (called after the g-code T{next_extruder}). You have access to {next_extruder} and {previous_extruder}."
        " next_extruder is the 'extruder number' of the new milling tool, it's equal to the index (begining at 0) of the milling tool plus the number of extruders."
        " previous_extruder is the 'extruder number' of the previous tool, it may be a normal extruder, if it's below the number of extruders."
        " The number of extruder is available at {extruder} and the number of milling tool is available at {milling_cutter}.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("milling_toolchange_end_gcode", coStrings);
    def->label = L("G-Code to switch from this toolhead");
    def->category = OptionCategory::milling_extruders;
    def->tooltip = L("Enter here the gcode to end the toolhead action, like stopping the spindle. You have access to {next_extruder} and {previous_extruder}."
        " previous_extruder is the 'extruder number' of the current milling tool, it's equal to the index (begining at 0) of the milling tool plus the number of extruders."
        " next_extruder is the 'extruder number' of the next tool, it may be a normal extruder, if it's below the number of extruders."
        " The number of extruder is available at {extruder}and the number of milling tool is available at {milling_cutter}.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("milling_post_process", coBool);
    def->label = L("Milling post-processing");
    def->category = OptionCategory::milling;
    def->tooltip = L("If activated, at the end of each layer, the printer will switch to a milling head and mill the external perimeters."
        "\nYou should set the 'Milling extra XY size' to a value high enough to have enough plastic to mill. Also, be sure that your piece is firmly glued to the bed.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("milling_extra_size", coFloatOrPercent);
    def->label = L("Milling extra XY size");
    def->category = OptionCategory::milling;
    def->tooltip = L("This increases the size of the object by a certain amount to have enough plastic to mill."
        " You can set a number of mm or a percentage of the calculated optimal extra width (from flow calculation).");
    def->sidetext = L("mm or %");
    def->ratio_over = "computed_on_the_fly";
    def->max_literal = { 20, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(150, true));

    def = this->add("milling_after_z", coFloatOrPercent);
    def->label = L("Milling only after");
    def->category = OptionCategory::milling;
    def->tooltip = L("This setting restricts the post-process milling to a certain height, to avoid milling the bed. It can be a mm or a % of the first layer height (so it can depend on the object).");
    def->sidetext = L("mm or %");
    def->ratio_over = "first_layer_height";
    def->max_literal = { 10, false };
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloatOrPercent(200, true));

    def = this->add("milling_speed", coFloat);
    def->label = L("Milling Speed");
    def->category = OptionCategory::milling;
    def->tooltip = L("Speed for milling tool.");
    def->sidetext = L("mm/s");
    def->min = 0;
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(30));
}

void PrintConfigDef::init_laser_params()
{
    // ConfigOptionFloats, ConfigOptionPercents, ConfigOptionBools, ConfigOptionStrings
    m_laser_option_keys = {
        //"laser_diameter",
        "laser_toolchange_end_gcode",
        "laser_toolchange_start_gcode",
        "laser_offset",
        "laser_z_offset",
        "laser_power",
        "laser_disable_gcode",
        "laser_enable_gcode"

    };

    ConfigOptionDef* def;

    // laser Printer settings

    def = this->add("laser_head", coInt);
    def->label = L("laser cutter");
    def->category = OptionCategory::general;
    def->tooltip = L("The laser head to use (unless more specific tool settings are specified). ");
    def->min = 0;  // 0 = inherit defaults
    def->set_enum_labels(ConfigOptionDef::GUIType::i_enum_open, 
        { L("default"), "1", "2", "3", "4", "5", "6", "7", "8", "9" }); // override label for item 0

    //def = this->add("laser_diameter", coFloats);
    //def->label = L("laser diameter");
    //def->category = OptionCategory::laser_tools;
    //def->tooltip = L("This is the diameter of the laser on the focus point.");
    //def->sidetext = L("mm");
    //def->mode = comAdvancedE | comSuSi;
    //def->is_vector_extruder = true;
    //def->set_default_value(new ConfigOptionFloats(3.14));

    def = this->add("laser_offset", coPoints);
    def->label = L("Tool offset");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("If your firmware doesn't handle the tool displacement you need the G-code "
        "to take it into account. This option lets you specify the displacement of each tool "
        "with respect to the first one. It expects positive coordinates (they will be subtracted "
        "from the XY coordinate).");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionPoints{ Vec2d(0,0) });

    def = this->add("laser_z_offset", coFloats);
    def->label = L("Tool z offset");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L(".");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionFloats(0));

    def = this->add("laser_power", coFloats);
    def->label = L("Laser power");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("Power in watt from the laser head when activated.");
    def->sidetext = L("W");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloats(2));

    def = this->add("laser_toolchange_start_gcode", coStrings);
    def->label = L("G-Code to switch to this toolhead");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("Put here the gcode to change the toolhead (called after the g-code T{next_extruder}). You have access to {next_extruder} and {previous_extruder}."
        " next_extruder is the 'extruder number' of the new laser tool, it's equal to the index (begining at 0) of the laser tool plus the number of extruders."
        " previous_extruder is the 'extruder number' of the previous tool, it may be a normal extruder, if it's below the number of extruders."
        " The number of extruder is available at {extruder} and the number of laser tool is available at {laser_cutter}.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("laser_toolchange_end_gcode", coStrings);
    def->label = L("G-Code to switch from this toolhead");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("Enter here the gcode to end the toolhead action, like stopping the spindle. You have access to {next_extruder} and {previous_extruder}."
        " previous_extruder is the 'extruder number' of the current laser tool, it's equal to the index (begining at 0) of the laser tool plus the number of extruders."
        " next_extruder is the 'extruder number' of the next tool, it may be a normal extruder, if it's below the number of extruders."
        " The number of extruder is available at {extruder}and the number of laser tool is available at {laser_cutter}.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("laser_enable_gcode", coStrings);
    def->label = L("G-Code to fire the laser");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("The gcode to activate the laser .");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    def = this->add("laser_disable_gcode", coStrings);
    def->label = L("G-Code to turn off the laser");
    def->category = OptionCategory::laser_tools;
    def->tooltip = L("The gcode to disable the laser.");
    def->multiline = true;
    def->full_width = true;
    def->height = 12;
    def->mode = comAdvancedE | comSuSi;
    def->is_vector_extruder = true;
    def->set_default_value(new ConfigOptionStrings(""));

    // laser Print settings

    def = this->add("laser_support_interface_pp", coBool);
    def->label = L("laser post-processing");
    def->category = OptionCategory::laser;
    def->tooltip = L("If activated, at the end of each layer with support interfec, the printer will switch to a laser head and will pass over the interface to charcoal them.");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("laser_energy", coFloat);
    def->label = L("laser enrgy");
    def->category = OptionCategory::laser;
    def->tooltip = L("How much energy to use per mm of travel? Increase it to make the laser move slower to tranfert more enrgy. If your laser has a power of 1W, set this to 0.1 so it takes 1 second to do 10mm of post-processing.");
    def->sidetext = L("Ws/mm");
    def->mode = comAdvancedE | comSuSi;
    def->set_default_value(new ConfigOptionFloat(1.f));
}

void PrintConfigDef::init_sla_support_params(const std::string &prefix)
{
    ConfigOptionDef* def;

    def = this->add(prefix + "support_head_front_diameter", coFloat);
    def->label = L("Pinhead front diameter");
    def->category = OptionCategory::support;
    def->tooltip = L("Diameter of the pointing side of the head");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.4));

    def = this->add(prefix + "support_head_penetration", coFloat);
    def->label = L("Head penetration");
    def->category = OptionCategory::support;
    def->tooltip = L("How much the pinhead has to penetrate the model surface");
    def->sidetext = L("mm");
    def->mode = comAdvancedE | comPrusa;
    def->min = 0;
    def->set_default_value(new ConfigOptionFloat(0.2));

    def = this->add(prefix + "support_head_width", coFloat);
    def->label = L("Pinhead width");
    def->category = OptionCategory::support;
    def->tooltip = L("Width from the back sphere center to the front sphere center");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 20;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add(prefix + "support_pillar_diameter", coFloat);
    def->label = L("Pillar diameter");
    def->category = OptionCategory::support;
    def->tooltip = L("Diameter in mm of the support pillars");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 15;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add(prefix + "support_small_pillar_diameter_percent", coPercent);
    def->label = L("Small pillar diameter percent");
    def->category = OptionCategory::support;
    def->tooltip = L("The percentage of smaller pillars compared to the normal pillar diameter "
                      "which are used in problematic areas where a normal pilla cannot fit.");
    def->sidetext = L("%");
    def->min = 1;
    def->max = 100;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionPercent(50));

    def = this->add(prefix + "support_max_bridges_on_pillar", coInt);
    def->label = L("Max bridges on a pillar");
    def->tooltip = L(
        "Maximum number of bridges that can be placed on a pillar. Bridges "
        "hold support point pinheads and connect to pillars as small branches.");
    def->min = 0;
    def->max = 50;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionInt(prefix == "branching" ? 2 : 3));

    def = this->add(prefix + "support_max_weight_on_model", coFloat);
    def->label = L("Max weight on model");
    def->category = OptionCategory::support;
    def->tooltip  = L(
        "Maximum weight of sub-trees that terminate on the model instead of the print bed. The weight is the sum of the lenghts of all "
        "branches emanating from the endpoint.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(10.));

    def = this->add(prefix + "support_pillar_connection_mode", coEnum);
    def->label = L("Pillar connection mode");
    def->tooltip = L("Controls the bridge type between two neighboring pillars."
                            " Can be zig-zag, cross (double zig-zag) or dynamic which"
                            " will automatically switch between the first two depending"
                            " on the distance of the two pillars.");
    def->set_enum<SLAPillarConnectionMode>(
        ConfigOptionEnum<SLAPillarConnectionMode>::get_enum_names(),
        { L("Zig-Zag"), L("Cross"), L("Dynamic") });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum(SLAPillarConnectionMode::dynamic));

    def = this->add(prefix + "support_buildplate_only", coBool);
    def->label = L("Support on build plate only");
    def->category = OptionCategory::support;
    def->tooltip = L("Only create support if it lies on a build plate. Don't create support on a print.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add(prefix + "support_pillar_widening_factor", coFloat);
    def->label = L("Pillar widening factor");
    def->category = OptionCategory::support;

    def->tooltip  = 
        L("Merging bridges or pillars into another pillars can "
        "increase the radius. Zero means no increase, one means "
        "full increase. The exact amount of increase is unspecified and can "
        "change in the future.");

    def->min = 0;
    def->max = 1;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.5));

    def = this->add(prefix + "support_base_diameter", coFloat);
    def->label = L("Support base diameter");
    def->category = OptionCategory::support;
    def->tooltip = L("Diameter in mm of the pillar base");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 30;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(4.0));

    def = this->add(prefix + "support_base_height", coFloat);
    def->label = L("Support base height");
    def->category = OptionCategory::support;
    def->tooltip = L("The height of the pillar base cone");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add(prefix + "support_base_safety_distance", coFloat);
    def->label = L("Support base safety distance");
    def->category = OptionCategory::support;
    def->tooltip  = L(
        "The minimum distance of the pillar base from the model in mm. "
        "Makes sense in zero elevation mode where a gap according "
        "to this parameter is inserted between the model and the pad.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 10;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1));

    def = this->add(prefix + "support_critical_angle", coFloat);
    def->label = L("Critical angle");
    def->category = OptionCategory::support;
    def->tooltip = L("The default angle for connecting support sticks and junctions.");
    def->sidetext = L("°");
                    def->min = 0;
    def->max = 90;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(45));

    def = this->add(prefix + "support_max_bridge_length", coFloat);
    def->label = L("Max bridge length");
    def->category = OptionCategory::support;
    def->tooltip = L("The max length of a bridge");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;

    double default_val = 15.0;
    if (prefix == "branching")
        default_val = 5.0;

    def->set_default_value(new ConfigOptionFloat(default_val));

    def = this->add(prefix + "support_max_pillar_link_distance", coFloat);
    def->label = L("Max pillar linking distance");
    def->category = OptionCategory::support;
    def->tooltip = L("The max distance of two pillars to get linked with each other."
                               " A zero value will prohibit pillar cascading.");
    def->sidetext = L("mm");
    def->min = 0;   // 0 means no linking
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(10.0));

    def = this->add(prefix + "support_object_elevation", coFloat);
    def->label = L("Object elevation");
    def->category = OptionCategory::support;
    def->tooltip = L("How much the supports should lift up the supported object. "
                      "If \"Pad around object\" is enabled, this value is ignored.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 150; // This is the max height of print on SL1
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(5.0));
}

void PrintConfigDef::init_sla_params()
{
    m_material_overrides_option_keys = {
        "branchingsupport_head_front_diameter",
        "branchingsupport_head_penetration",
        "branchingsupport_head_width",
        "branchingsupport_pillar_diameter",
        "first_layer_size_compensation",
        "relative_correction_x",
        "relative_correction_y",
        "relative_correction_z",
        "support_head_front_diameter",
        "support_head_penetration",
        "support_head_width",
        "support_pillar_diameter",
        "support_points_density_relative",
    };

    ConfigOptionDef* def;

    // SLA Printer settings

    def = this->add("display_width", coFloat);
    def->label = L("Display width");
    def->tooltip = L("Width of the display");
    def->min = 1;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(120.));

    def = this->add("display_height", coFloat);
    def->label = L("Display height");
    def->tooltip = L("Height of the display");
    def->min = 1;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(68.));

    def = this->add("display_pixels_x", coInt);
    def->full_label = L("Number of pixels in");
    def->label = L("X");
    def->tooltip = L("Number of pixels in X");
    def->min = 100;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(2560));

    def = this->add("display_pixels_y", coInt);
    def->label = L("Y");
    def->tooltip = L("Number of pixels in Y");
    def->min = 100;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(1440));

    def = this->add("display_mirror_x", coBool);
    def->full_label = L("Display horizontal mirroring");
    def->label = L("Mirror horizontally");
    def->tooltip = L("Enable horizontal mirroring of output images");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("display_mirror_y", coBool);
    def->full_label = L("Display vertical mirroring");
    def->label = L("Mirror vertically");
    def->tooltip = L("Enable vertical mirroring of output images");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("display_orientation", coEnum);
    def->label = L("Display orientation");
    def->tooltip = L("Set the actual LCD display orientation inside the SLA printer."
                     " Portrait mode will flip the meaning of display width and height parameters"
                     " and the output images will be rotated by 90 degrees.");
    def->set_enum<SLADisplayOrientation>({
        { "landscape",  L("Landscape") },
        { "portrait",   L("Portrait") }
    });
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionEnum<SLADisplayOrientation>(sladoPortrait));

    def = this->add("fast_tilt_time", coFloat);
    def->label = L("Fast");
    def->full_label = L("Fast tilt");
    def->tooltip = L("Time of the fast tilt");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(5.));

    def = this->add("slow_tilt_time", coFloat);
    def->label = L("Slow");
    def->full_label = L("Slow tilt");
    def->tooltip = L("Time of the slow tilt");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(8.));

    def = this->add("high_viscosity_tilt_time", coFloat);
    def->label = L("High viscosity");
    def->full_label = L("Tilt for high viscosity resin");
    def->tooltip = L("Time of the super slow tilt");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert;
    def->set_default_value(new ConfigOptionFloat(10.));

    def = this->add("area_fill", coFloat);
    def->label = L("Area fill");
    def->tooltip = L("The percentage of the bed area. \nIf the print area exceeds the specified value, \nthen a slow tilt will be used, otherwise - a fast tilt");
    def->sidetext = L("%");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(50.));

    def = this->add("relative_correction", coFloats);
    def->label = L("Printer scaling correction");
    def->full_label = L("Printer scaling correction");
    def->tooltip  = L("Printer scaling correction");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloats( { 1., 1.} ));

    def = this->add("relative_correction_x", coFloat);
    def->label = L("Printer scaling correction in X axis");
    def->full_label = L("Printer scaling X axis correction");
    def->tooltip = L("Printer scaling correction in X axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("relative_correction_y", coFloat);
    def->label = L("Printer scaling correction in Y axis");
    def->full_label = L("Printer scaling Y axis correction");
    def->tooltip = L("Printer scaling correction in Y axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("relative_correction_z", coFloat);
    def->label = L("Printer scaling correction in Z axis");
    def->full_label = L("Printer scaling Z axis correction");
    def->tooltip = L("Printer scaling correction in Z axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("absolute_correction", coFloat);
    def->label = L("Printer absolute correction");
    def->full_label = L("Printer absolute correction");
    def->tooltip  = L("Will inflate or deflate the sliced 2D polygons according "
                      "to the sign of the correction.");
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.0));
    
    def = this->add("elephant_foot_min_width", coFloat);
    def->label = L("minimum width");
    def->category = OptionCategory::slicing;
    def->tooltip = L("Minimum width of features to maintain when doing the first layer compensation.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.2));

    def = this->add("gamma_correction", coFloat);
    def->label = L("Printer gamma correction");
    def->full_label = L("Printer gamma correction");
    def->tooltip  = L("This will apply a gamma correction to the rasterized 2D "
                      "polygons. A gamma value of zero means thresholding with "
                      "the threshold in the middle. This behaviour eliminates "
                      "antialiasing without losing holes in polygons.");
    def->min = 0;
    def->max = 1;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));


    // SLA Material settings.

    def = this->add("material_colour", coString);
    def->label = L("Color");
    def->tooltip = L("This is only used in the Slic3r interface as a visual help.");
    def->gui_type = ConfigOptionDef::GUIType::color;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionString("#29B2B2"));

    def = this->add("material_type", coString);
    def->label = L("SLA material type");
    def->tooltip = L("SLA material type");
    def->gui_flags = "show_value";
    def->set_enum_values(ConfigOptionDef::GUIType::select_open,
        { "Tough", "Flexible", "Casting", "Dental", "Heat-resistant" });
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionString("Tough"));

    def = this->add("initial_layer_height", coFloat);
    def->label = L("Initial layer height");
    def->tooltip = L("Initial layer height");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.3));

    def = this->add("idle_temperature", coInts);
    def->label = L("Idle temperature");
    def->tooltip = L("Nozzle temperature when the tool is currently not used in multi-tool setups."
                     "\nThis is only used when 'Ooze prevention' is active in Print Settings.");
    def->sidetext = L("°C");
    def->min = 0;
    def->max = max_temp;
    def->can_be_disabled = true;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(disable_defaultoption(new ConfigOptionInts{30}));

    def = this->add("bottle_volume", coFloat);
    def->label = L("Bottle volume");
    def->tooltip = L("Bottle volume");
    def->sidetext = L("ml");
    def->min = 50;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1000.0));

    def = this->add("bottle_weight", coFloat);
    def->label = L("Bottle weight");
    def->tooltip = L("Bottle weight");
    def->sidetext = L("kg");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add("material_density", coFloat);
    def->label = L("Density");
    def->tooltip = L("Density");
    def->sidetext = L("g/ml");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add("bottle_cost", coFloat);
    def->label = L("Cost");
    def->tooltip = L("Cost");
    def->sidetext = L("money/bottle");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.0));

    def = this->add("faded_layers", coInt);
    def->label = L("Faded layers");
    def->tooltip = L("Number of the layers needed for the exposure time fade from initial exposure time to the exposure time");
    def->min = 3;
    def->max = 20;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionInt(10));

    def = this->add("min_exposure_time", coFloat);
    def->label = L("Minimum exposure time");
    def->tooltip = L("Minimum exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("max_exposure_time", coFloat);
    def->label = L("Maximum exposure time");
    def->tooltip = L("Maximum exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(100));

    def = this->add("exposure_time", coFloat);
    def->label = L("Exposure time");
    def->tooltip = L("Exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(10));

    def = this->add("min_initial_exposure_time", coFloat);
    def->label = L("Minimum initial exposure time");
    def->tooltip = L("Minimum initial exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("max_initial_exposure_time", coFloat);
    def->label = L("Maximum initial exposure time");
    def->tooltip = L("Maximum initial exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(150));

    def = this->add("initial_exposure_time", coFloat);
    def->label = L("Initial exposure time");
    def->tooltip = L("Initial exposure time");
    def->sidetext = L("s");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(15));

    def = this->add("material_correction", coFloats);
    def->label = L("Correction for expansion");
    def->tooltip  = L("Correction for expansion");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloats({ 1., 1., 1. }));

    def = this->add("material_correction_x", coFloat);
    def->full_label = L("Correction for expansion in X axis");
    def->tooltip = L("Correction for expansion in X axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("material_correction_y", coFloat);
    def->full_label = L("Correction for expansion in Y axis");
    def->tooltip = L("Correction for expansion in Y axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("material_correction_z", coFloat);
    def->full_label = L("Correction for expansion in Z axis");
    def->tooltip = L("Correction for expansion in Z axis");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("material_notes", coString);
    def->label = L("SLA print material notes");
    def->tooltip = L("You can put your notes regarding the SLA print material here.");
    def->multiline = true;
    def->full_width = true;
    def->height = 13;
    // TODO currently notes are the only way to pass data
    // for non-PrusaResearch printers. We therefore need to always show them 
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionString(""));

    def = this->add("material_vendor", coString);
    def->set_default_value(new ConfigOptionString(L("(Unknown)")));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("default_sla_material_profile", coString);
    def->label = L("Default SLA material profile");
    def->tooltip = L("Default print profile associated with the current printer profile. "
                   "On selection of the current printer profile, this print profile will be activated.");
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("sla_material_settings_id", coString);
    def->set_default_value(new ConfigOptionString(""));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("sla_material_settings_modified", coBool);
    def->set_default_value(new ConfigOptionBool(false));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("default_sla_print_profile", coString);
    def->label = L("Default SLA material profile");
    def->tooltip = L("Default print profile associated with the current printer profile. "
                   "On selection of the current printer profile, this print profile will be activated.");
    def->set_default_value(new ConfigOptionString());
    def->cli = ConfigOptionDef::nocli;

    def = this->add("sla_print_settings_id", coBool);
    def->set_default_value(new ConfigOptionBool(false));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("sla_print_settings_modified", coString);
    def->set_default_value(new ConfigOptionString(""));
    def->cli = ConfigOptionDef::nocli;

    def = this->add("supports_enable", coBool);
    def->label = L("Generate supports");
    def->category = OptionCategory::support;
    def->tooltip = L("Generate supports for the models");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("support_tree_type", coEnum);
    def->label = L("Support tree type");
    def->tooltip = L("Support tree building strategy");
    def->set_enum<sla::SupportTreeType>(
        ConfigOptionEnum<sla::SupportTreeType>::get_enum_names(),
        { L("Default"),
    // TRN One of the "Support tree type"s on SLAPrintSettings : Supports
            L("Branching (experimental)") });
    // TODO: def->enum_def->labels[2] = L("Organic");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionEnum(sla::SupportTreeType::Default));

    init_sla_support_params("");
    init_sla_support_params("branching");

    def = this->add("support_enforcers_only", coBool);
    def->label = L("Support only in enforced regions");
    def->category = OptionCategory::support;
    def->tooltip = L("Only create support if it lies in a support enforcer.");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("support_points_density_relative", coInt);
    def->label = L("Support points density");
    def->category = OptionCategory::support;
    def->tooltip = L("This is a relative measure of support points density.");
    def->sidetext = L("%");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionInt(100));

    def = this->add("support_points_minimal_distance", coFloat);
    def->label = L("Minimal distance of the support points");
    def->category = OptionCategory::support;
    def->tooltip = L("No support points will be placed closer than this threshold.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.));

    def = this->add("pad_enable", coBool);
    def->label = L("Use pad");
    def->category = OptionCategory::pad;
    def->tooltip = L("Add a pad underneath the supported model");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("pad_wall_thickness", coFloat);
    def->label = L("Pad wall thickness");
    def->category = OptionCategory::pad;
     def->tooltip = L("The thickness of the pad and its optional cavity walls.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 30;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(2.0));

    def = this->add("pad_wall_height", coFloat);
    def->label = L("Pad wall height");
    def->tooltip = L("Defines the pad cavity depth. Set zero to disable the cavity. "
                     "Be careful when enabling this feature, as some resins may "
                     "produce an extreme suction effect inside the cavity, "
                     "which makes peeling the print off the vat foil difficult.");
    def->category = OptionCategory::pad;
//     def->tooltip = L("");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 30;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.));
    
    def = this->add("pad_brim_size", coFloat);
    def->label = L("Pad brim size");
    def->tooltip = L("How far should the pad extend around the contained geometry");
    def->category = OptionCategory::pad;
    //     def->tooltip = L("");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 30;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1.6));

    def = this->add("pad_max_merge_distance", coFloat);
    def->label = L("Max merge distance");
    def->category = OptionCategory::pad;
     def->tooltip = L("Some objects can get along with a few smaller pads "
                      "instead of a single big one. This parameter defines "
                      "how far the center of two smaller pads should be. If they"
                      "are closer, they will get merged into one pad.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(50.0));

    // This is disabled on the UI. I hope it will never be enabled.
//    def = this->add("pad_edge_radius", coFloat);
//    def->label = L("Pad edge radius");
//    def->category = OptionCategory::pad;
////     def->tooltip = L("");
//    def->sidetext = L("mm");
//    def->min = 0;
//    def->mode = comAdvancedE | comPrusa;
//    def->set_default_value(new ConfigOptionFloat(1.0));

    def = this->add("pad_wall_slope", coFloat);
    def->label = L("Pad wall slope");
    def->category = OptionCategory::pad;
    def->tooltip = L("The slope of the pad wall relative to the bed plane. "
                     "90 degrees means straight walls.");
    def->sidetext = L("°");
    def->min = 45;
    def->max = 90;
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(90.0));

    def = this->add("pad_around_object", coBool);
    def->label = L("Pad around object");
    def->category = OptionCategory::pad;
    def->tooltip = L("Create pad around object and ignore the support elevation");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));
    
    def = this->add("pad_around_object_everywhere", coBool);
    def->label = L("Pad around object everywhere");
    def->category = OptionCategory::pad;
    def->tooltip = L("Force pad around object everywhere");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("pad_object_gap", coFloat);
    def->label = L("Pad object gap");
    def->category = OptionCategory::pad;
    def->tooltip  = L("The gap between the object bottom and the generated "
                      "pad in zero elevation mode.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 10;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(1));

    def = this->add("pad_object_connector_stride", coFloat);
    def->label = L("Pad object connector stride");
    def->category = OptionCategory::pad;
    def->tooltip = L("Distance between two connector sticks which connect the object and the generated pad.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(10));

    def = this->add("pad_object_connector_width", coFloat);
    def->label = L("Pad object connector width");
    def->category = OptionCategory::pad;
    def->tooltip  = L("Width of the connector sticks which connect the object and the generated pad.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.5));

    def = this->add("pad_object_connector_penetration", coFloat);
    def->label = L("Pad object connector penetration");
    def->category = OptionCategory::pad;
    def->tooltip  = L(
        "How much should the tiny connectors penetrate into the model body.");
    def->sidetext = L("mm");
    def->min = 0;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.3));
    
    def = this->add("hollowing_enable", coBool);
    def->label = L("Enable hollowing");
    def->category = OptionCategory::hollowing;
    def->tooltip = L("Hollow out a model to have an empty interior");
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionBool(false));
    
    def = this->add("hollowing_min_thickness", coFloat);
    def->label = L("Wall thickness");
    def->category = OptionCategory::hollowing;
    def->tooltip  = L("Minimum wall thickness of a hollowed model.");
    def->sidetext = L("mm");
    def->min = 1;
    def->max = 10;
    def->mode = comSimpleAE | comPrusa;
    def->set_default_value(new ConfigOptionFloat(3.));
    
    def = this->add("hollowing_quality", coFloat);
    def->label = L("Accuracy");
    def->category = OptionCategory::hollowing;
    def->tooltip  = L("Performance vs accuracy of calculation. Lower values may produce unwanted artifacts.");
    def->min = 0;
    def->max = 1;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.5));
    
    def = this->add("hollowing_closing_distance", coFloat);
    def->label = L("Closing distance");
    def->category = OptionCategory::hollowing;
    def->tooltip  = L(
        "Hollowing is done in two steps: first, an imaginary interior is "
        "calculated deeper (offset plus the closing distance) in the object and "
        "then it's inflated back to the specified offset. A greater closing "
        "distance makes the interior more rounded. At zero, the interior will "
        "resemble the exterior the most.");
    def->sidetext = L("mm");
    def->min = 0;
    def->max = 10;
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(2.0));

    def = this->add("material_print_speed", coEnum);
    def->label = L("Print speed");
    def->tooltip = L(
        "A slower printing profile might be necessary when using materials with higher viscosity "
        "or with some hollowed parts. It slows down the tilt movement and adds a delay before exposure.");
    def->set_enum<SLAMaterialSpeed>({
        { "slow",           L("Slow") },
        { "fast",           L("Fast") },
        { "high_viscosity", L("High viscosity") }
    });
    def->mode = comAdvancedE | comPrusa;
    def->set_default_value(new ConfigOptionEnum<SLAMaterialSpeed>(slamsFast));

    //def = this->add("sla_archive_format", coString);
    //def->label = L("Format of the output SLA archive");
    //def->mode = comAdvancedE | comPrusa;
    //def->set_default_value(new ConfigOptionString("SL1"));

    def = this->add("output_format", coEnum);
    def->label = L("Output Format");
    def->tooltip = L("Select the output format for this printer.");
    def->set_enum<OutputFormat>({
        {"SL1", L("Prusa SL1")},
        {"SL1_SVG", L("Prusa SL1 with SVG")},
        {"mCWS", L("Masked CWS")},
        {"AnyMono", L("Anycubic Mono")},
        {"AnyMonoX", L("Anycubic Mono X")},
        {"AnyMonoSE", L("Anycubic Mono SE")},
    });
    def->mode = comAdvancedE | comSuSi; // output_format should be preconfigured in profiles;
    def->set_default_value(new ConfigOptionEnum<OutputFormat>(ofSL1));

    def = this->add("sla_output_precision", coFloat);
    def->label = L("SLA output precision");
    def->tooltip = L("Minimum resolution in nanometers");
    def->sidetext = L("mm");
    def->min = float(SCALING_FACTOR);
    def->mode = comExpert | comPrusa;
    def->set_default_value(new ConfigOptionFloat(0.001));

    // Declare retract values for material profile, overriding the print and printer profiles.
    for (const std::string &opt_key : m_material_overrides_option_keys) {
        auto it_opt = options.find(opt_key);
        assert(it_opt != options.end());
        def = this->add(std::string("material_ow_") + opt_key, it_opt->second.type);
        def->can_be_disabled = true;
        def->is_optional = true;
        def->label = it_opt->second.label;
        def->full_label = it_opt->second.full_label;
        def->tooltip = it_opt->second.tooltip;
        def->sidetext = it_opt->second.sidetext;
        def->min  = it_opt->second.min;
        def->max  = it_opt->second.max;
        def->mode = it_opt->second.mode;
        ConfigOption *default_opt = it_opt->second.default_value->clone();
        default_opt->set_can_be_disabled(true);
        def->set_default_value(default_opt);
        assert(!def->default_value->is_enabled());
    }
}

// Ignore the following obsolete configuration keys:
static std::set<std::string> PrintConfigDef_ignore = {
    "clip_multipart_objects",
    "duplicate_x", "duplicate_y", "gcode_arcs", "multiply_x", "multiply_y",
    "support_material_tool", "acceleration", "adjust_overhang_flow",
    "standby_temperature", "scale", "rotate", "duplicate", "duplicate_grid",
    "start_perimeters_at_concave_points", "start_perimeters_at_non_overhang", "randomize_start",
    "seal_position", "vibration_limit", "bed_size",
    "print_center", "g0", "threads", "pressure_advance", "wipe_tower_per_color_wipe",
    "cooling", "serial_port", "serial_speed",
    "exact_last_layer_height",
    // Introduced in some PrusaSlicer 2.3.1 alpha, later renamed or removed.
    "fuzzy_skin_perimeter_mode", "fuzzy_skin_shape",
    //replaced by brim_per_object, but can't translate the value as the old one is only used for complete_objects (and the default are differents).
    "complete_objects_one_brim",
    // Introduced in PrusaSlicer 2.3.0-alpha2, later replaced by automatic calculation based on extrusion width.
    "wall_add_middle_threshold", "wall_split_middle_threshold",
    // Replaced by new concentric ensuring in 2.6.0-alpha5
    "ensure_vertical_shell_thickness",
    // Disabled in 2.6.0-alpha6, this option is problematic
//    "infill_only_where_needed", <- ignore only if deactivated
    "gcode_binary", // Introduced in 2.7.0-alpha1, removed in 2.7.1 (replaced by binary_gcode).
    "gcode_resolution", // now in printer config.
    "enable_dynamic_fan_speeds", "overhang_fan_speed_0","overhang_fan_speed_1","overhang_fan_speed_2","overhang_fan_speed_3", // converted in composite_legacy
    "enable_dynamic_overhang_speeds", "overhang_speed_0", "overhang_speed_1", "overhang_speed_2", "overhang_speed_3", // converted in composite_legacy
    "travel_max_lift", "filament_travel_max_lift" // removed, using retract_lift also for rampping lift instead.
};

void PrintConfigDef::handle_legacy(t_config_option_key &opt_key, std::string &value, bool remove_unkown_keys)
{
    // handle legacy options (other than aliases)
    if (opt_key == "extrusion_width_ratio" || opt_key == "bottom_layer_speed_ratio"
        || opt_key == "first_layer_height_ratio") {
        boost::replace_first(opt_key, "_ratio", "");
        if (opt_key == "bottom_layer_speed") opt_key = "first_layer_speed";
        try {
            float v = boost::lexical_cast<float>(value);
            if (v != 0)
                value = boost::lexical_cast<std::string>(v*100) + "%";
        } catch (boost::bad_lexical_cast &) {
            value = "0";
        }
    }
    if ("infill_only_where_needed" == opt_key && "0" == value) {
        opt_key = "";
        value = "";
    }
    if (opt_key == "gcode_flavor") {
        if (value == "makerbot")
            value = "makerware";
        else if (value == "marlinfirmware")
            // the "new" marlin firmware flavor used to be called "marlinfirmware" for some time during PrusaSlicer 2.4.0-alpha development.
            value = "marlin2";
    } else if (opt_key == "host_type" && value == "mainsail") {
        // the "mainsail" key (introduced in 2.6.0-alpha6) was renamed to "moonraker" (in 2.6.0-rc1).
        value = "moonraker";
    } else if (opt_key == "fill_density" && value.find("%") == std::string::npos) {
        try {
            // fill_density was turned into a percent value
            float v = boost::lexical_cast<float>(value);
            value = boost::lexical_cast<std::string>(v*100) + "%";
        } catch (boost::bad_lexical_cast &) {}
    }
    if (opt_key == "randomize_start" && value == "1") {
        opt_key = "seam_position";
        value = "random";
    }
    if (opt_key == "bed_size" && !value.empty()) {
        opt_key = "bed_shape";
        ConfigOptionPoint p;
        p.deserialize(value, ForwardCompatibilitySubstitutionRule::Disable);
        std::ostringstream oss;
        oss << "0x0," << p.value(0) << "x0," << p.value(0) << "x" << p.value(1) << ",0x" << p.value(1);
        value = oss.str();
    }
    //if ((opt_key == "perimeter_acceleration" && value == "25")
    //    || (opt_key == "infill_acceleration" && value == "50")) {
    //    /*  For historical reasons, the world's full of configs having these very low values;
    //        to avoid unexpected behavior we need to ignore them. Banning these two hard-coded
    //        values is a dirty hack and will need to be removed sometime in the future, but it
    //        will avoid lots of complaints for now. */
    //    value = "0";
    //} // i think it's time now.
    //TODO: change defautl from 0 (no forbidden) to 100% for all speed & acceleration
    //if ("0" == value && ("infill_acceleration" == opt_key || "solid_infill_acceleration" == opt_key ||
    //    "top_solid_infill_acceleration" == opt_key || "bridge_acceleration" == opt_key ||
    //    "default_acceleration" == opt_key || "perimeter_acceleration" == opt_key || "overhangs_speed" == opt_key ||
    //    "ironing_speed" == opt_key || "perimeter_speed" == opt_key || "infill_speed" == opt_key ||
    //    "bridge_speed" == opt_key || "support_material_speed" == opt_key || "max_print_speed" == opt_key)) {
    //    // 100% is the same, and easier to understand.
    //    value = "100%";
    //}
    if (opt_key == "support_material_pattern" && value == "pillars") {
        // Slic3r PE does not support the pillars. They never worked well.
        value = "rectilinear";
    }
    if (opt_key == "skirt_height" && value == "-1") {
        // PrusaSlicer no more accepts skirt_height == -1 to print a draft shield to the top of the highest object.
        // A new "draft_shield" enum config value is used instead.
        opt_key = "draft_shield";
        value = "enabled";
    } else if (opt_key == "draft_shield" && (value == "1" || value == "0")) {
        // draft_shield used to be a bool, it was turned into an enum in PrusaSlicer 2.4.0.
        value = value == "1" ? "enabled" : "disabled";
    } else if (("label_printed_objects" == opt_key || "gcode_label_objects" == opt_key) && (value == "1" || value == "0")) {
        // gcode_label_objects used to be a bool (the behavior was nothing or "octoprint"), it is
        // and enum since PrusaSlicer 2.6.2.
        value = value == "1" ? "octoprint" : "disabled";
    } else if (opt_key == "octoprint_host") {
        opt_key = "print_host";
    }
    if (opt_key == "octoprint_cafile") {
        opt_key = "printhost_cafile";
    }
    if (opt_key == "octoprint_apikey") {
        opt_key = "printhost_apikey";
    }
    if (opt_key == "elefant_foot_compensation") {
        opt_key = "first_layer_size_compensation";
        float v = boost::lexical_cast<float>(value);
        if (v > 0)
            value = boost::lexical_cast<std::string>(-v);
    }
    if ("elefant_foot_min_width" == opt_key) {
        opt_key = "elephant_foot_min_width";
    }
    if (opt_key == "thumbnails") {
        if (value.empty())
            value = "0x0,0x0";
    }
    if (opt_key == "z_steps_per_mm") {
        opt_key = "z_step";
        float v = boost::lexical_cast<float>(value);
        if(v > 0)
            value = boost::lexical_cast<std::string>(1/v);
    }
    if (opt_key == "infill_not_connected") {
        opt_key = "infill_connection";
        if (value == "1")
            value = "notconnected";
        else
            value = "connected";
    }
    if (opt_key == "preset_name") {
        opt_key = "preset_names";
    }
    if (opt_key == "seam_travel") {
        if (value == "1") {
            opt_key = "seam_travel_cost";
            value = "200%";
        } else {
            opt_key = "";
        }
    }
    if (opt_key == "seam_position") {
        if (value == "hidden") {
            value = "cost";
        } else if ("near" == value || "nearest" == value) {
            value = "cost";
            //FIXME can we change the cost?
        }
    }
    if (opt_key == "perimeter_loop_seam") {
        if (value == "hidden") {
            value = "nearest";
        }
    }
    if (opt_key == "overhangs") {
        opt_key = "overhangs_width_speed";
        if (value == "1")
            value = "50%";
        else
            value = "0";
    }
    if (opt_key == "print_machine_envelope") {
        opt_key = "machine_limits_usage";
        if (value == "1")
            value = "emit_to_gcode";
        else
            value = "time_estimate_only";
    }
    if (opt_key == "retract_lift_not_last_layer") {
        opt_key = "retract_lift_top";
        if (value == "1")
            value = "Not on top";
        else
            value = "All surfaces";
    }
    if ("gcode_precision_e" == opt_key) {
        if (value.find(",") != std::string::npos)
            value = value.substr(0, value.find(","));
        try {
            int val = boost::lexical_cast<int>(value);
            if (val > 0)
                value = boost::lexical_cast<std::string>(val);
        }
        catch (boost::bad_lexical_cast&) {
            value = "5";
        }
    }
    if ("first_layer_min_speed" == opt_key && value.back() == '%')
        value = value.substr(0, value.length() - 1); //no percent.
    if (!value.empty() && value.back() != '%' && std::set<std::string>{"bridge_flow_ratio", "bridge_flow_ratio", "over_bridge_flow_ratio", "fill_top_flow_ratio", "first_layer_flow_ratio"}.count(opt_key) > 0 ) {
        //need percent
        try {
            float val = boost::lexical_cast<float>(value);
            if (val < 2)
                value = boost::lexical_cast<std::string>(val*100) + "%";
            else
                value = "100%";
        }
        catch (boost::bad_lexical_cast&) {
            value = "100%";
        }
    }
    if("thick_bridges" == opt_key) {
        opt_key = "bridge_type";
        if (value == "1")
            value = "nozzle";
        else
            value = "flow";
    }
    if ("sla_archive_format" == opt_key) {
        opt_key = "output_format";
    }

    // In PrusaSlicer 2.3.0-alpha0 the "monotonic" infill was introduced, which was later renamed to "monotonous".
    if (value == "monotonous" && (opt_key == "top_fill_pattern" || opt_key == "bottom_fill_pattern" || opt_key == "fill_pattern"
            || opt_key == "solid_fill_pattern" || opt_key == "bridge_fill_pattern" || opt_key == "support_material_interface_pattern")) {
        value = "monotonic";
    }
    // some changes has occurs between rectilineargapfill and monotonicgapfill. Set them at the right value for each type
    if (value == "rectilineargapfill" && (opt_key == "top_fill_pattern" || opt_key == "bottom_fill_pattern") )
        value = "monotonicgapfill";
    if (opt_key == "fill_pattern" || opt_key == "support_material_interface_pattern" || opt_key == "support_material_top_interface_pattern" || opt_key == "support_material_bottom_interface_pattern") {
        if (value == "rectilineargapfill") {
            value = "rectilinear";
        } else if (value == "monotonicgapfill") {
            value = "monotonic";
        }
    }
    //in ps 2.4, the raft_first_layer_density is now more powerful than the support_material_solid_first_layer, also it always does the perimeter.
    if ("support_material_solid_first_layer" == opt_key) {
        opt_key = "raft_first_layer_density";
        value = "100";
    }
    if (boost::starts_with(opt_key, "thin_perimeters") && value == "1") {
        value = "100%";
    }

    //prusa
    if ("gcode_flavor" == opt_key) {
        if ("reprap" == value)
            value = "sprinter";
    }

    if (PrintConfigDef_ignore.find(opt_key) != PrintConfigDef_ignore.end()) {
        opt_key = "";
        return;
    }

    if ("fan_always_on" == opt_key) {
        if (value != "1") {
            //min_fan_speed is already converted to default_fan_speed, just has to deactivate it if not always_on
            opt_key = "default_fan_speed"; // note: maybe this doesn't works, as default_fan_speed can also get its value from min_fan_speed
            value = "0";
        } else {
            opt_key = "";
        }
        return;
    }
    if ("arc_fitting" == opt_key) {
        if (value == "1")
            value = "bambu";
        else if (value == "0")
            value = "disabled";
    }

    if (!print_config_def.has(opt_key)) {
        //check the aliases
        for(const auto& entry : print_config_def.options) {
            for (const std::string& alias : entry.second.aliases) {
                if (alias == opt_key) {
                    // translate
                    opt_key = entry.first;
                    goto use_alias;
                }
            }
        }
        if (remove_unkown_keys) {
            opt_key = "";
        }
    use_alias:;
        if (!print_config_def.has(opt_key)) {
            return;
        }
    }

    //fan speed: activate disable.
    if (opt_key.find("_fan_speed") != std::string::npos) {
        if ("max_fan_speed" != opt_key && "filament_toolchange_part_fan_speed" != opt_key && "min_fan_speed" != opt_key
            && "overhangs_dynamic_fan_speed" != opt_key && "enable_dynamic_fan_speeds" != opt_key && opt_key.find("overhang_fan_speed_") == std::string::npos) {
            assert(print_config_def.get(opt_key) && print_config_def.get(opt_key)->type == coInts);
            //if vector, split it.
            ConfigOptionInts opt_decoder;
            opt_decoder.set_can_be_disabled();
            opt_decoder.deserialize(value);
            for (size_t idx = 0; idx < opt_decoder.size(); ++idx) {
                if (opt_decoder.is_enabled(idx)) {
                    if (opt_decoder.get_at(idx) < 0) {
                        opt_decoder.set_at(100, idx);
                        opt_decoder.set_enabled(false, idx);
                    } else if (opt_decoder.get_at(idx) <= 1) {
                        // for now, still consider "1" as a "0", to be able to import old config where the 1 means 0
                        // (and 0 was disable).
                        opt_decoder.set_at(0, idx);
                    }
                }
            }
            value = opt_decoder.serialize();
        }
    }
    if ("max_layer_height" == opt_key && "0" == value) {
        value = "!75%";
    }
    if (value == "-1") {
        if ("overhangs_bridge_threshold" == opt_key) {value = "!0";}
        if ("overhangs_bridge_upper_layers" == opt_key) {value = "!2";}
        if ("perimeters_hole" == opt_key) {value = "!0";}
        if ("support_material_bottom_interface_layers" == opt_key) {value = "!0";}
    }
    //nil-> disabled
    if (value.find("e+") != std::string::npos) {
        const ConfigOptionDef *def = print_config_def.get(opt_key);
        if (def && def->can_be_disabled) {
            ConfigOption *default_opt = def->default_value->clone();
            default_opt->deserialize(value);
            float max_value = std::numeric_limits<int32_t>::max() / 2;
            switch (default_opt->type()) {
            case coInt:
            case coPercent:
            case coFloat:
            case coFloatOrPercent:
            case coInts:
            case coPercents:
            case coFloats:
            case coFloatsOrPercents: {
                for (size_t idx = 0; idx < default_opt->size(); idx++) {
                    if (std::abs(default_opt->get_float(idx)) > std::numeric_limits<int>::max() / 2) {
                        default_opt->set(def->default_value.get(), idx);
                        default_opt->set_enabled(false, idx);
                    }
                }
            }
            break;
            default:;
            }
            value = default_opt->serialize();
            delete default_opt;
        }
    }
    //nil-> disabled
    if (value.find("nil") != std::string::npos) {
        const ConfigOptionDef *def = print_config_def.get(opt_key);
        if (def->type != coString && def->type != coStrings) {
            assert(def && def->can_be_disabled);
            if (def && def->can_be_disabled) {
                ConfigOption *default_opt = def->default_value->clone();
                default_opt->set_enabled(false);
                value = default_opt->serialize();
                delete default_opt;
            }
        }
    }
}

// Called after a config is loaded as a whole.
// Perform composite conversions, for example merging multiple keys into one key.
// Don't convert single options here, implement such conversion in PrintConfigDef::handle_legacy() instead.
void PrintConfigDef::handle_legacy_composite(DynamicPrintConfig &config, std::vector<std::pair<t_config_option_key, std::string>> &opt_deleted)
{
    std::map<t_config_option_key, std::string> useful_items;
    for (auto& opt_pair : opt_deleted) {
        t_config_option_key &opt_key = opt_pair.first;
        std::string &value = opt_pair.second;
        if (opt_key.find("overhang_fan_speed_") != std::string::npos) {
            useful_items[opt_key] = value;
            opt_key = "";
        }
        if ("enable_dynamic_fan_speeds" == opt_key) {
            useful_items[opt_key] = value;
            opt_key = "";
        }
        if (opt_key.find("overhang_speed_") != std::string::npos) {
            useful_items[opt_key] = value;
            opt_key = "";
        }
        if ("enable_dynamic_overhang_speeds" == opt_key) {
            useful_items[opt_key] = value;
            opt_key = "";
        }
    }
    bool old = true;
    if (config.has("print_version")) {
        std::string str_version = config.option<ConfigOptionString>("print_version")->value;
        old = str_version.size() < 4+1+7;
        old = old || str_version.substr(0,4) != "SUSI";
        assert(old || str_version[4] == '_');
        if (!old) {
            std::optional<Semver> version = Semver::parse(str_version.substr(5));
            if (version) {
                if (version->maj() <= 2 && version->min() <= 6) {
                    old = true;
                }
            } else {
                old = true;
            }
        }
    }
    if (old && config.has("bridge_angle") && config.get_float("bridge_angle") == 0 && config.is_enabled("bridge_angle")) {
        config.option("bridge_angle")->set_enabled(false);
    }
    if (old && config.has("overhangs_width_speed") && config.get_float("overhangs_width_speed") == 0 && config.is_enabled("overhangs_width_speed")) {
        config.option("overhangs_width_speed")->set_enabled(false);
    }
    if (old && config.has("overhangs_width") && config.get_float("overhangs_width") == 0 && config.is_enabled("overhangs_width")) {
        config.option("overhangs_width")->set_enabled(false);
    }
    if (useful_items.find("enable_dynamic_overhang_speeds") != useful_items.end()) {
        ConfigOptionBool enable_dynamic_overhang_speeds;
        enable_dynamic_overhang_speeds.deserialize(useful_items["enable_dynamic_overhang_speeds"]);
        std::vector<ConfigOptionFloatOrPercent> values;
        values.resize(4);
        values[0].deserialize(useful_items["overhang_speed_0"]);
        values[1].deserialize(useful_items["overhang_speed_1"]);
        values[2].deserialize(useful_items["overhang_speed_2"]);
        values[3].deserialize(useful_items["overhang_speed_3"]);
        double external_perimeter_speed = config.get_computed_value("external_perimeter_speed");
        double max = external_perimeter_speed;
        double min = external_perimeter_speed;
        for (int x = 0; x < values.size(); ++x) {
            if (values[x].percent) {
                min = std::min(min, values[x].get_abs_value(external_perimeter_speed));
                max = std::max(max, values[x].get_abs_value(external_perimeter_speed));
            } else {
                min = std::min(min, values[x].value);
                max = std::max(max, values[x].value);
            }
        }
        //can't have both (min < external_perimeter_speed) & (max > external_perimeter_speed) at same time.
        if (min < external_perimeter_speed) {
            config.set_key_value("overhangs_speed", new ConfigOptionFloatOrPercent(min, false));
            max = external_perimeter_speed;
        } else if (max > external_perimeter_speed) {
            config.set_key_value("overhangs_speed", new ConfigOptionFloatOrPercent(max, false));
            min = external_perimeter_speed;
        } else {
            assert(min == max && min ==external_perimeter_speed);
        }
        ConfigOptionGraph opt;
        opt.set_can_be_disabled();
        // extract values
        Pointfs graph_curve;
        for (int x = 0; x < values.size(); ++x) {
            double speed = std::clamp(values[x].value, min, max);
            if (values[x].percent) {
                speed = values[x].get_abs_value(external_perimeter_speed);
            }
            double percent = (speed - min) / (max - min);
            if (min == external_perimeter_speed) {
                percent = 1 - percent;
            }
            graph_curve.push_back(Vec2d(x * 25, int(percent * 100)));
        }
        if (min == external_perimeter_speed) {
            graph_curve.push_back(Vec2d(100, 0));
        } else {
            graph_curve.push_back(Vec2d(100, 100));
        }
        opt.value = GraphData(graph_curve);
        opt.set_enabled(enable_dynamic_overhang_speeds.value);
        config.set_key_value("overhangs_dynamic_speed", opt.clone());
    }
    if (useful_items.find("enable_dynamic_fan_speeds") != useful_items.end()) {
        ConfigOptionBools enable_dynamic_fan_speeds;
        enable_dynamic_fan_speeds.deserialize(useful_items["enable_dynamic_fan_speeds"]);
        auto *external_perimeter_fan_speed = config.option<ConfigOptionInts>("external_perimeter_fan_speed");
        auto *perimeter_fan_speed = config.option<ConfigOptionInts>("perimeter_fan_speed");
        auto *default_fan_speed = config.option<ConfigOptionInts>("default_fan_speed");
        std::vector<ConfigOptionFloats> values;
        values.resize(4);
        values[0].deserialize(useful_items["overhang_fan_speed_0"]);
        values[1].deserialize(useful_items["overhang_fan_speed_1"]);
        values[2].deserialize(useful_items["overhang_fan_speed_2"]);
        values[3].deserialize(useful_items["overhang_fan_speed_3"]);
        ConfigOptionGraphs opt;
        opt.set_can_be_disabled();
        std::vector<GraphData> graph_data;
        // while there is a value
        assert(enable_dynamic_fan_speeds.size() == values[0].size());
        assert(values[0].size() == values[1].size());
        assert(values[0].size() == values[2].size());
        assert(values[0].size() == values[3].size());
        for(int idx = 0 ;idx < enable_dynamic_fan_speeds.size(); ++idx) {
            // extract values
            Pointfs graph_curve;
            for (int x = 0; x < values.size(); ++x) {
                graph_curve.push_back(Vec2d(x*25, values[x].get_at(idx)));
            }
            if (external_perimeter_fan_speed && external_perimeter_fan_speed->is_enabled(idx)) {
                graph_curve.push_back(Vec2d(100, external_perimeter_fan_speed->get_at(idx)));
            } else if (perimeter_fan_speed && perimeter_fan_speed->is_enabled(idx)) {
                graph_curve.push_back(Vec2d(100, perimeter_fan_speed->get_at(idx)));
            } else if (default_fan_speed && default_fan_speed->is_enabled(idx)) {
                graph_curve.push_back(Vec2d(100, default_fan_speed->get_at(idx)));
            } else {
                graph_curve.push_back(Vec2d(100, graph_curve.back().y()));
            }
            graph_data.emplace_back(graph_curve);
        }
        // recreate fan speed graph
        opt.set(graph_data);
        for (int idx = 0; idx < enable_dynamic_fan_speeds.size(); ++idx) {
            opt.set_enabled(enable_dynamic_fan_speeds.get_at(idx), idx);
        }
        config.set_key_value("overhangs_dynamic_fan_speed", opt.clone());
    }

    //if (config.has("thumbnails")) {
    //    std::string extention;
    //    if (config.has("thumbnails_format")) {
    //        if (const ConfigOptionDef* opt = config.def()->get("thumbnails_format")) {
    //            auto label = opt->enum_def->enum_to_label(config.option("thumbnails_format")->getInt());
    //            if (label.has_value())
    //                extention = *label;
    //        }
    //    }

    //    std::string thumbnails_str = config.opt_string("thumbnails");
    //    auto [thumbnails_list, errors] = GCodeThumbnails::make_and_check_thumbnail_list(thumbnails_str, extention);

    //    if (errors != enum_bitmask<ThumbnailError>()) {
    //        std::string error_str = "\n" + format("Invalid value provided for parameter %1%: %2%", "thumbnails", thumbnails_str);
    //        error_str += GCodeThumbnails::get_error_string(errors);
    //        throw BadOptionValueException(error_str);
    //    }

    //    if (!thumbnails_list.empty()) {
    //        const auto& extentions = ConfigOptionEnum<GCodeThumbnailsFormat>::get_enum_names();
    //        thumbnails_str.clear();
    //        for (const auto& [ext, size] : thumbnails_list)
    //            thumbnails_str += format("%1%x%2%/%3%, ", size.x(), size.y(), extentions[int(ext)]);
    //        thumbnails_str.resize(thumbnails_str.length() - 2);

    //        config.set_key_value("thumbnails", new ConfigOptionString(thumbnails_str));
    //    }
    //}
}

bool PrintConfigDef::is_defined(t_config_option_key &opt_key) { return print_config_def.has(opt_key); }

// this is for extra things to add / modify from prusa that can't be handled otherwise.
// after handle_legacy
std::map<std::string,std::string> PrintConfigDef::from_prusa(t_config_option_key& opt_key, std::string& value, const DynamicConfig& all_conf) {
    std::map<std::string, std::string> output;
    if ("toolchange_gcode" == opt_key) {
        if (!value.empty() && value.find("T") == std::string::npos) {
            value = "T{next_extruder}\n" + value;
        }
    }
    if ("xy_size_compensation" == opt_key) {
        output["xy_inner_size_compensation"] = value;
    }
    if ("infill_anchor_max" == opt_key) {
        if(value == "0")
            output["infill_connection"] = "notconnected";
    }
    if ("first_layer_speed" == opt_key) {
        output["first_layer_min_speed"] = value;
        output["first_layer_infill_speed"] = value;
    }
    //dep between solid infill & infill accel are inverted for prusa.
    if ("infill_acceleration" == opt_key) {
        if (value == "0" && all_conf.get_float("solid_infill_acceleration") != 0) {
            value = std::to_string(all_conf.get_computed_value("default_acceleration"));
        }
    }
    if ("solid_infill_acceleration" == opt_key) {
        if (value == "0" && all_conf.get_float("infill_acceleration") != 0) {
            value = all_conf.get_float("infill_acceleration");
        }
    }
    if ("brim_type" == opt_key) {
        opt_key = "";
        if ("no_brim" == value) {
            output["brim_width"] = "0";
            output["brim_width_interior"] = "0";
        } else if ("outer_only" == value) {
            output["brim_width_interior"] = "0";
        } else if ("inner_only" == value) {
            output["brim_width_interior"] = all_conf.get_computed_value("brim_width");
            output["brim_width"] = "0";
        } else if ("outer_and_inner" == value) {
            output["brim_width_interior"] = all_conf.get_computed_value("brim_width");
        } else if ("auto_brim" == value) { //orca
            output["brim_width"] = "2";
            output["brim_width_interior"] = "0";
        } else if ("brim_ears" == value) { //orca
            output["brim_ears"] = "1";
        }
    }
    if ("support_material_contact_distance" == 0) {
        output["support_material_contact_distance_type"] = "none";
    }
    if (opt_key == "seam_position") {
        if ("cost" == value ) { // eqauls to "near" == value || "nearest" == value
            output["seam_angle_cost"] = "50%";
            output["seam_travel_cost"] = "50%";
        }
    }
    if ("bridge_type" == opt_key) { // seems like thick_bridge to 0
        if (value == "flow") {
            output["bridge_overlap_min"] = "60%";
            output["bridge_overlap"] = "75%";
        }
    }
    if ("first_layer_height" == opt_key) {
        if (!value.empty() && value.back() == '%') {
            // A first_layer_height isn't a % of layer_height but from nozzle_diameter now!
            // can't really convert right now, so put it at a safe value liek 50%.
            value = "50%";
        }
    }
    if ("resolution" == opt_key && value == "0") {
        value = "0.0125";
    }
    if ("gcode_resolution" == opt_key) {
        output["gcode_min_resolution"] = value;
    }
    if (("brim_width" == opt_key || "brim_width_interior" == opt_key) && all_conf.option("brim_separation") ) {
        // add brim_separation to brim_width & brim_width_interior
        float val = boost::lexical_cast<float>(value);
        if (val > 0) {
            val += all_conf.option("brim_separation")->get_float();
            value = boost::lexical_cast<std::string>(val);
        }
    }
    if ("fill_pattern" == opt_key && "alignedrectilinear" == value) {
        value = "rectilinear";
        output["fill_angle_increment"] = "90";
    } else if ("alignedrectilinear" == value) {
        value = "rectilinear";
    }
    if ("fan_always_on" == opt_key) {
        //min_fan_speed is already converted to default_fan_speed, just has to deactivate it if not always_on
        if (value != "1")
            output["default_fan_speed"] = "0";
    }
    if ("bridge_angle" == opt_key && "0" == value) {
        value = "!0";
    }
    if ("thumbnails_format" == opt_key) {
        // by default, no thumbnails_tag_format for png output
        if (value == "PNG")
            output["thumbnails_tag_format"] = "0";
    }
    if ("thumbnails" == opt_key && value.find('/') != std::string::npos) {
        // new (string from 2.7) .x./. , not old (Points before) type .x.
        auto [thumbnails_list, errors] = GCodeThumbnails::make_and_check_thumbnail_list_from_prusa(value);
        std::vector<Vec2d> pts;
        ConfigOptionEnum<GCodeThumbnailsFormat> opt_format;
        opt_format.value = thumbnails_list.empty() ? GCodeThumbnailsFormat::PNG : thumbnails_list.front().first;
        for (auto [format, pt] : thumbnails_list) {
            pts.push_back(pt);
        }
        value = ConfigOptionPoints(pts).serialize();
        // format (the first) is still set by prusa, no need to parse it.
        //output["thumbnails_format"] = opt_format.serialize();
    }
/*
    if ("thumbnails" == opt_key) {
        //check if their format is inside the size
        if (value.find('/') != std::string::npos) {
            std::vector<std::string> sizes;
            boost::split(sizes, value, boost::is_any_of(","), boost::token_compress_off);
            value = "";
            std::string coma = "";
            size_t added = 0;
            for (std::string &size : sizes) {
                size_t pos = size.find('/');
                assert(pos != std::string::npos);
                if (pos != std::string::npos) {
                    assert(size.find('/', pos + 1) == std::string::npos);
                    value = value + coma + size.substr(0, pos);
                } else {
                    value = value + coma + size;
                }
                coma  = ",";
                added++;
                if (added >= 2)
                    break;
            }
            //if less than 2: add 0X0 until two.
            while (added < 2) {
                value = value + coma + "0x0";
                coma  = ",";
                added++;
            }
            // format (the first) is still set by prusa, no need to parse it.
        }
    }
*/
    
    // ---- custom gcode: ----
    static const std::vector<std::pair<std::string, std::string>> custom_gcode_replace =
        {{"[temperature]", "{temperature+extruder_temperature_offset}"},
         {"{temperature}", "{temperature+extruder_temperature_offset}"},
         {"{temperature[initial_extruder]}", "{temperature[initial_extruder]+extruder_temperature_offset[initial_extruder]}"},
         {"[first_layer_temperature]", "{first_layer_temperature+extruder_temperature_offset}"},
         {"{first_layer_temperature}", "{first_layer_temperature+extruder_temperature_offset}"},
         {"{first_layer_temperature[initial_extruder]}", "{first_layer_temperature[initial_extruder]+extruder_temperature_offset[initial_extruder]}"}};

    static const std::set<t_config_option_key> custom_gcode_keys =
        {"template_custom_gcode", "toolchange_gcode", "before_layer_gcode",
         "between_objects_gcode", "end_gcode",        "layer_gcode",
         "feature_gcode",         "start_gcode",      "color_change_gcode",
         "pause_print_gcode",     "toolchange_gcode", "end_filament_gcode",
         "start_filament_gcode"};
    if (custom_gcode_keys.find(opt_key) != custom_gcode_keys.end()) {
        // check & replace
        for (auto &entry : custom_gcode_replace) {
            boost::replace_all(value, entry.first, entry.second);
        }
    }

    return output;
}

//keys that needs to go through from_prusa before beeing deserialized.
const std::unordered_set<std::string> prusa_import_to_review_keys =
{
    "thumbnails"
};

const std::vector<std::pair<t_config_option_key, t_config_option_key>> prusa_import_widths_2_spacings_for_phony_fix =
    {{"extrusion_width", "extrusion_spacing"},
     {"perimeter_extrusion_width", "perimeter_extrusion_spacing"},
     {"external_perimeter_extrusion_width", "external_perimeter_extrusion_spacing"},
     {"first_layer_extrusion_width", "first_layer_extrusion_spacing"},
     {"infill_extrusion_width", "infill_extrusion_spacing"},
     {"solid_infill_extrusion_width", "solid_infill_extrusion_spacing"},
     {"top_infill_extrusion_width", "top_infill_extrusion_spacing"}};

template<typename CONFIG_CLASS>
void _convert_from_prusa(CONFIG_CLASS& conf, const DynamicPrintConfig& global_config, bool with_phony) {
    //void convert_from_prusa(DynamicPrintConfig& conf, const DynamicPrintConfig & global_config) {
    //void convert_from_prusa(ModelConfigObject& conf, const DynamicPrintConfig& global_config) {
    std::map<std::string, std::string> results;
    for (const t_config_option_key& opt_key : conf.keys()) {
        const ConfigOption* opt = conf.option(opt_key);
        std::string serialized = opt->serialize();
        std::string key = opt_key;
        std::map<std::string, std::string> result = PrintConfigDef::from_prusa(key, serialized, global_config);
        if (key != opt_key) {
            conf.erase(opt_key);
        }
        if (!key.empty() && serialized != opt->serialize()) {
            ConfigOption* opt_new = opt->clone();
            opt_new->deserialize(serialized);
            conf.set_key_value(key, opt_new);
        }
        results.insert(result.begin(), result.end());
    }
    for (auto entry : results) {
        const ConfigOptionDef* def = print_config_def.get(entry.first);
        if (def) {
            ConfigOption* opt_new = def->default_value.get()->clone();
            opt_new->deserialize(entry.second); // note: deserialize don't set phony, only the ConfigBase::set_deserialize*
            conf.set_key_value(entry.first, opt_new);
        }
    }

    // set phony entries
    if (with_phony) {
        for (auto & [opt_key_width, opt_key_spacing] : prusa_import_widths_2_spacings_for_phony_fix) {
            // if prusa has defined a width, or if the conf has a default spacing that need to be overwritten
            if (conf.option(opt_key_width) != nullptr || conf.option(opt_key_spacing) != nullptr) {
                ConfigOption *opt_new = print_config_def.get(opt_key_spacing)->default_value.get()->clone();
                opt_new->deserialize(""); // note: deserialize don't set phony, only the ConfigBase::set_deserialize*
                opt_new->set_phony(true);
                conf.set_key_value(opt_key_spacing, opt_new);
            }
        }
    }
}

void DynamicPrintConfig::convert_from_prusa(bool with_phony) {
    _convert_from_prusa<DynamicPrintConfig>(*this, *this, with_phony);
}
void ModelConfig::convert_from_prusa(const DynamicPrintConfig& global_config, bool with_phony) {
    _convert_from_prusa<ModelConfig>(*this, global_config, with_phony);
}


template<typename CONFIG_CLASS>
void _deserialize_maybe_from_prusa(const std::map<t_config_option_key, std::string> settings,
                                           CONFIG_CLASS &                             config,
                                           const DynamicPrintConfig &                 global_config,
                                           ConfigSubstitutionContext &                config_substitutions,
                                           bool                                       with_phony,
                                           bool                                       check_prusa)
{
    std::vector<std::pair<t_config_option_key, std::string>> deleted_keys;
    std::vector<std::pair<t_config_option_key, std::string>> unknown_keys;
    const ConfigDef *def = config.def();
    for (const auto &[key, value] : settings) {
        try {
            t_config_option_key opt_key = key;
            std::string opt_value = value;
            PrintConfigDef::handle_legacy(opt_key, opt_value, false);
            if (!opt_key.empty()) {
                if (!def->has(opt_key) ||
                    (check_prusa && prusa_import_to_review_keys.find(opt_key) != prusa_import_to_review_keys.end())) {
                    unknown_keys.emplace_back(key, value);
                } else {
                    config.set_deserialize(opt_key, opt_value, config_substitutions);
                }
            } else {
                deleted_keys.emplace_back(key, value);
            }
        } catch (UnknownOptionException & /* e */) {
            // log & ignore
            if (config_substitutions.rule != ForwardCompatibilitySubstitutionRule::Disable)
                config_substitutions.add(ConfigSubstitution(key, value));
            assert(false);
        } catch (BadOptionValueException &e) {
            if (config_substitutions.rule == ForwardCompatibilitySubstitutionRule::Disable)
                throw e;
            // log the error
            const ConfigDef *def = config.def();
            if (def == nullptr)
                throw e;
            const ConfigOptionDef *optdef = def->get(key);
            config_substitutions.emplace(optdef,std::string(value), ConfigOptionUniquePtr(optdef->default_value->clone()));
        }
    }
    config.handle_legacy_composite(deleted_keys);
    // from prusa: try again with from_prusa before handle_legacy
    if (check_prusa) {
        std::map<t_config_option_key, std::string> settings_to_change;
        for (auto& [key, value] : unknown_keys) {
            t_config_option_key                        opt_key = key;
            std::map<t_config_option_key, std::string> result  = PrintConfigDef::from_prusa(opt_key, value, global_config);
            settings_to_change.insert(result.begin(), result.end());
            if (!opt_key.empty())
                //check if good this time
                PrintConfigDef::handle_legacy(opt_key, value, false);
            if (!opt_key.empty()) {
                if (!def->has(opt_key)) {
                    if (config_substitutions.rule != ForwardCompatibilitySubstitutionRule::Disable) {
                        config_substitutions.add(ConfigSubstitution(key, value));
                    }
                } else {
                    try {
                        config.set_deserialize(opt_key, value, config_substitutions);
                    } catch (BadOptionValueException &e) {
                        if (config_substitutions.rule == ForwardCompatibilitySubstitutionRule::Disable)
                            throw e;
                        // log the error
                        if (def == nullptr)
                            throw e;
                        const ConfigOptionDef *optdef = def->get(key);
                        config_substitutions.emplace(optdef, std::string(value), ConfigOptionUniquePtr(optdef->default_value->clone()));
                    }
                }
            }
        }
        for (const auto &entry : settings_to_change)
            config.set_deserialize(entry.first, entry.second, config_substitutions);
    } else {
        for (const auto& [key, value] : unknown_keys) {
            if (config_substitutions.rule != ForwardCompatibilitySubstitutionRule::Disable) {
                config_substitutions.add(ConfigSubstitution(key, value));
            }
        }
    }

    // set phony entries
    if (with_phony) {
        const ConfigDef *def = config.def();
        for (auto & [opt_key_width, opt_key_spacing] : prusa_import_widths_2_spacings_for_phony_fix) {
            const ConfigOption *opt_width = config.option(opt_key_width);
            const ConfigOption *opt_spacing = config.option(opt_key_spacing);
            if (opt_width && opt_spacing) {
                // if the config has a default spacing that need to be overwritten (if the width wasn't deserialized as phony)
                if (settings.find(opt_key_spacing) == settings.end()) {
                    if (opt_width->is_phony()) {
                        if (opt_spacing->is_phony()) {
                            ConfigOption *opt_new = opt_spacing->clone();
                            opt_new->set_phony(false);
                            config.set_key_value(opt_key_spacing, opt_new);
                        }
                    } else {
                        if (!opt_spacing->is_phony()) {
                            ConfigOption *opt_new = opt_spacing->clone();
                            opt_new->set_phony(true);
                            config.set_key_value(opt_key_spacing, opt_new);
                        }
                    }
                } else {
                    //spacing exist in the config, make sure one if phony
                    if (opt_spacing->is_phony() && opt_width->is_phony()) {
                        ConfigOption *opt_new = opt_width->clone();
                        opt_new->set_phony(false);
                        config.set_key_value(opt_key_width, opt_new);
                    }
                    if (!opt_spacing->is_phony() && !opt_width->is_phony()) {
                        ConfigOption *opt_new = opt_width->clone();
                        opt_new->set_phony(true);
                        config.set_key_value(opt_key_width, opt_new);
                    }
                }
                assert(config.option(opt_key_width)->is_phony() != config.option(opt_key_spacing)->is_phony());
            }
        }
    }
}
void deserialize_maybe_from_prusa(std::map<t_config_option_key, std::string> settings,
                                  ModelConfig &                              config,
                                  const DynamicPrintConfig &                 global_config,
                                  ConfigSubstitutionContext &                config_substitutions,
                                  bool                                       with_phony,
                                  bool                                       check_prusa)
{
    _deserialize_maybe_from_prusa(settings, config, global_config, config_substitutions, with_phony, check_prusa);
}
void deserialize_maybe_from_prusa(std::map<t_config_option_key, std::string> settings,
                                  DynamicPrintConfig &                       config,
                                  ConfigSubstitutionContext &                config_substitutions,
                                  bool                                       with_phony,
                                  bool                                       check_prusa)
{
    _deserialize_maybe_from_prusa(settings, config, config, config_substitutions, with_phony, check_prusa);
}


std::unordered_set<std::string> prusa_export_to_remove_keys = {
"allow_empty_layers",
"arc_fitting_tolerance",
"avoid_crossing_not_first_layer",
"avoid_crossing_top",
"bridge_fill_pattern",
"bridge_precision",
"bridge_overlap",
"bridge_overlap_min",
"bridge_type",
"bridged_infill_margin",
"brim_acceleration",
"brim_ears_detection_length",
"brim_ears_max_angle",
"brim_ears_pattern",
"brim_ears",
"brim_inside_holes",
"brim_per_object",
"brim_speed",
"brim_width_interior",
"chamber_temperature",
"complete_objects_one_skirt",
"complete_objects_sort",
"curve_smoothing_angle_concave",
"curve_smoothing_angle_convex",
"curve_smoothing_cutoff_dist",
"curve_smoothing_precision",
"default_speed",
"enforce_full_fill_volume",
// "exact_last_layer_height",
"external_infill_margin",
"external_perimeter_cut_corners",
"external_perimeter_extrusion_spacing",
"external_perimeter_extrusion_change_odd_layers",
"external_perimeter_fan_speed",
"external_perimeter_overlap",
"external_perimeters_hole",
"external_perimeters_nothole",
"external_perimeters_vase",
"extra_perimeters_odd_layers",
"extruder_extrusion_multiplier_speed",
"extruder_fan_offset",
"extruder_temperature_offset",
"extrusion_spacing",
"fan_kickstart",
"fan_name",
"fan_percentage",
"fan_printer_min_speed",
"fan_speedup_overhangs",
"fan_speedup_time",
"feature_gcode",
"filament_cooling_zone_pause",
"filament_custom_variables",
"filament_dip_extraction_speed",
"filament_dip_insertion_speed",
"filament_enable_toolchange_part_fan",
"filament_enable_toolchange_temp",
"filament_fill_top_flow_ratio",
"filament_first_layer_flow_ratio",
"filament_max_speed",
"filament_max_wipe_tower_speed",
"filament_melt_zone_pause",
"filament_max_overlap",
"filament_retract_lift_before_travel",
"filament_shrink",
"filament_skinnydip_distance",
"filament_toolchange_part_fan_speed",
"filament_toolchange_temp",
"filament_use_fast_skinnydip",
"filament_use_skinnydip",
"filament_wipe_advanced_pigment",
//pa
"filament_bridge_internal_pa",
"filament_bridge_pa",
"filament_brim_pa",
"filament_default_pa",
"filament_external_perimeter_pa",
"filament_first_layer_pa",
"filament_first_layer_pa_over_raft",
"filament_gap_fill_pa",
"filament_infill_pa",
"filament_ironing_pa",
"filament_overhangs_pa",
"filament_perimeter_pa",
"filament_solid_infill_pa",
"filament_support_material_pa",
"filament_support_material_interface_pa",
"filament_thin_walls_pa",
"filament_top_solid_infill_pa",
"filament_travel_pa",
//pa end
"fill_aligned_z",
"fill_angle_increment",
"fill_angle_cross",
"fill_angle_follow_model",
"fill_angle_template",
"fill_smooth_distribution",
"fill_smooth_width",
"fill_top_flow_ratio",
"fill_top_flow_ratio",
"first_layer_extrusion_spacing",
"first_layer_infill_extrusion_width",
"first_layer_infill_extrusion_spacing",
"first_layer_flow_ratio",
"first_layer_infill_speed",
"first_layer_min_speed",
"first_layer_size_compensation_layers",
"gcode_ascii",
"gcode_command_buffer",
"gcode_min_length",
"gcode_min_resolution",
"gap_fill_acceleration",
"gap_fill_extension",
"gap_fill_fan_speed",
"gap_fill_flow_match_perimeter",
"gap_fill_last",
"gap_fill_infill",
"gap_fill_min_area",
"gap_fill_max_width",
"gap_fill_min_length",
"gap_fill_min_width",
"gap_fill_overlap",
"gcode_filename_illegal_char",
"gcode_precision_e",
"gcode_precision_xyz",
"hole_size_compensation",
"hole_size_threshold",
"hole_to_polyhole_threshold",
"hole_to_polyhole_twisted",
"hole_to_polyhole",
"infill_connection",
"infill_connection_bottom",
"infill_connection_bridge",
"infill_connection_solid",
"infill_connection_top",
"infill_dense_algo",
"infill_dense",
"infill_extrusion_change_odd_layers",
"infill_extrusion_spacing",
"infill_fan_speed",
"init_z_rotate",
"internal_bridge_acceleration",
"internal_bridge_fan_speed",
"internal_bridge_speed",
"ironing_acceleration",
"ironing_angle",
"lift_min",
"max_gcode_per_second",
"max_speed_reduction",
"milling_after_z",
"milling_cutter",
"milling_diameter",
"milling_extra_size",
"milling_offset",
"milling_post_process",
"milling_speed",
"milling_toolchange_end_gcode",
"milling_toolchange_start_gcode",
"milling_z_lift",
"milling_z_offset",
"laser_activate_gcode",
"laser_diameter",
"laser_disable_gcode",
"laser_enable_gcode",
"laser_energy",
"laser_head",
"laser_support_interface_pp",
"laser_power",
"laser_offset",
"laser_toolchange_end_gcode",
"laser_toolchange_start_gcode",
"laser_z_offset",
"min_width_top_surface",
"model_precision",
"no_perimeter_unsupported_algo",
"object_gcode",
"only_one_perimeter_top_other_algo",
"only_one_perimeter_top",
"only_one_perimeter_first_layer",
"over_bridge_flow_ratio",
"overhangs_acceleration",
"overhangs_fan_speed",
"overhangs_max_slope",
"overhangs_bridge_threshold",
"overhangs_bridge_upper_layers",
"overhangs_reverse_threshold",
"overhangs_reverse",
"overhangs_speed_enforce",
"overhangs_width_speed",
"parallel_objects_step",
"perimeter_bonding",
"perimeter_extrusion_change_odd_layers",
"perimeter_extrusion_spacing",
"perimeter_direction",
"perimeter_fan_speed",
"perimeter_loop_seam",
"perimeter_loop",
"perimeter_overlap",
"perimeter_reverse",
"perimeter_round_corners",
"perimeters_hole",
"priming_position",
"print_extrusion_multiplier",
"print_first_layer_temperature",
"print_custom_variables",
"print_retract_length",
"print_retract_lift",
"print_temperature",
"printer_custom_variables",
"printhost_client_cert",
"printhost_client_cert_password",
"raft_layer_height",
"raft_interface_layer_height",
"region_gcode",
"remaining_times_type",
"resolution_internal",
"retract_lift_first_layer",
"retract_lift_top",
"retract_lift_before_travel",
"seam_angle_cost",
"seam_gap",
"seam_gap_external",
"solid_over_perimeters",
"filament_seam_gap", // filament override
"filament_seam_gap_external", // filament override
"seam_notch_all",
"seam_notch_angle",
"seam_notch_inner",
"seam_notch_outer",
"seam_travel_cost",
"seam_visibility",
"skirt_brim",
"skirt_distance_from_brim",
"skirt_extrusion_width",
"small_area_infill_flow_compensation",
"small_area_infill_flow_compensation_model",
"small_perimeter_max_length",
"small_perimeter_min_length",
"solid_fill_pattern",
"solid_infill_extrusion_change_odd_layers",
"solid_infill_extrusion_spacing",
"solid_infill_fan_speed",
"solid_infill_overlap",
"start_gcode_manual",
"solid_infill_below_layer_area",
"solid_infill_below_width",
"support_material_angle_height",
"support_material_acceleration",
"support_material_contact_distance_type",
"support_material_fan_speed",
"support_material_interface_acceleration",
"support_material_interface_angle",
"support_material_interface_angle_increment",
"support_material_interface_fan_speed",
"support_material_interface_layer_height",
"support_material_bottom_interface_pattern",
"support_material_layer_height",
"thin_perimeters_all",
"thin_perimeters",
"thin_walls_acceleration",
"thin_walls_merge",
"thin_walls_min_width",
"thin_walls_overlap",
"thin_walls_speed",
"thumbnails_color",
"thumbnails_custom_color",
"thumbnails_end_file",
"thumbnails_tag_format",
"thumbnails_with_bed",
"thumbnails_with_support",
"time_cost",
"time_estimation_compensation",
"time_start_gcode",
"time_toolchange",
"tool_name",
"top_fan_speed",
"top_infill_extrusion_spacing",
"top_solid_infill_overlap",
"travel_acceleration",
"travel_deceleration_use_target",
"wipe_advanced_algo",
"wipe_advanced_multiplier",
"wipe_advanced_nozzle_melted_volume",
"wipe_advanced",
"wipe_extra_perimeter",
"wipe_inside_depth",
"wipe_inside_end",
"wipe_inside_start",
"wipe_only_crossing",
"wipe_speed",
"filament_wipe_extra_perimeter", // filament override
"filament_wipe_inside_depth", // filament override
"filament_wipe_inside_end", // filament override
"filament_wipe_inside_start", // filament override
"filament_wipe_only_crossing", // filament override
"filament_wipe_speed", // filament override
"wipe_tower_speed",
"wipe_tower_wipe_starting_speed",
"xy_size_compensation",
"xy_inner_size_compensation",
"z_step",

"print_version",
};

std::unordered_set<std::string> prusa_export_to_change_keys =
{
"default_fan_speed", // used to convert to min_fan_speed & fan_always_on
"overhangs_width",
"overhangs_speed",
};

std::map<std::string, std::string> PrintConfigDef::to_prusa(t_config_option_key& opt_key, std::string& value, const DynamicConfig& all_conf) {
    std::map<std::string, std::string> new_entries;
    //looks if it's to be removed, or have to be transformed
    if (prusa_export_to_remove_keys.find(opt_key) != prusa_export_to_remove_keys.end()) {
        assert((all_conf.def()->get(opt_key)->mode & comPrusa) == 0);
        opt_key = "";
        value = "";
        return new_entries;
    }
    if (!opt_key.empty() && prusa_export_to_change_keys.find(opt_key) == prusa_export_to_change_keys.end()) {
        auto mode = all_conf.def()->get(opt_key)->mode;
        assert( (mode & comPrusa) == comPrusa);
    }
    if (opt_key.find("_pattern") != std::string::npos) {
        if ("smooth" == value || "smoothtriple" == value || "smoothhilbert" == value || "rectiwithperimeter" == value || "scatteredrectilinear" == value || "rectilineargapfill" == value || "sawtooth" == value) {
            value = "rectilinear";
        } else if ("concentricgapfill" == value) {
            value = "concentric";
        } else if ("monotonicgapfill" == value) {
            value = "monotonic";
        }
        if (all_conf.has("fill_angle_increment") && ((int(all_conf.option("fill_angle_increment")->get_float())-90)%180) == 0 && "rectilinear" == value
            && ("fill_pattern" == opt_key || "top_fill_pattern" == opt_key)) {
            value = "alignedrectilinear";
        }
        if ("support_material_top_interface_pattern" == opt_key) {
            opt_key = "support_material_interface_pattern";
        }
    } else if ("seam_position" == opt_key) {
        if ("cost" == value) {
            value = "nearest";
        }else if ("allrandom" == value) {
            value = "random";
        }else if ("contiguous" == value) {
            value = "aligned";
        }
    } else if ("first_layer_size_compensation" == opt_key) {
        opt_key = "elefant_foot_compensation";
        if (!value.empty()) {
            if (value[0] == '-') {
                value = value.substr(1);
            } else {
                value = "0";
            }
        }
    } else if ("elephant_foot_min_width" == opt_key) {
        opt_key = "elefant_foot_min_width";
    } else if("first_layer_acceleration" == opt_key || "first_layer_acceleration_over_raft" == opt_key) {
        if (value.find("%") != std::string::npos) {
            // can't support %, so we uese the default accel a baseline for half-assed conversion
            value = std::to_string(all_conf.get_abs_value(opt_key, all_conf.get_computed_value("default_acceleration")));
        }
    } else if ("infill_acceleration" == opt_key || "solid_infill_acceleration" == opt_key || "top_solid_infill_acceleration" == opt_key
        || "bridge_acceleration" == opt_key || "default_acceleration" == opt_key || "perimeter_acceleration" == opt_key
        || "overhangs_speed" == opt_key || "ironing_speed" == opt_key || "perimeter_speed" == opt_key 
        || "infill_speed" == opt_key || "bridge_speed" == opt_key || "support_material_speed" == opt_key
        || "max_print_speed" == opt_key
        ) {
        // remove '%', change 0 with different meanings
        if (value.find("%") != std::string::npos) {
            if (value == "100%")
                value = "0";
            else
                value = std::to_string(all_conf.get_computed_value(opt_key));
        }
        // infill_acceleration & solid_infill_acceleration dep are inverted
        if ("infill_acceleration" == opt_key && value == "0") {
            value = std::to_string(all_conf.get_computed_value("solid_infill_acceleration"));
        } else if ("solid_infill_acceleration" == opt_key && value == "0") {
            value = std::to_string(all_conf.get_computed_value("default_acceleration"));
        }
    } else if ("gap_fill_speed" == opt_key && all_conf.has("gap_fill_enabled") && !all_conf.option<ConfigOptionBool>("gap_fill_enabled")->value) {
        value = "0";
    } else if ("bridge_flow_ratio" == opt_key && all_conf.has("bridge_flow_ratio")) {
        value = boost::lexical_cast<std::string>(all_conf.option<ConfigOptionPercent>("bridge_flow_ratio")->get_abs_value(1));
    } else if ("overhangs_width" == opt_key) {
        opt_key = "overhangs";
        if ((!value.empty() && value.front() == '!') || !all_conf.is_enabled("overhangs_width_speed")) {
            value = "0";
        } else {
            value = "1";
        }
    } else if ("support_material_contact_distance_top" == opt_key) {
        opt_key = "support_material_contact_distance";
        //default : get the top value or 0.2 if a %
        if (value.find("%") != std::string::npos)
            value = "0.2";
        try { //avoid most useless cheks and multiple corners cases with this try catch
            SupportZDistanceType dist_type = all_conf.option<ConfigOptionEnum<SupportZDistanceType>>("support_material_contact_distance_type")->value;
            if (SupportZDistanceType::zdNone == dist_type) {
                value = "0";
            } else {
                double val = all_conf.option<ConfigOptionFloatOrPercent>("support_material_contact_distance_top")->get_abs_value(all_conf.option<ConfigOptionFloats>("nozzle_diameter")->get_at(0));
                if (SupportZDistanceType::zdFilament == dist_type) { // not exact but good enough effort
                    val += all_conf.option<ConfigOptionFloats>("nozzle_diameter")->get_at(0);
                    val -= all_conf.get_computed_value("layer_height", 0);
                }
                value = boost::lexical_cast<std::string>(val);
            }
        }
        catch (...) {
        }
    } else if ("support_material_contact_distance_bottom" == opt_key) {
            opt_key = "support_material_bottom_contact_distance";
            //default : get the top value or 0.2 if a %
            if (value.find("%") != std::string::npos)
                value = "0.2";
            try { //avoid most useless cheks and multiple corners cases with this try catch
                SupportZDistanceType dist_type = all_conf.option<ConfigOptionEnum<SupportZDistanceType>>("support_material_contact_distance_type")->value;
                if (SupportZDistanceType::zdNone == dist_type) {
                    value = "0";
                } else {
                    double val = all_conf.option<ConfigOptionFloatOrPercent>("support_material_contact_distance_bottom")->get_abs_value(all_conf.option<ConfigOptionFloats>("nozzle_diameter")->get_at(0));
                    if (SupportZDistanceType::zdFilament == dist_type) { // not exact but good enough effort
                        val += all_conf.option<ConfigOptionFloats>("nozzle_diameter")->get_at(0);
                        val -= all_conf.get_computed_value("layer_height", 0);
                    }
                    value = boost::lexical_cast<std::string>(val);
                }
            }
            catch (...) {
            }
        } else if ("gcode_flavor" == opt_key) {
        if ("sprinter" == value)
            value = "reprap";
        else if ("lerdge" == value)
            value = "marlin";
        else if ("klipper" == value)
            value = "reprap";
    } else if ("host_type" == opt_key) {
        if ("klipper" == value)
            value = "octoprint";
    } else if (opt_key.find("extrusion_width") != std::string::npos) {
        if (std::set<std::string>{"extrusion_width", "first_layer_extrusion_width", "perimeter_extrusion_width", "external_perimeter_extrusion_width", 
            "infill_extrusion_width", "solid_infill_extrusion_width", "top_infill_extrusion_width", "support_material_extrusion_width"}.count(opt_key) > 0) {
            const ConfigOptionFloatOrPercent* opt = all_conf.option<ConfigOptionFloatOrPercent>(opt_key);
            if (opt->is_phony() || opt->percent) {
                if (opt->percent) {
                    ConfigOptionFloat opt_temp{ opt->get_abs_value(all_conf.option<ConfigOptionFloats>("nozzle_diameter")->get_at(0)) };
                    value = opt_temp.serialize();
                } else {
                    //bypass the phony kill switch from Config::opt_serialize
                    value = opt->serialize();
                }
            }
        }
    }
    if ("infill_anchor_max" == opt_key) {
        //it's infill_anchor == 0 that disable it for prusa
        if (all_conf.opt_serialize("infill_connection") == "notconnected") {
            value = "0";
        }
    }
    if ("brim_width" == opt_key) {
        double brim_width_interior = all_conf.get_computed_value("brim_width_interior");
        if (value == "0" && brim_width_interior == 0)
            new_entries["brim_type"] = "no_brim";
        else if (value != "0" && brim_width_interior == 0)
            new_entries["brim_type"] = "outer_only";
        else if (value == "0" && brim_width_interior != 0)
            new_entries["brim_type"] = "inner_only";
        else if (value != "0" && brim_width_interior != 0)
            new_entries["brim_type"] = "outer_and_inner";
    }
    if ("output_format" == opt_key) {
        opt_key = "sla_archive_format";
    }
    if ("host_type" == opt_key) {
        if ("klipper" == value || "mpmdv2" == value || "monoprice" == value) value = "octoprint";
    }
    if ("fan_below_layer_time" == opt_key) {
        if (value.find('.') != std::string::npos)
            value = value.substr(0, value.find('.'));
    }
    if ("bed_custom_texture" == opt_key || "Bed custom texture" == opt_key) {
        value = Slic3r::find_full_path(value, value).generic_string();
    }
    if ("default_fan_speed" == opt_key) {
        if (!value.empty() && value.front() == '!') {
            value = "1";
        }
        if (value == "0") {
            opt_key = "min_fan_speed";
            value = std::to_string(all_conf.option("fan_printer_min_speed")->get_float());
            new_entries["fan_always_on"] = "0";
        } else {
            opt_key = "min_fan_speed";
            new_entries["fan_always_on"] = "1";
        }
    }
    if ("bridge_fan_speed" == opt_key) {
        if (!value.empty() && value.front() == '!') {
            value = std::to_string(all_conf.option("fan_printer_min_speed")->get_float());
        }
    }
    if ("travel_ramping_lift" == opt_key && "1" == value) {
        // also add travel_max_lift from retract_lift & same from filament
        new_entries["travel_ramping_lift"] = all_conf.option("retract_lift")->serialize();
        new_entries["filament_travel_ramping_lift"] = all_conf.option("filament_retract_lift")->serialize();
    }

    // compute max & min height from % to flat value
    if ("min_layer_height" == opt_key || "max_layer_height" == opt_key) {
        if ("max_layer_height" == opt_key && !value.empty() && value.front() == '!') {
            value = "0";
        } else {
            ConfigOptionFloats computed_opt;
            const ConfigOptionFloatsOrPercents *current_opt = all_conf.option<ConfigOptionFloatsOrPercents>(opt_key);
            const ConfigOptionFloats *nozzle_diameters = all_conf.option<ConfigOptionFloats>("nozzle_diameter");
            assert(current_opt && nozzle_diameters);
            assert(current_opt->size() == nozzle_diameters->size());
            for (int i = 0; i < current_opt->size(); i++) {
                computed_opt.set_at(current_opt->get_abs_value(i, nozzle_diameters->get_at(i)), i);
            }
            assert(computed_opt.size() == nozzle_diameters->size());
            value = computed_opt.serialize();
        }
    }
    if ("arc_fitting" == opt_key && "bambu" == value) {
        value = "emit_center";
    }
    
    if ("thumbnails" == opt_key) {
    // add format to thumbnails
        const ConfigOptionEnum<GCodeThumbnailsFormat> *format_opt = all_conf.option<ConfigOptionEnum<GCodeThumbnailsFormat>>("thumbnails_format");
        std::string format = format_opt->serialize();
        std::vector<std::string> sizes;
        boost::split(sizes, value, boost::is_any_of(","), boost::token_compress_off);
        value = "";
        std::string coma = "";
        for (std::string &size : sizes) {
            //if first or second dimension is 0: ignore.
            if (size.find("0x") == 0 || size.find("x0") + 2 == size.size())
                continue;
            assert(size.find('/') == std::string::npos);
            value = value + coma + size + std::string("/") + format;
            coma = ",";
        }
    }

    // disabled
    if (!value.empty() && value.front() == '!') {
        // ----- -1 -> disabled -----
        static const std::set<t_config_option_key> minus_1_is_disabled = {"overhangs_bridge_threshold",
              "overhangs_bridge_upper_layers", "perimeters_hole", "support_material_bottom_interface_layers"};
        // ---- filament override ------
        if (boost::starts_with(opt_key, "filament_")) {
            std::string extruder_key = opt_key.substr(strlen("filament_"));
            if (print_config_def.filament_override_option_keys().find(extruder_key) !=
                print_config_def.filament_override_option_keys().end()) {
                value = "nil";
            }
        }
        if (boost::starts_with(opt_key, "material_ow_")) {
            std::string normal_key = opt_key.substr(strlen("material_ow_"));
            if (print_config_def.material_overrides_option_keys().find(opt_key) !=
                print_config_def.material_overrides_option_keys().end()) {
                value = "nil";
            }
        }
        if ("idle temperature" == opt_key) {
            value = "nil";
        } else if (minus_1_is_disabled.find(opt_key) != minus_1_is_disabled.end()) {
            value = "-1";
        }
    }
    
    if ("bridge_angle" == opt_key && !value.empty() && value.front() == '!') {
        value = "0";
    }

    // ---- custom gcode: ----
    static const std::vector<std::pair<std::string, std::string>> custom_gcode_replace =
        {{"extruder_temperature_offset[initial_extruder]", "0"}, {"extruder_temperature_offset", "0"}};
    static const std::set<t_config_option_key> custom_gcode_keys =
        {"template_custom_gcode", "toolchange_gcode", "before_layer_gcode",
         "between_objects_gcode", "end_gcode",        "layer_gcode",
         "feature_gcode",         "start_gcode",      "color_change_gcode",
         "pause_print_gcode",     "toolchange_gcode", "end_filament_gcode",
         "start_filament_gcode"};
    if (custom_gcode_keys.find(opt_key) != custom_gcode_keys.end()) {
        // check & replace
        for (auto &entry : custom_gcode_replace) {
            boost::replace_all(value, entry.first, entry.second);
        }
    }
    // --end-- custom gcode: --end--


    return new_entries;
}

const PrintConfigDef print_config_def;

DynamicPrintConfig DynamicPrintConfig::full_print_config()
{
	return DynamicPrintConfig((const PrintRegionConfig&)FullPrintConfig::defaults());
}

DynamicPrintConfig::DynamicPrintConfig(const StaticPrintConfig& rhs) : DynamicConfig(rhs, rhs.keys_ref())
{
}

DynamicPrintConfig* DynamicPrintConfig::new_from_defaults_keys(const std::vector<std::string> &keys)
{
    auto *out = new DynamicPrintConfig();
    out->apply_only(FullPrintConfig::defaults(), keys);
    return out;
}

const ConfigOption *MultiPtrPrintConfig::optptr(const t_config_option_key &opt_key) const
{
    for (ConfigBase *conf : storages) {
#ifdef _DEBUG
        assert(conf->exists());
#endif
        const ConfigOption *opt = conf->optptr(opt_key);
        if (opt)
            return opt;
    }
    return nullptr;
}
ConfigOption *MultiPtrPrintConfig::optptr(const t_config_option_key &opt_key, bool create)
{
    assert(!create);
    for (ConfigBase *conf : storages) {
#ifdef _DEBUG
        assert(conf->exists());
#endif
        ConfigOption *opt = conf->optptr(opt_key);
        if (opt)
            return opt;
    }
    return nullptr;
}
t_config_option_keys MultiPtrPrintConfig::keys() const
{
    assert(false);
    // shouldn't need ot call that
    t_config_option_keys keys;
    for (ConfigBase *conf : storages) {
#ifdef _DEBUG
        assert(conf->exists());
#endif
        append(keys, conf->keys());
    }
    return keys;
}

PrinterTechnology printer_technology(const ConfigBase& cfg)
{
    const ConfigOptionEnum<PrinterTechnology>* opt = cfg.option<ConfigOptionEnum<PrinterTechnology>>("printer_technology");

    if (opt) return opt->value;

    const ConfigOptionBool* export_opt = cfg.option<ConfigOptionBool>("export_sla");
    if (export_opt && export_opt->get_bool()) return ptSLA;

    export_opt = cfg.option<ConfigOptionBool>("export_gcode");
    if (export_opt && export_opt->get_bool()) return ptFFF;

    return ptUnknown;
}

OutputFormat output_format(const ConfigBase& cfg)
{
    std::cerr << "Detected technology " << printer_technology(cfg) << "\n";
    if (printer_technology(cfg) == ptFFF) return ofGCode;
    const ConfigOptionEnum<OutputFormat>* opt = cfg.option<ConfigOptionEnum<OutputFormat>>("output_format");
    if (opt) return opt->value;

    return ofUnknown;
}

/*
double min_object_distance(const ConfigBase &cfg)
{
    const ConfigOptionEnum<PrinterTechnology> *opt_printer_technology = cfg.option<ConfigOptionEnum<PrinterTechnology>>("printer_technology");
    auto printer_technology = opt_printer_technology ? opt_printer_technology->value : ptUnknown;
    double ret = 0.;

    if (printer_technology == ptSLA)
        ret = 6.;
    else {
        auto ecr_opt = cfg.option<ConfigOptionFloat>("extruder_clearance_radius");
        auto dd_opt  = cfg.option<ConfigOptionFloat>("duplicate_distance");
        auto co_opt  = cfg.option<ConfigOptionBool>("complete_objects");

        if (!ecr_opt || !dd_opt || !co_opt) ret = 0.;
        else {
            // min object distance is max(duplicate_distance, clearance_radius)
            ret = (co_opt->value && ecr_opt->value > dd_opt->value) ?
                      ecr_opt->value : dd_opt->value;
        }
    }

    return ret;
}*/

double min_object_distance(const PrintConfig& config)
{
    return min_object_distance(static_cast<const ConfigBase*>(&config));
}

double min_object_distance(const ConfigBase *config, double ref_height /* = 0*/)
{
    if (printer_technology(*config) == ptSLA) return 6.;
    
    const ConfigOptionFloat* dd_opt = config->option<ConfigOptionFloat>("duplicate_distance");
    //test if called from usaslicer::l240 where it's called on an empty config...
    if (dd_opt == nullptr) return 0;

    double base_dist = 0;
    //std::cout << "START min_object_distance =>" << base_dist << "\n";
    const ConfigOptionBool* co_opt = config->option<ConfigOptionBool>("complete_objects");
    if ((config->option("parallel_objects_step")->get_float() > 0) || (co_opt && co_opt->value)) {
        double skirt_dist = 0;
        try {
            std::vector<double> vals = dynamic_cast<const ConfigOptionFloats*>(config->option("nozzle_diameter"))->get_values();
            double max_nozzle_diam = 0;
            for (double val : vals) max_nozzle_diam = std::fmax(max_nozzle_diam, val);

            // min object distance is max(duplicate_distance, clearance_radius)
            // /2 becasue we only count the grawing for the current object
            //add 1 as safety offset.
            double extruder_clearance_radius = config->option("extruder_clearance_radius")->get_float() / 2;
            if (extruder_clearance_radius > base_dist) {
                base_dist = extruder_clearance_radius;
            }

            // we use the max nozzle, just to be on the safe side
            //ideally, we should use print::first_layer_height()
            const double first_layer_height = dynamic_cast<const ConfigOptionFloatOrPercent*>(config->option("first_layer_height"))->get_abs_value(max_nozzle_diam);
            //add the skirt
            int skirts = config->option("skirts")->get_int();
            if (skirts > 0 && ref_height == 0)
                skirts += config->option("skirt_brim")->get_int();
            if (skirts > 0 && config->option("skirt_height")->get_int() >= 1 && !config->option("complete_objects_one_skirt")->get_bool()) {
                float overlap_ratio = 1;
                //can't know the extruder, so we settle on the worst: 100%
                //if (config->option<ConfigOptionPercents>("filament_max_overlap")) overlap_ratio = config->get_computed_value("filament_max_overlap");
                if (ref_height == 0) {
                    skirt_dist = config->option("skirt_distance")->get_float();
                    Flow skirt_flow = Flow::new_from_config_width(
                        frPerimeter,
                        *Flow::extrusion_width_option("skirt", *config),
                        *Flow::extrusion_spacing_option("skirt", *config),
                        (float)max_nozzle_diam,
                        (float)first_layer_height,
                        overlap_ratio,
                        0
                    );
                    skirt_dist += skirt_flow.width() + (skirt_flow.spacing() * ((double)skirts - 1));
                    base_dist = std::max(base_dist, skirt_dist + 1);
                    //set to 0 becasue it's incorporated into the base_dist, so we don't want to be added in to it again.
                    skirt_dist = 0;
                } else {
                    double skirt_height = ((double)config->option("skirt_height")->get_int() - 1) * config->get_computed_value("layer_height") + first_layer_height;
                    if (ref_height <= skirt_height) {
                        skirt_dist = config->option("skirt_distance")->get_float();
                        Flow skirt_flow = Flow::new_from_config_width(
                            frPerimeter,
                            *Flow::extrusion_width_option("skirt", *config),
                            *Flow::extrusion_spacing_option("skirt", *config),
                            (float)max_nozzle_diam,
                            (float)first_layer_height,
                            overlap_ratio,
                            0
                        );
                        skirt_dist += skirt_flow.width() + (skirt_flow.spacing() * ((double)skirts - 1));
                    }
                }
            }
        }
        catch (const std::exception & ex) {
            boost::nowide::cerr << ex.what() << std::endl;
        }
        return base_dist + skirt_dist;
    }
    return base_dist;
}

void DynamicPrintConfig::normalize_fdm()
{
    if (this->has("extruder")) {
        int extruder = this->option("extruder")->get_int();
        this->erase("extruder");
        if (extruder != 0) {
            if (!this->has("infill_extruder"))
                this->option<ConfigOptionInt>("infill_extruder", true)->value = (extruder);
            if (!this->has("perimeter_extruder"))
                this->option<ConfigOptionInt>("perimeter_extruder", true)->value = (extruder);
            // Don't propagate the current extruder to support.
            // For non-soluble supports, the default "0" extruder means to use the active extruder,
            // for soluble supports one certainly does not want to set the extruder to non-soluble.
            // if (!this->has("support_material_extruder"))
            //     this->option("support_material_extruder", true)->setInt(extruder);
            // if (!this->has("support_material_interface_extruder"))
            //     this->option("support_material_interface_extruder", true)->setInt(extruder);
        }
    }
    if (this->has("first_layer_extruder"))
        this->erase("first_layer_extruder");

    if (this->has("wipe_tower_extruder")) {
        // If invalid, replace with 0.
        int extruder = this->opt<ConfigOptionInt>("wipe_tower_extruder")->value;
        int num_extruders = this->opt<ConfigOptionFloats>("nozzle_diameter")->size();
        if (extruder < 0 || extruder > num_extruders)
            this->opt<ConfigOptionInt>("wipe_tower_extruder")->value = 0;
    }

    if (!this->has("solid_infill_extruder") && this->has("infill_extruder"))
        this->option<ConfigOptionInt>("solid_infill_extruder", true)->value = (this->option("infill_extruder")->get_int());

    if (this->has("spiral_vase") && this->opt<ConfigOptionBool>("spiral_vase", true)->value) {
        {
            // this should be actually done only on the spiral layers instead of all
            auto* opt = this->opt<ConfigOptionBools>("retract_layer_change", true);
            opt->set(std::vector<uint8_t>(opt->size(), false));  // set all values to false
            // Disable retract on layer change also for filament overrides.
            auto* opt_n = this->opt<ConfigOptionBools>("filament_retract_layer_change", true);
            opt_n->set(std::vector<uint8_t>(opt_n->size(), false));  // Set all values to false.
        }
        {
            this->opt<ConfigOptionInt>("top_solid_layers", true)->value = 0;
            this->opt<ConfigOptionPercent>("fill_density", true)->value = 0;
            this->opt<ConfigOptionEnumGeneric>("perimeter_generator", true)->value = (int)PerimeterGeneratorType::Classic;
            this->opt<ConfigOptionBool>("support_material", true)->value = false;
            this->opt<ConfigOptionInt>("solid_over_perimeters")->value = 0;
            this->opt<ConfigOptionInt>("support_material_enforce_layers")->value = 0;
            // this->opt<ConfigOptionBool>("exact_last_layer_height", true)->value = false;
            this->opt<ConfigOptionBool>("infill_dense", true)->value = false;
            this->opt<ConfigOptionBool>("extra_perimeters", true)->value = false;
            this->opt<ConfigOptionBool>("extra_perimeters_on_overhangs", true)->value = false;
            this->opt<ConfigOptionBool>("extra_perimeters_odd_layers", true)->value = false;
            this->opt<ConfigOptionBool>("overhangs_reverse", true)->value = false; 
            this->opt<ConfigOptionBool>("perimeter_reverse", true)->value = false; 
        }
    }

// merill : why?
//    if (auto* opt_gcode_resolution = this->opt<ConfigOptionFloat>("gcode_resolution", false); opt_gcode_resolution)
//        // Resolution will be above 1um.
//        opt_gcode_resolution->value = std::max(opt_gcode_resolution->value, 0.001);

    if (auto *opt_min_bead_width = this->opt<ConfigOptionFloat>("min_bead_width", false); opt_min_bead_width)
        opt_min_bead_width->value = std::max(opt_min_bead_width->value, 0.001);
    if (auto *opt_wall_transition_length = this->opt<ConfigOptionFloat>("wall_transition_length", false); opt_wall_transition_length)
        opt_wall_transition_length->value = std::max(opt_wall_transition_length->value, 0.001);

    if (auto *opt_min_bead_width = this->opt<ConfigOptionFloat>("min_bead_width", false); opt_min_bead_width)
        opt_min_bead_width->value = std::max(opt_min_bead_width->value, 0.001);
    if (auto *opt_wall_transition_length = this->opt<ConfigOptionFloat>("wall_transition_length", false); opt_wall_transition_length)
        opt_wall_transition_length->value = std::max(opt_wall_transition_length->value, 0.001);
}

void  handle_legacy_sla(DynamicPrintConfig& config)
{
    for (std::string corr : {"relative_correction", "material_correction"}) {
        if (config.has(corr)) {
            if (std::string corr_x = corr + "_x"; !config.has(corr_x)) {
                auto* opt = config.opt<ConfigOptionFloat>(corr_x, true);
                opt->value = config.opt<ConfigOptionFloats>(corr)->get_at(0);
            }

            if (std::string corr_y = corr + "_y"; !config.has(corr_y)) {
                auto* opt = config.opt<ConfigOptionFloat>(corr_y, true);
                opt->value = config.opt<ConfigOptionFloats>(corr)->get_at(0);
            }

            if (std::string corr_z = corr + "_z"; !config.has(corr_z)) {
                auto* opt = config.opt<ConfigOptionFloat>(corr_z, true);
                opt->value = config.opt<ConfigOptionFloats>(corr)->get_at(1);
            }
        }
    }
}

void DynamicPrintConfig::set_num_extruders(unsigned int num_extruders)
{
    const auto &defaults = FullPrintConfig::defaults();
    for (const std::string &key : print_config_def.extruder_option_keys()) {
        if (key == "default_filament_profile")
            // Don't resize this field, as it is presented to the user at the "Dependencies" page of the Printer profile and we don't want to present
            // empty fields there, if not defined by the system profile.
            continue;
        auto *opt = this->option(key, false);
        assert(opt != nullptr && opt->is_vector());
        if (opt != nullptr && opt->is_vector())
            static_cast<ConfigOptionVectorBase*>(opt)->resize(num_extruders, defaults.option(key));
    }
}

void DynamicPrintConfig::set_num_milling(unsigned int num_milling)
{
    const auto& defaults = FullPrintConfig::defaults();
    for (const std::string& key : print_config_def.milling_option_keys()) {
        auto* opt = this->option(key, false);
        assert(opt != nullptr);
        assert(opt->is_vector());
        if (opt != nullptr && opt->is_vector())
            static_cast<ConfigOptionVectorBase*>(opt)->resize(num_milling, defaults.option(key));
    }
}

void DynamicPrintConfig::set_num_laser(unsigned int num_laser)
{
    const auto& defaults = FullPrintConfig::defaults();
    for (const std::string& key : print_config_def.laser_option_keys()) {
        auto* opt = this->option(key, false);
        assert(opt != nullptr);
        assert(opt->is_vector());
        if (opt != nullptr && opt->is_vector())
            static_cast<ConfigOptionVectorBase*>(opt)->resize(num_laser, defaults.option(key));
    }
}

std::string DynamicPrintConfig::validate()
{
    // Full print config is initialized from the defaults.
    const ConfigOption *opt = this->option("printer_technology", false);
    auto printer_technology = (opt == nullptr) ? ptFFF : static_cast<PrinterTechnology>(dynamic_cast<const ConfigOptionEnumGeneric*>(opt)->value);
    switch (printer_technology) {
    case ptFFF:
    {
        FullPrintConfig fpc;
        fpc.apply(*this, true);
        // Verify this print options through the FullPrintConfig.
        return Slic3r::validate(fpc);
    }
    default:
        //FIXME no validation on SLA data?
        return std::string();
    }
}

template<typename TYPE>
const TYPE* find_option(const t_config_option_key &opt_key,const  DynamicPrintConfig* default_config, const std::vector<const DynamicPrintConfig*> &other_config) {
    const TYPE* option = default_config->option<TYPE>(opt_key);
    if (option)
        return option;
    for (const DynamicPrintConfig* conf : other_config) {
        option = conf->option<TYPE>(opt_key);
        if (option)
            return option;
    }
    return nullptr;
}

const DynamicPrintConfig* DynamicPrintConfig::update_phony(const std::vector<const DynamicPrintConfig*> config_collection, bool exclude_default_extrusion /*= false*/) {
    const DynamicPrintConfig* something_changed = nullptr;
    //update width/spacing links
    const char* widths[] = { "", "external_perimeter_", "perimeter_", "infill_", "solid_infill_", "top_infill_", "support_material_", "first_layer_", "first_layer_infill_", "skirt_" };
    for (size_t i = exclude_default_extrusion?1:0; i < sizeof(widths) / sizeof(widths[i]); ++i) {
        std::string key_width(widths[i]);
        key_width += "extrusion_width";
        std::string key_spacing(widths[i]);
        key_spacing += "extrusion_spacing";
        ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>(key_width);
        ConfigOptionFloatOrPercent* spacing_option = this->option<ConfigOptionFloatOrPercent>(key_spacing);
        if (width_option && spacing_option){
            const DynamicPrintConfig* returned_value;
            if (!spacing_option->is_phony() && width_option->is_phony())
                returned_value = value_changed(key_spacing, config_collection);
             else
                returned_value = value_changed(key_width, config_collection);
             if (something_changed == nullptr)
                something_changed = returned_value;
        }
    }

    return something_changed;
}

//note: width<-> spacing conversion is done via float, so max 6-7 digit of precision.
const DynamicPrintConfig* DynamicPrintConfig::value_changed(const t_config_option_key& opt_key, const std::vector<const DynamicPrintConfig*> config_collection) {
    if (opt_key == "layer_height") {
        const ConfigOptionFloat* layer_height_option = find_option<ConfigOptionFloat>("layer_height", this, config_collection);
        //if bad layer height, slip to be able to go to the check part without outputing exceptions.
        if (layer_height_option && layer_height_option->value < EPSILON)
            return nullptr;
        if (this->update_phony(config_collection) != nullptr)
            return this;
        return nullptr;
    }
    if (opt_key == "filament_max_overlap" || opt_key == "perimeter_overlap" ||
        opt_key == "external_perimeter_overlap" || opt_key == "solid_infill_overlap" ||
        opt_key == "top_solid_infill_overlap") {
        if (this->option("extrusion_width")) {
            if (this->update_phony(config_collection) != nullptr) {
                return this;
            }
        }
        return nullptr;
    }

    bool something_changed = false;
    // width -> spacing
    if (opt_key.find("extrusion_spacing") != std::string::npos) {
        const ConfigOptionFloats* nozzle_diameter_option = find_option<ConfigOptionFloats>("nozzle_diameter", this, config_collection);
        const ConfigOptionFloat* layer_height_option = find_option<ConfigOptionFloat>("layer_height", this, config_collection);
        ConfigOptionFloatOrPercent* spacing_option = this->option<ConfigOptionFloatOrPercent>(opt_key);
        if (layer_height_option && spacing_option && nozzle_diameter_option) {
            //compute spacing with current height and change the width
            double max_nozzle_diameter = 0;
            for (double dmr : nozzle_diameter_option->get_values())
                max_nozzle_diameter = std::max(max_nozzle_diameter, dmr);
            double spacing_value = spacing_option->get_abs_value(max_nozzle_diameter);
            float overlap_ratio = 1;
            const ConfigOptionPercents* filament_max_overlap_option = find_option<ConfigOptionPercents>("filament_max_overlap", this, config_collection);
            if (filament_max_overlap_option) overlap_ratio = filament_max_overlap_option->get_abs_value(0, 1.);
            Flow flow = Flow::new_from_spacing(spacing_value, max_nozzle_diameter,layer_height_option->value, overlap_ratio, false);
            //test for valid height. If too high, revert to round shape
            if (spacing_value > 0 && flow.height() > spacing_value / (1 - (1. - 0.25 * PI) * flow.spacing_ratio())) {
                flow = flow.with_width(spacing_value / (1 - (1. - 0.25 * PI) * flow.spacing_ratio()));
                flow = flow.with_height(flow.width());
            }
            if (opt_key == "extrusion_spacing") {
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "first_layer_extrusion_spacing") {
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("first_layer_extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "first_layer_infill_extrusion_spacing") {
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("first_layer_infill_extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "perimeter_extrusion_spacing") {
                const ConfigOptionPercent* perimeter_overlap_option = find_option<ConfigOptionPercent>("perimeter_overlap", this, config_collection);
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("perimeter_extrusion_width");
                if (width_option && perimeter_overlap_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if(spacing_value == 0)
                        width_option->value = 0;
                    else {
                        float spacing_ratio = (std::min(flow.spacing_ratio(), float(perimeter_overlap_option->get_abs_value(1))));
                        flow = flow.with_width( spacing_option->get_abs_value(max_nozzle_diameter) + layer_height_option->value * (1. - 0.25 * PI) * spacing_ratio);
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    }
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "external_perimeter_extrusion_spacing") {
                const ConfigOptionPercent* external_perimeter_overlap_option = find_option<ConfigOptionPercent>("external_perimeter_overlap", this, config_collection);
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("external_perimeter_extrusion_width");
                if (width_option && external_perimeter_overlap_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else {
                        float spacing_ratio = (std::min(flow.spacing_ratio() / 2, float(external_perimeter_overlap_option->get_abs_value(0.5))));
                        flow = flow.with_width(spacing_option->get_abs_value(max_nozzle_diameter) + layer_height_option->value * (1. - 0.25 * PI) * spacing_ratio);
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    }
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "infill_extrusion_spacing") {
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("infill_extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "solid_infill_extrusion_spacing") {
                const ConfigOptionPercent* solid_infill_overlap_option = find_option<ConfigOptionPercent>("solid_infill_overlap", this, config_collection);
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("solid_infill_extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else {
                        float spacing_ratio = (std::min(flow.spacing_ratio(), float(solid_infill_overlap_option->get_abs_value(1))));
                        flow = flow.with_width(spacing_option->get_abs_value(max_nozzle_diameter) + layer_height_option->value * (1. - 0.25 * PI) * spacing_ratio);
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    }
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            if (opt_key == "top_infill_extrusion_spacing") {
                const ConfigOptionPercent* top_solid_infill_overlap_option = find_option<ConfigOptionPercent>("top_solid_infill_overlap", this, config_collection);
                ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>("top_infill_extrusion_width");
                if (width_option) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    if (spacing_value == 0)
                        width_option->value = 0;
                    else {
                        float spacing_ratio = (std::min(flow.spacing_ratio(), float(top_solid_infill_overlap_option->get_abs_value(1))));
                        flow = flow.with_width(spacing_option->get_abs_value(max_nozzle_diameter) + layer_height_option->value * (1. - 0.25 * PI) * spacing_ratio);
                        width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    }
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                }
            }
            /*if (opt_key == "support_material_extrusion_spacing") {
                if (spacing_option->percent)
                    this->set_key_value("support_material_extrusion_width", new ConfigOptionFloatOrPercent(std::round(100 * flow.width / max_nozzle_diameter), true));
                else
                    this->set_key_value("support_material_extrusion_width", new ConfigOptionFloatOrPercent(std::round(flow.width * 10000) / 10000, false));
                something_changed = true;
            }
            if (opt_key == "skirt_extrusion_spacing") {
                if (spacing_option->percent)
                    this->set_key_value("skirt_extrusion_width", new ConfigOptionFloatOrPercent(std::round(100 * flow.width / max_nozzle_diameter), true));
                else
                    this->set_key_value("skirt_extrusion_width", new ConfigOptionFloatOrPercent(std::round(flow.width * 10000) / 10000, false));
                something_changed = true;
            }*/
        }
    }
    if (opt_key.find("extrusion_width") != std::string::npos) {
        const ConfigOptionFloats* nozzle_diameter_option = find_option<ConfigOptionFloats>("nozzle_diameter", this, config_collection);
        const ConfigOptionFloat* layer_height_option = find_option<ConfigOptionFloat>("layer_height", this, config_collection);
        ConfigOptionFloatOrPercent* default_width_option = this->option<ConfigOptionFloatOrPercent>("extrusion_width");
        ConfigOptionFloatOrPercent* width_option = this->option<ConfigOptionFloatOrPercent>(opt_key);
        float overlap_ratio = 1;
        const ConfigOptionPercents* filament_max_overlap_option = find_option<ConfigOptionPercents>("filament_max_overlap", this, config_collection);
        if (filament_max_overlap_option) overlap_ratio = filament_max_overlap_option->get_abs_value(0, 1.);
        if (layer_height_option && width_option && nozzle_diameter_option) {
            //compute spacing with current height and change the width
            float max_nozzle_diameter = 0;
            for (double dmr : nozzle_diameter_option->get_values())
                max_nozzle_diameter = std::max(max_nozzle_diameter, (float)dmr);
            ConfigOptionFloatOrPercent* spacing_option = nullptr;
            try {
                if (opt_key == "extrusion_width") {
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("extrusion_spacing");
                    if (width_option) {
                            width_option->set_phony(false);
                            spacing_option->set_phony(true);
                            if (width_option->value == 0)
                                spacing_option->value = 0;
                            else {
                                Flow flow = Flow::new_from_config_width(FlowRole::frPerimeter, width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, max_nozzle_diameter, layer_height_option->value, overlap_ratio, 0);
                                if (flow.width() < flow.height()) flow.with_height(flow.width());
                                spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                            }
                            spacing_option->percent = width_option->percent;
                            something_changed = true;
                    }
                }
                if (opt_key == "first_layer_extrusion_width") {
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("first_layer_extrusion_spacing");
                    if (width_option) {
                            width_option->set_phony(false);
                            spacing_option->set_phony(true);
                            if (width_option->value == 0)
                                spacing_option->value = 0;
                            else {
                                Flow flow = Flow::new_from_config_width(FlowRole::frPerimeter, 
                                    width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, 
                                    max_nozzle_diameter, layer_height_option->value, overlap_ratio, 0);
                                if (flow.width() < flow.height()) flow.with_height(flow.width());
                                spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                            }
                            spacing_option->percent = width_option->percent;
                            something_changed = true;
                    }
                }
                if (opt_key == "first_layer_infill_extrusion_width") {
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("first_layer_infill_extrusion_spacing");
                    if (width_option) {
                            width_option->set_phony(false);
                            spacing_option->set_phony(true);
                            if (width_option->value == 0)
                                spacing_option->value = 0;
                            else {
                                Flow flow = Flow::new_from_config_width(FlowRole::frPerimeter, 
                                    width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, 
                                    max_nozzle_diameter, layer_height_option->value, overlap_ratio, 0);
                                if (flow.width() < flow.height()) flow.with_height(flow.width());
                                spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                            }
                            spacing_option->percent = width_option->percent;
                            something_changed = true;
                    }
                }
                if (opt_key == "perimeter_extrusion_width") {
                    const ConfigOptionPercent* perimeter_overlap_option = find_option<ConfigOptionPercent>("perimeter_overlap", this, config_collection);
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("perimeter_extrusion_spacing");
                    if (width_option && perimeter_overlap_option) {
                        width_option->set_phony(false);
                        spacing_option->set_phony(true);
                        if (width_option->value == 0)
                            spacing_option->value = 0;
                        else {
                            Flow flow = Flow::new_from_config_width(FlowRole::frExternalPerimeter, 
                                width_option->value == 0 ? *default_width_option : *width_option,  *spacing_option, 
                                max_nozzle_diameter, layer_height_option->value, 
                                std::min(overlap_ratio, (float)perimeter_overlap_option->get_abs_value(1)), 0);
                            if (flow.width() < flow.height()) flow = flow.with_height(flow.width());
                            spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                        }
                        spacing_option->percent = width_option->percent;
                        something_changed = true;
                    }
                }
                if (opt_key == "external_perimeter_extrusion_width") {
                    const ConfigOptionPercent* external_perimeter_overlap_option = find_option<ConfigOptionPercent>("external_perimeter_overlap", this, config_collection);
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("external_perimeter_extrusion_spacing");
                    if (width_option && external_perimeter_overlap_option) {
                        width_option->set_phony(false);
                        spacing_option->set_phony(true);
                        if (width_option->value == 0)
                            spacing_option->value = 0;
                        else {
                            Flow ext_perimeter_flow = Flow::new_from_config_width(FlowRole::frPerimeter, 
                                width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, 
                                max_nozzle_diameter, layer_height_option->value, 
                                std::min(overlap_ratio * 0.5f, float(external_perimeter_overlap_option->get_abs_value(0.5))), 0);
                            if (ext_perimeter_flow.width() < ext_perimeter_flow.height()) ext_perimeter_flow = ext_perimeter_flow.with_height(ext_perimeter_flow.width());
                            spacing_option->value = (width_option->percent) ? std::round(100 * ext_perimeter_flow.spacing() / max_nozzle_diameter) : (std::round(ext_perimeter_flow.spacing() * 10000) / 10000);
                        }
                        spacing_option->percent = width_option->percent;
                        something_changed = true;
                    }
                }
                if (opt_key == "infill_extrusion_width") {
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("infill_extrusion_spacing");
                    if (width_option) {
                        width_option->set_phony(false);
                        spacing_option->set_phony(true);
                        if (width_option->value == 0)
                            spacing_option->value = 0;
                        else {
                            Flow flow = Flow::new_from_config_width(FlowRole::frInfill, width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, max_nozzle_diameter, layer_height_option->value, overlap_ratio, 0);
                            if (flow.width() < flow.height()) flow = flow.with_height(flow.width());
                            spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                        }
                        spacing_option->percent = width_option->percent;
                        something_changed = true;
                    }
                }
                if (opt_key == "solid_infill_extrusion_width") {
                    const ConfigOptionPercent* solid_infill_overlap_option = find_option<ConfigOptionPercent>("solid_infill_overlap", this, config_collection);
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("solid_infill_extrusion_spacing");
                    if (width_option) {
                        width_option->set_phony(false);
                        spacing_option->set_phony(true);
                        if (width_option->value == 0)
                            spacing_option->value = 0;
                        else {
                            Flow flow = Flow::new_from_config_width(FlowRole::frSolidInfill, 
                                width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, 
                                max_nozzle_diameter, layer_height_option->value, 
                                std::min(overlap_ratio, float(solid_infill_overlap_option->get_abs_value(1.))), 0);
                            if (flow.width() < flow.height()) flow = flow.with_height(flow.width());
                            spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                        }
                        spacing_option->percent = width_option->percent;
                        something_changed = true;
                    }
                }
                if (opt_key == "top_infill_extrusion_width") {
                    const ConfigOptionPercent* top_solid_infill_overlap_option = find_option<ConfigOptionPercent>("top_solid_infill_overlap", this, config_collection);
                    spacing_option = this->option<ConfigOptionFloatOrPercent>("top_infill_extrusion_spacing");
                    if (width_option) {
                        width_option->set_phony(false);
                        spacing_option->set_phony(true);
                        if (width_option->value == 0)
                            spacing_option->value = 0;
                        else {
                            Flow flow = Flow::new_from_config_width(FlowRole::frTopSolidInfill, 
                                width_option->value == 0 ? *default_width_option : *width_option, *spacing_option, 
                                max_nozzle_diameter, layer_height_option->value,
                                std::min(overlap_ratio, float(top_solid_infill_overlap_option->get_abs_value(1.))), 0);
                            if (flow.width() < flow.height()) flow = flow.with_height(flow.width());
                            spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                        }
                        spacing_option->percent = width_option->percent;
                        something_changed = true;
                    }
                }
                //if (opt_key == "support_material_extrusion_width") {
                //    Flow flow = Flow::new_from_config_width(FlowRole::frSupportMaterial, width_option->value == 0 ? *default_width_option : *width_option, max_nozzle_diameter, layer_height_option->value, 0);
                //    if (width_option->percent)
                //        this->set_key_value("support_material_extrusion_spacing", new ConfigOptionFloatOrPercent(std::round(100 * flow.spacing() / max_nozzle_diameter), true));
                //    else
                //        this->set_key_value("support_material_extrusion_spacing", new ConfigOptionFloatOrPercent(std::round(flow.spacing() * 10000) / 10000, false));
                //    something_changed = true;
                //}
                //if (opt_key == "skirt_extrusion_width") {
                //    Flow flow = Flow::new_from_config_width(FlowRole::frPerimeter, width_option->value == 0 ? *default_width_option : *width_option, max_nozzle_diameter, layer_height_option->value, 0);
                //    if (width_option->percent)
                //        this->set_key_value("skirt_extrusion_spacing", new ConfigOptionFloatOrPercent(std::round(100 * flow.spacing() / max_nozzle_diameter), true));
                //    else
                //        this->set_key_value("skirt_extrusion_spacing", new ConfigOptionFloatOrPercent(std::round(flow.spacing() * 10000) / 10000, false));
                //    something_changed = true;
                //}
            } catch (FlowErrorNegativeSpacing) {
                if (spacing_option != nullptr) {
                    width_option->set_phony(true);
                    spacing_option->set_phony(false);
                    spacing_option->value = 100;
                    spacing_option->percent = true;
                    Flow flow = Flow::new_from_spacing(spacing_option->get_abs_value(max_nozzle_diameter), max_nozzle_diameter, layer_height_option->value, overlap_ratio, false);
                    width_option->value = (spacing_option->percent) ? std::round(100 * flow.width() / max_nozzle_diameter) : (std::round(flow.width() * 10000) / 10000);
                    width_option->percent = spacing_option->percent;
                    something_changed = true;
                } else {
                    width_option->value = 100;
                    width_option->percent = true;
                    width_option->set_phony(false);
                    spacing_option->set_phony(true);
                    Flow flow = Flow::new_from_config_width(FlowRole::frPerimeter, width_option->value == 0 ? *width_option : *default_width_option, *spacing_option, max_nozzle_diameter, layer_height_option->value, overlap_ratio, 0);
                    spacing_option->value = (width_option->percent) ? std::round(100 * flow.spacing() / max_nozzle_diameter) : (std::round(flow.spacing() * 10000) / 10000);
                    spacing_option->percent = width_option->percent;
                    something_changed = true;
                }
            }
        }
    }
    //update phony counterpark of 0-set fields
    // now they show 0, no need to update them
    //if (opt_key == "extrusion_width" || opt_key == "extrusion_spacing") {
    //    for (auto conf : config_collection) {
    //        if (conf->option("extrusion_width"))
    //            if (!conf->update_phony(config_collection, true).empty())
    //                return { conf };
    //    }
    //    return {};
    //}
    if(something_changed)
        return this;
    return nullptr;
}

//FIXME localize this function.
//note: seems only called for config export & command line. Most of the validation work for the gui is done elsewhere... So this function may be a bit out-of-sync
std::string validate(const FullPrintConfig& cfg)
{
    // --layer-height
    if (cfg.get_computed_value("layer_height") <= 0)
        return "Invalid value for --layer-height";
    if (fabs(fmod(cfg.get_computed_value("layer_height"), SCALING_FACTOR)) > 1e-4)
        return "--layer-height must be a multiple of print resolution";

    // --first-layer-height
    //if (cfg.get_abs_value("first_layer_height") <= 0) //can't do that, as the extruder isn't defined
    if(cfg.first_layer_height.value <= 0)
        return "Invalid value for --first-layer-height";

    // --filament-diameter
    for (double fd : cfg.filament_diameter.get_values())
        if (fd < 1)
            return "Invalid value for --filament-diameter";

    // --nozzle-diameter
    for (double nd : cfg.nozzle_diameter.get_values())
        if (nd < 0.005)
            return "Invalid value for --nozzle-diameter";

    // --perimeters
    if (cfg.perimeters.value < 0)
        return "Invalid value for --perimeters";

    // --solid-layers
    if (cfg.top_solid_layers < 0)
        return "Invalid value for --top-solid-layers";
    if (cfg.bottom_solid_layers < 0)
        return "Invalid value for --bottom-solid-layers";

    if (cfg.use_firmware_retraction.value &&
        cfg.gcode_flavor.value != gcfSmoothie &&
        cfg.gcode_flavor.value != gcfSprinter &&
        cfg.gcode_flavor.value != gcfRepRap &&
        cfg.gcode_flavor.value != gcfMarlinLegacy &&
        cfg.gcode_flavor.value != gcfMarlinFirmware &&
        cfg.gcode_flavor.value != gcfMachinekit &&
        cfg.gcode_flavor.value != gcfRepetier &&
        cfg.gcode_flavor.value != gcfKlipper)
        return "--use-firmware-retraction is only supported by Marlin 1&2, Smoothie, Sprinter, Reprap, Repetier, Machinekit, Repetier, Klipper? and Lerdge firmware";

    if (cfg.use_firmware_retraction.value)
        for (unsigned char wipe : cfg.wipe.get_values())
             if (wipe)
                return "--use-firmware-retraction is not compatible with --wipe";

    // --gcode-flavor
    if (! print_config_def.get("gcode_flavor")->has_enum_value(cfg.gcode_flavor.serialize()))
        return "Invalid value for --gcode-flavor";

    // --fill-pattern
    if (! print_config_def.get("fill_pattern")->has_enum_value(cfg.fill_pattern.serialize()))
        return "Invalid value for --fill-pattern";

    // --top-fill-pattern
    if (!print_config_def.get("top_fill_pattern")->has_enum_value(cfg.top_fill_pattern.serialize()))
        return "Invalid value for --top-fill-pattern";

    // --bottom-fill-pattern
    if (! print_config_def.get("bottom_fill_pattern")->has_enum_value(cfg.bottom_fill_pattern.serialize()))
        return "Invalid value for --bottom-fill-pattern";

    // --solid-fill-pattern
    if (!print_config_def.get("solid_fill_pattern")->has_enum_value(cfg.solid_fill_pattern.serialize()))
        return "Invalid value for --solid-fill-pattern";

    // --brim-ears-pattern
    if (!print_config_def.get("brim_ears_pattern")->has_enum_value(cfg.brim_ears_pattern.serialize()))
        return "Invalid value for --brim-ears-pattern";

    // --fill-density
    if (fabs(cfg.fill_density.value - 100.) < EPSILON &&
        (! print_config_def.get("top_fill_pattern")->has_enum_value(cfg.fill_pattern.serialize())
        && ! print_config_def.get("bottom_fill_pattern")->has_enum_value(cfg.fill_pattern.serialize())
        ))
        return "The selected fill pattern is not supposed to work at 100% density";

    // --infill-every-layers
    if (cfg.infill_every_layers < 1)
        return "Invalid value for --infill-every-layers";

    // --skirt-height
    if (cfg.skirt_height < 0)
        return "Invalid value for --skirt-height";
    
    // extruder clearance
    if (cfg.extruder_clearance_radius <= 0)
        return "Invalid value for --extruder-clearance-radius";
    if (cfg.extruder_clearance_height <= 0)
        return "Invalid value for --extruder-clearance-height";

    // --extrusion-multiplier
    for (double em : cfg.extrusion_multiplier.get_values())
        if (em <= 0)
            return "Invalid value for --extrusion-multiplier";

    // --spiral-vase
    if (cfg.spiral_vase) {
        // Note that we might want to have more than one perimeter on the bottom
        // solid layers.
        if (cfg.perimeters > 1)
            return "Can't make more than one perimeter when spiral vase mode is enabled";
        else if (cfg.perimeters < 1)
            return "Can't make less than one perimeter when spiral vase mode is enabled";
        if (cfg.fill_density > 0)
            return "Spiral vase mode can only print hollow objects, so you need to set Fill density to 0";
        if (cfg.top_solid_layers > 0)
            return "Spiral vase mode is not compatible with top solid layers";
        if (cfg.support_material || cfg.support_material_enforce_layers > 0)
            return "Spiral vase mode is not compatible with support material";
        if (cfg.infill_dense)
            return "Spiral vase mode can only print hollow objects and have no top surface, so you don't need any dense infill";
        if (cfg.extra_perimeters || cfg.extra_perimeters_on_overhangs || cfg.extra_perimeters_odd_layers)
            return "Can't make more than one perimeter when spiral vase mode is enabled";
        if (cfg.overhangs_reverse)
            return "Can't reverse the direction of the overhangs every layer when spiral vase mode is enabled";
        if (cfg.perimeter_reverse)
            return "Can't reverse the direction of the perimeters every layer when spiral vase mode is enabled";
    }

    // extrusion widths
    {
        double max_nozzle_diameter = 0.;
        for (double dmr : cfg.nozzle_diameter.get_values())
            max_nozzle_diameter = std::max(max_nozzle_diameter, dmr);
        const char *widths[] = { "", "external_perimeter_", "perimeter_", "infill_", "solid_infill_", "top_infill_", "support_material_", "first_layer_", "first_layer_infill_", "skirt_" };
        for (size_t i = 0; i < sizeof(widths) / sizeof(widths[i]); ++ i) {
            std::string key(widths[i]);
            key += "extrusion_width";
            if (cfg.get_abs_value(key, max_nozzle_diameter) > 10. * max_nozzle_diameter)
                return std::string("Invalid extrusion width (too large): ") + key;
        }
    }

    // Out of range validation of numeric values.
    for (const std::string &opt_key : cfg.keys()) {
        const ConfigOption      *opt    = cfg.optptr(opt_key);
        assert(opt != nullptr);
        const ConfigOptionDef   *optdef = print_config_def.get(opt_key);
        assert(optdef != nullptr);

        if (!opt->is_enabled()) {
            // Do not check disabled values
            continue;
        }

        bool out_of_range = false;
        switch (opt->type()) {
        case coFloat:
        case coPercent:
        {
            auto *fopt = static_cast<const ConfigOptionFloat*>(opt);
            out_of_range = fopt->value < optdef->min || fopt->value > optdef->max;
            break;
        }
        case coFloatOrPercent:
        {
            auto *fopt = static_cast<const ConfigOptionFloatOrPercent*>(opt);
            out_of_range = fopt->get_abs_value(1) < optdef->min || fopt->get_abs_value(1) > optdef->max;
            break;
        }
        case coPercents:
        case coFloats:
        {
            const auto* vec = static_cast<const ConfigOptionVector<double>*>(opt);
            for (size_t i = 0; i < vec->size(); ++i) {
                if (!vec->is_enabled(i))
                    continue;
                double v = vec->get_at(i);
                if (v < optdef->min || v > optdef->max) {
                    out_of_range = true;
                    break;
                }
            }
            break;
        }
        case coFloatsOrPercents:
        {
            const auto* vec = static_cast<const ConfigOptionVector<FloatOrPercent>*>(opt);
            for (size_t i = 0; i < vec->size(); ++i) {
                if (!vec->is_enabled(i))
                    continue;
                const FloatOrPercent &v = vec->get_at(i);
                if (v.value < optdef->min || v.value > optdef->max) {
                    out_of_range = true;
                    break;
                }
            }
            break;
        }
        case coInt:
        {
            auto *iopt = static_cast<const ConfigOptionInt*>(opt);
            out_of_range = iopt->value < optdef->min || iopt->value > optdef->max;
            break;
        }
        case coInts:
        {
            const auto* vec = static_cast<const ConfigOptionVector<int32_t>*>(opt);
            for (size_t i = 0; i < vec->size(); ++i) {
                if (!vec->is_enabled(i))
                    continue;
                int v = vec->get_at(i);
                if (v < optdef->min || v > optdef->max) {
                    out_of_range = true;
                    break;
                }
            }
            break;
        }
        default:;
        }
        if (out_of_range)
            return std::string("Value out of range: " + opt_key);
    }

    // The configuration is valid.
    return "";
}

// Declare and initialize static caches of StaticPrintConfig derived classes.
#define PRINT_CONFIG_CACHE_ELEMENT_DEFINITION(r, data, CLASS_NAME) StaticPrintConfig::StaticCache<class Slic3r::CLASS_NAME> BOOST_PP_CAT(CLASS_NAME::s_cache_, CLASS_NAME);
#define PRINT_CONFIG_CACHE_ELEMENT_INITIALIZATION(r, data, CLASS_NAME) Slic3r::CLASS_NAME::initialize_cache();
#define PRINT_CONFIG_CACHE_INITIALIZE(CLASSES_SEQ) \
    BOOST_PP_SEQ_FOR_EACH(PRINT_CONFIG_CACHE_ELEMENT_DEFINITION, _, BOOST_PP_TUPLE_TO_SEQ(CLASSES_SEQ)) \
    int print_config_static_initializer() { \
        /* Putting a trace here to avoid the compiler to optimize out this function. */ \
        /*BOOST_LOG_TRIVIAL(trace) << "Initializing StaticPrintConfigs";*/ \
        /* Tamas: alternative solution through a static volatile int. Boost log pollutes stdout and prevents tests from generating clean output */ \
        static volatile int ret = 1; \
        BOOST_PP_SEQ_FOR_EACH(PRINT_CONFIG_CACHE_ELEMENT_INITIALIZATION, _, BOOST_PP_TUPLE_TO_SEQ(CLASSES_SEQ)) \
        return ret; \
    }
PRINT_CONFIG_CACHE_INITIALIZE((
    PrintObjectConfig, PrintRegionConfig, MachineEnvelopeConfig, GCodeConfig, PrintConfig, FullPrintConfig, 
    SLAMaterialConfig, SLAPrintConfig, SLAPrintObjectConfig, SLAPrinterConfig, SLAFullPrintConfig))
static int print_config_static_initialized = print_config_static_initializer();

CLIActionsConfigDef::CLIActionsConfigDef()
{
    ConfigOptionDef* def;

    // Actions:
    def = this->add("export_obj", coBool);
    def->label = L("Export OBJ");
    def->tooltip = L("Export the model(s) as OBJ.");
    def->set_default_value(new ConfigOptionBool(false));

/*
    def = this->add("export_svg", coBool);
    def->label = L("Export SVG");
    def->tooltip = L("Slice the model and export solid slices as SVG.");
    def->set_default_value(new ConfigOptionBool(false));
*/

    def = this->add("export_sla", coBool);
    def->label = L("Export SLA");
    def->tooltip = L("Slice the model and export SLA printing layers as PNG.");
    def->cli = "export-sla|sla";
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("export_3mf", coBool);
    def->label = L("Export 3MF");
    def->tooltip = L("Export the model(s) as 3MF.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("export_amf", coBool);
    def->label = L("Export AMF");
    def->tooltip = L("Export the model(s) as AMF.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("export_stl", coBool);
    def->label = L("Export STL");
    def->tooltip = L("Export the model(s) as STL.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("export_gcode", coBool);
    def->label = L("Export G-code");
    def->tooltip = L("Slice the model and export toolpaths as G-code.");
    def->cli = "export-gcode|gcode|g";
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("gcodeviewer", coBool);
    def->label = L("G-code viewer");
    def->tooltip = L("Visualize an already sliced and saved G-code");
    def->cli = "gcodeviewer";
    def->set_default_value(new ConfigOptionBool(false));

#if ENABLE_GL_CORE_PROFILE
    def = this->add("opengl-version", coString);
    def->label = L("OpenGL version");
    def->tooltip = L("Select a specific version of OpenGL");
    def->cli = "opengl-version";
    def->set_default_value(new ConfigOptionString());

    def = this->add("opengl-compatibility", coBool);
    def->label = L("OpenGL compatibility profile");
    def->tooltip = L("Enable OpenGL compatibility profile");
    def->cli = "opengl-compatibility";
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("opengl-debug", coBool);
    def->label = L("OpenGL debug output");
    def->tooltip = L("Activate OpenGL debug output on graphic cards which support it (OpenGL 4.3 or higher)");
    def->cli = "opengl-debug";
    def->set_default_value(new ConfigOptionBool(false));
#endif // ENABLE_GL_CORE_PROFILE

    def = this->add("slice", coBool);
    def->label = L("Slice");
    def->tooltip = L("Slice the model as FFF or SLA based on the printer_technology configuration value.");
    def->cli = "slice|s";
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("help", coBool);
    def->label = L("Help");
    def->tooltip = L("Show this help.");
    def->cli = "help|h";
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("help_fff", coBool);
    def->label = L("Help (FFF options)");
    def->tooltip = L("Show the full list of print/G-code configuration options.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("help_sla", coBool);
    def->label = L("Help (SLA options)");
    def->tooltip = L("Show the full list of SLA print configuration options.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("info", coBool);
    def->label = L("Output Model Info");
    def->tooltip = L("Write information about the model to the console.");
    def->set_default_value(new ConfigOptionBool(false));

    def = this->add("save", coString);
    def->label = L("Save config file");
    def->tooltip = L("Save configuration to the specified file.");
    def->set_default_value(new ConfigOptionString());
}

CLITransformConfigDef::CLITransformConfigDef()
{
    ConfigOptionDef* def;

    // Transform options:
    def = this->add("align_xy", coPoint);
    def->label = L("Align XY");
    def->tooltip = L("Align the model to the given point.");
    def->set_default_value(new ConfigOptionPoint(Vec2d(100,100)));

    def = this->add("cut", coFloat);
    def->label = L("Cut");
    def->tooltip = L("Cut model at the given Z.");
    def->set_default_value(new ConfigOptionFloat(0));

/*
    def = this->add("cut_grid", coFloat);
    def->label = L("Cut");
    def->tooltip = L("Cut model in the XY plane into tiles of the specified max size.");
    def->set_default_value(new ConfigOptionPoint());

    def = this->add("cut_x", coFloat);
    def->label = L("Cut");
    def->tooltip = L("Cut model at the given X.");
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("cut_y", coFloat);
    def->label = L("Cut");
    def->tooltip = L("Cut model at the given Y.");
    def->set_default_value(new ConfigOptionFloat(0));
*/

    def = this->add("center", coPoint);
    def->label = L("Center");
    def->tooltip = L("Center the print around the given center.");
    def->set_default_value(new ConfigOptionPoint(Vec2d(100,100)));

    def = this->add("dont_arrange", coBool);
    def->label = L("Don't arrange");
    def->tooltip = L("Do not rearrange the given models before merging and keep their original XY coordinates.");

    def = this->add("ensure_on_bed", coBool);
    def->label = L("Ensure on bed");
    def->tooltip = L("Lift the object above the bed when it is partially below. Enabled by default, use --no-ensure-on-bed to disable.");
    def->set_default_value(new ConfigOptionBool(true));

    def = this->add("duplicate", coInt);
    def->label = L("Duplicate");
    def->tooltip =L("Multiply copies by this factor.");
    def->min = 1;

    def = this->add("duplicate_grid", coPoint);
    def->label = L("Duplicate by grid");
    def->tooltip = L("Multiply copies by creating a grid.");

    def = this->add("merge", coBool);
    def->label = L("Merge");
    def->tooltip = L("Arrange the supplied models in a plate and merge them in a single model in order to perform actions once.");
    def->cli = "merge|m";

    def = this->add("repair", coBool);
    def->label = L("Repair");
    def->tooltip = L("Try to repair any non-manifold meshes (this option is implicitly added whenever we need to slice the model to perform the requested action).");

    def = this->add("rotate", coFloat);
    def->label = L("Rotate");
    def->tooltip = L("Rotation angle around the Z axis in degrees.");
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("rotate_x", coFloat);
    def->label = L("Rotate around X");
    def->tooltip = L("Rotation angle around the X axis in degrees.");
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("rotate_y", coFloat);
    def->label = L("Rotate around Y");
    def->tooltip = L("Rotation angle around the Y axis in degrees.");
    def->set_default_value(new ConfigOptionFloat(0));

    def = this->add("scale", coFloatOrPercent);
    def->label = L("Scale");
    def->tooltip = L("Scaling factor or percentage.");
    def->set_default_value(new ConfigOptionFloatOrPercent(1, false));

    def = this->add("split", coBool);
    def->label = L("Split");
    def->tooltip = L("Detect unconnected parts in the given model(s) and split them into separate objects.");

    def = this->add("scale_to_fit", coPoint3);
    def->label = L("Scale to Fit");
    def->tooltip = L("Scale to fit the given volume.");
    def->set_default_value(new ConfigOptionPoint3(Vec3d(0,0,0)));

    def = this->add("delete-after-load", coString);
    def->label = L("Delete files after loading");
    def->tooltip = L("Delete files after loading.");
}

CLIMiscConfigDef::CLIMiscConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("ignore_nonexistent_config", coBool);
    def->label = L("Ignore non-existent config files");
    def->tooltip = L("Do not fail if a file supplied to --load does not exist.");

    def = this->add("config_compatibility", coEnum);
    def->label = L("Forward-compatibility rule when loading configurations from config files and project files (3MF, AMF).");
    def->tooltip = L("This version of Slic3r may not understand configurations produced by the newest Slic3r versions. "
                     "For example, newer Slic3r may extend the list of supported firmware flavors. One may decide to "
                     "bail out or to substitute an unknown value with a default silently or verbosely.");
    def->set_enum<ForwardCompatibilitySubstitutionRule>({
        { "disable",        L("Bail out on unknown configuration values") },
        { "enable",         L("Enable reading unknown configuration values by verbosely substituting them with defaults.") },
        { "enable_silent",  L("Enable reading unknown configuration values by silently substituting them with defaults.") }
    });
    def->set_default_value(new ConfigOptionEnum<ForwardCompatibilitySubstitutionRule>(ForwardCompatibilitySubstitutionRule::Enable));

    def = this->add("load", coStrings);
    def->label = L("Load config file");
    def->tooltip = L("Load configuration from the specified file. It can be used more than once to load options from multiple files.");

    def = this->add("output", coString);
    def->label = L("Output File");
    def->tooltip = L("The file where the output will be written (if not specified, it will be based on the input file).");
    def->cli = "output|o";

    def = this->add("single_instance", coBool);
    def->label = L("Single instance mode");
    def->tooltip = L("If enabled, the command line arguments are sent to an existing instance of GUI Slic3r, "
                     "or an existing Slic3r window is activated. "
                     "Overrides the \"single_instance\" configuration value from application preferences.");

    def = this->add("datadir", coString);
    def->label = L("Data directory");
    def->tooltip = L("Load and store settings at the given directory. This is useful for maintaining different profiles or including configurations from a network storage.");

    def = this->add("threads", coInt);
    def->label = L("Maximum number of threads");
    def->tooltip = L("Sets the maximum number of threads the slicing process will use. If not defined, it will be decided automatically.");
    def->min = 1;

    def = this->add("loglevel", coInt);
    def->label = L("Logging level");
    def->tooltip = L("Sets logging sensitivity. 0:fatal, 1:error, 2:warning, 3:info, 4:debug, 5:trace\n"
                     "For example. loglevel=2 logs fatal, error and warning level messages.");
    def->min = 0;

#if (defined(_MSC_VER) || defined(__MINGW32__)) && defined(SLIC3R_GUI)
    def = this->add("sw_renderer", coBool);
    def->label = L("Render with a software renderer");
    def->tooltip = L("Render with a software renderer. The bundled MESA software renderer is loaded instead of the default OpenGL driver.");
    def->min = 0;
#endif /* _MSC_VER */
}

const CLIActionsConfigDef    cli_actions_config_def;
const CLITransformConfigDef  cli_transform_config_def;
const CLIMiscConfigDef       cli_misc_config_def;

DynamicPrintAndCLIConfig::PrintAndCLIConfigDef DynamicPrintAndCLIConfig::s_def;

void DynamicPrintAndCLIConfig::handle_legacy(t_config_option_key &opt_key, std::string &value) const
{
    if (cli_actions_config_def  .options.find(opt_key) == cli_actions_config_def  .options.end() &&
        cli_transform_config_def.options.find(opt_key) == cli_transform_config_def.options.end() &&
        cli_misc_config_def     .options.find(opt_key) == cli_misc_config_def     .options.end()) {
        PrintConfigDef::handle_legacy(opt_key, value);
    }
}

// SlicingStatesConfigDefs

ReadOnlySlicingStatesConfigDef::ReadOnlySlicingStatesConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("zhop", coFloat);
    def->label = L("Current z-hop");
    def->tooltip = L("Contains z-hop present at the beginning of the custom G-code block.");
}

ReadWriteSlicingStatesConfigDef::ReadWriteSlicingStatesConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("position", coFloats);
    def->label = L("Position");
    def->tooltip = L("Position of the extruder at the beginning of the custom G-code block. If the custom G-code travels somewhere else, "
                     "it should write to this variable so PrusaSlicer knows where it travels from when it gets control back.");

    def = this->add("e_retracted", coFloats);
    def->label = L("Retraction");
    def->tooltip = L("Retraction state at the beginning of the custom G-code block. If the custom G-code moves the extruder axis, "
                     "it should write to this variable so PrusaSlicer deretracts correctly when it gets control back.");

    def = this->add("e_restart_extra", coFloats);
    def->label = L("Extra deretraction");
    def->tooltip = L("Currently planned extra extruder priming after deretraction.");

    def = this->add("e_position", coFloats);
    def->label = L("Absolute E position");
    def->tooltip = L("Current position of the extruder axis. Only used with absolute extruder addressing.");
}

OtherSlicingStatesConfigDef::OtherSlicingStatesConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("current_extruder", coInt);
    def->label = L("Current extruder");
    def->tooltip = L("Zero-based index of currently used extruder.");

    def = this->add("current_object_idx", coInt);
    def->label = L("Current object index");
    def->tooltip = L("Specific for sequential printing. Zero-based index of currently printed object.");

    def = this->add("has_single_extruder_multi_material_priming", coBool);
    def->label = L("Has single extruder MM priming");
    def->tooltip = L("Are the extra multi-material priming regions used in this print?");

    def = this->add("has_wipe_tower", coBool);
    def->label = L("Has wipe tower");
    def->tooltip = L("Whether or not wipe tower is being generated in the print.");

    def = this->add("initial_extruder", coInt);
    def->label = L("Initial extruder");
    def->tooltip = L("Zero-based index of the first extruder used in the print. Same as initial_tool.");

    def = this->add("initial_filament_type", coString);
    // TRN: Meaning 'filament type of the initial filament'
    def->label = L("Initial filament type");
    def->tooltip = L("String containing filament type of the first used extruder.");

    def = this->add("initial_tool", coInt);
    def->label = L("Initial tool");
    def->tooltip = L("Zero-based index of the first extruder used in the print. Same as initial_extruder.");

    def = this->add("is_extruder_used", coBools);
    def->label = L("Is extruder used?");
    def->tooltip = L("Vector of booleans stating whether a given extruder is used in the print.");
}

PrintStatisticsConfigDef::PrintStatisticsConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("extruded_volume", coFloats);
    def->label = L("Volume per extruder");
    def->tooltip = L("Total filament volume extruded per extruder during the entire print.");

    def = this->add("normal_print_time", coString);
    def->label = L("Print time (normal mode)");
    def->tooltip = L("Estimated print time when printed in normal mode (i.e. not in silent mode). Same as print_time.");
    
    def = this->add("num_printing_extruders", coInt);
    def->label = L("Number of printing extruders");
    def->tooltip = L("Number of extruders used during the print.");

    def = this->add("print_time", coString);
    def->label = L("Print time (normal mode)");
    def->tooltip = L("Estimated print time when printed in normal mode (i.e. not in silent mode). Same as normal_print_time.");

    def = this->add("printing_filament_types", coString);
    def->label = L("Used filament types");
    def->tooltip = L("Comma-separated list of all filament types used during the print.");

    def = this->add("silent_print_time", coString);
    def->label = L("Print time (silent mode)");
    def->tooltip = L("Estimated print time when printed in silent mode.");

    def = this->add("total_cost", coFloat);
    def->label = L("Total cost");
    def->tooltip = L("Total cost of all material used in the print. Calculated from cost in Filament Settings.");

    def = this->add("total_weight", coFloat);
    def->label = L("Total weight");
    def->tooltip = L("Total weight of the print. Calculated from density in Filament Settings.");

    def = this->add("total_wipe_tower_cost", coFloat);
    def->label = L("Total wipe tower cost");
    def->tooltip = L("Total cost of the material wasted on the wipe tower. Calculated from cost in Filament Settings.");

    def = this->add("total_wipe_tower_filament", coFloat);
    def->label = L("Wipe tower volume");
    def->tooltip = L("Total filament volume extruded on the wipe tower.");

    def = this->add("used_filament", coFloat);
    def->label = L("Used filament");
    def->tooltip = L("Total length of filament used in the print.");

    def = this->add("total_toolchanges", coInt);
    def->label = L("Total number of toolchanges");
    def->tooltip = L("Number of toolchanges during the print.");

    def = this->add("extruded_volume_total", coFloat);
    def->label = L("Total volume");
    def->tooltip = L("Total volume of filament used during the entire print.");

    def = this->add("extruded_weight", coFloats);
    def->label = L("Weight per extruder");
    def->tooltip = L("Weight per extruder extruded during the entire print. Calculated from density in Filament Settings.");

    def = this->add("extruded_weight_total", coFloat);
    def->label = L("Total weight");
    def->tooltip = L("Total weight of the print. Calculated from density in Filament Settings.");

    def = this->add("total_layer_count", coInt);
    def->label = L("Total layer count");
    def->tooltip = L("Number of layers in the entire print.");
}

ObjectsInfoConfigDef::ObjectsInfoConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("num_objects", coInt);
    def->label = L("Number of objects");
    def->tooltip = L("Total number of objects in the print.");

    def = this->add("num_instances", coInt);
    def->label = L("Number of instances");
    def->tooltip = L("Total number of object instances in the print, summed over all objects.");

    def = this->add("scale", coStrings);
    def->label = L("Scale per object");
    def->tooltip = L("Contains a string with the information about what scaling was applied to the individual objects. "
                     "Indexing of the objects is zero-based (first object has index 0).\n"
                     "Example: 'x:100% y:50% z:100%'.");

    def = this->add("input_filename_base", coString);
    def->label = L("Input filename without extension");
    def->tooltip = L("Source filename of the first object, without extension.");
}

DimensionsConfigDef::DimensionsConfigDef()
{
    ConfigOptionDef* def;

    const std::string point_tooltip   = L("The vector has two elements: x and y coordinate of the point. Values in mm.");
    const std::string bb_size_tooltip = L("The vector has two elements: x and y dimension of the bounding box. Values in mm.");

    def = this->add("first_layer_print_convex_hull", coPoints);
    def->label = L("First layer convex hull");
    def->tooltip = L("Vector of points of the first layer convex hull. Each element has the following format: "
                     "'[x, y]' (x and y are floating-point numbers in mm).");

    def = this->add("first_layer_print_min", coFloats);
    def->label = L("Bottom-left corner of first layer bounding box");
    def->tooltip = point_tooltip;

    def = this->add("first_layer_print_max", coFloats);
    def->label = L("Top-right corner of first layer bounding box");
    def->tooltip = point_tooltip;

    def = this->add("first_layer_print_size", coFloats);
    def->label = L("Size of the first layer bounding box");
    def->tooltip = bb_size_tooltip;

    def = this->add("print_bed_min", coFloats);
    def->label = L("Bottom-left corner of print bed bounding box");
    def->tooltip = point_tooltip;

    def = this->add("print_bed_max", coFloats);
    def->label = L("Top-right corner of print bed bounding box");
    def->tooltip = point_tooltip;

    def = this->add("print_bed_size", coFloats);
    def->label = L("Size of the print bed bounding box");
    def->tooltip = bb_size_tooltip;
}

TimestampsConfigDef::TimestampsConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("timestamp", coString);
    def->label = L("Timestamp");
    def->tooltip = L("String containing current time in yyyyMMdd-hhmmss format.");

    def = this->add("year", coInt);
    def->label = L("Year");

    def = this->add("month", coInt);
    def->label = L("Month");

    def = this->add("day", coInt);
    def->label = L("Day");

    def = this->add("hour", coInt);
    def->label = L("Hour");

    def = this->add("minute", coInt);
    def->label = L("Minute");

    def = this->add("second", coInt);
    def->label = L("Second");
}

OtherPresetsConfigDef::OtherPresetsConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("num_extruders", coInt);
    def->label = L("Number of extruders");
    def->tooltip = L("Total number of extruders, regardless of whether they are used in the current print.");

    def = this->add("num_milling", coInt);
    def->label = L("Number of mills");
    def->tooltip = L("Total number of mills, regardless of whether they are used in the current print.");

    def = this->add("print_preset", coString);
    def->label = L("Print preset name");
    def->tooltip = L("Name of the print preset used for slicing.");

    def = this->add("filament_preset", coStrings);
    def->label = L("Filament preset name");
    def->tooltip = L("Names of the filament presets used for slicing. The variable is a vector "
                     "containing one name for each extruder.");
    def->is_vector_extruder = true;

    def = this->add("printer_preset", coString);
    def->label = L("Printer preset name");
    def->tooltip = L("Name of the printer preset used for slicing.");

    def = this->add("physical_printer_preset", coString);
    def->label = L("Physical printer name");
    def->tooltip = L("Name of the physical printer used for slicing.");
}


static std::map<t_custom_gcode_key, t_config_option_keys> s_CustomGcodeSpecificPlaceholders{
    {"start_filament_gcode",    {"layer_num", "layer_z", "max_layer_z", "filament_extruder_id", "previous_extruder", "next_extruder"}},
    {"end_filament_gcode",      {"layer_num", "layer_z", "max_layer_z", "filament_extruder_id", "previous_extruder", "next_extruder"}},
    {"milling_toolchange_start_gcode", {"layer_num", "layer_z", "previous_layer_z", "max_layer_z", "previous_extruder", "next_extruder"}},
    {"milling_toolchange_end_gcode",   {"layer_num", "layer_z", "previous_layer_z", "max_layer_z", "previous_extruder", "next_extruder"}},
    {"end_gcode",               {"layer_num", "layer_z", "max_layer_z", "filament_extruder_id", "previous_extruder", "next_extruder"}},
    {"before_layer_gcode",      {"layer_num", "layer_z", "previous_layer_z", "max_layer_z"}},
    {"layer_gcode",             {"layer_num", "layer_z", "previous_layer_z", "max_layer_z"}},
    {"feature_gcode",           {"layer_num", "layer_z", "max_layer_z", "previous_extrusion_role", "next_extrusion_role", /*deprecated*/"extrusion_role", "last_extrusion_role" /*deprecated*/}},
    {"toolchange_gcode",        {"layer_num", "layer_z", "max_layer_z", "previous_extruder", "next_extruder", "toolchange_z"}},
    {"color_change_gcode",      {"color_change_extruder"}},
    {"pause_print_gcode",       {"color_change_extruder"}},
    {"between_objects_gcode",   {"layer_num", "layer_z"}},
};

const std::map<t_custom_gcode_key, t_config_option_keys>& custom_gcode_specific_placeholders()
{
    return s_CustomGcodeSpecificPlaceholders;
}

CustomGcodeSpecificConfigDef::CustomGcodeSpecificConfigDef()
{
    ConfigOptionDef* def;

    def = this->add("layer_num", coInt);
    def->label = L("Layer number");
    def->tooltip = L("Zero-based index of the current layer (i.e. first layer is number 0).");

    def = this->add("layer_z", coFloat);
    def->label = L("Layer Z");
    def->tooltip = L("Height of the current layer above the print bed, measured to the top of the layer.");

    def = this->add("previous_layer_z", coFloat);
    def->label = L("Previous Layer Z");
    def->tooltip = L("Height of the previous layer.");

    def = this->add("max_layer_z", coFloat);
    def->label = L("Maximal layer Z");
    def->tooltip = L("Height of the last layer above the print bed.");

    def = this->add("filament_extruder_id", coInt);
    def->label = L("Current extruder index");
    def->tooltip = L("Zero-based index of currently used extruder (i.e. first extruder has index 0)."
        "\nWith multiple extruders, an extra 'end_filament_gcode' is applied at the end of the print for each extruder. In this case, 'filament_extruder_id' will be the index of the extruder to 'finalize'.");

    def = this->add("previous_extruder", coInt);
    def->label = L("Previous extruder");
    def->tooltip = L("Index of the extruder that is being unloaded. The index is zero based (first extruder has index 0). -1 if there is no extruder before (like in the start gcode)."
        "\nWith multiple extruders, an extra 'end_filament_gcode' is applied at the end of the print for each extruder. In this case, 'previous_extruder' will be the index of the last extruder used.");

    def = this->add("next_extruder", coInt);
    def->label = L("Next extruder");
    def->tooltip = L("Index of the extruder that is being loaded. The index is zero based (first extruder has index 0). -1 if there is no extruder after (like in the end gcode)."
        "\nWith multiple extruders, an extra 'end_filament_gcode' is applied at the end of the print for each extruder. In this case, 'next_extruder' will be -1.");

    def = this->add("toolchange_z", coFloat);
    def->label = L("Toolchange Z");
    def->tooltip = L("Height above the print bed when the toolchange takes place. Usually the same as layer_z, but can be different.");

    def = this->add("color_change_extruder", coInt);
    // TRN: This is a label in custom g-code editor dialog, belonging to color_change_extruder. Denoted index of the extruder for which color change is performed.
    def->label = L("Color change extruder");
    def->tooltip = L("Index of the extruder for which color change will be performed. The index is zero based (first extruder has index 0).");
}

const CustomGcodeSpecificConfigDef custom_gcode_specific_config_def;

uint64_t ModelConfig::s_last_timestamp = 1;

static Points to_points(const std::vector<Vec2d> &dpts)
{
    Points pts; pts.reserve(dpts.size());
    for (auto &v : dpts)
        pts.emplace_back( coord_t(scale_(v.x())), coord_t(scale_(v.y())) );
    return pts;    
}

Points get_bed_shape(const DynamicPrintConfig &config)
{
    const auto *bed_shape_opt = config.opt<ConfigOptionPoints>("bed_shape");
    if (!bed_shape_opt) {
        
        // Here, it is certain that the bed shape is missing, so an infinite one
        // has to be used, but still, the center of bed can be queried
        if (auto center_opt = config.opt<ConfigOptionPoint>("center"))
            return { scaled(center_opt->value) };
        
        return {};
    }
    
    return to_points(bed_shape_opt->get_values());
}

Points get_bed_shape(const PrintConfig &cfg)
{
    return to_points(cfg.bed_shape.get_values());
}

Points get_bed_shape(const SLAPrinterConfig &cfg) { return to_points(cfg.bed_shape.get_values()); }

std::string get_sla_suptree_prefix(const DynamicPrintConfig &config)
{
    const auto *suptreetype = config.option<ConfigOptionEnum<sla::SupportTreeType>>("support_tree_type");
    std::string slatree = "";
    if (suptreetype) {
        auto ttype = static_cast<sla::SupportTreeType>(suptreetype->get_int());
        switch (ttype) {
        case sla::SupportTreeType::Branching: slatree = "branching"; break;
        case sla::SupportTreeType::Organic: slatree = "organic"; break;
        default:
            ;
        }
    }

    return slatree;
}

static bool is_XL_printer(const std::string& printer_notes)
{
    return boost::algorithm::contains(printer_notes, "PRINTER_VENDOR_PRUSA3D")
        && boost::algorithm::contains(printer_notes, "PRINTER_MODEL_XL");
}

bool is_XL_printer(const DynamicPrintConfig &cfg)
{
    auto *printer_notes = cfg.opt<ConfigOptionString>("printer_notes");
    return printer_notes && is_XL_printer(printer_notes->value);
}

bool is_XL_printer(const PrintConfig &cfg)
{
    return is_XL_printer(cfg.printer_notes.value);
}

} // namespace Slic3r

#include <cereal/types/polymorphic.hpp>
CEREAL_REGISTER_TYPE(Slic3r::DynamicPrintConfig)
CEREAL_REGISTER_POLYMORPHIC_RELATION(Slic3r::DynamicConfig, Slic3r::DynamicPrintConfig)
