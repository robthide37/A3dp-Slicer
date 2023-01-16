#include "ExtrusionRole.hpp"

#include <string>
#include <string_view>

#define L(s) (s)

namespace Slic3r {

// Convert a rich bitmask based ExtrusionRole to a less expressive ordinal GCodeExtrusionRole.
// GCodeExtrusionRole is to be serialized into G-code and deserialized by G-code viewer,
GCodeExtrusionRole extrusion_role_to_gcode_extrusion_role(ExtrusionRole role)
{
    if (role == ExtrusionRole::None)                return erNone;
    if (role == ExtrusionRole::Perimeter)           return erPerimeter;
    if (role == ExtrusionRole::ExternalPerimeter)   return erExternalPerimeter;
    if (role == ExtrusionRole::OverhangPerimeter)   return erOverhangPerimeter;
    if (role == ExtrusionRole::InternalInfill)      return erInternalInfill;
    if (role == ExtrusionRole::SolidInfill)         return erSolidInfill;
    if (role == ExtrusionRole::TopSolidInfill)      return erTopSolidInfill;
    if (role == ExtrusionRole::Ironing)             return erIroning;
    if (role == ExtrusionRole::BridgeInfill)        return erBridgeInfill;
    if (role == ExtrusionRole::GapFill)             return erGapFill;
    if (role == ExtrusionRole::Skirt)               return erSkirt;
    if (role == ExtrusionRole::SupportMaterial)     return erSupportMaterial;
    if (role == ExtrusionRole::SupportMaterialInterface) return erSupportMaterialInterface;
    if (role == ExtrusionRole::WipeTower)           return erWipeTower;
    assert(false);
    return erNone;
}

std::string gcode_extrusion_role_to_string(GCodeExtrusionRole role)
{
    switch (role) {
        case erNone                         : return L("Unknown");
        case erPerimeter                    : return L("Perimeter");
        case erExternalPerimeter            : return L("External perimeter");
        case erOverhangPerimeter            : return L("Overhang perimeter");
        case erInternalInfill               : return L("Internal infill");
        case erSolidInfill                  : return L("Solid infill");
        case erTopSolidInfill               : return L("Top solid infill");
        case erIroning                      : return L("Ironing");
        case erBridgeInfill                 : return L("Bridge infill");
        case erGapFill                      : return L("Gap fill");
        case erSkirt                        : return L("Skirt/Brim");
        case erSupportMaterial              : return L("Support material");
        case erSupportMaterialInterface     : return L("Support material interface");
        case erWipeTower                    : return L("Wipe tower");
        case erCustom                       : return L("Custom");
        default                             : assert(false);
    }
    return {};
}

GCodeExtrusionRole string_to_gcode_extrusion_role(const std::string_view role)
{
    if (role == L("Perimeter"))
        return erPerimeter;
    else if (role == L("External perimeter"))
        return erExternalPerimeter;
    else if (role == L("Overhang perimeter"))
        return erOverhangPerimeter;
    else if (role == L("Internal infill"))
        return erInternalInfill;
    else if (role == L("Solid infill"))
        return erSolidInfill;
    else if (role == L("Top solid infill"))
        return erTopSolidInfill;
    else if (role == L("Ironing"))
        return erIroning;
    else if (role == L("Bridge infill"))
        return erBridgeInfill;
    else if (role == L("Gap fill"))
        return erGapFill;
    else if (role == L("Skirt") || role == L("Skirt/Brim")) // "Skirt" is for backward compatibility with 2.3.1 and earlier
        return erSkirt;
    else if (role == L("Support material"))
        return erSupportMaterial;
    else if (role == L("Support material interface"))
        return erSupportMaterialInterface;
    else if (role == L("Wipe tower"))
        return erWipeTower;
    else if (role == L("Custom"))
        return erCustom;
    else
        return erNone;
}

}
