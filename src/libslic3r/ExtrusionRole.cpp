///|/ Copyright (c) Prusa Research 2023 Pavel Mikuš @Godrak, Oleksandra Iushchenko @YuSanka, Vojtěch Bubník @bubnikv
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "ExtrusionRole.hpp"
#include "I18N.hpp"

#include <string>
#include <string_view>
#include <cassert>


namespace Slic3r {

// Convert a rich bitmask based ExtrusionRole to a less expressive ordinal GCodeExtrusionRole.
// GCodeExtrusionRole is to be serialized into G-code and deserialized by G-code viewer,
GCodeExtrusionRole extrusion_role_to_gcode_extrusion_role(ExtrusionRole role)
{
    assert(role != ExtrusionRole::Mixed);
    if (role == ExtrusionRole::None)                return GCodeExtrusionRole::None;
    if (role.is_perimeter()) {
        return role.is_bridge() ? GCodeExtrusionRole::OverhangPerimeter :
               role.is_external() ? GCodeExtrusionRole::ExternalPerimeter : GCodeExtrusionRole::Perimeter;
    }
    if (role == ExtrusionRole::InternalInfill)      return GCodeExtrusionRole::InternalInfill;
    if (role == ExtrusionRole::SolidInfill)         return GCodeExtrusionRole::SolidInfill;
    if (role == ExtrusionRole::TopSolidInfill)      return GCodeExtrusionRole::TopSolidInfill;
    if (role == ExtrusionRole::Ironing)             return GCodeExtrusionRole::Ironing;
    if (role == ExtrusionRole::BridgeInfill)        return GCodeExtrusionRole::BridgeInfill;
    if (role == ExtrusionRole::InternalBridgeInfill)return GCodeExtrusionRole::InternalBridgeInfill;
    if (role == ExtrusionRole::ThinWall)            return GCodeExtrusionRole::ThinWall;
    if (role == ExtrusionRole::GapFill)             return GCodeExtrusionRole::GapFill;
    if (role == ExtrusionRole::Skirt)               return GCodeExtrusionRole::Skirt;
    if (role == ExtrusionRole::SupportMaterial)     return GCodeExtrusionRole::SupportMaterial;
    if (role == ExtrusionRole::SupportMaterialInterface) return GCodeExtrusionRole::SupportMaterialInterface;
    if (role == ExtrusionRole::WipeTower)           return GCodeExtrusionRole::WipeTower;
    if (role == ExtrusionRole::Milling)             return GCodeExtrusionRole::Milling;
    assert(false);
    return GCodeExtrusionRole::None;
}

std::string gcode_extrusion_role_to_string(GCodeExtrusionRole role)
{
    switch (role) {
        case GCodeExtrusionRole::None                         : return L("Unknown");
        case GCodeExtrusionRole::Perimeter                    : return L("Perimeter");
        case GCodeExtrusionRole::ExternalPerimeter            : return L("External perimeter");
        case GCodeExtrusionRole::OverhangPerimeter            : return L("Overhang perimeter");
        case GCodeExtrusionRole::InternalInfill               : return L("Internal infill");
        case GCodeExtrusionRole::SolidInfill                  : return L("Solid infill");
        case GCodeExtrusionRole::TopSolidInfill               : return L("Top solid infill");
        case GCodeExtrusionRole::Ironing                      : return L("Ironing");
        case GCodeExtrusionRole::BridgeInfill                 : return L("Bridge infill");
        case GCodeExtrusionRole::InternalBridgeInfill         : return L("Internal bridge infill");
        case GCodeExtrusionRole::ThinWall                     : return L("Thin wall");
        case GCodeExtrusionRole::GapFill                      : return L("Gap fill");
        case GCodeExtrusionRole::Skirt                        : return L("Skirt/Brim");
        case GCodeExtrusionRole::SupportMaterial              : return L("Support material");
        case GCodeExtrusionRole::SupportMaterialInterface     : return L("Support material interface");
        case GCodeExtrusionRole::WipeTower                    : return L("Wipe tower");
        case GCodeExtrusionRole::Milling                      : return L("Milling");
        case GCodeExtrusionRole::Custom                       : return L("Custom");
        case GCodeExtrusionRole::Travel                       : return L("Travel");
        default                             : assert(false);
    }
    return {};
}

GCodeExtrusionRole string_to_gcode_extrusion_role(const std::string_view role)
{
    if (role == L("Perimeter"))
        return GCodeExtrusionRole::Perimeter;
    else if (role == L("External perimeter"))
        return GCodeExtrusionRole::ExternalPerimeter;
    else if (role == L("Overhang perimeter"))
        return GCodeExtrusionRole::OverhangPerimeter;
    else if (role == L("Internal infill"))
        return GCodeExtrusionRole::InternalInfill;
    else if (role == L("Solid infill"))
        return GCodeExtrusionRole::SolidInfill;
    else if (role == L("Top solid infill"))
        return GCodeExtrusionRole::TopSolidInfill;
    else if (role == L("Ironing"))
        return GCodeExtrusionRole::Ironing;
    else if (role == L("Bridge infill"))
        return GCodeExtrusionRole::BridgeInfill;
    else if (role == L("Internal bridge infill"))
        return GCodeExtrusionRole::InternalBridgeInfill;
    else if (role == L("Thin wall"))
        return GCodeExtrusionRole::ThinWall;
    else if (role == L("Gap fill"))
        return GCodeExtrusionRole::GapFill;
    else if (role == L("Skirt") || role == L("Skirt/Brim")) // "Skirt" is for backward compatibility with 2.3.1 and earlier
        return GCodeExtrusionRole::Skirt;
    else if (role == L("Support material"))
        return GCodeExtrusionRole::SupportMaterial;
    else if (role == L("Support material interface"))
        return GCodeExtrusionRole::SupportMaterialInterface;
    else if (role == L("Wipe tower"))
        return GCodeExtrusionRole::WipeTower;
    else if (role == L("Milling"))
        return GCodeExtrusionRole::Milling;
    else if (role == L("Custom"))
        return GCodeExtrusionRole::Custom;
    else if (role == L("Trabel"))
        return GCodeExtrusionRole::Travel;
    else
        assert(false);
    return GCodeExtrusionRole::None;
}

std::string role_to_code(ExtrusionRole role)
{
    if (ExtrusionRole::None == role)
        return L("None");
    else if (ExtrusionRole::Perimeter == role)
        return L("IPeri");
    else if (ExtrusionRole::ExternalPerimeter == role)
        return L("EPeri");
    else if (ExtrusionRole::OverhangPerimeter == role)
        return L("OPeri");
    else if (ExtrusionRole::InternalInfill == role)
        return L("IFill");
    else if (ExtrusionRole::SolidInfill == role)
        return L("SFill");
    else if (ExtrusionRole::TopSolidInfill == role)
        return L("TFill");
    else if (ExtrusionRole::Ironing == role)
        return L("Iron");
    else if (ExtrusionRole::BridgeInfill == role)
        return L("EBridge");
    else if (ExtrusionRole::InternalBridgeInfill == role)
        return L("IBridge");
    else if (ExtrusionRole::ThinWall == role)
        return L("ThinW");
    else if (ExtrusionRole::GapFill == role)
        return L("GFill");
    else if (ExtrusionRole::Skirt == role)
        return L("Skirt");
    else if (ExtrusionRole::SupportMaterial == role)
        return L("Supp");
    else if (ExtrusionRole::SupportMaterialInterface == role)
        return L("SuppI");
    else if (ExtrusionRole::WipeTower == role)
        return L("WTower");
    else if (ExtrusionRole::Milling == role)
        return L("Mill");
    //else if (ExtrusionRole::Custom == role)
    //    return L("Custom");
    else if (ExtrusionRole::Mixed == role)
        return L("Mixed");
    else if (ExtrusionRole::Travel == role)
        return L("Travel");
    else
        assert(false);
    return "";
}

std::string looprole_to_code(ExtrusionLoopRole looprole)
{
    std::string code;
    if(elrDefault == (looprole & elrDefault))
        code += std::string("D");
    if(elrInternal == (looprole & elrInternal))
        code += std::string("Int");
    if(elrSkirt == (looprole & elrSkirt))
        code += std::string("Skirt");
    if(elrHole == (looprole & elrHole))
        code += std::string("Hole");
    if(elrVase == (looprole & elrVase))
        code += std::string("Vase");
    if(elrFirstLoop == (looprole & elrFirstLoop))
        code += std::string("First");

    return code;
}

}
