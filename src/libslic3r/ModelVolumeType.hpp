#ifndef slic3r_ModelVolumeType_hpp_
#define slic3r_ModelVolumeType_hpp_

namespace Slic3r {

enum class ModelVolumeType : int {
    INVALID    = -1,
    MODEL_PART = 0,
    NEGATIVE_VOLUME,
    PARAMETER_MODIFIER,
    SUPPORT_BLOCKER,
    SUPPORT_ENFORCER,
};

} // namespace Slic3r
#endif /* slic3r_ModelVolumeType_hpp_ */
