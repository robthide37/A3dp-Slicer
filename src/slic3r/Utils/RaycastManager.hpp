#ifndef slic3r_RaycastManager_hpp_
#define slic3r_RaycastManager_hpp_

#include <memory> // unique_ptr
#include <optional> // unique_ptr
#include <map>
#include "slic3r/GUI/MeshUtils.hpp" // MeshRaycaster
#include "libslic3r/Point.hpp" // Transform3d
#include "libslic3r/ObjectID.hpp"
#include "libslic3r/Model.hpp" // ModelObjectPtrs, ModelObject, ModelInstance, ModelVolume

namespace Slic3r::GUI{

/// <summary>
/// Cast rays from camera to scene
/// Used for find hit point on model volume under mouse cursor
/// </summary>
class RaycastManager
{
    //               ModelVolume
    std::map<size_t, std::unique_ptr<MeshRaycaster>> m_raycasters;

    // Key for transformation consist of unique volume and instance
    //                 ModelInstance, ModelVolume
    using TrKey = std::pair<size_t, size_t>;
    std::map<TrKey, Transform3d> m_transformations;

    // should contain shared pointer to camera but it is not shared pointer so it need it every time when casts rays

public:
    class ISkip{
    public:
        virtual ~ISkip() = default;
        /// <summary>
        /// Condition to not process specific transformation
        /// </summary>
        /// <param name="key">Transformation key</param>
        /// <returns>True on skip otherwise false</returns>
        //virtual bool skip(const TrKey &key) const { return false; }

        /// <summary>
        /// Condition to not process model volume
        /// </summary>
        /// <param name="model_volume_id">ObjectID of model volume to not process</param>
        /// <returns>True on skip otherwise false</returns>
        virtual bool skip(const size_t &model_volume_id) const { return false; }
    };

    /// <summary>
    /// Actualize raycasters + transformation
    /// Detection of removed object
    /// Detection of removed instance
    /// Detection of removed volume
    /// </summary>
    /// <param name="skip">Condifiton for skip actualization</param>
    /// <param name="objects">Model representation</param>
    void actualize(const ModelObjectPtrs &objects,
                   const ISkip *          skip = nullptr);

    void actualize(const ModelObject *object, const ISkip *skip = nullptr);

    // TODO: it is more general object move outside of this class
    struct SurfacePoint
    {
        Vec3f position = Vec3f::Zero();
        Vec3f normal   = Vec3f::UnitZ();
        SurfacePoint() = default;
        SurfacePoint(Vec3f position, Vec3f normal)
            : position(position), normal(normal)
        {}
    };

    struct Hit: public SurfacePoint
    {
        TrKey tr_key;
        Hit(TrKey tr_key, Vec3f position, Vec3f normal)
            : SurfacePoint(position, normal), tr_key(tr_key)
        {} 
    };    

    class SkipVolume: public ISkip
    {
        size_t volume_id;
    public:
        SkipVolume(size_t volume_id) : volume_id(volume_id) {}
        bool skip(const size_t &model_volume_id) const override { return model_volume_id == volume_id; }
    };

    class AllowVolumes: public ISkip
    {
        std::vector<size_t> allowed_id;
    public:
        AllowVolumes(std::vector<size_t> allowed_id) : allowed_id(allowed_id) {}
        bool skip(const size_t &model_volume_id) const override {
            auto it = std::find(allowed_id.begin(), allowed_id.end(), model_volume_id);
            return it == allowed_id.end(); 
        }
    };

    /// <summary>
    /// Unproject on mesh by Mesh raycasters
    /// Note: Function use current camera position from wxGetApp()
    /// </summary>
    /// <param name="mouse_pos">Position of mouse on screen</param>
    /// <param name="camera">Projection params</param>
    /// <param name="skip">Define which caster will be skipped, null mean no skip</param>
    /// <returns>Position on surface, normal direction and transformation key, which define hitted object instance</returns>
    std::optional<Hit> unproject(const Vec2d &mouse_pos, 
                                 const Camera &camera,
                                 const ISkip *skip = nullptr) const;

    /// <summary>
    /// Getter on transformation
    /// </summary>
    /// <param name="tr_key">Define transformation</param>
    /// <returns>Transformation for key</returns>
    Transform3d get_transformation(const TrKey &tr_key) const;
};

} // namespace Slic3r::GUI

#endif // slic3r_RaycastManager_hpp_
