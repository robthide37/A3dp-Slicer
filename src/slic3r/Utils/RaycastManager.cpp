#include "RaycastManager.hpp"
#include <utility>

// include for earn camera
#include "slic3r/GUI/GUI_App.hpp" 
#include "slic3r/GUI/Plater.hpp"  
#include "slic3r/GUI/Camera.hpp"

using namespace Slic3r::GUI;   

void RaycastManager::actualize(const ModelObject *object, const ISkip *skip)
{
    // check if volume was removed
    std::vector<bool> removed_casters(m_raycasters.size(), {true});
    // check if inscance was removed
    std::vector<bool> removed_transf(m_transformations.size(), {true});
    
    // actualize MeshRaycaster
    for (const ModelVolume *volume : object->volumes) {
        size_t oid = volume->id().id;
        if (skip != nullptr && skip->skip(oid)) 
            continue;
        auto item = std::find_if(m_raycasters.begin(), m_raycasters.end(),
                              [oid](const RaycastManager::Raycaster &it)->bool {
                                  return oid == it.first;
                              });
        if (item == m_raycasters.end()) {
            // add new raycaster
#if ENABLE_RAYCAST_PICKING
            auto raycaster = std::make_unique<MeshRaycaster>(volume->get_mesh_shared_ptr());
#else // !ENABLE_RAYCAST_PICKING
            auto raycaster = std::make_unique<MeshRaycaster>(volume->mesh());
#endif // ENABLE_RAYCAST_PICKING
            m_raycasters.emplace_back(std::make_pair(oid, std::move(raycaster)));
        } else {
            size_t index = item - m_raycasters.begin();
            removed_casters[index] = false;
        }
    }

    // actualize transformation matrices
    for (const ModelVolume *volume : object->volumes) {
        if (skip != nullptr && skip->skip(volume->id().id)) continue;
        const Transform3d &volume_tr = volume->get_matrix();
        for (const ModelInstance *instance : object->instances) {
            const Transform3d &instrance_tr   = instance->get_matrix();
            Transform3d        transformation = instrance_tr * volume_tr;
            // TODO: add SLA shift Z
            // transformation.translation()(2) += m_sla_shift_z;
            TrKey tr_key = std::make_pair(instance->id().id, volume->id().id);
            auto  item   = std::find_if(m_transformations.begin(),
                                        m_transformations.end(),
                                        [&tr_key](const TrItem &it) -> bool {
                                         return it.first == tr_key;
                                        });
            if (item != m_transformations.end()) {
                // actualize transformation all the time
                item->second = transformation;
                size_t index = item - m_transformations.begin();
                removed_transf[index] = false;
            } else {
                // add new transformation
                m_transformations.emplace_back(
                    std::make_pair(tr_key, transformation));
            }
        }
    }

    // clean other raycasters
    for (int i = removed_casters.size() - 1; i >= 0; --i) 
        if (removed_casters[i]) 
            m_raycasters.erase(m_raycasters.begin() + i);
    
    // clean other transformation
    for (int i = removed_transf.size() - 1; i >= 0; --i)
        if (removed_transf[i])
            m_transformations.erase(m_transformations.begin() + i);
}

std::optional<RaycastManager::Hit> RaycastManager::unproject(
    const Vec2d &mouse_pos, const Camera &camera, const ISkip *skip) const
{
    struct HitWithDistance: public Hit
    {
        double squared_distance;
        HitWithDistance(double              squared_distance,
                        const TrKey &       key,
                        const SurfacePoint &surface_point)
            : Hit(key, surface_point.position, surface_point.normal)
            , squared_distance(squared_distance)
        {}
    };
    std::optional<HitWithDistance> closest;
    for (const auto &item : m_transformations) { 
        const TrKey &key = item.first;
        size_t       volume_id = key.second;
        if (skip != nullptr && skip->skip(volume_id)) continue;
        const Transform3d &transformation = item.second;
        auto raycaster_it =
            std::find_if(m_raycasters.begin(), m_raycasters.end(),
                         [volume_id](const RaycastManager::Raycaster &it)
                             -> bool { return volume_id == it.first; });
        if (raycaster_it == m_raycasters.end()) continue;
        const MeshRaycaster &raycaster = *(raycaster_it->second);
        SurfacePoint surface_point;
        bool success = raycaster.unproject_on_mesh(
            mouse_pos, transformation, camera,
            surface_point.position, surface_point.normal);
        if (!success) continue;

        Vec3d act_hit_tr = transformation * surface_point.position.cast<double>();
        double squared_distance = (camera.get_position() - act_hit_tr).squaredNorm();
        if (closest.has_value() &&
            closest->squared_distance < squared_distance)
            continue;
        closest = HitWithDistance(squared_distance, key, surface_point);
    }

    //if (!closest.has_value()) return {};
    return closest;
}

Slic3r::Transform3d RaycastManager::get_transformation(const TrKey &tr_key) const {
    auto item = std::find_if(m_transformations.begin(),
                             m_transformations.end(),
                             [&tr_key](const TrItem &it) -> bool {
                                 return it.first == tr_key;
                             });
    if (item == m_transformations.end()) return Transform3d::Identity();
    return item->second;
}