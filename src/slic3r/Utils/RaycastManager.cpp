#include "RaycastManager.hpp"

// include for earn camera
#include "slic3r/GUI/GUI_App.hpp" 
#include "slic3r/GUI/Plater.hpp"  
#include "slic3r/GUI/Camera.hpp"

using namespace Slic3r::GUI;   

void RaycastManager::actualize(const ModelObjectPtrs &objects,
                               const ISkip *          skip)
{
    // check if volume was removed
    std::set<size_t> removed_casters;
    for (const auto &raycaster_item : m_raycasters)
        removed_casters.insert(raycaster_item.first);
    
    // check if inscance was removed
    std::set<TrKey> removed_transformation;
    for (const auto &item : m_transformations)
        removed_transformation.insert(item.first);

    for (const ModelObject *object : objects) {
        // actualize MeshRaycaster
        for (const ModelVolume *volume : object->volumes) {
            size_t oid = volume->id().id;
            if (skip != nullptr && skip->skip(oid)) { 
                removed_casters.erase(oid);
                continue;
            }
            auto item = m_raycasters.find(oid);
            if (item != m_raycasters.end()) {
                removed_casters.erase(oid);
                // alredy in list only actualize
                // TODO: check triangles when change than actualize MeshRaycaster
            } else {
                // add new raycaster
                auto raycaster = std::make_unique<MeshRaycaster>(
                    volume->mesh());
                m_raycasters.insert(std::make_pair(oid, std::move(raycaster)));
            }
        }

        // actualize transformation matrices
        for (const ModelVolume *volume : object->volumes) {
            if (skip != nullptr && skip->skip(volume->id().id)) continue;
            const Transform3d &volume_tr = volume->get_matrix();
            for (const ModelInstance *instance : object->instances) {
                const Transform3d& instrance_tr = instance->get_matrix();
                Transform3d        transformation = instrance_tr * volume_tr;
                // TODO: add SLA shift Z
                // transformation.translation()(2) += m_sla_shift_z;

                TrKey tr_key = std::make_pair(instance->id().id, volume->id().id);
                auto  item   = m_transformations.find(tr_key);
                if (item != m_transformations.end()) {
                    // actualize transformation all the time
                    item->second = transformation;
                    removed_transformation.erase(tr_key);                
                } else {
                    // add new transformation
                    m_transformations.insert(
                        std::make_pair(tr_key, transformation));
                }
            }
        }
    }

    // remove non existing volumes
    for (size_t volume_oid : removed_casters) m_raycasters.erase(volume_oid);
    // remove non existing transformations
    for (const TrKey& transformation_key : removed_transformation)
        m_transformations.erase(transformation_key);
}

void RaycastManager::actualize(const ModelObject *object, const ISkip *skip)
{
    // actualize MeshRaycaster
    for (const ModelVolume *volume : object->volumes) {
        size_t oid = volume->id().id;
        if (skip != nullptr && skip->skip(oid)) 
            continue;
        auto item = m_raycasters.find(oid);
        if (item == m_raycasters.end()) {
            // add new raycaster
            auto raycaster = std::make_unique<MeshRaycaster>(volume->mesh());
            m_raycasters.insert(std::make_pair(oid, std::move(raycaster)));
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
            auto  item   = m_transformations.find(tr_key);
            if (item != m_transformations.end()) {
                // actualize transformation all the time
                item->second = transformation;
            } else {
                // add new transformation
                m_transformations.insert(
                    std::make_pair(tr_key, transformation));
            }
        }
    }
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
        auto raycaster_it = m_raycasters.find(volume_id);
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
    auto item = m_transformations.find(tr_key);
    if (item == m_transformations.end()) return Transform3d::Identity();
    return item->second;
}