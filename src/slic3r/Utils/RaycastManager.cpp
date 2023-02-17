#include "RaycastManager.hpp"
#include <utility>

// include for earn camera
#include "slic3r/GUI/GUI_App.hpp" 
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/CameraUtils.hpp"

using namespace Slic3r::GUI;

namespace priv {
using namespace Slic3r;
static void actualize(RaycastManager::Meshes &meshes, const ModelVolumePtrs &volumes, const RaycastManager::ISkip *skip);
static const AABBMesh * get_mesh(const RaycastManager::Meshes &meshes, size_t volume_id);
static RaycastManager::TrKey create_key(const ModelVolume* volume, const ModelInstance* instance){ 
    return std::make_pair(instance->id().id, volume->id().id); }
static RaycastManager::TrItems::iterator find(RaycastManager::TrItems &items, const RaycastManager::TrKey &key);
static bool is_lower_key(const RaycastManager::TrKey &k1, const RaycastManager::TrKey &k2) {
    return k1.first < k2.first || k1.first == k2.first && k1.second < k2.second; }
static bool is_lower(const RaycastManager::TrItem &i1, const RaycastManager::TrItem &i2) {
    return is_lower_key(i1.first, i2.first); };
}

void RaycastManager::actualize(const ModelObject *object, const ISkip *skip)
{
    // actualize MeshRaycaster
    priv::actualize(m_meshes, object->volumes, skip);

    // check if inscance was removed
    std::vector<bool> removed_transf(m_transformations.size(), {true});

    bool need_sort = false;
    // actualize transformation matrices
    for (const ModelVolume *volume : object->volumes) {
        if (skip != nullptr && skip->skip(volume->id().id)) continue;
        const Transform3d &volume_tr = volume->get_matrix();
        for (const ModelInstance *instance : object->instances) {
            const Transform3d &instrance_tr = instance->get_matrix();
            Transform3d transformation = instrance_tr * volume_tr;
            TrKey key = priv::create_key(volume, instance);
            auto item = priv::find(m_transformations, key);
            if (item != m_transformations.end()) {
                // actualize transformation all the time
                item->second = transformation;
                size_t index = item - m_transformations.begin();
                removed_transf[index] = false;
            } else {
                // add new transformation
                m_transformations.emplace_back(std::make_pair(key, transformation));
                need_sort = true;
            }
        }
    }

    // clean other transformation
    for (int i = removed_transf.size() - 1; i >= 0; --i)
        if (removed_transf[i])
            m_transformations.erase(m_transformations.begin() + i);

    if (need_sort)
        std::sort(m_transformations.begin(), m_transformations.end(), priv::is_lower);
}

void RaycastManager::actualize(const ModelInstance *instance, const ISkip *skip) {
    const ModelVolumePtrs &volumes = instance->get_object()->volumes;

    // actualize MeshRaycaster
    priv::actualize(m_meshes, volumes, skip);

    // check if inscance was removed
    std::vector<bool> removed_transf(m_transformations.size(), {true});

    bool need_sort = false;
    // actualize transformation matrices
    for (const ModelVolume *volume : volumes) {
        if (skip != nullptr && skip->skip(volume->id().id))
            continue;
        const Transform3d &volume_tr = volume->get_matrix();
        const Transform3d &instrance_tr   = instance->get_matrix();
        Transform3d        transformation = instrance_tr * volume_tr;
        TrKey key = priv::create_key(volume, instance);
        auto item = priv::find(m_transformations, key);
        if (item != m_transformations.end()) {
            // actualize transformation all the time
            item->second          = transformation;
            size_t index          = item - m_transformations.begin();
            removed_transf[index] = false;
        } else {
            // add new transformation
            m_transformations.emplace_back(std::make_pair(key, transformation));
            need_sort = true;
        }        
    }

    // clean other transformation
    for (int i = removed_transf.size() - 1; i >= 0; --i)
        if (removed_transf[i])
            m_transformations.erase(m_transformations.begin() + i);

    if (need_sort)
        std::sort(m_transformations.begin(), m_transformations.end(), priv::is_lower);
}
 
std::optional<RaycastManager::Hit> RaycastManager::ray_from_camera(
    const Vec2d &mouse_pos, const Camera &camera, const ISkip *skip) const
{
    // Improve it is not neccessaru to use AABBMesh and calc normal in 

    struct Result
    {
        const AABBMesh *mesh = nullptr;
        double squared_distance;
        int face;
        Vec3d hit_world;
        const Transform3d *tramsformation;
        const TrKey *key;
    }result;
    for (const auto &item : m_transformations) { 
        const TrKey &key = item.first;
        size_t volume_id = key.second;
        if (skip != nullptr && skip->skip(volume_id)) continue;
        const AABBMesh *mesh = priv::get_mesh(m_meshes, volume_id);
        if (mesh == nullptr) continue;
        const Transform3d &transformation = item.second;

        Vec3d point;
        Vec3d direction;
        CameraUtils::ray_from_screen_pos(camera, mouse_pos, point, direction);
        Transform3d inv = transformation.inverse();
        point           = inv * point;
        direction       = inv.linear() * direction;
        std::vector<AABBMesh::hit_result> hits = mesh->query_ray_hits(point, direction);
        if (hits.empty()) continue; // no intersection found

        const AABBMesh::hit_result &hit = hits.front();

        // convert to world
        Vec3d hit_world = transformation * hit.position();
        double squared_distance = (camera.get_position() - hit_world).squaredNorm();
        if (result.mesh != nullptr &&
            result.squared_distance < squared_distance)
            continue;

        result.mesh = mesh;
        result.squared_distance = squared_distance;
        result.face = hit.face();
        result.hit_world = hit_world;
        result.tramsformation = &transformation;
        result.key = &key;
    }

    if (result.mesh == nullptr)
        return {};

    const Vec3i tri = result.mesh->indices(result.face);
    Vec3d pts[3];
    auto tr = result.tramsformation->linear();
    for (int i = 0; i < 3; ++i)
        pts[i] = tr * result.mesh->vertices(tri[i]).cast<double>();
    Vec3d normal_world = (pts[1] - pts[0]).cross(pts[2] - pts[1]);
    normal_world.normalize();

    SurfacePoint<double> point_world{result.hit_world, normal_world};
    return RaycastManager::Hit{point_world, *result.key, result.squared_distance};
}

std::optional<RaycastManager::Hit> RaycastManager::unproject(const Vec3d &point, const Vec3d &direction, const ISkip *skip) const
{
    std::optional<Hit> closest;
    for (const auto &item : m_transformations) { 
        const TrKey &key = item.first;
        size_t volume_id = key.second;
        if (skip != nullptr && skip->skip(volume_id)) continue;
        const Transform3d &transformation = item.second;
        const AABBMesh *mesh = priv::get_mesh(m_meshes, volume_id);
        if (mesh == nullptr) continue;
        Transform3d tr_inv = transformation.inverse();
        Vec3d mesh_point = tr_inv * point;
        Vec3d mesh_direction = tr_inv.linear() * direction;

        // Need for detect that actual point position is on correct place
        Vec3d point_positive = mesh_point - mesh_direction;
        Vec3d point_negative = mesh_point + mesh_direction;

        // Throw ray to both directions of ray
        std::vector<AABBMesh::hit_result> hits = mesh->query_ray_hits(point_positive, mesh_direction);
        std::vector<AABBMesh::hit_result> hits_neg = mesh->query_ray_hits(point_negative, -mesh_direction);
        hits.insert(hits.end(), std::make_move_iterator(hits_neg.begin()), std::make_move_iterator(hits_neg.end()));
        for (const AABBMesh::hit_result &hit : hits) { 
            double squared_distance = (mesh_point - hit.position()).squaredNorm();
            if (closest.has_value() &&
                closest->squared_distance < squared_distance)
                continue;
            closest = Hit{{hit.position(), hit.normal()}, key, squared_distance};
        }
    }
    return closest;
}

std::optional<RaycastManager::ClosePoint> RaycastManager::closest(const Vec3d &point, const ISkip *skip) const
{
    std::optional<ClosePoint> closest;
    for (const auto &item : m_transformations) {
        const TrKey &key       = item.first;
        size_t       volume_id = key.second;
        if (skip != nullptr && skip->skip(volume_id))
            continue;
        const AABBMesh *mesh = priv::get_mesh(m_meshes, volume_id);
        if (mesh == nullptr) continue;
        const Transform3d &transformation = item.second;
        Transform3d tr_inv = transformation.inverse();
        Vec3d mesh_point = tr_inv * point;
                
        int   face_idx = 0;
        Vec3d closest_point;
        Vec3d pointd = point.cast<double>();
        mesh->squared_distance(pointd, face_idx, closest_point);

        double squared_distance = (mesh_point - closest_point).squaredNorm();
        if (closest.has_value() && closest->squared_distance < squared_distance)
            continue;

        closest = ClosePoint{key, closest_point, squared_distance};
    }
    return closest;
}

Slic3r::Transform3d RaycastManager::get_transformation(const TrKey &tr_key) const {
    // TODO: transformations are sorted use lower bound
    auto item = std::find_if(m_transformations.begin(),
                             m_transformations.end(),
                             [&tr_key](const TrItem &it) -> bool {
                                 return it.first == tr_key;
                             });
    if (item == m_transformations.end()) return Transform3d::Identity();
    return item->second;
}

void priv::actualize(RaycastManager::Meshes &meshes, const ModelVolumePtrs &volumes, const RaycastManager::ISkip *skip)
{
    // check if volume was removed
    std::vector<bool> removed_meshes(meshes.size(), {true});
    bool need_sort = false;
    // actualize MeshRaycaster
    for (const ModelVolume *volume : volumes) {
        size_t oid = volume->id().id;
        if (skip != nullptr && skip->skip(oid))
            continue;
        auto item = std::find_if(meshes.begin(), meshes.end(), [oid](const RaycastManager::Mesh &it) -> bool { return oid == it.first; });
        if (item == meshes.end()) {
            // add new raycaster
            bool calculate_epsilon = true;
            auto mesh = std::make_unique<AABBMesh>(volume->mesh(), calculate_epsilon);
            meshes.emplace_back(std::make_pair(oid, std::move(mesh)));
            need_sort = true;
        } else {
            size_t index = item - meshes.begin();
            removed_meshes[index] = false;
        }
    }

    // clean other raycasters
    for (int i = removed_meshes.size() - 1; i >= 0; --i)
        if (removed_meshes[i])
            meshes.erase(meshes.begin() + i);

    // All the time meshes must be sorted by volume id - for faster search
    if (need_sort) {
        auto is_lower = [](const RaycastManager::Mesh &m1, const RaycastManager::Mesh &m2) { return m1.first < m2.first; };
        std::sort(meshes.begin(), meshes.end(), is_lower);
    }
}

const Slic3r::AABBMesh *priv::get_mesh(const RaycastManager::Meshes &meshes, size_t volume_id)
{
    auto is_lower_index = [](const RaycastManager::Mesh &m, size_t i) -> bool { return m.first < i; };
    auto it = std::lower_bound(meshes.begin(), meshes.end(), volume_id, is_lower_index);
    if (it == meshes.end() || it->first != volume_id)
        return nullptr;
    return &(*(it->second));
}

RaycastManager::TrItems::iterator priv::find(RaycastManager::TrItems &items, const RaycastManager::TrKey &key) {
    auto fnc = [](const RaycastManager::TrItem &it, const RaycastManager::TrKey &key)->bool { 
        return priv::is_lower_key(it.first, key);
    };
    auto it = std::lower_bound(items.begin(), items.end(), key, fnc);
    if (it == items.end() || it->first != key)
        return items.end();
    return it;
}
