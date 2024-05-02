#include "ExtrusionEntityCollection.hpp"
#include "ShortestPath.hpp"
#include <algorithm>
#include <cmath>
#include <map>

namespace Slic3r {

void filter_by_extrusion_role_in_place(ExtrusionEntitiesPtr &extrusion_entities, ExtrusionRole role)
{
	if (role != erMixed) {
		auto first  = extrusion_entities.begin();
		auto last   = extrusion_entities.end();
        extrusion_entities.erase(
            std::remove_if(first, last, [&role](const ExtrusionEntity* ee) {
                return ee->role() != role; }),
            last);
	}
}

ExtrusionEntityCollection::ExtrusionEntityCollection(const ExtrusionPaths &paths)
    : m_no_sort(false), ExtrusionEntity(true)
{
    this->append(paths);
}

ExtrusionEntityCollection& ExtrusionEntityCollection::operator= (const ExtrusionEntityCollection &other)
{
    this->m_no_sort = other.m_no_sort;
    this->m_can_reverse = other.m_can_reverse;
    clear();
    this->append(other.m_entities);
    return *this;
}

void ExtrusionEntityCollection::swap(ExtrusionEntityCollection &c)
{
    std::swap(this->m_entities, c.m_entities);
    std::swap(this->m_no_sort, c.m_no_sort);
    std::swap(this->m_can_reverse, c.m_can_reverse);
}

void ExtrusionEntityCollection::clear()
{
	for (size_t i = 0; i < this->m_entities.size(); ++i)
		delete this->m_entities[i];
    this->m_entities.clear();
}

ExtrusionEntityCollection::operator ExtrusionPaths() const
{
    ExtrusionPaths paths;
    for (const ExtrusionEntity *ptr : this->entities()) {
        if (const ExtrusionPath *path = dynamic_cast<const ExtrusionPath*>(ptr))
            paths.push_back(*path);
    }
    return paths;
}

void ExtrusionEntityCollection::reverse()
{
    for (ExtrusionEntity *ptr : this->m_entities)
    {
        // Don't reverse it if it's a loop, as it doesn't change anything in terms of elements ordering
        // and caller might rely on winding order
        if (ptr->can_reverse() && !ptr->is_loop())
            ptr->reverse();
    }
    std::reverse(this->m_entities.begin(), this->m_entities.end());
}

void ExtrusionEntityCollection::replace(size_t i, const ExtrusionEntity &entity)
{
    delete this->m_entities[i];
    this->m_entities[i] = entity.clone();
}

void ExtrusionEntityCollection::remove(size_t i)
{
    delete this->m_entities[i];
    this->m_entities.erase(this->m_entities.begin() + i);
}

void ExtrusionEntityCollection::chained_path_from(const Point &start_near)
{
    if (this->m_no_sort)
        return;
    chain_and_reorder_extrusion_entities(this->m_entities, &start_near);
}

void ExtrusionEntityCollection::polygons_covered_by_width(Polygons &out, const float scaled_epsilon) const
{
    for (const ExtrusionEntity *entity : this->entities())
        entity->polygons_covered_by_width(out, scaled_epsilon);
}

void ExtrusionEntityCollection::polygons_covered_by_spacing(Polygons &out, const float spacing_ratio, const float scaled_epsilon) const
{
    for (const ExtrusionEntity *entity : this->entities())
        entity->polygons_covered_by_spacing(out, spacing_ratio, scaled_epsilon);
}

// Recursively count paths and loops contained in this collection.
size_t ExtrusionEntityCollection::items_count() const
{
    return CountEntities().count(*this);
}

void
CountEntities::use(const ExtrusionEntityCollection &coll) {
    for (const ExtrusionEntity* entity : coll.entities()) {
        entity->visit(*this);
    }
}

// Returns a single vector of pointers to all non-collection items contained in this one.
ExtrusionEntityCollection ExtrusionEntityCollection::flatten(bool preserve_ordering) const
{
    //ExtrusionEntityCollection coll;
    //this->flatten(&coll, preserve_ordering);
    //return coll;
    return FlatenEntities(preserve_ordering).flatten(*this);

}
void
FlatenEntities::use(const ExtrusionEntityCollection &coll) {
    if ((!coll.can_sort() || !this->to_fill.can_sort()) && preserve_ordering) {
        FlatenEntities unsortable(coll, preserve_ordering);
        for (const ExtrusionEntity* entity : coll.entities()) {
            entity->visit(unsortable);
        }
        to_fill.append(std::move(unsortable.to_fill));
    } else {
        for (const ExtrusionEntity* entity : coll.entities()) {
            entity->visit(*this);
        }
    }
}

ExtrusionEntityCollection&&
FlatenEntities::flatten(const ExtrusionEntityCollection &to_flatten) && {
    use(to_flatten);
    return std::move(to_fill);
}

}
