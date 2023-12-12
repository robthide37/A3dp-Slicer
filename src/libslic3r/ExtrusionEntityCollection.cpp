///|/ Copyright (c) Prusa Research 2016 - 2023 Vojtěch Bubník @bubnikv, Lukáš Hejl @hejllukas, Lukáš Matěna @lukasmatena
///|/ Copyright (c) SuperSlicer 2023 Remi Durand @supermerill
///|/ Copyright (c) Slic3r 2013 - 2016 Alessandro Ranellucci @alranel
///|/ Copyright (c) 2015 Maksim Derbasov @ntfshard
///|/ Copyright (c) 2014 Petr Ledvina @ledvinap
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "ExtrusionEntityCollection.hpp"
#include "ShortestPath.hpp"
#include <algorithm>
#include <cmath>
#include <map>

namespace Slic3r {

#if 0
void filter_by_extrusion_role_in_place(ExtrusionEntitiesPtr &extrusion_entities, ExtrusionRole role)
{
	if (role != ExtrusionRole::Mixed) {
		auto first  = extrusion_entities.begin();
		auto last   = extrusion_entities.end();
        extrusion_entities.erase(
            std::remove_if(first, last, [&role](const ExtrusionEntity* ee) {
                return ee->role() != role; }),
            last);
	}
}
#endif

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

void ExtrusionEntityCollection::flatten(bool preserve_ordering = false, ExtrusionEntityCollection& out) const
{
    if ((!coll.can_sort() || !this->to_fill.can_sort()) && preserve_ordering) {
        out.push_back(this->flatten(preserve_ordering))
    }else{
        FlatenEntities(preserve_ordering) flattener;
        flattener.use(*this);
        //tranfert owner of entities.
        out.entities.insert(out.entities.begin(), flattener.get().entities.begin(), flattener.get().entities.end());
        flattener.set().entities.clear();
    }
}

void FlatenEntities::use(const ExtrusionEntityCollection &coll) {
    if ((!coll.can_sort() || !this->to_fill.can_sort()) && preserve_ordering) {
        ExtrusionEntityCollection copy{pattern.can_sort(), pattern.can_reverse()};
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

ExtrusionEntityCollection&& FlatenEntities::flatten(const ExtrusionEntityCollection &to_flatten) && {
    use(to_flatten);
    return std::move(to_fill);
}

}
