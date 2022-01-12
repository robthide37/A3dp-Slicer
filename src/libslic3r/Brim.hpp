#ifndef slic3r_Brim_hpp_
#define slic3r_Brim_hpp_

namespace Slic3r {

class Print;

#if 0
// Produce brim lines around those objects, that have the brim enabled.
// Collect islands_area to be merged into the final 1st layer convex hull.
ExtrusionEntityCollection make_brim(const Print &print, PrintTryCancel try_cancel, Polygons &islands_area);
#endif
void make_brim(const Print& print, const Flow& flow, const PrintObjectPtrs& objects, ExPolygons& unbrimmable, ExtrusionEntityCollection& out);
void make_brim_ears(const Print& print, const Flow& flow, const PrintObjectPtrs& objects, ExPolygons& unbrimmable, ExtrusionEntityCollection& out);
void make_brim_interior(const Print& print, const Flow& flow, const PrintObjectPtrs& objects, ExPolygons& unbrimmable_areas, ExtrusionEntityCollection& out);

} // Slic3r

#endif // slic3r_Brim_hpp_
