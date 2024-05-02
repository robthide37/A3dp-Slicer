#ifndef slic3r_ExtrusionEntity_hpp_
#define slic3r_ExtrusionEntity_hpp_

#include "libslic3r.h"
#include "Polygon.hpp"
#include "Polyline.hpp"

#include <assert.h>
#include <string_view>
#include <numeric>

namespace Slic3r {

class ExPolygonCollection;
class ExtrusionEntityCollection;
class Extruder;

// Each ExtrusionRole value identifies a distinct set of { extruder, speed }
/*
enum ExtrusionRoleModifier : uint16_t {
    ermPerimeter = (1 << 0),
    ermInfill = (2 << 1),
    ermThin = (2 << 2),
    ermSkirt = (2 << 3),
    ermOther = (2 << 4),
    ermInternal = (1 << 10),
    ermExternal = (1 << 11),
    ermSolid = (1 << 12),
    ermBridge = (1 << 13),
    ermSupport = (1 << 13)
};
enum ExtrusionRole : uint16_t {
    erNone = 0,
    erPerimeter = ermPerimeter | ermInternal,
    erExternalPerimeter = ermPerimeter | ermExternal,
    erOverhangPerimeter = ermPerimeter | ermBridge,
    erInternalInfill = ermInfill | ermInternal,
    erSolidInfill = ermInfill | ermSolid | ermInternal,
    erTopSolidInfill = ermInfill | ermSolid | ermExternal,
    erBridgeInfill = ermInfill | ermSolid | ermBridge | ermExternal,
    erInternalBridgeInfill = ermInfill | ermSolid | ermBridge,
    erThinWall = ermThin | ermInternal,
    erGapFill = ermThin | ermExternal,
    erSkirt = ermSkirt,
    erSupportMaterial = ermInfill | ermSupport | ermInternal,
    erSupportMaterialInterface = ermInfill | ermSupport | ermExternal,
    erWipeTower = ermSkirt | ermSupport,
    erMilling = ermOther | ermPerimeter,
    erLaser = ermOther | erSupportMaterialInterface,
    erCustom = ermOther | ermSkirt,
    // Extrusion role for a collection with multiple extrusion roles.
    erMixed = ermOther
};
*/
enum ExtrusionRole : uint8_t {
    erNone, // erNone needs to be 0 (for cooling and everything that use the role as index). It must not be used. Please use erTravel for not-extrusion role.
    erPerimeter,
    erExternalPerimeter,
    erOverhangPerimeter,
    erInternalInfill,
    erSolidInfill,
    erTopSolidInfill,
    erIroning,
    erBridgeInfill,
    erInternalBridgeInfill,
    erThinWall,
    erGapFill,
    erSkirt,
    erSupportMaterial,
    erSupportMaterialInterface,
    erWipeTower,
    erMilling,
    erLaser,
    erCustom,
    // Extrusion role for a collection with multiple extrusion roles.
    erMixed,
    //extrusion role when there is no extrusion
    erTravel,
    erCount
};

// perimeter / infill / support / skirt / gapfill / wipetower / custom / mixed
// side / internal / top / bottom
// bridge

// Special flags describing loop
enum ExtrusionLoopRole : uint16_t {
    // useless
    elrDefault = 1 << 0, //1
    // doesn't contains more contour: it's the most internal one
    elrInternal = 1 << 1, //2
    elrSkirt =    1 << 2,  //4
    //it's a modifier that indicate that the loop is around a hole, not around the infill
    elrHole = 1 << 3, // 8
    //it's a modifier that indicate that the loop should be printed as vase
    elrVase = 1 << 4, //16
    //it's a modifier that indicate that the loop does not contains an inner loop, used for random seam
    elrFirstLoop = 1 << 5, //32
};


inline bool is_perimeter(ExtrusionRole role)
{
    return role == erPerimeter
        || role == erExternalPerimeter
        || role == erThinWall
        || role == erOverhangPerimeter;
}

inline bool is_infill(ExtrusionRole role)
{
    return role == erBridgeInfill
        || role == erInternalBridgeInfill
        || role == erInternalInfill
        || role == erSolidInfill
        || role == erTopSolidInfill
        || role == erIroning;
}

inline bool is_solid_infill(ExtrusionRole role)
{
    return role == erBridgeInfill
        || role == erInternalBridgeInfill
        || role == erSolidInfill
        || role == erTopSolidInfill
        || role == erIroning;
}

inline bool is_bridge(ExtrusionRole role) {
    return role == erBridgeInfill
        || role == erInternalBridgeInfill
        || role == erOverhangPerimeter;
}


class ExtrusionEntity;
class ExtrusionPath;
class ExtrusionPath3D;
class ExtrusionMultiPath;
class ExtrusionMultiPath3D;
class ExtrusionLoop;
//
//class ExtrusionVisitor {
//public:
//    virtual void default_use(ExtrusionEntity &entity) { assert(false); };
//    virtual void use(ExtrusionPath &path) { ExtrusionEntity &entity = path; default_use(entity); };
//    virtual void use(ExtrusionPath3D &path3D) { ExtrusionPath &path = path3D;  use(path); };
//    virtual void use(ExtrusionMultiPath &multipath) { ExtrusionEntity &entity = multipath; default_use(entity); };
//    virtual void use(ExtrusionMultiPath3D &multipath3D) { ExtrusionEntity &entity = multipath3D; default_use(entity); };
//    virtual void use(ExtrusionLoop &loop) { ExtrusionEntity &entity = loop; default_use(entity); };
//    virtual void use(ExtrusionEntityCollection &collection) { ExtrusionEntity &entity = collection;  default_use(entity); };
//};
//class ExtrusionVisitorConst {
//public:
//    virtual void default_use(const ExtrusionEntity &entity) { assert(false); };
//    virtual void use(const ExtrusionPath &path) { const ExtrusionEntity &entity = path; default_use(entity); };
//    virtual void use(const ExtrusionPath3D &path3D) { const ExtrusionPath &path = path3D;  use(path); };
//    virtual void use(const ExtrusionMultiPath &multipath) { const ExtrusionEntity &entity = multipath; default_use(entity); };
//    virtual void use(const ExtrusionMultiPath3D &multipath3D) { const ExtrusionEntity &entity = multipath3D; default_use(entity); };
//    virtual void use(const ExtrusionLoop &loop) { const ExtrusionEntity &entity = loop; default_use(entity); };
//    virtual void use(const ExtrusionEntityCollection &collection) { const ExtrusionEntity &entity = collection;  default_use(entity); };
//};

class ExtrusionVisitor {
public:
    virtual void default_use(ExtrusionEntity &entity) { assert(false); };
    virtual void use(ExtrusionPath &path);
    virtual void use(ExtrusionPath3D &path3D);
    virtual void use(ExtrusionMultiPath &multipath);
    virtual void use(ExtrusionMultiPath3D &multipath3D);
    virtual void use(ExtrusionLoop &loop);
    virtual void use(ExtrusionEntityCollection &collection);
};
class ExtrusionVisitorConst {
public:
    virtual void default_use(const ExtrusionEntity &entity) { assert(false); };
    virtual void use(const ExtrusionPath &path);
    virtual void use(const ExtrusionPath3D &path3D);
    virtual void use(const ExtrusionMultiPath &multipath);
    virtual void use(const ExtrusionMultiPath3D &multipath3D);
    virtual void use(const ExtrusionLoop &loop);
    virtual void use(const ExtrusionEntityCollection &collection);
};

class ExtrusionEntity
{
protected:
    // even if no_sort, allow to reverse() us (and our entities if they allow it, but they should) 
    bool m_can_reverse;
    ExtrusionEntity(bool can_reverse) : m_can_reverse(can_reverse) {}
public:
    virtual ExtrusionRole role() const = 0;
    virtual bool is_collection() const { return false; }
    virtual bool is_loop() const { return false; }
    virtual bool can_reverse() const { return m_can_reverse; }
    virtual ExtrusionEntity* clone() const = 0;
    // Create a new object, initialize it with this object using the move semantics.
    virtual ExtrusionEntity* clone_move() = 0;
    virtual ~ExtrusionEntity() = default;
    virtual void reverse() = 0;
    virtual const Point& first_point() const = 0;
    virtual const Point& last_point() const = 0;
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion width.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    virtual void polygons_covered_by_width(Polygons &out, const float scaled_epsilon) const = 0;
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion spacing.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    // Useful to calculate area of an infill, which has been really filled in by a 100% rectilinear infill.
    virtual void polygons_covered_by_spacing(Polygons &out, const float spacing_ratio, const float scaled_epsilon) const = 0;
    virtual Polygons polygons_covered_by_width(const float scaled_epsilon = 0.f) const
        { Polygons out; this->polygons_covered_by_width(out, scaled_epsilon); return out; }
    virtual Polygons polygons_covered_by_spacing(const float spacing_ratio, const float scaled_epsilon) const
        { Polygons out; this->polygons_covered_by_spacing(out, spacing_ratio, scaled_epsilon); return out; }
    virtual PolylineOrArc as_polyline() const = 0;
    virtual void   collect_polylines(PolylinesOrArcs &dst) const = 0;
    virtual void   collect_points(Points &dst) const = 0;
    virtual PolylinesOrArcs as_polylines() const { PolylinesOrArcs dst; this->collect_polylines(dst); return dst; }
    virtual double length() const = 0;
    virtual double total_volume() const = 0;
    virtual void visit(ExtrusionVisitor &visitor) = 0;
    virtual void visit(ExtrusionVisitorConst &visitor) const = 0;

    static std::string role_to_string(ExtrusionRole role);
    static ExtrusionRole string_to_role(const std::string_view role);
};

//FIXME: don't use that unsafe container. use a vector of sharedpointer or a ExtrusionEntityCollection.
typedef std::vector<ExtrusionEntity*> ExtrusionEntitiesPtr;

class ExtrusionPath : public ExtrusionEntity
{
public:
    PolylineOrArc polyline;
    // Volumetric velocity. mm^3 of plastic per mm of linear head motion. Used by the G-code generator.
    double mm3_per_mm;
    // Width of the extrusion, used for visualization purposes & for seam notch %. Unscaled
    float width;
    // Height of the extrusion, used for visualization purposes. Unscaled
    float height;

    ExtrusionPath(ExtrusionRole role) : mm3_per_mm(-1), width(-1), height(-1), m_role(role), ExtrusionEntity(true) {}
    ExtrusionPath(ExtrusionRole role, bool can_reverse) : mm3_per_mm(-1), width(-1), height(-1), m_role(role), ExtrusionEntity(can_reverse) {}
    ExtrusionPath(ExtrusionRole role, double mm3_per_mm, float width, float height, bool can_reverse) : mm3_per_mm(mm3_per_mm), width(width), height(height), m_role(role), ExtrusionEntity(can_reverse) { assert(mm3_per_mm == mm3_per_mm); assert(width == width); assert(height == height); }
    ExtrusionPath(const ExtrusionPath& rhs) : polyline(rhs.polyline), mm3_per_mm(rhs.mm3_per_mm), width(rhs.width), height(rhs.height), m_role(rhs.m_role), ExtrusionEntity(rhs.can_reverse()) { assert(mm3_per_mm == mm3_per_mm); assert(width == width); assert(height == height); }
    ExtrusionPath(ExtrusionPath&& rhs) : polyline(std::move(rhs.polyline)), mm3_per_mm(rhs.mm3_per_mm), width(rhs.width), height(rhs.height), m_role(rhs.m_role), ExtrusionEntity(rhs.can_reverse()) { assert(mm3_per_mm == mm3_per_mm); assert(width == width); assert(height == height); }
    ExtrusionPath(const PolylineOrArc& polyline, const ExtrusionPath& rhs) : polyline(polyline), mm3_per_mm(rhs.mm3_per_mm), width(rhs.width), height(rhs.height), m_role(rhs.m_role), ExtrusionEntity(rhs.can_reverse()) { assert(mm3_per_mm == mm3_per_mm); assert(width == width); assert(height == height); }
    ExtrusionPath(PolylineOrArc &&polyline, const ExtrusionPath &rhs) : polyline(std::move(polyline)), mm3_per_mm(rhs.mm3_per_mm), width(rhs.width), height(rhs.height), m_role(rhs.m_role), ExtrusionEntity(rhs.can_reverse()) { assert(mm3_per_mm == mm3_per_mm); assert(width == width); assert(height == height); }

    ExtrusionPath& operator=(const ExtrusionPath& rhs) { m_role = rhs.m_role; this->mm3_per_mm = rhs.mm3_per_mm; this->width = rhs.width; this->height = rhs.height; this->polyline = rhs.polyline; this->m_can_reverse = rhs.m_can_reverse; return *this; }
    ExtrusionPath& operator=(ExtrusionPath&& rhs) { m_role = rhs.m_role; this->mm3_per_mm = rhs.mm3_per_mm; this->width = rhs.width; this->height = rhs.height; this->polyline = std::move(rhs.polyline); this->m_can_reverse = rhs.m_can_reverse; return *this; }

    virtual ExtrusionPath* clone() const override { return new ExtrusionPath(*this); }
    // Create a new object, initialize it with this object using the move semantics.
    virtual ExtrusionPath* clone_move() override { return new ExtrusionPath(std::move(*this)); }
    void reverse() override { this->polyline.reverse(); }
    void set_can_reverse(bool can_reverse) { this->m_can_reverse = can_reverse; }
    const Point& first_point() const override { return this->polyline.front(); }
    const Point& last_point() const override { return this->polyline.back(); }
    size_t size() const { return this->polyline.size(); }
    bool empty() const { return this->polyline.empty(); }
    bool is_closed() const { return ! this->empty() && this->polyline.front() == this->polyline.back(); }
    // Produce a list of extrusion paths into retval by clipping this path by ExPolygonCollection.
    // Currently not used.
    void intersect_expolygons(const ExPolygonCollection &collection, ExtrusionEntityCollection* retval) const;
    // Produce a list of extrusion paths into retval by removing parts of this path by ExPolygonCollection.
    // Currently not used.
    void subtract_expolygons(const ExPolygonCollection &collection, ExtrusionEntityCollection* retval) const;
    void clip_end(coordf_t distance);
    virtual void simplify(coordf_t tolerance, bool with_fitting_arc, double fitting_arc_tolerance);
    double length() const override;
    ExtrusionRole role() const override { return m_role; }
    void set_role(ExtrusionRole new_role) { m_role = new_role; }
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion width.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    void polygons_covered_by_width(Polygons &out, const float scaled_epsilon) const override;
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion spacing.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    // Useful to calculate area of an infill, which has been really filled in by a 100% rectilinear infill.
    void polygons_covered_by_spacing(Polygons &out, const float spacing_ratio, const float scaled_epsilon) const override;
    virtual Polygons polygons_covered_by_width(const float scaled_epsilon = 0.f) const
        { Polygons out; this->polygons_covered_by_width(out, scaled_epsilon); return out; }
    virtual Polygons polygons_covered_by_spacing(const float spacing_ratio, const float scaled_epsilon) const
        { Polygons out; this->polygons_covered_by_spacing(out, spacing_ratio, scaled_epsilon); return out; }
    PolylineOrArc as_polyline() const override { return this->polyline; }
    void   collect_polylines(PolylinesOrArcs &dst) const override { if (! this->polyline.empty()) dst.emplace_back(this->polyline); }
    void   collect_points(Points &dst) const override { append(dst, this->polyline.get_points()); }
    double total_volume() const override { return mm3_per_mm * unscale<double>(length()); }
    virtual void visit(ExtrusionVisitor &visitor) override { visitor.use(*this); };
    virtual void visit(ExtrusionVisitorConst &visitor) const override { visitor.use(*this); };

protected:
    void _inflate_collection(const Polylines &polylines, ExtrusionEntityCollection* collection) const;

    ExtrusionRole m_role;
};
typedef std::vector<ExtrusionPath> ExtrusionPaths;
ExtrusionPaths clip_end(ExtrusionPaths& paths, coordf_t distance);

class ExtrusionPath3D : public ExtrusionPath {
public:
    std::vector<coord_t> z_offsets;

    ExtrusionPath3D(ExtrusionRole role) : ExtrusionPath(role) { /*std::cout << "new path3D\n"; */};
    ExtrusionPath3D(ExtrusionRole role, double mm3_per_mm, float width, float height, bool can_reverse) : ExtrusionPath(role, mm3_per_mm, width, height, can_reverse) { /*std::cout << "new path3D++\n";*/ };
    ExtrusionPath3D(const ExtrusionPath &rhs) : ExtrusionPath(rhs) { /*std::cout << "new path3D from path "<<size()<<"?"<<z_offsets.size()<<"\n";*/ }
    ExtrusionPath3D(ExtrusionPath &&rhs) : ExtrusionPath(rhs) { /*std::cout << "new path3D from path " << size() << "?" << z_offsets.size()<<"\n";*/ }
    ExtrusionPath3D(const ExtrusionPath3D &rhs) : ExtrusionPath(rhs), z_offsets(rhs.z_offsets) { /*std::cout << "new path3D from path3D " << size() << "?" << z_offsets.size()<<"\n";*/ }
    ExtrusionPath3D(ExtrusionPath3D &&rhs) : ExtrusionPath(rhs), z_offsets(std::move(rhs.z_offsets)) { /*std::cout << "new2 path3D from path3D " << size() << "?" << z_offsets.size()<<"\n";*/ }
    //    ExtrusionPath(ExtrusionRole role, const Flow &flow) : m_role(role), mm3_per_mm(flow.mm3_per_mm()), width(flow.width), height(flow.height), feedrate(0.0f), extruder_id(0) {};

    ExtrusionPath3D& operator=(const ExtrusionPath3D &rhs) { m_role = rhs.m_role; this->mm3_per_mm = rhs.mm3_per_mm; this->width = rhs.width; this->height = rhs.height; 
        this->polyline = rhs.polyline; z_offsets = rhs.z_offsets; return *this;
    }
    ExtrusionPath3D& operator=(ExtrusionPath3D &&rhs) { m_role = rhs.m_role; this->mm3_per_mm = rhs.mm3_per_mm; this->width = rhs.width; this->height = rhs.height; 
        this->polyline = std::move(rhs.polyline); z_offsets = std::move(rhs.z_offsets); return *this;
    }
    virtual ExtrusionPath3D* clone() const override { return new ExtrusionPath3D(*this); }
    virtual ExtrusionPath3D* clone_move() override { return new ExtrusionPath3D(std::move(*this)); }
    virtual void visit(ExtrusionVisitor &visitor) override { visitor.use(*this); };
    virtual void visit(ExtrusionVisitorConst &visitor) const override { visitor.use(*this); };

    void push_back(Point p, coord_t z_offset) { 
        assert(!polyline.has_arc());
        polyline.set_points().push_back(p);
        z_offsets.push_back(z_offset);
        polyline.reset_arc();
    }

    //TODO: simplify only for points that have the same z-offset
    void simplify(double tolerance, bool use_arc_fitting, double fitting_arc_tolerance) override;
};
typedef std::vector<ExtrusionPath3D> ExtrusionPaths3D;

// Single continuous extrusion path, possibly with varying extrusion thickness, extrusion height or bridging / non bridging.
template <typename THING = ExtrusionEntity>
class ExtrusionMultiEntity : public ExtrusionEntity {
public:
    std::vector<THING> paths;

    ExtrusionMultiEntity(): ExtrusionEntity(false) {};
    ExtrusionMultiEntity(const ExtrusionMultiEntity &rhs) : paths(rhs.paths), ExtrusionEntity(false) {}
    ExtrusionMultiEntity(ExtrusionMultiEntity &&rhs) : paths(std::move(rhs.paths)), ExtrusionEntity(false) {}
    ExtrusionMultiEntity(const std::vector<THING> &paths) : paths(paths), ExtrusionEntity(false) {};
    ExtrusionMultiEntity(const THING &path): ExtrusionEntity(false) { this->paths.push_back(path); }

    ExtrusionMultiEntity& operator=(const ExtrusionMultiEntity &rhs) { this->paths = rhs.paths; return *this; }
    ExtrusionMultiEntity& operator=(ExtrusionMultiEntity &&rhs) { this->paths = std::move(rhs.paths); return *this; }

    bool is_loop() const override { return false; }
    ExtrusionRole role() const override
    {
        if (this->paths.empty())
            return erNone;
        ExtrusionRole role = this->paths.front().role();
        for (const ExtrusionPath &path : this->paths)
            if (role != path.role()) {
                return erMixed;
            }
        return role;
    }
    virtual const Point& first_point() const override { return this->paths.front().polyline.as_polyline().front(); }
    virtual const Point& last_point() const override { return this->paths.back().polyline.as_polyline().back(); }

    virtual void reverse() override {
        for (THING &entity : this->paths)
            entity.reverse();
        std::reverse(this->paths.begin(), this->paths.end());
    }

    size_t size() const { return this->paths.size(); }
    bool empty() const { return this->paths.empty(); }
    double length() const override {
        double len = 0;
        for (const THING &entity : this->paths)
            len += entity.length();
        return len;
    }

    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion width.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    void polygons_covered_by_width(Polygons &out, const float scaled_epsilon) const override {
        for (const THING &entity : this->paths)
            entity.polygons_covered_by_width(out, scaled_epsilon);
    }

    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion spacing.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    // Useful to calculate area of an infill, which has been really filled in by a 100% rectilinear infill.
    void polygons_covered_by_spacing(Polygons &out, const float spacing_ratio, const float scaled_epsilon) const override {
        for (const THING &entity : this->paths)
            entity.polygons_covered_by_spacing(out, spacing_ratio, scaled_epsilon);
    }

    PolylineOrArc as_polyline() const override {
        PolylineOrArc out;
        if (!paths.empty()) {
            for (const ExtrusionPath& path : paths) {
                out.append(path.as_polyline());
            }
        }
        return out;
    }
    Polygons polygons_covered_by_width(const float scaled_epsilon = 0.f) const override{ Polygons out; this->polygons_covered_by_width(out, scaled_epsilon); return out; }
    Polygons polygons_covered_by_spacing(const float spacing_ratio, const float scaled_epsilon) const override { Polygons out; this->polygons_covered_by_spacing(out, spacing_ratio,  scaled_epsilon); return out; }
    void collect_polylines(PolylinesOrArcs &dst) const override { PolylineOrArc pl = this->as_polyline(); if (!pl.empty()) dst.emplace_back(std::move(pl)); }
    void collect_points(Points &dst) const override { 
        size_t n = std::accumulate(paths.begin(), paths.end(), 0, [](const size_t n, const ExtrusionPath &p){ return n + p.polyline.size(); });
        dst.reserve(dst.size() + n);
        for (const ExtrusionPath &p : this->paths)
            append(dst, p.polyline.get_points());
    }
    double total_volume() const override { double volume = 0.; for (const auto& path : paths) volume += path.total_volume(); return volume; }
};

// Single continuous extrusion path, possibly with varying extrusion thickness, extrusion height or bridging / non bridging.
class ExtrusionMultiPath : public ExtrusionMultiEntity<ExtrusionPath> {
public:

    ExtrusionMultiPath() {};
    ExtrusionMultiPath(const ExtrusionMultiPath &rhs) : ExtrusionMultiEntity(rhs) {}
    ExtrusionMultiPath(ExtrusionMultiPath &&rhs) : ExtrusionMultiEntity(rhs) {}
    ExtrusionMultiPath(const ExtrusionPaths &paths) : ExtrusionMultiEntity(paths) {};
    ExtrusionMultiPath(const ExtrusionPath &path) :ExtrusionMultiEntity(path) {}

    ExtrusionMultiPath& operator=(const ExtrusionMultiPath& rhs) { this->paths = rhs.paths; return *this; }
    ExtrusionMultiPath& operator=(ExtrusionMultiPath&& rhs) { this->paths = std::move(rhs.paths); return *this; }

    void set_can_reverse(bool can_reverse) { m_can_reverse = can_reverse; }

    virtual ExtrusionMultiPath* clone() const override { return new ExtrusionMultiPath(*this); }
    virtual ExtrusionMultiPath* clone_move() override { return new ExtrusionMultiPath(std::move(*this)); }

    virtual void visit(ExtrusionVisitor &visitor) override { visitor.use(*this); };
    virtual void visit(ExtrusionVisitorConst &visitor) const override { visitor.use(*this); };
};
// Single continuous extrusion path, possibly with varying extrusion thickness, extrusion height or bridging / non bridging.
class ExtrusionMultiPath3D : public ExtrusionMultiEntity<ExtrusionPath3D> {
public:

    ExtrusionMultiPath3D() {};
    ExtrusionMultiPath3D(const ExtrusionMultiPath3D &rhs) : ExtrusionMultiEntity(rhs) {}
    ExtrusionMultiPath3D(ExtrusionMultiPath3D &&rhs) : ExtrusionMultiEntity(rhs) {}
    ExtrusionMultiPath3D(const ExtrusionPaths3D &paths) : ExtrusionMultiEntity(paths) {};
    ExtrusionMultiPath3D(const ExtrusionPath3D &path) :ExtrusionMultiEntity(path) {}

    ExtrusionMultiPath3D& operator=(const ExtrusionMultiPath3D& rhs) { this->paths = rhs.paths; return *this; }
    ExtrusionMultiPath3D& operator=(ExtrusionMultiPath3D&& rhs) { this->paths = std::move(rhs.paths); return *this; }

    virtual ExtrusionMultiPath3D* clone() const override { return new ExtrusionMultiPath3D(*this); }
    virtual ExtrusionMultiPath3D* clone_move() override { return new ExtrusionMultiPath3D(std::move(*this)); }

    virtual void visit(ExtrusionVisitor &visitor) override { visitor.use(*this); };
    virtual void visit(ExtrusionVisitorConst &visitor) const override { visitor.use(*this); };

    virtual void reverse() override {
        std::cout << "I SAID NO REVERSE!!!FFFS\n";
    }
};

// Single continuous extrusion loop, possibly with varying extrusion thickness, extrusion height or bridging / non bridging.
class ExtrusionLoop : public ExtrusionEntity
{
public:
    ExtrusionPaths paths;
    
    ExtrusionLoop(ExtrusionLoopRole role = elrDefault) : m_loop_role(role) , ExtrusionEntity(false) {}
    ExtrusionLoop(const ExtrusionPaths &paths, ExtrusionLoopRole role = elrDefault) : paths(paths), m_loop_role(role), ExtrusionEntity(false) { 
        assert(!this->paths.empty());
        assert(this->first_point().coincides_with_epsilon(this->paths.back().polyline.back()));
    }
    ExtrusionLoop(ExtrusionPaths &&paths, ExtrusionLoopRole role = elrDefault) : paths(std::move(paths)), m_loop_role(role), ExtrusionEntity(false) {
        assert(!this->paths.empty());
        assert(this->first_point().coincides_with_epsilon(this->paths.back().polyline.back()));
    }
    ExtrusionLoop(const ExtrusionPath &path, ExtrusionLoopRole role = elrDefault) : m_loop_role(role), ExtrusionEntity(false) {
        this->paths.push_back(path);
        assert(!this->paths.empty());
        assert(this->first_point().coincides_with_epsilon(this->paths.back().polyline.back()));
    }
    ExtrusionLoop(ExtrusionPath &&path, ExtrusionLoopRole role = elrDefault) : m_loop_role(role), ExtrusionEntity(false) {
        this->paths.emplace_back(std::move(path));
        assert(!this->paths.empty());
        assert(this->first_point().coincides_with_epsilon(this->paths.back().polyline.back()));
    }
    virtual bool is_loop() const override{ return true; }
    virtual ExtrusionEntity* clone() const override{ return new ExtrusionLoop (*this); }
    // Create a new object, initialize it with this object using the move semantics.
    virtual ExtrusionEntity* clone_move() override { return new ExtrusionLoop(std::move(*this)); }
    bool make_clockwise();
    bool make_counter_clockwise();
    virtual void reverse() override;
    const Point& first_point() const override { return this->paths.front().polyline.front(); }
    const Point& last_point() const override { assert(this->first_point() == this->paths.back().polyline.back()); return this->first_point(); }
    Polygon polygon() const;
    double length() const override;
    bool split_at_vertex(const Point &point, const coordf_t scaled_epsilon = scale_d(0.001));
    void split_at(const Point &point, bool prefer_non_overhang, const coordf_t scaled_epsilon = scale_d(0.001));
    struct ClosestPathPoint {
        size_t path_idx;
        size_t segment_idx;
        Point  foot_pt;
    };
    ClosestPathPoint get_closest_path_and_point(const Point& point, bool prefer_non_overhang) const;
    // Test, whether the point is extruded by a bridging flow.
    // This used to be used to avoid placing seams on overhangs, but now the EdgeGrid is used instead.
    bool has_overhang_point(const Point &point) const;
    ExtrusionRole role() const override;
    ExtrusionLoopRole loop_role() const { return m_loop_role; }
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion width.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    void polygons_covered_by_width(Polygons &out, const float scaled_epsilon) const override;
    // Produce a list of 2D polygons covered by the extruded paths, offsetted by the extrusion spacing.
    // Increase the offset by scaled_epsilon to achieve an overlap, so a union will produce no gaps.
    // Useful to calculate area of an infill, which has been really filled in by a 100% rectilinear infill.
    void polygons_covered_by_spacing(Polygons &out, const float spacing_ratio, const float scaled_epsilon) const  override;
    Polygons polygons_covered_by_width(const float scaled_epsilon = 0.f) const
        { Polygons out; this->polygons_covered_by_width(out, scaled_epsilon); return out; }
    Polygons polygons_covered_by_spacing(const float spacing_ratio, const float scaled_epsilon) const
        { Polygons out; this->polygons_covered_by_spacing(out, spacing_ratio, scaled_epsilon); return out; }
    PolylineOrArc as_polyline() const override;
    void   collect_polylines(PolylinesOrArcs &dst) const override { PolylineOrArc pl = this->as_polyline(); if (! pl.empty()) dst.emplace_back(std::move(pl)); }
    void   collect_points(Points &dst) const override { 
        size_t n = std::accumulate(paths.begin(), paths.end(), 0, [](const size_t n, const ExtrusionPath &p){ return n + p.polyline.size(); });
        dst.reserve(dst.size() + n);
        for (const ExtrusionPath &p : this->paths)
            append(dst, p.polyline.get_points());
    }
    double total_volume() const override { double volume =0.; for (const auto& path : paths) volume += path.total_volume(); return volume; }

    virtual void visit(ExtrusionVisitor &visitor) override { visitor.use(*this); };
    virtual void visit(ExtrusionVisitorConst &visitor) const override { visitor.use(*this); };

    //static inline std::string role_to_string(ExtrusionLoopRole role);

#ifndef NDEBUG
	bool validate() const {
		assert(this->first_point() == this->paths.back().polyline.back());
		for (size_t i = 1; i < paths.size(); ++ i)
			assert(this->paths[i - 1].polyline.back() == this->paths[i].polyline.front());
		return true;
	}
#endif /* NDEBUG */

private:
    ExtrusionLoopRole m_loop_role{ elrDefault };
};

inline void extrusion_paths_append(ExtrusionPaths &dst, Polylines &polylines, ExtrusionRole role, double mm3_per_mm, float width, float height, bool can_reverse = true)
{
    dst.reserve(dst.size() + polylines.size());
    for (Polyline &polyline : polylines)
        if (polyline.is_valid()) {
            dst.push_back(ExtrusionPath(role, mm3_per_mm, width, height, can_reverse));
            dst.back().polyline = polyline;
        }
}

inline void extrusion_paths_append(ExtrusionPaths &dst, Polylines &&polylines, ExtrusionRole role, double mm3_per_mm, float width, float height, bool can_reverse = true)
{
    dst.reserve(dst.size() + polylines.size());
    for (Polyline &polyline : polylines)
        if (polyline.is_valid()) {
            dst.push_back(ExtrusionPath(role, mm3_per_mm, width, height, can_reverse));
            dst.back().polyline = std::move(polyline);
        }
    polylines.clear();
}

class ExtrusionPrinter : public ExtrusionVisitorConst {
    std::stringstream ss;
    double mult;
    int trunc;
    bool json;
public:
    ExtrusionPrinter(double mult = 0.000001, int trunc = 0, bool json = false) : mult(mult), trunc(trunc), json(json) { }
    virtual void use(const ExtrusionPath& path) override;
    virtual void use(const ExtrusionPath3D& path3D) override;
    virtual void use(const ExtrusionMultiPath& multipath) override;
    virtual void use(const ExtrusionMultiPath3D& multipath) override;
    virtual void use(const ExtrusionLoop& loop) override;
    virtual void use(const ExtrusionEntityCollection& collection) override;
    std::string str() { return ss.str(); }
    std::string print(const ExtrusionEntity& entity)&& {
        entity.visit(*this);
        return ss.str();
    }
};

class ExtrusionLength : public ExtrusionVisitorConst {
    coordf_t dist;
public:
    ExtrusionLength() : dist(0){ }
    virtual void default_use(const ExtrusionEntity& path) override;
    virtual void use(const ExtrusionEntityCollection& collection) override;
    double get() { return dist; }
    double length(const ExtrusionEntity& entity)&& {
        entity.visit(*this);
        return get();
    }
};

class ExtrusionVisitorRecursiveConst : public ExtrusionVisitorConst {
public:
    virtual void use(const ExtrusionMultiPath& multipath) override;
    virtual void use(const ExtrusionMultiPath3D& multipath) override;
    virtual void use(const ExtrusionLoop& loop) override;
    virtual void use(const ExtrusionEntityCollection& collection) override;
};

class ExtrusionVisitorRecursive : public ExtrusionVisitor {
public:
    virtual void use(ExtrusionMultiPath& multipath) override;
    virtual void use(ExtrusionMultiPath3D& multipath) override;
    virtual void use(ExtrusionLoop& loop) override;
    virtual void use(ExtrusionEntityCollection& collection) override;
};

class HasRoleVisitor : public ExtrusionVisitorConst{
public:
    bool found = false;
    void use(const ExtrusionMultiPath& multipath) override;
    void use(const ExtrusionMultiPath3D& multipath3D) override;
    void use(const ExtrusionLoop& loop) override;
    void use(const ExtrusionEntityCollection& collection) override;
    static bool search(const ExtrusionEntity &entity, HasRoleVisitor&& visitor);
    static bool search(const ExtrusionEntitiesPtr &entities, HasRoleVisitor&& visitor);
};
struct HasInfillVisitor : public HasRoleVisitor{
    void default_use(const ExtrusionEntity &entity) override { found = is_infill(entity.role()); };
};
struct HasSolidInfillVisitor : public HasRoleVisitor{
    void default_use(const ExtrusionEntity &entity) override { found = is_solid_infill(entity.role()); };
};
struct HasThisRoleVisitor : public HasRoleVisitor{
    ExtrusionRole role_to_find;
    HasThisRoleVisitor(ExtrusionRole role) : role_to_find(role) {}
    void default_use(const ExtrusionEntity &entity) override { found = entity.role() == role_to_find; };
};


//call simplify for all paths.
class SimplifyVisitor : public ExtrusionVisitorRecursive {
    bool m_use_arc_fitting;
    coordf_t m_scaled_resolution;
    const ConfigOptionFloatOrPercent* m_arc_fitting_tolearance;
public:
    SimplifyVisitor(coordf_t scaled_resolution, bool use_arc_fitting, const ConfigOptionFloatOrPercent* arc_fitting_tolearance) : m_scaled_resolution(scaled_resolution), m_use_arc_fitting(use_arc_fitting), m_arc_fitting_tolearance(arc_fitting_tolearance){}
    virtual void use(ExtrusionPath& path) override {
        path.simplify(m_scaled_resolution, m_use_arc_fitting, scale_d(m_arc_fitting_tolearance->get_abs_value(path.width)));
    }
    virtual void use(ExtrusionPath3D& path3D) override {
        path3D.simplify(m_scaled_resolution, m_use_arc_fitting, scale_d(m_arc_fitting_tolearance->get_abs_value(path3D.width)));
    }
};
class GetPathsVisitor : public ExtrusionVisitorRecursive {
public:
    std::vector<ExtrusionPath*> paths;
    std::vector<ExtrusionPath3D*> paths3D;
    virtual void use(ExtrusionPath& path) override {
        paths.push_back(&path);
    }
    virtual void use(ExtrusionPath3D& path3D) override {
        paths3D.push_back(&path3D);
    }
};

class ExtrusionVolume : public ExtrusionVisitorRecursiveConst {
    bool _with_gap_fill = true;
public:
    double volume = 0; //unscaled
    ExtrusionVolume(bool with_gap_fill = true) : _with_gap_fill(with_gap_fill) {}
    void use(const ExtrusionPath &path) override {
        if(path.role() == erGapFill && !_with_gap_fill) return;
        volume += unscaled(path.length()) * path.mm3_per_mm; }
    void use(const ExtrusionPath3D &path3D) override { volume += unscaled(path3D.length()) * path3D.mm3_per_mm; }
    double get(const ExtrusionEntityCollection &coll);
};

class ExtrusionModifyFlow : public ExtrusionVisitorRecursive {
    double _flow_mult = 1.;
public:
    ExtrusionModifyFlow(double flow_mult) : _flow_mult(flow_mult) {}
    void use(ExtrusionPath &path) override { path.mm3_per_mm *= _flow_mult; path.width *= _flow_mult; }
    void use(ExtrusionPath3D &path3D) override { path3D.mm3_per_mm *= _flow_mult; path3D.width *= _flow_mult; }
    void set(ExtrusionEntityCollection &coll);
};

    
#if _DEBUG
struct LoopAssertVisitor : public ExtrusionVisitorRecursiveConst {
    virtual void default_use(const ExtrusionEntity& entity) override {};
    virtual void use(const ExtrusionPath &path) override { assert(path.length() > SCALED_EPSILON); }
    virtual void use(const ExtrusionLoop &loop) override {
        for (auto it = std::next(loop.paths.begin()); it != loop.paths.end(); ++it) {
            assert(it->polyline.size() >= 2);
            assert(std::prev(it)->polyline.back() == it->polyline.front());
        }
        for (auto it = loop.paths.begin(); it != loop.paths.end(); ++it) {
            assert(it->length() > SCALED_EPSILON);
        }
        assert(loop.paths.front().first_point() == loop.paths.back().last_point());
    }
};
#endif

}

#endif
