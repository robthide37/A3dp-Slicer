#ifndef Slic3r_Measure_hpp_
#define Slic3r_Measure_hpp_

#include <memory>

#include "Point.hpp"


struct indexed_triangle_set;



namespace Slic3r {
namespace Measure {


enum class SurfaceFeatureType {
        Edge      = 1 << 0,
        Circle    = 1 << 1,
        Plane     = 1 << 2
    };

class SurfaceFeature {
public:        
    virtual SurfaceFeatureType get_type() const = 0;
};

class Edge : public SurfaceFeature {
public:
    Edge(const Vec3d& start, const Vec3d& end) : m_start{start}, m_end{end} {}
    Edge(const Vec3d& start, const Vec3d& end, const Vec3d& pin) : m_start{start}, m_end{end},
        m_pin{std::unique_ptr<Vec3d>(new Vec3d(pin))} {}
    SurfaceFeatureType      get_type() const override { return SurfaceFeatureType::Edge; }
    std::pair<Vec3d, Vec3d> get_edge() const { return std::make_pair(m_start, m_end); }
    const Vec3d*            get_point_of_interest() const { return m_pin.get(); }
private:
    Vec3d m_start;
    Vec3d m_end;
    std::unique_ptr<Vec3d> m_pin;
};

class Circle : public SurfaceFeature {
public:
    Circle(const Vec3d& center, double radius, const Vec3d& normal)
    : m_center{center}, m_radius{radius}, m_normal{normal} {}
    SurfaceFeatureType   get_type()   const override { return SurfaceFeatureType::Circle; }
    Vec3d  get_center() const { return m_center; }
    double get_radius() const { return m_radius; }
    Vec3d  get_normal() const { return m_normal; }
private:
    Vec3d m_center;
    double m_radius;
    Vec3d m_normal;
};

class Plane : public SurfaceFeature {
public:
    Plane(int idx) : m_idx(idx) {}
    SurfaceFeatureType get_type() const override { return SurfaceFeatureType::Plane; }
    int get_plane_idx() const { return m_idx; } // index into vector provided by Measuring::get_plane_triangle_indices

private:
    int m_idx;
};


class MeasuringImpl;


class Measuring {
public:
    // Construct the measurement object on a given its. The its must remain
    // valid and unchanged during the whole lifetime of the object.
    explicit Measuring(const indexed_triangle_set& its);
    ~Measuring();
    
    // Return a reference to a list of all features identified on the its.
    [[deprecated]]const std::vector<const SurfaceFeature*>& get_features() const;

    // Given a face_idx where the mouse cursor points, return a feature that
    // should be highlighted or nullptr.
    const SurfaceFeature* get_feature(size_t face_idx, const Vec3d& point) const;

    // Returns a list of triangle indices for each identified plane. Each
    // Plane object contains an index into this vector.
    const std::vector<std::vector<int>> get_planes_triangle_indices() const;



    // Returns distance between two SurfaceFeatures.
    static double get_distance(const SurfaceFeature* a, const SurfaceFeature* b);

    // Returns distance between a SurfaceFeature and a point.
    static double get_distance(const SurfaceFeature* a, const Vec3d* pt);

    // Returns true if measuring angles between features makes sense.
    // If so, result contains the angle in radians.
    static bool   get_angle(const SurfaceFeature* a, const SurfaceFeature* b, double& result);


private: 
    std::unique_ptr<MeasuringImpl> priv;  
};


} // namespace Measure
} // namespace Slic3r

#endif // Slic3r_Measure_hpp_
