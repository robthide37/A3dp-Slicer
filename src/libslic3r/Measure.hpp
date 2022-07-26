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
    SurfaceFeatureType      get_type() const override { return SurfaceFeatureType::Edge; }
    std::pair<Vec3d, Vec3d> get_edge() const { return std::make_pair(m_start, m_end); }
private:
    Vec3d m_start;
    Vec3d m_end;
};

class Circle : public SurfaceFeature {
public:
    Circle(const Vec3d& center, double radius) : m_center{center}, m_radius{radius} {}
    SurfaceFeatureType   get_type()   const override { return SurfaceFeatureType::Circle; }
    Vec3d  get_center() const { return m_center; }
    double get_radius() const { return m_radius; }
private:
    Vec3d m_center;
    double m_radius;
};

class Plane : public SurfaceFeature {
public:
    SurfaceFeatureType get_type() const override { return SurfaceFeatureType::Plane; }

};


class MeasuringImpl;


class Measuring {
public:
    // Construct the measurement object on a given its. The its must remain
    // valid and unchanged during the whole lifetime of the object.
    explicit Measuring(const indexed_triangle_set& its);
    ~Measuring();
    
    // Return a reference to a list of all features identified on the its.
    const std::vector<SurfaceFeature*>& get_features() const;

    // Given a face_idx where the mouse cursor points, return a feature that
    // should be highlighted or nullptr.
    const SurfaceFeature* get_feature(size_t face_idx, const Vec3d& point) const;

    // Returns distance between two SurfaceFeatures.
    static double get_distance(const SurfaceFeature* a, const SurfaceFeature* b);

    // Returns true if an x/y/z distance between features makes sense.
    // If so, result contains the distances.
    static bool   get_distances(const SurfaceFeature* a, const SurfaceFeature* b, std::array<double, 3>& result);

    // Returns true if an x/y/z distance between feature and a point makes sense.
    // If so, result contains the distances.
    static bool   get_axis_aligned_distances(const SurfaceFeature* feature, const Vec3d* pt, std::array<double, 3>& result);

    // Returns true if measuring angles between features makes sense.
    // If so, result contains the angle in radians.
    static bool   get_angle(const SurfaceFeature* a, const SurfaceFeature* b, double& result);


private:
    
    std::unique_ptr<MeasuringImpl> priv;  
};


} // namespace Measure
} // namespace Slic3r

#endif // Slic3r_Measure_hpp_
