#ifndef SRC_LIBSLIC3R_AABBTREELINES_HPP_
#define SRC_LIBSLIC3R_AABBTREELINES_HPP_

#include "Utils.hpp"
#include "libslic3r.h"
#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/Line.hpp"
#include <type_traits>
#include <vector>

namespace Slic3r {

namespace AABBTreeLines {

namespace detail {

template<typename ALineType, typename ATreeType, typename AVectorType>
struct IndexedLinesDistancer {
    using LineType = ALineType;
    using TreeType = ATreeType;
    using VectorType = AVectorType;
    using ScalarType = typename VectorType::Scalar;

    const std::vector<LineType> &lines;
    const TreeType &tree;

    const VectorType origin;

    inline VectorType closest_point_to_origin(size_t primitive_index,
            ScalarType &squared_distance) const {
        Vec<LineType::Dim, typename LineType::Scalar> nearest_point;
        const LineType &line = lines[primitive_index];
        squared_distance = line_alg::distance_to_squared(line, origin.template cast<typename LineType::Scalar>(), &nearest_point);
        return nearest_point.template cast<ScalarType>();
    }
};

}

// Build a balanced AABB Tree over a vector of lines, balancing the tree
// on centroids of the lines.
// Epsilon is applied to the bounding boxes of the AABB Tree to cope with numeric inaccuracies
// during tree traversal.
template<typename LineType>
inline AABBTreeIndirect::Tree<2, typename LineType::Scalar> build_aabb_tree_over_indexed_lines(
        const std::vector<LineType> &lines,
        //FIXME do we want to apply an epsilon?
        const double eps = 0)
        {
    using TreeType = AABBTreeIndirect::Tree<2, typename LineType::Scalar>;
//    using              CoordType      = typename TreeType::CoordType;
    using VectorType = typename TreeType::VectorType;
    using BoundingBox = typename TreeType::BoundingBox;

    struct InputType {
        size_t idx() const {
            return m_idx;
        }
        const BoundingBox& bbox() const {
            return m_bbox;
        }
        const VectorType& centroid() const {
            return m_centroid;
        }

        size_t m_idx;
        BoundingBox m_bbox;
        VectorType m_centroid;
    };

    std::vector<InputType> input;
    input.reserve(lines.size());
    const VectorType veps(eps, eps);
    for (size_t i = 0; i < lines.size(); ++i) {
        const LineType &line = lines[i];
        InputType n;
        n.m_idx = i;
        n.m_centroid = (line.a + line.b) * 0.5;
        n.m_bbox = BoundingBox(line.a, line.a);
        n.m_bbox.extend(line.b);
        n.m_bbox.min() -= veps;
        n.m_bbox.max() += veps;
        input.emplace_back(n);
    }

    TreeType out;
    out.build(std::move(input));
    return out;
}

// Finding a closest line, its closest point and squared distance to the closest point
// Returns squared distance to the closest point or -1 if the input is empty.
// or no closer point than max_sq_dist
template<typename LineType, typename TreeType, typename VectorType, typename Scalar = typename VectorType::Scalar>
inline typename Scalar squared_distance_to_indexed_lines(const std::vector<LineType>        &lines,
                                                         const TreeType                     &tree,
                                                         const VectorType                   &point,
                                                         size_t                             &hit_idx_out,
                                                         Eigen::PlainObjectBase<VectorType> &hit_point_out,
                                                         Scalar max_sqr_dist = std::numeric_limits<Scalar>::infinity())
{
    if (tree.empty()) return Scalar(-1);
    auto distancer = detail::IndexedLinesDistancer<LineType, TreeType, VectorType>{lines, tree, point};
    return AABBTreeIndirect::detail::squared_distance_to_indexed_primitives_recursive(
        distancer, size_t(0), Scalar(0), max_sqr_dist, hit_idx_out, hit_point_out);
}

// Returns all lines within the given radius limit
template<typename LineType, typename TreeType, typename VectorType>
inline std::vector<size_t> all_lines_in_radius(const std::vector<LineType> &lines,
                                               const TreeType              &tree,
                                               const VectorType            &point,
                                               typename VectorType::Scalar max_distance_squared)
{
    auto distancer = detail::IndexedLinesDistancer<LineType, TreeType, VectorType>{lines, tree, point};

    if (tree.empty()) { return {}; }

    std::vector<size_t> found_lines{};
    AABBTreeIndirect::detail::indexed_primitives_within_distance_squared_recurisve(distancer, size_t(0), max_distance_squared, found_lines);
    return found_lines;
}

template<typename LineType> class LinesDistancer
{
private:
    std::vector<LineType> lines;
    using Scalar   = typename LineType::Scalar;
    using Floating = typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type;
    AABBTreeIndirect::Tree<2, Scalar> tree;

public:
    explicit LinesDistancer(const std::vector<LineType> &lines) : lines(lines)
    {
        tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(this->lines);
    }

    // negative sign means inside
    std::tuple<Floating, size_t, Vec<2, Floating>> signed_distance_from_lines_extra(const Vec<2, Scalar> &point) const
    {
        size_t           nearest_line_index_out = size_t(-1);
        Vec<2, Floating> nearest_point_out      = Vec<2, Floating>::Zero();
        Vec<2, Floating> p                      = point.template cast<Floating>();
        auto distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, p, nearest_line_index_out, nearest_point_out);

        if (distance < 0) { return {std::numeric_limits<Floating>::infinity(), nearest_line_index_out, nearest_point_out}; }
        distance              = sqrt(distance);
        const LineType  &line = lines[nearest_line_index_out];
        Vec<2, Floating> v1   = (line.b - line.a).template cast<Floating>();
        Vec<2, Floating> v2   = (point - line.a).template cast<Floating>();
        auto d1 = (v1.x() * v2.y()) - (v1.y() * v2.x()); 
        
        LineType second_line = line;
        if ((line.a.template cast<Floating>() - nearest_point_out).squaredNorm() < SCALED_EPSILON) {
            second_line = lines[prev_idx_modulo(nearest_line_index_out, lines.size())];
        } else {
            second_line = lines[next_idx_modulo(nearest_line_index_out, lines.size())];
        }
        v1   = (second_line.b - second_line.a).template cast<Floating>();
        v2   = (point - second_line.a).template cast<Floating>();
        auto d2 = (v1.x() * v2.y()) - (v1.y() * v2.x()); 
        
        auto d = abs(d1) > abs(d2) ? d1 : d2;
        
        if (d > 0.0) { distance *= -1.0; }
        
        return {distance, nearest_line_index_out, nearest_point_out};
    }

    Floating signed_distance_from_lines(const Vec<2, typename LineType::Scalar> &point) const
    {
        auto [dist, idx, np] = signed_distance_from_lines_extra(point);
        return dist;
    }

    std::vector<size_t> all_lines_in_radius(const Vec<2, typename LineType::Scalar> &point, Floating radius)
    {
        return all_lines_in_radius(this->lines, this->tree, point, radius * radius);
    }

    const LineType &get_line(size_t line_idx) const { return lines[line_idx]; }

    const std::vector<LineType> &get_lines() const { return lines; }
};

}} // namespace Slic3r::AABBTreeLines

#endif /* SRC_LIBSLIC3R_AABBTREELINES_HPP_ */
