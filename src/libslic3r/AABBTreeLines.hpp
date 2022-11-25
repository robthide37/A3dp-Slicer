#ifndef SRC_LIBSLIC3R_AABBTREELINES_HPP_
#define SRC_LIBSLIC3R_AABBTREELINES_HPP_

#include "Utils.hpp"
#include "libslic3r.h"
#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/Line.hpp"
#include <type_traits>
#include <vector>

namespace Slic3r { namespace AABBTreeLines {

namespace detail {

template<typename ALineType, typename ATreeType, typename AVectorType> struct IndexedLinesDistancer
{
    using LineType   = ALineType;
    using TreeType   = ATreeType;
    using VectorType = AVectorType;
    using ScalarType = typename VectorType::Scalar;

    const std::vector<LineType> &lines;
    const TreeType              &tree;

    const VectorType origin;

    inline VectorType closest_point_to_origin(size_t primitive_index, ScalarType &squared_distance) const
    {
        Vec<LineType::Dim, typename LineType::Scalar> nearest_point;
        const LineType                               &line = lines[primitive_index];
        squared_distance = line_alg::distance_to_squared(line, origin.template cast<typename LineType::Scalar>(), &nearest_point);
        return nearest_point.template cast<ScalarType>();
    }
};

// !!! NOTE: algorithm expects the BoundingBoxes to be snug, no epsilon is allowed
template<typename LineType, typename TreeType, typename VectorType, int coordinate>
inline std::tuple<size_t, size_t> coordinate_aligned_ray_hit_count(size_t node_idx, const TreeType &tree, const VectorType &ray_origin)
{
    static constexpr int other_coordinate = (coordinate + 1) % 2;
    using Scalar                          = typename LineType::Scalar;
    using Floating                        = typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type;
    const auto &node                      = tree.node(node_idx);
    assert(node.is_valid());
    if (node.is_leaf()) {
        if (ray_origin[coordinate] > node.bbox.max()[coordinate]) {
            return {1, 0};
        } else if (ray_origin[coordinate] < node.bbox.min()[coordinate]) {
            return {0, 1};
        } else {
            auto     sizes = node.bbox.sizes();
            Floating t     = (ray_origin[other_coordinate] - node.bbox.min()[other_coordinate]) /
                         (sizes[other_coordinate] > 0 ? sizes[other_coordinate] : 1);
            auto intersection = node.bbox.min()[coordinate] + t * sizes[coordinate];
            if (ray_origin[coordinate] > intersection) {
                return {1, 0};
            } else if (ray_origin[coordinate] < intersection) {
                return {0, 1};
            } else { // ray origin is on boundary
                return {1, 1};
            }
        }
    } else {
        size_t      intersections_above = 0;
        size_t      intersections_below = 0;
        size_t      left_node_idx       = node_idx * 2 + 1;
        size_t      right_node_idx      = left_node_idx + 1;
        const auto &node_left           = tree.node(left_node_idx);
        const auto &node_right          = tree.node(right_node_idx);
        assert(node_left.is_valid());
        assert(node_right.is_valid());

        if (node_left.bbox.min()[other_coordinate] <= ray_origin[other_coordinate] &&
            node_left.bbox.max()[other_coordinate] > ray_origin[other_coordinate]) { // sharp inequality, beacuse we do not want to count point common to two lines more than once
            auto [above, below] = coordinate_aligned_ray_hit_count<LineType, TreeType, VectorType, coordinate>(left_node_idx, tree,
                                                                                                               ray_origin);
            intersections_above += above;
            intersections_below += below;
        }

        if (node_right.bbox.min()[other_coordinate] <= ray_origin[other_coordinate] &&
            node_right.bbox.max()[other_coordinate] > ray_origin[other_coordinate]) {
            auto [above, below] = coordinate_aligned_ray_hit_count<LineType, TreeType, VectorType, coordinate>(right_node_idx, tree,
                                                                                                               ray_origin);
            intersections_above += above;
            intersections_below += below;
        }
        return {intersections_above, intersections_below};
    }
}

} // namespace detail

// Build a balanced AABB Tree over a vector of lines, balancing the tree
// on centroids of the lines.
// Epsilon is applied to the bounding boxes of the AABB Tree to cope with numeric inaccuracies
// during tree traversal.
template<typename LineType>
inline AABBTreeIndirect::Tree<2, typename LineType::Scalar> build_aabb_tree_over_indexed_lines(const std::vector<LineType> &lines)
{
    using TreeType = AABBTreeIndirect::Tree<2, typename LineType::Scalar>;
    //    using              CoordType      = typename TreeType::CoordType;
    using VectorType  = typename TreeType::VectorType;
    using BoundingBox = typename TreeType::BoundingBox;

    struct InputType
    {
        size_t             idx() const { return m_idx; }
        const BoundingBox &bbox() const { return m_bbox; }
        const VectorType  &centroid() const { return m_centroid; }

        size_t      m_idx;
        BoundingBox m_bbox;
        VectorType  m_centroid;
    };

    std::vector<InputType> input;
    input.reserve(lines.size());
    for (size_t i = 0; i < lines.size(); ++i) {
        const LineType &line = lines[i];
        InputType       n;
        n.m_idx      = i;
        n.m_centroid = (line.a + line.b) * 0.5;
        n.m_bbox     = BoundingBox(line.a, line.a);
        n.m_bbox.extend(line.b);
        input.emplace_back(n);
    }

    TreeType out;
    out.build(std::move(input));
    return out;
}

// Finding a closest line, its closest point and squared distance to the closest point
// Returns squared distance to the closest point or -1 if the input is empty.
// or no closer point than max_sq_dist
template<typename LineType, typename TreeType, typename VectorType>
inline typename VectorType::Scalar squared_distance_to_indexed_lines(
    const std::vector<LineType>        &lines,
    const TreeType                     &tree,
    const VectorType                   &point,
    size_t                             &hit_idx_out,
    Eigen::PlainObjectBase<VectorType> &hit_point_out,
    typename VectorType::Scalar         max_sqr_dist = std::numeric_limits<typename VectorType::Scalar>::infinity())
{
    using Scalar = typename VectorType::Scalar;
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
                                               typename VectorType::Scalar  max_distance_squared)
{
    auto distancer = detail::IndexedLinesDistancer<LineType, TreeType, VectorType>{lines, tree, point};

    if (tree.empty()) { return {}; }

    std::vector<size_t> found_lines{};
    AABBTreeIndirect::detail::indexed_primitives_within_distance_squared_recurisve(distancer, size_t(0), max_distance_squared, found_lines);
    return found_lines;
}

// return 1 if true, -1 if false, 0 if cannot be determined
template<typename LineType, typename TreeType, typename VectorType>
inline int point_outside_closed_contours(const std::vector<LineType> &lines, const TreeType &tree, const VectorType &point)
{
    if (tree.empty()) { return 1; }

    auto [hits_above, hits_below] = detail::coordinate_aligned_ray_hit_count<LineType, TreeType, VectorType, 0>(0, tree, point);
    std::cout << "hits_above:  " << hits_above << "   hits_below: " << hits_below << std::endl; 
    if (hits_above % 2 == 1 && hits_below % 2 == 1) {
        return -1;
    } else if (hits_above % 2 == 0 && hits_below % 2 == 0) {
        return 1;
    } else { // this should not happen with closed contours. lets check it in Y direction
        auto [hits_above, hits_below] = detail::coordinate_aligned_ray_hit_count<LineType, TreeType, VectorType, 1>(0, tree, point);
        if (hits_above % 2 == 1 && hits_below % 2 == 1) {
            return -1;
        } else if (hits_above % 2 == 0 && hits_below % 2 == 0) {
            return 1;
        } else { // both results were unclear
            return 0;
        }
    }
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

    explicit LinesDistancer(std::vector<LineType> &&lines) : lines(lines)
    {
        tree = AABBTreeLines::build_aabb_tree_over_indexed_lines(this->lines);
    }

    LinesDistancer() = default;

    // 1 true, -1 false, 0 cannot determine
    int outside(const Vec<2, Scalar> &point) const { return point_outside_closed_contours(lines, tree, point); }

    // negative sign means inside
    std::tuple<Floating, size_t, Vec<2, Floating>> signed_distance_from_lines_extra(const Vec<2, Scalar> &point) const
    {
        size_t           nearest_line_index_out = size_t(-1);
        Vec<2, Floating> nearest_point_out      = Vec<2, Floating>::Zero();
        Vec<2, Floating> p                      = point.template cast<Floating>();
        auto distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, p, nearest_line_index_out, nearest_point_out);

        if (distance < 0) { return {std::numeric_limits<Floating>::infinity(), nearest_line_index_out, nearest_point_out}; }
        distance = sqrt(distance);
        distance *= outside(point);

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
