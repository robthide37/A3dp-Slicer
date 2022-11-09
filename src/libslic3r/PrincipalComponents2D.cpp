#include "PrincipalComponents2D.hpp"
#include "Point.hpp"

namespace Slic3r {

// returns two eigenvectors of the area covered by given polygons. The vectors are sorted by their corresponding eigenvalue, largest first
std::tuple<Vec2d, Vec2d> compute_principal_components(const Polygons &polys, const double unscaled_resolution)
{
    // USING UNSCALED VALUES
    const Vec2d pixel_size  = Vec2d(unscaled_resolution, unscaled_resolution);
    const auto  bb          = get_extents(polys);
    const Vec2i pixel_count = unscaled(bb.size()).cwiseQuotient(pixel_size).cast<int>() + Vec2i::Ones();

    std::vector<Linef> lines{};
    for (Line l : to_lines(polys)) { lines.emplace_back(unscaled(l.a), unscaled(l.b)); }
    AABBTreeIndirect::Tree<2, double> tree      = AABBTreeLines::build_aabb_tree_over_indexed_lines(lines);
    auto                              is_inside = [&](const Vec2d &point) {
        size_t nearest_line_index_out = 0;
        Vec2d  nearest_point_out      = Vec2d::Zero();
        auto   distance = AABBTreeLines::squared_distance_to_indexed_lines(lines, tree, point, nearest_line_index_out, nearest_point_out);
        if (distance < 0) return false;
        const Linef &line = lines[nearest_line_index_out];
        Vec2d        v1   = line.b - line.a;
        Vec2d        v2   = point - line.a;
        if ((v1.x() * v2.y()) - (v1.y() * v2.x()) > 0.0) { return true; }
        return false;
    };

    double pixel_area                                   = pixel_size.x() * pixel_size.y();
    Vec2d  centroid_accumulator                         = Vec2d::Zero();
    Vec2d  second_moment_of_area_accumulator            = Vec2d::Zero();
    double second_moment_of_area_covariance_accumulator = 0.0;
    double area                                         = 0.0;

    for (int x = 0; x < pixel_count.x(); x++) {
        for (int y = 0; y < pixel_count.y(); y++) {
            Vec2d position = unscaled(bb.min) + pixel_size.cwiseProduct(Vec2d{x, y});
            if (is_inside(position)) {
                area += pixel_area;
                centroid_accumulator += pixel_area * position;
                second_moment_of_area_accumulator += pixel_area * position.cwiseProduct(position);
                second_moment_of_area_covariance_accumulator += pixel_area * position.x() * position.y();
            }
        }
    }

    if (area <= 0.0) {
        return {Vec2d::Zero(), Vec2d::Zero()};
    }

    Vec2d  centroid   = centroid_accumulator / area;
    Vec2d  variance   = second_moment_of_area_accumulator / area - centroid.cwiseProduct(centroid);
    double covariance = second_moment_of_area_covariance_accumulator / area - centroid.x() * centroid.y();
#if 0
        std::cout << "area : " << area << std::endl;
        std::cout << "variancex : " << variance.x() << std::endl;
        std::cout << "variancey : " << variance.y() << std::endl;
        std::cout << "covariance : " << covariance << std::endl;
#endif
    if (abs(covariance) < EPSILON) {
        std::tuple<Vec2d, Vec2d> result{Vec2d{variance.x(), 0.0}, Vec2d{0.0, variance.y()}};
        if (variance.y() > variance.x()) {
            return {std::get<1>(result), std::get<0>(result)};
        } else
            return result;
    }

    // now we find the first principal component of the covered area by computing max eigenvalue and the correspoding eigenvector of
    // covariance matrix
    //  covaraince matrix C is :  | VarX  Cov  |
    //                            | Cov   VarY |
    // Eigenvalues are solutions to det(C - lI) = 0, where l is the eigenvalue and I unit matrix
    // Eigenvector for eigenvalue l is any vector v such that Cv = lv

    double eigenvalue_a = 0.5 * (variance.x() + variance.y() +
                                 sqrt((variance.x() - variance.y()) * (variance.x() - variance.y()) + 4 * covariance * covariance));
    double eigenvalue_b = 0.5 * (variance.x() + variance.y() -
                                 sqrt((variance.x() - variance.y()) * (variance.x() - variance.y()) + 4 * covariance * covariance));
    Vec2d  eigenvector_a{(eigenvalue_a - variance.y()) / covariance, 1.0};
    Vec2d  eigenvector_b{(eigenvalue_b - variance.y()) / covariance, 1.0};

#if 0
        std::cout << "eigenvalue_a: " << eigenvalue_a << std::endl;
        std::cout << "eigenvalue_b: " << eigenvalue_b << std::endl;
        std::cout << "eigenvectorA: " << eigenvector_a.x() <<  "  " << eigenvector_a.y() << std::endl;
        std::cout << "eigenvectorB: " << eigenvector_b.x() <<  "  " << eigenvector_b.y() << std::endl;
#endif

    if (eigenvalue_a > eigenvalue_b) {
        return {eigenvector_a, eigenvector_b};
    } else {
        return {eigenvector_b, eigenvector_a};
    }
}

}