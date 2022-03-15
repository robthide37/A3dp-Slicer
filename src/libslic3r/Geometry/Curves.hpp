#ifndef SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_
#define SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_

#include "libslic3r/Point.hpp"

#include <iostream>

namespace Slic3r {
namespace Geometry {

struct PolynomialCurve {
    std::vector<Vec2f> coefficients;

    explicit PolynomialCurve(const std::vector<Vec2f> &coefficients) :
            coefficients(coefficients) {
    }

    Vec3f get_fitted_point(float z) const;
};

//https://towardsdatascience.com/least-square-polynomial-CURVES-using-c-eigen-package-c0673728bd01
// interpolates points in z (treats z coordinates as time) and returns coefficients for axis x and y
PolynomialCurve fit_polynomial(const std::vector<Vec3f> &points, const std::vector<float> &weights, size_t order);

namespace CurveSmoothingKernels {
//NOTE: Kernel functions are used in range [-1,1].

// konstant kernel is used mainly for tests
template<typename NumberType>
struct ConstantKernel {
    float operator()(NumberType w) const {
        return NumberType(1);
    }
};

//https://en.wikipedia.org/wiki/Kernel_(statistics)
template<typename NumberType>
struct EpanechnikovKernel {
    float operator()(NumberType w) const {
        return NumberType(0.25) * (NumberType(1) - w * w);
    }
};

}

template<typename NumberType>
using VecX = Eigen::Matrix<NumberType, Eigen::Dynamic, 1>;

template<int Dimension, typename NumberType, typename Kernel>
struct PiecewiseFittedCurve {
    std::vector<VecX<NumberType>> coefficients;
    Kernel kernel;
    NumberType normalized_kernel_bandwidth;
    NumberType length;
    NumberType segment_size;

    NumberType get_segment_center(size_t segment_index) const {
        return segment_size * segment_index;
    }

    //t should be in normalized range [0,1] with respect to the length of the original polygon line
    Vec<Dimension, NumberType> get_fitted_point(const NumberType &t) const {
        Vec<Dimension, NumberType> result { 0 };
        for (int coeff_index = 0; coeff_index < coefficients[0].size(); ++coeff_index) {
            NumberType segment_center = this->get_segment_center(coeff_index);
            NumberType normalized_center_dis = (segment_center - t)
                    / (NumberType(0.5) * normalized_kernel_bandwidth);
            if (normalized_center_dis >= NumberType(-1) && normalized_center_dis <= NumberType(1)) {
                for (size_t dim = 0; dim < Dimension; ++dim) {
                    result[dim] += kernel(normalized_center_dis) * coefficients[dim][coeff_index];
                }
            }
        }
        return result;
    }
};

// number_of_segments: how many curve segments (kernel applications) are used. Must be at least 2, because we are centering the segments on the first and last point
// normalized_kernel_bandwidth (0,..): how spread the kernel is over the points in normalized coordinates (points are mapped to parametric range [0,1])
// for example, 0.5 normalized_kernel_bandwidth means that one curve segment covers half of the points
template<int Dimension, typename NumberType, typename Kernel>
PiecewiseFittedCurve<Dimension, NumberType, Kernel> fit_curve(const std::vector<Vec<Dimension, NumberType>> &points,
        size_t number_of_segments, const NumberType &normalized_kernel_bandwidth,
        Kernel kernel) {

    // check to make sure inputs are correct
    assert(normalized_kernel_bandwidth > 0);
    assert(number_of_segments >= 2);
    assert(number_of_segments <= points.size());
    NumberType length { };
    std::vector<NumberType> knots { };
    for (size_t point_index = 0; point_index < points.size() - 1; ++point_index) {
        knots.push_back(length);
        length += (points[point_index + 1] - points[point_index]).norm();
    }
    //last point
    knots.push_back(length);

    //normalize knots
    for (NumberType &knot : knots) {
        knot = knot / length;
    }

    PiecewiseFittedCurve<Dimension, NumberType, Kernel> result { };

    result.kernel = kernel;
    result.normalized_kernel_bandwidth = normalized_kernel_bandwidth;
    result.length = length;
    result.segment_size = NumberType(1) / (number_of_segments - NumberType(1));

    std::vector<VecX<NumberType>> data_points(Dimension);
    for (size_t dim = 0; dim < Dimension; ++dim) {
        data_points[dim] = Eigen::Matrix<NumberType, Eigen::Dynamic, 1>(points.size());
    }

    for (size_t index = 0; index < points.size(); index++) {
        for (size_t dim = 0; dim < Dimension; ++dim) {
            data_points[dim](index) = points[index](dim);
        }
    }

    Eigen::MatrixXf T(knots.size(), number_of_segments);
    for (size_t i = 0; i < knots.size(); ++i) {
        for (size_t j = 0; j < number_of_segments; ++j) {
            NumberType knot_val = knots[i];
            NumberType segment_center = result.get_segment_center(j);
            NumberType normalized_center_dis = (segment_center - knot_val)
                    / (NumberType(0.5) * normalized_kernel_bandwidth);
            if (normalized_center_dis >= NumberType(-1) && normalized_center_dis <= NumberType(1)) {
                T(i, j) = kernel(normalized_center_dis);
            } else {
                T(i, j) = 0;
            }
        }
    }

    // Solve for linear least square fit
    std::vector<VecX<NumberType>> coefficients(Dimension);
    const auto QR = T.colPivHouseholderQr();
    for (size_t dim = 0; dim < Dimension; ++dim) {
        coefficients[dim] = QR.solve(data_points[dim]);
    }

    result.coefficients = coefficients;

    return result;
}

template<int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, CurveSmoothingKernels::EpanechnikovKernel<NumberType>>
fit_epanechnikov_curve(const std::vector<Vec<Dimension, NumberType>> &points,
        size_t number_of_segments, const NumberType &normalized_kernel_bandwidth) {
    return fit_curve(points, number_of_segments, normalized_kernel_bandwidth,
            CurveSmoothingKernels::EpanechnikovKernel<NumberType> { });
}

}
}

#endif /* SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_ */
