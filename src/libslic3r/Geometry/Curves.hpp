#ifndef SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_
#define SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_

#include "libslic3r/Point.hpp"
#include "Bicubic.hpp"

#include <iostream>

namespace Slic3r {
namespace Geometry {

template<int Dimension, typename NumberType>
struct PolynomialCurve {
    std::vector<DynVec<NumberType>> coefficients;

    explicit PolynomialCurve(std::vector<DynVec<NumberType>> coefficients) :
            coefficients(coefficients) {
    }

    Vec3f get_fitted_value(const NumberType value) const {
        Vec<Dimension, NumberType> result = Vec<Dimension, NumberType>::Zero();
        size_t order = this->coefficients.size() - 1;
        for (size_t index = 0; index < order + 1; ++index) {
            float powered = pow(value, index);
            for (size_t dim = 0; dim < Dimension; ++dim) {
                result(dim) += powered * this->coefficients[dim](index);
            }
        }
        return result;
    }
};

//https://towardsdatascience.com/least-square-polynomial-CURVES-using-c-eigen-package-c0673728bd01
// interpolates points in z (treats z coordinates as time) and returns coefficients for axis x and y
template<int Dimension, typename NumberType>
PolynomialCurve<Dimension, NumberType> fit_polynomial(const std::vector<Vec<Dimension, NumberType>> &observations,
        const std::vector<NumberType> &observation_points,
        const std::vector<NumberType> &weights, size_t order) {
    // check to make sure inputs are correct
    assert(observation_points.size() >= order + 1);
    assert(observation_points.size() == weights.size());
    assert(observations.size() == weights.size());

    std::vector<float> squared_weights(weights.size());
    for (size_t index = 0; index < weights.size(); ++index) {
        squared_weights[index] = sqrt(weights[index]);
    }

    std::vector<DynVec<NumberType>> data_points(Dimension);
    for (size_t dim = 0; dim < Dimension; ++dim) {
        data_points[dim] = Eigen::Matrix<NumberType, Eigen::Dynamic, 1>(
                observations.size());
    }
    for (size_t index = 0; index < observations.size(); index++) {
        for (size_t dim = 0; dim < Dimension; ++dim) {
            data_points[dim](index) = observations[index](dim) * squared_weights[index];
        }
    }

    Eigen::MatrixXf T(observation_points.size(), order + 1);
    // Populate the matrix
    for (size_t i = 0; i < observation_points.size(); ++i) {
        for (size_t j = 0; j < order + 1; ++j) {
            T(i, j) = pow(observation_points[i], j) * squared_weights[i];
        }
    }

    const auto QR = T.householderQr();
    std::vector<DynVec<NumberType>> coefficients(Dimension);
    // Solve for linear least square fit
    for (size_t dim = 0; dim < Dimension; ++dim) {
        coefficients[dim] = QR.solve(data_points[dim]);
    }

    return PolynomialCurve<Dimension, NumberType>(coefficients);
}

template<size_t Dimension, typename NumberType, typename Kernel>
struct PiecewiseFittedCurve {
    std::vector<DynVec<NumberType>> coefficients;
    Kernel kernel;
    NumberType start;
    NumberType length;
    NumberType n_segment_size;
    size_t segments_count;

    NumberType get_n_segment_start(int segment_index) const {
        return n_segment_size * segment_index;
    }

    NumberType normalize(const NumberType &observation_point) const {
        return (observation_point - start) / length;
    }

    Vec<Dimension, NumberType> get_fitted_value(const NumberType &observation_point) const {
        Vec<Dimension, NumberType> result = Vec<Dimension, NumberType>::Zero();
        NumberType t = normalize(observation_point);

        int start_segment_idx = int(floor(t / this->n_segment_size)) - int(Kernel::kernel_span * 0.5f - 1.0f);
        for (int segment_index = start_segment_idx; segment_index < int(start_segment_idx + Kernel::kernel_span);
                segment_index++) {
            if (segment_index < 0 || segment_index >= int(this->segments_count)) {
                continue;
            }
            NumberType segment_start = this->get_n_segment_start(segment_index);
            NumberType normalized_segment_distance = (segment_start - t) / this->n_segment_size;

            for (size_t dim = 0; dim < Dimension; ++dim) {
                result(dim) += kernel.kernel(normalized_segment_distance) * coefficients[dim](segment_index);
            }
        }
        return result;
    }
}
;

// observations: data to be fitted by the curve
// observation points: growing sequence of points where the observations were made.
//      In other words, for function f(x) = y, observations are y0...yn, and observation points are x0...xn
// weights: how important the observation is
//  number_of_inner_splines: how many full inner splines are fit into the normalized valid range 0,1;
//          final number of knots is Kernel::kernel_span times larger + additional segments on start and end
// Kernel: model used for the curve fitting
template<int Dimension, typename NumberType, typename Kernel>
PiecewiseFittedCurve<Dimension, NumberType, Kernel> fit_curve(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        const std::vector<NumberType> &observation_points,
        const std::vector<NumberType> &weights,
        size_t number_of_inner_splines,
        Kernel kernel) {

    // check to make sure inputs are correct
    assert(number_of_inner_splines >= 1);
    assert(observation_points.size() == observations.size());
    assert(observation_points.size() == weights.size());
    assert(number_of_inner_splines <= observations.size());
    assert(observations.size() >= 4);

    size_t extremes_repetition = Kernel::kernel_span - 1; //how many (additional) times is the first and last point repeated

    //prepare sqrt of weights, which will then be applied to both matrix T and observed data: https://en.wikipedia.org/wiki/Weighted_least_squares
    std::vector<float> sqrt_weights(weights.size() + extremes_repetition * 2);
    for (size_t index = 0; index < weights.size(); ++index) {
        assert(weights[index] > 0);
        sqrt_weights[index + extremes_repetition] = sqrt(weights[index]);
    }
    //repeat weights for addtional segments
    for (int index = 0; index < int(extremes_repetition); ++index) {
        sqrt_weights[index] = sqrt(weights.front());
        sqrt_weights[sqrt_weights.size() - index - 1] = sqrt(weights.back());
    }

    // prepare result and compute metadata
    PiecewiseFittedCurve<Dimension, NumberType, Kernel> result { };

    NumberType original_len = observation_points.back() - observation_points.front();
    NumberType orig_segment_size = original_len / NumberType(number_of_inner_splines * Kernel::kernel_span);
    result.kernel = kernel;
    result.start = observation_points.front() - extremes_repetition * orig_segment_size;
    result.length = observation_points.back() + extremes_repetition * orig_segment_size - result.start;
    result.segments_count = number_of_inner_splines * Kernel::kernel_span + extremes_repetition * 2;
    result.n_segment_size = NumberType(1) / NumberType(result.segments_count - 1);

    //normalize observations points by start and length
    std::vector<NumberType> normalized_obs_points(observation_points.size() + extremes_repetition * 2);
    for (size_t index = 0; index < observation_points.size(); ++index) {
        normalized_obs_points[index + extremes_repetition] = result.normalize(observation_points[index]);
    }
    // create artificial values at the extremes for constant curve fitting
    for (int index = 0; index < int(extremes_repetition); ++index) {
        normalized_obs_points[extremes_repetition - 1 - index] = result.normalize(observation_points.front()
                - index * orig_segment_size - NumberType(0.5) * orig_segment_size);

        normalized_obs_points[normalized_obs_points.size() - extremes_repetition + index] = result.normalize(
                observation_points.back() + index * orig_segment_size + NumberType(0.5) * orig_segment_size);
    }

    // prepare observed data
    std::vector<DynVec<NumberType>> data_points(Dimension);
    for (size_t dim = 0; dim < Dimension; ++dim) {
        data_points[dim] = Eigen::Matrix<NumberType, Eigen::Dynamic, 1>(
                observations.size() + extremes_repetition * 2);
    }
    for (size_t index = 0; index < observations.size(); index++) {
        for (size_t dim = 0; dim < Dimension; ++dim) {
            data_points[dim](index + extremes_repetition) = observations[index](dim)
                    * sqrt_weights[index + extremes_repetition];
        }
    }
    //duplicate observed data at the start and end
    for (int index = 0; index < int(extremes_repetition); index++) {
        for (size_t dim = 0; dim < Dimension; ++dim) {
            data_points[dim](index) = observations.front()(dim) * sqrt_weights[index];
            data_points[dim](data_points[dim].size() - index - 1) = observations.back()(dim)
                    * sqrt_weights[data_points[dim].size() - index - 1];
        }
    }

    //Create weight matrix for each point and each segment.
    Eigen::MatrixXf T(normalized_obs_points.size(), result.segments_count);
    for (size_t i = 0; i < normalized_obs_points.size(); ++i) {
        for (size_t j = 0; j < result.segments_count; ++j) {
            T(i, j) = NumberType(0);
        }
    }

    for (size_t i = 0; i < normalized_obs_points.size(); ++i) {
        NumberType knot_val = normalized_obs_points[i];
        int start_segment_idx = int(floor(knot_val / result.n_segment_size)) - int(Kernel::kernel_span * 0.5f - 1.0f);
        for (int segment_index = start_segment_idx; segment_index < int(start_segment_idx + Kernel::kernel_span);
                segment_index++) {
            if (segment_index < 0 || segment_index >= int(result.segments_count)) {
                continue;
            }
            NumberType segment_start = result.get_n_segment_start(segment_index);
            NumberType normalized_segment_distance = (segment_start - knot_val) / result.n_segment_size;
            // fill in kernel value with weight applied

            T(i, segment_index) += kernel.kernel(normalized_segment_distance) * sqrt_weights[i];
        }
    }

    // Solve for linear least square fit
    std::vector<DynVec<NumberType>> coefficients(Dimension);
    const auto QR = T.fullPivHouseholderQr();
    for (size_t dim = 0; dim < Dimension; ++dim) {
        coefficients[dim] = QR.solve(data_points[dim]);
    }

    // store coefficients in result
    result.coefficients = coefficients;

    return result;
}

template<int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, CubicBSplineKernel<NumberType>>
fit_cubic_bspline(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        std::vector<NumberType> observation_points,
        std::vector<NumberType> weights,
        size_t number_of_segments) {
    return fit_curve(observations, observation_points, weights, number_of_segments,
            CubicBSplineKernel<NumberType> { });
}

template<int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, CubicCatmulRomKernel<NumberType>>
fit_catmul_rom_spline(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        std::vector<NumberType> observation_points,
        std::vector<NumberType> weights,
        size_t number_of_segments) {
    return fit_curve(observations, observation_points, weights, number_of_segments,
            CubicCatmulRomKernel<NumberType> { });
}

}
}

#endif /* SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_ */
