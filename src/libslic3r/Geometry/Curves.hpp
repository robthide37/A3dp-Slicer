#ifndef SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_
#define SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_

#include "libslic3r/Point.hpp"
#include "Bicubic.hpp"

#include <iostream>

namespace Slic3r {
namespace Geometry {

template<int Dimension, typename NumberType>
struct PolynomialCurve {
    Eigen::MatrixXf coefficients;

    Vec3f get_fitted_value(const NumberType value) const {
        auto result = Vec<Dimension, NumberType>::Zero();
        size_t order = this->coefficients.rows() - 1;
        auto x = NumberType(1.);
        for (size_t index = 0; index < order + 1; ++index, x *= value)
            result += x * this->coefficients.col(index);
        return result;
    }
};

//https://towardsdatascience.com/least-square-polynomial-CURVES-using-c-eigen-package-c0673728bd01
template<int Dimension, typename NumberType>
PolynomialCurve<Dimension, NumberType> fit_polynomial(const std::vector<Vec<Dimension, NumberType>> &observations,
        const std::vector<NumberType> &observation_points,
        const std::vector<NumberType> &weights, size_t order) {
    // check to make sure inputs are correct
    size_t cols = order + 1;
    assert(observation_points.size() >= cols);
    assert(observation_points.size() == weights.size());
    assert(observations.size() == weights.size());

    Eigen::MatrixXf data_points(Dimension, observations.size());
    Eigen::MatrixXf T(observations.size(), cols);
    for (size_t i = 0; i < weights.size(); ++i) {
        auto squared_weight = sqrt(weights[i]);
        data_points.col(i) = observations[i] * squared_weight;
        // Populate the matrix
        auto x = squared_weight;
        auto c = observation_points[i];
        for (size_t j = 0; j < cols; ++j, x *= c)
            T(i, j) = x;
    }

    const auto QR = T.householderQr();
    Eigen::MatrixXf coefficients(Dimension, cols);
    // Solve for linear least square fit
    for (size_t dim = 0; dim < Dimension; ++dim) {
        coefficients.row(dim) = QR.solve(data_points.row(dim).transpose());
    }

    return { std::move(coefficients) };
}

template<size_t Dimension, typename NumberType, typename KernelType>
struct PiecewiseFittedCurve {
    using Kernel = KernelType;

    Eigen::MatrixXf coefficients;
    NumberType start;
    NumberType length;
    NumberType n_segment_size;

    size_t     segments() const { return this->coefficients.cols(); }

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
            if (segment_index < 0 || segment_index >= int(this->coefficients.cols())) {
                continue;
            }
            NumberType segment_start = this->get_n_segment_start(segment_index);
            NumberType normalized_segment_distance = (segment_start - t) / this->n_segment_size;

            result += Kernel::kernel(normalized_segment_distance) * coefficients.col(segment_index);
        }
        return result;
    }
};

// observations: data to be fitted by the curve
// observation points: growing sequence of points where the observations were made.
//      In other words, for function f(x) = y, observations are y0...yn, and observation points are x0...xn
// weights: how important the observation is
//  number_of_inner_splines: how many full inner splines are fit into the normalized valid range 0,1;
//          final number of knots is Kernel::kernel_span times number_of_inner_splines + additional segments on start and end
// Kernel: model used for the curve fitting
template<typename Kernel, int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, Kernel> fit_curve(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        const std::vector<NumberType> &observation_points,
        const std::vector<NumberType> &weights,
        size_t number_of_inner_splines) {

    // check to make sure inputs are correct
    assert(number_of_inner_splines >= 1);
    assert(observation_points.size() == observations.size());
    assert(observation_points.size() == weights.size());
    assert(number_of_inner_splines <= observations.size());
    assert(observations.size() >= 4);

    size_t extremes_repetition = Kernel::kernel_span - 1; //how many (additional) times is the first and last point repeated

    //prepare sqrt of weights, which will then be applied to both matrix T and observed data: https://en.wikipedia.org/wiki/Weighted_least_squares
    std::vector<NumberType> sqrt_weights(weights.size() + extremes_repetition * 2);
    for (size_t index = 0; index < weights.size(); ++index) {
        assert(weights[index] > 0);
        sqrt_weights[index + extremes_repetition] = sqrt(weights[index]);
    }
    //repeat weights for addtional extreme segments
    {
        auto first = sqrt_weights[extremes_repetition];
        auto last  = sqrt_weights[extremes_repetition + weights.size() - 1];
        for (int index = 0; index < int(extremes_repetition); ++index) {
            sqrt_weights[index] = first;
            sqrt_weights[sqrt_weights.size() - index - 1] = last;
        }
    }

    // prepare result and compute metadata
    PiecewiseFittedCurve<Dimension, NumberType, Kernel> result { };

    NumberType orig_len = observation_points.back() - observation_points.front();
    NumberType orig_segment_size = orig_len / NumberType(number_of_inner_splines * Kernel::kernel_span);
    result.start = observation_points.front() - extremes_repetition * orig_segment_size;
    result.length = observation_points.back() + extremes_repetition * orig_segment_size - result.start;
    size_t segments_count = number_of_inner_splines * Kernel::kernel_span + extremes_repetition * 2;
    result.n_segment_size = NumberType(1) / NumberType(segments_count - 1);

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
    // Eigen defaults to column major memory layout.
    Eigen::MatrixXf data_points(Dimension, observations.size() + extremes_repetition * 2);
    for (size_t index = 0; index < observations.size(); ++ index) {
        data_points.col(index + extremes_repetition) = observations[index]
                * sqrt_weights[index + extremes_repetition];
    }
    //duplicate observed data at the extremes
    for (int index = 0; index < int(extremes_repetition); index++) {
        data_points.col(index) = observations.front() * sqrt_weights[index];
        data_points.col(data_points.cols() - index - 1) = observations.back()
                * sqrt_weights[data_points.cols() - index - 1];
    }

    //Create weight matrix T for each point and each segment;
    Eigen::MatrixXf T(normalized_obs_points.size(), segments_count);
    T.setZero();

    //Fill the weight matrix
    for (size_t i = 0; i < normalized_obs_points.size(); ++i) {
        NumberType knot_val = normalized_obs_points[i];
        //find index of first segment that is affected by the point i; this can be deduced from kernel_span
        int start_segment_idx = int(floor(knot_val / result.n_segment_size)) - int(Kernel::kernel_span * 0.5f - 1.0f);
        for (int segment_index = start_segment_idx; segment_index < int(start_segment_idx + Kernel::kernel_span);
                segment_index++) {
            // skip if we overshoot segment_index - happens at the extremes
            if (segment_index < 0 || segment_index >= int(segments_count)) {
                continue;
            }
            NumberType segment_start = result.get_n_segment_start(segment_index);
            NumberType normalized_segment_distance = (segment_start - knot_val) / result.n_segment_size;
            // fill in kernel value with weight applied

            T(i, segment_index) += Kernel::kernel(normalized_segment_distance) * sqrt_weights[i];
        }
    }

    // Solve for linear least square fit
    result.coefficients.resize(Dimension, segments_count);
    const auto QR = T.fullPivHouseholderQr();
    for (size_t dim = 0; dim < Dimension; ++dim) {
         result.coefficients.row(dim) = QR.solve(data_points.row(dim).transpose());
    }

    return result;
}

template<int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, CubicBSplineKernel<NumberType>>
fit_cubic_bspline(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        std::vector<NumberType> observation_points,
        std::vector<NumberType> weights,
        size_t number_of_segments) {
    return fit_curve<CubicBSplineKernel<NumberType>>(observations, observation_points, weights, number_of_segments);
}

template<int Dimension, typename NumberType>
PiecewiseFittedCurve<Dimension, NumberType, CubicCatmulRomKernel<NumberType>>
fit_catmul_rom_spline(
        const std::vector<Vec<Dimension, NumberType>> &observations,
        std::vector<NumberType> observation_points,
        std::vector<NumberType> weights,
        size_t number_of_segments) {
    return fit_curveCubicCatmulRomKernel<NumberType>(observations, observation_points, weights, number_of_segments);
}

}
}

#endif /* SRC_LIBSLIC3R_GEOMETRY_CURVES_HPP_ */
