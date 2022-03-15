#include <Eigen/Dense>
#include <Eigen/QR>
#include "Curves.hpp"


namespace Slic3r {
namespace Geometry {

PolynomialCurve fit_polynomial(const std::vector<Vec3f> &points, const std::vector<float> &weights, size_t order) {
    // check to make sure inputs are correct
    assert(points.size() >= order + 1);
    assert(points.size() == weights.size());

    std::vector<float> squared_weights(weights.size());
    for (size_t index = 0; index < weights.size(); ++index) {
        squared_weights[index] = sqrt(weights[index]);
    }

    Eigen::VectorXf V0(points.size());
    Eigen::VectorXf V1(points.size());
    Eigen::VectorXf V2(points.size());
    for (size_t index = 0; index < points.size(); index++) {
        V0(index) = points[index].x() * squared_weights[index];
        V1(index) = points[index].y() * squared_weights[index];
        V2(index) = points[index].z();
    }

    // Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
    Eigen::MatrixXf T(points.size(), order + 1);
    // Populate the matrix
    for (size_t i = 0; i < points.size(); ++i)
            {
        for (size_t j = 0; j < order + 1; ++j)
                {
            T(i, j) = pow(V2(i), j) * squared_weights[i];
        }
    }

    // Solve for linear least square fit
    const auto QR = T.householderQr();
    Eigen::VectorXf result0 = QR.solve(V0);
    Eigen::VectorXf result1 = QR.solve(V1);
    std::vector<Vec2f> coeff { order + 1 };
    for (size_t k = 0; k < order + 1; k++) {
        coeff[k] = Vec2f { result0[k], result1[k] };
    }

    return PolynomialCurve{coeff};
}

Vec3f PolynomialCurve::get_fitted_point(float z) const {
    size_t order = this->coefficients.size() - 1;
    float fitted_x = 0;
    float fitted_y = 0;
    for (size_t index = 0; index < order + 1; ++index) {
        float z_pow = pow(z, index);
        fitted_x += this->coefficients[index].x() * z_pow;
        fitted_y += this->coefficients[index].y() * z_pow;
    }

    return Vec3f { fitted_x, fitted_y, z };
}

}
}
