#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/Geometry/Curves.hpp>

TEST_CASE("Curves: constant kernel fitting", "[Curves]") {
    using namespace Slic3r;
    using namespace Slic3r::Geometry;

    std::vector<Vec<1, float>> points { Vec<1, float> { 0 }, Vec<1, float> { 2 }, Vec<1, float> { 7 } };
    size_t number_of_segments = 2;
    float normalized_kernel_bandwidth = 1.0f;
    auto kernel = CurveSmoothingKernels::ConstantKernel<float> { };
    auto curve = fit_curve(points, number_of_segments, normalized_kernel_bandwidth, kernel);

    REQUIRE(curve.length == Approx(7.0f));
    REQUIRE(curve.coefficients[0].size() == number_of_segments);
    REQUIRE(curve.coefficients[0][0] == Approx(1.0f));
    REQUIRE(curve.coefficients[0][1] == Approx(7.0f));

    REQUIRE(curve.get_fitted_point(0.33)[0] == Approx(1.0f));
}

TEST_CASE("Curves: constant kernel fitting 2", "[Curves]") {
    using namespace Slic3r;
    using namespace Slic3r::Geometry;

    std::vector<Vec<1, float>> points { Vec<1, float> { 0 }, Vec<1, float> { 2 }, Vec<1, float> { 2 },
            Vec<1, float> { 4 } };
    size_t number_of_segments = 2;
    float normalized_kernel_bandwidth = 2.0f;
    auto kernel = CurveSmoothingKernels::ConstantKernel<float> { };
    auto curve = fit_curve(points, number_of_segments, normalized_kernel_bandwidth, kernel);

    REQUIRE(curve.length == Approx(4.0f));
    REQUIRE(curve.coefficients[0].size() == number_of_segments);
    REQUIRE(curve.get_fitted_point(0.33)[0] == Approx(2.0f));
}

TEST_CASE("Curves: 2D constant kernel fitting", "[Curves]") {
    using namespace Slic3r;
    using namespace Slic3r::Geometry;

    std::vector<Vec2f> points { Vec2f { 0, 0 }, Vec2f { 2, 1 }, Vec2f { 4, 2 }, Vec2f { 6, 3 } };
    size_t number_of_segments = 4;
    float normalized_kernel_bandwidth = 0.49f;
    auto kernel = CurveSmoothingKernels::ConstantKernel<float> { };
    auto curve = fit_curve(points, number_of_segments, normalized_kernel_bandwidth, kernel);

    REQUIRE(curve.length == Approx(sqrt(6 * 6 + 3 * 3)));
    REQUIRE(curve.coefficients.size() == 2);
    REQUIRE(curve.coefficients[0].size() == number_of_segments);
    REQUIRE(curve.coefficients[0][0] == Approx(0.0f));
    REQUIRE(curve.coefficients[0][1] == Approx(2.0f));
    REQUIRE(curve.coefficients[0][2] == Approx(4.0f));
    REQUIRE(curve.coefficients[0][3] == Approx(6.0f));

    REQUIRE(curve.coefficients[1][0] == Approx(0.0f));
    REQUIRE(curve.coefficients[1][1] == Approx(1.0f));
    REQUIRE(curve.coefficients[1][2] == Approx(2.0f));
    REQUIRE(curve.coefficients[1][3] == Approx(3.0f));
}
