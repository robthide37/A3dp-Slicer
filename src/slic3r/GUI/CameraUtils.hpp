#ifndef slic3r_CameraUtils_hpp_
#define slic3r_CameraUtils_hpp_

#include "Camera.hpp"
#include "libslic3r/Point.hpp"

namespace Slic3r::GUI {

/// <summary>
/// Help divide camera data and camera functions
/// This utility work with camera data by static funtions
/// </summary>
class CameraUtils
{
public:
    CameraUtils() = delete; // only static functions

	/// <summary>
	/// project point throw camera to 2d coordinate into imgui window
	/// </summary>
	/// <param name="its">IN/OUT triangle mesh to be simplified.</param>
	/// <param name="its">IN/OUT triangle mesh to be simplified.</param>
    /// <returns>projected points by camera into coordinate of camera.
    /// x(from left to right), y(from top to bottom)</returns>
	static Points project(const Camera& camera, const std::vector<Vec3d> &points);

};
} // Slic3r::GUI

#endif /* slic3r_CameraUtils_hpp_ */
