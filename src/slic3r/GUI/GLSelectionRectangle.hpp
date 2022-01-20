#ifndef slic3r_GLSelectionRectangle_hpp_
#define slic3r_GLSelectionRectangle_hpp_

#include "libslic3r/Point.hpp"
#if ENABLE_GLBEGIN_GLEND_REMOVAL
#include "GLModel.hpp"
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

namespace Slic3r {
namespace GUI {
    
struct Camera;
class GLCanvas3D;

class GLSelectionRectangle {
public:
    enum EState {
            Off,
            Select,
            Deselect
    };

    // Initiates the rectangle. 
    void start_dragging(const Vec2d& mouse_position, EState state);

    // To be called on mouse move.
    void dragging(const Vec2d& mouse_position);

    // Given a vector of points in world coordinates, the function returns indices of those
    // that are in the rectangle. It then disables the rectangle.
    std::vector<unsigned int> stop_dragging(const GLCanvas3D& canvas, const std::vector<Vec3d>& points);

    // Disables the rectangle.
    void stop_dragging();

#if ENABLE_GLBEGIN_GLEND_REMOVAL
    void render(const GLCanvas3D& canvas);
#else
    void render(const GLCanvas3D& canvas) const;
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL

    bool is_dragging() const { return m_state != Off; }
    EState get_state() const { return m_state; }

    float get_width() const { return std::abs(m_start_corner(0) - m_end_corner(0)); }
    float get_height() const { return std::abs(m_start_corner(1) - m_end_corner(1)); }
    float get_left() const { return std::min(m_start_corner(0), m_end_corner(0)); }
    float get_right() const { return std::max(m_start_corner(0), m_end_corner(0)); }
    float get_top() const { return std::max(m_start_corner(1), m_end_corner(1)); }
    float get_bottom() const { return std::min(m_start_corner(1), m_end_corner(1)); }

private:
    EState m_state{ Off };
    Vec2d m_start_corner{ Vec2d::Zero() };
    Vec2d m_end_corner{ Vec2d::Zero() };
#if ENABLE_GLBEGIN_GLEND_REMOVAL
    GLModel m_rectangle;
    Vec2d m_old_start_corner{ Vec2d::Zero() };
    Vec2d m_old_end_corner{ Vec2d::Zero() };
#endif // ENABLE_GLBEGIN_GLEND_REMOVAL
};

    
} // namespace GUI
} // namespace Slic3r


#endif // slic3r_GLGizmoSlaSupports_hpp_
