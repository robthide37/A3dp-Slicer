#include <catch2/catch.hpp>

#include <libslic3r/Emboss.hpp>
#include <libslic3r/SVG.hpp> // only debug visualization

#include <optional>
#include <libslic3r/AABBTreeIndirect.hpp>

using namespace Slic3r;

namespace Private{
        
// calculate multiplication of ray dir to intersect - inspired by
// segment_segment_intersection when ray dir is normalized retur distance from
// ray point to intersection No value mean no intersection
std::optional<double> ray_segment_intersection(const Vec2d &r_point,
                                               const Vec2d &r_dir,
                                               const Vec2d &s0,
                                               const Vec2d &s1)
{
    auto denominate = [](const Vec2d &v0, const Vec2d &v1) -> double {
        return v0.x() * v1.y() - v1.x() * v0.y();
    };

    Vec2d  segment_dir = s1 - s0;
    double d           = denominate(segment_dir, r_dir);
    if (std::abs(d) < std::numeric_limits<double>::epsilon())
        // Line and ray are collinear.
        return {};

    Vec2d  s12         = s0 - r_point;
    double s_number    = denominate(r_dir, s12);
    bool   change_sign = false;
    if (d < 0.) {
        change_sign = true;
        d           = -d;
        s_number    = -s_number;
    }

    if (s_number < 0. || s_number > d)
        // Intersection outside of segment.
        return {};

    double r_number = denominate(segment_dir, s12);
    if (change_sign) r_number = -r_number;

    if (r_number < 0.)
        // Intersection before ray start.
        return {};

    return r_number / d;
}

Vec2d get_intersection(const Vec2d &               point,
                       const Vec2d &               dir,
                       const std::array<Vec2d, 3> &triangle)
{
    std::optional<double> t;
    for (size_t i = 0; i < 3; ++i) {
        size_t i2 = i + 1;
        if (i2 == 3) i2 = 0;
        if (!t.has_value()) {
            t = ray_segment_intersection(point, dir, triangle[i],
                                         triangle[i2]);
            continue;
        }

        // small distance could be preccission inconsistance
        std::optional<double> t2 = ray_segment_intersection(point, dir,
                                                            triangle[i],
                                                            triangle[i2]);
        if (t2.has_value() && *t2 > *t) t = t2;
    }
    assert(t.has_value()); // Not found intersection.
    return point + dir * (*t);
}

Vec3d calc_hit_point(const igl::Hit &          h,
                     const Vec3i &             triangle,
                     const std::vector<Vec3f> &vertices)
{
    double c1 = h.u;
    double c2 = h.v;
    double c0 = 1.0 - c1 - c2;
    Vec3d  v0 = vertices[triangle[0]].cast<double>();
    Vec3d  v1 = vertices[triangle[1]].cast<double>();
    Vec3d  v2 = vertices[triangle[2]].cast<double>();
    return v0 * c0 + v1 * c1 + v2 * c2;
}

Vec3d calc_hit_point(const igl::Hit &h, indexed_triangle_set &its)
{
    return calc_hit_point(h, its.indices[h.id], its.vertices);
}
} // namespace Private

std::string get_font_filepath() {
    std::string resource_dir = 
        std::string(TEST_DATA_DIR) + "/../../resources/";
    return resource_dir + "fonts/NotoSans-Regular.ttf";
}

#include "imgui/imstb_truetype.h"
TEST_CASE("Read glyph C shape from font, stb library calls ONLY", "[Emboss]") {
    std::string font_path = get_font_filepath();
    char  letter   = 'C';
    
    // Read  font file
    FILE *file = fopen(font_path.c_str(), "rb");
    REQUIRE(file != nullptr);
    // find size of file
    REQUIRE(fseek(file, 0L, SEEK_END) == 0);
    size_t size = ftell(file);
    REQUIRE(size != 0);
    rewind(file);
    std::vector<unsigned char> buffer(size);
    size_t count_loaded_bytes = fread((void *) &buffer.front(), 1, size, file);
    REQUIRE(count_loaded_bytes == size);

    // Use stb true type library
    int font_offset = stbtt_GetFontOffsetForIndex(buffer.data(), 0);
    REQUIRE(font_offset >= 0);
    stbtt_fontinfo font_info;
    REQUIRE(stbtt_InitFont(&font_info, buffer.data(), font_offset) != 0);    
    int unicode_letter = (int) letter;
    int glyph_index = stbtt_FindGlyphIndex(&font_info, unicode_letter);
    REQUIRE(glyph_index != 0);
    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&font_info, glyph_index, &vertices);
    CHECK(num_verts > 0);
}

#include <libslic3r/Utils.hpp>
TEST_CASE("Convert glyph % to model", "[Emboss]") 
{
    std::string font_path = get_font_filepath();
    char  letter   = '%';
    float flatness = 2.;

    auto font = Emboss::create_font_file(font_path.c_str());
    REQUIRE(font != nullptr);

    std::optional<Emboss::Glyph> glyph = Emboss::letter2glyph(*font, letter, flatness);
    REQUIRE(glyph.has_value());

    ExPolygons shape = glyph->shape;    
    REQUIRE(!shape.empty());

    float z_depth = 1.f;
    Emboss::ProjectZ projection(z_depth);
    indexed_triangle_set its = Emboss::polygons2model(shape, projection);

    CHECK(!its.indices.empty());    
}

TEST_CASE("Test hit point", "[AABBTreeIndirect]")
{
    indexed_triangle_set its;
    its.vertices = {
        Vec3f(1, 1, 1),
        Vec3f(2, 10, 2),
        Vec3f(10, 0, 2),
    };
    its.indices = {Vec3i(0, 2, 1)};
    auto tree   = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
        its.vertices, its.indices);

    Vec3d    ray_point(8, 1, 0);
    Vec3d    ray_dir(0, 0, 1);
    igl::Hit hit;
    AABBTreeIndirect::intersect_ray_first_hit(its.vertices, its.indices, tree,
                                              ray_point, ray_dir, hit);
    Vec3d hp = Private::calc_hit_point(hit, its);
    CHECK(abs(hp.x() - ray_point.x()) < .1);
    CHECK(abs(hp.y() - ray_point.y()) < .1);
}

TEST_CASE("ray segment intersection", "[MeshBoolean]")
{
    Vec2d r_point(1, 1);
    Vec2d r_dir(1, 0);

    // colinear
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 0), Vec2d(2, 0)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(0, 0)).has_value());

    // before ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 0), Vec2d(0, 2)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(0, 2), Vec2d(0, 0)).has_value());

    // above ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 2), Vec2d(2, 3)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 3), Vec2d(2, 2)).has_value());

    // belove ray
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(2, -1)).has_value());
    CHECK(!Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, -1), Vec2d(2, 0)).has_value());

    // intersection at [2,1] distance 1
    auto t1 = Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 0), Vec2d(2, 2));
    REQUIRE(t1.has_value());
    auto t2 = Private::ray_segment_intersection(r_point, r_dir, Vec2d(2, 2), Vec2d(2, 0));
    REQUIRE(t2.has_value());

    CHECK(abs(*t1 - *t2) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("triangle intersection", "[]")
{
    Vec2d                point(1, 1);
    Vec2d                dir(-1, 0);
    std::array<Vec2d, 3> triangle = {Vec2d(0, 0), Vec2d(5, 0), Vec2d(0, 5)};
    Vec2d                i = Private::get_intersection(point, dir, triangle);
    CHECK(abs(i.x()) < std::numeric_limits<double>::epsilon());
    CHECK(abs(i.y() - 1.) < std::numeric_limits<double>::epsilon());
}

#ifndef __APPLE__
#include <string>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
// Check function Emboss::is_italic that exist some italic and some non-italic font.
TEST_CASE("Italic check", "[Emboss]") 
{  
    std::queue<std::string> dir_paths;
#ifdef _WIN32
    dir_paths.push("C:/Windows/Fonts");
#elif defined(__linux__)
    dir_paths.push("/usr/share/fonts");
//#elif defined(__APPLE__)
//    dir_paths.push("//System/Library/Fonts");
#endif
    bool exist_italic = false;
    bool exist_non_italic = false;
    while (!dir_paths.empty()) {
        std::string dir_path = dir_paths.front();
        dir_paths.pop();
        for (const auto &entry : fs::directory_iterator(dir_path)) {
            const fs::path &act_path = entry.path();
            if (entry.is_directory()) {
                dir_paths.push(act_path.u8string());
                continue;
            }
            std::string ext = act_path.extension().u8string();
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (ext != ".ttf") continue;
            std::string path_str = act_path.u8string();
            auto        font_opt = Emboss::create_font_file(path_str.c_str());
            if (font_opt == nullptr) continue;

            unsigned int collection_number = 0;
            if (Emboss::is_italic(*font_opt, collection_number))
                exist_italic = true;
            else
                exist_non_italic = true;

            if (exist_italic && exist_non_italic) break;
            //std::cout << ((Emboss::is_italic(*font_opt)) ? "[yes] " : "[no ] ") << entry.path() << std::endl;
        }
    }
    CHECK(exist_italic);
    CHECK(exist_non_italic);
}
#endif // not __APPLE__
