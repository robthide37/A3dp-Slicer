#ifndef SLIC3R_TEST_UTILS
#define SLIC3R_TEST_UTILS

#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/Format/OBJ.hpp>
#include <random>

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR R"(\)"
#else
#define PATH_SEPARATOR R"(/)"
#endif

inline Slic3r::TriangleMesh load_model(const std::string &obj_filename)
{
    Slic3r::TriangleMesh mesh;
    auto fpath = TEST_DATA_DIR PATH_SEPARATOR + obj_filename;
    Slic3r::load_obj(fpath.c_str(), &mesh);
    return mesh;
}

inline std::string get_model_path(const std::string &obj_filename)
{
    std::string fpath = TEST_DATA_DIR PATH_SEPARATOR + obj_filename;
    return fpath;
}
template<class T>
Slic3r::FloatingOnly<T> random_value(T minv, T maxv)
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<T> dist(minv, maxv);

    return dist(rng);
}

template<class T>
Slic3r::IntegerOnly<T> random_value(T minv, T maxv)
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<T> dist(minv, maxv);

    return dist(rng);
}

#endif // SLIC3R_TEST_UTILS
