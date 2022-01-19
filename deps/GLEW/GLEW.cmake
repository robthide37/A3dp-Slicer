# We have to check for OpenGL to compile GLEW
set(OpenGL_GL_PREFERENCE "LEGACY") # to prevent a nasty warning by cmake
find_package(OpenGL QUIET REQUIRED)

prusaslicer_add_cmake_project(
  GLEW
  URL https://sourceforge.net/projects/glew/files/glew/2.1.0/glew-2.1.0.zip
  URL_HASH SHA256=2700383d4de2455f06114fbaf872684f15529d4bdc5cdea69b5fb0e9aa7763f1
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS
    -DBUILD_UTILS=OFF
)

if (MSVC)
    add_debug_dep(dep_GLEW)
endif ()
