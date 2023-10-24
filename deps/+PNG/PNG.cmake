if (APPLE)
    # Only disable NEON extension for Apple ARM builds, leave it enabled for Raspberry PI.
    set(_disable_neon_extension "-DPNG_ARM_NEON:STRING=off")
else ()
    set(_disable_neon_extension "")
endif ()

add_cmake_project(PNG 
    URL https://github.com/glennrp/libpng/archive/refs/tags/v1.6.40.zip
    URL_HASH SHA256=ab3f88779f0661bbb07c60e778fda782216bff70355d86848fbf6a327084563a
    CMAKE_ARGS
        -DPNG_SHARED=OFF
        -DPNG_STATIC=ON
        -DPNG_PREFIX=prusaslicer_
        -DPNG_TESTS=OFF
        -DPNG_EXECUTABLES=OFF
        ${_disable_neon_extension}
)

set(DEP_PNG_DEPENDS ZLIB)
