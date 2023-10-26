if (APPLE)
    # Only disable NEON extension for Apple ARM builds, leave it enabled for Raspberry PI.
    set(_disable_neon_extension "-DPNG_ARM_NEON:STRING=off")
else ()
    set(_disable_neon_extension "")
endif ()

set(_patch_cmd PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/genout.cmake.in scripts/genout.cmake.in)

if (APPLE)
    set(_patch_cmd ${_patch_cmd} && ${PATCH_CMD} ${CMAKE_CURRENT_LIST_DIR}/PNG.patch)
endif ()

add_cmake_project(PNG 
    URL https://github.com/glennrp/libpng/archive/refs/tags/v1.6.40.zip
    URL_HASH SHA256=ab3f88779f0661bbb07c60e778fda782216bff70355d86848fbf6a327084563a
    #PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/genout.cmake.in scripts/genout.cmake.in 
    PATCH_COMMAND "${_patch_cmd}"
    # SOURCE_DIR /home/quarky/Workspace/prusa3d/PrusaSlicer/prusaslicer-src-master/deps/build-default/dep_PNG-prefix/src/dep_PNG/
    CMAKE_ARGS
        -DPNG_SHARED=OFF
        -DPNG_STATIC=ON
        -DPNG_PREFIX=prusaslicer_
        -DPNG_TESTS=OFF
        -DPNG_EXECUTABLES=OFF
        ${_disable_neon_extension}
)

set(DEP_PNG_DEPENDS ZLIB)
