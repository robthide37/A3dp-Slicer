set(LibBGCode_SOURCE_DIR "" CACHE PATH "Optionally specify local LibBGCode source directory")

set(_source_dir_line
    URL https://github.com/prusa3d/libbgcode/archive/767a0970a0df277c366a08f32b72e39b6e628b54.zip
    URL_HASH SHA256=60025397e9ce301ac925b59a12cfced26034b6124866b618c4771842d109847f
)

if (LibBGCode_SOURCE_DIR)
    set(_source_dir_line "SOURCE_DIR;${LibBGCode_SOURCE_DIR};BUILD_ALWAYS;ON")
endif ()

prusaslicer_add_cmake_project(LibBGCode_deps
    ${_source_dir_line}
    SOURCE_SUBDIR deps
    DEPENDS dep_Boost ${ZLIB_PKG}
    CMAKE_ARGS
        -DDEP_DOWNLOAD_DIR:PATH=${DEP_DOWNLOAD_DIR}
        -DLibBGCode_Deps_SELECT_ALL:BOOL=OFF
        -DLibBGCode_Deps_SELECT_heatshrink:BOOL=ON
        -DDESTDIR=${DESTDIR}
)

prusaslicer_add_cmake_project(LibBGCode
    ${_source_dir_line}
    DEPENDS dep_LibBGCode_deps
    CMAKE_ARGS
        -DLibBGCode_BUILD_TESTS:BOOL=OFF
        -DLibBGCode_BUILD_CMD_TOOL:BOOL=OFF
)

if (MSVC)
    add_debug_dep(dep_LibBGCode)
endif ()