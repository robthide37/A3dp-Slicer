set(LibBGCode_SOURCE_DIR "" CACHE PATH "Optionally specify local LibBGCode source directory")

set(_source_dir_line
    URL https://github.com/prusa3d/libbgcode/archive/6b4e4f9fd1270f58f9202a674c7f8e7e00d9cfb0.zip
	URL_HASH SHA256=1c613f8340657f87a3179711cb14de346c8ba3b160ac97b7b0d8bfa81d00b98a
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
        -DDEP_CMAKE_OPTS:STRING=-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
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