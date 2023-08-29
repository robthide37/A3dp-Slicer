set(LibBGCode_SOURCE_DIR "" CACHE PATH "Optionally specify local LibBGCode source directory")

set(_source_dir_line
    URL https://github.com/prusa3d/libbgcode/archive/49e7b47b7a9a7e9365613d3f90ad49d71a0e93d6.zip
	URL_HASH SHA256=b9513ab2bbaf06a79c35c4d69226d211f1226ca625ce07c04f54b89b93e9bc75
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
        -DDEP_CMAKE_OPTS:STRING="-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
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