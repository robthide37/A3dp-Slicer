set(LibBGCode_SOURCE_DIR "" CACHE PATH "Optionally specify local LibBGCode source directory")

set(_source_dir_line "URL;https://github.com/prusa3d/libbgcode/archive/refs/heads/main.zip")

if (LibBGCode_SOURCE_DIR)
    set(_source_dir_line "SOURCE_DIR;${LibBGCode_SOURCE_DIR};BUILD_ALWAYS;ON")
endif ()

prusaslicer_add_cmake_project(LibBGCode_deps
    ${_source_dir_line}
    SOURCE_SUBDIR deps
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
)