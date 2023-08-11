include(GNUInstallDirs)
# WARNING: Qhull does this nice thing to remove previous configs when installing.
# It does that only on Unix. If you build release, then install, then build
# debug, then install, it will delete the previous release config....
add_cmake_project(Qhull
    URL "https://github.com/qhull/qhull/archive/v8.0.1.zip"
    URL_HASH SHA256=5287f5edd6a0372588f5d6640799086a4033d89d19711023ef8229dd9301d69b
    CMAKE_ARGS 
        -DINCLUDE_INSTALL_DIR=${CMAKE_INSTALL_INCLUDEDIR}
)
