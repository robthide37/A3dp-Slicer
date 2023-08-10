#/|/ Copyright (c) Prusa Research 2018 - 2022 Lukáš Matěna @lukasmatena, Tomáš Mészáros @tamasmeszaros, Filip Sykala @Jony01, Vojtěch Bubník @bubnikv, Vojtěch Král @vojtechkral
#/|/
#/|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
#/|/
#
# This CMake project downloads, configures and builds PrusaSlicer dependencies on Unix and Windows.
#
# When using this script, it's recommended to perform an out-of-source build using CMake.
#
# All the dependencies are installed in a `destdir` directory in the root of the build directory,
# in a traditional Unix-style prefix structure. The destdir can be used directly by CMake
# when building PrusaSlicer - to do this, set the CMAKE_PREFIX_PATH to ${destdir}/usr/local.
# Warning: On UNIX/Linux, you also need to set -DSLIC3R_STATIC=1 when building PrusaSlicer.
#
# For better clarity of console output, it's recommended to _not_ use a parallelized build
# for the top-level command, ie. use `make -j 1` or `ninja -j 1` to force single-threaded top-level
# build. This doesn't degrade performance as individual dependencies are built in parallel fashion
# if supported by the dependency.
#
# On Windows, architecture (64 vs 32 bits) is judged based on the compiler variant.
# To build dependencies for either 64 or 32 bit OS, use the respective compiler command line.
#
# WARNING: On UNIX platforms wxWidgets hardcode the destdir path into its `wx-conffig` utility,
# therefore, unfortunatelly, the installation cannot be copied/moved elsewhere without re-installing wxWidgets.
#
cmake_minimum_required(VERSION 3.12)
project(PrusaSlicer_deps)

# Always ON options:
if (NOT MSVC)
    set(DEP_CMAKE_OPTS "-DCMAKE_POSITION_INDEPENDENT_CODE=ON")
endif ()

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    cmake_policy(SET CMP0135 NEW)
endif ()

include(${PROJECT_SOURCE_DIR}/../cmake/modules/AddCMakeProject.cmake)

macro(list_projects result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
      string(REGEX MATCH "^\\+([a-zA-Z0-9]+)" is_package ${child})
      if(is_package) 
        list(APPEND dirlist ${CMAKE_MATCH_1})
      endif()
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

if (APPLE)
    # This ensures dependencies don't use SDK features which are not available in the version specified by Deployment target
    # That can happen when one uses a recent SDK but specifies an older Deployment target
    set(DEP_WERRORS_SDK "-Werror=partial-availability -Werror=unguarded-availability -Werror=unguarded-availability-new")
    
    set(DEP_CMAKE_OPTS
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
        "-DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}"
        "-DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}"
        "-DCMAKE_CXX_FLAGS=${DEP_WERRORS_SDK}"
        "-DCMAKE_C_FLAGS=${DEP_WERRORS_SDK}"
        "-DCMAKE_FIND_FRAMEWORK=LAST"
        "-DCMAKE_FIND_APPBUNDLE=LAST"
    )
endif ()

list_projects(SUPPORTED_PACKAGES ${CMAKE_CURRENT_LIST_DIR})
message(STATUS "Found external package definitions: ${SUPPORTED_PACKAGES}")

set(REQUIRED_PACKAGES
    Boost
    Catch2
    Cereal
    CURL
    EXPAT
    wxWidgets
)

set(SYSTEM_PROVIDED_PACKAGES OpenGL)

if (MSVC)
    list(APPEND REQUIRED_PACKAGES ZLIB)
elseif (UNIX)
    # On UNIX systems (including Apple) ZLIB should be available
    list(APPEND SYSTEM_PROVIDED_PACKAGES ZLIB)
endif ()

include(CMakeDependentOption)
option(${PROJECT_NAME}_SELECT_ALL "Choose all external projects to be built." ON)

find_package(Git REQUIRED)

# The default command line for patching. Only works for newer 
set(PATCH_CMD ${GIT_EXECUTABLE} apply --verbose --ignore-space-change --whitespace=fix)

# all required package targets that have existing definitions will be gathered here
set(_dep_list "")
set(_build_list "")

foreach (pkg ${SUPPORTED_PACKAGES})
    cmake_dependent_option(${PROJECT_NAME}_SELECT_${pkg} "Select package ${pkg} to be built." OFF  "NOT ${PROJECT_NAME}_SELECT_ALL" OFF)
    if (${PROJECT_NAME}_SELECT_ALL OR ${PROJECT_NAME}_SELECT_${pkg})
        include(+${pkg}/${pkg}.cmake)

        if (${pkg} IN_LIST REQUIRED_PACKAGES)
            list(APPEND _dep_list dep_${pkg})
        endif ()
    endif ()
endforeach()

# Establish dependency graph
foreach (pkg ${REQUIRED_PACKAGES})
  foreach(deppkg ${DEP_${pkg}_DEPENDS})
    if (${deppkg} IN_LIST SYSTEM_PROVIDED_PACKAGES)
      find_package(${deppkg} QUIET)
      if (NOT ${deppkg}_FOUND)
        message(WARNING "No ${deppkg} found in system altough marked as system provided. This might cause trouble building the dependencies on this platform")
      endif ()
    else ()
      message(STATUS "Mapping dep_${deppkg} => dep_${pkg}")
      add_dependencies(dep_${pkg} dep_${deppkg})
      list(APPEND _build_list dep_${deppkg})
    endif ()
  endforeach()
endforeach()

list(APPEND _build_list ${_dep_list})
list(REMOVE_DUPLICATES _build_list)
message(STATUS "Building dep targets: ${_build_list}")
add_custom_target(deps ALL DEPENDS ${_dep_list})

install(CODE "message(STATUS \"Built packages succesfully.\")")