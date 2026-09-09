cmake_minimum_required(VERSION 3.20)
project(oILAB VERSION 1.0 LANGUAGES C CXX)

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------- Build type ----------
# Default to Release, but only when the caller has not asked for something else.  Setting
# CMAKE_BUILD_TYPE unconditionally creates a normal variable that shadows the cache for the rest of
# the configure, so -DCMAKE_BUILD_TYPE=Debug (or an IDE's Debug profile) would be silently ignored
# and every build would come out optimised with NDEBUG defined -- which disables every assert in
# the code base.
# FORCE is required: project() already put an empty CMAKE_BUILD_TYPE in the cache, and a plain
# set(... CACHE ...) will not overwrite an existing entry.  It is safe here because the guard is
# false whenever the caller supplied a build type, so an explicit choice is never overridden.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type: Debug, Release, RelWithDebInfo or MinSizeRel" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)

if(MSVC)
    add_compile_options(/bigobj /MP)
endif()

# ---------- Dependencies ----------
find_package(Eigen3 5.0 REQUIRED NO_MODULE)
find_package(CGAL REQUIRED CONFIG)
find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)
# CGAL_NDEBUG is deliberately NOT set.  It compiles out CGAL's preconditions and internal
# invariant checks, and the periodic Delaunay triangulation relies on those to reject point sets
# it cannot handle -- degenerate or near-degenerate arrangements on the torus, which the mesostate
# enumeration produces in quantity.  With the checks compiled out those inputs corrupted the
# triangulation instead of being rejected, and the run died later inside the covering-space
# conversion.  Keeping them costs a little speed and turns a segmentation fault into a catchable
# error.
# add_compile_definitions(CGAL_NDEBUG)
find_package(OpenMP REQUIRED)

# ---------- Subdirectories ----------
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(userExamples)

# ---------- Testing (optional) ----------
option(ENABLE_TESTING "Build tests" ON)
if(ENABLE_TESTING)
    add_subdirectory(tests)
endif()

 # ---------- Python bindings (optional) ----------
option(BUILD_PYBINDINGS "Build Python bindings (pyOILAB)" ON)
if(BUILD_PYBINDINGS)
    add_subdirectory(bindings)
endif()

# ---------- Installation ----------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install header files (preserves subdirectory structure under oILAB/)
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oILAB
    FILES_MATCHING PATTERN "*.h"
)

# Generate package config files
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/oILABConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install CMake config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/oILABConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)

# Export targets
install(EXPORT oILABTargets
    FILE oILABTargets.cmake
    NAMESPACE oILAB::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oILAB
)
