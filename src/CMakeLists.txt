add_library(oILAB SHARED    Lattices/BiCrystal.cpp
                            Lattices/Gb.cpp
                            Lattices/Lattice.cpp 
                            Lattices/LatticeVector.cpp
                            Math/RationalMatrix.cpp
                            Lattices/LatticeCore.cpp
                            Math/RLLL.cpp
                            Lattices/ReciprocalLatticeDirection.cpp
                            Math/LLL.cpp
                            Lattices/LatticeDirection.cpp
                            Lattices/RationalLatticeDirection.cpp
                            Lattices/RationalReciprocalLatticeDirection.cpp
                            Lattices/ReciprocalLatticeVector.cpp
                            Math/Farey.cpp
                            Lattices/GbShifts.cpp
                            Lattices/GbFacet.cpp)


# Conditionally apply the export property for MSVC on Windows
if(MSVC)
    set_property(TARGET oILAB PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_include_directories(oILAB
    PUBLIC
        # For build tree: point to include/ and its subdirectories for current includes to work
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        # For install tree: point to include/ so users can do #include <oILAB/Lattices/LatticeModule.h>
        $<INSTALL_INTERFACE:include>
)

# GbFacet.h is a public header that includes CGAL, so CGAL must be PUBLIC.
#
# OpenMP is PUBLIC for a reason beyond GbMesoState's use of omp_get_thread_num(): CGAL decides
# whether it may use unguarded static variables from whether _OPENMP (or BOOST_HAS_THREADS) is
# defined when its headers are parsed.  Compiling this library without OpenMP while the callers
# have it left CGAL's triangulation believing it was single threaded, and calling it from many
# threads then crashed sporadically inside the periodic Delaunay code.  Both sides must agree.
target_link_libraries(oILAB PUBLIC Eigen3::Eigen CGAL::CGAL OpenMP::OpenMP_CXX)

# Install library
install(TARGETS oILAB
    EXPORT oILABTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

