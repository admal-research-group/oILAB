configure_file(sigma123.txt sigma123.txt)
configure_file(sigma5.txt sigma5.txt)
# The EAM potential the energies are evaluated with.  Shared with tests/testGbMesoState rather
# than duplicated, so both use the same file.
configure_file(${PROJECT_SOURCE_DIR}/tests/testGbMesoState/Cu_mishin1.eam.alloy
               Cu_mishin1.eam.alloy COPYONLY)
add_executable(GbFacetMeshing GbFacetMeshing.cpp)

target_include_directories(GbFacetMeshing PRIVATE
        ${PROJECT_SOURCE_DIR}/include/IO
        ${PROJECT_SOURCE_DIR}/include/Lattices
)

target_link_libraries(GbFacetMeshing PRIVATE oILAB OpenMP::OpenMP_CXX)

# ---------------------------------------------------------------------------------------------
# tcmalloc, and it is not a micro-optimisation: it is what stops a long sweep decaying.
#
# A faceted sweep under glibc malloc gets slower the longer it runs, and badly.  Measured on
# sigma5 (310), axisScaling 1, slabHalfThickness 0.40, asymmetric displacements, tethered
# relaxation only, over 89 matched progress intervals of a 1,045,533-state run:
#
#                        first interval      last interval
#     LAMMPS   glibc        35.7                177.1      core-ms/state   +396%
#              tcmalloc     36.5                 37.8                        +4%
#     construct glibc        2.06                 4.22
#              tcmalloc      2.02                 2.16
#     coinc    glibc         0.09                 0.18
#              tcmalloc      0.09                 0.09
#
# Same floor at the start -- tcmalloc costs nothing to adopt -- and 4.7x apart by the end.  The
# effect grows with the length of the run, so the sweeps that matter most suffer most: the
# 10.7M-state production sweep was heading for 35 hours with its rate still falling.
#
# It is the ALLOCATOR, not anything in this code.  The phases that decay are the allocation-heavy
# ones, and countCoincidences decays too although it never touches LAMMPS, which is how we know
# no LAMMPS-side explanation could be right.  Ruled out first, each by measurement: memo
# saturation, memo hit rate, memo capacity, shard lock contention, the omp critical section,
# minimiser iterations (flat at ~9), LAMMPS-internal memory (constant at 4.236 MB all run) and
# the atom-ID space (regenerated 1..N per state).  See performanceVsOldFormulation.txt.
#
# Results are unaffected: 150,908 states byte-identical across the allocator change.
#
# Optional, because a machine without it should still build -- but a sweep run without it will
# be several times slower than it needs to be, and nothing in the output will say so.
# The versioned sonames are listed too: the bare libtcmalloc_minimal.so symlink lives in
# gperftools-devel, which a machine that merely RUNS the sweep has no reason to have installed,
# and linking straight against the .so.4 works because the soname is what ends up in the binary
# either way.
find_library(TCMALLOC_LIBRARY
             NAMES tcmalloc_minimal tcmalloc
                   libtcmalloc_minimal.so.4 libtcmalloc.so.4)
if (TCMALLOC_LIBRARY)
    target_link_libraries(GbFacetMeshing PRIVATE ${TCMALLOC_LIBRARY})
    message(STATUS "GbFacetMeshing: linking ${TCMALLOC_LIBRARY} (see the note in this file --
    without it a long sweep slows by a factor of several as it runs)")
else()
    message(WARNING "GbFacetMeshing: tcmalloc not found; falling back to the system allocator.
    Long sweeps will decay -- LAMMPS cost per state rose 396% over one measured run under glibc
    against 4% under tcmalloc.  Install gperftools, or launch with
    LD_PRELOAD=/path/to/libtcmalloc_minimal.so")
endif()

# Sweep output goes to <project>/runs, not to the build tree.  Baked in as an absolute path so
# that the results land in the same place however the example is started -- from CLion, from the
# build directory, or from anywhere else -- and so that clearing or regenerating the build tree
# cannot reach them.
target_compile_definitions(GbFacetMeshing PRIVATE
        GBFACETMESHING_RUNS_DIR="${PROJECT_SOURCE_DIR}/runs")
