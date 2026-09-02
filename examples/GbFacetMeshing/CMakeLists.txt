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
