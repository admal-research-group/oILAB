add_subdirectory(bicrystal)
add_subdirectory(GbFacetMeshing)

# Copy python examples to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/python
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
