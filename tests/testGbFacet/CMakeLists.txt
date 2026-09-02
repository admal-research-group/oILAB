# add the executable
add_executable(testGbFacet testGbFacet.cpp)
target_link_libraries(testGbFacet oILAB)

add_test(TestGbFacet testGbFacet)
