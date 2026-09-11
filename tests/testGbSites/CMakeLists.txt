# add the executable
add_executable(testGbSites testGbSites.cpp)
target_link_libraries(testGbSites oILAB OpenMP::OpenMP_CXX)

# arguments: <slab half-thickness in A> <dMax in A>; the source carries the defaults
add_test(TestGbSites testGbSites)
