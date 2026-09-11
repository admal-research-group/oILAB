# add the executable
add_executable(testGbShiftsNonFlat testGbShiftsNonFlat.cpp)
target_link_libraries(testGbShiftsNonFlat oILAB OpenMP::OpenMP_CXX)

# the default (Full) run; edit the USER CHOICE block in the source, or pass arguments, to vary it
add_test(TestGbShiftsNonFlat testGbShiftsNonFlat)
add_test(TestGbShiftsFlatPath testGbShiftsNonFlat flat 1 1)
