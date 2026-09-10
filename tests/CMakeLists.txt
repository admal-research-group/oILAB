enable_testing()
add_subdirectory(testSNF)
add_subdirectory(testSNF_cubic)
add_subdirectory(testPerpendicularDirections)
add_subdirectory(testPlaneParallelLatticeDirections)
add_subdirectory(testBezout)
add_subdirectory(testCCUM)
add_subdirectory(testCoincidentRotations)
add_subdirectory(testGb)
add_subdirectory(testGB3d)
add_subdirectory(testLattice)
add_subdirectory(testGenerateGBs)
add_subdirectory(testGbMesoState)
add_subdirectory(testGbShifts)
add_subdirectory(testGbFacet)
add_subdirectory(testGbShiftsNonFlat)
add_subdirectory(testGbSites)
add_subdirectory(testMoire)

# The three Monte Carlo tests are out of the build until they are brought up to the current API.
# They do not compile, and have not for some time: each sets GbMaterialTensors::lambda and ::mu
# and calls GbContinuum<3>::reset(), and neither the class nor the function exists any more.
# Their LAMMPS plumbing is current -- they no longer name an `lmp` binary -- so what is left is
# to replace those three calls with whatever now carries the elastic constants.
#add_subdirectory(testMonteCarlo)
#add_subdirectory(testLandauWang)
#add_subdirectory(testLandauWang2D)
