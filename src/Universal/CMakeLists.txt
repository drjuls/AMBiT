set(MODS_UNIVERSAL  Eigensolver.cpp
                    ExpLattice.cpp
                    FornbergDifferentiator.cpp
                    Include.cpp
                    Interpolator.cpp
                    Lattice.cpp
                    MathConstant.cpp
                    PhysicalConstant.cpp
                    ScalapackMatrix.cpp
                    SpinorFunction.cpp
                    Davidson.f
                    CACHE INTERNAL "")

add_library(universal "${MODS_UNIVERSAL}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(universal PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(universal PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(universal PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(universal PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(universal PUBLIC Eigen3::Eigen)
target_link_libraries(universal PUBLIC GSL::gsl)
target_include_directories(universal PRIVATE ${Boost_INCLUDE_DIRS})
