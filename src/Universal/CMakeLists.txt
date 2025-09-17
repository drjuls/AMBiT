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

if(USE_SCALAPACK)
  set(MODS_UNIVERSAL ${MODS_UNIVERSAL} "ScalapackMatrix.cpp")
endif()

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

if(USE_SCALAPACK)
  target_compile_definitions(universal PRIVATE -DAMBIT_USE_SCALAPACK)
  if(USE_MKL)
    target_link_libraries(universal PRIVATE MKL::MKL_SCALAPACK)
  else()
    find_package(scalapack REQUIRED)
    target_link_libraries(universal PUBLIC scalapack)
  endif()
endif()

target_link_libraries(universal PUBLIC Eigen3::Eigen)
target_link_libraries(universal PUBLIC GSL::gsl)
target_link_libraries(universal PUBLIC absl::flat_hash_map)
target_include_directories(universal PRIVATE ${Boost_INCLUDE_DIRS})
