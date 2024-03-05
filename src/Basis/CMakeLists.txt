set(MODS_BASIS BasisGenerator.cpp
               BSplineBasis.cpp
               TimesRBasis.cpp
               HartreeFockBasis.cpp
               OrbitalManager.cpp
               CACHE INTERNAL "")

add_library(basis "${MODS_BASIS}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(basis PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(basis PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(basis PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(basis PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(basis PUBLIC Eigen3::Eigen)
target_include_directories(basis PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(basis PUBLIC hartreefock)
