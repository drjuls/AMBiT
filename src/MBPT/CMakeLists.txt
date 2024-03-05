set(MODS_MBPT BruecknerDecorator.cpp
              BruecknerSigmaCalculator.cpp
              CoreMBPTCalculator.cpp
              MBPTCalculator.cpp
              OneElectronMBPT.cpp
              Sigma3Calculator.cpp
              SigmaPotential.cpp
              TwoElectronCoulombOperator.cpp
              ValenceMBPTCalculator.cpp
              CACHE INTERNAL "")

add_library(mbpt "${MODS_MBPT}")

target_link_libraries(mbpt PUBLIC unordered_dense::unordered_dense)

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(mbpt PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(mbpt PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(mbpt PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(mbpt PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(mbpt PUBLIC basis)
