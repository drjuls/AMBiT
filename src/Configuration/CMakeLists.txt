set(MODS_CONFIGURATION  AngularData.cpp
                        ConfigGenerator.cpp
                        ElectronInfo.cpp
                        HamiltonianID.cpp
                        HamiltonianMatrix.cpp
                        LevelMap.cpp
                        NonRelConfiguration.cpp
                        Projection.cpp
                        RelativisticConfiguration.cpp
                        RelativisticConfigList.cpp
                        CACHE INTERNAL "")

add_library(configuration "${MODS_CONFIGURATION}")

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(configuration PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(configuration PUBLIC OpenMP::OpenMP_CXX)
endif()

if(USE_SCALAPACK)
  target_compile_definitions(configuration PRIVATE -DAMBIT_USE_SCALAPACK)
  if(USE_MKL)
    target_link_libraries(configuration PRIVATE MKL::MKL_SCALAPACK)
  else()
    find_package(scalapack REQUIRED)
    target_link_libraries(configuration PUBLIC scalapack)
  endif()
endif()

target_link_libraries(configuration PUBLIC Eigen3::Eigen)
target_include_directories(configuration PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(configuration PUBLIC mbpt)
