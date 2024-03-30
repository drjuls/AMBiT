# First, generate the gitinfo header file and dummy target
execute_process(COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/getGitInfo.sh" 
                OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/gitInfo.h")
              add_custom_target(generate_gitinfo ALL DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/gitInfo.h")

# Main ambit exe
add_executable(ambit ambit.cpp ambit.h)
# Need to explicitly add the dependency on gitInfo.h to ensure it gets built before we do anything
# else
add_dependencies(ambit generate_gitinfo)

# External libraries
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(absl CONFIG REQUIRED)

# LAPACK and BLAS
# First, check to see if we need to use MKL
set(USE_MKL OFF CACHE BOOL "Use Intel MKL for BLAS/LAPACK subroutines")
if(USE_MKL)
  # MKL 2021.3 and later has a nice CMake support, so make sure we use that
  # Documentation can be found here:
  # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/cmake-config-for-onemkl.html
  # and
  # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/using-the-cmake-config-file.html
  # Get the path from the shell or CMake command-line options (e.g. set -DMKL_ROOT=/path/to/MKL)
  find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})  
  target_link_libraries(ambit PRIVATE MKL::MKL)
  # Also tell Eigen to use MKL as its backend
  target_compile_definitions(ambit PRIVATE -DEIGEN_USE_MKL_ALL)
else() # Use vanilla BLAS/LAPACK
  find_package(LAPACK REQUIRED)

  target_link_libraries(ambit PRIVATE LAPACK::LAPACK)
  find_package(BLAS REQUIRED)
  target_link_libraries(ambit PRIVATE BLAS::BLAS)
endif()

# MPI options. Note that we can't automatically find the MPI compiler wrappers (e.g. mpiCC), so the
# user will have to manually set CMAKE_CXX_COMPILER=mpiCC or whatever (TODO EVK: put this in the
# instructions)

# Important note: our user-facing "CACHE" variables (i.e. we can set them with cmake -DUSE_MPI=yes)
# have to be different from the internal AMBiT preprocessor definitions. AMBiT just checks if
# AMBIT_USE_(OPENMP|MPI) are defined, regardless of its actual value. If we expose these variables
# in the CMake cache then they'll get set to FALSE, but the preprocessor will see that they exist
# and thus will always generate the OpenMP/MPI code, even though that's not what the user has
# requested. Consequently, we need "PRIVATE" compile definitions that shadow the "CACHE" ones in
# order to make this work properly.
set(USE_MPI OFF CACHE BOOL "Enable MPI support in AMBiT")
if(USE_MPI)
  target_compile_definitions(ambit PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(ambit PUBLIC MPI::MPI_CXX)
endif()

# OpenMP options. CMake should be able to automatically infer the correct flags for the current
# compiler
set(USE_OPENMP OFF CACHE BOOL "Enable OpenMP support in AMBiT")
if(USE_OPENMP)
  target_compile_definitions(ambit PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(ambit PUBLIC OpenMP::OpenMP_CXX)
endif()

# Now include the directories for the different ambit modules and link against the static libs
add_subdirectory(Atom)
add_subdirectory(Basis)
add_subdirectory(Configuration)
add_subdirectory(ExternalField)
add_subdirectory(HartreeFock)
add_subdirectory(MBPT)
add_subdirectory(Universal)

target_link_libraries(ambit PUBLIC atom basis configuration externalfield hartreefock mbpt universal)

# Sparsehash
# EVK TODO: Sparsehash doesn't have support for CMake, which makes automatically finding it kind of
# annoying. I can think of a few options:
# 1.  Make the user specify a root directory for Sparsehash and just bail out if we can't find it
#     (this is the current behaviour with SCons)
# 1a. Optional: download Sparsehash if we can't find it on the host system
# 2.  There's a C++11 version of sparsehash which has CMake support: 
#     https://github.com/sparsehash/sparsehash-c11/. We could switch to using this, although I
#     don't know if there will be API incompatibilities or something
# 3.  Switch to something like Abseil: https://github.com/abseil/abseil-cpp. This is the successor
#     library to Sparsehash (also from Google) and (allegedly) has better performance than
#     Sparsehash but would require some code changes for the new API. We could also automatically
#     download Abseil if the user doesn't supply a location (which would be easier to do since we
#     can add it as a submodule using CMake)
#
# I'm going to use the first option for now, but I'm leaning towards replacing this with Abseil
# once I have the CMake build system done.
set(SPARSEHASH_ROOT "" CACHE STRING "Sparsehash root installation directory")
target_include_directories(ambit PRIVATE "${SPARSEHASH_ROOT}/include")

# Install the "ambit" exe into "./bin".
install(TARGETS ambit DESTINATION "bin")
