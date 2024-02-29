cmake_minimum_required(VERSION 3.18)

project(ambit VERSION 3.1.1
        DESCRIPTION "Relativistic CI+MBPT atomic structure code"
        LANGUAGES C CXX Fortran)

enable_language(CXX)
enable_language(Fortran)

# Set default C++ and Fortran standards
set(CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FVERSION "-std=gnu") # Modern Fortran with GNU extensions
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} -cpp")

# Set the Fortran compiler. Note that this is case-sensitive: "Fortran" not "FORTRAN"
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


###### AMBiT stuff ###### 

# First, generate the gitinfo header file and dummy target
execute_process(COMMAND bash "${CMAKE_SOURCE_DIR}/getGitInfo.sh" 
                OUTPUT_FILE "${CMAKE_SOURCE_DIR}/gitInfo.h")
add_custom_target(generate_gitinfo ALL DEPENDS "${CMAKE_SOURCE_DIR}/gitInfo.h")

# Set directory for angular data caching
set(ANGULAR_DATA_DIR "${CMAKE_SOURCE_DIR}/AngularData" CACHE STRING "Directory to store angular data")
add_definitions(-DANGULAR_DATA_DIRECTORY=${ANGULAR_DATA_DIR})


# Main ambit exe
add_executable(ambit ambit.cpp ambit.h)
# Need to explicitly add the dependency on gitInfo.h to ensure it gets built before we do anything
# else
add_dependencies(ambit generate_gitinfo)
include_directories(${CMAKE_SOURCE_DIR})

# Eigen headers. This needs to be here so the sub-modules can find it
# TODO EVK: this seems weird, is this a problem with Eigen?
find_package(Eigen3 REQUIRED NO_MODULE)
if(EIGEN3_FOUND)
  target_link_libraries(ambit PUBLIC Eigen3::Eigen)
  include_directories(${EIGEN_INCLUDE_DIR} PUBLIC)
else()
  message(FATAL_ERROR "Could not find Eigen. Set EIGEN_INCLUDE_DIR to the Eigen installation location")
endif()

# Now include the directories for the different ambit modules and link against the static libs
add_subdirectory(Atom)
add_subdirectory(Basis)
add_subdirectory(Configuration)
add_subdirectory(ExternalField)
add_subdirectory(HartreeFock)
add_subdirectory(MBPT)
add_subdirectory(Universal)

# TODO EVK: Should probably transition these to using CMake namespaces
target_link_libraries(ambit PRIVATE atom)
target_link_libraries(ambit PRIVATE basis)
target_link_libraries(ambit PRIVATE configuration)
target_link_libraries(ambit PRIVATE externalfield)
target_link_libraries(ambit PRIVATE hartreefock)
target_link_libraries(ambit PRIVATE mbpt)
target_link_libraries(ambit PRIVATE universal)

# External libraries
# GSL
find_package(GSL REQUIRED)
target_link_libraries(ambit PUBLIC GSL::gsl)

# Boost
find_package(Boost REQUIRED)
           target_include_directories(ambit PUBLIC ${Boost_INCLUDE_DIRS})

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
  target_link_libraries(ambit PUBLIC MKL::MKL)
  # Also tell Eigen to use MKL as its backend
  target_compile_definitions(ambit PRIVATE -DEIGEN_USE_MKL_ALL)
else() # Use vanilla BLAS/LAPACK
  find_package(LAPACK REQUIRED)

  if(LAPACK_FOUND)
    target_link_libraries(ambit PUBLIC LAPACK::LAPACK)
  else()
    message(FATAL_ERROR "Could not find LAPACK libraries. Set LAPACK_LIBRARIES CMake variable to LAPACK installation directory.")
  endif()

  find_package(BLAS REQUIRED)
  if(BLAS_FOUND)
    target_link_libraries(ambit PUBLIC BLAS::BLAS)
  else()
    message(FATAL_ERROR "Could not find BLAS libraries. Set BLAS_LIBRARIES CMake variable to BLAS installation directory.")
  endif()
endif()
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
#add_definitions(-DSPARSEHASH_ROOT=${SPARSEHASH_ROOT})
target_include_directories(ambit PUBLIC "${SPARSEHASH_ROOT}/include")

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


# Install the "ambit" exe into "./bin".
# TODO EVK: check what Julian wants to do here
install(TARGETS ambit DESTINATION "bin")

############# Unit tests ############# 
# Now do unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # Download Googletest if we don't have it already
  enable_testing()
  add_subdirectory(UnitTest)

  # Copy the template input file (used for testing MultirunOptions) to the build directory so CTest
  # can find it
  configure_file("template.input" "${CMAKE_CURRENT_BINARY_DIR}/UnitTest/template.input" COPYONLY)

endif()
