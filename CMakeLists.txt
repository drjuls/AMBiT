cmake_minimum_required(VERSION 3.18)

project(ambit VERSION 3.2
        DESCRIPTION "Relativistic CI+MBPT atomic structure code"
        LANGUAGES CXX Fortran)

enable_language(CXX)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  message(STATUS "Build type set to Release")
endif()

# Set default C++ and Fortran standards
set(CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(FVERSION "-std=f2003") # Modern Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} -cpp")

# Set the Fortran compiler. Note that this is case-sensitive: "Fortran" not "FORTRAN"
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

###### AMBiT stuff ###### 

# Set directory for angular data caching
set(ANGULAR_DATA_DIR "${CMAKE_SOURCE_DIR}/AngularData" CACHE STRING "Directory to store angular data")
add_definitions(-DANGULAR_DATA_DIRECTORY=${ANGULAR_DATA_DIR})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory(src)

############# Unit tests ############# 
# Now do unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # Unit tests don't (currently) work with MPI, so abort if we're building with MPI
  if(USE_MPI)
    message(FATAL_ERROR "AMBiT test suite does not support MPI. Disable MPI support to build tests")
  endif()

  enable_testing()
  add_subdirectory(test)

  # Copy the template input file (used for testing MultirunOptions) to the build directory so CTest
  # can find it
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/template.input" "${CMAKE_CURRENT_BINARY_DIR}/test/template.input" COPYONLY)
endif()

