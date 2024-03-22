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

# Set directory for angular data caching
set(ANGULAR_DATA_DIR "${CMAKE_SOURCE_DIR}/AngularData" CACHE STRING "Directory to store angular data")
add_definitions(-DANGULAR_DATA_DIRECTORY=${ANGULAR_DATA_DIR})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory(src)

############# Unit tests ############# 
# Now do unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # Download Googletest if we don't have it already
  enable_testing()
  add_subdirectory(test)

  # Copy the template input file (used for testing MultirunOptions) to the build directory so CTest
  # can find it
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/template.input" "${CMAKE_CURRENT_BINARY_DIR}/test/template.input" COPYONLY)
endif()

