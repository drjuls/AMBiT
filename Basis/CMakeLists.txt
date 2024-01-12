set(MODS_BASIS BasisGenerator.cpp
               BSplineBasis.cpp
               TimesRBasis.cpp
               HartreeFockBasis.cpp
               OrbitalManager.cpp
               CACHE INTERNAL "")

add_library(basis "${MODS_BASIS}")

#find_package(Eigen3 REQUIRED NO_MODULE)
#target_link_libraries(basis PUBLIC Eigen3::Eigen)
#include_directories(${EIGEN_INCLUDE_DIR} PUBLIC)

