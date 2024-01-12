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

#find_package(Eigen3 REQUIRED NO_MODULE)
#target_link_libraries(universal PUBLIC Eigen3::Eigen)
#include_directories(${EIGEN_INCLUDE_DIR} PUBLIC)

