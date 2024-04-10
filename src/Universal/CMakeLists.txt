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

target_link_libraries(universal PUBLIC Eigen3::Eigen)
target_link_libraries(universal PUBLIC GSL::gsl)
target_link_libraries(universal PUBLIC absl::flat_hash_map)
target_include_directories(universal PRIVATE ${Boost_INCLUDE_DIRS})
