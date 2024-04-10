set(MODS_BASIS BasisGenerator.cpp
               BSplineBasis.cpp
               TimesRBasis.cpp
               HartreeFockBasis.cpp
               OrbitalManager.cpp
               CACHE INTERNAL "")

add_library(basis "${MODS_BASIS}")

target_link_libraries(basis PUBLIC Eigen3::Eigen)
target_include_directories(basis PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(basis PUBLIC hartreefock)
