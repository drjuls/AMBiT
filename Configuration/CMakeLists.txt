set(MODS_CONFIGURATION  AngularData.cpp
                        ConfigGenerator.cpp
                        ElectronInfo.cpp
                        HamiltonianMatrix.cpp
                        Level.cpp
                        LevelMap.cpp
                        NonRelConfiguration.cpp
                        Projection.cpp
                        RelativisticConfiguration.cpp
                        RelativisticConfigList.cpp
                        CACHE INTERNAL "")
add_library(configuration "${MODS_CONFIGURATION}")

#find_package(Eigen3 REQUIRED NO_MODULE)
#target_link_libraries(configuration PUBLIC Eigen3::Eigen)
#include_directories(${EIGEN_INCLUDE_DIR} PUBLIC)

