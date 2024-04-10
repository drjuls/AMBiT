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

target_link_libraries(configuration PUBLIC Eigen3::Eigen)
target_include_directories(configuration PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(configuration PUBLIC mbpt)
