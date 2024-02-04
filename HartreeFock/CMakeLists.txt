set(MODS_HARTREEFOCK  ConfigurationParser.cpp
                      Core.cpp
                      CoulombOperator.cpp
                      ExchangeDecorator.cpp
                      GreensMethodODE.cpp
                      HartreeFocker.cpp
                      HartreeY.cpp
                      HFOperator.cpp
                      LocalPotentialDecorator.cpp
                      NucleusDecorator.cpp
                      Integrator.cpp
                      Orbital.cpp
                      OrbitalInfo.cpp
                      OrbitalMap.cpp
                      ODESolver.cpp
                      SpinorODE.cpp
                      ThomasFermiDecorator.cpp
                      CACHE INTERNAL "")

add_library(hartreefock "${MODS_HARTREEFOCK}")

if(USE_MPI)
  target_compile_definitions(hartreefock PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(hartreefock PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(hartreefock PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(hartreefock PUBLIC OpenMP::OpenMP_CXX)
endif()
