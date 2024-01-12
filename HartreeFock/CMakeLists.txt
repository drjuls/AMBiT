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
