set(MODS_MBPT BruecknerDecorator.cpp
              BruecknerSigmaCalculator.cpp
              CoreMBPTCalculator.cpp
              MBPTCalculator.cpp
              OneElectronMBPT.cpp
              Sigma3Calculator.cpp
              SigmaPotential.cpp
              TwoElectronCoulombOperator.cpp
              ValenceMBPTCalculator.cpp
              CACHE INTERNAL "")

add_library(mbpt "${MODS_MBPT}")

target_link_libraries(mbpt PUBLIC absl::flat_hash_map)

target_link_libraries(mbpt PUBLIC basis)
