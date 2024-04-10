set(MODS_EXTERNALFIELD  BreitHFDecorator.cpp
                        BreitZero.cpp
                        EJOperator.cpp
                        FieldShift.cpp
                        Hyperfine.cpp
                        KineticEnergy.cpp
                        LorentzInvarianceT2.cpp
                        NormalMassShiftDecorator.cpp
		                    NuclearPolarisability.cpp
                        RadiativePotential.cpp
                        RPAOperator.cpp
                        RPAOrbital.cpp
                        RPASolver.cpp
                        SpecificMassShiftDecorator.cpp
                        Transitions.cpp
                        TwoBodySMSOperator.cpp
                        YukawaPotential.cpp
                        CACHE INTERNAL "")

add_library(externalfield "${MODS_EXTERNALFIELD}")

target_link_libraries(externalfield PRIVATE GSL::gsl)
target_include_directories(externalfield PRIVATE ${Boost_INCLUDE_DIRS})


target_link_libraries(externalfield PUBLIC configuration)
