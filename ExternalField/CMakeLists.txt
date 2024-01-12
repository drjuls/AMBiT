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
