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

# Set the internal MPI and OpenMP #ifdefs for this library
if(USE_MPI)
  target_compile_definitions(externalfield PRIVATE -DAMBIT_USE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(externalfield PUBLIC MPI::MPI_CXX)
endif()

if(USE_OPENMP)
  target_compile_definitions(externalfield PRIVATE -DAMBIT_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(externalfield PUBLIC OpenMP::OpenMP_CXX)
endif()
