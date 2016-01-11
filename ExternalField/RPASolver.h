#ifndef RPA_SOLVER_H
#define RPA_SOLVER_H

#include "HartreeFock/Core.h"
#include "HartreeFock/HFOperator.h"
#include "Hyperfine.h"
#include "RPAOrbital.h"

class RPASolver
{
public:
    RPASolver(): remove_spin_orbit_partner(false) {}
    ~RPASolver() {}

    /** Populate the DeltaOrbitals self-consistently.
        PRE: core should be self-consistent solution of hf.
     */
    void SolveRPACore(pCore core, pHFOperator hf, pHyperfineRPAOperator rpa, bool include_negative_basis = true);

    /** Return RPA energy correction to excited state. */
    double CalculateRPAExcited(pRPAOrbital orbital, pHyperfineRPAOperator rpa);

    /** Iterate a single deltaOrbital.
        Return change in its deltaEnergy.
     */
    double IterateDeltaOrbital(pDeltaOrbital orbital, pHyperfineRPAOperator rpa);

    double EnergyTolerance = 1.e-14;

protected:
    std::map<int, pOrbitalMap> basis;   //!< DeltaOrbital basis for each kappa
    bool remove_spin_orbit_partner;     //!< Remove spin-orbit partner from basis (small energy denominators can be fatal)
    unsigned int MaxRPAIterations = 100;
};

#endif
