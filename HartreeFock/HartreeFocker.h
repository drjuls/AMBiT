#ifndef HARTREE_FOCKER_H
#define HARTREE_FOCKER_H

#include "HFOperator.h"

class HartreeFocker
{
public:
    HartreeFocker(ODESolver* ode_solver): odesolver(ode_solver) {}

    /** Iterate all orbitals in core until self-consistency is reached. */
    void SolveCore(Core* core, pHFOperator hf);

    /** Find self-consistent solution to hf operator, including exchange.
        Return change in energy.
     */
    double SolveOrbital(pOrbital orbital, pHFOperator hf);

    unsigned int CalculateExcitedState(pOrbital orbital, pHFOperator hf);

    /** Find energy eigenvalue for orbital with a given exchange potential.
        If exchange is not given, generate from hf.
        Note: this function does not iterate/update the exchange potential, 
        so the final orbital is not an eigenvalue of the hf operator.
        Return change in energy.
     */
    double IterateOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunction exchange = pSpinorFunction());

    /** Find energy eigenvalue for orbital using tail matching method.
        Unlike the Greens method used in IterateOrbital(), this method doesn't require a
        source term (exchange term), so it can solve without exchange or with a local exchange approximation.
        Return number of loops.
     */
    unsigned int IterateOrbitalTailMatching(pOrbital orbital, pHFOperator hf);

protected:
    ODESolver* odesolver;

    unsigned int MaxHFIterations = 1000;
    double WavefunctionTolerance = 1.E-11;
    double EnergyTolerance = 1.E-14;
    double TailMatchingEnergyTolerance = 1.E-8;
};

#endif
