#ifndef HARTREE_FOCKER_H
#define HARTREE_FOCKER_H

#include "HFOperator.h"

class HartreeFocker
{
public:
    HartreeFocker() {}

    /** Iterate all orbitals in core until self-consistency is reached. */
    void SolveCore(Core* core, HFOperator* hf, SpinorODE* hf_decorated = NULL);

    /** Find self-consistent solution to hf operator, including exchange.
        Return change in energy.
     */
    double SolveOrbital(Orbital* orbital, SpinorODE* hf);

    /** Find energy eigenvalue for orbital with a given exchange potential.
        If exchange is not given, generate from hf.
        Note: this function does not iterate/update the exchange potential, 
        so the final orbital is not an eigenvalue of the hf operator.
        Return change in energy.
     */
    double IterateOrbital(Orbital* orbital, SpinorODE* hf, SpinorFunction* exchange = NULL);

protected:
    unsigned int MaxHFIterations = 1000;
    double WavefunctionTolerance = 1.E-11;
    double WavefunctionEnergyTolerance = 1.E-14;
};

#endif
