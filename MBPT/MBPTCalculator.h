#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "HartreeFock/SigmaPotential.h"

class MBPTCalculator
{
public:
    MBPTCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    ~MBPTCalculator(void) {}

    /** Create a sigma operator for the given state to second order.
     */
    void GetSecondOrderSigma(int kappa, SigmaPotential* sigma) const;

    /** Create a sigma operator for the given state to second order.
        if(sigma == NULL)
            no sigma potential is created, the function just calculates the matrix element.
        Return value is the energy of the state including second order matrix element.
     */
    double GetSecondOrderSigma(const State* s, SigmaPotential* sigma = NULL) const;

    /** Return value is the matrix element <s1 | Sigma | s2>. */
    double GetSecondOrderSigma(const State* s1, const State* s2, SigmaPotential* sigma = NULL) const;

    /** Use Brillouin-Wigner perturbation theory, where the energy of external lines is
        kept constant in the energy denominator (this ensures that the operator is hermitian).
        The energy of the external single-particle is supplied by the user, and should
        usually correspond to the leading configuration.
     */
    void UseBrillouinWignerPT(double single_particle_energy)
    {   ValenceEnergy1 = single_particle_energy;
        BrillouinWignerPT = true;
    }

    /** Use Rayleigh-Shrodinger perturbation theory. */
    void UseRayleighShrodingerPT()
    {   BrillouinWignerPT = false;
    }

protected:
    /** Calculate diagrams of second order (shown here 1 through 4).
        ->>------>------>>--  ->>------>------<---  ->>------<------>>--  ->>------<------>---
         i   |   3   |   f     i   |   3   |   2     i   |   3   |   f     i   |   3   |   4  
             |       |             |       |             |       |             |       |      
        --<------>------<---  --<------>------>>--  --<------>------<---  -->------<------>>--
         2       4       2     2       4       f     2       4       2     4       2       f  

         PRE: si.kappa == sf.kappa
     */
    double CalculateCorrelation1and3(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateCorrelation2(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateCorrelation4(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;

protected:
    Lattice* lattice;
    const Core* core;
    const ExcitedStates* excited;

    bool BrillouinWignerPT;
    double ValenceEnergy1;
};

#endif
