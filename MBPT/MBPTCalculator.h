#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "SigmaPotential.h"

class MBPTCalculator
{
public:
    MBPTCalculator(Lattice* lattice, Core* atom_core, ExcitedStates* excited_states);
    ~MBPTCalculator(void) {}

    /** Create a sigma operator for the given state to second order.
        if(sigma == NULL)
            no sigma potential is created, the function just calculates the second order correlations
        Return value is the energy of the state including second order matrix elements.
     */
    double GetSecondOrderSigma(const State* s, SigmaPotential* sigma = NULL);

protected:
    /** Calculate diagrams of second order (shown here 1 through 4).
        ->>------>------>>--  ->>------>------<---  ->>------<------>>--  ->>------<------>---
         1   |   3   |   1     1   |   3   |   2     1   |   3   |   1     1   |   3   |   4  
             |       |             |       |             |       |             |       |      
        --<------>------<---  --<------>------>>--  --<------>------<---  -->------<------>>--
         2       4       2     2       4       1     2       4       2     4       2       1  
     */
    double CalculateCorrelation1and3(const State* s, SigmaPotential* sigma = NULL);
    double CalculateCorrelation2(const State* s, SigmaPotential* sigma = NULL);
    double CalculateCorrelation4(const State* s, SigmaPotential* sigma = NULL);

protected:
    Lattice* lattice;
    Core* core;
    ExcitedStates* excited;
};

#endif
