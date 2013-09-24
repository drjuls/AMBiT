#ifndef HARTREE_FOCKER_H
#define HARTREE_FOCKER_H

#include "HFOperator.h"

class HartreeFocker
{
public:
    HartreeFocker() {}
    
    void IterateCore(Core* core, HFOperator* hf, SpinorODE* hf_decorated = NULL);
    void IterateOrbital(Orbital* orbital, SpinorODE* hf);

protected:
    unsigned int MaxHFIterations = 1000;
    double WavefunctionTolerance = 1.E-11;
    double WavefunctionEnergyTolerance = 1.E-14;
};

#endif
