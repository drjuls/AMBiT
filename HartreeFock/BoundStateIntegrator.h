#ifndef BOUND_STATE_INTEGRATOR_H
#define BOUND_STATE_INTEGRATOR_H

#include "StateIntegrator.h"

class BoundStateIntegrator :
    public StateIntegrator
{
public:
    BoundStateIntegrator(Lattice& lat): StateIntegrator(lat) {}
    virtual ~BoundStateIntegrator(void) {}

protected:
    virtual void SetUpForwardsIntegral(State& s, const std::vector<double>& HFPotential, double nuclear_charge);
    virtual void SetUpBackwardsIntegral(State& s, const std::vector<double>& HFPotential);

protected:
    void SolveDiracBoundary(State& s, const std::vector<double>& HFPotential, unsigned int start_point, bool forwards = true);
};

#endif
