#ifndef RADIATIVE_POTENTIAL_H
#define RADIATIVE_POTENTIAL_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/LocalPotentialDecorator.h"

/** Uehling potential decorator.
    Uses step potential from Ginges & Berengut.
 */
class UehlingDecorator: public LocalPotentialDecorator<UehlingDecorator>
{
public:
    UehlingDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pOPIntegrator integration_strategy = pOPIntegrator());

protected:
    /** Set local direct_potential to Uehling potential. */
    void GenerateStepUehling(double nuclear_rms_radius);
};

#endif
