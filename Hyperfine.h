#ifndef HYPERFINE_H
#define HYPERFINE_H

#include "HartreeFock/SpinorOperator.h"

class HyperfineDipoleOperator : public SpinorOperator
{
public:
    /** Point-like nucleus. */
    HyperfineDipoleOperator(pLattice lattice, pOPIntegrator integration_strategy);

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pLattice lattice;
};

#endif
