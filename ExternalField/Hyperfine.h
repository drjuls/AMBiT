#ifndef HYPERFINE_H
#define HYPERFINE_H

#include "HartreeFock/SpinorOperator.h"

/** Hyperfine operator based on relativistic Johnson formula. */
class HyperfineDipoleOperator : public SpinorOperator
{
public:
    /** Use nuclear_magnetic_radius_fm to assume uniformly magnetised spherical nucleus (in fm).
        Default nuclear_magnetic_radius = 0, corresponding to a point-like nucleus.
     */
    HyperfineDipoleOperator(pOPIntegrator integration_strategy, double nuclear_magnetic_radius_fm = 0.0);

public:
    /** Hyperfine shift with nuclear magneton in atomic units. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    double nuclear_radius;
    unsigned int nuclear_radius_lattice;
    pLattice lattice;
};

#endif
