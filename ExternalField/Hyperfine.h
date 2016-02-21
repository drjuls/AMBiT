#ifndef HYPERFINE_H
#define HYPERFINE_H

#include "HartreeFock/SpinorOperator.h"
#include "Transitions.h"

/** Hyperfine operator based on relativistic Johnson formula. */
class HyperfineDipoleOperator : public SpinorOperator
{
public:
    /** Use nuclear_magnetic_radius_fm to assume uniformly magnetised spherical nucleus (in fm).
        Default nuclear_magnetic_radius = 0, corresponding to a point-like nucleus.
     */
    HyperfineDipoleOperator(pIntegrator integration_strategy, double nuclear_magnetic_radius_fm = 0.0);

public:
    /** Hyperfine shift with nuclear magneton in atomic units. */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    double nuclear_radius;
    unsigned int nuclear_radius_lattice;
    pLattice lattice;
};

class HyperfineDipoleCalculator : public TransitionCalculator
{
public:
    HyperfineDipoleCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double g_I;
};

/** Hyperfine electric quadrupole operator based on relativistic Johnson formula.
        t = -C^(2)/r^3
 */
class HyperfineQuadrupoleOperator: public SpinorOperator
{
public:
    HyperfineQuadrupoleOperator(pIntegrator integration_strategy);

public:
    /** Hyperfine shift for Q = 1 barn in atomic units. */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pLattice lattice;
};

class HyperfineQuadrupoleCalculator : public TransitionCalculator
{
public:
    HyperfineQuadrupoleCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double Q;
};

#endif
