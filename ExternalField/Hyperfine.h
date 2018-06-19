#ifndef HYPERFINE_H
#define HYPERFINE_H

#include "HartreeFock/SpinorOperator.h"
#include "Transitions.h"

namespace Ambit
{
/** Hyperfine electric operator based on relativistic formula
        t = -C^(k)/r^(k+1)
 */
class HyperfineEJOperator : public SpinorOperator
{
public:
    /** Use nuclear_magnetic_radius_fm to assume uniformly magnetised spherical nucleus (in fm).
        Default nuclear_radius_fm = 0, corresponding to a point-like nucleus.
     */
    HyperfineEJOperator(int J, pIntegrator integration_strategy, double nuclear_radius_fm = 0.0);

public:
    /** Hyperfine shift in atomic units. */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    double nuclear_radius;
    unsigned int nuclear_radius_lattice;
    pLattice lattice;
};

/** Hyperfine magnetic operator based on relativistic formula
        t = -i/r^(k+1) sqrt((k+1)/k) alpha.C^k (r)
 */
class HyperfineMJOperator : public SpinorOperator
{
public:
    /** Use nuclear_radius_fm to assume uniformly magnetised spherical nucleus (in fm).
        Default nuclear_radius_fm = 0, corresponding to a point-like nucleus.
     */
    HyperfineMJOperator(int J, pIntegrator integration_strategy, double nuclear_radius_fm = 0.0);

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

class HyperfineQuadrupoleCalculator : public TransitionCalculator
{
public:
    HyperfineQuadrupoleCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double Q;
};

class GeneralisedHyperfineCalculator : public TransitionCalculator
{
public:
    GeneralisedHyperfineCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    MultipolarityType EorM;
    double units;
};

}
#endif
