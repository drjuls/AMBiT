#ifndef HYPERFINE_H
#define HYPERFINE_H

#include "HartreeFock/SpinorOperator.h"
#include "HartreeFock/HFOperator.h"
#include "RPAOrbital.h"
#include "HartreeFock/HartreeY.h"

/** Hyperfine operator based on relativistic Johnson formula. */
class HyperfineDipoleOperator : public SpinorOperator
{
public:
    /** Point-like nucleus. */
    HyperfineDipoleOperator(pOPIntegrator integration_strategy);

public:
    /** Hyperfine shift with nuclear magneton in atomic units. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pLattice lattice;
};

class HyperfineRPAOperator : public SpinorOperator
{
public:
    HyperfineRPAOperator(pCoreConst core, pHartreeY hartreeY):
        SpinorOperator(1, hartreeY->GetOPIntegrator()), hartreeY(hartreeY), core(core),
        hyperfine(hartreeY->GetOPIntegrator())
    {}

    virtual void SetCore(pCoreConst rpa_core) { core = rpa_core; }

    /** Return (f + deltaVhf)|a> */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    HyperfineDipoleOperator hyperfine;
    pCoreConst core;
    pHartreeY hartreeY;
};

typedef std::shared_ptr<HyperfineRPAOperator> pHyperfineRPAOperator;

/** Add external dipole field and change in HF operator to RHS of RPA equations.
    These are added only to DeltaOrbital;
    see Dzuba, Flambaum, Sushkov, J. Phys. B 17, 1953 (1984).
 */
class HyperfineDipoleRPADecorator : public HFOperatorDecorator<HFBasicDecorator, HyperfineDipoleRPADecorator>
{
public:
    HyperfineDipoleRPADecorator(pHFOperator wrapped, pHartreeY hartreeY, pOPIntegrator integration_strategy = pOPIntegrator());

    void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation) const override;

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    HyperfineDipoleOperator hyperfine;
    pHartreeY hartreeY;

    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const;
};

#endif
