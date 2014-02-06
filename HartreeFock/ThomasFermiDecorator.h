#ifndef THOMAS_FERMI_DECORATOR_H
#define THOMAS_FERMI_DECORATOR_H

#include "HFOperator.h"

/** Replace part of Hartree-Fock potential with simple Thomas-Fermi type potential,
    where direct part \f$ V = Z_{eff}/r \f$ and exchange is always zero
    (it may be approximated by wrapping LocalExchangeApproximation).
    The user can supply a change factor to SetCore(), to slowly transition to a Hartree-Fock potential.
    The effective charge for the Thomas-Fermi potential is a monotonically decreasing function that is
    Z at the origin and charge at infinity (for neutral atoms this charge = 1).
    The final (mixed with HF) potential is forced to have a minimum value of charge/r at all points.
 */
class ThomasFermiDecorator: public HFOperatorDecorator
{
public:
    ThomasFermiDecorator(pHFOperator decorated_object, pOPIntegrator integration_strategy = pOPIntegrator());

    virtual void SetCore(pCoreConst hf_core, double hf_mixing = 0.0);
    virtual RadialFunction GetDirectPotential() const;
    virtual void ExtendPotential();

    virtual SpinorFunction GetExchange(pSingleParticleWavefunctionConst approximation = pSingleParticleWavefunctionConst()) const;
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const;
};

typedef boost::shared_ptr<ThomasFermiDecorator> pThomasFermiDecorator;
typedef boost::shared_ptr<const ThomasFermiDecorator> pThomasFermiDecoratorConst;

#endif
