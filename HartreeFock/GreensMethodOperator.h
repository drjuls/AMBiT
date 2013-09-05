#ifndef GREENS_METHOD_OPERATOR_H
#define GREENS_METHOD_OPERATOR_H

#include "Operator.h"
#include "SpinorODE.h"
#include "ODESolver.h"

/** Operator for solving equations of the form
        L psi = source(r)
    where L is a first order differential operator.
    GreensMethodOperator is supplied with solutions to the homogenous equation
        L psi = 0
    which are regular at the origin and at infinity, s0 and sInf.
    Then one can express the solution of the original equation as a linear combination
        psi = s0 * A(r) + sInf * B(r).
    A(r) and B(r) will be proportional to functions of the form
        d(GInf)/dr = source(r) * sInf(r) / W
        d(G0)/dr = source(r) * s0(r) / W
    where W is the Wronskian = 1/(s0.f * sInf.g - sInf.f * s0.g).

    For the operator form, rather than the ODE form,
        GInf(r) = Integral[ source(r') * sInf(r') / W , {r' = r, Infinity}]
        G0(r)   = Integral[ source(r') * s0(r') / W , {r' = 0, r}]
    which is easier to obtain from the ODE form. An ODESolver may be supplied to replace the default.
 */
class GreensMethodOperator : public OneDimensionalODE
{
public:
    GreensMethodOperator(Lattice* lattice);

public:
    virtual void SetHomogenousSolutions(const SpinorFunction* fromOrigin, const SpinorFunction* fromInfinity);
    virtual void SetSourceTerm(const SpinorFunction* sourceTerm, bool fromOrigin);

public:
    /** Get df/dr = w[0] given point r, f.
     PRE: w should be an allocated double.
     */
    virtual void GetODEFunction(unsigned int latticepoint, const std::vector<double>& f, double* w) const;
    
    /** Get numerical coefficients of the ODE at the point r, f.
     PRE: w_f and w_const should be allocated 2 dimensional arrays.
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const std::vector<double>& f, double* w_f, double* w_const) const;
    
    /** Get Jacobian dw[i]/df and dw[i]/dr at a point r, f.
     PRE: jacobian and dwdr should allocated doubles.
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const std::vector<double>& f, double* jacobian, double* dwdr) const;
    
    /** Get approximation to solution for first numpoints near the origin. */
    virtual void EstimateSolutionNearOrigin(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const;
    
    /** Get approximation to solution for last numpoints far from the origin. */
    virtual void EstimateSolutionNearInfinity(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const;

protected:
    const SpinorFunction* s0;
    const SpinorFunction* sInf;
    const SpinorFunction* source;
    bool solutionRegularAtOrigin;

    /** The Wronskian is supposed to be a constant, but here we include it as a function of r.
            wronskian(r) = f0 gInf - fInf g0
     */
    std::vector<double> wronskian;
};

#endif
