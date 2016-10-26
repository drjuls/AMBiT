#ifndef SPINOR_ODE_H
#define SPINOR_ODE_H

#include "Universal/SpinorFunction.h"
#include "Orbital.h"
#include "Core.h"

/** SpinorODE is an abstract class for numerical integration of coupled linear
    ordinary differential equations (ODEs) of the form
        df/dr = w[0] = w_f[0] f + w_g[0] g + w_const[0]
        dg/dr = w[1] = w_f[1] f + w_g[1] g + w_const[1]
    where w is a linear function of f and g (e.g. Hartree-Fock).
    w_const is the "nonlocal" part (the exchange part).
    It follows the Decorator (Wrapper) pattern, so it is recursively extensive.
    It is also a LatticeObserver, and is guaranteed to provide correct differential equations over
    the entire lattice.
 */
class SpinorODE : public LatticeObserver
{
public:
    SpinorODE(pLattice lattice);
    virtual ~SpinorODE() {}

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, const SpinorFunction* exchange = nullptr) = 0;
    
    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) = 0;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation = nullptr) const = 0;

    /** Tell SpinorODE whether to include the nonlocal (w_const) terms in GetODEFunction, GetODECoefficients, and GetODEJacobian. */
    virtual void IncludeExchange(bool include_exchange);

    /** Interrogate whether to include nonlocal (w_const) terms in GetODEFunction, GetODECoefficients, and GetODEJacobian. */
    virtual bool IncludeExchange() const;

    /** Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g).
        PRE: w should be an allocated 2 dimensional array;
             latticepoint < size().
     */
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const = 0;

    /** Get numerical coefficients of the ODE at the point r, (f,g).
        w_f and w_g are coefficients of f and g in w; w_const is the constant term of w (not proportional to f or g).
        PRE: w_f, w_g, and w_const should be allocated 2 dimensional arrays;
             latticepoint < size().
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const = 0;
    
    /** Get Jacobian (dw[i]/df and dw[i]/dg), dw[i]/dr at a point r, (f, g).
        PRE: jacobian should be an allocated 2x2 matrix;
             dwdr and w_const should be allocated 2 dimensional arrays;
             latticepoint < size().
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const = 0;
    
    /** Get approximation to solution for first numpoints near the origin. */
    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const = 0;
    
    /** Get approximation to solution for last numpoints far from the origin.
        This routine can change the size of the orbital.
     */
    virtual void EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const = 0;

    /** Get df/dr and dg/dr given (f, g).
        If set_parameters is true, this function will call SetODEParameters(), changing exchange and include_exchange.
     */
    virtual void GetDerivative(Orbital& fg, bool set_parameters);
    virtual void GetDerivative(Orbital& fg) const;

protected:
    bool include_nonlocal;
};

typedef std::shared_ptr<SpinorODE> pSpinorODE;
typedef std::shared_ptr<const SpinorODE> pSpinorODEConst;

#endif
