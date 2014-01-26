#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include "SpinorODE.h"

/** OneDimensionalODE is an abstract class for numerical integration of a simple one-dimensional
    linear ordinary differential equations (ODE) of the form
        df/dr = w[0] = w_f(r) f + w_const(r)
 */
class OneDimensionalODE
{
public:
    OneDimensionalODE(pLattice lattice): lattice(lattice) {}

    /** Get df/dr = w[0] given point r, f.
        PRE: w should be an allocated double.
     */
    virtual void GetODEFunction(unsigned int latticepoint, const RadialFunction& f, double* w) const = 0;
    
    /** Get numerical coefficients of the ODE at the point r, f.
        PRE: w_f and w_const should be allocated 2 dimensional arrays.
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const RadialFunction& f, double* w_f, double* w_const) const = 0;
    
    /** Get Jacobian dw[i]/df and dw[i]/dr at a point r, f.
        PRE: jacobian and dwdr should allocated doubles.
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const RadialFunction& f, double* jacobian, double* dwdr) const = 0;
    
    /** Get approximation to solution for first numpoints near the origin. */
    virtual void EstimateSolutionNearOrigin(unsigned int numpoints, RadialFunction& f) const = 0;
    
    /** Get approximation to solution for last numpoints far from the origin. */
    virtual void EstimateSolutionNearInfinity(unsigned int numpoints, RadialFunction& f) const = 0;

protected:
    pLattice lattice;
};

typedef boost::shared_ptr<OneDimensionalODE> pOneDimensionalODE;
typedef boost::shared_ptr<const OneDimensionalODE> pOneDimensionalODEConst;

class ODESolver
{
public:
    ODESolver(pLattice lat): lattice(lat) {}
    virtual ~ODESolver() {}

public:
    /** Get eigenstate of the operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const OneDimensionalODE* op, RadialFunction* solution) = 0;
    
    /** Get eigenstate of the operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const OneDimensionalODE* op, RadialFunction* solution) = 0;
    
public:
    /** Get solution of the ODE operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(pSpinorODEConst op, SpinorFunction* solution) = 0;

    /** Get solution of the ODE operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(pSpinorODEConst op, Orbital* solution) = 0;

protected:
    pLattice lattice;
};

typedef boost::shared_ptr<ODESolver> pODESolver;
typedef boost::shared_ptr<const ODESolver> pODESolverConst;

class AdamsSolver : public ODESolver
{
public:
    AdamsSolver(pLattice lat);
    virtual ~AdamsSolver();

public:
    /** Get eigenstate of the operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const OneDimensionalODE* op, RadialFunction* solution);
    
    /** Get eigenstate of the operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const OneDimensionalODE* op, RadialFunction* solution);
    
public:
    /** Get solution of the ODE operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(pSpinorODEConst op, SpinorFunction* solution);
    
    /** Get solution of the ODE operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(pSpinorODEConst op, Orbital* solution);

public:
    /** Get solution from r=rmax backwards to the first maximum.
        Return lattice position of peak.
     */
    virtual unsigned int IntegrateBackwardsUntilPeak(pSpinorODEConst op, Orbital* solution);
    
protected:
    unsigned int order;
    double* adams_coeff;
};

#endif
