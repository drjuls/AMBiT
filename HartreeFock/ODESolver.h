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
    OneDimensionalODE(Lattice* lattice): lattice(lattice) {}

    /** Get df/dr = w[0] given point r, f.
        PRE: w should be an allocated double.
     */
    virtual void GetODEFunction(unsigned int latticepoint, const std::vector<double>& f, double* w) const = 0;
    
    /** Get numerical coefficients of the ODE at the point r, f.
        PRE: w_f and w_const should be allocated 2 dimensional arrays.
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const std::vector<double>& f, double* w_f, double* w_const) const = 0;
    
    /** Get Jacobian dw[i]/df and dw[i]/dr at a point r, f.
        PRE: jacobian and dwdr should allocated doubles.
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const std::vector<double>& f, double* jacobian, double* dwdr) const = 0;
    
    /** Get approximation to solution for first numpoints near the origin. */
    virtual void EstimateSolutionNearOrigin(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const = 0;
    
    /** Get approximation to solution for last numpoints far from the origin. */
    virtual void EstimateSolutionNearInfinity(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const = 0;

protected:
    Lattice* lattice;
};

class ODESolver
{
public:
    ODESolver(Lattice* lat): lattice(lat) {}
    virtual ~ODESolver() {}

public:
    /** Get eigenstate of the operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const OneDimensionalODE* op, std::vector<double>* sol, std::vector<double>* dsoldr) = 0;
    
    /** Get eigenstate of the operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const OneDimensionalODE* op, std::vector<double>* sol, std::vector<double>* dsoldr) = 0;
    
public:
    /** Get solution of the ODE operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const SpinorODE* op, SpinorFunction* solution) = 0;

    /** Get solution of the ODE operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const SpinorODE* op, Orbital* solution) = 0;

protected:
    Lattice* lattice;
};


class AdamsSolver : public ODESolver
{
public:
    AdamsSolver(Lattice* lat);
    virtual ~AdamsSolver();

public:
    /** Get eigenstate of the operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const OneDimensionalODE* op, std::vector<double>* sol, std::vector<double>* dsoldr);
    
    /** Get eigenstate of the operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const OneDimensionalODE* op, std::vector<double>* sol, std::vector<double>* dsoldr);
    
public:
    /** Get solution of the ODE operator by integrating from r=0 to rmax. */
    virtual void IntegrateForwards(const SpinorODE* op, SpinorFunction* solution);
    
    /** Get solution of the ODE operator by integrating from r=rmax to 0. */
    virtual void IntegrateBackwards(const SpinorODE* op, Orbital* solution);

protected:
    unsigned int order;
    double* adams_coeff;
};

#endif
