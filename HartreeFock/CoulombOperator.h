#ifndef COULOMB_OPERATOR
#define COULOMB_OPERATOR

#include "ODESolver.h"

/** Find the potential due to a charge density
                           r< ^k
        I(r) = Integral[ ---------.density(r').dr' ]
                         r> ^(k+1)
    Split into two parts,
                   1
        I1(r) = ------- . Integral[ r' ^k .density(r').dr' ]
                r ^(k+1)    0->r
                                      1
        I2(r) = r ^k . Integral  [ ------- .density(r').dr' ]
                      r->infinity  r'^(k+1)
 
        dI1/dr = -(k+1)/r .I1 + density/r
        dI2/dr =    k/r .I2   - density/r
 */
class CoulombOperator : public OneDimensionalODE
{
public:
    CoulombOperator(pLattice lattice, pODESolver ode = pODESolver());

    void SetK(int multipole_k) { k = multipole_k; }
    void SetDensity(const RadialFunction& density);

    unsigned int GetK() const { return k; }
    void GetPotential(int k, const RadialFunction& density, RadialFunction& pot, pODESolver ode = pODESolver());

    /** Get zero-multipole potential, but renormalise density so that potential function goes as charge/r at infinity. */
    void GetPotential(RadialFunction& density, RadialFunction& pot, double charge, pODESolver ode = pODESolver());

public:
    /** Get df/dr = w[0] given point r, f.
     PRE: w should be an allocated double.
     */
    virtual void GetODEFunction(unsigned int latticepoint, const RadialFunction& f, double* w) const;
    
    /** Get numerical coefficients of the ODE at the point r, f.
     PRE: w_f and w_const should be allocated 2 dimensional arrays.
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const RadialFunction& f, double* w_f, double* w_const) const;
    
    /** Get Jacobian dw[i]/df and dw[i]/dr at a point r, f.
     PRE: jacobian and dwdr should allocated doubles.
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const RadialFunction& f, double* jacobian, double* dwdr) const;
    
    /** Get approximation to solution for first numpoints near the origin. */
    virtual void EstimateSolutionNearOrigin(unsigned int numpoints, RadialFunction& f) const;
    
    /** Get approximation to solution for last numpoints far from the origin. */
    virtual void EstimateSolutionNearInfinity(unsigned int numpoints, RadialFunction& f) const;

protected:
    int k;
    RadialFunction rho;
    pODESolver ode_solver;
    bool fwd_direction;
};

typedef boost::shared_ptr<CoulombOperator> pCoulombOperator;
typedef boost::shared_ptr<const CoulombOperator> pCoulombOperatorConst;

#endif
