#include "GreensMethodOperator.h"
#include "Include.h"

GreensMethodOperator::GreensMethodOperator(Lattice* lattice): OneDimensionalODE(lattice)
{
    s0 = NULL;
    sInf = NULL;
    source = NULL;
    solutionRegularAtOrigin = true;
}

/** Functions for ODE solving. */
void GreensMethodOperator::SetHomogenousSolutions(const SpinorFunction* fromOrigin, const SpinorFunction* fromInfinity)
{
    s0 = fromOrigin;
    sInf = fromInfinity;

    wronskian.resize(mmin(s0->Size(),sInf->Size()));
    unsigned int i;
    for(i = 0; i < wronskian.size(); i++)
    {
        wronskian[i] = s0->f[i] * sInf->g[i] - sInf->f[i] * s0->g[i];
    }
}

void GreensMethodOperator::SetSourceTerm(const SpinorFunction* sourceTerm, bool fromOrigin)
{
    source = sourceTerm;
    solutionRegularAtOrigin = fromOrigin;
}

/** Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g).
 PRE: w should be an allocated 2 dimensional array.
 */
void GreensMethodOperator::GetODEFunction(unsigned int latticepoint, const std::vector<double>& f, double* w) const
{
    double integrand = 0;

    if(latticepoint < wronskian.size() && latticepoint < source->Size())
    {
        if(solutionRegularAtOrigin)
            integrand = s0->f[latticepoint] * source->f[latticepoint] + s0->g[latticepoint] * source->g[latticepoint];
        else
            integrand = sInf->f[latticepoint] * source->f[latticepoint] + sInf->g[latticepoint] * source->g[latticepoint];

        integrand = integrand/wronskian[latticepoint];
    }

    *w = integrand;
}

/** Get numerical coefficients of the ODE at the point r, (f,g).
 PRE: w_f, w_g, and w_const should be allocated 2 dimensional arrays.
 */
void GreensMethodOperator::GetODECoefficients(unsigned int latticepoint, const std::vector<double>& f, double* w_f, double* w_const) const
{
    *w_f = 0.;

    GetODEFunction(latticepoint, f, w_const);
}

/** Get Jacobian (dw[i]/df and dw[i]/dg), dw[i]/dr, and w_const at a point r, (f, g).
 w_const is the constant term of w (not proportional to f or g).
 PRE: jacobian should be an allocated 2x2 matrix,
 dwdr and w_const should be allocated 2 dimensional arrays.
 */
void GreensMethodOperator::GetODEJacobian(unsigned int latticepoint, const std::vector<double>& f, double* jacobian, double* dwdr) const
{
    *jacobian = 0.;

    double w;
    const SpinorFunction* sregular;
    if(solutionRegularAtOrigin)
        sregular = s0;
    else
        sregular = sInf;

    w = sregular->dfdr[latticepoint] * source->f[latticepoint] + sregular->f[latticepoint] * source->dfdr[latticepoint]
        + sregular->dgdr[latticepoint] * source->g[latticepoint] + sregular->g[latticepoint] * source->dgdr[latticepoint];

    *dwdr = w;
}

/** Get approximation to eigenfunction for first numpoints near the origin. */
void GreensMethodOperator::EstimateSolutionNearOrigin(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const
{
    double w_i;
    const double* dR = lattice->dR();

    unsigned int i = 0;

    GetODEFunction(i, f, &w_i);
    f[i] = 0.5 * w_i * dR[i];
    dfdr[i] = w_i;

    for(i=1; i < numpoints; i++)
    {
        GetODEFunction(i, f, &w_i);
        dfdr[i] = w_i;

        f[i] = f[i-1] + 0.5 * (dfdr[i-1] * dR[i-1] + dfdr[i] * dR[i]);
    }
}

/** Get approximation to eigenfunction for last numpoints far from the origin.
 This routine can change the size of the orbital.
 */
void GreensMethodOperator::EstimateSolutionNearInfinity(unsigned int numpoints, std::vector<double>& f, std::vector<double>& dfdr) const
{
    double w_i;
    const double* dR = lattice->dR();
    
    unsigned int i = f.size()-1;
    
    GetODEFunction(i, f, &w_i);
    f[i] = 0.5 * w_i * dR[i];
    dfdr[i] = w_i;
    
    for(i=f.size()-2; i >= f.size()-numpoints; i--)
    {
        GetODEFunction(i, f, &w_i);
        dfdr[i] = w_i;
        
        f[i] = f[i+1] - 0.5 * (dfdr[i+1] * dR[i+1] + dfdr[i] * dR[i]);
    }
}
