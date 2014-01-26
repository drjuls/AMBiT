#include "CoulombOperator.h"

CoulombOperator::CoulombOperator(pLattice lattice, pODESolver ode):
    OneDimensionalODE(lattice), fwd_direction(true), ode_solver(ode)
{   SetK(0);
}

void CoulombOperator::SetDensity(const RadialFunction& density)
{
    rho = density;
}

void CoulombOperator::GetPotential(unsigned int k, const RadialFunction& density, RadialFunction& pot, pODESolver ode)
{
    SetK(k);
    SetDensity(density);

    pODESolver ode_to_use;
    if(ode)
        ode_to_use = ode;
    else if (ode_solver)
        ode_to_use = ode_solver;
    else
        ode_to_use = pODESolver(new AdamsSolver(lattice));

    if(pot.Size() < density.Size())
        pot.ReSize(density.Size());

    // Integrate forwards to obtain I1
    fwd_direction = true;
    ode_to_use->IntegrateForwards(this, &pot);

    // Integrate backwards to obtain I2
    fwd_direction = false;
    RadialFunction I2(pot.Size());

    ode_to_use->IntegrateBackwards(this, &I2);

    pot += I2;
}

void CoulombOperator::GetODEFunction(unsigned int latticepoint, const RadialFunction& f, double* w) const
{
    double r = lattice->R(latticepoint);
    double density = 0.;
    if(latticepoint < rho.Size())
        density = rho.f[latticepoint];

    if(fwd_direction)
        *w = - double(k + 1)/r * f.f[latticepoint] + density/r;
    else
        *w = double(k)/r * f.f[latticepoint] - density/r;
}

void CoulombOperator::GetODECoefficients(unsigned int latticepoint, const RadialFunction& f, double* w_f, double* w_const) const
{
    double r = lattice->R(latticepoint);
    double density = 0.;
    if(latticepoint < rho.Size())
        density = rho.f[latticepoint];

    if(fwd_direction)
    {   *w_f = - double(k + 1)/r;
        *w_const = density/r;
    }
    else
    {   *w_f = double(k)/r;
        *w_const = - density/r;
    }
}

void CoulombOperator::GetODEJacobian(unsigned int latticepoint, const RadialFunction& f, double* jacobian, double* dwdr) const
{
    // dI1/dr = -(k+1)/r .I1 + rho/r = w1
    // dI2/dr =    k/r .I2   - rho/r = w2

    // dw1/dr = (k+1)(k+2)/r^2 .I1 - (k+2)/r^2 .rho + 1/r drho/dr
    // dw2/dr = k(k-1)/r^2 .I2 - (k-1)/r^2 .rho - 1/r drho/dr

    double r = lattice->R(latticepoint);
    double density = 0.;
    double drhodr = 0.;
    if(latticepoint < rho.Size())
    {   density = rho.f[latticepoint];
        drhodr = rho.dfdr[latticepoint];
    }

    if(fwd_direction)
    {   *jacobian = - double(k + 1)/r;
        *dwdr = double(k + 1)*(k + 2)/(r*r) * f.f[latticepoint]
                - double(k + 2)/(r*r) * density
                + drhodr/r;
    }
    else
    {   *jacobian = double(k)/r;
        *dwdr = double(k)*(k - 1)/(r*r) * f.f[latticepoint]
                - double(k - 1)/(r*r) * density
                - drhodr/r;
    }
}

void CoulombOperator::EstimateSolutionNearOrigin(unsigned int numpoints, RadialFunction& f) const
{
    // Assume fwd_direction == true and estimate I1.
    double integrand = 0.0;
    double w;
    if(f.Size() < numpoints)
        f.ReSize(numpoints);

    // Use trapezoidal rule for integrand
    for(unsigned int i = 0; i < numpoints; i++)
    {
        integrand += 0.5 * pow(lattice->R(i), k) * rho.f[i] * lattice->dR(i);
        f.f[i] = 1./pow(lattice->R(i), k+1) * integrand;

        GetODEFunction(i, f, &w);
        f.dfdr[i] = w;
        integrand +=  0.5 * pow(lattice->R(i), k) * rho.f[i] * lattice->dR(i);
    }
}

void CoulombOperator::EstimateSolutionNearInfinity(unsigned int numpoints, RadialFunction& f) const
{
    // Assume fwd_direction == false and estimate I2
    double integrand = 0.0;
    double w;

    for(unsigned int i = f.Size() - 1; i >= f.Size() - numpoints; i--)
    {
        integrand += 0.5 * 1./pow(lattice->R(i), k+1) * rho.f[i] * lattice->dR(i);
        f.f[i] = pow(lattice->R(i), k) * integrand;

        GetODEFunction(i, f, &w);
        f.dfdr[i] = (-double(k+1) * f.f[i] + rho.f[i])/lattice->R(i);
        integrand += 0.5 * 1./pow(lattice->R(i), k+1) * rho.f[i] * lattice->dR(i);
    }
}
