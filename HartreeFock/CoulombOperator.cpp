#include "CoulombOperator.h"
#include <gsl/gsl_math.h>

CoulombOperator::CoulombOperator(pLattice lattice, pODESolver ode):
    OneDimensionalODE(lattice), fwd_direction(true), ode_solver(ode)
{   SetK(0);
}

void CoulombOperator::SetDensity(const RadialFunction& density)
{
    rho = density;
}

void CoulombOperator::GetPotential(int k, const RadialFunction& density, RadialFunction& pot, pODESolver ode)
{
    SetK(k);
    SetDensity(density);

    pODESolver ode_to_use;
    if(ode)
        ode_to_use = ode;
    else if (ode_solver)
        ode_to_use = ode_solver;
    else
    {   pIntegrator integrator(new SimpsonsIntegrator(lattice));
        ode_to_use = pODESolver(new AdamsSolver(integrator));
    }

    if(pot.size() < density.size())
        pot.resize(density.size());

    // Integrate forwards to obtain I1
    fwd_direction = true;
    ode_to_use->IntegrateForwards(this, &pot);

    // Integrate backwards to obtain I2
    fwd_direction = false;
    RadialFunction I2(pot.size());

    ode_to_use->IntegrateBackwards(this, &I2);

    pot += I2;
}

/** Get zero-multipole potential, but renormalise density so that potential function goes as charge/r at infinity. */
void CoulombOperator::GetPotential(RadialFunction& density, RadialFunction& pot, double charge, pODESolver ode)
{
    SetK(0);
    SetDensity(density);

    pODESolver ode_to_use;
    if(ode)
        ode_to_use = ode;
    else if (ode_solver)
        ode_to_use = ode_solver;
    else
    {   pIntegrator integrator(new SimpsonsIntegrator(lattice));
        ode_to_use = pODESolver(new AdamsSolver(integrator));
    }

    if(pot.size() < density.size())
        pot.resize(density.size());

    // Integrate forwards to obtain I1
    fwd_direction = true;
    ode_to_use->IntegrateForwards(this, &pot);

    // Renormalise.
    // To generalise this function use potential = charge/R^(k+1) here rather than just R;
    if(fabs(charge) > 1.e-6)
    {
        double norm = pot.f[pot.size()-1] * lattice->R(pot.size()-1);
        norm = charge/norm;

        pot *= norm;
        density *= norm;
    }

    // Integrate backwards to obtain I2
    fwd_direction = false;
    RadialFunction I2(pot.size());

    ode_to_use->IntegrateBackwards(this, &I2);
    
    pot += I2;
}

void CoulombOperator::GetForwardPotential(int k, const RadialFunction& density, RadialFunction& pot, pODESolver ode)
{
    SetK(k);
    SetDensity(density);

    pODESolver ode_to_use;
    if(ode)
        ode_to_use = ode;
    else if (ode_solver)
        ode_to_use = ode_solver;
    else
    {   pIntegrator integrator(new SimpsonsIntegrator(lattice));
        ode_to_use = pODESolver(new AdamsSolver(integrator));
    }

    if(pot.size() < density.size())
        pot.resize(density.size());

    // Integrate forwards to obtain I1
    fwd_direction = true;
    ode_to_use->IntegrateForwards(this, &pot);
}

void CoulombOperator::GetBackwardPotential(int k, const RadialFunction& density, RadialFunction& pot, pODESolver ode)
{
    SetK(k);
    SetDensity(density);

    pODESolver ode_to_use;
    if(ode)
        ode_to_use = ode;
    else if (ode_solver)
        ode_to_use = ode_solver;
    else
    {   pIntegrator integrator(new SimpsonsIntegrator(lattice));
        ode_to_use = pODESolver(new AdamsSolver(integrator));
    }

    if(pot.size() < density.size())
        pot.resize(density.size());

    // Integrate backwards to obtain I2
    fwd_direction = false;
    ode_to_use->IntegrateBackwards(this, &pot);
}

void CoulombOperator::GetODEFunction(unsigned int latticepoint, const RadialFunction& f, double* w) const
{
    double r = lattice->R(latticepoint);
    double density = 0.;
    if(latticepoint < rho.size())
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
    if(latticepoint < rho.size())
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
    if(latticepoint < rho.size())
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
    if(f.size() < numpoints)
        f.resize(numpoints);

    // Use trapezoidal rule for integrand
    for(unsigned int i = 0; i < numpoints; i++)
    {
        double rhof = 0.;
        if(i < rho.size())
            rhof = rho.f[i];
        integrand += 0.5 * gsl_pow_int(lattice->R(i), k) * rhof * lattice->dR(i);
        f.f[i] = 1./gsl_pow_int(lattice->R(i), k+1) * integrand;

        GetODEFunction(i, f, &w);
        f.dfdr[i] = w;
        integrand +=  0.5 * gsl_pow_int(lattice->R(i), k) * rhof * lattice->dR(i);
    }
}

void CoulombOperator::EstimateSolutionNearInfinity(unsigned int numpoints, RadialFunction& f) const
{
    // Assume fwd_direction == false and estimate I2
    double integrand = 0.0;
    double w;

    for(unsigned int i = f.size() - 1; i >= f.size() - numpoints; i--)
    {
        double rhof = 0.;
        if(i < rho.size())
            rhof = rho.f[i];
        integrand += 0.5 * 1./gsl_pow_int(lattice->R(i), k+1) * rhof * lattice->dR(i);
        f.f[i] = gsl_pow_int(lattice->R(i), k) * integrand;

        GetODEFunction(i, f, &w);
        f.dfdr[i] = (-double(k+1) * f.f[i] + rhof)/lattice->R(i);
        integrand += 0.5 * 1./gsl_pow_int(lattice->R(i), k+1) * rhof * lattice->dR(i);
    }
}
