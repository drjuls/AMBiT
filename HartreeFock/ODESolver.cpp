#include "ODESolver.h"
#include "Include.h"

AdamsSolver::AdamsSolver(pOPIntegrator integrator): ODESolver(integrator)
{
    order = 10;
    adams_coeff
        = {2082753./7257600., 9449717./7257600., -11271304./7257600., 16002320./7257600., -17283646./7257600.,
           13510082./7257600., -7394032./7257600., 2687864./7257600., -583435./7257600., 57281./7257600.};
}

void AdamsSolver::IntegrateForwards(const OneDimensionalODE* op, RadialFunction* solution)
{
    RadialFunction& s = *solution;
    op->EstimateSolutionNearOrigin(order, s);

    int start_point = order;
    int end_point = s.size();

    double w_f, w_const;

    const double* dR = lattice->dR();

    for(int i = start_point; i < end_point; i++)
    {
        op->GetODECoefficients(i, s, &w_f, &w_const);
        
        double f_next = s.f[i-1] + adams_coeff[0] * w_const * dR[i];
        
        for(unsigned int j=1; j<order; j++)
        {
            f_next = f_next + adams_coeff[j] * s.dfdr[i-j] * dR[i-j];
        }
        
        double D = 1. - adams_coeff[0] * w_f * dR[i];
        s.f[i] = f_next / D;
        s.dfdr[i] = w_f * s.f[i] + w_const;
    }
}

void AdamsSolver::IntegrateBackwards(const OneDimensionalODE* op, RadialFunction* solution)
{
    RadialFunction& s = *solution;
    op->EstimateSolutionNearInfinity(order, s);
    
    int start_point = s.size() - order;
    int end_point = 0;
    
    double w_f, w_const;
    
    const double* dR = lattice->dR();
    
    for(int i = start_point; i >= end_point; i--)
    {
        op->GetODECoefficients(i, s, &w_f, &w_const);
        
        double f_next = s.f[i+1] - adams_coeff[0] * w_const * dR[i];
        
        for(unsigned int j=1; j<order; j++)
        {
            f_next = f_next - adams_coeff[j] * s.dfdr[i+j] * dR[i+j];
        }
        
        double D = 1. + adams_coeff[0] * w_f * dR[i];
        s.f[i] = f_next / D;
        s.dfdr[i] = w_f * s.f[i] + w_const;
    }
}

void AdamsSolver::IntegrateForwards(pSpinorODEConst op, SpinorFunction* solution)
{
    SpinorFunction& s = *solution;
    op->EstimateOrbitalNearOrigin(order, s);

    int start_point = order;
    int end_point = s.size();

    double w_f[2];
    double w_g[2];
    double w_const[2];

    const double* dR = lattice->dR();

    for(int i = start_point; i < end_point; i++)
    {
        op->GetODECoefficients(i, s, w_f, w_g, w_const);

        double f_next = s.f[i-1] + adams_coeff[0] * w_const[0] * dR[i];
        double g_next = s.g[i-1] + adams_coeff[0] * w_const[1] * dR[i];

        for(unsigned int j=1; j<order; j++)
        {
            f_next = f_next + adams_coeff[j] * s.dfdr[i-j] * dR[i-j];
            g_next = g_next + adams_coeff[j] * s.dgdr[i-j] * dR[i-j];
        }
        
        double D = (1. - adams_coeff[0] * w_f[0] * dR[i])
                    * (1. - adams_coeff[0] * w_g[1]  * dR[i])
                  -(adams_coeff[0] * w_g[0] * dR[i])
                    * (adams_coeff[0] * w_f[1] * dR[i]);
        
        s.f[i] = (f_next * (1. - adams_coeff[0] * w_g[1] * dR[i])
                  + g_next * (adams_coeff[0] * w_g[0] * dR[i])) / D;
        s.g[i] = (g_next * (1. - adams_coeff[0] * w_f[0] * dR[i])
                  + f_next * (adams_coeff[0] * w_f[1] * dR[i])) / D;
        s.dfdr[i] = w_f[0] * s.f[i] + w_g[0] * s.g[i] + w_const[0];
        s.dgdr[i] = w_f[1] * s.f[i] + w_g[1] * s.g[i] + w_const[1];
    }
}

void AdamsSolver::IntegrateBackwards(pSpinorODEConst op, Orbital* solution)
{
    Orbital& s = *solution;
    op->EstimateOrbitalNearInfinity(order, s);

    int start_point = s.size()-order;
    int end_point = 0;
    
    double w_f[2];
    double w_g[2];
    double w_const[2];

    const double* dR = lattice->dR();

    for(int i = start_point; i >= end_point; i--)
    {
        op->GetODECoefficients(i, s, w_f, w_g, w_const);

        double f_next = s.f[i+1] - adams_coeff[0] * w_const[0] * dR[i];
        double g_next = s.g[i+1] - adams_coeff[0] * w_const[1] * dR[i];
        
        for(unsigned int j=1; j<order; j++)
        {
            f_next = f_next - adams_coeff[j] * s.dfdr[i+j] * dR[i+j];
            g_next = g_next - adams_coeff[j] * s.dgdr[i+j] * dR[i+j];
        }
        
        double D = (1. + adams_coeff[0] * w_f[0] * dR[i])
                    * (1. + adams_coeff[0] * w_g[1] * dR[i])
                  -(adams_coeff[0] * w_g[0] * dR[i])
                    * (adams_coeff[0] * w_f[1] * dR[i]);
        
        s.f[i] = (f_next * (1. + adams_coeff[0] * w_g[1] * dR[i])
                  - g_next * (adams_coeff[0] * w_g[0] * dR[i])) / D;
        s.g[i] = (g_next * (1. + adams_coeff[0] * w_f[0] * dR[i])
                  - f_next * (adams_coeff[0] * w_f[1] * dR[i])) / D;
        s.dfdr[i] = w_f[0] * s.f[i] + w_g[0] * s.g[i] + w_const[0];
        s.dgdr[i] = w_f[1] * s.f[i] + w_g[1] * s.g[i] + w_const[1];
    }
}

unsigned int AdamsSolver::IntegrateBackwardsUntilPeak(pSpinorODEConst op, Orbital* solution, int classical_turning_point)
{
    Orbital& s = *solution;
    op->EstimateOrbitalNearInfinity(order, s);

    int start_point = s.size()-order;
    int end_point = 0;
    unsigned int peak = 0;

    double w_f[2];
    double w_g[2];
    double w_const[2];

    const double* dR = lattice->dR();

    for(int i = start_point; i >= end_point; i--)
    {
        op->GetODECoefficients(i, s, w_f, w_g, w_const);
        
        double f_next = s.f[i+1] - adams_coeff[0] * w_const[0] * dR[i];
        double g_next = s.g[i+1] - adams_coeff[0] * w_const[1] * dR[i];
        
        for(unsigned int j=1; j<order; j++)
        {
            f_next = f_next - adams_coeff[j] * s.dfdr[i+j] * dR[i+j];
            g_next = g_next - adams_coeff[j] * s.dgdr[i+j] * dR[i+j];
        }
        
        double D = (1. + adams_coeff[0] * w_f[0] * dR[i])
                    * (1. + adams_coeff[0] * w_g[1] * dR[i])
                  -(adams_coeff[0] * w_g[0] * dR[i])
                    * (adams_coeff[0] * w_f[1] * dR[i]);
        
        s.f[i] = (f_next * (1. + adams_coeff[0] * w_g[1] * dR[i])
                  - g_next * (adams_coeff[0] * w_g[0] * dR[i])) / D;
        s.g[i] = (g_next * (1. + adams_coeff[0] * w_f[0] * dR[i])
                  - f_next * (adams_coeff[0] * w_f[1] * dR[i])) / D;
        s.dfdr[i] = w_f[0] * s.f[i] + w_g[0] * s.g[i] + w_const[0];
        s.dgdr[i] = w_f[1] * s.f[i] + w_g[1] * s.g[i] + w_const[1];

        // Break when peak is reached
        if(s.dfdr[i]/s.dfdr[i+1] <= 0.0 && i <= classical_turning_point)
        {   peak = i;
            break;
        }
    }

    return peak;
}
