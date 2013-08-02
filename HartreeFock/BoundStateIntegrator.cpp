#include "Include.h"
#include "BoundStateIntegrator.h"
#include "Universal/Eigensolver.h"
#include "Universal/PhysicalConstant.h"

void BoundStateIntegrator::SetUpForwardsIntegral(SingleParticleWavefunction& s, const std::vector<double>& HFPotential, double nuclear_charge = 1.)
{
    unsigned int i;

    double correction = s.f[0];
    if(s.Kappa() < 0)
    {   s.f[0] = pow(lattice->R(0), -s.Kappa());
        s.g[0] = s.f[0] * lattice->R(0) * HFPotential[0] / (2 * s.Kappa() - 1);
    }
    else
    {   s.g[0] = pow(lattice->R(0), s.Kappa());
        s.f[0] = s.g[0] * lattice->R(0) * PhysicalConstant::Instance()->GetAlphaSquared() * HFPotential[0] / (2 * s.Kappa() + 1);
    }

    // Determine an appropriate scaling to make the norm close to unit.
    if(s.f[0] && correction)
        correction = correction/s.f[0];
    else
        correction = nuclear_charge * nuclear_charge;

    SolveDiracBoundary(s, HFPotential, 0);

    for(i=0; i<adams_N; i++)
    {   s.f[i] = s.f[i] * correction;
        s.g[i] = s.g[i] * correction;
        s.df[i] = s.df[i] * correction;
        s.dg[i] = s.dg[i] * correction;
    }
}

void BoundStateIntegrator::SetUpBackwardsIntegral(SingleParticleWavefunction& s, const std::vector<double>& HFPotential)
{
    // Get start point
    unsigned int start_point = s.Size() - 1;
    unsigned int i = start_point - (adams_N-2);
    double P;
    while(i < HFPotential.size())
    {   P = -2.*(HFPotential[i] + s.Energy()) + double(s.Kappa()*(s.Kappa() + 1))/pow(lattice->R(i),2.);
        if(P > 0.)
            break;
        i++;
    }

    start_point = i + (adams_N-2);
    s.ReSize(start_point+1);

    double correction = s.f[start_point];
    double S = -9.;

    P = -2*(HFPotential[start_point] + s.Energy()) + s.Kappa()*(s.Kappa() + 1)/pow(lattice->R(start_point),2.);
    //assert(P>0);
    P = sqrt(P);
    S = S + 0.5 * P * lattice->dR(start_point);

    s.f[start_point] = exp(S)/sqrt(P);
    s.g[start_point] = s.f[start_point] * (s.Kappa()/lattice->R(start_point) - P) * 0.5;

    SolveDiracBoundary(s, HFPotential, start_point, false);

    if(correction)
    {   correction = correction/s.f[start_point];
        for(unsigned int i=start_point; i>start_point-(adams_N-1); i--)
        {   s.f[i] = s.f[i] * correction;
            s.g[i] = s.g[i] * correction;
            s.df[i] = s.df[i] * correction;
            s.dg[i] = s.dg[i] * correction;
        }
    }
}

void BoundStateIntegrator::SolveDiracBoundary(SingleParticleWavefunction& s, const std::vector<double>& HFPotential, unsigned int start_point, bool forwards)
{
    unsigned int i, j;

    const double f0 = s.f[start_point];
    const double g0 = s.g[start_point];

    // Get first DM points using exact Adam's method calculation
    // coeff[i][j]/coeff_denom gives the coefficient of y_j at x = x_i
    unsigned int DM = 9;    // Num points - 1
    double coeff_denom = 362880;   // DM!

    double coeff[5][10]
    = {{-1026576, 3265920, -6531840, 10160640, -11430720,  9144576, -5080320, 1866240, -408240, 40320},
       {  -40320, -623376,  1451520, -1693440,   1693440, -1270080,   677376, -241920,   51840, -5040},
       {    5040,  -90720,  -396576,   846720,   -635040,   423360,  -211680,   72576,  -15120,  1440},
       {   -1440,   19440,  -155520,  -223776,    544320,  -272160,   120960,  -38880,    7776,  -720},
       {     720,   -8640,    51840,  -241920,    -72576,   362880,  -120960,   34560,   -6480,   576}};
    
    //double coeff[3][6] = {{-274.,  600., -600.,  400., -150., 24.},
    //                      { -24., -130.,  240., -120.,   40., -6.},
    //                      {   6.,  -60.,  -40.,  120.,  -30.,  4.}};

    double* A = new double[4*DM*DM];
    double* B = new double[2*DM];

    for(i=0; i<4*DM*DM; i++)
        A[i] = 0.;

    // Set up upper left quarter
    for(i=0; i<DM; i++)
    {   for(j=0; j<DM; j++)
        {   if(i < DM/2)
                A[i*2*DM + j] = coeff[i+1][j+1]/coeff_denom;
            else
                A[i*2*DM + j] = -coeff[DM-i-1][DM-j-1]/coeff_denom;
            
            if(i == j)
            {   if(forwards)
                    A[i*2*DM + j] += s.Kappa()/lattice->R(start_point + i+1) * lattice->dR(start_point + i+1);
                else
                    A[i*2*DM + j] -= s.Kappa()/lattice->R(start_point - (i+1)) * lattice->dR(start_point - (i+1));
            }
        }
    }
    
    // Set up upper right quarter
    for(i=0; i<DM; i++)
    {
        if(forwards)
            A[i*2*DM + i+DM] = -(2. + PhysicalConstant::Instance()->GetAlphaSquared() * HFPotential[start_point + i+1])*lattice->dR(start_point + i+1);
        else
            A[i*2*DM + i+DM] = (2. + PhysicalConstant::Instance()->GetAlphaSquared() * HFPotential[start_point -(i+1)])*lattice->dR(start_point - (i+1));
    }

    // Set up lower left quarter
    for(i=0; i<DM; i++)
    {
        if(forwards)
            A[(i+DM)*2*DM + i] = HFPotential[start_point + i+1]*lattice->dR(start_point + i+1);
        else
            A[(i+DM)*2*DM + i] = -HFPotential[start_point - (i+1)]*lattice->dR(start_point - (i+1));
    }

    // Set up lower right quarter
    for(i=0; i<DM; i++)
    {   for(j=0; j<DM; j++)
        {   if(i < DM/2)
                A[(i+DM)*2*DM + j+DM] = coeff[i+1][j+1]/coeff_denom;
            else
                A[(i+DM)*2*DM + j+DM] = -coeff[DM-i-1][DM-j-1]/coeff_denom;
            
            if(i == j)
            {   if(forwards)
                    A[(i+DM)*2*DM + j+DM] -= s.Kappa()/lattice->R(start_point + i+1) * lattice->dR(start_point + i+1);
                else
                    A[(i+DM)*2*DM + j+DM] += s.Kappa()/lattice->R(start_point - (i+1)) * lattice->dR(start_point - (i+1));
            }
        }
    }

    // Set up B
    for(i=0; i<DM; i++)
    {
        if(i < DM/2)
        {   B[i] = -coeff[i+1][0]/coeff_denom * f0;
            B[i+DM] = -coeff[i+1][0]/coeff_denom * g0;
        }
        else
        {   B[i] = coeff[DM-i-1][DM]/coeff_denom * f0;
            B[i+DM] = coeff[DM-i-1][DM]/coeff_denom * g0;
        }
    }

    // Solve A*x = B
    Eigensolver E;
    E.SolveSimultaneousEquations(A, B, 2*DM);

    if(forwards)
        for(i=0; i<DM; i++)
        {   s.f[start_point + i+1] = B[i];
            s.g[start_point + i+1] = B[DM+i];
        }
    else
        for(i=0; i<DM; i++)
        {   s.f[start_point - (i+1)] = B[i];
            s.g[start_point - (i+1)] = B[DM+i];
        }

    // Calculate derivatives
    if(forwards)
        for(i=start_point; i<start_point+DM; i++)
        {   s.df[i] = (-s.Kappa()/lattice->R(i)*s.f[i] + (2. + PhysicalConstant::Instance()->GetAlphaSquared() * HFPotential[i])*s.g[i])*lattice->dR(i);
            s.dg[i] = (-HFPotential[i] * s.f[i] + s.Kappa()/lattice->R(i)*s.g[i])*lattice->dR(i);
        }
    else
        for(i=start_point; i>start_point-DM; i--)
        {   s.df[i] = (-s.Kappa()/lattice->R(i)*s.f[i] + (2. + PhysicalConstant::Instance()->GetAlphaSquared() * HFPotential[i])*s.g[i])*lattice->dR(i);
            s.dg[i] = (-HFPotential[i] * s.f[i] + s.Kappa()/lattice->R(i)*s.g[i])*lattice->dR(i);
        }

    // Calculate other points up to adams_N
    if(adams_N > DM+1)
    {   unsigned int old_N = adams_N;
        SetAdamsOrder(5);

        CoupledFunction exchange(s.Size());
        StateFunction Af(lattice);
        Af.SetHFPotential(HFPotential);
        Af.SetState(&s);
        Af.SetExchange(&exchange);

        if(forwards)
            Integrate2(Af, s, start_point+DM, start_point+old_N);
        else
            Integrate2(Af, s, start_point-DM, start_point-old_N);

        SetAdamsOrder(old_N);
    }

    delete[] A;
    delete[] B;
}
