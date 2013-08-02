#include "GreensIntegrator.h"
#include "Universal/PhysicalConstant.h"
#include "Include.h"

std::vector<double> GreensIntegrator::GetGreensInfinity()
{
    std::vector<double> Ginf;

    if(!G || !s_0 || !s_inf)
    {   return Ginf;
    }

    unsigned int limit = mmax(s_inf->Size(), G->Size());
    Ginf.resize(limit);

    // Make integrand
    unsigned int i;
    double W;
    double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();

    std::vector<double> integrand(limit);
    for(i = 0; i<mmin(s_inf->Size(), G->Size()); i++)
    {   W = s_0->f[i] * s_inf->g[i] - s_inf->f[i] * s_0->g[i];
        integrand[i] = - (G->f[i] * s_inf->f[i] + alphasquared * G->g[i] * s_inf->g[i])/W;
    }

    const double* dR = lattice->dR();

    i = limit-1;
    Ginf[i] = - 0.5 * integrand[i] * dR[i];

    for(i = limit-2; i > limit - adams_N; i--)
    {   Ginf[i] = Ginf[i+1] - 0.5 * integrand[i+1] * dR[i+1]
                            - 0.5 * integrand[i] * dR[i];
    }

    Integrate0(integrand, Ginf, limit - adams_N, -1);

    return Ginf;
}

std::vector<double> GreensIntegrator::GetGreensOrigin()
{
    std::vector<double> G0;

    if(!G || !s_0 || !s_inf)
    {   return G0;
    }

    unsigned int limit = mmax(s_0->Size(), G->Size());
    G0.resize(limit);
    double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();

    // Make integrand
    unsigned int i;
    double W;
    std::vector<double> integrand(limit);
    for(i = 0; i<mmin(s_0->Size(), G->Size()); i++)
    {   W = s_0->f[i] * s_inf->g[i] - s_inf->f[i] * s_0->g[i];
        integrand[i] = (G->f[i] * s_0->f[i] + alphasquared * G->g[i] * s_0->g[i])/W;
    }

    const double* dR = lattice->dR();

    i = 0;
    G0[i] = 0.5 * integrand[i] * dR[i];

    for(i = 1; i < adams_N - 1; i++)
    {   G0[i] = G0[i-1] + 0.5 * integrand[i-1] * dR[i-1]
                        + 0.5 * integrand[i] * dR[i];
    }

    Integrate0(integrand, G0, adams_N-1, limit);

    return G0;
}
