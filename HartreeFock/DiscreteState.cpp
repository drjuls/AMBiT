#include "Include.h"
#include "DiscreteState.h"
#include "Universal/Constant.h"
#include <math.h>

DiscreteState::DiscreteState(Lattice* lat, unsigned int num_points): 
    State(lat, num_points)
{}

DiscreteState::DiscreteState(Lattice* lat, unsigned int PrincipalQN, int Kappa, unsigned int num_points):
    State(lat, Kappa, num_points), pqn(PrincipalQN)
{
    occupancy = 2*abs(kappa);
}

DiscreteState::DiscreteState(const DiscreteState& other):
    State(other), pqn(other.pqn), occupancy(other.occupancy)
{}

double DiscreteState::Energy() const
{   if(nu > 0.)
        return -1./(2.*nu*nu);
    else if(nu < 0.)
        return 1./(2.*nu*nu);
    else
        return 0.;
}

double DiscreteState::Norm() const
{
    double norm = 0.;
    unsigned int i;
    const double* dR = lattice->dR();

    norm = f[0]*f[0] * Constant::AlphaSquared * g[0]*g[0] * dR[0];
    for(i=0; i<Size()-2; i+=2)
    {   norm = norm + 4. * (f[i]*f[i] + Constant::AlphaSquared * g[i]*g[i]) * dR[i];
        norm = norm + 2. * (f[i+1]*f[i+1] + Constant::AlphaSquared * g[i+1]*g[i+1]) * dR[i+1];
    }

    norm = norm/3.;

    while(i<Size())
    {   norm = norm + (f[i]*f[i] + Constant::AlphaSquared * g[i]*g[i]) * dR[i];
        i++;
    }

    return norm;
}

std::string DiscreteState::Name() const
{
    unsigned int level;
    level = RequiredPQN();

    char buffer[20];
    sprintf(buffer, "%d", level);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    if(kappa < -1)
        ret.append(1, '+');

    return ret;
}

bool DiscreteState::CheckSize(double tolerance)
{
    double maximum = 0.;
    unsigned int i = 0;
    unsigned int max_point = 0;
    while(i < f.size())
    {   if(maximum < fabs(f[i]))
        {   maximum = fabs(f[i]);
            max_point = i;
        }
        i++;
    }
    
    i = max_point;
    while((i<f.size()) && (fabs(f[i])/maximum > tolerance))
        i++;
    
    if(i < f.size())
    {
        if(i == f.size() - 1)
            return true;

        ReSize(i+1);
        return false;
    }
    else
    {   // add points to wavefunction
        unsigned int max = f.size()-1;
        unsigned int old_size = f.size();
        double f_max = fabs(f[max]);
        double f_ratio = f[max]/f[max-1],
               g_ratio = g[max]/g[max-1];

        if(f_ratio > 1.)
            return true;

        while(f_max/maximum >= tolerance)
        {   max++;
            f_max = f_max * f_ratio;
        }

        ReSize(max+1);
        for(unsigned int i=old_size; i<max+1; i++)
        {   f[i] = f[i-1]*f_ratio;
            g[i] = g[i-1]*g_ratio;
        }

        // Check lattice is okay
        lattice->R(max);
        return false;
    }
}

void DiscreteState::ReNormalise(double norm)
{
    double scaling = sqrt(norm/Norm());
    Scale(scaling);
}

void DiscreteState::Write(FILE* fp) const
{
    // As well as the CoupledFunction vectors, we need to output some other things
    fwrite(&kappa, sizeof(int), 1, fp);
    fwrite(&nu, sizeof(double), 1, fp);
    fwrite(&pqn, sizeof(unsigned int), 1, fp);
    fwrite(&occupancy, sizeof(double), 1, fp);

    State::Write(fp);
}

void DiscreteState::Read(FILE* fp)
{
    fread(&kappa, sizeof(int), 1, fp);
    fread(&nu, sizeof(double), 1, fp);
    fread(&pqn, sizeof(unsigned int), 1, fp);
    fread(&occupancy, sizeof(double), 1, fp);

    State::Read(fp);
}
