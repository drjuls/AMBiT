#include "Include.h"
#include "State.h"
#include "Universal/Constant.h"

State::State(Lattice* lat, unsigned int num_points):
    CoupledFunction(num_points), lattice(lat), nu(0.), kappa(0)
{}

State::State(Lattice* lat, int Kappa, unsigned int num_points):
    CoupledFunction(num_points), lattice(lat), nu(0.), kappa(Kappa)
{}

State::State(const State& other):
    CoupledFunction(other), kappa(other.kappa), nu(other.nu), lattice(other.lattice)
{}

unsigned int State::NumZeroes() const
{
    // Count number of zeroes
    unsigned int zeroes = 0;
    std::vector<double>::const_iterator i = f.begin();
    double prev = *i++;
    while(prev == 0.)
        prev = *i++;
    while(i != f.end())
    {   if((*i) * prev < 0.)
            zeroes++;
        prev = *i;
        i++;
    }
    return zeroes;
}

void State::Read(FILE* fp)
{
    CoupledFunction::Read(fp);
    lattice->R(Size()-1);
}

double State::Overlap(const State& other) const
{
    double total = 0.;

    if(kappa == other.kappa)
    {
        const double* dR = lattice->dR();

        for(unsigned int i=0; i<mmin(Size(), other.Size()); i++)
        {
            total += (f[i] * other.f[i] + Constant::AlphaSquared * g[i] * other.g[i])*dR[i];
        }
    }

    return total;
}
