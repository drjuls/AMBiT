#include "Include.h"
#include "State.h"
#include "Universal/Constant.h"

State::State():
    CoupledFunction(), nu(0.), kappa(0)
{}

State::State(int Kappa):
    CoupledFunction(), nu(0.), kappa(Kappa)
{}

State::State(const State& other):
    CoupledFunction(other), kappa(other.kappa), nu(other.nu)
{}

unsigned int State::NumZeroes() const
{
    // Count number of zeroes
    unsigned int zeroes = 0;

    // Get maximum point
    double fmax = 0.;
    unsigned int i_fmax = 0;
    unsigned int i = 0;
    while(i < f.size())
    {   if(fabs(f[i]) > fmax)
        {   fmax = fabs(f[i]);
            i_fmax = i;
        }
        i++;
    }
    
    // Get effective end point
    unsigned int end = f.size() - 1;
    while(fabs(f[end]) < 1.e-7 * fmax)
        end--;

    i = 0;
    double prev = f[i];
    while(fabs(prev) < 1.e-7 * fmax)
        prev = f[i++];
    while(i < end)
    {   if(f[i] * prev < 0.)
            zeroes++;
        prev = f[i];
        i++;
    }
    return zeroes;
}

void State::Read(FILE* fp)
{
    CoupledFunction::Read(fp);
}

double State::Overlap(const State& other, const Lattice* lattice) const
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

bool State::Print(Lattice* lattice) const
{
    return Print(Name() + "_orbital.txt", lattice);
}

bool State::Print(const std::string& filename, Lattice* lattice) const
{
    FILE* fp = fopen(filename.c_str(), "wt");

    if(fp)
    {   bool success = Print(fp, lattice);
        fclose(fp);
        return success;
    }
    else
        return false;
}

bool State::Print(FILE* fp, Lattice* lattice) const
{
    if(lattice)
        for(unsigned int i = 0; i < Size(); i++)
        {
            fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e\n", lattice->R(i),
                f[i], g[i]*Constant::Alpha, df[i], dg[i]*Constant::Alpha);
        }
    else
        for(unsigned int i = 0; i < Size(); i++)
        {
            fprintf(fp, "%d %12.6e %12.6e %12.6e %12.6e\n", i,
                f[i], g[i]*Constant::Alpha, df[i], dg[i]*Constant::Alpha);
        }
    
    return true;
}
