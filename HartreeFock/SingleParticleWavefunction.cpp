#include "Include.h"
#include "SingleParticleWavefunction.h"
#include "Universal/Constant.h"

SingleParticleWavefunction::SingleParticleWavefunction():
    CoupledFunction(), nu(0.), kappa(0)
{}

SingleParticleWavefunction::SingleParticleWavefunction(int Kappa):
    CoupledFunction(), nu(0.), kappa(Kappa)
{}

SingleParticleWavefunction::SingleParticleWavefunction(const SingleParticleWavefunction& other):
    CoupledFunction(other), kappa(other.kappa), nu(other.nu)
{}

void SingleParticleWavefunction::Read(FILE* fp)
{
    CoupledFunction::Read(fp);
}

double SingleParticleWavefunction::Overlap(const SingleParticleWavefunction& other, const Lattice* lattice) const
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

bool SingleParticleWavefunction::Print(Lattice* lattice) const
{
    return Print(Name() + "_orbital.txt", lattice);
}

bool SingleParticleWavefunction::Print(const std::string& filename, Lattice* lattice) const
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

bool SingleParticleWavefunction::Print(FILE* fp, Lattice* lattice) const
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
