#include "Include.h"
#include "SingleParticleWavefunction.h"
#include "OrbitalInfo.h"
#include "Universal/MathConstant.h"

SingleParticleWavefunction::SingleParticleWavefunction(const OrbitalInfo& info):
    SpinorFunction(info.Kappa()), pqn(info.PQN()), energy(0.0)
{}

SingleParticleWavefunction::SingleParticleWavefunction(int kappa, unsigned int pqn, double energy, unsigned int size):
    SpinorFunction(kappa, size), pqn(pqn), energy(energy)
{}

SingleParticleWavefunction::SingleParticleWavefunction(const SingleParticleWavefunction& other):
    SpinorFunction(other), pqn(other.pqn), energy(other.energy)
{}

const SingleParticleWavefunction& SingleParticleWavefunction::operator=(const SingleParticleWavefunction& other)
{
    SpinorFunction::operator=(other);
    pqn = other.pqn;
    energy = other.energy;
    return *this;
}

SingleParticleWavefunction& SingleParticleWavefunction::operator=(SingleParticleWavefunction&& other)
{
    SpinorFunction::operator=(other);
    pqn = other.pqn;
    energy = other.energy;
    return *this;
}

double SingleParticleWavefunction::Energy() const
{
    return energy;
}

/** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
double SingleParticleWavefunction::Nu() const
{
    double nu;
    if(std::fabs(energy) < 1.e-10)
        nu = nan("energy is zero");
    else
    {   nu = std::sqrt(std::fabs(-0.5/energy));
        if(energy > 0.)
            nu = -nu;
    }

    return nu;
}

int SingleParticleWavefunction::PQN() const
{
    return pqn;
}

void SingleParticleWavefunction::SetEnergy(double Energy)
{
    energy = Energy;
}

void SingleParticleWavefunction::SetNu(double Nu)
{
    if(fabs(Nu) > 1.e-5)
        energy = -0.5/(Nu * Nu);
    else
        energy = -1.e10;

    if(Nu < 0)
        energy = -energy;
}

void SingleParticleWavefunction::SetPQN(int PQN)
{
    pqn = PQN;
}

double SingleParticleWavefunction::Overlap(const SingleParticleWavefunction &other, pLatticeConst lattice) const
{
    double total = 0.;
    
    if(kappa == other.kappa)
    {
        const double* dR = lattice->dR();

        for(unsigned int i=0; i<mmin(size(), other.size()); i++)
        {
            total += (f[i] * other.f[i] + g[i] * other.g[i])*dR[i];
        }
    }

    return total;
}

std::string SingleParticleWavefunction::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);
    
    ret.append(1, MathConstant::Instance()->GetSpectroscopicNotation(L()));
    
#ifdef USE_ALT_STATE_NOTATION
    if(kappa > 0)
        ret.append(1, '-');
    else
        ret.append(1, ' ');
#else
    if(kappa < -1)
        ret.append(1, '+');
#endif
    
    return ret;
}

const SingleParticleWavefunction& SingleParticleWavefunction::operator*=(double scale_factor)
{
    SpinorFunction::operator*=(scale_factor);
    return *this;
}

SingleParticleWavefunction SingleParticleWavefunction::operator*(double scale_factor) const
{
    SingleParticleWavefunction ret(*this);
    ret *= scale_factor;
    return ret;
}

const SingleParticleWavefunction& SingleParticleWavefunction::operator+=(const SingleParticleWavefunction& other)
{
    SpinorFunction::operator+=(other);
    return *this;
}

const SingleParticleWavefunction& SingleParticleWavefunction::operator-=(const SingleParticleWavefunction& other)
{
    SpinorFunction::operator-=(other);
    return *this;
}

SingleParticleWavefunction SingleParticleWavefunction::operator+(const SingleParticleWavefunction& other) const
{
    SingleParticleWavefunction ret(*this);
    ret += other;
    return ret;
}

SingleParticleWavefunction SingleParticleWavefunction::operator-(const SingleParticleWavefunction& other) const
{
    SingleParticleWavefunction ret(*this);
    ret -= other;
    return ret;    
}

const SingleParticleWavefunction& SingleParticleWavefunction::operator*=(const RadialFunction& chi)
{
    SpinorFunction::operator*=(chi);
    return *this;
}

SingleParticleWavefunction SingleParticleWavefunction::operator*(const RadialFunction& chi) const
{
    SingleParticleWavefunction ret(*this);
    ret *= chi;
    return ret;
}

void SingleParticleWavefunction::Write(FILE* fp) const
{
    // As well as the SpinorFunction vectors, we need to output some other things
    fwrite(&pqn, sizeof(int), 1, fp);
    fwrite(&energy, sizeof(double), 1, fp);

    SpinorFunction::Write(fp);
}

void SingleParticleWavefunction::Read(FILE* fp)
{
    fread(&pqn, sizeof(int), 1, fp);
    fread(&energy, sizeof(double), 1, fp);

    SpinorFunction::Read(fp);
}

bool SingleParticleWavefunction::Print(pLattice lattice) const
{
    return Print(Name() + "_orbital.txt", lattice);
}

bool SingleParticleWavefunction::Print(const std::string& filename, pLattice lattice) const
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

bool SingleParticleWavefunction::Print(FILE* fp, pLattice lattice) const
{
    unsigned int i;
    if(lattice)
        for(i = 0; i < size(); i++)
        {
            fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e\n", lattice->R(i),
                f[i], g[i], dfdr[i], dgdr[i]);
        }
    else
        for(i = 0; i < size(); i++)
        {
            fprintf(fp, "%d %12.6e %12.6e %12.6e %12.6e\n", i,
                f[i], g[i], dfdr[i], dgdr[i]);
        }
    
    return true;
}
