#include "Include.h"
#include "Orbital.h"
#include "Universal/PhysicalConstant.h"
#include "OrbitalInfo.h"
#include <math.h>

Orbital::Orbital(const OrbitalInfo& info):
    SpinorFunction(info.Kappa()), pqn(info.PQN()), energy(0.0)
{}

Orbital::Orbital(int kappa, int pqn, double energy, unsigned int size):
    SpinorFunction(kappa, size), pqn(pqn), energy(energy)
{}

Orbital::Orbital(const Orbital& other):
    SpinorFunction(other), pqn(other.pqn), energy(other.energy)
{}

Orbital::Orbital(Orbital&& other):
    SpinorFunction(other), pqn(other.pqn), energy(other.energy)
{}

const Orbital& Orbital::operator=(const Orbital& other)
{
    SpinorFunction::operator=(other);
    return *this;
}

Orbital& Orbital::operator=(Orbital&& other)
{
    SpinorFunction::operator=(other);
    return *this;
}

const Orbital& Orbital::operator*=(double scale_factor)
{
    SpinorFunction::operator*=(scale_factor);
    return *this;
}

Orbital Orbital::operator*(double scale_factor) const
{
    Orbital ret(*this);
    ret *= scale_factor;
    return ret;
}

const Orbital& Orbital::operator+=(const Orbital& other)
{
    SpinorFunction::operator+=(other);
    return *this;
}

const Orbital& Orbital::operator-=(const Orbital& other)
{
    SpinorFunction::operator-=(other);
    return *this;
}

Orbital Orbital::operator+(const Orbital& other) const
{
    Orbital ret(*this);
    ret += other;
    return ret;
}

Orbital Orbital::operator-(const Orbital& other) const
{
    Orbital ret(*this);
    ret -= other;
    return ret;
}

const Orbital& Orbital::operator*=(const RadialFunction& chi)
{
    SpinorFunction::operator*=(chi);
    return *this;
}

Orbital Orbital::operator*(const RadialFunction& chi) const
{
    Orbital ret(*this);
    ret *= chi;
    return ret;
}

double Orbital::Energy() const
{
    return energy;
}

/** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
double Orbital::Nu() const
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

int Orbital::PQN() const
{
    return pqn;
}

void Orbital::SetEnergy(double Energy)
{
    energy = Energy;
}

void Orbital::SetNu(double Nu)
{
    if(fabs(Nu) > 1.e-5)
        energy = -0.5/(Nu * Nu);
    else
        energy = -1.e10;
    
    if(Nu < 0)
        energy = -energy;
}

void Orbital::SetPQN(int PQN)
{
    pqn = PQN;
}

std::string Orbital::Name() const
{
    return OrbitalInfo(this).Name();
}

void Orbital::Write(FILE* fp) const
{
    // As well as the SpinorFunction vectors, we need to output some other things
    fwrite(&pqn, sizeof(int), 1, fp);
    fwrite(&energy, sizeof(double), 1, fp);
    
    SpinorFunction::Write(fp);
}

void Orbital::Read(FILE* fp)
{
    fread(&pqn, sizeof(int), 1, fp);
    fread(&energy, sizeof(double), 1, fp);
    
    SpinorFunction::Read(fp);
}

bool Orbital::Print(pLattice lattice) const
{
    return Print(Name() + "_orbital.txt", lattice);
}

bool Orbital::Print(const std::string& filename, pLattice lattice) const
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

bool Orbital::Print(FILE* fp, pLattice lattice) const
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

bool Orbital::CheckSize(pLattice lattice, double tolerance)
{
    double maximum = 0.;
    unsigned int i = 0;
    unsigned int max_point = 0;
    while(i < f.size())
    {   if(fabs(f[i]) >= maximum)
        {   maximum = fabs(f[i]);
            max_point = i;
        }
        i++;
    }

    if(maximum < tolerance*100)
    {   *errstream << "Orbital::Checksize: Zero function. " << std::endl;
        PAUSE
        exit(1);
    }

    i = f.size() - 1;
    while(fabs(f[i])/maximum < tolerance)
        i--;

    if(i == f.size() - 1)
    {   // add points to wavefunction
        unsigned int max = f.size();
        double f_max, f_ratio, g_ratio;
        double log_f_ratio, log_g_ratio, dr_max;

        // Strip off any nearby node
        do
        {   max--;
            f_max = fabs(f[max]);
            f_ratio = f[max]/f[max-1];
            g_ratio = g[max]/g[max-1];
        }while(f_ratio < 0. || g_ratio < 0.);

        // Make sure we are tailing off
        if(f_ratio > 0.96)
            f_ratio = 0.96;
        if(g_ratio > 0.96)
            g_ratio = 0.96;

        log_f_ratio = log(f_ratio);
        log_g_ratio = log(g_ratio);
        dr_max = lattice->R(max)-lattice->R(max-1);

        // Resize the state (this is a slight overestimate assuming dr is constant).
        unsigned int old_size = max;
        while(f_max/maximum >= tolerance)
        {   max++;
            f_max = f_max * f_ratio;
        }
        resize(max+1);

        if(lattice->size() < max+1)
            lattice->resize(max+1);

        // Exponential decay (assumes dr changes slowly).
        unsigned int i = old_size;
        while((i < max) && (fabs(f[i])/maximum > tolerance))
        {
            double d2r = (lattice->R(i+1) - lattice->R(i))/dr_max -1.;

            f[i+1] = f[i] * f_ratio * (1. + log_f_ratio * d2r);
            g[i+1] = g[i] * g_ratio * (1. + log_g_ratio * d2r);

            i++;
        }
        resize(i+1);

        return false;
    }
    else if(i+2 < f.size())
    {   // Reduce size
        resize(i+2);
        return false;
    }
    else return true;
}

double Orbital::Norm(pOPIntegrator integrator) const
{
    return integrator->GetNorm(*this);
}

void Orbital::ReNormalise(pOPIntegrator integrator, double norm)
{
    if(norm > 0.)
        (*this) *= sqrt(norm/Norm(integrator));
    else
        Clear();
}

unsigned int Orbital::NumNodes() const
{
    // Count number of zeros
    unsigned int zeros = 0;

    // Get maximum point. This is generally past the last node.
    // We want to ignore small oscillations at the tail of the wavefunction,
    // these are due to the exchange interaction.
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
    while(fabs(f[end]) < 1.e-2 * fmax)
        end--;

    // Get effective start point
    i = 0;
    double prev = f[i];
    while(fabs(prev) < 1.e-7 * fmax)
        prev = f[i++];

    while(i < end)
    {   if(f[i] * prev < 0.)
            zeros++;
        prev = f[i];
        i++;
    }
    return zeros;
}
