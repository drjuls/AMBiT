#include "Include.h"
#include "Orbital.h"
#include "Universal/PhysicalConstant.h"
#include "OrbitalInfo.h"
#include <math.h>

namespace Ambit
{
OrbitalBase::OrbitalBase(const OrbitalInfo& info):
    SpinorFunction(info.Kappa()), pqn(info.PQN()), energy(0.0)
{}

OrbitalBase::OrbitalBase(int kappa, int pqn, double energy, unsigned int size):
    SpinorFunction(kappa, size), pqn(pqn), energy(energy)
{}

double OrbitalBase::Energy() const
{
    return energy;
}

/** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
double OrbitalBase::Nu() const
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

int OrbitalBase::PQN() const
{
    return pqn;
}

void OrbitalBase::SetEnergy(double Energy)
{
    energy = Energy;
}

void OrbitalBase::SetNu(double Nu)
{
    if(fabs(Nu) > 1.e-5)
        energy = -0.5/(Nu * Nu);
    else
        energy = -1.e10;
    
    if(Nu < 0)
        energy = -energy;
}

void OrbitalBase::SetPQN(int PQN)
{
    pqn = PQN;
}

std::string OrbitalBase::Name() const
{
    return OrbitalInfo(std::static_pointer_cast<const Orbital>(shared_from_this())).Name();
}

void OrbitalBase::Write(FILE* fp) const
{
    // As well as the SpinorFunction vectors, we need to output some other things
    file_err_handler->fwrite(&pqn, sizeof(int), 1, fp);
    file_err_handler->fwrite(&energy, sizeof(double), 1, fp);
    
    SpinorFunction::Write(fp);
}

void OrbitalBase::Read(FILE* fp)
{
    file_err_handler->fread(&pqn, sizeof(int), 1, fp);
    file_err_handler->fread(&energy, sizeof(double), 1, fp);
    
    SpinorFunction::Read(fp);
}

bool OrbitalBase::Print(pLattice lattice) const
{
    return Print(Name() + "_orbital.txt", lattice);
}

bool OrbitalBase::Print(const std::string& filename, pLattice lattice) const
{
    FILE* fp = file_err_handler->fopen(filename.c_str(), "wt");
    
    if(fp)
    {   bool success = Print(fp, lattice);
        file_err_handler->fclose(fp);
        return success;
    }
    else
        return false;
}

bool OrbitalBase::Print(FILE* fp, pLattice lattice, bool jacobian) const
{
    unsigned int i;
    if(lattice)
    {
        if(jacobian)
        {   for(i = 0; i < size(); i++)
            {
                fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",
                        lattice->R(i), lattice->dR(i),
                        f[i], g[i], dfdr[i], dgdr[i]);
            }
        }
        else
        {   for(i = 0; i < size(); i++)
            {
                fprintf(fp, "%12.6e %12.6e %12.6e %12.6e %12.6e\n", lattice->R(i),
                        f[i], g[i], dfdr[i], dgdr[i]);
            }
        }
    }
    else
    {   for(i = 0; i < size(); i++)
        {
            fprintf(fp, "%d %12.6e %12.6e %12.6e %12.6e\n", i,
                    f[i], g[i], dfdr[i], dgdr[i]);
        }
    }

    return true;
}

bool OrbitalBase::CheckSize(pLattice lattice, double tolerance)
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
    {   *errstream << "Orbital::Checksize: Zero function in " << Name() << std::endl;
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

        if(lattice->size() <= max+1)
            lattice->resize(max+2);

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

double OrbitalBase::Norm(pIntegrator integrator) const
{
    return integrator->GetNorm(*this);
}

void OrbitalBase::ReNormalise(pIntegrator integrator, double norm)
{
    if(norm > 0.)
        (*this) *= sqrt(norm/Norm(integrator));
    else
        Clear();
}

unsigned int OrbitalBase::NumNodes() const
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
}
