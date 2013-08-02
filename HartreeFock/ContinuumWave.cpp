#include "Include.h"
#include "ContinuumWave.h"
#include "Universal/MathConstant.h"

ContinuumWave::ContinuumWave(const ContinuumWave& other):
    SingleParticleWavefunction(other)
{}

ContinuumWave::ContinuumWave(double energy, int Kappa):
    SingleParticleWavefunction(Kappa)
{
    SetEnergy(energy);
}

std::string ContinuumWave::Name() const
{
    char buffer[20];
    sprintf(buffer, "%5.4f", nu);
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

void ContinuumWave::SetEnergy(double energy)
{
    if(energy > 0)
        nu = sqrt(0.5/energy);
    else
    {   *errstream << "ContinuumWave: energy <= 0" << std::endl;
        exit(1);
    }
}


void ContinuumWave::Write(FILE* fp) const
{
    // As well as the CoupledFunction vectors, we need to output some other things
    fwrite(&kappa, sizeof(int), 1, fp);
    fwrite(&nu, sizeof(double), 1, fp);

    SingleParticleWavefunction::Write(fp);
}

void ContinuumWave::Read(FILE* fp)
{
    fread(&kappa, sizeof(int), 1, fp);
    fread(&nu, sizeof(double), 1, fp);

    SingleParticleWavefunction::Read(fp);
}
