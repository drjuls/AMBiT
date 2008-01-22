#include "Include.h"
#include "ContinuumState.h"
#include "Universal/Constant.h"

ContinuumState::ContinuumState(const ContinuumState& other):
    State(other)
{}

ContinuumState::ContinuumState(double energy, int Kappa):
    State(Kappa)
{
    SetEnergy(energy);
}

std::string ContinuumState::Name() const
{
    char buffer[20];
    sprintf(buffer, "%5.4f", nu);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);

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

void ContinuumState::SetEnergy(double energy)
{
    if(energy > 0)
        nu = sqrt(0.5/energy);
    else
    {   *errstream << "ContinuumState: energy <= 0" << std::endl;
        exit(1);
    }
}


void ContinuumState::Write(FILE* fp) const
{
    // As well as the CoupledFunction vectors, we need to output some other things
    fwrite(&kappa, sizeof(int), 1, fp);
    fwrite(&nu, sizeof(double), 1, fp);

    State::Write(fp);
}

void ContinuumState::Read(FILE* fp)
{
    fread(&kappa, sizeof(int), 1, fp);
    fread(&nu, sizeof(double), 1, fp);

    State::Read(fp);
}
