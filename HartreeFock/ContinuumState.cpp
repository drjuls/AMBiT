#include "Include.h"
#include "ContinuumState.h"
#include "Universal/Constant.h"

ContinuumState::ContinuumState(Lattice* lat, unsigned int num_points):
    State(lat, num_points)
{}

ContinuumState::ContinuumState(Lattice* lat, double Nu, int Kappa, unsigned int num_points):
    State(lat, Kappa, num_points)
{
    nu = Nu;
}

ContinuumState::ContinuumState(const ContinuumState& other):
    State(other)
{}

std::string ContinuumState::Name() const
{
    char buffer[20];
    sprintf(buffer, "%5.4f", nu);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    if(kappa < -1)
        ret.append(1, '+');

    return ret;
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
