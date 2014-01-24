#include "SpinorODE.h"

SpinorODE::SpinorODE(pLattice lattice): lattice(lattice), core(NULL), include_nonlocal(true)
{}

SpinorODE::SpinorODE(const Core* core): lattice(core->GetLattice()), include_nonlocal(true)
{
    SetCore(core);
}

void SpinorODE::SetCore(const Core* hf_core)
{
    core = hf_core;
}

const Core* SpinorODE::GetCore() const
{
    return core;
}

void SpinorODE::IncludeExchangeInODE(bool include_exchange)
{
    include_nonlocal = include_exchange;
}

void SpinorODE::GetDerivative(SingleParticleWavefunction& fg)
{
    SetODEParameters(fg);

    double w[2];
    for(unsigned int i = 0; i < fg.Size(); i++)
    {
        GetODEFunction(i, fg, w);
        fg.dfdr[i] = w[0];
        fg.dgdr[i] = w[1];
    }
}
