#include "SpinorODE.h"

SpinorODE::SpinorODE(Lattice* lattice): lattice(lattice), core(NULL), include_nonlocal(true)
{}

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
