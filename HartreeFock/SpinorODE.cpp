#include "SpinorODE.h"

namespace Ambit
{
SpinorODE::SpinorODE(pLattice lattice): LatticeObserver(lattice), include_nonlocal(true)
{}

void SpinorODE::IncludeExchange(bool include_exchange)
{
    include_nonlocal = include_exchange;
}

bool SpinorODE::IncludeExchange() const
{
    return include_nonlocal;
}

void SpinorODE::GetDerivative(Orbital& fg, bool set_parameters)
{
    if(set_parameters)
        SetODEParameters(fg);

    GetDerivative(fg);
}

void SpinorODE::GetDerivative(Orbital& fg) const
{
    double w[2];
    for(unsigned int i = 0; i < fg.size(); i++)
    {
        GetODEFunction(i, fg, w);
        fg.dfdr[i] = w[0];
        fg.dgdr[i] = w[1];
    }
}
}
