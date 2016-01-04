#include "RPAOrbital.h"

double DeltaOrbital::Energy() const
{
    return parent.lock()->Energy();
}

RPAOrbital::RPAOrbital(const RPAOrbital& other):
Orbital(other)
{
    for(const auto& pair: other.deltapsi)
    {
        pDeltaOrbital clone(new DeltaOrbital(*pair.second));
        deltapsi[pair.first] = clone;
    }
}
