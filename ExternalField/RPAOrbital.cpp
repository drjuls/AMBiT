#include "RPAOrbital.h"

double DeltaOrbital::Energy() const
{
    return parent.lock()->Energy();
}

RPAOrbital::RPAOrbital(const RPAOrbital& other):
    BaseClass(other)
{
    for(const auto& pair: other.deltapsi)
    {
        pDeltaOrbital clone = std::static_pointer_cast<DeltaOrbital>(pair.second->Clone());
        clone->parent = std::static_pointer_cast<RPAOrbital>(shared_from_this());
        deltapsi[pair.first] = clone;
    }
}
