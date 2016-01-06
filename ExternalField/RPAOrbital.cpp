#include "RPAOrbital.h"

double DeltaOrbital::Energy() const
{
    return parent.lock()->Energy();
}

pDeltaOrbital DeltaOrbital::CloneWithNewParent(std::shared_ptr<RPAOrbital> new_parent) const
{
    pDeltaOrbital ret = std::make_shared<DeltaOrbital>(static_cast<DeltaOrbital const &>(*this));
    ret->parent = new_parent;
    return ret;
}

RPAOrbital::RPAOrbital(const RPAOrbital& other):
    BaseClass(other)
{
    for(const auto& pair: other.deltapsi)
    {
        pDeltaOrbital clone = pair.second->CloneWithNewParent(std::static_pointer_cast<RPAOrbital>(shared_from_this()));
        deltapsi[pair.first] = clone;
    }
}
