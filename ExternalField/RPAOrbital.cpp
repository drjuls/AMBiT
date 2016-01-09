#include "RPAOrbital.h"

double DeltaOrbital::Energy() const
{
    pOrbitalConst sp = parent.lock();
    if(sp)
        return sp->Energy();
    else
        return 0.;
}

pDeltaOrbital DeltaOrbital::CloneWithNewParent(std::shared_ptr<RPAOrbital> new_parent) const
{
    pDeltaOrbital ret = std::make_shared<DeltaOrbital>(static_cast<DeltaOrbital const &>(*this));
    ret->parent = new_parent;
    return ret;
}
