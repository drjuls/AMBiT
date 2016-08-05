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

auto RPAOrbital::GetDeltaPsi(int kappa) -> decltype(deltapsi)::iterator
{
    auto compare_kappas = [&kappa](const std::pair<pDeltaOrbital, pDeltaOrbital>& item)
        {   return item.first->Kappa() == kappa;
        };

    return std::find_if(deltapsi.begin(), deltapsi.end(), compare_kappas);
}

auto RPAOrbital::GetDeltaPsi(int kappa) const -> decltype(deltapsi)::const_iterator
{
    auto compare_kappas = [&kappa](const std::pair<pDeltaOrbital, pDeltaOrbital>& item)
    {   return item.first->Kappa() == kappa;
    };

    return std::find_if(deltapsi.begin(), deltapsi.end(), compare_kappas);
}
