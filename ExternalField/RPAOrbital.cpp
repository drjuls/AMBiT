#include "RPAOrbital.h"

namespace Ambit
{
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

const RPAOrbital& RPAOrbital::operator*=(double scale_factor)
{
    for(auto& deltaorbitals: deltapsi)
    {
        if(deltaorbitals.first)
            *(deltaorbitals.first) *= scale_factor;
        if(deltaorbitals.second)
            *(deltaorbitals.second) *= scale_factor;
    }

    BaseClass::operator*=(scale_factor);
    return *this;
}

/** Adding or subtracting two spinor functions can only occur if both have same angular part. */
const RPAOrbital& RPAOrbital::operator+=(const RPAOrbital& other)
{
    for(auto& deltaorbitals: deltapsi)
    {
        auto other_it = other.GetDeltaPsi(deltaorbitals.first->Kappa());
        if(other_it != other.deltapsi.end())
        {
            *(deltaorbitals.first) += *(other_it->first);
            if(other_it->second)
            {
                if(deltaorbitals.second)
                    *(deltaorbitals.second) += *(other_it->second);
                else
                    deltaorbitals.second = other_it->second->CloneWithNewParent(std::static_pointer_cast<RPAOrbital>(shared_from_this()));
            }
        }
    }

    BaseClass::operator+=(other);
    return *this;
}

const RPAOrbital& RPAOrbital::operator-=(const RPAOrbital& other)
{
    (*this) += other * (-1.0);
    return *this;
}

RPAOrbital RPAOrbital::operator+(const RPAOrbital& other) const
{
    RPAOrbital ret(*this);
    return ret += other;
}

RPAOrbital RPAOrbital::operator-(const RPAOrbital& other) const
{
    RPAOrbital ret(*this);
    return ret -= other;
}

const RPAOrbital& RPAOrbital::operator*=(const RadialFunction& chi)
{
    for(auto& deltaorbitals: deltapsi)
    {
        if(deltaorbitals.first)
            *(deltaorbitals.first) *= chi;
        if(deltaorbitals.second)
            *(deltaorbitals.second) *= chi;
    }

    BaseClass::operator*=(chi);
    return *this;
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
}
