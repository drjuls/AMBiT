#ifndef RPA_ORBITAL_H
#define RPA_ORBITAL_H

#include "HartreeFock/Orbital.h"
#include <map>

class RPAOrbital;

/** DeltaOrbital is a perturbation to an existing orbital (which is an RPAOrbital).
    Member DeltaEnergy() is the shift in energy of parent due to this DeltaOrbital.
    It has a definite kappa.
    The PQN is zero: it cannot be uniquely identified from OrbitalInfo.
 */
class DeltaOrbital : public OrbitalTemplate<Orbital, DeltaOrbital>
{
public:
    DeltaOrbital(int kappa, std::shared_ptr<RPAOrbital> parent):
        BaseClass(kappa), deltaEnergy(0.), parent(parent)
    {}
    DeltaOrbital(const Orbital& other):
        BaseClass(other), deltaEnergy(0.), parent()
    {}

    /** Get an identical Clone() except with a new parent. */
    std::shared_ptr<DeltaOrbital> CloneWithNewParent(std::shared_ptr<RPAOrbital> new_parent) const;

    /** Return parent energy. */
    virtual double Energy() const override;

    /** Get parent. */
    virtual std::shared_ptr<RPAOrbital> GetParent() const { return parent.lock(); }

    /** Set parent. */
    virtual void SetParent(std::shared_ptr<RPAOrbital> new_parent) { parent = new_parent; }

    /** Get DeltaEnergy: shift in energy of parent due to this DeltaOrbital. */
    virtual double DeltaEnergy() const { return deltaEnergy; }
    virtual void SetDeltaEnergy(double energy) { deltaEnergy = energy; }

protected:
    double deltaEnergy;
    std::weak_ptr<RPAOrbital> parent;
};

typedef std::shared_ptr<DeltaOrbital> pDeltaOrbital;
typedef std::shared_ptr<const DeltaOrbital> pDeltaOrbitalConst;

/** RPAOrbital can have first-order perturbations with different kappa: deltapsi.
 */
class RPAOrbital : public OrbitalTemplate<Orbital, RPAOrbital>
{
public:
    RPAOrbital(const OrbitalInfo& info): BaseClass(info) {}
    RPAOrbital(int kappa, int pqn = 0, double energy = 0.0, unsigned int size = 0):
        BaseClass(kappa, pqn, energy, size) {}
    RPAOrbital(const Orbital& other): BaseClass(other) {}
    virtual ~RPAOrbital() {}

    /** Clone makes deep copy of deltapsi. */
    virtual pOrbital Clone() const override
    {
        std::shared_ptr<RPAOrbital> ret = std::make_shared<RPAOrbital>(static_cast<const Orbital&>(*this));

        for(const auto& pair: deltapsi)
        {
            pDeltaOrbital clone = pair.second->CloneWithNewParent(ret);
            ret->deltapsi[pair.first] = clone;
        }
        return ret;
    }

    /** Map from kappa to orbital contribution. */
    std::map<int, pDeltaOrbital> deltapsi;
};

typedef std::shared_ptr<RPAOrbital> pRPAOrbital;
typedef std::shared_ptr<const RPAOrbital> pRPAOrbitalConst;

#endif
