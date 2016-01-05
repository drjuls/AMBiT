#ifndef RPA_ORBITAL_H
#define RPA_ORBITAL_H

#include "HartreeFock/Orbital.h"
#include <map>

class RPAOrbital;

/** DeltaOrbital is a perturbation to an existing orbital (which is an RPAOrbital).
 It has a definite kappa.
 The PQN is zero: it cannot be uniquely identified from OrbitalInfo.
 */
class DeltaOrbital : public OrbitalTemplate<Orbital, DeltaOrbital>
{
public:
    DeltaOrbital(int kappa, std::shared_ptr<RPAOrbital> parent):
        BaseClass(kappa), deltaEnergy(0.), parent(parent)
    {}

    /** Energy() should return the parent energy. */
    virtual double Energy() const override;

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
    RPAOrbital(const RPAOrbital& other);    //!< Deep copy of deltapsi
    virtual ~RPAOrbital() {}

    /** Clone makes deep copy of deltapsi. */
    virtual pOrbital Clone() const
    {
        pOrbital ret = std::make_shared<RPAOrbital>(static_cast<const RPAOrbital&>(*this));
        return ret;
    }

    /** Map from kappa to orbital contribution. */
    std::map<int, pDeltaOrbital> deltapsi;
};

typedef std::shared_ptr<RPAOrbital> pRPAOrbital;
typedef std::shared_ptr<const RPAOrbital> pRPAOrbitalConst;

#endif
