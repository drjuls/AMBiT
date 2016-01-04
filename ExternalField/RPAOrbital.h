#ifndef RPA_ORBITAL_H
#define RPA_ORBITAL_H

#include "HartreeFock/Orbital.h"
#include <map>

class RPAOrbital;

/** DeltaOrbital is a perturbation to an existing orbital (which is an RPAOrbital).
 It has a definite kappa.
 The PQN is zero: it cannot be uniquely identified from OrbitalInfo.
 */
class DeltaOrbital : public Orbital
{
public:
    DeltaOrbital(int kappa, std::shared_ptr<RPAOrbital> parent):
    Orbital(kappa), deltaEnergy(0.), parent(parent)
    {}
    DeltaOrbital(const DeltaOrbital& other):
    Orbital(other), deltaEnergy(other.deltaEnergy), parent(other.parent)
    {}
    DeltaOrbital(DeltaOrbital&& other):
    Orbital(other), deltaEnergy(other.deltaEnergy), parent(other.parent)
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
class RPAOrbital : public Orbital
{
public:
    RPAOrbital(const OrbitalInfo& info): Orbital(info) {}
    RPAOrbital(int kappa, int pqn = 0, double energy = 0.0, unsigned int size = 0): Orbital(kappa, pqn, energy, size) {}
    RPAOrbital(const Orbital& other): Orbital(other) {}
    RPAOrbital(const RPAOrbital& other);
    RPAOrbital(RPAOrbital&& other): Orbital(other), deltapsi(other.deltapsi) {}
    virtual ~RPAOrbital() {}

    /** Map from kappa to orbital contribution. */
    std::map<int, pDeltaOrbital> deltapsi;
};

typedef std::shared_ptr<RPAOrbital> pRPAOrbital;
typedef std::shared_ptr<const RPAOrbital> pRPAOrbitalConst;

#endif
