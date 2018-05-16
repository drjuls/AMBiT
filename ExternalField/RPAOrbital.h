#ifndef RPA_ORBITAL_H
#define RPA_ORBITAL_H

#include "HartreeFock/Orbital.h"

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

    /** Get DeltaEnergy: the shift in energy of parent due to this DeltaOrbital.
        Warning: This is the "reduced" energy, < a || f + deltaV || a > where a is the parent.
     */
    virtual double DeltaEnergy() const { return deltaEnergy; }
    virtual void SetDeltaEnergy(double energy) { deltaEnergy = energy; }

protected:
    double deltaEnergy = 0.;
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
    RPAOrbital(const RPAOrbital& other) = default;
    RPAOrbital(RPAOrbital&& other) = default;
    virtual ~RPAOrbital() = default;

    /** Clone makes deep copy of deltapsi. */
    virtual pOrbital Clone() const override
    {
        std::shared_ptr<RPAOrbital> ret = std::make_shared<RPAOrbital>(static_cast<const Orbital&>(*this));

        for(const auto& pair: deltapsi)
        {
            pDeltaOrbital first_clone = (pair.first? pair.first->CloneWithNewParent(ret): nullptr);
            pDeltaOrbital second_clone = (pair.second? pair.second->CloneWithNewParent(ret): nullptr);
            ret->deltapsi.emplace_back(std::make_pair(first_clone, second_clone));
        }
        return ret;
    }

    /** Multiply orbital and deltapsi by the scale factor. */
    const RPAOrbital& operator*=(double scale_factor);

    /** Adding or subtracting two RPAOrbitals can only occur if both have same angular part. */
    const RPAOrbital& operator+=(const RPAOrbital& other);
    const RPAOrbital& operator-=(const RPAOrbital& other);
    RPAOrbital operator+(const RPAOrbital& other) const;
    RPAOrbital operator-(const RPAOrbital& other) const;

    /** Multiply spinor function by another radial function (assumed zero outside range). */
    const RPAOrbital& operator*=(const RadialFunction& chi);

    /** Orbital contributions for different kappas: pair<dPsi(+), dPsi(-)>. */
    std::vector<std::pair<pDeltaOrbital, pDeltaOrbital>> deltapsi;

    /** Search deltaPsi for pair with given kappa. Returns deltapsi.end() if none found. */
    auto GetDeltaPsi(int kappa) -> decltype(deltapsi)::iterator;
    auto GetDeltaPsi(int kappa) const -> decltype(deltapsi)::const_iterator;
};

typedef std::shared_ptr<RPAOrbital> pRPAOrbital;
typedef std::shared_ptr<const RPAOrbital> pRPAOrbitalConst;

#endif
