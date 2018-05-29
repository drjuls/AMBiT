#ifndef CORE_H
#define CORE_H

#include "OrbitalMap.h"
#include "Atom/Debug.h"
#include "Configuration.h"

namespace Ambit
{
typedef Configuration<OrbitalInfo, double> OccupationMap;

/** Core is a OrbitalMap container class with orbital occupancies added. */
class Core : public OrbitalMap
{
public:
    Core(pLattice lat, const std::string& filling = "");
    Core(const OrbitalMap& other);  //!< Create closed-shell core with fully occupied orbitals from other.
    virtual ~Core() {}

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    virtual Core* Clone(pLattice new_lattice = nullptr) const override;

    virtual double GetOccupancy(OrbitalInfo info) const;
    virtual double GetOccupancy(pOrbital info) const;

    /** Set occupancies for orbitals.
        Remove orbitals with no occupancy and create orbitals if they don't exist already.
     */
    virtual void SetOccupancies(const OccupationMap& occupancies);
    virtual const OccupationMap& GetOccupancies() const;

    virtual void Write(FILE* fp) const override;
    virtual void Read(FILE* fp) override;

    /** Delete all currently stored states. */
    virtual void clear() override;
    virtual iterator erase(const_iterator position) override;

    virtual int NumElectrons() const;

protected:
    OccupationMap occupancy;
};

typedef std::shared_ptr<Core> pCore;
typedef std::shared_ptr<const Core> pCoreConst;

}
#endif
