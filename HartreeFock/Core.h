#ifndef CORE_H
#define CORE_H

#include "OrbitalMap.h"
#include "Atom/Debug.h"
#include "SigmaPotential.h"
#include "Configuration.h"

typedef Configuration<OrbitalInfo, double> OccupationMap;

/** Core is a OrbitalMap container class with orbital occupancies added. */
class Core : public OrbitalMap
{
public:
    Core(pLattice lat, const std::string& filling = "");
    Core(const Core& other);
    Core(Core&& other);
    Core(const OrbitalMap& other);  //!< Create closed-shell core with fully occupied orbitals from other.
    virtual ~Core() {}

    const Core& operator=(const Core& other);

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    const Core& Clone(const Core& other, pLattice new_lattice = pLattice());
    Core Clone(pLattice new_lattice = pLattice()) const;

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

    virtual int NumElectrons() const;

protected:
    OccupationMap occupancy;
};

typedef boost::shared_ptr<Core> pCore;
typedef boost::shared_ptr<const Core> pCoreConst;

#endif
