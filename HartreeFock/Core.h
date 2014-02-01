#ifndef CORE_H
#define CORE_H

#include "StateManager.h"
#include "StateIterator.h"
#include "Atom/Debug.h"
#include "SigmaPotential.h"
#include "ConfigurationParser.h"

/** Core is a StateManager container class with orbital occupancies added. */
class Core : public StateManager
{
public:
    Core(pLattice lat, const std::string& filling = "");
    Core(const Core& other);
    virtual ~Core() {}

    const Core& operator=(const Core& other);

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    const Core& Copy(const Core& other, pLattice new_lattice = pLattice());
    Core Copy(pLattice new_lattice = pLattice()) const;

    virtual double GetOccupancy(OrbitalInfo info) const;
    virtual double GetOccupancy(pOrbital info) const;

    /** Set occupancies for orbitals.
        Remove orbitals with no occupancy and create orbitals if they don't exist already.
     */
    virtual void SetOccupancies(const OccupationMap& occupancies);
    virtual const OccupationMap& GetOccupancies() const;

    virtual void Write(FILE* fp) const;
    virtual void Read(FILE* fp);

    /** Delete all currently stored states. */
    virtual void Clear();

    virtual int NumElectrons() const;

protected:
    OccupationMap occupancy;
};

#endif
