#ifndef STATE_MANAGER_H
#define STATE_MANAGER_H

#include "Orbital.h"
#include "OrbitalInfo.h"
#include <set>
#include <map>

namespace Ambit
{
/** OrbitalMap is functionally a map from OrbitalInfo to pOrbital.
 */
class OrbitalMap
{
protected:
    typedef std::map<OrbitalInfo, pOrbital> BaseMap;

public:
    OrbitalMap(pLattice lat): lattice(lat) {}
    virtual ~OrbitalMap() {}

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    virtual OrbitalMap* Clone(pLattice new_lattice = nullptr) const;

    typedef BaseMap::iterator iterator;
    typedef BaseMap::const_iterator const_iterator;

    virtual iterator begin() { return m_orbitals.begin(); }
    virtual const_iterator begin() const { return m_orbitals.begin(); }
    virtual void clear() { m_orbitals.clear(); }
    virtual unsigned int count(const OrbitalInfo& key) const { return m_orbitals.count(key); }
    virtual bool empty() const { return m_orbitals.empty(); }
    virtual iterator end() { return m_orbitals.end(); }
    virtual const_iterator end() const { return m_orbitals.end(); }
    virtual iterator erase(const_iterator position) { return m_orbitals.erase(position); }
    virtual iterator find(const OrbitalInfo& key) { return m_orbitals.find(key); }
    virtual const_iterator find(const OrbitalInfo& key) const { return m_orbitals.find(key); }
    virtual unsigned int size() const { return m_orbitals.size(); }

    /** Get pointer to discrete state.
        Return null if no such state exists.
     */
    virtual pOrbitalConst GetState(const OrbitalInfo& info) const;
    virtual pOrbital GetState(const OrbitalInfo& info);

    /** Replaces electron states with those previously stored on disk. Ignores states
        that are not in both sets.
        Manager needs to know whether they're discrete or continuum states.
     */
    virtual void Read(FILE* fp);

    /** Write all electron states to a file. */
    virtual void Write(FILE* fp) const;

    /** Print summary of orbitals to outstream. */
    virtual void Print() const;

    pLattice GetLattice() const { return lattice; }

    /** Add or replace state. */
    virtual void AddState(pOrbital s);

    /** Add all states from another OrbitalMap. Does not replace duplicates. */
    virtual void AddStates(const OrbitalMap& other);

    /** Get size of largest orbital. */
    virtual unsigned int LargestOrbitalSize() const;

protected:
    BaseMap m_orbitals;
    pLattice lattice;
};

typedef std::shared_ptr<OrbitalMap> pOrbitalMap;
typedef std::shared_ptr<const OrbitalMap> pOrbitalMapConst;

}
#endif
