#ifndef STATE_MANAGER_H
#define STATE_MANAGER_H

#include "Orbital.h"
#include "OrbitalInfo.h"
#include <set>
#include <map>

typedef std::map<OrbitalInfo, pOrbital> OrbitalMap;

/** StateManager is functionally a map from OrbitalInfo to pOrbital.
 */
class StateManager
{
public:
    StateManager(pLattice lat): lattice(lat) {}
    StateManager(const StateManager& other);
    virtual ~StateManager(void);

    const StateManager& operator=(const StateManager& other);

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    const StateManager& Copy(const StateManager& other, pLattice new_lattice = pLattice());
    StateManager Copy(pLattice new_lattice = pLattice()) const;

    typedef OrbitalMap::iterator iterator;
    typedef OrbitalMap::const_iterator const_iterator;

    virtual iterator begin() { return AllStates.begin(); }
    virtual const_iterator begin() const { return AllStates.begin(); }
    virtual void clear() { AllStates.clear(); }
    virtual bool empty() const { return AllStates.empty(); }
    virtual iterator end() { return AllStates.end(); }
    virtual const_iterator end() const { return AllStates.end(); }
    virtual unsigned int size() const { return static_cast<unsigned int>(AllStates.size()); }

    /** Get pointer to discrete state.
        Return null if no such state exists.
     */
    virtual pOrbitalConst GetState(const OrbitalInfo& info) const;
    virtual pOrbital GetState(const OrbitalInfo& info);

    /** Write all electron states to a file. */
    virtual void Write(FILE* fp) const;

    /** Replaces electron states with those previously stored on disk. Ignores states
        that are not in both sets.
        Manager needs to know whether they're discrete or continuum states.
     */
    virtual void Read(FILE* fp);

    pLattice GetLattice() const { return lattice; }

    /** Add or replace state. */
    virtual void AddState(pOrbital s);

    /** Add all states from another StateManager. */
    virtual void AddStates(StateManager& other);

protected:
    OrbitalMap AllStates;
    pLattice lattice;
};

typedef boost::shared_ptr<StateManager> pStateManager;
typedef boost::shared_ptr<const StateManager> pStateManagerConst;

#endif
