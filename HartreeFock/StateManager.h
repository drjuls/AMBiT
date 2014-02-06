#ifndef STATE_MANAGER_H
#define STATE_MANAGER_H

#include "Orbital.h"
#include "OrbitalInfo.h"
#include <set>
#include <map>

typedef std::map<OrbitalInfo, pOrbital> StateSet;

class StateIterator;
class ConstStateIterator;

class StateManager
{
public:
    friend class StateIterator;
    friend class ConstStateIterator;

public:
    StateManager(pLattice lat);
    StateManager(const StateManager& other);
    virtual ~StateManager(void);

    const StateManager& operator=(const StateManager& other);

    /** Deep copy of all orbitals, interpolating from new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice.
     */
    const StateManager& Copy(const StateManager& other, pLattice new_lattice = pLattice());
    StateManager Copy(pLattice new_lattice = pLattice()) const;

    virtual bool Empty() const { return AllStates.empty(); }
    virtual unsigned int NumStates() const { return static_cast<unsigned int>(AllStates.size()); }

    /** Get pointer to discrete state.
        Return null if no such state exists.
     */
    virtual pOrbitalConst GetState(const OrbitalInfo& info) const;
    virtual pOrbital GetState(const OrbitalInfo& info);

    virtual StateIterator GetStateIterator();
    virtual ConstStateIterator GetConstStateIterator() const;

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

    /** Delete all currently stored states. */
    virtual void Clear();

protected:
    StateSet AllStates;
    pLattice lattice;
};

typedef boost::shared_ptr<StateManager> pStateManager;
typedef boost::shared_ptr<const StateManager> pStateManagerConst;

#endif
