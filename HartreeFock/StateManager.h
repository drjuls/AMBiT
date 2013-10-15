#ifndef STATE_MANAGER_H
#define STATE_MANAGER_H

#include "Orbital.h"
#include "StatePointer.h"
#include "OrbitalInfo.h"
#include <set>
#include <map>

typedef std::map<OrbitalInfo, StatePointer> StateSet;

class StateIterator;
class ConstStateIterator;

class StateManager
{
public:
    friend class StateIterator;
    friend class ConstStateIterator;

public:
    StateManager(Lattice* lat, unsigned int atomic_number, int ion_charge);
    /** Copy all states, interpolating onto new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice (same object, no copy made).
     */
    StateManager(const StateManager& other, Lattice* new_lattice = NULL);
    virtual ~StateManager(void);

    const StateManager& operator=(const StateManager& other);

    virtual bool Empty() const { return AllStates.empty(); }
    virtual unsigned int NumStates() const { return static_cast<unsigned int>(AllStates.size()); }

    /** Get pointer to discrete state.
        Return null if no such state exists.
     */
    virtual const Orbital* GetState(const OrbitalInfo& info) const;
    virtual Orbital* GetState(const OrbitalInfo& info);

    virtual StateIterator GetStateIterator();
    virtual ConstStateIterator GetConstStateIterator() const;

    /** Write all electron states to a file. */
    virtual void Write(FILE* fp) const;

    /** Replaces electron states with those previously stored on disk. Ignores states
        that are not in both sets.
        Manager needs to know whether they're discrete or continuum states.
     */
    virtual void Read(FILE* fp);

    Lattice* GetLattice() const { return lattice; }

    /** Test for orthogonality of states.
        Return largest overlap.
     */
    double TestOrthogonality() const;

    virtual void AddState(Orbital* s);

protected:
    /** Delete all currently stored states. */
    virtual void Clear();

    double Z, Charge;
    StateSet AllStates;

    Lattice* lattice;
};

#endif
