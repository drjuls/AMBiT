#ifndef DISCRETE_EXCITED_STATES_H
#define DISCRETE_EXCITED_STATES_H

#include "ExcitedStates.h"
#include "HartreeFock/DiscreteStateIterator.h"

class DiscreteExcitedStates : public ExcitedStates
{
public:
    friend class DiscreteStateIterator;
    friend class ConstDiscreteStateIterator;

public:
    DiscreteExcitedStates(Lattice* lattice, Core* atom_core): ExcitedStates(lattice, atom_core) {}
    virtual ~DiscreteExcitedStates() {}

    /** These override the methods in StateManager to return a DiscreteState. */
    virtual DiscreteState* GetState(const StateInfo& info);
    virtual const DiscreteState* GetState(const StateInfo& info) const;
    
    virtual DiscreteStateIterator GetDiscreteStateIterator();
    virtual ConstDiscreteStateIterator GetConstDiscreteStateIterator() const;

    virtual void Read(FILE* fp);
    virtual void Write(FILE* fp) const;

protected:
    /** current.f = r * previous.f
        PRE: previous.kappa = current.kappa
     */
    void MultiplyByR(const DiscreteState* previous, DiscreteState* current) const;

    /** current.f = sin(kr) * previous.f, where k = Pi/R_max
        PRE: previous.kappa = current.kappa
     */
    void MultiplyBySinR(const DiscreteState* previous, DiscreteState* current) const;

    /** current.f = r * sin(kr) * previous.f
        PRE: current.l = previous.l + 1
     */
    void MultiplyByRSinR(const DiscreteState* previous, DiscreteState* current) const;

    void Orthogonalise(DiscreteState* current) const;
};

#endif
