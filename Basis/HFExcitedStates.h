#ifndef HF_EXCITED_STATES_H
#define HF_EXCITED_STATES_H

#include "ExcitedStates.h"

class HFExcitedIterator;
class ConstHFExcitedIterator;

class HFExcitedStates : public ExcitedStates
{
public:
    friend class HFExcitedIterator;
    friend class ConstHFExcitedIterator;

public:
    HFExcitedStates(Lattice* lattice, Core* atom_core): ExcitedStates(lattice, atom_core) {}
    virtual ~HFExcitedStates() {}

    virtual void AddState(State* s);

    virtual void Write(FILE* fp) const;
    virtual void Read(FILE* fp);

public:
    /** Create Hartree-Fock discrete virtual states. */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);
    virtual void Update();

    /** Build the continuum states.
        Creates num_nu * num_kappa states, 
          spaced quadratically in nu from start_nu to end_nu
          with kappa from -1, 1, -2, 2, ... up to num_kappa of them
        PRE: end_nu >= num_states * start_nu
             num_kappa should be odd to have all states with any given L.
     */
    virtual void CreateContinuum(double start_nu, double end_nu, unsigned int num_nu, unsigned int num_kappa);

protected:
    std::map<int, unsigned int> LastPQNForKappa;
    std::set<double> ContinuumQNs;
};

class HFExcitedIterator : public StateIterator
{
public:
    HFExcitedIterator(HFExcitedStates* manager):
      StateIterator(manager), excited_manager(manager) {}
    virtual ~HFExcitedIterator(void) {}

    virtual double Weight();

protected:
    HFExcitedStates* excited_manager;
};

class ConstHFExcitedIterator : public ConstStateIterator
{
public:
    ConstHFExcitedIterator(const HFExcitedStates* manager):
      ConstStateIterator(manager), excited_manager(manager) {}
    virtual ~ConstHFExcitedIterator(void) {}

    virtual double Weight();

protected:
    const HFExcitedStates* excited_manager;
};

#endif
