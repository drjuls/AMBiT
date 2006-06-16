#ifndef CUSTOM_BASIS_H
#define CUSTOM_BASIS_H

#include "ExcitedStates.h"
#include "HartreeFock/NonRelInfo.h"

class CustomBasis : public ExcitedStates
{
public:
    CustomBasis(Lattice* lattice, Core* atom_core): ExcitedStates(lattice, atom_core) {}
    virtual ~CustomBasis() {}

    /** Create excited states by following the prescription in CustomBasis.txt
        Only creates up to num_states_per_l states in each wave. If num_states_per_l is zero length,
        it creates all states in CustomBasis.txt.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();

protected:
    virtual NonRelInfo ReadNonRelInfo(char* buffer, unsigned int& num_read);
    std::vector<unsigned int> NumStatesPerL;
};

#endif
