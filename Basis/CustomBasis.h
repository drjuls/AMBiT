#ifndef CUSTOM_BASIS_H
#define CUSTOM_BASIS_H

#include "DiscreteExcitedStates.h"
#include "Configuration/NonRelInfo.h"

class CustomBasis : public DiscreteExcitedStates
{
public:
    CustomBasis(Lattice* lattice, Core* atom_core): DiscreteExcitedStates(lattice, atom_core) {}
    virtual ~CustomBasis() {}

    /** Create excited states by following the prescription in CustomBasis.txt
        Only creates up to num_states_per_l states in each wave, and adds more
        states than specified in CustomBasis.txt if necessary by multiplication by sin(kr).
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();

protected:
    virtual NonRelInfo ReadNonRelInfo(char* buffer, unsigned int& num_read);
};

#endif