#ifndef HARTREE_FOCK_BASIS_H
#define HARTREE_FOCK_BASIS_H

#include "ExcitedStates.h"

class HartreeFockBasis : public ExcitedStates
{
public:
    HartreeFockBasis(Lattice* lattice, Core* atom_core): ExcitedStates(lattice, atom_core) {}
    virtual ~HartreeFockBasis() {}

    /** Create excited states by following the prescription in CustomBasis.txt
        Only creates up to num_states_per_l states in each wave. If num_states_per_l is zero length,
        it creates all states in CustomBasis.txt.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();
};

#endif
