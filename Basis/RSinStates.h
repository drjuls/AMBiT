#ifndef R_SIN_STATES_H
#define R_SIN_STATES_H

#include "RStates.h"

class RSinStates : public RStates
{
public:
    RSinStates(Lattice* lattice, Core* atom_core): RStates(lattice, atom_core) {}
    virtual ~RSinStates(void) {}

    /** Create excited states by multiplying states by R and orthogonalising.
        The first excited state is made by HF iterations.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);
};

#endif