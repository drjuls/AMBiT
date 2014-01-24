#ifndef RSTATES_H
#define RSTATES_H

#include "ExcitedStates.h"

class RStates : public ExcitedStates
{
public:
    RStates(pLattice lattice, Core* atom_core): ExcitedStates(lattice, atom_core) {}
    virtual ~RStates(void) {}

    /** Create excited states by multiplying states by R and orthogonalising.
        The first excited state is made by HF iterations.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();

protected:
    std::vector<unsigned int> NumStatesPerL;
};

#endif
