#ifndef SUBSET_BASIS_H
#define SUBSET_BASIS_H

#include "ExcitedStates.h"

class SubsetBasis : public ExcitedStates
{
    /** Simple basis formed by taking pointers to a pre-existing basis set (the superset).
        SubsetBasis is on best behaviour not to modify the superset basis.
     */
public:
    SubsetBasis(pLattice lattice, ExcitedStates* superset_basis);
    virtual ~SubsetBasis();
     
    /** Create excited states by taking a subset of the superset_basis.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** This doesn't change the states because we don't own them. */
    virtual void Update() { ClearSigmas(); }

    /** Do not delete states: they are owned by superset basis. */
    virtual void Clear();

    virtual const ExcitedStates* GetSuperset() const;
    virtual void SetSuperset(ExcitedStates* superset_basis);

protected:
    ExcitedStates* superset;
};

#endif
