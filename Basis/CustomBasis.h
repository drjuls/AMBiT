#ifndef CUSTOM_BASIS_H
#define CUSTOM_BASIS_H

#include "ExcitedStates.h"
#include "HartreeFock/NonRelInfo.h"

class CustomBasis : public ExcitedStates
{
public:
    CustomBasis(pLattice lattice, pCoreConst atom_core):
        ExcitedStates(lattice, atom_core), filename("CustomBasis.txt")
    {}
    virtual ~CustomBasis() {}

    /** Create excited states by following the prescription in filename (CustomBasis.txt).
        Only creates up to num_states_per_l states in each wave. If num_states_per_l is zero length,
        it creates all states in CustomBasis.txt.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l, const std::string& file);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();

    virtual void SetFile(const std::string& file);

protected:
    virtual NonRelInfo ReadNonRelInfo(char* buffer, unsigned int& num_read);
    std::vector<unsigned int> NumStatesPerL;
    std::string filename;
};

#endif
