#ifndef READ_BASIS_H
#define READ_BASIS_H

#include "Include.h"
#include "ExcitedStates.h"

class ReadBasis : public ExcitedStates
{
    /** Read basis functions from a GRASP file.
     */
public:
    ReadBasis(pLattice lattice, Core* atom_core, const std::string& input_file):
        ExcitedStates(lattice, atom_core), filename(input_file) {}
    virtual ~ReadBasis() {}

    /** Create excited states from the file.
        Only creates up to num_states_per_l states in each wave.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Create excited states from the file, but also update any core orbitals from the file.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l, Core* atom_core);

    /** Update all of the excited states because the core has changed.
        In the case of externally read basis sets, this does nothing.
     */
    virtual void Update() {}

protected:
    std::string filename;
};

#endif
