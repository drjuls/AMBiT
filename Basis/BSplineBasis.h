#ifndef BSPLINE_BASIS_H
#define BSPLINE_BASIS_H

#include "DiscreteExcitedStates.h"

class BSplineBasis : public DiscreteExcitedStates
{
public:
    BSplineBasis(Lattice* lattice, Core* atom_core):
        DiscreteExcitedStates(lattice, atom_core), n(40), k(7), rmax(50.) {}
    virtual ~BSplineBasis() {}

    /** Create virtual states above the core.
        num_states_per_l is simply a vector to say how many states
        should be included above the core for each angular momentum, L.
        eg: {3, 2, 1} would make 3 s-wave, 2 p-wave, 1 d-wave, etc.
        If a state already exists, then this routine just skips it.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l);

    /** Update all of the excited states because the core has changed. */
    virtual void Update();

    /** Set the number of splines and their order, as well as the cavity radius. */
    inline void SetParameters(int num_splines, int order, double r = 50.);

protected:
    std::vector<unsigned int> NumStatesPerL;

    int n;     // Number of splines
    int k;      // Order of splines. Maximum of k is 15.
    double rmax;
};

inline void BSplineBasis::SetParameters(int num_splines, int order, double r)
{   n = num_splines;
    k = order;
    rmax = r;
}

#endif
