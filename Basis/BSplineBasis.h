#ifndef BSPLINE_BASIS_H
#define BSPLINE_BASIS_H

#include "ExcitedStates.h"

class BSplineBasis : public ExcitedStates
{
protected:
    enum SplineType {NotreDame, Reno, Vanderbilt};

public:
    BSplineBasis(Lattice* lattice, Core* atom_core):
        ExcitedStates(lattice, atom_core), n(40), k(7), rmax(50.), spline_type(Reno) {}
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

protected:
    class BSpline : public Orbital
    {
        /** Not really an Orbital. This is abuse of the class. But it should work well for us. */
    public:
        BSpline(int Kappa, unsigned int size = 0): Orbital()
        {
            SetKappa(Kappa);
            ReSize(size);
        }
    };

    SplineType spline_type;
};

inline void BSplineBasis::SetParameters(int num_splines, int order, double r)
{   n = num_splines;
    k = order;
    rmax = r;
}

#endif
