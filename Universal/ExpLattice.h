#ifndef EXP_LATTICE_H
#define EXP_LATTICE_H

#include "Lattice.h"

class ExpLattice : public Lattice
{
    /** Simple exponential conversion between a lattice with even spacing (x),
        and a "real" space (r):
            r = rmin * exp(x)
            x = ln(r/rmin)
        beta (from Lattice class) is not used here.
     */
public:
    ExpLattice(const ExpLattice& other);
    ExpLattice(unsigned int numpoints, double r_min, double H);
    virtual ~ExpLattice() {}

    /** Equality does not check size of lattice. */
    virtual bool operator==(const ExpLattice& other) const;
    virtual bool operator==(const Lattice& other) const override { return false; }

    /** Calculate the value that r[i] should be. */
    virtual double lattice_to_real(unsigned int i) const override;

    /** Calculate the closest point i to r_point. */
    virtual unsigned int real_to_lattice(double r_point) const override;

    /** Calculate the lattice spacing at a point. */
    virtual double calculate_dr(double r_point) const override;
};

#endif
