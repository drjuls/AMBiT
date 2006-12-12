#ifndef EXP_LATTICE_H
#define EXP_LATTICE_H

#include "Lattice.h"

/**
 * The lattice class provides a conversion between a lattice with even spacing (x),
 * and a "real" space which may not (r).
 *   x = r + beta*ln(r/rmin)
 */
class ExpLattice : public Lattice
{
public:
    ExpLattice(unsigned int numpoints, double r_min, double H);
    virtual ~ExpLattice(void) {}

    /** Return the first lattice point greater than or equal to "r_point".
        If no such point exists, create it.
     */
    //virtual unsigned int real_to_lattice(double r_point);

protected:
    /** Calculate the value that r[i] should be. */
    virtual double lattice_to_real(unsigned int i) const;

    /** Calculate the lattice spacing at a point. */
    virtual double calculate_dr(double r_point) const;
};

#endif
