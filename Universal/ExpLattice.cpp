#include "Include.h"
#include "ExpLattice.h"

ExpLattice::ExpLattice(const ExpLattice& other):
    Lattice(other)
{}

ExpLattice::ExpLattice(unsigned int numpoints, double r_min, double H):
    Lattice()
{
    beta = 0.;
    NumPoints = numpoints;
    rmin = r_min;
    h = H;

    r.resize(numpoints);
    dr.resize(numpoints);

    for(unsigned int i=0; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }
}

bool ExpLattice::operator==(const ExpLattice& other) const
{
    return ((h == other.h) && (rmin == other.rmin));
}

double ExpLattice::lattice_to_real(unsigned int i) const
{
    double x = h*double(i);
    double r;

    r = rmin * exp(x);

    return r;
}

double ExpLattice::calculate_dr(double r_point) const
{
    return r_point*h;
}
