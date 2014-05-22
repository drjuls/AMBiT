#include "Include.h"
#include "ExpLattice.h"

ExpLattice::ExpLattice(const ExpLattice& other):
    Lattice(other)
{}

ExpLattice::ExpLattice(unsigned int numpoints, double r_min, double H):
    Lattice()
{
    beta = 0.;
    num_points = numpoints;
    original_size = numpoints;
    rmin = r_min;
    h = H;

    r.resize(num_points);
    dr.resize(num_points);

    for(unsigned int i=0; i<num_points; i++)
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
    double r = rmin * exp(x);

    return r;
}

unsigned int ExpLattice::real_to_lattice(double r_point) const
{
    if(r_point <= rmin)
        return 0;

    // x = i h = ln(r/rmin)
    double i_point = log(r_point/rmin)/h;
    return (unsigned int)(i_point + 0.5);
}

double ExpLattice::calculate_dr(double r_point) const
{
    return r_point*h;
}
