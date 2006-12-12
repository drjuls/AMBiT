#include "Include.h"
#include "ExpLattice.h"

ExpLattice::ExpLattice(unsigned int numpoints, double r_min, double H) :
    Lattice()
{
    //Lattice(), NumPoints(numpoints), rmin(r_min), h(H)
    NumPoints = numpoints;
    rmin = r_min;
    h = H;

    r = (double*)malloc(NumPoints * sizeof(double));
    dr = (double*)malloc(NumPoints * sizeof(double));

    for(unsigned int i=0; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }
}

/*
unsigned int ExpLattice::real_to_lattice(double r_point)
{
    while(r_point > r[NumPoints-1])
    {   ReSize(NumPoints*2-1);
    }
    
    double x = log(r_point/rmin);
    unsigned int i = (unsigned int)(x/h);
    
    if(i>0)
        i--;

    while(r_point > r[i])
       i++;

    return i;
}
*/

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
