#include "Include.h"
#include "Lattice.h"
#include <fstream>

Lattice::Lattice(const Lattice& other):
    r(other.dr), dr(other.dr), beta(other.beta), NumPoints(other.NumPoints), rmin(other.rmin), h(other.h), r_power(other.r_power)
{}

Lattice::Lattice(unsigned int numpoints, double r_min, double r_max):
    beta(4.0), NumPoints(numpoints), rmin(r_min)
{
    r.resize(NumPoints);
    dr.resize(NumPoints);
    
    h =(r_max - rmin + beta*log(r_max/rmin))/(NumPoints-1);
    
    for(unsigned int i=0; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }
}

Lattice::Lattice(FILE* binary_infile)
{
    fread(&beta, sizeof(double), 1, binary_infile);
    fread(&h, sizeof(double), 1, binary_infile);
    fread(&rmin, sizeof(double), 1, binary_infile);
    fread(&NumPoints, sizeof(unsigned int), 1, binary_infile);

    r.resize(NumPoints);
    dr.resize(NumPoints);

    fread(r.data(), sizeof(double), NumPoints, binary_infile);
    fread(dr.data(), sizeof(double), NumPoints, binary_infile);
}

Lattice::~Lattice(void)
{}

void Lattice::resize(double min_size)
{
    if(min_size <= r[NumPoints-1])
        return;

    unsigned int old_size = NumPoints;

    do
    {   NumPoints *= 2;
    }while(lattice_to_real(NumPoints - 1) <= min_size);

    r.resize(NumPoints);
    dr.resize(NumPoints);

    for(unsigned int i=old_size; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }

    // Do powers of r
    for(unsigned int k = 0; k < r_power.size(); k++)
    {
        std::vector<double>& previous = (k == 0)? r : r_power[k-1];
        std::vector<double>& current = r_power[k];

        current.resize(NumPoints);

        for(unsigned int i=old_size; i<NumPoints; i++)
            current[i] = previous[i] * r[i];
    }
}

double Lattice::lattice_to_real(unsigned int i) const
{
    // Lattice coordinate x = rmin + ih,
    //      y = x + beta * log(rmin)
    double y = rmin + beta*log(rmin) + h*double(i);
    double r, rold;

    if(beta <= y)
        r = y;
    else
        r = exp(y/beta);

    // Solve using Newton's method f(r) == 0, where
    //    f(r) = r + beta*log(r/rmin) - x
    //         = r + beta*log(r) - y
    do
    {   rold = r;
        r = r - (r + beta*log(r) - y)/(1. + beta/r);
        if(r <= 0.0)
            r = rold/2.0;
    } while(fabs(r - rold)/r > 1.e-13);

    return r;
}

unsigned int Lattice::real_to_lattice(double r_point)
{
    if(r_point > MaxRealDistance())
        resize(r_point);
    else if (r_point <= r[0])
        return 0;

    unsigned int i_min = 0;
    unsigned int i_max = NumPoints;
    unsigned int i_mid;

    // Bisection search
    while(i_max - i_min > 1)
    {   i_mid = (i_min + i_max)/2;
        (r_point >= r[i_mid]) ? i_min = i_mid : i_max = i_mid;
    }

    if(r_point == r[i_min])
        return i_min;
    else
        return i_min + 1;
}

bool Lattice::operator==(const Lattice& other) const
{
    return ((beta == other.beta) && (h == other.h) && (rmin == other.rmin));
}

double Lattice::calculate_dr(double r_point) const
{
    return r_point/(beta + r_point) * h;
}

const double* Lattice::Calculate_Rpower(unsigned int k)
{
    unsigned int kminustwo = k - 2;

    unsigned int old_size = r_power.size();
    if(kminustwo >= old_size)
        r_power.resize(kminustwo+1);

    // Create all powers up to and including k
    for(unsigned int new_k = old_size; new_k <= kminustwo; new_k++)
    {
        std::vector<double>& previous = (new_k == 0)? r : r_power[new_k-1];
        std::vector<double>& current = r_power[new_k];

        current.resize(NumPoints);

        for(unsigned int i = 0; i < NumPoints; i++)
            current[i] = previous[i] * r[i];
    }

    return r_power[kminustwo].data();
}

void Lattice::Write(FILE* binary_outfile) const
{
    fwrite(&beta, sizeof(double), 1, binary_outfile);
    fwrite(&h, sizeof(double), 1, binary_outfile);
    fwrite(&rmin, sizeof(double), 1, binary_outfile);
    fwrite(&NumPoints, sizeof(unsigned int), 1, binary_outfile);

    fwrite(r.data(), sizeof(double), NumPoints, binary_outfile);
    fwrite(dr.data(), sizeof(double), NumPoints, binary_outfile);
}
