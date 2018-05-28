#include "Include.h"
#include "Lattice.h"
#include <fstream>

namespace Ambit
{
Lattice::Lattice(unsigned int numpoints, double r_min, double r_max):
    beta(4.0), num_points(numpoints), original_size(numpoints), rmin(r_min)
{
    r.resize(num_points);
    dr.resize(num_points);
    
    h =(r_max - rmin + beta*log(r_max/rmin))/(num_points-1);
    
    for(unsigned int i=0; i<num_points; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }
}

Lattice::Lattice(FILE* binary_infile)
{
    fread(&beta, sizeof(double), 1, binary_infile);
    fread(&h, sizeof(double), 1, binary_infile);
    fread(&rmin, sizeof(double), 1, binary_infile);
    fread(&num_points, sizeof(unsigned int), 1, binary_infile);

    r.resize(num_points);
    dr.resize(num_points);

    fread(r.data(), sizeof(double), num_points, binary_infile);
    fread(dr.data(), sizeof(double), num_points, binary_infile);

    original_size = num_points;
}

unsigned int Lattice::resize(unsigned int new_size)
{
    unsigned int old_size = size();
    new_size = mmax(new_size, original_size);

    if(old_size != new_size)
    {
        r.resize(new_size);
        dr.resize(new_size);
        for(auto& r_k: r_power)
            r_k.resize(new_size);

        for(unsigned int i = old_size; i < new_size; i++)
        {   r[i] = lattice_to_real(i);
            dr[i] = calculate_dr(r[i]);
        }

        // Do powers of r
        for(unsigned int k = 0; k < r_power.size(); k++)
        {
            std::vector<double>& previous = (k == 0)? r : r_power[k-1];
            std::vector<double>& current = r_power[k];

            current.resize(new_size);

            for(unsigned int i = old_size; i < new_size; i++)
                current[i] = previous[i] * r[i];
        }

        num_points = new_size;
        Notify();
    }

    return num_points;
}

unsigned int Lattice::resize_to_r(double r_max)
{
    return resize(real_to_lattice(r_max));
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

unsigned int Lattice::real_to_lattice(double r_point) const
{
    if(r_point <= rmin)
        return 0;

    double i_point = (r_point - rmin + beta * log(r_point/rmin))/h;
    return (unsigned int)(i_point + 0.5);
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

        current.resize(num_points);

        for(unsigned int i = 0; i < num_points; i++)
            current[i] = previous[i] * r[i];
    }

    return r_power[kminustwo].data();
}

void Lattice::Subscribe(LatticeObserver* observer)
{
    // Often (but not always) last to subscribe is first to unsubscribe.
    // So we add to front so that it is found more quickly in list.
    observers.push_front(observer);
}

void Lattice::Unsubscribe(LatticeObserver* observer)
{
    auto it = observers.begin();
    while(it != observers.end())
    {
        if(*it == observer)
        {   observers.erase(it);
            break;
        }
        else
            it++;
    }
}

void Lattice::Write(FILE* binary_outfile) const
{
    file_err_handler->fwrite(&beta, sizeof(double), 1, binary_outfile);
    file_err_handler->fwrite(&h, sizeof(double), 1, binary_outfile);
    file_err_handler->fwrite(&rmin, sizeof(double), 1, binary_outfile);
    file_err_handler->fwrite(&num_points, sizeof(unsigned int), 1, binary_outfile);

    file_err_handler->fwrite(r.data(), sizeof(double), num_points, binary_outfile);
    file_err_handler->fwrite(dr.data(), sizeof(double), num_points, binary_outfile);
}
}
