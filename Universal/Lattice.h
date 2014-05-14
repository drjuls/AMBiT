#ifndef LATTICE_H
#define LATTICE_H

#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

class Lattice
{
    /** The lattice class provides a conversion between a lattice with even spacing (x),
        and a "real" space which may not (r).
            x = r + beta*ln(r/rmin)
        Lattice equality is defined by type, and in this case, beta, h, and rmin.
        None of the public functions can change these (constructor aside), so it is safe
        to pass non-const lattice pointers around, which is useful since then the lattice
        can expand if necessary.
     */
public:
    Lattice(const Lattice& other);
    Lattice(unsigned int numpoints, double r_min, double r_max);
    Lattice(FILE* binary_infile);
    virtual ~Lattice(void);

    unsigned int Size() const { return NumPoints; }
    double MaxRealDistance() const { return r[NumPoints-1]; }

    double R(unsigned int i);
    double dR(unsigned int i);

    /** Return all points of R, from 0 to size()-1 */
    const double* R() const { return r.data(); }
    const double* dR() const { return dr.data(); }

    /** Return all points of R^k, from 0 to size()-1.
        PRE: k > 0
      */
    const double* Rpower(unsigned int k);

    double H() const { return h; }

    /** Return the first lattice point greater than or equal to "r_point".
        If no such point exists, create it.
     */
    virtual unsigned int real_to_lattice(double r_point);

    /** Equality does not check size of lattice. */
    virtual bool operator==(const Lattice& other) const;

    virtual void Write(FILE* binary_outfile) const;

protected:
    Lattice(): NumPoints(0), beta(0.), h(0.), rmin(0.) {}

    /** Calculate the value that r[i] should be. */
    virtual double lattice_to_real(unsigned int i) const;

    /** Calculate the lattice spacing at a point. */
    virtual double calculate_dr(double r_point) const;

    /** Calculate R^power and store, for all powers up to k.
        PRE: k >= 2
     */
    const double* Calculate_Rpower(unsigned int k);
    
    /** Resizes the lattice such that NumPoints > min_size. */
    virtual void resize(double min_size);

protected:
    unsigned int NumPoints;
    std::vector<double> r, dr;

    // r_power[k-2] = R^k, defined for k >= 2.
    std::vector<std::vector<double>> r_power;

    double beta, h, rmin;
};

typedef boost::shared_ptr<Lattice> pLattice;
typedef boost::shared_ptr<const Lattice> pLatticeConst;

inline double Lattice::R(unsigned int i)
{   if(i >= NumPoints)
        resize(i);
    return r[i];
}

inline double Lattice::dR(unsigned int i)
{   if(i >= NumPoints)
        resize(i);
    return dr[i];
}

inline const double* Lattice::Rpower(unsigned int k)
{
    if(k == 0)
        return NULL;
    if(k == 1)
        return r.data();
    if(k <= r_power.size()+1)
        return r_power[k - 2].data();
    else
        return Calculate_Rpower(k);
}

#endif
