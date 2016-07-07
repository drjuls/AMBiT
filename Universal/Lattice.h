#ifndef LATTICE_H
#define LATTICE_H

#include <string>
#include <vector>
#include <list>
#include <memory>
#include <boost/weak_ptr.hpp>

class LatticeObserver;

/** The lattice class provides a conversion between a lattice with even spacing (x),
    and a "real" space which may not (r).
        x = r + beta*ln(r/rmin)
    Lattice equality is defined by type, and in this case, beta, h, and rmin.
    None of the public functions can change these (constructor aside), so it is safe
    to pass non-const lattice pointers around, which is useful since then the lattice
    can expand if necessary.
    When the lattice size changes the lattice notifies observers that have registered
    using Attach().
 */
class Lattice
{
public:
    Lattice(unsigned int numpoints, double r_min, double r_max);
    Lattice(FILE* binary_infile);
    virtual ~Lattice() {}

    /** Equality does not check size of lattice. */
    virtual bool operator==(const Lattice& other) const;

    /** Current size of lattice. */
    unsigned int size() const { return num_points; }

    /** Resize the lattice to max(new_size, original_size).
        That is, the lattice size is never smaller than original_size.
        Notifies observers if size changes.
        Returns new lattice size.
     */
    unsigned int resize(unsigned int new_size);

    /** Resize the lattice to radius r_max (however the lattice size is never smaller than original_size).
        Notifies observers if size changes.
        Returns new lattice size.
     */
    unsigned int resize_to_r(double r_max);

    double MaxRealDistance() const { return r[num_points-1]; }

    /** PRE: i < size() */
    inline double R(unsigned int i) const { return r[i]; }

    /** PRE: i < size() */
    inline double dR(unsigned int i) const { return dr[i]; }

    /** Return all points of R, from 0 to size()-1 */
    const double* R() const { return r.data(); }
    const double* dR() const { return dr.data(); }

    /** Return all points of R^k, from 0 to size()-1.
        PRE: k > 0
      */
    const double* Rpower(unsigned int k);

    /** Add to observer list. */
    void Subscribe(LatticeObserver* observer);

    /** Find and remove from observer list. */
    void Unsubscribe(LatticeObserver* observer);

    /** Calculate the value that r[i] should be. */
    virtual double lattice_to_real(unsigned int i) const;

    /** Calculate the closest point i to r_point. */
    virtual unsigned int real_to_lattice(double r_point) const;

    /** Calculate the lattice spacing at a point. */
    virtual double calculate_dr(double r_point) const;
    
    /** Return spacing of transformed lattice x = rmin + i h,
        where i is the lattice point.
     */
    double H() const { return h; }

    virtual void Write(FILE* binary_outfile) const;

protected:
    Lattice(): original_size(0), num_points(0), beta(0.), h(0.), rmin(0.) {}

    /** Calculate R^power and store, for all powers up to k.
        PRE: k >= 2
     */
    const double* Calculate_Rpower(unsigned int k);

    void Notify();

protected:
    unsigned int original_size;
    double beta, h, rmin;

    // Current size and points
    unsigned int num_points;
    std::vector<double> r, dr;

    // r_power[k-2] = R^k, defined for k >= 2.
    std::vector<std::vector<double>> r_power;

    // Observers
    std::list<LatticeObserver*> observers;
};

typedef std::shared_ptr<Lattice> pLattice;
typedef std::shared_ptr<const Lattice> pLatticeConst;

/** LatticeObserver is an abstract interface for any class that needs to observe changes in the size of a lattice.
    The lattice that is under observation is passed in the Alert() function, since in some cases the observer
    may not need to keep a pointer to the lattice, but just be registered with it.
 */
class LatticeObserver
{
public:
    LatticeObserver(pLattice lattice): lattice(lattice)
    {   lattice->Subscribe(this);
    }
    /** Need to include non-default copy constructor: the new LatticeObserver needs to subscribe. */
    LatticeObserver(const LatticeObserver& other):
        LatticeObserver(other.lattice)
    {}
    virtual ~LatticeObserver()
    {   lattice->Unsubscribe(this);
    }

    /** lattice has changed size; react accordingly. */
    virtual void Alert() = 0;

    /** All LatticeObservers can share the lattice, even when const.
        The lattice can be considered a "mutable" member of classes since it can expand and in any case
        all observers will be alerted if it changes.
        Hence it is reasonable to return pLattice rather than pLatticeConst.
     */
    virtual pLattice GetLattice() const { return lattice; }

protected:
    pLattice lattice;
};

inline void Lattice::Notify()
{
    for(auto& observer: observers)
        observer->Alert();
}

inline const double* Lattice::Rpower(unsigned int k)
{
    if(k == 1)
        return r.data();
    if(k <= r_power.size()+1)
        return r_power[k - 2].data();
    else
        return Calculate_Rpower(k);
}

#endif
