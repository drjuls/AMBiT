#ifndef ORBITAL_H
#define ORBITAL_H

#include "SingleParticleWavefunction.h"
#include <boost/shared_ptr.hpp>

/** Orbital is really no different to SingleParticleWavefunction,
    however it is explicitly bounded and hence normalisable, integrable, etc.
 */
class Orbital : public SingleParticleWavefunction  
{
public:
    Orbital(const OrbitalInfo& info);
    Orbital(int kappa, double energy = 0., unsigned int pqn = 0, unsigned int size = 0);
    Orbital(const Orbital& other);
    virtual ~Orbital() {}

    double Norm(const Lattice* lattice) const; // Deprecate

    inline void SetOccupancy(double number_electrons); // Deprecate
    double Occupancy() const { return occupancy; } // Deprecate

    virtual std::string Name() const;

    const Orbital& operator=(const Orbital& other);
    const Orbital& operator*=(double scale_factor);
    Orbital operator*(double scale_factor) const;
    const Orbital& operator+=(const Orbital& other);
    const Orbital& operator-=(const Orbital& other);
    Orbital operator+(const Orbital& other) const;
    Orbital operator-(const Orbital& other) const;
    const Orbital& operator*=(const RadialFunction& chi);
    Orbital operator*(const RadialFunction& chi) const;

    /** Check that the ratio
          f[Size()-1]/f_max < tolerance
        and
          f[Size()-2]/f_max > tolerance
        and if they are not, enlarge or reduce the size of the wavefunctions.
        Return true if the size was correct, false if there was a change.
      */
    bool CheckSize(Lattice* lattice, double tolerance);

    /** Scale the state so that it is normalised to "norm". */
    void ReNormalise(const Lattice* lattice, double norm = 1.); // Deprecate

    /** Count the number of nodes of the wavefunction. */
    unsigned int NumNodes() const;

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

protected:
    double occupancy;
};

typedef boost::shared_ptr<Orbital> pOrbital;
typedef boost::shared_ptr<const Orbital> pOrbitalConst;

inline void Orbital::SetOccupancy(double number_electrons)
{   if(number_electrons <= (double)(2*abs(kappa)))
        occupancy = number_electrons;
}

#endif

