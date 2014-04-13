#ifndef ORBITAL_H
#define ORBITAL_H

#include "SingleParticleWavefunction.h"
#include "OpIntegrator.h"
#include <boost/shared_ptr.hpp>

/** Orbital is really no different to SingleParticleWavefunction,
    however it is explicitly bounded and hence normalisable, integrable, etc.
 */
class Orbital : public SingleParticleWavefunction  
{
public:
    Orbital(const OrbitalInfo& info): SingleParticleWavefunction(info) {}
    Orbital(int kappa, unsigned int pqn = 0, double energy = 0.0, unsigned int size = 0):
        SingleParticleWavefunction(kappa, pqn, energy, size) {}
    Orbital(const Orbital& other): SingleParticleWavefunction(other) {}
    Orbital(Orbital&& other): SingleParticleWavefunction(other) {}
    virtual ~Orbital() {}

    const Orbital& operator=(const Orbital& other);
    Orbital& operator=(Orbital&& other);

    virtual std::string Name() const;

    const Orbital& operator*=(double scale_factor);
    Orbital operator*(double scale_factor) const;
    const Orbital& operator+=(const Orbital& other);
    const Orbital& operator-=(const Orbital& other);
    Orbital operator+(const Orbital& other) const;
    Orbital operator-(const Orbital& other) const;
    const Orbital& operator*=(const RadialFunction& chi);
    Orbital operator*(const RadialFunction& chi) const;

    /** Check that the ratio
          f[size()-1]/f_max < tolerance
        and
          f[size()-2]/f_max > tolerance
        and if they are not, enlarge or reduce the size of the wavefunctions.
        Return true if the size was correct, false if there was a change.
      */
    bool CheckSize(pLattice lattice, double tolerance);

    /** Get current normalisation of the orbital. */
    double Norm(pOPIntegrator integrator) const;

    /** Scale the orbital so that it is normalised to "norm". */
    void ReNormalise(pOPIntegrator integrator, double norm = 1.);

    /** Count the number of nodes of the wavefunction. */
    unsigned int NumNodes() const;

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);
};

typedef boost::shared_ptr<Orbital> pOrbital;
typedef boost::shared_ptr<const Orbital> pOrbitalConst;

#endif

