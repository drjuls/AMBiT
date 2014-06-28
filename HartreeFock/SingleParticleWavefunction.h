#ifndef SINGLE_PARTICLE_WAVEFUNCTION_H
#define SINGLE_PARTICLE_WAVEFUNCTION_H

#include "Universal/SpinorFunction.h"
#include "Universal/Lattice.h"
#include <boost/shared_ptr.hpp>
#include <string>
#include <math.h>
#include <stdlib.h>

class OrbitalInfo;

/** Single particle wavefunction is for solutions of a single-particle Hamiltonian.
    That is, it is a SpinorFunction with an energy eigenvalue E:
        \f[ h|\psi> = E|\psi> \f]
    It also contains a principal quantum number, which may or may not correspond
    to the spectroscopic pqn.
 */
class SingleParticleWavefunction : public SpinorFunction
{
public:
    SingleParticleWavefunction(const OrbitalInfo& info);
    SingleParticleWavefunction(int kappa, int pqn = 0, double energy = 0.0, unsigned int size = 0);
    SingleParticleWavefunction(const SingleParticleWavefunction& other);
    SingleParticleWavefunction(SingleParticleWavefunction&& other);
    virtual ~SingleParticleWavefunction() {}

    const SingleParticleWavefunction& operator=(const SingleParticleWavefunction& other);
    SingleParticleWavefunction& operator=(SingleParticleWavefunction&& other);

    virtual double Energy() const;

    /** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
    virtual double Nu() const;
    virtual int PQN() const;

    virtual void SetEnergy(double Energy);
    virtual void SetNu(double Nu);
    virtual void SetPQN(int PQN);

    virtual std::string Name() const;

    const SingleParticleWavefunction& operator*=(double scale_factor);
    SingleParticleWavefunction operator*(double scale_factor) const;
    const SingleParticleWavefunction& operator+=(const SingleParticleWavefunction& other);
    const SingleParticleWavefunction& operator-=(const SingleParticleWavefunction& other);
    SingleParticleWavefunction operator+(const SingleParticleWavefunction& other) const;
    SingleParticleWavefunction operator-(const SingleParticleWavefunction& other) const;
    const SingleParticleWavefunction& operator*=(const RadialFunction& chi);
    SingleParticleWavefunction operator*(const RadialFunction& chi) const;

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

    /** Print state to file, optionally printing lattice. Return success. */
    bool Print(pLattice lattice = pLattice()) const;
    bool Print(const std::string& filename, pLattice lattice = pLattice()) const;
    bool Print(FILE* fp, pLattice lattice = pLattice()) const;

protected:
    int pqn;
    double energy;
};

typedef boost::shared_ptr<SingleParticleWavefunction> pSingleParticleWavefunction;
typedef boost::shared_ptr<const SingleParticleWavefunction> pSingleParticleWavefunctionConst;

typedef SingleParticleWavefunction ContinuumWave;
typedef pSingleParticleWavefunction pContinuumWave;
typedef pSingleParticleWavefunctionConst pContinuumWaveConst;

#endif
