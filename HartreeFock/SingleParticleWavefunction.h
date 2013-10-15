#ifndef SINGLE_PARTICLE_WAVEFUNCTION_H
#define SINGLE_PARTICLE_WAVEFUNCTION_H

//#define USE_ALT_STATE_NOTATION

#include "Universal/SpinorFunction.h"
#include "Universal/Lattice.h"
#include <string>
#include <math.h>
#include <stdlib.h>

class OrbitalInfo;

/** Single particle wavefunction is for solutions of a single-particle Hamiltonian.
    That is, it is a SpinorFunction with an energy eigenvalue E:
        h|psi> = E|psi>
    It also contains a principal quantum number, which may or may not correspond
    to the spectroscopic pqn.
 */
class SingleParticleWavefunction : public SpinorFunction
{
public:
    SingleParticleWavefunction(const OrbitalInfo& info);
    SingleParticleWavefunction(int kappa, double energy = 0.0, unsigned int pqn = 0, unsigned int size = 0);
    SingleParticleWavefunction(const SingleParticleWavefunction& other);
    virtual ~SingleParticleWavefunction() {}

    virtual double GetEnergy() const;

    /** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
    virtual double GetNu() const;
    virtual unsigned int GetPQN() const;

    virtual void SetEnergy(double Energy);
    virtual void SetNu(double Nu);
    virtual void SetPQN(int PQN);

    virtual std::string Name() const;

    const SingleParticleWavefunction& operator=(const SingleParticleWavefunction& other);
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

    double Overlap(const SingleParticleWavefunction& other, const Lattice* lattice) const; // Deprecate

    /** Print state to file, optionally printing lattice. Return success. */
    bool Print(Lattice* lattice = NULL) const;
    bool Print(const std::string& filename, Lattice* lattice = NULL) const;
    bool Print(FILE* fp, Lattice* lattice = NULL) const;

protected:
    int pqn;
    double energy;
};

typedef SingleParticleWavefunction ContinuumWave;

#endif
