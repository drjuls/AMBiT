#ifndef CONTINUUM_WAVE_H
#define CONTINUUM_WAVE_H

#include "SingleParticleWavefunction.h"

class ContinuumWave: public SingleParticleWavefunction
{
public:
    ContinuumWave(): SingleParticleWavefunction() {}
    ContinuumWave(const ContinuumWave& other);
    ContinuumWave(double energy, int Kappa);
    virtual ~ContinuumWave(void) {}

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

    virtual void SetEnergy(double energy);

public:
    virtual double Energy() const { return 0.5/(nu*nu); }
    virtual std::string Name() const;
};

#endif
