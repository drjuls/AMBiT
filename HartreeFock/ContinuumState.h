#ifndef CONTINUUM_STATE_H
#define CONTINUUM_STATE_H

#include "State.h"

class ContinuumState: public State
{
public:
    ContinuumState(Lattice* lat, unsigned int num_points = 0);
    ContinuumState(Lattice* lat, double Nu, int Kappa, unsigned int num_points = 0);
    ContinuumState(const ContinuumState& other);
    ContinuumState(double energy, int Kappa);
    virtual ~ContinuumState(void) {}

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

public:
    virtual double Energy() const { return 0.5/(nu*nu); }
    virtual std::string Name() const;
};

#endif
