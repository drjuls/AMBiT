#ifndef STATE_H
#define STATE_H

#include "CoupledFunction.h"
#include "Universal/Lattice.h"
#include <string>
#include <math.h>

/**
 * Electron state. Psi = (        f/r. Y_jlm)
 *                       (i.alpha.g/r.~Y_jlm)
 */
class State : public CoupledFunction
{
public:
    State(Lattice* lat, unsigned int num_points = 0);
    State(Lattice* lat, int Kappa, unsigned int num_points = 0);
    State(const State& other);
    virtual ~State(void) {}

    virtual double Energy() const = 0;
    /** Nu is the effective principal quantum number */
    double Nu() const { return nu; }
    unsigned int NumZeroes() const;
    
    int Kappa() const { return kappa; }
    unsigned int L() const;
    double J() const;
    virtual std::string Name() const = 0;

    inline void SetKappa(int kappa) { this->kappa = kappa; }
    inline void SetNu(double nu) { this->nu = nu; }

    const State& operator=(const State& other)
    {   CoupledFunction::operator=(other);
        kappa = other.kappa;
        nu = other.nu;

        lattice = other.lattice;
        return *this;
    }

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const { CoupledFunction::Write(fp); }

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

    /** Get overlap of this state with another. */
    double Overlap(const State& other) const;
protected:
    int kappa;
    double nu;          // effective principal quantum number - determines energy

    Lattice* lattice;
};

inline unsigned int State::L() const
{   if (kappa > 0)
        return (unsigned int)(kappa);
    else
        return (unsigned int)(-kappa-1);
}

inline double State::J() const
{   return (double)abs(kappa) - 0.5;
}

#endif
