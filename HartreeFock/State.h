#ifndef STATE_H
#define STATE_H

//#define USE_ALT_STATE_NOTATION

#include "Universal/CoupledFunction.h"
#include "Universal/Lattice.h"
#include <string>
#include <math.h>

class State : public CoupledFunction
{
    /** Electron state. Psi = (        f/r. Y_jlm)
                              (i.alpha.g/r.~Y_jlm)
     */
public:
    State();
    State(int Kappa);
    State(const State& other);
    virtual ~State(void) {}

    virtual double Energy() const = 0;
    /** Nu is the effective principal quantum number */
    double Nu() const { return nu; }
    
    int Kappa() const { return kappa; }
    unsigned int L() const;
    double J() const;
    unsigned int TwoJ() const;
    virtual std::string Name() const = 0;

    inline void SetKappa(int kappa) { this->kappa = kappa; }
    inline void SetNu(double nu) { this->nu = nu; }

    const State& operator=(const State& other)
    {   CoupledFunction::operator=(other);
        kappa = other.kappa;
        nu = other.nu;
        return *this;
    }

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const { CoupledFunction::Write(fp); }

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

    /** Get overlap of this state with another. */
    double Overlap(const State& other, const Lattice* lattice) const;

    /** Print state to file, optionally printing lattice. Return success. */
    bool Print(Lattice* lattice = NULL) const;
    bool Print(const std::string& filename, Lattice* lattice = NULL) const;
    bool Print(FILE* fp, Lattice* lattice = NULL) const;

protected:
    int kappa;
    double nu;          // effective principal quantum number - determines energy
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

inline unsigned int State::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

#endif
