#ifndef ORBITAL_H
#define ORBITAL_H

#include "State.h"

class Orbital : public State  
{
public:
    Orbital(): State() {}
    Orbital(unsigned int PrincipalQN, int Kappa);
    Orbital(const Orbital& other);
    virtual ~Orbital(void) {}

    virtual double Energy() const;
    double Norm(const Lattice* lattice) const;

    inline void SetRequiredPrincipalQN(unsigned int pqn) { this->pqn = pqn; }
    inline void SetOccupancy(double number_electrons);
    inline void SetEnergy(double energy);

    virtual std::string Name() const;
    unsigned int RequiredPQN() const { return pqn; }
    double Occupancy() const { return occupancy; }

    /** Check that the ratio
          f[Size()-1]/f_max < tolerance
        and
          f[Size()-2]/f_max > tolerance
        and if they are not, enlarge or reduce the size of the wavefunctions.
        Return true if the size was correct, false if there was a change.
      */
    bool CheckSize(Lattice* lattice, double tolerance);

    /** Scale the state so that it is normalised to "norm". */
    void ReNormalise(const Lattice* lattice, double norm = 1.);

    /** Count the number of nodes of the wavefunction. */
    unsigned int NumNodes() const;

    const Orbital& operator=(const Orbital& other)
    {   State::operator=(other);
        pqn = other.pqn;
        occupancy = other.occupancy;
        return *this;
    }

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;

    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);

protected:
    unsigned int pqn;   // principal quantum number (eg: the 4 in "4s")
    double occupancy;
};

inline void Orbital::SetOccupancy(double number_electrons)
{   if(number_electrons <= (double)(2*abs(kappa)))
        occupancy = number_electrons;
}

inline void Orbital::SetEnergy(double energy)
{
    if (energy < 0.)
        nu = 1./sqrt(-2.*energy);
    else if(energy > 0.)
        nu = -1./sqrt(2.*energy);
    else
        nu = 0.;
}

#endif

