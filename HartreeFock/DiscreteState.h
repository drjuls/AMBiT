#ifndef DISCRETE_STATE_H
#define DISCRETE_STATE_H

#include "State.h"

class DiscreteState : public State  
{
public:
    DiscreteState(): State() {}
    DiscreteState(unsigned int PrincipalQN, int Kappa);
    DiscreteState(const DiscreteState& other);
    virtual ~DiscreteState(void) {}

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
    void ReNormalise(const Lattice* lattice, double norm = 1.);

    const DiscreteState& operator=(const DiscreteState& other)
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

inline void DiscreteState::SetOccupancy(double number_electrons)
{   if(number_electrons <= (double)(2*abs(kappa)))
        occupancy = number_electrons;
}

inline void DiscreteState::SetEnergy(double energy)
{
    if (energy < 0.)
        nu = 1./sqrt(-2.*energy);
    else if(energy > 0.)
        nu = -1./sqrt(2.*energy);
    else
        nu = 0.;
}

#endif

