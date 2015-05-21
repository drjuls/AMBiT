#ifndef ORBITAL_INFO_H
#define ORBITAL_INFO_H

#include "Orbital.h"
#include <string>

class OrbitalInfo
{
    /** Single particle state information.
        Stores pqn and kappa (L). Has an inbuilt ordering.
     */
public:
    OrbitalInfo(int pqn, int kappa);
    OrbitalInfo(int pqn, int two_j, Parity p);
    OrbitalInfo(pOrbitalConst s);
    OrbitalInfo(const Orbital* s);
    OrbitalInfo(const OrbitalInfo& other);
    virtual ~OrbitalInfo(void) {}

    const OrbitalInfo& operator=(const OrbitalInfo& other);

    virtual bool operator<(const OrbitalInfo& other) const;
    virtual bool operator==(const OrbitalInfo& other) const;
    virtual bool operator!=(const OrbitalInfo& other) const;

    inline int PQN() const { return pqn; }
    inline int Kappa() const { return kappa; }
    inline int L() const;
    inline double J() const;
    inline int TwoJ() const;
    inline int L_Prime() const; //!< Return the value of L for the lower component of the wavefunction.
    inline Parity GetParity() const;

    virtual int MaxNumElectrons() const { return 2*abs(kappa); }
    virtual std::string Name() const;

protected:
    int pqn;
    int kappa;
};

inline OrbitalInfo::OrbitalInfo(int pqn, int two_j, Parity p):
    pqn(pqn)
{
    kappa = (two_j + 1)/2;
    int L = kappa;
    // If inconsistent, swap kappa
    if(!(L%2 == 0 && p == Parity::even) && !(L%2 == 1 && p == Parity::odd))
    {   kappa = -kappa;
    }
}

inline OrbitalInfo::OrbitalInfo(pOrbitalConst s)
{
    pqn = s->PQN();
    kappa = s->Kappa();
}

inline OrbitalInfo::OrbitalInfo(const Orbital* s)
{
    pqn = s->PQN();
    kappa = s->Kappa();
}

inline OrbitalInfo::OrbitalInfo(int principal_qn, int kap):
    pqn(principal_qn), kappa(kap)
{}

inline OrbitalInfo::OrbitalInfo(const OrbitalInfo& other):
    pqn(other.pqn), kappa(other.kappa)
{}

inline int OrbitalInfo::L() const
{
    if(kappa > 0)
        return kappa;
    else
        return (-kappa-1);
}

inline double OrbitalInfo::J() const
{
    return fabs(double(kappa)) - 0.5;
}

// Returns the value of L for the lower component of the wavefunction
inline int OrbitalInfo::L_Prime() const
{
    // The lower component has L corresponding to -Kappa
    if(kappa < 0)
        return -kappa;
    else
        return kappa-1;
}

inline int OrbitalInfo::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

inline Parity OrbitalInfo::GetParity() const
{
    if(L()%2 == 0)
        return Parity::even;
    else
        return Parity::odd;
}

#endif
