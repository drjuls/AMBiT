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
    OrbitalInfo(pOrbitalConst s);
    OrbitalInfo(const Orbital* s);
    OrbitalInfo(const OrbitalInfo& other);
    virtual ~OrbitalInfo(void) {}

    virtual bool operator<(const OrbitalInfo& other) const;
    virtual bool operator==(const OrbitalInfo& other) const;
    virtual bool operator!=(const OrbitalInfo& other) const;

    inline int PQN() const { return pqn; }
    inline int Kappa() const { return kappa; }
    inline int L() const;
    inline double J() const;
    inline unsigned int TwoJ() const;
    // Return the value of L for the lower component of the wavefunction
    inline int L_Prime() const;

    inline unsigned int MaxNumElectrons() const { return 2*abs(kappa); }
    virtual std::string Name() const;

protected:
    int pqn;
    int kappa;
};

inline OrbitalInfo::OrbitalInfo(pOrbitalConst s)
{
    pqn = s->GetPQN();
    kappa = s->Kappa();
}

inline OrbitalInfo::OrbitalInfo(const Orbital* s)
{
    pqn = s->GetPQN();
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

inline unsigned int OrbitalInfo::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

#endif
