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
    OrbitalInfo(unsigned int pqn, int kappa);
    OrbitalInfo(pOrbitalConst s);
    OrbitalInfo(const Orbital* s);
    OrbitalInfo(const OrbitalInfo& other);
    virtual ~OrbitalInfo(void) {}

    virtual bool operator<(const OrbitalInfo& other) const;
    virtual bool operator==(const OrbitalInfo& other) const;
    virtual bool operator!=(const OrbitalInfo& other) const;

    inline unsigned int PQN() const { return pqn; }
    inline int Kappa() const { return kappa; }
    inline unsigned int L() const { return l; }
    inline double J() const;
    inline unsigned int TwoJ() const;
    // Return the value of L for the lower component of the wavefunction
    inline unsigned int L_Prime() const;

    inline unsigned int MaxNumElectrons() const { return 2*abs(kappa); }
    virtual std::string Name() const;

protected:
    unsigned int pqn;
    int kappa;
    unsigned int l;
};

inline OrbitalInfo::OrbitalInfo(pOrbitalConst s)
{
    pqn = s->GetPQN();
    kappa = s->Kappa();
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline OrbitalInfo::OrbitalInfo(const Orbital* s)
{
    pqn = s->GetPQN();
    kappa = s->Kappa();
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline OrbitalInfo::OrbitalInfo(unsigned int principal_qn, int kap):
    pqn(principal_qn), kappa(kap)
{
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline OrbitalInfo::OrbitalInfo(const OrbitalInfo& other):
    pqn(other.pqn), kappa(other.kappa), l(other.l)
{}

inline double OrbitalInfo::J() const
{
    return fabs(double(kappa)) - 0.5;
}

// Returns the value of L for the lower component of the wavefunction
inline unsigned int OrbitalInfo::L_Prime() const
{
    // The lower component has L corresponding to -Kappa
    if(kappa < 0)
    {
        return (L() + 1);
    }
    else
    {
        return (L() - 1);
    }
}

inline unsigned int OrbitalInfo::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

#endif
