#ifndef SINGLE_PARTICLE_INFO_H
#define SINGLE_PARTICLE_INFO_H

#include "Include.h"

class SingleParticleInfo
{
    /** Just stores pqn and kappa (L). Has an inbuilt ordering. */
public:
    SingleParticleInfo(unsigned int principal_qn, int kap);
    SingleParticleInfo(const SingleParticleInfo& other);
    virtual ~SingleParticleInfo() {}

    virtual bool operator<(const SingleParticleInfo& other) const;
    virtual bool operator==(const SingleParticleInfo& other) const;
    virtual bool operator!=(const SingleParticleInfo& other) const;

    inline unsigned int PQN() const { return pqn; }
    inline unsigned int L() const { return l; }
    inline double J() const;
    inline unsigned int TwoJ() const;

protected:
    unsigned int pqn;
    int kappa;
    unsigned int l;
};

inline SingleParticleInfo::SingleParticleInfo(unsigned int principal_qn, int kap):
    pqn(principal_qn), kappa(kap)
{
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline SingleParticleInfo::SingleParticleInfo(const SingleParticleInfo& other):
    pqn(other.pqn), kappa(other.kappa), l(other.l)
{}

inline double SingleParticleInfo::J() const
{
    return fabs(double(kappa)) - 0.5;
}

inline unsigned int SingleParticleInfo::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

inline bool SingleParticleInfo::operator<(const SingleParticleInfo& other) const
{
    if(this->pqn < other.pqn)
       return true;
    else if(this->pqn > other.pqn)
        return false;

    if(this->l < other.l)
        return true;
    else if(this->l > other.l)
        return false;
    else return (this->kappa > other.kappa);
}

inline bool SingleParticleInfo::operator==(const SingleParticleInfo& other) const
{
    if((this->pqn == other.pqn) &&
       (this->kappa == other.kappa))
        return true;
    else
        return false;
}

inline bool SingleParticleInfo::operator!=(const SingleParticleInfo& other) const
{
    return !(*this == other);
}

#endif