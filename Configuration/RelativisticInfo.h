#ifndef RELATIVISTIC_INFO_H
#define RELATIVISTIC_INFO_H

#include "SingleParticleInfo.h"
#include "HartreeFock/StateInfo.h"
#include <string>

class NonRelInfo;

class RelativisticInfo : public SingleParticleInfo
{
public:
    RelativisticInfo(unsigned int principal_qn, int kappa):
        SingleParticleInfo(principal_qn, kappa)
    {}
    RelativisticInfo(const SingleParticleInfo& other):
        SingleParticleInfo(other)
    {}
    virtual ~RelativisticInfo(void) {}

    inline int Kappa() const { return kappa; }

    NonRelInfo GetNonRelInfo() const;
    StateInfo GetStateInfo() const;

    inline unsigned int MaxNumElectrons() const { return 2*abs(kappa); }
    virtual std::string Name() const;
};

#endif
