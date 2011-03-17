#ifndef VALENCE_CALCULATOR_H
#define VALENCE_CALCULATOR_H

#include "CoreMBPTCalculator.h"

class ValenceCalculator : public MBPTCalculator
{
    /** Calculate valence-valence diagrams of many-body perturbation theory.
     */
public:
    ValenceCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    virtual ~ValenceCalculator(void) {}

    /** No storage yet. */
    virtual unsigned int GetStorageSize(const ExcitedStates* valence_states)
    {   return 0;
    }
    virtual void UpdateIntegrals(const ExcitedStates* valence_states)
    {}

    /** Returns one-electron valence-valence diagrams. */
    double GetOneElectronValence(const SingleParticleWavefunction* s1, const SingleParticleWavefunction* s2);

    /** Returns two-electron valence-valence diagrams. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronValence(const SingleParticleWavefunction* s1, const SingleParticleWavefunction* s2, const SingleParticleWavefunction* s3, const SingleParticleWavefunction* s4, unsigned int k);

    /** Returns two-electron "box" valence-valence diagrams, which can have wrong parity.
        Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronBoxValence(const SingleParticleWavefunction* s1, const SingleParticleWavefunction* s2, const SingleParticleWavefunction* s3, const SingleParticleWavefunction* s4, unsigned int k);

protected:
    /** Calculate one-electron valence-valence diagram of second order.
             x       x
             |       |
             |       |
        ->>------>------>>--
          i      2       f

         PRE: si.kappa == sf.kappa
     */
    double CalculateOneElectronValence1(const SingleParticleWavefunction& si, const SingleParticleWavefunction& sf) const;
    
    /** Calculate two-electron valence-valence diagrams of second order.
                                           x
                                           |
                                           |
        ->>------>------>>--  ->>------>------>>--
          a  |   3   |   c      a  |   3       c
             |       |             |
        ->>------>------>>--  ->>------------->>--
          b      4       d      b              d

                              There are four diagrams, with the excited state
                              being connected to each of the four valence lines.
     */
    double CalculateTwoElectronValence1(const SingleParticleWavefunction& sa, const SingleParticleWavefunction& sb, const SingleParticleWavefunction& sc, const SingleParticleWavefunction& sd, unsigned int k) const;
    double CalculateTwoElectronValence2(const SingleParticleWavefunction& sa, const SingleParticleWavefunction& sb, const SingleParticleWavefunction& sc, const SingleParticleWavefunction& sd, unsigned int k) const;

protected:
    unsigned int MaxStateSize;
};

#endif
