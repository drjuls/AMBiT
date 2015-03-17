#ifndef VALENCE_CALCULATOR_H
#define VALENCE_CALCULATOR_H

#include "CoreMBPTCalculator.h"

class ValenceCalculator : public MBPTCalculator
{
    /** Calculate valence-valence diagrams of many-body perturbation theory.
     */
public:
    ValenceCalculator(pLattice lattice, pCoreConst atom_core, pExcitedStatesConst excited_states);
    virtual ~ValenceCalculator(void) {}

    /** No storage yet. */
    virtual unsigned int GetStorageSize(pExcitedStatesConst valence_states)
    {   return 0;
    }
    virtual void UpdateIntegrals(pExcitedStatesConst valence_states)
    {}

    /** Returns one-electron valence-valence diagrams. */
    double GetOneElectronValence(pOrbitalConst s1, pOrbitalConst s2);

    /** Returns two-electron valence-valence diagrams. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronValence(pOrbitalConst s1, pOrbitalConst s2, pOrbitalConst s3, pOrbitalConst s4, unsigned int k);

    /** Returns two-electron "box" valence-valence diagrams, which can have wrong parity.
        Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronBoxValence(pOrbitalConst s1, pOrbitalConst s2, pOrbitalConst s3, pOrbitalConst s4, unsigned int k);

protected:
    /** Calculate one-electron valence-valence diagram of second order.
             x       x
             |       |
             |       |
        ->>------>------>>--
          i      2       f

         PRE: si.kappa == sf.kappa
     */
    double CalculateOneElectronValence1(const Orbital& si, const Orbital& sf) const;
    
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
    double CalculateTwoElectronValence1(const Orbital& sa, const Orbital& sb, const Orbital& sc, const Orbital& sd, unsigned int k) const;
    double CalculateTwoElectronValence2(const Orbital& sa, const Orbital& sb, const Orbital& sc, const Orbital& sd, unsigned int k) const;

protected:
    unsigned int MaxStateSize;
};

#endif
