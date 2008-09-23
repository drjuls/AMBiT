#ifndef VALENCE_CALCULATOR_H
#define VALENCE_CALCULATOR_H

#include "MBPTCalculator.h"

class ValenceCalculator : public MBPTCalculator
{
    /** Calculate valence-valence diagrams of many-body perturbation theory.
     */
public:
    ValenceCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    virtual ~ValenceCalculator(void) {}

    /** Returns one-electron valence-valence diagrams. */
    double GetOneElectronValence(const State* s1, const State* s2);

    /** Returns two-electron valence-valence diagrams. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronValence(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k);

    /** Returns two-electron "box" valence-valence diagrams, which can have wrong parity.
        Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronBoxValence(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k);
           
protected:
    /** Calculate one-electron valence-valence diagram of second order.
             x       x
             |       |
             |       |
        ->>------>------>>--
          i      2       f

         PRE: si.kappa == sf.kappa
     */
    double CalculateOneElectronValence1(const State& si, const State& sf) const;
    
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
    double CalculateTwoElectronValence1(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;
    double CalculateTwoElectronValence2(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;

protected:
    unsigned int MaxStateSize;
};

#endif
