#ifndef SIGMA3_CALCULATOR_H
#define SIGMA3_CALCULATOR_H

#include "MBPTCalculator.h"
#include "Configuration/ElectronInfo.h"

class Sigma3Calculator : public MBPTCalculator
{
    /** Calculate effective three-body diagrams of many-body perturbation theory.
        Three-body diagrams have the form:
          a___ _______ ___d
              |       |
          b___|       |___e
              |       |
          c___|_______|___f
     */
public:
    Sigma3Calculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    virtual ~Sigma3Calculator(void) {}

    /** Returns <e1, e2, e3 | Sigma2(k) | e4, e5, e6>. Only calculates in Brillouin-Wigner PT.
        -->>------------->>--
          a           |  d
                      |
        -->>-------------<---
          b              n

        --<-------------->>--
          n   |          e
              |
        -->>------------->>--
          c              f

        Assumes momentum projections are okay:
            e1.TwoM() + e2.TwoM() + e3.TwoM() == e4.TwoM() + e5.TwoM() + e6.TwoM()
        and parity is okay:
            (e1.L() + e2.L() + e3.L() + e4.L() + e5.L() + e6.L())%2 == 0
      */
    double GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6);

protected:
    /** Update all integrals (on the assumption that the excited states have changed). */
    virtual void Update();

    virtual void UpdateStateIndexes();

    /** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
    virtual double GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4);

    /** GetSMSIntegral(i, j) = <i|p|j> */
    double GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const;

    /** Change ordering of states so that it corresponds to a stored integral.
        Returns false if SMS sign needs to be changed.
     */
    virtual bool TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const;

protected:
    unsigned int NumStates;
    // The ordering of states is not arbitrary; they should be ordered by pqn first.
    std::map<StateInfo, unsigned int> state_index;
    std::map<unsigned int, StateInfo> reverse_state_index;

    // SMSIntegrals(i, j) = <i|p|j>
    std::map<unsigned int, double> SMSIntegrals;

    // TwoElectronIntegrals(k, i, j, l, n) = R_k(ij, ln): i->l, j->n
    std::map<unsigned int, double> TwoElectronIntegrals;
};

#endif
