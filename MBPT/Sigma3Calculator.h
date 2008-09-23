#ifndef SIGMA3_CALCULATOR_H
#define SIGMA3_CALCULATOR_H

#include "CoreMBPTCalculator.h"
#include "SlaterIntegrals.h"
#include "Configuration/ElectronInfo.h"

class Sigma3Integrals : public SlaterIntegrals
{
    /** Needs to store two electron integrals of the form
          < n a | 1/r | b c >
        that is, three valence and one core line.
        One-electron integrals are not needed, but SMS integrals are.
     */
public:
    Sigma3Integrals(const ExcitedStates* excited_states):
        SlaterIntegrals(excited_states)
    {}
    virtual ~Sigma3Integrals() {}

    /** Return number of two-electron integrals that will be stored. */
    virtual unsigned int GetStorageSize(const ExcitedStates& valence);

protected:
    /** No one-electron Hamiltonian integrals, but SMS integrals are still needed. */
    virtual void UpdateOneElectronIntegrals();
    virtual void UpdateTwoElectronIntegrals();
};

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
    virtual ~Sigma3Calculator(void);

    virtual unsigned int GetStorageSize(const ExcitedStates* valence_states);
    virtual void UpdateIntegrals(const ExcitedStates* valence_states);

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
    Sigma3Integrals* integrals;
};

#endif
