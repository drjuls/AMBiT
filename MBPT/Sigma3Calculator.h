#ifndef SIGMA3_CALCULATOR_H
#define SIGMA3_CALCULATOR_H

#include "MBPTCalculator.h"
#include "SlaterIntegrals.h"
#include "Configuration/ElectronInfo.h"

namespace Ambit
{
/** Calculate effective three-body diagrams of many-body perturbation theory.
    Three-body diagrams have the form:
    a___ _______ ___d
        |       |
    b___|       |___e
        |       |
    c___|_______|___f

    Sigma3Calculator is different to other MBPT calculators since these diagrams are calculated on the fly,
    rather than being pre-calculated and stored.
    Sigma3Calculator has the GetMatrixElement() method that can be used directly in ManyBodyOperator.
 */
class Sigma3Calculator : public MBPTCalculator
{
public:
    Sigma3Calculator(pOrbitalManagerConst orbitals, pSlaterIntegrals two_body, const std::string& fermi_orbitals = "");
    virtual ~Sigma3Calculator();

    void IncludeCore(bool include_mbpt);
    void IncludeValence(bool include_mbpt);

    virtual unsigned int GetStorageSize() override;
    virtual void UpdateIntegrals() override;

    /** Returns <e1, e2, e3 | Sigma2(k) | e4, e5, e6>.
        -->>------------->>--        -->>------------->>--
          a           |  d             a   |          d
                      |                    |
        -->>-------------<---    +   -->>------>------>>--
          b              n             b     alpha |  e
                                                   |
        --<-------------->>--        -->>------------->>--
          n   |          e             c              f
              |
        -->>------------->>--        (if include_valence)
          c              f
        and all permutations of lines.
      */
    double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
                            const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;

protected:
    /** Get the single diagram above (not the line permutations).
        Assumes momentum projections are okay:
            e1.TwoM() + e2.TwoM() + e3.TwoM() == e4.TwoM() + e5.TwoM() + e6.TwoM()
        and parity is okay:
            (e1.L() + e2.L() + e3.L() + e4.L() + e5.L() + e6.L())%2 == 0
     */
    double GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
                                const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;

    bool include_core;
    bool include_valence;
    pSlaterIntegrals two_body;

    pOrbitalMapConst deep;
    pOrbitalMapConst high;
};

typedef std::shared_ptr<Sigma3Calculator> pSigma3Calculator;

}
#endif
