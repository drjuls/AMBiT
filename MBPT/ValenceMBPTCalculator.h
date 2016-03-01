#ifndef VALENCE_CALCULATOR_H
#define VALENCE_CALCULATOR_H

#include "CoreMBPTCalculator.h"

/** Calculate valence-valence diagrams of many-body perturbation theory.
    In these diagrams there are no (non-valence) holes.
 */
class ValenceMBPTCalculator : public MBPTCalculator
{
public:
    ValenceMBPTCalculator(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body);
    ValenceMBPTCalculator(pOrbitalManagerConst orbitals, pHFIntegrals one_body);
    virtual ~ValenceMBPTCalculator();

    /** No storage yet. */
    virtual unsigned int GetStorageSize() override;
    virtual void UpdateIntegrals() override;

    /** Returns one-electron valence-valence (subtraction) diagrams. */
    double GetOneElectronSubtraction(const OrbitalInfo& s1, const OrbitalInfo& s2);

    /** Returns two-electron valence-valence diagrams. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronValence(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4);
    double GetTwoElectronSubtraction(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4);

    /** Returns two-electron "box" valence-valence diagrams, which can have wrong parity.
        Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronBoxValence(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4);

protected:
    /** Calculate one-electron valence-valence diagram of second order.
             x       x
             |       |
             |       |
        ->>------>------>>--
          a    alpha     b

         PRE: si.kappa == sf.kappa
     */
    double CalculateOneElectronSub(const OrbitalInfo& sa, const OrbitalInfo& sb) const;
    
    /** Calculate two-electron valence-valence diagrams of second order.
                                           x
                                           |
                                           |
        ->>------>------>>--  ->>------>------>>--
          a  | alpha |   c      a  | alpha     c
             |       |             |
        ->>------>------>>--  ->>------------->>--
          b    beta      d      b              d

                              There are four diagrams, with the excited state
                              being connected to each of the four valence lines.
     */
    double CalculateTwoElectronValence(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const;
    double CalculateTwoElectronSub(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const;

protected:
    pHFIntegrals one_body;
    pSlaterIntegrals two_body;

    pOrbitalMapConst excited;
    pOrbitalMapConst high;
};

typedef std::shared_ptr<ValenceMBPTCalculator> pValenceMBPTCalculator;

#endif
