#ifndef ONE_ELECTRON_MBPT_H
#define ONE_ELECTRON_MBPT_H

#include "OneElectronIntegrals.h"
#include "CoreMBPTCalculator.h"

class OneElectronMBPT : public HFIntegrals
{
public:
    OneElectronMBPT(pOrbitalManagerConst orbitals, pHFIntegrals bare_one_body, pSlaterIntegrals bare_two_body);
    OneElectronMBPT(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator, pSpinorMatrixElementConst pOperator = nullptr);

    virtual ~OneElectronMBPT() {}

    /** Calculate one-electron integrals and return number of integrals stored.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
        PRE: OrbitalMaps should only include a subset of valence orbitals.
     */
    virtual unsigned int CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only = false);

    void IncludeCore(bool include_mbpt, bool include_subtraction);
    void IncludeValence(bool include_mbpt, bool include_subtraction);

protected:
    pCoreMBPTCalculator core_PT;
    //ValenceCalculator* valence_PT;

    /** Core MBPT effects. */
    bool include_core;
    bool include_core_subtraction;

    /** Valence-valence MBPT effects. */
    bool include_valence;
    bool include_valence_subtraction;
};

typedef std::shared_ptr<OneElectronMBPT> pOneElectronMBPT;

#endif
