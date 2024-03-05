#ifndef ONE_ELECTRON_MBPT_H
#define ONE_ELECTRON_MBPT_H

#include "OneElectronIntegrals.h"
#include "CoreMBPTCalculator.h"
#include "ValenceMBPTCalculator.h"

namespace Ambit
{
/** Valence-valence HFIntegrals including core and/or virtual correlations via MBPT.
    Default setting: include all core MBPT and no valence MBPT. Change using IncludeCore()/IncludeValence().
    Writing to file is overriden, and occurs automatically when calculating as a failsafe.
 */
class OneElectronMBPT : public HFIntegrals
{
public:
    OneElectronMBPT(pOrbitalManagerConst orbitals, pHFIntegrals bare_one_body, pSlaterIntegrals bare_two_body, const std::string& write_file);
    OneElectronMBPT(pOrbitalManagerConst orbitals, pSpinorMatrixElementConst pOperator, pCoreMBPTCalculator core_mbpt_calculator, pValenceMBPTCalculator valence_mbpt_calculator, const std::string& write_file);

    virtual ~OneElectronMBPT() {}

    /** Calculate one-electron integrals and return number of integrals stored.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
        PRE: OrbitalMaps should only include a subset of valence orbitals.
     */
    virtual unsigned int CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only = false);

#ifdef AMBIT_USE_MPI
    /** MPI version of Write() overridden to gather and write. */
    virtual void Write(const std::string& filename) const override;
#endif

    void IncludeCore(bool include_mbpt, bool include_subtraction);
    void IncludeValence(bool include_subtraction);

protected:
    pCoreMBPTCalculator core_PT;
    pValenceMBPTCalculator valence_PT;

    /** Core MBPT effects. */
    bool include_core;
    bool include_core_subtraction;

    /** Valence-valence MBPT effects. */
    bool include_valence_subtraction;

    std::vector<unsigned int> new_keys;
    std::vector<double> new_values;
    std::string write_file;
};

typedef std::shared_ptr<OneElectronMBPT> pOneElectronMBPT;

}
#endif
