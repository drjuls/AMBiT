#ifndef CORE_VALENCE_INTEGRALS_H
#define CORE_VALENCE_INTEGRALS_H

#include "SlaterIntegrals.h"
#include "CoreMBPTCalculator.h"

/** Valence-valence SlaterIntegrals (for CI) including core and/or virtual correlations via MBPT.
    CoreValenceIntegrals<MapType> both is a SlaterIntegrals<MapType> and has one for use in the MBPTCalculator.
    Default setting: include all core MBPT and no valence MBPT. Change using IncludeCore()/IncludeValence().
 */
template <class MapType>
class CoreValenceIntegrals : public SlaterIntegrals<MapType>
{
protected:
    typedef typename SlaterIntegrals<MapType>::KeyType KeyType;

public:
    /** Use bare integrals with the same MapType: SlaterIntegrals<MapType>. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pHartreeY hartreeY_op);

    /** Specify container object for bare integrals. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals bare_integrals);

    /** Specify CoreMBPTCalculator directly. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator);

    virtual ~CoreValenceIntegrals();

    /** Calculate two-electron Slater integrals including requested MBPT.
        PRE: OrbitalMaps should only include a subset of valence orbitals.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only = false) override;

    void IncludeCore(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams);
    void IncludeValence(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams);

protected:
    pCoreMBPTCalculator core_PT;
    //ValenceCalculator* valence_PT;

    /** Core MBPT effects. */
    bool include_core;
    bool include_core_subtraction;
    bool include_core_extra_box;

    /** Valence-valence MBPT effects. */
    bool include_valence;
    bool include_valence_subtraction;
    bool include_valence_extra_box;
};

typedef CoreValenceIntegrals<std::map<unsigned long long int, double>> CoreValenceIntegralsMap;
typedef std::shared_ptr<CoreValenceIntegralsMap> pCoreValenceIntegralsMap;

#include "CoreValenceIntegrals.cpp"

#endif
