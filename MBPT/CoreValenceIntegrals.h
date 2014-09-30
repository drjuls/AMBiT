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
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pOneElectronIntegrals one_body, pHartreeY hartreeY_op);

    /** Specify container object for bare integrals. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pOneElectronIntegrals one_body, pSlaterIntegrals bare_integrals);

    /** Specify CoreMBPTCalculator directly. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator);

    virtual ~CoreValenceIntegrals();

    /** Calculate two-electron Slater integrals including requested MBPT.
        PRE: OrbitalMaps should only include a subset of valence orbitals.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only = false) override;

    /** Calculate two electron integrals on valence orbitals with limits on the PQNs of the orbitals.
        Use max_pqn_1, max_pqn_2 and max_pqn_3 to keep size down.
        For two electron integrals:
            i1.pqn <= limit1
            (i2.pqn or i3.pqn) <= limit2
            (i2.pqn and i3.pqn) <= limit3
        For 'x'* 3 waves (spd) and 'y'* 4 waves (spdf) in basis set,
            N = 61 x^4       N = 279 y^4.
        After max_pqn_1 = 4,
            N ~ 502 x^3      N ~ 1858 y^3,
        and hopefully after max_pqn_2 and then max_pqn_3
            N ~ x^2 and then N ~ x, respectively.
     */
    //virtual unsigned int CalculateTwoElectronIntegrals(int limit1 = 0, int limit2 = 0, int limit3 = 0, bool check_size_only = false);

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
