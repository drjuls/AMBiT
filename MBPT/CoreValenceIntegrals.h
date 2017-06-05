#ifndef CORE_VALENCE_INTEGRALS_H
#define CORE_VALENCE_INTEGRALS_H

#include "SlaterIntegrals.h"
#include "CoreMBPTCalculator.h"
#include "ValenceMBPTCalculator.h"

/** Valence-valence SlaterIntegrals (for CI) including core and/or virtual correlations via MBPT.
    CoreValenceIntegrals<MapType> both is a SlaterIntegrals<MapType> and has one for use in the MBPTCalculator.
    Default setting: include all core MBPT and no valence MBPT. Change using IncludeCore()/IncludeValence().
    Writing to file is overriden, and occurs automatically when calculating as a failsafe.
 */
template <class MapType>
class CoreValenceIntegrals : public SlaterIntegrals<MapType>
{
protected:
    typedef typename SlaterIntegrals<MapType>::KeyType KeyType;

public:
    /** Use bare integrals with the same MapType: SlaterIntegrals<MapType>. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pHartreeY hartreeY_op, const std::string& write_file);

    /** Specify container object for bare integrals. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals bare_integrals, const std::string& write_file);

    /** Specify CoreMBPTCalculator and ValenceMBPTCalculator directly. */
    CoreValenceIntegrals(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator, pValenceMBPTCalculator valence_mbpt_calculator, const std::string& write_file);

    virtual ~CoreValenceIntegrals();

    virtual bool OffParityExists() const override { return include_core_extra_box || include_valence_extra_box; }

    /** Calculate two-electron requested MBPT. Write to file.
        PRE: OrbitalMaps should only include a subset of valence orbitals.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only = false) override;

#ifdef AMBIT_USE_MPI
    /** MPI version of Write() overridden to gather and write. */
    virtual void Write(const std::string& filename) const override;
#endif

    void IncludeCore(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams);
    void IncludeValence(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams);

protected:
    pCoreMBPTCalculator core_PT;
    pValenceMBPTCalculator valence_PT;

    /** Core MBPT effects. */
    bool include_core;
    bool include_core_subtraction;
    bool include_core_extra_box;

    /** Valence-valence MBPT effects. */
    bool include_valence;
    bool include_valence_subtraction;
    bool include_valence_extra_box;

    std::string write_file;
#ifdef AMBIT_USE_MPI
    bool my_calculations_done;
    mutable bool root_complete;
    std::vector<KeyType> new_keys;
    std::vector<double> new_values;
#endif
};

typedef CoreValenceIntegrals<std::map<unsigned long long int, double>> CoreValenceIntegralsMap;
typedef std::shared_ptr<CoreValenceIntegralsMap> pCoreValenceIntegralsMap;

#include "CoreValenceIntegrals.cpp"

#endif
