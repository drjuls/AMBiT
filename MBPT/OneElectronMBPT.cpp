#include "OneElectronMBPT.h"

OneElectronMBPT::OneElectronMBPT(pOrbitalManagerConst orbitals, pHFIntegrals bare_one_body, pSlaterIntegrals bare_two_body):
    OneElectronMBPT(orbitals, pCoreMBPTCalculator(new CoreMBPTCalculator(orbitals, bare_one_body, bare_two_body)), bare_one_body->GetOperator())
{}

OneElectronMBPT::OneElectronMBPT(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator, pSpinorMatrixElementConst pOperator):
    OneElectronIntegrals(orbitals, pOperator), core_PT(core_mbpt_calculator),
    include_core(true), include_core_subtraction(true), include_valence(false), include_valence_subtraction(false)
{}

unsigned int OneElectronMBPT::CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only)
{
    unsigned int i1, i2;

    std::set<unsigned int> found_keys;   // For check_size_only

    if(!check_size_only)
        core_PT->UpdateIntegrals();

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        const OrbitalInfo& s1 = it_1->first;

        auto it_2 = orbital_map_2->begin();
        while(it_2 != orbital_map_2->end())
        {
            if(it_1->first.Kappa() == it_2->first.Kappa())
            {
                i2 = orbitals->state_index.at(it_2->first);
                const OrbitalInfo& s2 = it_2->first;

                if(check_size_only)
                {
                    if(include_core || include_core_subtraction || include_valence || include_valence_subtraction)
                        found_keys.insert(GetKey(i1, i2));
                }
                else
                {   double value = 0.;
                    if(include_core)
                        value += core_PT->GetOneElectronDiagrams(s1, s2);
                    if(include_core_subtraction)
                        value += core_PT->GetOneElectronSubtraction(s1, s2);

                    integrals[GetKey(i1, i2)] = value;
                }
            }
            it_2++;
        }
        it_1++;
    }

    if(check_size_only)
        return found_keys.size();
    else
        return integrals.size();
}

void OneElectronMBPT::IncludeCore(bool include_mbpt, bool include_subtraction)
{
    include_core = include_mbpt;
    include_core_subtraction = include_subtraction;
}

void OneElectronMBPT::IncludeValence(bool include_mbpt, bool include_subtraction)
{
    include_valence = include_mbpt;
    include_valence_subtraction = include_subtraction;
}
