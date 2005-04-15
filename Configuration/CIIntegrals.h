#ifndef CI_INTEGRALS_H
#define CI_INTEGRALS_H

#include "Include.h"
#include "Basis/ExcitedStates.h"
#include "ElectronInfo.h"

class CIIntegrals
{
    /** Class to hold Coulomb integrals for use in CI calculation. */
public:
    CIIntegrals(const ExcitedStates& excited_states, unsigned int limit1 = 100, unsigned int limit2 = 100, unsigned int limit3 = 100):
        states(excited_states), include_valence_sms(false)
    {   max_pqn_1 = limit1;
        max_pqn_2 = limit2;
        max_pqn_3 = limit3;
    }
    virtual ~CIIntegrals() {}

    /** Calculate number of elements that will be stored.
        48 bytes are stored for each element, so this should be kept to around 2 million,
        assuming it can take up 100MB.
        Use max_pqn_1, max_pqn_2 and max_pqn_3 to keep size down.
        For two electron integrals:
            i1.pqn <= max_pqn_1
            (i2.pqn or i3.pqn) <= max_pqn_2
            (i2.pqn and i3.pqn) <= max_pqn_3
        For 'x'* 3 waves (spd) and 'y'* 4 waves (spdf) in basis set,
            N = 61 x^4       N = 279 y^4.
        After max_pqn_1 = 4,
            N ~ 502 x^3      N ~ 1858 y^3,
        and hopefully after max_pqn_2 and then max_pqn_3
            N ~ x^2 and then N ~ x, respectively.
     */
    unsigned int GetStorageSize() const;

    /** Update all integrals (on the assumption that the excited states have changed). */
    void Update();

    /** Include the scaled specific mass shift in the two electron integrals. */
    inline void IncludeValenceSMS(bool include)
    {   include_valence_sms = include;
    }

    /* GetOneElectronIntegral(i, j) = <i|H|j> */
    double GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const;

    /* GetSMSIntegral(i, j) = <i|p|j> */
    double GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const;
    
    /* GetOverlapIntegral(i, j) = <i|j> */
    double GetOverlapIntegral(const StateInfo& s1, const StateInfo& s2) const;
    
    /* GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
    double GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    inline double GetNuclearInverseMass() const
    {   return states.GetCore()->GetNuclearInverseMass();
    }

protected:
    const ExcitedStates& states;

    unsigned int NumStates;
    std::map<StateInfo, unsigned int> state_index;
    std::map<unsigned int, StateInfo> reverse_state_index;

    // Storage for one and two electron integrals.
    // If these are null, it means that there is not enough space in memory to store them,
    // they must be generated as needed.

    // OneElectronIntegrals(i, j) = <i|H|j>
    std::map<unsigned int, double> OneElectronIntegrals;

    // TwoElectronIntegrals(k, i, j, l, m) = R_k(ij, lm): i->l, j->m
    std::map<unsigned int, double> TwoElectronIntegrals;

    // SMSIntegrals(i, j) = <i|p|j>
    std::map<unsigned int, double> SMSIntegrals;

    // Overlap = <i|j>
    std::map<unsigned int, double> OverlapIntegrals;

    bool include_valence_sms;

    unsigned int max_pqn_1, max_pqn_2, max_pqn_3;
};

#endif
