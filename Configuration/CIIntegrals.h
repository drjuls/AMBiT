#ifndef CI_INTEGRALS_H
#define CI_INTEGRALS_H

#include "Include.h"
#include "Basis/ExcitedStates.h"
#include "ElectronInfo.h"
#include "HartreeFock/SigmaPotential.h"
#include "MBPT/MBPTCalculator.h"

class CIIntegrals
{
    /** Class to hold Coulomb integrals for use in CI calculation. */
public:
    /** If storage_id is used we assume that the user either wants to use stored integrals
        or output integrals once calculated.
     */
    CIIntegrals(const ExcitedStates& excited_states, const std::string& storage_id = ""):
        states(excited_states), id(storage_id), PT(NULL),
        include_valence_sms(false), include_sigma1(false)
    {   SetTwoElectronStorageLimits();
        SetValenceEnergies();
    }
    virtual ~CIIntegrals();

    /** Set limits on storage of two electron integrals (can check using GetStorageSize()).
        48 bytes are stored for each element, so this should be kept to around 2 million,
        assuming it can take up 100MB.
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
    void SetTwoElectronStorageLimits(unsigned int limit1 = 100, unsigned int limit2 = 100, unsigned int limit3 = 100)
    {   max_pqn_1 = limit1;
        max_pqn_2 = limit2;
        max_pqn_3 = limit3;
    }

    /** Calculate number of elements that will be stored. */
    unsigned int GetStorageSize() const;

    /** Update all integrals (on the assumption that the excited states have changed).
        The sigma_id string can be used to get Sigma operators from disk, or else
        the normal storage id is used.
      */
    void Update(const std::string& sigma_id = "");

    /** Include the scaled specific mass shift in the two electron integrals. */
    inline void IncludeValenceSMS(bool include)
    {   include_valence_sms = include;
    }

    /** Include the scaled specific mass shift in the two electron integrals. */
    inline void IncludeSigma1(bool include, MBPTCalculator* mbpt = NULL)
    {   include_sigma1 = include;
        if(mbpt)
            PT = mbpt;
    }

    /** GetOneElectronIntegral(i, j) = <i|H|j> */
    double GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const;

    /** GetSMSIntegral(i, j) = <i|p|j> */
    double GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const;
    
    /** GetOverlapIntegral(i, j) = <i|j>.
        PRE: i.L() == j.L()
     */
    double GetOverlapIntegral(const StateInfo& s1, const StateInfo& s2) const;
    
    /** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
    double GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    inline double GetNuclearInverseMass() const
    {   return states.GetCore()->GetNuclearInverseMass();
    }
 
    inline void SetIdentifier(const std::string& storage_id = "")
    {   id = storage_id;
    }

    /** Write out sigma potentials. */
    void WriteSigmaPotentials() const;

    /** Write single electron integrals to binary *.one.int file. */
    void WriteOneElectronIntegrals() const;

protected:
    /** Read single electron integrals from binary *.one.int file. */
    void ReadOneElectronIntegrals(FILE* fp);

protected:
    std::string id;
    const ExcitedStates& states;
    MBPTCalculator* PT;

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

    /** Include SMS in two-body integrals. */
    bool include_valence_sms;

    /** Limits on stored two-body integrals. */
    unsigned int max_pqn_1, max_pqn_2, max_pqn_3;

    /** Single-electron MBPT effects. */
    bool include_sigma1;
    std::map<int, SigmaPotential*> Sigma1;

    /** Valence energies for single-electron Brillouin-Wigner MBPT. */
    std::map<int, double> ValenceEnergies;
    void SetValenceEnergies();
};

#endif
