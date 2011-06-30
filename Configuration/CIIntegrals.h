#ifndef CI_INTEGRALS_H
#define CI_INTEGRALS_H

#include "Include.h"
#include "Basis/ExcitedStates.h"

class CIIntegrals
{
    /** Class to hold Coulomb integrals for use in CI calculation. */
public:
    /** If storage_id is used we assume that the user either wants to use stored integrals
        or output integrals once calculated.
     */
    CIIntegrals(const ExcitedStates& excited_states, const std::string& storage_id = ""):
        states(excited_states), include_valence_sms(false)
    {   SetTwoElectronStorageLimits();
        SetIdentifier(storage_id);
        UpdateStateIndexes();
    }
    virtual ~CIIntegrals() {}

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
    void SetTwoElectronStorageLimits(unsigned int limit1 = 0, unsigned int limit2 = 0, unsigned int limit3 = 0)
    {
        if(limit1) max_pqn_1 = limit1;
        else max_pqn_1 = 100;
            
        if(limit2) max_pqn_2 = limit2;
        else max_pqn_2 = 100;
        
        if(limit3) max_pqn_3 = limit3;
        else max_pqn_3 = 100;
    }

    /** Calculate number of one-electron and two-electron integrals that will be stored.
        Return total.
     */
    virtual unsigned int GetStorageSize() const;

    /** Clear all integrals. */
    virtual void Clear();

    /** Update all integrals (on the assumption that the excited states have changed). */
    virtual void Update();

    /** Include the scaled specific mass shift in the two electron integrals. */
    inline void IncludeValenceSMS(bool include)
    {   include_valence_sms = include;
    }

    /** GetOneElectronIntegral(i, j) = <i|H|j> */
    double GetOneElectronIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const;

    /** GetSMSIntegral(i, j) = <i|p|j> */
    double GetSMSIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const;
    
    /** GetOverlapIntegral(i, j) = <i|j>.
        PRE: i.L() == j.L()
     */
    double GetOverlapIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const;
    
    /** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
    virtual double GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const;

    inline double GetNuclearInverseMass() const
    {   return states.GetCore()->GetNuclearInverseMass();
    }

    /** The identifier is used to choose filenames for integrals. */
    virtual void SetIdentifier(const std::string& storage_id = "");

    /** Structure of *.one.int and *.two.int files:
         stores basis information (start and end pqn for each wave)
                        s-wave          p-wave        d-wave
           -------------------------------------------------------------
           | MAX_L | start |  end  | start |  end  | start |  end  | ...
           |       |  pqn  |  pqn  |  pqn  |  pqn  |  pqn  |  pqn  |
           -------------------------------------------------------------
         then integrals
           -------------------------------------------------------------
           | size  | index |    value      | index |    value      | ...
           |       |       |   (double)    |       |   (double)    |
           -------------------------------------------------------------
     */

    /** Write single electron integrals to binary *.one.int file. */
    void WriteOneElectronIntegrals(bool use_read_id = false) const;

    /** Write two-electron integrals to binary *.two.int file. */
    void WriteTwoElectronIntegrals(bool use_read_id = false) const;

protected:
    /** Change ordering of states so that it corresponds to a stored integral.
        Returns false if SMS sign needs to be changed.
     */
    virtual bool TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const;

    /** Read single electron integrals from binary *.one.int file. */
    void ReadOneElectronIntegrals(FILE* fp);

    /** Read two-electron integrals from binary *.two.int file. */
    void ReadTwoElectronIntegrals(FILE* fp);

    virtual void UpdateStateIndexes();
    virtual void UpdateOneElectronIntegrals();
    virtual void UpdateTwoElectronIntegrals();

protected:
    std::string read_id;    // Read and write files may be different (e.g. multiprocessor)
    std::string write_id;
    const ExcitedStates& states;

    unsigned int NumStates;
    // The ordering of states is not arbitrary; they should be ordered by pqn first.
    std::map<OrbitalInfo, unsigned int> state_index;
    std::map<unsigned int, OrbitalInfo> reverse_state_index;

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

    // Include SMS in two-body integrals.
    bool include_valence_sms;

    // Limits on stored two-body integrals.
    unsigned int max_pqn_1, max_pqn_2, max_pqn_3;
};

#endif
