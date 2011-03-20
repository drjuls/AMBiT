#ifndef CI_INTEGRALS_MBPT_H
#define CI_INTEGRALS_MBPT_H

#include "CIIntegrals.h"
#include "HartreeFock/SigmaPotential.h"
#include "MBPT/CoreMBPTCalculator.h"
#include "MBPT/ValenceCalculator.h"

class CIIntegralsMBPT: public CIIntegrals
{
    /** Can include MBPT effects in the integrals.
        Stores around twice as many two-electron integrals as CIIntegrals because the MBPT
        reduces their symmetry.
     */
public:
    CIIntegralsMBPT(const ExcitedStates& excited_states, const std::string& storage_id = ""):
        CIIntegrals(excited_states, storage_id), PT(NULL), ValencePT(NULL),
        include_sigma1(false), include_mbpt1(false), include_mbpt1_subtraction(false),
        include_mbpt2(false), include_mbpt2_subtraction(false), include_extra_box(false),
        include_valence_mbpt1(false), include_valence_mbpt2(false), include_valence_extra_box(false)
    {}
    virtual ~CIIntegralsMBPT();

    /** Calculate number of elements that will be stored. */
    virtual unsigned int GetStorageSize() const;

    /** Update all integrals (on the assumption that the excited states have changed).
        The sigma_id string can be used to get Sigma operators from disk, or else
        the normal storage read_id is used.
      */
    virtual void Update();
    virtual void Update(const std::string& sigma_id);

    /** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
    virtual double GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const;

    /** Include single-particle MBPT diagrams via sigma matrix.
        This option requires sigma matrices to exist, otherwise it will create them
        even if the integrals are already stored on disk.
        If just a few more integrals are required, consider using IncludeMBPT1.
     */
    inline void IncludeSigma1(bool include, CoreMBPTCalculator* mbpt = NULL)
    {   include_sigma1 = include;
        if(include)
            include_mbpt1 = false;
        if(states.GetCore()->IsOpenShellCore())
            include_mbpt1_subtraction = include;
        if(mbpt)
            PT = mbpt;
    }

    /** Include single-particle MBPT diagrams. */
    inline void IncludeMBPT1(bool include, CoreMBPTCalculator* mbpt = NULL)
    {   include_mbpt1 = include;
        if(include)
            include_sigma1 = false;
        if(states.GetCore()->IsOpenShellCore())
            include_mbpt1_subtraction = include;
        if(mbpt)
            PT = mbpt;
    }

    /** Include two-particle MBPT diagrams. */
    inline void IncludeMBPT2(bool include, CoreMBPTCalculator* mbpt = NULL)
    {   include_mbpt2 = include;
        if(states.GetCore()->IsOpenShellCore())
            include_mbpt2_subtraction = include;
        if(mbpt)
            PT = mbpt;
    }

    /** Include two-particle box diagrams. */
    inline void IncludeExtraBoxDiagrams(bool include)
    {   include_extra_box = include;
        if(include)
            SetExtraBoxDiagramLimits();
    }

    /** Include single-particle valence-valence MBPT diagrams. */
    inline void IncludeValenceMBPT1(bool include, ValenceCalculator* mbpt = NULL)
    {   include_valence_mbpt1 = include;
        if(mbpt)
            ValencePT = mbpt;
    }

    /** Include two-particle MBPT diagrams. */
    inline void IncludeValenceMBPT2(bool include, ValenceCalculator* mbpt = NULL)
    {   include_valence_mbpt2 = include;
        if(mbpt)
            ValencePT = mbpt;
    }

    /** Include two-particle box diagrams. */
    inline void IncludeValenceExtraBoxDiagrams(bool include)
    {   include_valence_extra_box = include;
        if(include)
            SetExtraBoxDiagramLimits();
    }

    /** Set limits on storage of extra box-diagram integrals (can check using GetStorageSize()).
        See CIIntegrals::SetTwoElectronStorageLimits() for details on structure.
     */
    void SetExtraBoxDiagramLimits(unsigned int limit1 = 0, unsigned int limit2 = 0, unsigned int limit3 = 0)
    {
        if(limit1) box_max_pqn_1 = limit1;
        else box_max_pqn_1 = 100;
            
        if(limit2) box_max_pqn_2 = limit2;
        else box_max_pqn_2 = 100;
        
        if(limit3) box_max_pqn_3 = limit3;
        else box_max_pqn_3 = 100;
    }

    /** Write out sigma potentials. */
    void WriteSigmaPotentials() const;

    /** Read multiple sets of one-electron integrals from binary *.one.int files.
        The number of files should be num_files, and they should be named
            name_0.two.int, name_1.two.int, ...
    */
    void ReadMultipleOneElectronIntegrals(const std::string& name, unsigned int num_files);

    /** Read multiple sets of two-electron integrals from binary *.two.int files.
        The number of files should be num_files, and they should be named
            name_0.two.int, name_1.two.int, ...
    */
    void ReadMultipleTwoElectronIntegrals(const std::string& name, unsigned int num_files);

    /** Temporary: Add SMS to stored integrals - this is just to upgrade old files. This function:
            UpdatesOneElectronIntegrals
            Reads in two electron integrals from "name"
            Adds SMS to stored integrals
            Writes two-electron integrals to "id"
    */
    void AddSMSToTwoElectronIntegrals(const std::string& name);

protected:
    virtual bool TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const;

    virtual void UpdateOneElectronIntegrals(const std::string& sigma_id);

    // Unlike CIIntegrals, the CIIntegralsMBPT class includes the SMS in the radial integral
    // directly. Therefore one must not update the two-electron integrals before the
    // one-electron integrals if the SMS != 0.
    virtual void UpdateTwoElectronIntegrals();

    /** Include a set of box diagrams of "wrong" parity with the two electron integrals.
        Should only be done after UpdateTwoElectronIntegals().
     */
    virtual void UpdateTwoElectronBoxDiagrams();

protected:
    CoreMBPTCalculator* PT;
    ValenceCalculator* ValencePT;

    /** Single-electron core MBPT effects. */
    bool include_sigma1;
    bool include_mbpt1;
    bool include_mbpt1_subtraction;
    std::map<int, SigmaPotential*> Sigma1;

    /** Two-electron core MBPT effects. */
    bool include_mbpt2;
    bool include_mbpt2_subtraction;
    bool include_extra_box;
    
    /** Valence-valence MBPT effects. */
    bool include_valence_mbpt1;
    bool include_valence_mbpt2;
    bool include_valence_extra_box;

    // Limits on extra two-body box diagram integrals.
    // Used for both core-valence and valence-valence MBPT diagrams.
    unsigned int box_max_pqn_1, box_max_pqn_2, box_max_pqn_3;
};

#endif
