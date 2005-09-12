#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "HartreeFock/SigmaPotential.h"

class MBPTCalculator
{
    /** Calculate diagrams of many-body perturbation theory.
        One-body diagrams (functions GetSecondOrderSigma()) are of the form
               _____
         s1___|     |___s2    (non-zero only if s1.kappa == s2.kappa)
              |_____|

        Two-body diagrams (function GetTwoElectronDiagrams()) are of the form
          a___ _____ ___b
              |     |
          c___|_____|___d
     */
public:
    MBPTCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    ~MBPTCalculator(void) {}

    /** Create a sigma operator for the given state to second order.
     */
    void GetSecondOrderSigma(int kappa, SigmaPotential* sigma);

    /** Create a sigma operator for the given state to second order.
        if(sigma == NULL)
            no sigma potential is created, the function just calculates the matrix element.
        Return value is the energy of the state including second order matrix element.
     */
    double GetSecondOrderSigma(const State* s, SigmaPotential* sigma = NULL);

    /** Return value is the matrix element <s1 | Sigma1 | s2>. */
    double GetSecondOrderSigma(const State* s1, const State* s2, SigmaPotential* sigma = NULL);

    /** Returns subtraction diagrams for the matrix element <s1 | Sigma1 | s2>. */
    double GetSigmaSubtraction(const State* s1, const State* s2);

    /** Returns <s1, s2 | Sigma2(k) | s3, s4>. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronDiagrams(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k);

    /** Returns subtraction diagrams for the matrix element <s1, s2 | Sigma2(k) | s3, s4>. */
    double GetTwoElectronSubtraction(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k);

    /** Returns "box" diagrams of <s1, s2 | Sigma2(k) | s3, s4> (numbers 4, 5, and 6).
        Only calculates in Brillouin-Wigner PT.
     */
    double GetTwoElectronBoxDiagrams(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k);

    /** Use Brillouin-Wigner perturbation theory, where the energy of external lines is
        kept constant in the energy denominator (this ensures that the operator is hermitian).
     */
    void UseBrillouinWignerPT()
    {   SetValenceEnergies();
        BrillouinWignerPT = true;
    }

    /** Use Rayleigh-Shrodinger perturbation theory. */
    void UseRayleighShrodingerPT()
    {   BrillouinWignerPT = false;
    }

protected:
    /** Calculate diagrams of second order (shown here 1 through 4).
        ->>------>------>>--  ->>------>------<---  ->>------<------>>--  ->>------<------>---
         i   |   3   |   f     i   |   3   |   2     i   |   3   |   f     i   |   3   |   4  
             |       |             |       |             |       |             |       |      
        --<------>------<---  --<------>------>>--  --<------>------<---  -->------<------>>--
         2       4       2     2       4       f     2       4       2     4       2       f  

         PRE: si.kappa == sf.kappa
     */
    double CalculateCorrelation1and3(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateCorrelation2(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateCorrelation4(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    
    /** Calculate subtraction diagrams of second order (shown here 1 through 3).
        ->>------------->>--  ->>------>------<---  ->>------<------>>--
         i           |   f     i       3   |   2     i   |   2   |   f  
                     |                     |             |       |
        --<------>------<---  --<------>------>>--       x       x
         2   |   4       2     2   |   4       f
             |                     |
             x                     x
         PRE: si.kappa == sf.kappa
         NOTE: diagrams 1 and 2 each represent two diagrams (swapping the interaction lines)
     */
    double CalculateSubtraction1(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateSubtraction2(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;
    double CalculateSubtraction3(const State& si, const State& sf, SigmaPotential* sigma = NULL) const;

    /** Calculate two-electron screening diagrams of second order (shown here 1 through 3).
        ->>------------->>--  ->>-------------<---  ->>------------->>--
          a          |   c      a          |   2      a          |   c  
                     |                     |                     |
        --<------>------<---  --<------>------>>--  ->>------>------<---
          2  |   4       2      2  |   4       c      b  |   4       2     
             |                     |                     |
        ->>------------->>--  ->>------------->>--  --<-----------------
          b              d      b              d      2              d
         NOTE: each diagram has a mirror (swapping a<->b and c<->d)
     */
    double CalculateTwoElectron1(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;
    double CalculateTwoElectron2(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;
    double CalculateTwoElectron3(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;

    /** Calculate two-electron screening "box" diagrams of second order (shown here 4 through 6).
        ->>-------------<---  ->>------>------>>--  ->>-------------<---
          a          |   2      a  |   4   |   c      a          |   1  
                     |             |       |                     |
        --<------------->>--  ->>-------------<---  --<------------->>--
          2  |       |   c      b  |           2      1  |       |   c     
             |       |             |                     |       |
        ->>------>------>>--  --<------------->>--  ->>-------------<---
          b      4       d      2              d      b  |           2
                                                         |              
                                                    --<------------->>--
                                                      2              d  
         NOTE: diagrams 4 and 5 are mirrors of each other (swapping a<->b and c<->d)
     */
    double CalculateTwoElectron4(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;
    double CalculateTwoElectron5(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;
    double CalculateTwoElectron6(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;

    /** Calculate two-electron subtraction diagram of second order.
                     x
                     |          There are four diagrams, with the hole state
                     |          being connected to each of the four valence lines.
        ->>-------------<---
          a              2  
                            
        --<------------->>--
          2  |           c  
             |              
        ->>------------->>--
          b              d  
     */
    double CalculateTwoElectronSub(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const;

protected:
    Lattice* lattice;
    const Core* core;
    const ExcitedStates* excited;

    unsigned int MaxStateSize;

    bool BrillouinWignerPT;

    /** Valence energies for Brillouin-Wigner MBPT. */
    std::map<int, double> ValenceEnergies;
    void SetValenceEnergies();
};

#endif
