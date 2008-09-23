#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "Include.h"
#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "HartreeFock/SigmaPotential.h"
#include "SlaterIntegrals.h"

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

        Use Brillouin-Wigner perturbation theory, where the energy of external lines is
        kept constant in the energy denominator (this ensures that the operator is hermitian).
     */
public:
    MBPTCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states);
    virtual ~MBPTCalculator(void);

    unsigned int GetStorageSize(const ExcitedStates* valence_states);
    void UpdateIntegrals(const ExcitedStates* valence_states);

    /** Create a second-order one-electron MBPT (sigma1) operator. */
    void GetSecondOrderSigma(int kappa, SigmaPotential* sigma);

    /** Return value is the matrix element < s1 | Sigma1 | s2 >. */
    double GetOneElectronDiagrams(const StateInfo& s1, const StateInfo& s2) const;

    /** Returns subtraction diagrams for the matrix element < s1 | Sigma1 | s2 >. */
    double GetOneElectronSubtraction(const StateInfo& s1, const StateInfo& s2) const;

    /** Returns <s1, s2 | Sigma2(k) | s3, s4>. Only calculates in Brillouin-Wigner PT. */
    double GetTwoElectronDiagrams(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    /** Returns subtraction diagrams for the matrix element <s1, s2 | Sigma2(k) | s3, s4>. */
    double GetTwoElectronSubtraction(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    /** Returns "box" diagrams of <s1, s2 | Sigma2(k) | s3, s4> (numbers 4, 5, and 6).
        Only calculates in Brillouin-Wigner PT.
     */
    double GetTwoElectronBoxDiagrams(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    /** Add a constant, delta, to the energy denominator in all diagrams.
        This corresponds to H_0 -> H_0 + delta, so V -> V - delta
        and the energy denominator becomes
            1/(E + delta - H_0)
        To choose delta to adjust E closer to the BW value, use
            delta = E_CI - E_HF
     */
    void SetEnergyShift(double energy_shift)
    {   delta = energy_shift;
    }

protected:
    /** Calculate diagrams of second order (shown here 1 through 4).
        ->>------>------>>--  ->>------>------<---  ->>------<------>>--  ->>------<------>---
         a   |  beta |   b     a   |  beta |   n     a   |   m   |   b     a   |   m   | alpha
             |       |             |       |             |       |             |       |      
        --<------>------<---  --<------>------>>--  --<------>------<---  -->------<------>>--
         n     alpha     n     n     alpha     b     n     alpha     n    alpha    n       b  

         PRE: si.kappa == sf.kappa
     */
    void CalculateCorrelation1and3(int kappa, SigmaPotential* sigma) const;
    void CalculateCorrelation2(int kappa, SigmaPotential* sigma) const;
    void CalculateCorrelation4(int kappa, SigmaPotential* sigma) const;

    double CalculateCorrelation1and3(const StateInfo& sa, const StateInfo& sb) const;
    double CalculateCorrelation2(const StateInfo& sa, const StateInfo& sb) const;
    double CalculateCorrelation4(const StateInfo& sa, const StateInfo& sb) const;
    
    /** Calculate subtraction diagrams of second order (shown here 1 through 3).
        ->>------------->>--  ->>-------------<---  ->>------<------>>--
         a           |   b     a           |   n     a   |   n   |   b  
                     |                     |             |       |
        --<------>------<---  --<------>------>>--       x       x
         n   | alpha     n     n   | alpha     b
             |                     |
             x                     x
         PRE: sa.kappa == sb.kappa
         NOTE: diagrams 1 and 2 each represent two diagrams (swapping the interaction lines)
     */
    double CalculateSubtraction1(const StateInfo& sa, const StateInfo& sb) const;
    double CalculateSubtraction2(const StateInfo& sa, const StateInfo& sb) const;
    double CalculateSubtraction3(const StateInfo& sa, const StateInfo& sb) const;

    /** Calculate two-electron screening diagrams of second order (shown here 1 through 3).
        ->>------------->>--  ->>-------------<---  ->>------------->>--
          a          |   c      a          |   n      a          |   c  
                     |                     |                     |
        --<------>------<---  --<------>------>>--  ->>------>------<---
          n  | alpha     n      n  | alpha     c      b  | alpha     n     
             |                     |                     |
        ->>------------->>--  ->>------------->>--  --<-----------------
          b              d      b              d      n              d
         NOTE: each diagram has a mirror (swapping a<->b and c<->d)
     */
    double CalculateTwoElectron1(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;
    double CalculateTwoElectron2(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;
    double CalculateTwoElectron3(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;

    /** Calculate two-electron screening "box" diagrams of second order (shown here 4 through 6).
        ->>-------------<---  ->>------>------>>--  ->>-------------<---
          a          |   n      a  | alpha |   c      a          |   m  
                     |             |       |                     |
        --<------------->>--  ->>-------------<---  --<------------->>--
          n  |       |   c      b  |           n      m  |       |   c     
             |       |             |                     |       |
        ->>------>------>>--  --<------------->>--  ->>-------------<---
          b    alpha     d      n              d      b  |           n
                                                         |              
                                                    --<------------->>--
                                                      n              d  
         NOTE: diagrams 4 and 5 are mirrors of each other (swapping a<->b and c<->d)
     */
    double CalculateTwoElectron4(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;
    double CalculateTwoElectron5(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;
    double CalculateTwoElectron6(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;

    /** Calculate two-electron subtraction diagram of second order.
                     x
                     |          There are four diagrams, with the hole state
                     |          being connected to each of the four valence lines.
        ->>-------------<---
          a              n  
                            
        --<------------->>--
          n  |           c  
             |              
        ->>------------->>--
          b              d  
     */
    double CalculateTwoElectronSub(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const;

    /** Return absolute value of the difference of two unsigned integers. */
    unsigned int absdiff(unsigned int i, unsigned int j) const;

    /** Return minimum value of k connecting two states a and b, satisfying:
          triangle(a.J(), b.J(), kmin)
          eps(a.L() + b.L() + kmin)
    */
    unsigned int kmin(const StateInfo& a, const StateInfo& b) const;

    /** Return minimum value of k such that
          triangle(a.J(), b.J(), kmin)
          eps(a.L() + b.L() + kmin)
        and
          triangle(c.J(), d.J(), kmin)
          eps(c.L() + d.L() + kmin)
    */
    unsigned int kmin(const StateInfo& a, const StateInfo& b, const StateInfo& c, const StateInfo& d) const;

    /** Return maximum value of k that satisfies
          triangle(a.J(), b.J(), kmax)
    */
    unsigned int kmax(const StateInfo& a, const StateInfo& b) const;

    /** Return maximum value of k that satisfies
          triangle(a.J(), b.J(), kmax)
          triangle(c.J(), d.J(), kmax)
    */
    unsigned int kmax(const StateInfo& a, const StateInfo& b, const StateInfo& c, const StateInfo& d) const;
    
protected:
    Lattice* lattice;
    const Core* core;
    const ExcitedStates* excited;

    SlaterIntegrals* integrals;

    /** Valence energies for Brillouin-Wigner MBPT. */
    std::map<int, double> ValenceEnergies;
    void SetValenceEnergies();

    double delta;   // Shift in the energy denominator.
};

inline unsigned int MBPTCalculator::absdiff(unsigned int i, unsigned int j) const
{
    return abs(int(i) - int(j));
}

inline unsigned int MBPTCalculator::kmin(const StateInfo& a, const StateInfo& b) const
{
    // ensure eps(a.L() + b.L() + kmin)
    unsigned int ret = absdiff(a.L(), b.L());

    // ensure triangle(a.J(), b.J(), kmin)
    if(absdiff(a.TwoJ(), b.TwoJ()) > 2 * ret)
        ret += 2;

    return ret;
}

inline unsigned int MBPTCalculator::kmin(const StateInfo& a, const StateInfo& b, const StateInfo& c, const StateInfo& d) const
{
    unsigned int opt1 = kmin(a, b);
    unsigned int opt2 = kmin(c, d);

    return mmax(opt1, opt2);
}

inline unsigned int MBPTCalculator::kmax(const StateInfo& a, const StateInfo& b) const
{
    return (a.TwoJ() + b.TwoJ())/2;
}

inline unsigned int MBPTCalculator::kmax(const StateInfo& a, const StateInfo& b, const StateInfo& c, const StateInfo& d) const
{
    return mmin(a.TwoJ() + b.TwoJ(), c.TwoJ() + d.TwoJ())/2;
}

#endif
