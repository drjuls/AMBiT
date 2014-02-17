#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "Include.h"
#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"

class MBPTCalculator
{
    /** Calculate diagrams of many-body perturbation theory.
        Use Brillouin-Wigner perturbation theory, where the energy of external lines is
        kept constant in the energy denominator (this ensures that the operator is hermitian).
     */
public:
    MBPTCalculator(pLattice lattice, pCoreConst atom_core, pExcitedStatesConst excited_states);
    virtual ~MBPTCalculator(void);

    virtual unsigned int GetStorageSize(pExcitedStatesConst valence_states) = 0;
    virtual void UpdateIntegrals(pExcitedStatesConst valence_states) = 0;

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
    /** Return absolute value of the difference of two unsigned integers. */
    inline int absdiff(int i, int j) const;

    /** Return minimum value of k connecting two states a and b, satisfying:
          triangle(a.J(), b.J(), kmin)
          eps(a.L() + b.L() + kmin)
    */
    inline int kmin(const OrbitalInfo& a, const OrbitalInfo& b) const;

    /** Return minimum value of k such that
          triangle(a.J(), b.J(), kmin)
          eps(a.L() + b.L() + kmin)
        and
          triangle(c.J(), d.J(), kmin)
          eps(c.L() + d.L() + kmin)
    */
    inline int kmin(const OrbitalInfo& a, const OrbitalInfo& b, const OrbitalInfo& c, const OrbitalInfo& d) const;

    /** Return maximum value of k that satisfies
          triangle(a.J(), b.J(), kmax)
    */
    inline int kmax(const OrbitalInfo& a, const OrbitalInfo& b) const;

    /** Return maximum value of k that satisfies
          triangle(a.J(), b.J(), kmax)
          triangle(c.J(), d.J(), kmax)
    */
    inline int kmax(const OrbitalInfo& a, const OrbitalInfo& b, const OrbitalInfo& c, const OrbitalInfo& d) const;
    
protected:
    pLattice lattice;
    pCoreConst core;
    pExcitedStatesConst excited;

    /** Valence energies for Brillouin-Wigner MBPT. */
    std::map<int, double> ValenceEnergies;
    void SetValenceEnergies();

    double delta;   // Shift in the energy denominator.
};

inline int MBPTCalculator::absdiff(int i, int j) const
{
    return abs(i - j);
}

inline int MBPTCalculator::kmin(const OrbitalInfo& a, const OrbitalInfo& b) const
{
    // ensure eps(a.L() + b.L() + kmin)
    int ret = absdiff(a.L(), b.L());

    // ensure triangle(a.J(), b.J(), kmin)
    if(absdiff(a.TwoJ(), b.TwoJ()) > 2 * ret)
        ret += 2;

    return ret;
}

inline int MBPTCalculator::kmin(const OrbitalInfo& a, const OrbitalInfo& b, const OrbitalInfo& c, const OrbitalInfo& d) const
{
    int opt1 = kmin(a, b);
    int opt2 = kmin(c, d);

    return mmax(opt1, opt2);
}

inline int MBPTCalculator::kmax(const OrbitalInfo& a, const OrbitalInfo& b) const
{
    return (a.TwoJ() + b.TwoJ())/2;
}

inline int MBPTCalculator::kmax(const OrbitalInfo& a, const OrbitalInfo& b, const OrbitalInfo& c, const OrbitalInfo& d) const
{
    return mmin(a.TwoJ() + b.TwoJ(), c.TwoJ() + d.TwoJ())/2;
}

#endif
