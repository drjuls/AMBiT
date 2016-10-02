#ifndef MBPT_CALCULATOR_H
#define MBPT_CALCULATOR_H

#include "Include.h"
#include "HartreeFock/Core.h"
#include "Basis/OrbitalManager.h"

class MBPTCalculator
{
    /** Calculate diagrams of many-body perturbation theory.
        Use Brillouin-Wigner perturbation theory, where the energy of external lines is
        kept constant in the energy denominator (this ensures that the operator is hermitian).
     */
public:
    /** Optional: fermi_orbitals are the orbitals used for energy denominators of MBPT.
            Denote using a basis-style string, e.g. 5sp4df.
            Any unspecified angular momentum defaults to first orbital above fermi level.
     */
    MBPTCalculator(pOrbitalManagerConst pOrbitals, const std::string& fermi_orbitals = "", bool include_off_parity = false);
    virtual ~MBPTCalculator(void);

    virtual unsigned int GetStorageSize() = 0;
    virtual void UpdateIntegrals() = 0;

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

    /** Return true if one of the orbital info arguments is not in the valence set. */
    template <typename... OrbitalTypes>
    inline bool InQSpace(const OrbitalInfo& orb, const OrbitalTypes&... orbitals) const
    {
        return ((valence->count(orb) == 0) || InQSpace(orbitals...));
    }

    inline bool InQSpace(const OrbitalInfo& orb) const
    {
        return (valence->count(orb) == 0);
    }

protected:
    /** Return absolute value of the difference of two unsigned integers. */
    inline int absdiff(int i, int j) const;

    /** Parity check returns true if (a.L() + b.L() + k)%2 == 0 or include_off_parity. */
    inline bool ParityCheck(const OrbitalInfo& a, const OrbitalInfo& b, const int& k) const;

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
    pOrbitalManagerConst orbitals;
    pOrbitalMapConst valence;

    /** Two-body radial integrals may include off-parity k: (l1 + l2 + k) is odd.
        Set kstep to one to include these.
     */
    const bool include_off_parity;
    const int kstep;

    /** Valence energies for Brillouin-Wigner MBPT. */
    std::map<int, double> ValenceEnergies;

    /** Choose orbitals used to create ValenceEnergies for energy denominators of MBPT. */
    std::string fermi_orbitals;
    void SetValenceEnergies();

    double delta;   // Shift in the energy denominator.
};

inline int MBPTCalculator::absdiff(int i, int j) const
{
    return abs(i - j);
}

inline bool MBPTCalculator::ParityCheck(const OrbitalInfo& a, const OrbitalInfo& b, const int& k) const
{
    return (include_off_parity || ((a.L() + b.L() + k)%2 == 0));
}

inline int MBPTCalculator::kmin(const OrbitalInfo& a, const OrbitalInfo& b) const
{
    // ensure triangle(a.J(), b.J(), kmin)
    int k = absdiff(a.TwoJ(), b.TwoJ())/2;

    // ensure eps(a.L() + b.L() + kmin)
    if(!ParityCheck(a, b, k))
        k++;

    return k;
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
