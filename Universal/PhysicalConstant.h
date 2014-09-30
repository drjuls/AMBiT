#ifndef PHYSICAL_CONSTANT_H
#define PHYSICAL_CONSTANT_H

#include <memory>

/** Fundamental physical constants.
    Note that AMBiT uses atomic units, hbar = m_e = e = 1, so
      c = 1/alpha ~ 137 depends on alpha.
 */
class PhysicalConstant
{
public:
    PhysicalConstant();
    ~PhysicalConstant() {}

    double GetAlpha() const;
    double GetAlphaSquared() const;

    /** Speed of light = 1/alpha. */
    double GetSpeedOfLight() const;

    /** Nucleon to electron mass ratio (remember m_e = 1). */
    double GetNucleonMass() const;

    /** Set alpha^2 = alpha_0^2 (1 + ratio), where alpha_0 is the Earth value. */
    void SetAlphaSquaredIncreaseRatio(double ratio);

    /** Set alpha = alpha_0, where alpha_0 is the Earth value. */
    void SetAlphaToEarthValue();
    
    /** Set alpha to any value. */
    void SetAlpha(double alpha);

protected:
    double AlphaSquared;
    double NucleonElectronMassRatio;
};

typedef std::shared_ptr<PhysicalConstant> pPhysicalConstant;
typedef std::shared_ptr<const PhysicalConstant> pPhysicalConstantConst;

#endif
