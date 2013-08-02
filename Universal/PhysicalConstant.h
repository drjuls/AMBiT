#ifndef PHYSICAL_CONSTANT_H
#define PHYSICAL_CONSTANT_H

/** Fundamental physical constants, following the Singleton pattern.
    Note that AMBiT uses atomic units, hbar = m_e = e = 1, so
      c = 1/alpha ~ 137 depends on alpha.
 */
class PhysicalConstant
{
public:
    static PhysicalConstant* Instance();

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
    PhysicalConstant();
    ~PhysicalConstant();

protected:
    static PhysicalConstant* instance;

    double AlphaSquared;
    double NucleonElectronMassRatio;
};

#endif
