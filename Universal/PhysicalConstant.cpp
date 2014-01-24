#include "Include.h"
#include "PhysicalConstant.h"

#define ALPHA_0        0.007297352932703
#define NUCLEON_MASS_0 1822.888

PhysicalConstant* PhysicalConstant::Instance()
{
    static PhysicalConstant instance;
    return &instance;
}

PhysicalConstant::PhysicalConstant()
{
    NucleonElectronMassRatio = 1822.888;
    SetAlphaToEarthValue();
}

PhysicalConstant::~PhysicalConstant()
{}

double PhysicalConstant::GetAlpha() const
{
    return sqrt(AlphaSquared);
}

double PhysicalConstant::GetAlphaSquared() const
{
    return AlphaSquared;
}

double PhysicalConstant::GetSpeedOfLight() const
{
    return 1.0/sqrt(AlphaSquared);
}

double PhysicalConstant::GetNucleonMass() const
{
    return NucleonElectronMassRatio;
}

void PhysicalConstant::SetAlphaSquaredIncreaseRatio(double ratio)
{
    AlphaSquared = ALPHA_0 * ALPHA_0 * (1.0 + ratio);
}

void PhysicalConstant::SetAlphaToEarthValue()
{
    SetAlphaSquaredIncreaseRatio(0.0);
}

void PhysicalConstant::SetAlpha(double alpha)
{
    AlphaSquared = alpha * alpha;
}
