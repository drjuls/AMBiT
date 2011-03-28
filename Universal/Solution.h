#ifndef SOLUTION_H
#define SOLUTION_H

#include <map>
#include "Configuration/Configuration.h"
#include "Universal/Constant.h"
#include "Universal/Enums.h"

// Effectively a hash function for converting
// J, Pi and solution position in the list
// To a unique ID
class SolutionID
{
public:
    SolutionID(double aJ, ParityType::Enum aParity, unsigned int aID);
    
    bool operator<(const SolutionID& other) const;
    inline bool operator==(const SolutionID& other)
    {   return ((mJ == other.mJ) && (mParity == other.mParity) && (mID == other.mID));
    }
    
    double GetJ() const { return mJ; }
    ParityType::Enum GetParity() const { return mParity; }
    unsigned int GetID() const { return mID; }

protected:
    double mJ;
    ParityType::Enum mParity;
    unsigned int mID;
};

class Solution
{
public:
    // Creates a Solution with energy aEnergy in atomic units
    Solution(double aEnergy, double agFactor = 0.0);
    Solution(double aEnergy, std::map<Configuration, double> aPercentagesMap, double agFactor = 0.0);

    inline double GetEnergyInversecm() const
    {   return mEnergy * Constant::HartreeEnergy_cm;
    }
    inline double GetEnergyAtomicUnits() const
    {   return mEnergy;
    }
    inline double GetgFactor() const
    {   return mgFactor;
    }
    inline ConfigurationSet* GetConfigurationSet()
    {   return &mConfigurationSet;
    }

protected:
    double mEnergy;
    double mgFactor;
    ConfigurationSet mConfigurationSet;
};

class SolutionMap : public std::map<SolutionID, Solution>
{
public:
    void Print();
protected:
};

#endif