#include "Include.h"
#include "Universal/Enums.h"
#include "Universal/Solution.h"

SolutionID::SolutionID(double aJ, ParityType::Enum aParity, unsigned int aID)
{
    mJ = aJ;
    mParity = aParity;
    mID = aID;
}

bool SolutionID::operator<(const SolutionID& other) const
{
    if(mJ < other.mJ)
    {
        return true;
    }
    else if(mJ == other.mJ)
    {
        if(mParity < other.mParity)
        {
            return true;
        }
        else if(mParity == other.mParity)
        {
            if(mID < other.mID)
            {
                return true;
            }
        }
    }

    return false;
}

Solution::Solution(double aEnergy, double agFactor)
{
    mEnergy = aEnergy;
    mgFactor = agFactor;
}

Solution::Solution(double aEnergy, std::map<Configuration, double> aPercentagesMap, double agFactor)
{
    mEnergy = aEnergy;
    mgFactor = agFactor;
    
    std::map<Configuration, double>::iterator cd_it;
    for(cd_it = aPercentagesMap.begin(); cd_it != aPercentagesMap.end(); cd_it++)
    {
        if(cd_it->second > 1)
        {
            mConfigurationSet.insert(ConfigurationPair(cd_it->first, cd_it->second));
        }
    }
}

void SolutionMap::Print()
{
    SolutionMap::iterator sm_it;
    for(sm_it = begin(); sm_it != end(); sm_it++)
    {
        *outstream << "J = " << sm_it->first.GetJ() << " Parity = " << ParityType::Name(sm_it->first.GetParity()) << " ID = " << sm_it->first.GetID() << std::endl;
        *outstream << "E = " << sm_it->second.GetEnergyAtomicUnits() << std::endl;
        sm_it->second.GetConfigurationSet()->Print();
    }
}
