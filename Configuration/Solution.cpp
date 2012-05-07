#include <map>
#include <sstream>

#include "Configuration/Solution.h"

SolutionID::SolutionID(double aJ, ParityType::Enum aParity, unsigned int aID)
{
    mJ = aJ;
    mParity = aParity;
    mID = aID;
}

SolutionID::SolutionID(Symmetry aSymmetry, unsigned int aID)
{
    mJ = aSymmetry.GetJ();
    if(aSymmetry.GetParity() == even)
    {
        mParity = ParityType::Even;
    }
    else
    {
        mParity = ParityType::Odd;
    }    
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

std::string SolutionID::GetIdentifier() const
{
    std::stringstream ss;

    ss << (int) 2 * mJ << ParityType::VeryShortName(mParity) <<  mID;

    return ss.str();
}

Symmetry SolutionID::GetSymmetry() const
{
    if(mParity == ParityType::Even)
    {
        return Symmetry(mJ * 2, even);
    }
    else
    {
        return Symmetry(mJ * 2, odd);
    }
}

Solution::Solution(double aEnergy, double agFactor)
{
    mEnergy = aEnergy;
    mgFactor = agFactor;
    mConfigurationSet = new ConfigurationSet();
    mTransitionSet = new TransitionSet();
}

Solution::Solution(double aEnergy, std::map<Configuration, double> aPercentagesMap, double agFactor)
{
    mEnergy = aEnergy;
    mgFactor = agFactor;
    mConfigurationSet = new ConfigurationSet();
    mTransitionSet = new TransitionSet();

    std::map<Configuration, double>::iterator cd_it;
    for(cd_it = aPercentagesMap.begin(); cd_it != aPercentagesMap.end(); cd_it++)
    {
        if(cd_it->second > 1)
        {
            mConfigurationSet->insert(ConfigurationPair(cd_it->first, cd_it->second));
        }
    }
}

void SolutionMap::Print(DisplayOutputType::Enum aDisplayOutputType)
{
    if(aDisplayOutputType == DisplayOutputType::Standard)
    {
        SolutionMap::iterator sm_it;
        SolutionID LastSolutionID(-1, ParityType::Even, 0);
        for(sm_it = begin(); sm_it != end(); sm_it++)
        {
            if(!(sm_it->first.GetSymmetry() == LastSolutionID.GetSymmetry()))
            {
                *outstream << "Solutions for J = " << sm_it->first.GetJ() << ", P = " << ParityType::LowerName(sm_it->first.GetParity()) << ":" << std::endl;
            }
            *outstream << sm_it->first.GetID() << ": " << std::setprecision(8) << sm_it->second.GetEnergyAtomicUnits() << "    " << std::setprecision(12) << sm_it->second.GetEnergyInversecm() << " /cm" << std::endl;
            sm_it->second.GetConfigurationSet()->Print();
            if(sm_it->first.GetJ() != 0)
            {
                *outstream << "    g-factor = " << std::setprecision(5) << sm_it->second.GetgFactor() << std::endl;
            }
            *outstream << std::endl;
            LastSolutionID = sm_it->first;
        }
    }
    else if(aDisplayOutputType == DisplayOutputType::SpaceSeparated)
    {
        *outstream << "J P ID E g" << std::endl;
        SolutionMap::iterator sm_it;
        for(sm_it = begin(); sm_it != end(); sm_it++)
        {
            *outstream << sm_it->first.GetJ() << " " << ParityType::ShortName(sm_it->first.GetParity()) << " " << sm_it->first.GetID() << " " << std::setprecision(12) << sm_it->second.GetEnergyInversecm();
            *outstream << " " << std::setprecision(5) << sm_it->second.GetgFactor();
            // Add error checking for GetLargestConfiguration!
            *outstream << " " << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.ShortName() << " " << std::setprecision(2) << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->second << "%" << std::endl;
        }
        *outstream << std::endl;
    }
    else if(aDisplayOutputType == DisplayOutputType::CommaSeparated)
    {
        *outstream << "J,P,ID,E,g" << std::endl;
        SolutionMap::iterator sm_it;
        for(sm_it = begin(); sm_it != end(); sm_it++)
        {
            *outstream << sm_it->first.GetJ() << "," << ParityType::ShortName(sm_it->first.GetParity()) << "," << sm_it->first.GetID() << "," << std::setprecision(12) << sm_it->second.GetEnergyInversecm();
            *outstream << "," << std::setprecision(5) << sm_it->second.GetgFactor();
            // Add error checking for GetLargestConfiguration!
            *outstream << "," << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.Name(false) << "," << std::setprecision(2) << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->second << "%" << std::endl;
        }
        *outstream << std::endl;
    }
    else if(aDisplayOutputType == DisplayOutputType::TabSeparated)
    {
        *outstream << "J\tP\tID\tE\tg" << std::endl;
        SolutionMap::iterator sm_it;
        for(sm_it = begin(); sm_it != end(); sm_it++)
        {
            *outstream << sm_it->first.GetJ() << "\t" << ParityType::ShortName(sm_it->first.GetParity()) << "\t" << sm_it->first.GetID() << "\t" << std::setprecision(12) << sm_it->second.GetEnergyInversecm();
            *outstream << "\t" << std::setprecision(5) << sm_it->second.GetgFactor();
            // Add error checking for GetLargestConfiguration!
            *outstream << "\t" << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.Name(false) << "\t" << std::setprecision(2) << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->second << "%" << std::endl;
        }
        *outstream << std::endl;
    }
}

void SolutionMap::PrintID()
{
    SolutionMap::iterator sm_it;
    for(sm_it = begin(); sm_it != end(); sm_it++)
    {
        *outstream << sm_it->first.GetIdentifier() << " " << std::setprecision(12) << sm_it->second.GetEnergyInversecm() << " " << sm_it->second.GetLeadingConfiguration().Name(false) << std::endl;
    }
}

SolutionMap::iterator SolutionMap::FindByIdentifier(const std::string& aIdentifier)
{
    const char* IdentifierString = aIdentifier.c_str();

    int i = 0;
    int TwoJ = atoi(IdentifierString);
    ParityType::Enum Parity;
    int ID = 0;

    while(isdigit(IdentifierString[i]) && IdentifierString[i])
    {
        i++;
    }
    if(IdentifierString[i] == 'e')
    {
        Parity = ParityType::Even;
    }
    else if(IdentifierString[i] == 'o')
    {
        Parity = ParityType::Odd;
    }
    else
    {
        *errstream << "Invalid solution identifier string received: " << IdentifierString << std::endl;
        exit(1);
    }
    i++;
    ID = atoi(IdentifierString + i);

    return find(SolutionID((double) TwoJ/2, Parity, ID));
}

void SolutionMap::PrintSolution(SolutionMap::iterator aSolutionIterator)
{
    if(aSolutionIterator == end())
    {
        *errstream << "Solution not found!" << std::endl;
        exit(1);
    }
    else
    {
        *outstream << aSolutionIterator->first.GetIdentifier() << " " << std::setprecision(12) << aSolutionIterator->second.GetEnergyInversecm() << " " << aSolutionIterator->second.GetLeadingConfiguration().Name(false) << std::endl;
    }
}

void SolutionMapMap::Print()
{
    *outstream << "Solution summary for " << size() << " runs." << std::endl;
    SolutionMapMap::iterator smm_it, smm_first;
    *outstream << "J,P,ID,";
    for(smm_it = begin(); smm_it != end(); smm_it++)
    {
        *outstream << "," << "E" << smm_it->first;
    }
    for(smm_it = begin(); smm_it != end(); smm_it++)
    {
        *outstream << "," << "g" << smm_it->first;
    }
    *outstream << ",Config" << std::endl;
    
    
    smm_first = begin();
    SolutionMap::iterator sm_it;
    // Use labelling from first run
    for(sm_it = smm_first->second.begin(); sm_it != smm_first->second.end(); sm_it++)
    {
        *outstream << sm_it->first.GetJ() << "," << ParityType::ShortName(sm_it->first.GetParity()) << "," << sm_it->first.GetID();
        for(smm_it = begin(); smm_it != end(); smm_it++)
        {
            *outstream  << "," << std::setprecision(12) << (smm_it->second.find(sm_it->first))->second.GetEnergyInversecm();
        }
        for(smm_it = begin(); smm_it != end(); smm_it++)
        {
            *outstream  << "," << std::setprecision(5) << (smm_it->second.find(sm_it->first))->second.GetgFactor();
        }
        *outstream << "," << sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.Name(false);
        bool difference_found = false;
        for(smm_it = begin(); smm_it != end(); smm_it++)
        {
            if(sm_it->second.GetgFactor() != 0)
            {
                if(fabs(((smm_it->second.find(sm_it->first))->second.GetgFactor() - sm_it->second.GetgFactor())/sm_it->second.GetgFactor()) > 0.1)
                {
                    difference_found = true;
                }
            }
            if((smm_it->second.find(sm_it->first))->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.Name(false) != sm_it->second.GetConfigurationSet()->GetLargestConfigurationPair()->first.Name(false))
            {
                difference_found = true;
            }
        }
        if(difference_found)
        {
            *outstream << "*";
        }
        
        *outstream << std::endl;
    }
    
    *outstream << std::endl;
}

