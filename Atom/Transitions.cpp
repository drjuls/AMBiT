#include <string>
#include <sstream>

#include "Include.h"
#include "Atom/Atom.h"
#include "Atom/Transitions.h"
#include "Configuration/Eigenstates.h"
#include "Configuration/Symmetry.h"
#include "HartreeFock/StateInfo.h"
#include "Universal/Enums.h"

TransitionType::TransitionType()
{
    mType = MultipolarityType::E;
    mMultipole = 1;
}

TransitionType::TransitionType(int aType, unsigned int aMultipole)
{
    if(aType == 0)
    {
        mType = MultipolarityType::E;
    } 
    else if(aType ==1)
    {
        mType = MultipolarityType::M;
    }
    else
    {
        *errstream << "ERROR: Unknown multipole type: " << aType << std::endl;
        exit(1);
    }
    
    mMultipole = aMultipole;
}

TransitionType::TransitionType(MultipolarityType::Enum aType, unsigned int aMultipole)
{
    mType = aType;
    mMultipole = aMultipole;
}

bool TransitionType::ChangesParity()
{
    bool parity_change = false;

    if(mType == MultipolarityType::M)
    {
        !parity_change;
    }
    if(mMultipole%2)
    {
        !parity_change;
    }

    return parity_change;
}

bool TransitionType::IsAllowedTransition(Symmetry aSymFrom, Symmetry aSymTo)
{
    double DeltaJ = aSymFrom.GetJ() - aSymTo.GetJ();

    if((aSymFrom.GetParity() == aSymTo.GetParity()) && !ChangesParity())
    {
        if((aSymFrom.GetJ() <= ((mMultipole - 1)/2)) && aSymTo.GetJ() <= (mMultipole - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;    
    }
    else if((aSymFrom.GetParity() != aSymTo.GetParity()) && ChangesParity())
    {
        if((aSymFrom.GetJ() <= ((mMultipole - 1)/2)) && aSymTo.GetJ() <= (mMultipole - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;
    }

    return false;
}
     /*
static bool TransitionType::ExistsBetweenStates(Symmetry aSymFrom, Symmetry aSymTo, TransitionType aType)
{
    double DeltaJ = aSymFrom.GetJ() - aSymTo.GetJ();

    if((aSymFrom.GetParity() == aSymTo.GetParity()) && !aType.ChangesParity())
    {
        if((aSymFrom.GetJ() <= ((aType.GetMultipole() - 1)/2)) && aSymTo.GetJ() <= (aType.GetMultipole() - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;
    }
    else if((aSymFrom.GetParity() != aSymTo.GetParity()) && aType.ChangesParity())
    {
        if((aSymFrom.GetJ() <= ((aType.GetMultipole() - 1)/2)) && aSymTo.GetJ() <= (aType.GetMultipole() - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;
    }

    return false;
}           */

std::string TransitionType::Name()
{
    std::string name = "";
    std::stringstream ss;

    ss << mMultipole;

    name += MultipolarityType::Name(mType);
    name += ss.str();

    return name;
}

Transition::Transition(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry* aSymmetryFrom, unsigned int aSolutionIDFrom, Configuration* aLeadingConfigurationTo, Symmetry* aSymmetryTo, unsigned int aSolutionIDTo, TransitionGaugeType::Enum aGauge)
{
    mAtom = aAtom;
    mTransitionType = aTransitionType;
    mLeadingConfigurationFrom = aLeadingConfigurationFrom;
    mSymmetryFrom = aSymmetryFrom;
    mSolutionIDFrom = aSolutionIDFrom;
    mLeadingConfigurationTo = aLeadingConfigurationTo;
    mSymmetryTo = aSymmetryTo;
    mSolutionIDTo = aSolutionIDTo;

    Solve(aGauge);
}

void Transition::Solve(TransitionGaugeType::Enum aGauge)
{
    // Do some integrals, put solution into class

}

TransitionSet::TransitionSet(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry* aSymmetryFrom, Configuration* aLeadingConfigurationTo, Symmetry* aSymmetryTo ,TransitionGaugeType::Enum aGauge)
{
    if(!aTransitionType.IsAllowedTransition(*aSymmetryFrom, *aSymmetryTo))
    {
        *errstream << aTransitionType.Name() <<" transitions are not allowed by selection rules between " << aSymmetryFrom->GetString() << " and " << aSymmetryTo->GetString() << std::endl;
    }
    else
    {
        int iterator;
        Eigenstates* EigenstatesFrom = aAtom->GetEigenstates(*aSymmetryFrom);
        Eigenstates* EigenstatesTo = aAtom->GetEigenstates(*aSymmetryTo);
        std::map<int, Configuration> FromMap;
        std::map<int, Configuration> ToMap;
        int i, j;
        double largest_percentage;
        Configuration largest_configuration;
        std::map<Configuration, double>::const_iterator it_largest_percentage, it;
        RelativisticConfigList::const_iterator list_it;
        std::map<Configuration, double> percentages;

        for(i = 0; i < EigenstatesFrom->GetNumEigenvalues(); i++) {
            list_it = EigenstatesFrom->GetConfigGenerator()->GetRelConfigs()->begin();

            j = 0;
            while(list_it != EigenstatesFrom->GetConfigGenerator()->GetRelConfigs()->end())
            {
                Configuration nrconfig(list_it->GetNonRelConfiguration());
                if(percentages.find(nrconfig) == percentages.end())
                    percentages[nrconfig] = 0.;

                for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
                {
                    double coeff = EigenstatesFrom->GetEigenvectors()[i*EigenstatesFrom->GetEigenvectorLength() + j];
                    coeff = coeff * coeff * 100;

                    percentages[nrconfig] += coeff;
                    j++;
                }
    
                list_it++;
            }

            it_largest_percentage = percentages.begin();
            largest_percentage = 0.0;
    
            it = percentages.begin();
            while(it != percentages.end())
            {
                if(it->second > largest_percentage)
                {   it_largest_percentage = it;
                    largest_percentage = it->second;
                }
    
            }

            FromMap.insert(std::pair<int, Configuration>(i, it_largest_percentage->first));
        }

        for(i = 0; i < EigenstatesTo->GetNumEigenvalues(); i++) {
            list_it = EigenstatesTo->GetConfigGenerator()->GetRelConfigs()->begin();
    
            j = 0;
            while(list_it != EigenstatesTo->GetConfigGenerator()->GetRelConfigs()->end())
            {
                Configuration nrconfig(list_it->GetNonRelConfiguration());
                if(percentages.find(nrconfig) == percentages.end())
                    percentages[nrconfig] = 0.;

                for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
                {
                    double coeff = EigenstatesTo->GetEigenvectors()[i*EigenstatesTo->GetEigenvectorLength() + j];
                    coeff = coeff * coeff * 100;

                    percentages[nrconfig] += coeff;
                    j++;
                }
    
                list_it++;
            }

            it_largest_percentage = percentages.begin();
            largest_percentage = 0.0;
    
            it = percentages.begin();
            while(it != percentages.end())
            {
                if(it->second > largest_percentage)
                {   it_largest_percentage = it;
                    largest_percentage = it->second;
                }
    
            }

            ToMap.insert(std::pair<int, Configuration>(i, it_largest_percentage->first));
        }
        
        std::map<int, Configuration>::iterator from_it, to_it;
        bool found_transitions = false;
        for(from_it = FromMap.begin(); from_it != FromMap.end(); from_it++)
        {
            /*
            if(aLeadingConfigurationFrom->Name() == from_it->second.Name())
            {
                for(to_it = ToMap.begin(); to_it != ToMap.end(); to_it++)
                {
                    if(aLeadingConfigurationTo->Name() == to_it->second.Name())
                    {
                        insert(Transition(aAtom, aTransitionType, aLeadingConfigurationFrom, aSymmetryFrom, i, aLeadingConfigurationTo, aSymmetryTo, j, aGauge));
                        found_transitions = true;
                    }
                }
            }
            */
        }
    }
}
