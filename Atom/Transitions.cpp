#include <string>
#include <sstream>

#include "Atom/Transitions.h"

TransitionType::TransitionType()
{
    first = MultipolarityType::E;
    second = 1;
}

TransitionType::TransitionType(int aType, unsigned int aMultipole)
{
    if(aType == 0)
    {
        first = MultipolarityType::E;
    } 
    else if(aType ==1)
    {
        first = MultipolarityType::M;
    }
    else
    {
        *errstream << "ERROR: Unknown multipole type: " << aType << std::endl;
        exit(1);
    }
    
    second = aMultipole;
}

TransitionType::TransitionType(MultipolarityType::Enum aType, unsigned int aMultipole)
{
    first = aType;
    second = aMultipole;
}

TransitionType::TransitionType(std::string astring)
{
    assert(StringSpecifiesTransitionType(astring));
    
    if(astring[0] == 'E' || astring[0] == 'e')
    {
        first = MultipolarityType::E;
    }
    if(astring[0] == 'M' || astring[0] == 'm')
    {
        first = MultipolarityType::M;
    }
    second = atoi(astring.substr(1).c_str());
}

bool TransitionType::ChangesParity()
{
    bool parity_change = false;

    if(GetType() == MultipolarityType::M)
    {
        !parity_change;
    }
    if(GetMultipole()%2)
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
        if((aSymFrom.GetJ() <= ((GetMultipole() - 1)/2)) && aSymTo.GetJ() <= (GetMultipole() - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;    
    }
    else if((aSymFrom.GetParity() != aSymTo.GetParity()) && ChangesParity())
    {
        if((aSymFrom.GetJ() <= ((GetMultipole() - 1)/2)) && aSymTo.GetJ() <= (GetMultipole() - aSymFrom.GetJ() - 1))
        {
            return false;
        }
        return true;
    }

    return false;
}

bool TransitionType::StringSpecifiesTransitionType(std::string astring)
{
    if(astring.size() < 2)
    {
        return false;
    }

    if(astring[0] == 'E' || astring[0] == 'M' || astring[0] == 'e' || astring[0] == 'm')
    {
        if(atoi(astring.substr(1).c_str()) > 0)
        {
            return true;
        }
    }
    
    return false;
}

bool TransitionType::operator<(TransitionType& other)
{
    if(GetType() == other.GetType())
    {    return GetMultipole()< other.GetMultipole();
    }
    return (GetType() < other.GetType());
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

    ss << GetMultipole();

    name += MultipolarityType::Name(GetType());
    name += ss.str();

    return name;
}

Transition::Transition(Atom* aAtom, TransitionType aTransitionType, Symmetry aSymmetryFrom, unsigned int aSolutionIDFrom, Symmetry aSymmetryTo, unsigned int aSolutionIDTo, TransitionGaugeType::Enum aGauge) :
mSymmetryFrom(0, even), mSymmetryTo(0, even)
{
    mAtom = aAtom;
    mTransitionType = aTransitionType;
    mSymmetryFrom = aSymmetryFrom;
    mSolutionIDFrom = aSolutionIDFrom;
    mSymmetryTo = aSymmetryTo;
    mSolutionIDTo = aSolutionIDTo;

    // Add code to determine the leading configurations...
    // mLeadingConfigurationFrom =
    // mLeadingConfigurationTo =

    Solve(aGauge);
}

Transition::Transition(Atom* aAtom, TransitionType aTransitionType, SolutionID aSolutionIDFrom, SolutionID aSolutionIDTo, TransitionGaugeType::Enum aGauge) :
mSymmetryFrom(0, even), mSymmetryTo(0, even)
{
    mAtom = aAtom;
    mTransitionType = aTransitionType;
    mSymmetryFrom = aSolutionIDFrom.GetSymmetry();
    mSolutionIDFrom = aSolutionIDFrom.GetID();
    mSymmetryTo = aSolutionIDTo.GetSymmetry();
    mSolutionIDTo = aSolutionIDTo.GetID();

    // Add code to determine the leading configurations...
    // mLeadingConfigurationFrom =
    // mLeadingConfigurationTo =

    Solve(aGauge);
}

Transition::Transition(Atom* aAtom, SolutionID aSolutionIDFrom, SolutionID aSolutionIDTo, TransitionGaugeType::Enum aGauge) :
mSymmetryFrom(0, even), mSymmetryTo(0, even)
{
    mAtom = aAtom;
    mSymmetryFrom = aSolutionIDFrom.GetSymmetry();
    mSolutionIDFrom = aSolutionIDFrom.GetID();
    mSymmetryTo = aSolutionIDTo.GetSymmetry();
    mSolutionIDTo = aSolutionIDTo.GetID();
    
    // Determine dominant transition type
    double JFrom = mSymmetryFrom.GetJ();
    double JTo = mSymmetryTo.GetJ();
    double ChangeJ = fabs(JFrom - JTo);
    
    bool ParityChange = (mSymmetryFrom.GetParity() != mSymmetryTo.GetParity());
    unsigned int multipole = ChangeJ;
    if(multipole == 0)
    {
        multipole = 1;
    }
    
    if(JTo == 0.0 && JFrom == 0.0)
    {
        *outstream << "Warning! Attempting to calculate transition from J = 0 to J = 0." << std::endl;
    }
    
    MultipolarityType::Enum mType;
    if((pow(-1.0, multipole) == 1.0 && !ParityChange) || (pow(-1.0, multipole) == -1.0 && ParityChange))
    {
        mType = MultipolarityType::E;
    }
    else
    {
        mType = MultipolarityType::M;
    }
    
    mTransitionType = TransitionType(mType, multipole);
    
    // Add code to determine the leading configurations...
    // mLeadingConfigurationFrom =
    // mLeadingConfigurationTo =
    
    Solve(aGauge);
}

void Transition::Solve(TransitionGaugeType::Enum aGauge)
{
    RateCalculator rcalc(GetAtom()->GetBasis());
    SetTransitionRate(rcalc.CalculateMultipoleStrength(GetTransitionType().GetType(), GetTransitionType().GetMultipole(), GetAtom(), GetSymmetryFrom(), GetSolutionIDFrom(), GetSymmetryTo(), GetSolutionIDTo()));
    
    if(GetTransitionType() == TransitionType(MultipolarityType::E, 1))
    {
        rcalc.CalculateDipoleStrength(GetAtom(), GetSymmetryFrom(), GetSolutionIDFrom(), GetSymmetryTo(), GetSolutionIDTo());
    }
    /*
    if(GetTransitionType() == TransitionType(MultipolarityType::E, 1))
    {
        RateCalculator rcalc(GetAtom()->GetBasis());
        SetTransitionRate(rcalc.CalculateMultipoleStrength(MultipolarityType::M, 1, GetAtom(), GetSymmetryFrom(), GetSolutionIDFrom(), GetSymmetryTo(), GetSolutionIDTo()));
        //rcalc.CalculateDipoleStrength(GetAtom(), GetSymmetryFrom(), GetSolutionIDFrom(), GetSymmetryTo(), GetSolutionIDTo());
    }
    else
    {
        SetTransitionRate(0);
    }*/
}

bool Transition::operator<(const Transition& other) const
{
    return (GetTransitionType() < other.GetTransitionType());
}

bool Transition::operator==(const Transition& other) const
{
    bool isEqual = (GetAtom() == other.GetAtom());
    isEqual = isEqual && (GetTransitionType() == other.GetTransitionType());
    isEqual = isEqual && (GetSymmetryFrom() == other.GetSymmetryFrom());
    isEqual = isEqual && (GetSolutionIDFrom() == other.GetSolutionIDFrom());
    isEqual = isEqual && (GetSymmetryTo() == other.GetSymmetryTo());
    isEqual = isEqual && (GetSolutionIDTo() == other.GetSolutionIDTo());
    
    return isEqual;
}

/*
TransitionSet::TransitionSet(Atom* aAtom, TransitionType aTransitionType, TransitionGaugeType::Enum aGauge)
{
    SymmetryEigenstatesMap::const_iterator from_it, to_it;
    Symmetry CurrentSymmetryFrom(0, even), CurrentSymmetryTo(0, even);
    int i, j;

    CurrentSymmetryFrom = aAtom->GetSymmetryEigenstatesMap()->begin()->first;
    CurrentSymmetryTo = CurrentSymmetryFrom;

    for(from_it = aAtom->GetSymmetryEigenstatesMap()->begin(), i = 0; from_it != aAtom->GetSymmetryEigenstatesMap()->end(); from_it++, i++)
    {
        if(!(from_it->first == CurrentSymmetryFrom))
        {
            i = 0;
            CurrentSymmetryFrom = from_it->first;
        }
        for(to_it = from_it, j = i, CurrentSymmetryTo = CurrentSymmetryFrom; to_it != aAtom->GetSymmetryEigenstatesMap()->end(); to_it++, j++)
        {
            if(!(to_it->first == CurrentSymmetryTo))
            {
                j = 0;
                CurrentSymmetryTo = to_it->first;
            }
            if(aTransitionType.IsAllowedTransition(from_it->first, to_it->first))
            {
                insert(Transition(aAtom, aTransitionType, Symmetry(from_it->first), (unsigned int) i, Symmetry(to_it->first), (unsigned int) j, aGauge));
            }
        }
    }
}

TransitionSet::TransitionSet(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry aSymmetryFrom, Configuration* aLeadingConfigurationTo, Symmetry aSymmetryTo ,TransitionGaugeType::Enum aGauge)
{
    if(!aTransitionType.IsAllowedTransition(aSymmetryFrom, aSymmetryTo))
    {
        *errstream << aTransitionType.Name() <<" transitions are not allowed by selection rules between " << aSymmetryFrom.GetString() << " and " << aSymmetryTo.GetString() << std::endl;
    }
    else
    {
        int iterator;
        Eigenstates* EigenstatesFrom = aAtom->GetEigenstates(aSymmetryFrom);
        Eigenstates* EigenstatesTo = aAtom->GetEigenstates(aSymmetryTo);
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
            if(*aLeadingConfigurationFrom == from_it->second)
            {
                for(to_it = ToMap.begin(); to_it != ToMap.end(); to_it++)
                {
                    if(*aLeadingConfigurationTo == to_it->second)
                    {
                        insert(Transition(aAtom, aTransitionType, aSymmetryFrom, (unsigned int) from_it->first, aSymmetryTo, (unsigned int) to_it->first, aGauge));
                        found_transitions = true;
                    }
                }
            }
        }
    }
}
*/