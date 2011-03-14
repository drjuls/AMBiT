#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include "Atom/Atom.h"
#include "Configuration/Eigenstates.h"
#include "Configuration/Symmetry.h"
#include "HartreeFock/StateInfo.h"
#include "Universal/Enums.h"

class TransitionType
{
public:
    TransitionType();
    TransitionType(int aType, unsigned int aMultipole);
    TransitionType(MultipolarityType::Enum aType, unsigned int aMultipole);
    
    bool ChangesParity();

    virtual std::string Name();
    bool IsAllowedTransition(Symmetry aSymFrom, Symmetry aSymTo);
    //static bool ExistsBetweenStates(Symmetry aSymFrom, Symmetry aSymTo, TransitionType aType);
    
    inline MultipolarityType::Enum GetType()
    {
        return mType;
    }
    inline unsigned int GetMultipole()
    {
        return mMultipole;
    }

private:
    MultipolarityType::Enum mType;
    unsigned int mMultipole;
};

class Transition
{
public:
    Transition(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry* aSymmetryFrom, unsigned int aSolutionIDFrom, Configuration* aLeadingConfigurationTo, Symmetry* aSymmetryTo, unsigned int aSolutionIDTo, TransitionGaugeType::Enum aGauge);

    inline Configuration* GetLeadingConfigurationFrom()
    {   return mLeadingConfigurationFrom;
    }
    inline Configuration* GetLeadingConfigurationTo()
    {   return mLeadingConfigurationTo;
    }
    inline TransitionType GetTransitionType()
    {   return mTransitionType;
    }
    inline double GetTransitionProbability()
    {   return mTransitionProbability;
    }
    inline double GetTransitionRate()
    {   return mTransitionRate;
    }
    inline Symmetry* GetSymmetryFrom()
    {   return mSymmetryFrom;
    }
    inline Symmetry* GetSymmetryTo()
    {   return mSymmetryTo;
    }
    inline unsigned int GetSolutionIDFrom()
    {   return mSolutionIDFrom;
    }
    inline unsigned int GetSolutionIDTo()
    {   return mSolutionIDTo;
    }
private:
    Atom* mAtom;

    Configuration* mLeadingConfigurationFrom;
    Configuration* mLeadingConfigurationTo;
    TransitionType mTransitionType;

    Symmetry* mSymmetryFrom;
    unsigned int mSolutionIDFrom;
    Symmetry* mSymmetryTo;
    unsigned int mSolutionIDTo;

    double mTransitionProbability;
    double mTransitionRate;

    void Solve(TransitionGaugeType::Enum aGauge);
};

class TransitionSet : std::set<Transition>
{
public:
    TransitionSet(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry* aSymmetryFrom, Configuration* aLeadingConfigurationTo, Symmetry* aSymmetryTo ,TransitionGaugeType::Enum aGauge);
private:

};

#endif
