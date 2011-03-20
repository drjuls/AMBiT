#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include "Atom/Atom.h"
#include "Atom/RateCalculator.h"
#include "Configuration/Eigenstates.h"
#include "Configuration/Symmetry.h"
#include "HartreeFock/OrbitalInfo.h"
#include "Universal/Enums.h"

// Class storage for a transition type (E1, M1, E2 etc.) based on std::pair
class TransitionType : public std::pair<MultipolarityType::Enum, unsigned int>
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
    {   return first;
    }
    inline unsigned int GetMultipole()
    {   return second;
    }
    
    bool operator<(TransitionType& other);
    inline bool operator==(TransitionType& other)
    {    return ((GetType() == other.GetType()) && (GetMultipole() == other.GetMultipole()));
    }
};


class Transition
{
public:
    Transition(Atom* aAtom, TransitionType aTransitionType, Symmetry aSymmetryFrom, unsigned int aSolutionIDFrom, Symmetry aSymmetryTo, unsigned int aSolutionIDTo, TransitionGaugeType::Enum aGauge = TransitionGaugeType::Length);
    Transition(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry aSymmetryFrom, unsigned int aSolutionIDFrom, Configuration* aLeadingConfigurationTo, Symmetry aSymmetryTo, unsigned int aSolutionIDTo, TransitionGaugeType::Enum aGauge = TransitionGaugeType::Length);

    inline Atom* GetAtom() const
    {   return mAtom;
    }
    inline Configuration* GetLeadingConfigurationFrom() const
    {   return mLeadingConfigurationFrom;
    }
    inline Configuration* GetLeadingConfigurationTo() const
    {   return mLeadingConfigurationTo;
    }
    inline TransitionType GetTransitionType() const
    {   return mTransitionType;
    }
    inline double GetTransitionProbability() const
    {   return mTransitionProbability;
    }
    inline double GetTransitionRate() const
    {   return mTransitionRate;
    }
    inline Symmetry GetSymmetryFrom() const
    {   return mSymmetryFrom;
    }
    inline Symmetry GetSymmetryTo() const
    {   return mSymmetryTo;
    }
    inline unsigned int GetSolutionIDFrom() const
    {   return mSolutionIDFrom;
    }
    inline unsigned int GetSolutionIDTo() const
    {   return mSolutionIDTo;
    }

    inline void SetTransitionRate(double aTransitionRate)
    {   mTransitionRate = aTransitionRate;
    }

    bool operator<(const Transition& other) const;
    bool operator==(const Transition& other) const;

protected:
    Atom* mAtom;

    Configuration* mLeadingConfigurationFrom;
    Configuration* mLeadingConfigurationTo;
    TransitionType mTransitionType;

    Symmetry mSymmetryFrom;
    unsigned int mSolutionIDFrom;
    Symmetry mSymmetryTo;
    unsigned int mSolutionIDTo;

    double mTransitionProbability;
    double mTransitionRate;

    void Solve(TransitionGaugeType::Enum aGauge);
};

class TransitionSet : std::set<Transition>
{
public:
    TransitionSet(Atom* aAtom, TransitionType aTransitionType, TransitionGaugeType::Enum aGauge = TransitionGaugeType::Length);
    // Creates a set of Transitions with a specified Type, Symmetry and Configurations from and to
    TransitionSet(Atom* aAtom, TransitionType aTransitionType, Configuration* aLeadingConfigurationFrom, Symmetry aSymmetryFrom, Configuration* aLeadingConfigurationTo, Symmetry aSymmetryTo ,TransitionGaugeType::Enum aGauge = TransitionGaugeType::Length);
protected:

};

#endif
