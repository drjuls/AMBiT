#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include "Atom/Atom.h"
#include "Configuration/Level.h"
#include "Universal/Enums.h"
#include <ctype.h>
#include <map>

/** Type of electromagnetic transition, e.g. E1, M1, E2, ...
    TransitionType is ordered and stream operators are provided.
 */
typedef std::pair<MultipolarityType, int> TransitionType;

inline bool operator<(const TransitionType& left, const TransitionType& right)
{
    if(left.second < right.second)
        return true;
    else if(left.second > right.second)
        return false;

    if(left.first != right.first)
        return (left.first == MultipolarityType::E);
    return false;
}

inline std::ostream& operator<<(std::ostream& outstream, const TransitionType& transition_type)
{
    outstream << Name(transition_type.first) << itoa(transition_type.second);
    return outstream;
}

inline std::istream& operator>>(std::istream& instream, TransitionType& transition_type)
{
    char c;
    do{
        instream.get(c);
    } while(std::isspace(c));

    if(c == 'e' || c == 'E')
        transition_type.first = MultipolarityType::E;
    else if(c == 'm' || c == 'M')
        transition_type.first = MultipolarityType::M;

    instream >> transition_type.second;
    return instream;
}

class TransitionTypeComparator
{
public:
    bool operator()(const TransitionType& left, const TransitionType& right) const { return (left < right); }
};

/** LevelID uniquely specifies a Level (HamiltonianID and index).
    The principle of don't be a jerk programming applies: for any calculation we should only be using one
    type of HamiltonianID in Atom, so no comparisons between pHamiltonianIDs of different classes occurs.
 */
typedef std::pair<pHamiltonianID, int> LevelID;
LevelID make_LevelID(const std::string& name);
std::string Name(const LevelID& levelid);

/** A transition includes a start level, final level, and transition type. */
typedef std::tuple<LevelID, LevelID, TransitionType> TransitionID;

class TransitionIDComparator
{
public:
    bool operator()(const TransitionID& first, const TransitionID& second) const;
};

/** Transitions calculates transition strengths and stores them in a mapping
    between TransitionID and Strength (S).
 */
class TransitionMap : protected std::map<TransitionID, double, TransitionIDComparator>
{
protected:
    typedef std::map<TransitionID, double, TransitionIDComparator> Parent;

public:
    TransitionMap(Atom& atom, TransitionGauge gauge = TransitionGauge::Length, TransitionType max = std::make_pair(MultipolarityType::E, 3));

    typedef Parent::iterator iterator;
    typedef Parent::const_iterator const_iterator;

    using Parent::at;
    using Parent::begin;
    using Parent::clear;
    using Parent::count;
    using Parent::empty;
    using Parent::end;
    using Parent::find;
    using Parent::insert;
    using Parent::size;

    /** Check if transition exists or is excluded by symmetry considerations. */
    bool TransitionExists(const Symmetry& left, const Symmetry& right, TransitionType type) const;
    bool TransitionExists(const LevelID& left, const LevelID& right, TransitionType type) const
    {   return TransitionExists(left.first->GetSymmetry(), right.first->GetSymmetry(), type);
    }

    /** Calculate transition strength for smallest TransitionType possible and add to map.
        If smallest TransitionType > max_type, return zero.
     */
    double CalculateTransition(const LevelID& left, const LevelID& right);

    /** Calculate transition strength if allowed by type and add to map.
        While we're at it, calculate all transitions with same symmetry, since this hardly costs any more.
        Print transition that was requested.
     */
    double CalculateTransition(const LevelID& left, const LevelID& right, TransitionType type);

    void Print() const;

protected:
    /** TransitionID with left <= right. */
    TransitionID make_transitionID(const LevelID& left, const LevelID& right, TransitionType type) const;

protected:
    Atom& atom;
    pOrbitalManagerConst orbitals;
    pLevelStore levels;

    TransitionGauge preferred_gauge;
    TransitionType max_type;

    std::map<TransitionType, pTransitionIntegrals> integrals;
};

inline bool TransitionIDComparator::operator()(const TransitionID& first, const TransitionID& second) const
{
    if(std::get<0>(first) < std::get<0>(second))
        return true;
    else if(std::get<0>(first) > std::get<0>(second))
        return false;

    if(std::get<1>(first) < std::get<1>(second))
        return true;
    else if(std::get<1>(first) > std::get<1>(second))
        return false;

    return (std::get<2>(first) < std::get<2>(second));
}

#endif
