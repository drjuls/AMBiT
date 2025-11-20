#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include "Atom/Atom.h"
#include "Configuration/LevelVector.h"
#include "Universal/Enums.h"
#include "RPAOperator.h"
#include <map>

namespace Ambit
{
/** LevelID uniquely specifies a Level (HamiltonianID and index).
    The principle of don't be a jerk programming applies: for any calculation we should only be using one
    type of HamiltonianID in Atom, so no comparisons between pHamiltonianIDs of different classes occurs.
 */
typedef std::pair<pHamiltonianID, int> LevelID;
std::string Name(const LevelID& levelid);

typedef std::pair<LevelID, LevelID> TransitionID;

/** Abstract class for generating an operator from user_input, calculating strengths, and printing.
    Calculates matrix_elements for the operator and stores them in a mapping between TransitionID and Strength (S).
    Derived classes should implement
        - operator creation (in the constructor),
        - PrintHeader(), and
        - PrintTransition()
    TransitionCalculator provides generic RPA routines for authors of subclasses.
    When constructing the operator use something like
        if(user_input.search("--rpa"))
            op = MakeRPA(op, hf, hartreeY);
 */
class TransitionCalculator
{
public:
    TransitionCalculator(MultirunOptions& user_input, pOrbitalManagerConst orbitals, pLevelStore levels):
        user_input(user_input), orbitals(orbitals), levels(levels), op(nullptr), scale(1.0)
    {}
    TransitionCalculator(MultirunOptions& user_input, Atom& atom):
        TransitionCalculator(user_input, atom.GetBasis(), atom.GetLevels())
    {}

    /** Parse user_input for transition requests, calculate and print.
        Usually this will not need to be overridden.
     */
    virtual void CalculateAndPrint();

    /** Print all stored transitions.
        Usually this will not need to be overridden.
     */
    virtual void PrintAll() const;

    /** Print one-electron transition integrals. */
    virtual void PrintIntegrals();

    pSpinorMatrixElementConst GetOperator() const { return op; }
protected:
    /** Print the header line to outstream, explaining transition type, units, etc. */
    virtual void PrintHeader() const = 0;

    /** Print a line for the transition to outstream. */
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const = 0;

protected:
    /** Parse user_input for RPA options. Return an RPAOperator.
        PRE: hf->GetCore() should be a self-consistent solution of hf.
     */
    virtual pRPAOperator MakeRPA(pSpinorOperator external, pHFOperatorConst hf, pHartreeY hartreeY);

    /** Check if transition exists or is excluded by symmetry considerations. */
    inline bool TransitionExists(const Symmetry& left, const Symmetry& right) const
    {
        return ((abs(left.GetTwoJ() - right.GetTwoJ()) <= 2 * op->GetK()) &&
                (left.GetTwoJ() + right.GetTwoJ() >= 2 * op->GetK()) &&
                (left.GetParity() * right.GetParity() == op->GetParity()));
    }

    inline bool TransitionExists(const LevelID& left, const LevelID& right) const
    {
        return TransitionExists(left.first->GetSymmetry(), right.first->GetSymmetry());
    }

    /** Calculate transition matrix element and add to stored transitions.
        While we're at it, calculate all transitions with same symmetry, since this hardly costs any more.
        Return transition that was requested.
     */
    double CalculateTransition(const LevelID& left, const LevelID& right);

    /** Parse transition string and calculate matrix element. Transition string is of the form
            "HamiltonianID:i -> HamiltonianID:j"
        or simply "HamiltonianID:i" (meaning to itself). Whitespace ignored.
     */
    double CalculateTransition(const std::string& transition);

    /** Convert string to levelID. If name is not well-formed, return
            std::pair<nullptr, 0>
     */
    LevelID make_LevelID(const std::string& name);

    /** TransitionID with left <= right. */
    TransitionID make_transitionID(const LevelID& left, const LevelID& right) const
    {
        if(right < left)
            return std::make_pair(right, left);
        else
            return std::make_pair(left, right);
    }

protected:
    MultirunOptions& user_input;
    pOrbitalManagerConst orbitals;
    pLevelStore levels;
    pSpinorMatrixElement op = nullptr;

    double scale;
    pTransitionIntegrals integrals;     // Scaled one-electron integrals
    std::map<TransitionID, double> matrix_elements; // Not scaled
};

}
#endif
