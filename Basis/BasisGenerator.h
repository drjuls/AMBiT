#ifndef BASIS_GENERATOR_H
#define BASIS_GENERATOR_H

#include "Atom/MultirunOptions.h"
#include "HartreeFock/Orbital.h"
#include "HartreeFock/HFOperator.h"
#include "HartreeFock/HartreeY.h"
#include "OrbitalManager.h"
#include <list>

/** Based on user input:
    - create HartreeFock operator (with relevant decorators)
    - create (or read) core orbitals, solving the Hartree Fock equations as required
    - create (or read) excited states
 */
class BasisGenerator
{
public:
    BasisGenerator(pLattice lat, MultirunOptions& userInput, pPhysicalConstant physical_constant = nullptr);
    virtual ~BasisGenerator();

    /** Generate core orbitals. If open_shell_core is supplied, then use this as a starting approximation.
        Return open shell Hartree-Fock core.
     */
    virtual pCore GenerateHFCore(pCoreConst open_shell_core = pCoreConst());

    /** Generate excited states.
        PRE: core must have been built using GenerateHFCore() already.
     */
    virtual pOrbitalManagerConst GenerateBasis();

    /** Build open-shell core, Hartree-Fock operator, and basis using read orbitals and the input file.
        This function will not change orbital_manager->all, but may change other OrbitalMaps.
        Returns open-shell HF operator.
     */
    virtual pHFOperator RecreateBasis(pOrbitalManager orbital_manager);

    /** Get open-shell Hartree-Fock operator. */
    virtual pHFOperatorConst GetOpenHFOperator() const { return hf; }

    /** Get open-shell Hartree-Fock operator.
        This can be converted to a closed-shell HF operator using hf->SetCore().
     */
    virtual pHFOperator GetOpenHFOperator() { return hf; }

    /** Return clone of Hartree-Fock operator with closed-shell core. */
    virtual pHFOperator GetClosedHFOperator()
    {   pHFOperator ret(hf->Clone());
        ret->SetCore(closed_core);
        return ret;
    }

    /** Get open-shell core. */
    virtual pCoreConst GetHFCore() const { return open_core; }
    virtual pCore GetHFCore() { return open_core; }

    virtual pOrbitalManagerConst GetBasis() const { return orbitals; }

    /** Get HartreeY operator. */
    virtual pHartreeY GetHartreeY() { return hartreeY; }
    virtual pHartreeYConst GetHartreeY() const { return hartreeY; }

    /** Get Physical constants. */
    virtual pPhysicalConstant GetPhysicalConstant() { return physical_constant; }
    virtual pPhysicalConstantConst GetPhysicalConstant() const { return physical_constant; }

protected:
    /** Create open-shell hf operator and set open_core occupancies. Used by GenerateHFCore() and RecreateBasis().
        POST: undressed_hf is base HF with finite nuclear radius;
              this->hf is dressed HF operator;
              this->hartreeY is dressed HartreeY operator;
              this->open_core has correct occupancies.
     */
    virtual void InitialiseHF(pHFOperator& undressed_hf);

    /** Set orbital maps core, valence, etc according to user input.
        PRE:  orbitals->all and indexes must be up to date and include all states.
        POST: No change in orbitals->all or the indexes.
     */
    virtual void SetOrbitalMaps();

    /** When generating excited states, assume that states below Fermi level in "orbitals" are good to go. */
    virtual pOrbitalMap GenerateBSplines(const std::vector<int>& max_pqn);

    /** When generating excited states, assume that states below Fermi level in "orbitals" are good to go. */
    virtual pOrbitalMap GenerateHFExcited(const std::vector<int>& max_pqn);

    /** Orthogonalise to all states that have the same kappa and principal quantum number
        less than current (both in core and in excited states).
        PRE: all states with smaller pqn must already be orthogonal
     */
    virtual void Orthogonalise(pOrbital orbital) const;

    /** Test for orthogonality of states.
        Return largest overlap and set max_i and max_j to states that give it.
     */
    double TestOrthogonality(OrbitalInfo& max_i, OrbitalInfo& max_j) const;

protected:
    pLattice lattice;
    MultirunOptions& user_input;
    pPhysicalConstant physical_constant;

    pCore open_core;
    pCore closed_core;
    pOrbitalManager orbitals;

    pHFOperator hf;     //!< Open-shell hf operator
    pHartreeY hartreeY; //!< Dressed HartreeY operator
};

#endif
