#ifndef BASIS_GENERATOR_H
#define BASIS_GENERATOR_H

#include "Atom/MultirunOptions.h"
#include "HartreeFock/Orbital.h"
#include "HartreeFock/HFOperator.h"
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
    BasisGenerator(pLattice lat, MultirunOptions& userInput);
    virtual ~BasisGenerator();

    /** Generate core orbitals. If open_shell_core is supplied, then use this as a starting approximation.
        Return open shell Hartree-Fock core.
     */
    virtual pCore GenerateHFCore(pCoreConst open_shell_core = pCoreConst());

    /** Create open shell core and Hartree-Fock operator using existing orbitals.
        This function will not change orbital_manager, but later non-const functions might.
     */
    virtual pHFOperator CreateHFOperator(pOrbitalManager orbital_manager);

    /** Generate excited states.
        PRE: core must have been built using GenerateHFCore() already.
     */
    virtual pOrbitalManagerConst GenerateBasis();

    virtual pHFOperatorConst GetHFOperator() const { return hf; }
    virtual pCoreConst GetHFCore() const { return open_core; }
    virtual pOrbitalManagerConst GetOrbitals() const { return orbitals; }

protected:
    /** Create hf operator and set open_core occupancies. Used by GenerateHFCore() and CreateHFOperator().
        POST: this->hf is dressed HF operator;
              undressed_hf is base HF with finite nuclear radius;
              this->open_core has correct occupancies.
     */
    virtual void InitialiseHF(pHFOperator& undressed_hf);

    /** When generating excited states, assume that states below Fermi level in "orbitals" are good to go. */
    virtual pOrbitalMap GenerateBSplines(const std::vector<int>& max_pqn);

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
    pHFOperator hf;
    pCore open_core;

    pOrbitalManager orbitals;
};

#endif
