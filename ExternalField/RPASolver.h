#ifndef RPA_SOLVER_H
#define RPA_SOLVER_H

#include "HartreeFock/Core.h"
#include "HartreeFock/HFOperator.h"
#include "Basis/BSplineBasis.h"
#include "RPAOrbital.h"

class RPAOperator;

/** Solve RPA equations self-consistently:
        (H_0 - E_a)|alpha> = -(f + deltaV + deltaE)|a>
    where |alpha> is a small perturbation to |a>.
    The RHS is evaluated using RPA operator.
 */
class RPASolver
{
public:
    RPASolver(pBSplineBasis basis_maker, bool include_negative_basis = true):
        basis_maker(basis_maker), remove_spin_orbit_partner(false), include_dirac_sea(include_negative_basis), hf0(nullptr)
    {}
    RPASolver(pLattice lattice, double rmax = 50., bool include_negative_basis = true):
        remove_spin_orbit_partner(false), include_dirac_sea(include_negative_basis), hf0(nullptr)
    {
        basis_maker = std::make_shared<BSplineBasis>(lattice, 40, 7, rmax);
    }
    RPASolver(pCoreConst core, bool include_negative_basis = true):
        RPASolver(core->GetLattice(), core->GetLattice()->R(core->LargestOrbitalSize()), include_negative_basis)
    {}

    ~RPASolver() {}

    /** Solve RPA equations for core by expanding deltaOrbitals over a set of B-splines and
        taking matrix elements of RPAOperator to determine coefficients.
        If include_negative_basis, include basis states in Dirac sea.
        PRE: core should be self-consistent solution of hf.
     */
    void SolveRPACore(pHFOperatorConst hf, std::shared_ptr<RPAOperator> rpa);

    /** Return RPA energy correction to excited state. */
    double CalculateRPAExcited(pRPAOrbital orbital, std::shared_ptr<const RPAOperator> rpa);

    /** Add all appropriate deltaOrbitals to RPA orbital. */
    void CreateDeltaOrbitals(pRPAOrbital orbital, std::shared_ptr<const RPAOperator> rpa) const;

    /** Iterate a single deltaOrbital.
        Return change in its deltaEnergy.
     */
    double IterateDeltaOrbital(pDeltaOrbital orbital, std::shared_ptr<const RPAOperator> rpa) const;

    double EnergyTolerance = 1.e-14;

protected:
    pBSplineBasis basis_maker;
    pHFOperatorConst hf0;               //!< Keep HF operator for making additional basis orbitals
    bool include_dirac_sea;
    std::map<int, pOrbitalMap> basis;   //!< DeltaOrbital basis for each kappa

    bool remove_spin_orbit_partner;     //!< Remove spin-orbit partner from basis (small energy denominators can be fatal)
    unsigned int MaxRPAIterations = 100;
};

typedef std::shared_ptr<RPASolver> pRPASolver;

#endif
