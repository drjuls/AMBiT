#ifndef RPA_SOLVER_H
#define RPA_SOLVER_H

#include "HartreeFock/Core.h"
#include "HartreeFock/HFOperator.h"
#include "Basis/BSplineBasis.h"
#include "RPAOrbital.h"

namespace Ambit
{
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
        basis_maker(basis_maker), include_dirac_sea(include_negative_basis), hf0(nullptr)
    {}
    RPASolver(pLattice lattice, double rmax = 50., bool include_negative_basis = true):
        include_dirac_sea(include_negative_basis), hf0(nullptr)
    {
        basis_maker = std::make_shared<BSplineBasis>(lattice, 40, 7, rmax);
    }
    RPASolver(pCoreConst core, bool include_negative_basis = true):
        RPASolver(core->GetLattice(), core->GetLattice()->R(core->LargestOrbitalSize()), include_negative_basis)
    {}
    virtual ~RPASolver() = default;

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
    void CreateDeltaOrbitals(pRPAOrbital orbital, int operator_K, Parity operator_P, bool operator_static) const;

    /** Iterate a single deltaOrbital.
        Return change in its deltaEnergy.
     */
    double IterateDeltaOrbital(pDeltaOrbital orbital, std::shared_ptr<const RPAOperator> rpa, double propnew) const;

    /** Iterate a pair of deltaOrbitals (for frequency-dependent RPA).
        Return change in deltaEnergy of first orbital.
     */
    double IterateDeltaOrbital(std::pair<pDeltaOrbital, pDeltaOrbital>& orbitals, std::shared_ptr<const RPAOperator> rpa, double propnew) const;

    /** Set the weighting of successive TDHF loops.
        PRE: 0 < propnew <= 1
     */
    void SetTDHFWeighting(double propnew) { TDHF_propnew = propnew; }

    /** Set limit on the maximum angular momentum difference between each DeltaOrbital and its parent.
        PRE: limit_deltaK >= 0
     */
    void SetLimitDeltaK(int deltaK) { limit_deltaK = deltaK; }

    /** Set limit on the maximum number of iterations to RPA convergence.
        PRE: max_iterations > 0
     */
    void SetMaxRPAIterations(int max_iterations) { MaxRPAIterations = max_iterations; }

protected:
    pBSplineBasis basis_maker;
    pHFOperatorConst hf0;               //!< Keep HF operator for making additional basis orbitals
    bool include_dirac_sea;
    std::map<int, pOrbitalMap> basis;   //!< DeltaOrbital basis for each kappa
    int limit_deltaK = -1;              //!< Limits on maximum angular momentum difference between each DeltaOrbital and its parent.

    double TDHF_propnew = 0.5;          //!< Weighting to apply to each iteration
    double EnergyTolerance = 1.e-14;
    unsigned int MaxRPAIterations = 100;
};

typedef std::shared_ptr<RPASolver> pRPASolver;

}
#endif
