#ifndef LORENTZ_INVARIANCE_T2_H
#define LORENTZ_INVARIANCE_T2_H

#include "HartreeFock/SpinorOperator.h"
#include "Transitions.h"

namespace Ambit
{
/** Lorentz invariance violation operator T2 based on relativistic formula
        \f$ T^{(2)} = c\,(\vec{\alpha}\cdot\vec{p} - \alpha_3 p_3)\f$.
    Detailed equations may be found in Hohensee et al. PRL 111, 050401 (2013), Supplemental Material.
 */
class LorentzInvarianceT2Operator: public SpinorOperator
{
public:
    LorentzInvarianceT2Operator(pIntegrator integration_strategy):
        SpinorOperator(2, Parity::even, integration_strategy), lattice(integration_strategy->GetLattice())
    {
        differentiator = std::make_shared<FornbergDifferentiator>(lattice, 7, false);
    }

public:
    /** LLI T2 in atomic units. */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pLattice lattice;
    pFornbergDifferentiator differentiator;
};

class LorentzInvarianceT2Calculator : public TransitionCalculator
{
public:
    LorentzInvarianceT2Calculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;
};

}

#endif
