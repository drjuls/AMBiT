#include "Hyperfine.h"
#include "Include.h"
#include "Universal/MathConstant.h"

HyperfineDipoleOperator::HyperfineDipoleOperator(pOPIntegrator integration_strategy):
    SpinorOperator(1, integration_strategy), lattice(integration_strategy->GetLattice())
{}

SpinorFunction HyperfineDipoleOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b, a.size());
    MathConstant* math = MathConstant::Instance();

    double prefactor = -(a.Kappa() + kappa_b) * math->SphericalTensorReducedMatrixElement(-kappa_b, a.Kappa(), 1);
    if(!prefactor)
        return ret;
    prefactor *= math->NuclearMagneton();

    const double* R2 = lattice->Rpower(2);
    const double* R3 = lattice->Rpower(3);

    for(unsigned int i = 0; i < ret.size(); i++)
    {
        ret.f[i] = a.g[i]/R2[i];
        ret.dfdr[i] = a.dgdr[i]/R2[i] - 2. * a.g[i]/R3[i];
        ret.g[i] = a.f[i]/R2[i];
        ret.dgdr[i] = a.dfdr[i]/R2[i] - 2. * a.g[i]/R3[i];
    }

    return ret * prefactor;
}

SpinorFunction HyperfineRPAOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b);
    MathConstant* math = MathConstant::Instance();

    // To save confusion: kappa_b -> alph
    int kappa_alph = kappa_b;
    int alph_L = (kappa_b > 0? kappa_b: -kappa_b-1);
    int alph_TwoJ = 2 * abs(kappa_alph) - 1;
    double alph_J = double(alph_TwoJ)/2;

    pSpinorFunctionConst pa = a.shared_from_this();

    // First term: <kappa || f || parent>
    ret = hyperfine.ApplyTo(a, kappa_alph);

    // Second term: <kappa || deltaV || parent>
    // Sum deltaV for all core states. Core states are `b', delta is `beta'.
    for(const auto& cs: *core)
    {
        pRPAOrbitalConst b(std::dynamic_pointer_cast<const RPAOrbital>(cs.second));
        if(b == nullptr)
            continue;

        // Sum over beta
        for(const auto& cs_delta: b->deltapsi)
        {
            pDeltaOrbitalConst beta = cs_delta.second;
            double occupancy_factor = (a.TwoJ()+1) * (b->TwoJ()+1) * (alph_TwoJ+1) * (beta->TwoJ()+1);
            occupancy_factor = sqrt(occupancy_factor);
            // TODO: Open shells need to be scaled

            // C: alph -> b; beta -> a
            // Sum over all k
            if((alph_L + b->L())%2 == (a.L() + beta->L())%2)
            {
                int mink = mmax(abs(a.L() - beta->L()), abs(alph_L - b->L()));
                int maxk = mmin(a.L() + beta->L(), alph_L + b->L());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(1., alph_J, a.J(), k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k + 1);
                    coeff *= math->Electron3j(a.TwoJ(), beta->TwoJ(), k) * math->Electron3j(alph_TwoJ, b->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, beta, pa);
                        ret -= hartreeY->ApplyTo(*b) * coeff;
                    }
                }
            }

            // D: alph->beta; b -> a
            // Sum over all k
            if((a.L() + b->L())%2 == (alph_L + beta->L())%2)
            {
                int mink = mmax(abs(a.L() - b->L()), abs(alph_L - beta->L()));
                int maxk = mmin(a.L() + b->L(), alph_L + beta->L());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(1., a.J(), alph_J, k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k + 1);
                    coeff *= math->Electron3j(a.TwoJ(), b->TwoJ(), k) * math->Electron3j(alph_TwoJ, beta->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, b, pa);
                        ret -= hartreeY->ApplyTo(*beta) * coeff;
                    }
                }
            }
        }
    }

    return ret;
}

HyperfineDipoleRPADecorator::HyperfineDipoleRPADecorator(pHFOperator wrapped, pHartreeY hartreeY, pOPIntegrator integration_strategy):
    BaseDecorator(wrapped, integration_strategy), hyperfine(integration_strategy), hartreeY(hartreeY)
{}

void HyperfineDipoleRPADecorator::Alert()
{
    if(currentExchangePotential.size() > lattice->size())
        currentExchangePotential.resize(lattice->size());
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void HyperfineDipoleRPADecorator::SetODEParameters(const Orbital& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction HyperfineDipoleRPADecorator::GetExchange(pOrbitalConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == NULL)
        ret += currentExchangePotential;
    else
        ret += CalculateExtraExchange(*approximation);

    return ret;
}

void HyperfineDipoleRPADecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void HyperfineDipoleRPADecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void HyperfineDipoleRPADecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        dwdr[0] += alpha * currentExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * currentExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction HyperfineDipoleRPADecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= CalculateExtraExchange(a);

    return ta;
}

SpinorFunction HyperfineDipoleRPADecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa());

    // Only add additional RHS to deltaOrbital. Call deltaOrbital `alph', parent is `a'.
    const DeltaOrbital* alph = dynamic_cast<const DeltaOrbital*>(&s);
    if(alph == nullptr)
        return exchange;

    MathConstant* math = MathConstant::Instance();
    pRPAOrbital a(alph->GetParent());

    // First term: - <kappa || f || parent>
    exchange = hyperfine.ApplyTo(*a, alph->Kappa()) * (-1.);

    // Second term: <kappa || deltaV || parent>
    // Sum deltaV for all core states. Core states are `b', delta is `beta'.
    for(const auto& cs: *core)
    {
        pRPAOrbitalConst b(std::dynamic_pointer_cast<const RPAOrbital>(cs.second));
        if(b == nullptr)
            continue;

        // Sum over beta
        for(const auto& cs_delta: b->deltapsi)
        {
            pDeltaOrbitalConst beta = cs_delta.second;
            double occupancy_factor = (a->TwoJ()+1) * (b->TwoJ()+1) * (alph->TwoJ()+1) * (beta->TwoJ()+1);
            occupancy_factor = sqrt(occupancy_factor);
            // TODO: Open shells need to be scaled

            // C: alph -> b; beta -> a
            // Sum over all k
            if((alph->L() + b->L())%2 == (a->L() + beta->L())%2)
            {
                int mink = mmax(abs(a->L() - beta->L()), abs(alph->L() - b->L()));
                int maxk = mmin(a->L() + beta->L(), alph->L() + b->L());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(1., alph->J(), a->J(), k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a->TwoJ() + b->TwoJ())/2 + k + 1);
                    coeff *= math->Electron3j(a->TwoJ(), beta->TwoJ(), k) * math->Electron3j(alph->TwoJ(), b->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, beta, a);
                        exchange += hartreeY->ApplyTo(*b) * coeff;
                    }
                }
            }

            // D: alph->beta; b -> a
            // Sum over all k
            if((a->L() + b->L())%2 == (alph->L() + beta->L())%2)
            {
                int mink = mmax(abs(a->L() - b->L()), abs(alph->L() - beta->L()));
                int maxk = mmin(a->L() + b->L(), alph->L() + beta->L());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(1., a->J(), alph->J(), k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a->TwoJ() + b->TwoJ())/2 + k + 1);
                    coeff *= math->Electron3j(a->TwoJ(), b->TwoJ(), k) * math->Electron3j(alph->TwoJ(), beta->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, b, a);
                        exchange += hartreeY->ApplyTo(*beta) * coeff;
                    }
                }
            }
        }
    }

    // Third term: deltaE |parent >
    if(a->Kappa() == alph->Kappa())
        exchange += (*a) * alph->DeltaEnergy();

    return exchange;
}
