#ifndef RPA_OPERATOR_H
#define RPA_OPERATOR_H

#include "HartreeFock/SpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "RPAOrbital.h"

/** Add RPA corrections to external SpinorOperator. */
class RPAOperator : public SpinorOperator
{
public:
    RPAOperator(pSpinorOperator external, pHartreeY hartreeY):
        SpinorOperator(external->GetK(), external->GetParity(), external->GetIntegrator()),
        external(external), hartreeY(hartreeY), core(nullptr)
    {}

    /** Set core RPA orbitals. */
    virtual void SetCore(pCoreConst rpa_core) { core = rpa_core; }

    /** Return (f + deltaVhf)||a> */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pSpinorOperator external;
    pCoreConst core;
    pHartreeY hartreeY;
};

typedef std::shared_ptr<RPAOperator> pRPAOperator;
typedef std::shared_ptr<const RPAOperator> pRPAOperatorConst;

inline SpinorFunction RPAOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
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
    ret = external->ReducedApplyTo(a, kappa_alph);

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
                    double coeff = math->Wigner6j(K, alph_J, a.J(), k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k);
                    coeff *= math->Electron3j(a.TwoJ(), beta->TwoJ(), k) * math->Electron3j(alph_TwoJ, b->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, beta, pa);
                        ret += hartreeY->ApplyTo(*b, ret.Kappa()) * coeff;
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
                    double coeff = math->Wigner6j(K, a.J(), alph_J, k, beta->J(), b->J()) *
                                   math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k);
                    coeff *= math->Electron3j(a.TwoJ(), b->TwoJ(), k) * math->Electron3j(alph_TwoJ, beta->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, b, pa);
                        ret += hartreeY->ApplyTo(*beta, ret.Kappa()) * coeff;
                    }
                }
            }
        }
    }
    
    return ret;
}

#endif
